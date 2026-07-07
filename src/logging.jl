using Dates

# =============================================================================
# Abstract Logger Interface
# =============================================================================

"""
    AbstractProcessingLogger

Abstract base type for processing loggers. Enables type-stable dispatch
and zero-cost abstraction when logging is disabled.

Subtypes:
- `ProcessingLogger`: Active logger that records events
- `NoOpLogger`: Singleton that compiles away to no-ops
"""
abstract type AbstractProcessingLogger end

# =============================================================================
# NoOpLogger - Zero-cost disabled logging
# =============================================================================

"""
    NoOpLogger()

Singleton logger that performs no operations. All logging calls compile
to no-ops, providing zero runtime overhead when logging is disabled.

# Example
```julia
logger = NoOpLogger()
log_event!(logger, :qc, :bounds)  # Compiles to nothing
```
"""
struct NoOpLogger <: AbstractProcessingLogger end

@inline log_event!(::NoOpLogger, ::Symbol, ::Symbol; kwargs...) = nothing
@inline record_stage_time!(::NoOpLogger, ::Symbol, ::Real) = nothing
@inline write_processing_log(::NoOpLogger, ::AbstractString; append::Bool=false) = nothing
@inline log_index_runs!(::NoOpLogger, ::Symbol, ::Symbol, ::Union{Symbol,Nothing},
                        ::AbstractVector, ::AbstractVector{Int}; kwargs...) = nothing
@inline log_mask_runs!(::NoOpLogger, ::Symbol, ::Symbol, ::Union{Symbol,Nothing},
                       ::AbstractVector, ::AbstractVector{Bool}; kwargs...) = nothing
"""
    is_logging_enabled(logger::AbstractProcessingLogger) -> Bool

Return `true` if `logger` is an active logger that records events, and `false` for
`NoOpLogger`.
"""
function is_logging_enabled end

@inline is_logging_enabled(::NoOpLogger) = false

# =============================================================================
# ProcessingLogEntry - Immutable log record
# =============================================================================

"""
    ProcessingLogEntry

Immutable record of a single processing event.

# Fields
- `stage::Symbol`: Pipeline stage (e.g., `:quality_control`, `:despiking`)
- `category::Symbol`: Event category (e.g., `:bounds`, `:spike`)
- `variable::Union{Symbol,Nothing}`: Affected variable, if applicable
- `start_time::Union{DateTime,Nothing}`: Event start timestamp
- `end_time::Union{DateTime,Nothing}`: Event end timestamp
- `details::Dict{Symbol,Any}`: Additional metadata
"""
struct ProcessingLogEntry
    stage::Symbol
    category::Symbol
    variable::Union{Symbol,Nothing}
    start_time::Union{DateTime,Nothing}
    end_time::Union{DateTime,Nothing}
    details::Dict{Symbol,Any}
end

# =============================================================================
# ProcessingLogger - Active event logger
# =============================================================================

"""
    ProcessingLogger()

Mutable logger that accumulates processing events and stage durations.
Events are stored in memory and can be written to CSV via `write_processing_log`.

# Example
```julia
logger = ProcessingLogger()
log_event!(logger, :qc, :bounds; variable=:Ux, start_time=now())
record_stage_time!(logger, :qc, 1.5)
write_processing_log(logger, "processing.csv")
```
"""
mutable struct ProcessingLogger <: AbstractProcessingLogger
    entries::Vector{ProcessingLogEntry}
    stage_durations::Dict{Symbol,Float64}
    metadata::Dict{String,String}
end

ProcessingLogger() = ProcessingLogger(ProcessingLogEntry[], Dict{Symbol,Float64}(), Dict{String,String}())

@inline is_logging_enabled(::ProcessingLogger) = true

"""
    log_metadata!(logger, key, value)

Store a metadata key-value pair. Written as `# key = value` comment lines by `write_processing_log`.
"""
function log_metadata!(logger::ProcessingLogger, key::String, value::String)
    logger.metadata[key] = value
    return nothing
end
@inline log_metadata!(::NoOpLogger, ::String, ::String) = nothing

"""
    reset!(logger::ProcessingLogger)

Clear entries and stage durations while preserving metadata. Useful for reusing
a logger across multiple days within the same month.
"""
function reset!(logger::ProcessingLogger)
    empty!(logger.entries)
    empty!(logger.stage_durations)
    return nothing
end

# =============================================================================
# Core logging functions
# =============================================================================

"""
    log_event!(logger, stage, category; variable=nothing, start_time=nothing, end_time=nothing, kwargs...)

Record a processing event. Duration is automatically computed if both timestamps are provided.
"""
function log_event!(logger::ProcessingLogger, stage::Symbol, category::Symbol;
                    variable::Union{Symbol,Nothing}=nothing,
                    start_time::Union{DateTime,Nothing}=nothing,
                    end_time::Union{DateTime,Nothing}=nothing,
                    kwargs...)
    # Normalize end_time
    t_end = _normalize_end_time(start_time, end_time)
    
    # Build details dict and compute duration
    details = _build_details(start_time, t_end, kwargs)
    
    entry = ProcessingLogEntry(stage, category, variable, start_time, t_end, details)
    push!(logger.entries, entry)
    return nothing
end

"""
    record_stage_time!(logger, stage, seconds)

Accumulate runtime for a pipeline stage. Multiple calls for the same stage are summed.
"""
function record_stage_time!(logger::ProcessingLogger, stage::Symbol, seconds::Real)
    prev = get(logger.stage_durations, stage, 0.0)
    logger.stage_durations[stage] = prev + Float64(seconds)
    return nothing
end

"""
    write_processing_log(logger, filepath; append=false)

Write all logged events and stage durations to a CSV file.
When `append=true`, entries are appended without re-writing the metadata header.
"""
function write_processing_log(logger::ProcessingLogger, filepath::AbstractString; append::Bool=false)
    open(filepath, append ? "a" : "w") do io
        if !append
            _write_metadata(io, logger.metadata)
            _write_header(io)
        end
        _write_entries(io, logger.entries)
        _write_stage_durations(io, logger.stage_durations)
    end
    return nothing
end

# =============================================================================
# Batch logging utilities
# =============================================================================

"""
    log_index_runs!(logger, stage, category, variable, timestamps, indices; include_run_length=false, kwargs...)

Log contiguous runs of indices as separate events. Useful for logging flagged data points.
"""
function log_index_runs!(logger::ProcessingLogger, stage::Symbol, category::Symbol,
                         variable::Union{Symbol,Nothing}, timestamps::AbstractVector,
                         indices::AbstractVector{Int}; include_run_length::Bool=false,
                         kwargs...)
    isempty(indices) && return nothing
    
    sorted = sort(indices)
    n = length(sorted)
    run_start_idx = 1
    
    for i in 2:n
        if sorted[i] != sorted[i-1] + 1
            _log_run!(logger, stage, category, variable, timestamps,
                      sorted[run_start_idx], sorted[i-1], include_run_length; kwargs...)
            run_start_idx = i
        end
    end
    
    _log_run!(logger, stage, category, variable, timestamps,
              sorted[run_start_idx], sorted[n], include_run_length; kwargs...)
    
    return nothing
end

"""
    log_mask_runs!(logger, stage, category, variable, timestamps, mask; kwargs...)

Log contiguous runs of `true` values in a boolean mask as separate events.
"""
function log_mask_runs!(logger::ProcessingLogger, stage::Symbol, category::Symbol,
                        variable::Union{Symbol,Nothing}, timestamps::AbstractVector,
                        mask::AbstractVector{Bool}; kwargs...)
    n = min(length(timestamps), length(mask))
    n == 0 && return nothing
    
    run_start = 0  # 0 = not in a run (type-stable)
    
    for i in 1:n
        if mask[i]
            run_start == 0 && (run_start = i)
        elseif run_start != 0
            log_event!(logger, stage, category; variable=variable,
                       start_time=timestamps[run_start], end_time=timestamps[i-1], kwargs...)
            run_start = 0
        end
    end
    
    if run_start != 0
        log_event!(logger, stage, category; variable=variable,
                   start_time=timestamps[run_start], end_time=timestamps[n], kwargs...)
    end
    
    return nothing
end

# =============================================================================
# Internal helpers - Event building
# =============================================================================

@inline function _normalize_end_time(start_time::Union{DateTime,Nothing},
                                     end_time::Union{DateTime,Nothing})::Union{DateTime,Nothing}
    end_time !== nothing && return end_time
    return start_time  # Default to start_time if end_time not provided
end

function _build_details(start_time::Union{DateTime,Nothing},
                        end_time::Union{DateTime,Nothing},
                        kwargs)::Dict{Symbol,Any}
    details = Dict{Symbol,Any}(kwargs)
    
    if start_time !== nothing && end_time !== nothing && !haskey(details, :duration_seconds)
        details[:duration_seconds] = Dates.value(end_time - start_time) / 1000.0
    end
    
    return details
end

@inline function _log_run!(logger::ProcessingLogger, stage::Symbol, category::Symbol,
                           variable::Union{Symbol,Nothing}, timestamps::AbstractVector,
                           first_data_idx::Int, last_data_idx::Int,
                           include_run_length::Bool; kwargs...)
    if include_run_length
        run_len = last_data_idx - first_data_idx + 1
        log_event!(logger, stage, category; variable=variable,
                   start_time=timestamps[first_data_idx], end_time=timestamps[last_data_idx],
                   samples_in_run=run_len, kwargs...)
    else
        log_event!(logger, stage, category; variable=variable,
                   start_time=timestamps[first_data_idx], end_time=timestamps[last_data_idx],
                   kwargs...)
    end
    return nothing
end

# =============================================================================
# Internal helpers - CSV formatting
# =============================================================================

const _CSV_HEADER = "stage,category,variable,start_timestamp,end_timestamp,duration_seconds,details"
const _TIME_FORMAT = dateformat"yyyy-mm-ddTHH:MM:SS.s"

function _write_metadata(io::IO, metadata::Dict{String,String})
    for (k, v) in sort(collect(metadata))
        println(io, "# ", k, " = ", v)
    end
end

@inline _write_header(io::IO) = println(io, _CSV_HEADER)

function _write_entries(io::IO, entries::Vector{ProcessingLogEntry})
    for entry in entries
        _write_entry(io, entry)
    end
end

function _write_stage_durations(io::IO, durations::Dict{Symbol,Float64})
    for (stage, seconds) in sort!(collect(durations); by=first)
        println(io, stage, ",runtime,,,,", seconds, ",duration_seconds=", seconds)
    end
end

function _write_entry(io::IO, entry::ProcessingLogEntry)
    print(io, entry.stage, ",")
    print(io, entry.category, ",")
    print(io, _format_variable(entry.variable), ",")
    print(io, _format_time(entry.start_time), ",")
    print(io, _format_time(entry.end_time), ",")
    
    duration_str, details_str = _format_details(entry.details)
    print(io, duration_str, ",")
    println(io, details_str)
end

@inline _format_variable(v::Nothing)::String = ""
@inline _format_variable(v::Symbol)::String = string(v)

@inline _format_time(t::Nothing)::String = ""
@inline _format_time(t::DateTime)::String = Dates.format(t, _TIME_FORMAT)

function _format_details(details::Dict{Symbol,Any})::Tuple{String,String}
    isempty(details) && return "", ""
    
    duration_str = ""
    if haskey(details, :duration_seconds)
        duration_str = string(details[:duration_seconds])
    end
    
    other_keys = filter(k -> k !== :duration_seconds, collect(keys(details)))
    if isempty(other_keys)
        return duration_str, ""
    end
    
    sort!(other_keys; by=string)
    parts = [string(k, "=", details[k]) for k in other_keys]
    return duration_str, join(parts, ";")
end
