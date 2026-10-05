export EddyPipeline
export process!
export check_data

using ProgressMeter

"""
    EddyPipeline(; sensor, quality_control, gas_analyzer, despiking, make_continuous, gap_filling,
                  double_rotation, mrd, output, logger)

High-level orchestrator for the Peddy.jl processing pipeline.

Each step is a pluggable component that implements its respective abstract
interface. Any step can be set to `nothing` to skip it. The typical order is:

1. Quality control (`quality_control!`)
2. Gas analyzer correction (`correct_gas_analyzer!`)
3. Despiking (`despike!`)
4. Make continuous (`make_continuous!`) – insert missing timestamps with NaNs
5. Gap filling / interpolation (`fill_gaps!`)
6. Double rotation (`rotate!`)
7. MRD decomposition (`decompose!`)
8. Output writing (`write_data`)

Fields
- `sensor::AbstractSensor`: Sensor providing metadata (and coefficients)
- `quality_control::OptionalPipelineStep`
- `gas_analyzer::OptionalPipelineStep`
- `despiking::OptionalPipelineStep`
- `make_continuous::OptionalPipelineStep`
- `gap_filling::OptionalPipelineStep`
- `double_rotation::OptionalPipelineStep`
- `mrd::OptionalPipelineStep`
- `output::AbstractOutput`: Writer implementation
- `logger::AbstractProcessingLogger`: Logger for events and timing (default: `NoOpLogger()`)
"""
@kwdef struct EddyPipeline{SI<:AbstractSensor,QC<:OptionalPipelineStep,
                           D<:OptionalPipelineStep,MC<:OptionalPipelineStep,G<:OptionalPipelineStep,
                           GA<:OptionalPipelineStep,DR<:OptionalPipelineStep,
                           MRD<:OptionalPipelineStep,O<:AbstractOutput,
                           L<:AbstractProcessingLogger}
    sensor::SI
    quality_control::QC = nothing
    gas_analyzer::GA = nothing
    despiking::D = nothing
    make_continuous::MC = nothing
    gap_filling::G = nothing
    double_rotation::DR = nothing
    mrd::MRD = nothing
    output::O
    logger::L = NoOpLogger()
end

# Default to no operation
"""No-op quality control when `quality_control` is `nothing`."""
function quality_control!(qc::Nothing, high_frequency_data, low_frequency_data, sensor;
                          kwargs...)
    return nothing
end
"""No-op gas analyzer correction when `gas_analyzer` is `nothing`."""
function correct_gas_analyzer!(gas_analyzer::Nothing, high_frequency_data,
                               low_frequency_data, sensor::AbstractSensor; kwargs...)
    return nothing
end
"""No-op despiking when `despiking` is `nothing`."""
function despike!(despiking::Nothing, high_frequency_data, low_frequency_data; kwargs...)
    return nothing
end
"""No-op gap filling when `gap_filling` is `nothing`."""
function fill_gaps!(gap_filling::Nothing, high_frequency_data, low_frequency_data;
                    kwargs...)
    return nothing
end
"""No-op make_continuous when `make_continuous` is `nothing`."""
function make_continuous!(make_continuous::Nothing, high_frequency_data, low_frequency_data;
                          kwargs...)
    return high_frequency_data
end
"""No-op double rotation when `double_rotation` is `nothing`."""
function rotate!(double_rotation::Nothing, high_frequency_data, low_frequency_data;
                 kwargs...)
    return nothing
end
decompose!(mrd::Nothing, high_frequency_data, low_frequency_data; kwargs...) = nothing

# Data should be in the correct format
"""
    process!(pipeline::EddyPipeline, high_frequency_data::DimArray,
             low_frequency_data::Union{Nothing,DimArray}; kwargs...) -> Nothing

Run the configured pipeline over the provided data. Steps that are `nothing`
are skipped. Progress is shown via a spinner with status messages.

Arguments
- `pipeline`: An `EddyPipeline` instance
- `high_frequency_data`: DimArray with fast measurements (must have `Var` and `Ti` dims)
- `low_frequency_data`: DimArray with slow data or `nothing`
- `kwargs...`: Forwarded to step implementations
"""
function process!(pipeline::EddyPipeline, high_frequency_data::DimArray, low_frequency_data::Union{Nothing, DimArray}; kwargs...)
    logger = pipeline.logger

    # Log which pipeline steps are active
    for (step_name, step_field) in [
        (:quality_control, pipeline.quality_control),
        (:gas_analyzer,    pipeline.gas_analyzer),
        (:despiking,       pipeline.despiking),
        (:make_continuous, pipeline.make_continuous),
        (:gap_filling,     pipeline.gap_filling),
        (:double_rotation, pipeline.double_rotation),
        (:mrd,             pipeline.mrd),
    ]
        log_metadata!(logger, "step_$step_name", step_field === nothing ? "disabled" : string(typeof(step_field)))
    end

    # Log calibration coefficients if present
    coeffs = get_calibration_coefficients(pipeline.sensor)
    if coeffs !== nothing
        log_metadata!(logger, "calibration_A", string(coeffs.A))
        log_metadata!(logger, "calibration_B", string(coeffs.B))
        log_metadata!(logger, "calibration_C", string(coeffs.C))
        log_metadata!(logger, "calibration_H2O_Zero", string(coeffs.H2O_Zero))
        log_metadata!(logger, "calibration_H20_Span", string(coeffs.H20_Span))
    end

    log_metadata!(logger, "sensor_type", string(typeof(pipeline.sensor)))

    result = nothing
    total_seconds = @elapsed begin
        prog = ProgressUnknown(;desc="Peddy is cleaning your data...", spinner=true)

        check_data(high_frequency_data, low_frequency_data, pipeline.sensor)

        next!(prog; showvalues=[("Status", "Performing Quality Control...")], spinner="🔬")
        _run_stage(logger, :quality_control) do
            quality_control!(pipeline.quality_control, high_frequency_data, low_frequency_data,
                             pipeline.sensor; logger=logger, kwargs...)
        end

        next!(prog; showvalues=[("Status", "Correcting Gas Analyzer...")], spinner="🧹")
        _run_stage(logger, :gas_analyzer) do
            correct_gas_analyzer!(pipeline.gas_analyzer, high_frequency_data, low_frequency_data,
                                  pipeline.sensor; logger=logger, kwargs...)
        end

        next!(prog; showvalues=[("Status", "Removing Spikes...")], spinner="🦔")
        _run_stage(logger, :despiking) do
            despike!(pipeline.despiking, high_frequency_data, low_frequency_data; logger=logger, kwargs...)
        end

        next!(prog; showvalues=[("Status", "Making Time Continuous...")], spinner="⏱️")
        high_frequency_data = _run_stage(logger, :make_continuous) do
            make_continuous!(pipeline.make_continuous, high_frequency_data,
                             low_frequency_data; logger=logger, kwargs...)
        end

        next!(prog; showvalues=[("Status", "Filling Gaps...")], spinner="🧩")
        _run_stage(logger, :gap_filling) do
            fill_gaps!(pipeline.gap_filling, high_frequency_data, low_frequency_data; logger=logger, kwargs...)
        end

        next!(prog; showvalues=[("Status", "Applying Double Rotation...")], spinner="🌀")
        _run_stage(logger, :double_rotation) do
            rotate!(pipeline.double_rotation, high_frequency_data, low_frequency_data; logger=logger, kwargs...)
        end

        next!(prog; showvalues=[("Status", "Decomposing MRD...")], spinner="〰️")
        _run_stage(logger, :mrd) do
            decompose!(pipeline.mrd, high_frequency_data, low_frequency_data; logger=logger, kwargs...)
        end

        next!(prog; showvalues=[("Status", "Writing Data...")], spinner="💾")
        _run_stage(logger, :output) do
            write_data(pipeline.output, high_frequency_data, low_frequency_data; logger=logger, kwargs...)
        end

        result = finish!(prog; desc="Peddy is done cleaning your data!", spinner="🎉")
    end
    record_stage_time!(logger, :pipeline_total, total_seconds)
    return result
end

@inline function _run_stage(f::F, logger::ProcessingLogger, stage::Symbol) where {F}
    local result
    seconds = @elapsed result = f()
    record_stage_time!(logger, stage, seconds)
    return result
end

@inline _run_stage(f::F, ::NoOpLogger, ::Symbol) where {F} = f()

"""
    check_data(high_frequency_data::DimArray, low_frequency_data::Union{Nothing,DimArray},
               sensor::AbstractSensor) -> Nothing

Validate that the high-frequency data has the required dimensions and variables
for the specified `sensor`. Throws an `ArgumentError` if a requirement is not met.
Currently only the high-frequency data is validated.
"""
function check_data(high_frequency_data::DimArray, low_frequency_data::Union{Nothing,DimArray}, sensor::AbstractSensor)
    # FAQ: Should we check the low frequency data?

    needed_cols = needs_data_cols(sensor)
    if !(:Var in DimensionalData.name.(dims(high_frequency_data)))
        throw(ArgumentError("High frequency data must have a Var dimension"))
    end
    if !(:Ti in DimensionalData.name.(dims(high_frequency_data)))
        throw(ArgumentError("High frequency data must have a Time dimension"))
    end
    var_names = val(dims(high_frequency_data, :Var))
    for col in needed_cols
        if !(col in var_names)
            throw(ArgumentError("Var dimension must have a $col variable"))
        end
    end
    @debug "High frequency data checked, no checks performed on low frequency data"
end
