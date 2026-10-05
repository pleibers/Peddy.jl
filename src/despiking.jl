export SimpleSigmundDespiking, VariableGroup

using Statistics

"""
    VariableGroup(name::String, variables::Vector{Symbol}; spike_threshold::Real=6.0)

Defines a group of variables that are combined for spike detection.

# Parameters
- `name`: Descriptive name for the variable group
- `variables`: Vector of variable symbols to include in this group
- `spike_threshold`: Spike detection threshold for this group (default: 6.0, normalized by 0.6745)

# Examples
```julia
# Wind components group
wind_group = VariableGroup("Wind Components", [:Ux, :Uy, :Uz], spike_threshold=6.0)

# Gas analyzer group with different threshold
gas_group = VariableGroup("Gas Analyzer", [:CO2, :H2O], spike_threshold=6.0)
```
"""
struct VariableGroup{N<:Real}
    name::String
    variables::Vector{Symbol}
    spike_threshold::N
    
    function VariableGroup(name::String, variables::Vector{Symbol}; spike_threshold=6.0)
        new{typeof(spike_threshold)}(name, variables, spike_threshold)
    end
end

"""
    SimpleSigmundDespiking(; window_minutes=5.0, variable_groups=VariableGroup[])

Implements the modified MAD (Median Absolute Deviation) filter for spike detection
based on Sigmund et al. (2022). This despiking method:

1. Calculates rolling median and MAD over a specified window
2. Computes deviation patterns using neighboring points
3. Identifies spikes based on variable group thresholds (each group can have its own threshold)
4. Sets detected spikes to NaN

# Parameters
- `window_minutes`: Window size in minutes for rolling statistics (default: 5.0)
- `variable_groups`: Vector of VariableGroup structs, each with its own threshold
- `use_mad_floor`: When true, clamp MAD values to the per-variable floors in `mad_floor`
- `mad_floor`: Dict of variable => minimum MAD (default 0.02 for Ux, Uy, Uz, Ts)
- `max_ts_deviation`: Optional absolute limit (K). Before the MAD filter, every sample whose
  Ts deviates from the rolling median (same window) by more than this is removed for all
  sonic variables (Ux, Uy, Uz, Ts), since Ts and wind come from the same transit times.
  The limit is fixed, so a long sonic fault cannot raise its own threshold as the MAD does.
  Default `nothing` (off).
- `couple_sonic`: When true, a sample where despiking removes any of Ux, Uy, Uz or Ts loses all
  four. The sonic derives wind and temperature from the same transit times, so a spike in one
  means the others are unreliable too. Default `false`.

# Examples
```julia
# Default: wind components and temperature combined
SimpleSigmundDespiking()

# Custom groups with different thresholds
SimpleSigmundDespiking(
    variable_groups=[
        VariableGroup("Wind Components", [:Ux, :Uy, :Uz], spike_threshold=6.0),
        VariableGroup("Sonic Temperature", [:Ts], spike_threshold=6.0),
        VariableGroup("Gas Analyzer", [:CO2, :H2O], spike_threshold=6.0)
    ]
)

# Single variables with individual thresholds
SimpleSigmundDespiking(
    variable_groups=[
        VariableGroup("Wind U", [:Ux], spike_threshold=6.0),
        VariableGroup("Wind V", [:Uy], spike_threshold=6.0),
        VariableGroup("Wind W", [:Uz], spike_threshold=7.0)
    ]
)
```
"""
struct SimpleSigmundDespiking{N<:Real} <: AbstractDespiking
    window_minutes::N
    variable_groups::Vector{<:VariableGroup}
    use_mad_floor::Bool
    mad_floor::Dict{Symbol,Float64}
    max_ts_deviation::Union{Nothing,Float64}
    couple_sonic::Bool
    
    function SimpleSigmundDespiking(; window_minutes=5.0,
                                    number_type::Type{N}=Float64,
                                    variable_groups=VariableGroup[],
                                    use_mad_floor::Bool=false,
                                    mad_floor::Dict{Symbol,Float64}=Dict(
                                        :Ux => 0.02,
                                        :Uy => 0.02,
                                        :Uz => 0.02,
                                        :Ts => 0.02,
                                    ),
                                    max_ts_deviation::Union{Nothing,Real}=nothing,
                                    couple_sonic::Bool=false) where {N <: Real}
        # Default group if none provided
        if isempty(variable_groups)
            default_group = VariableGroup("Default Sonic", [:Ux, :Uy, :Uz, :Ts], spike_threshold=N(6.0))
            variable_groups = [default_group]
        end
        new{N}(N(window_minutes), variable_groups, use_mad_floor, mad_floor,
               max_ts_deviation === nothing ? nothing : Float64(max_ts_deviation), couple_sonic)
    end
end

Base.eltype(::SimpleSigmundDespiking{N}) where N = N

"""
    despike!(despiking::SimpleSigmundDespiking, high_frequency_data::DimArray, low_frequency_data; kwargs...)

Apply the Sigmund et al. (2022) modified MAD filter for spike detection.

The algorithm:
1. Calculates rolling median over specified window
2. Computes absolute deviations from rolling median
3. Calculates MAD (Median Absolute Deviation) over the same window
4. Applies pattern detection using neighboring points
5. Identifies spikes when combined normalized deviations exceed threshold
6. Sets detected spikes to NaN in-place
"""
function despike!(despiking::SimpleSigmundDespiking, high_frequency_data::DimArray, low_frequency_data; kwargs...)
    # Calculate window size from time dimension
    window_size_points = _calculate_window_size(high_frequency_data, despiking.window_minutes)
    logger = get(kwargs, :logger, nothing)
    timestamps = logger === nothing ? nothing : collect(dims(high_frequency_data, Ti))
    
    println("Processing spike detection with window size: $window_size_points points")

    sonic_vars = [v for v in (:Ux, :Uy, :Uz, :Ts) if v in dims(high_frequency_data, Var)]
    valid_before = despiking.couple_sonic ?
        [.!isnan.(collect(high_frequency_data[Var=At(v)])) for v in sonic_vars] : nothing

    if despiking.max_ts_deviation !== nothing
        _remove_ts_median_outliers!(high_frequency_data, window_size_points,
                                    despiking.max_ts_deviation, logger, timestamps)
    end
    println("Processing $(length(despiking.variable_groups)) variable groups")
    
    # Process each variable group independently
    for (group_idx, var_group) in enumerate(despiking.variable_groups)
        println("  Processing group $group_idx: '$(var_group.name)' with threshold $(var_group.spike_threshold)")
        println("    Variables: $(var_group.variables)")
        
        # Pre-allocate combined deviation array for this group
        num_time_points = length(dims(high_frequency_data, Ti))
        group_combined_deviations = zeros(eltype(high_frequency_data), num_time_points)
        num_variables_in_group = 0
        
        # Process each variable in this group
        for variable_name in var_group.variables
            if variable_name ∉ dims(high_frequency_data, Var)
                continue  # Skip if variable doesn't exist in data
            end
            
            # Extract time series for this variable
            variable_data = @view high_frequency_data[Var=At(variable_name)]
            
            if length(variable_data) < window_size_points
                @warn "Data length ($(length(variable_data))) smaller than window size ($window_size_points) for variable $variable_name"
                continue
            end
            
            # Calculate normalized deviations using MAD filter
            mad_floor_value = (despiking.use_mad_floor ? get(despiking.mad_floor, variable_name, 0.0) : 0.0)
            normalized_deviations = _calculate_mad_normalized_deviations(variable_data, window_size_points; mad_floor=mad_floor_value)
            
            # Accumulate deviations within this group
            num_variables_in_group += 1
            group_combined_deviations .+= normalized_deviations
        end
        
        # Apply spike detection for this group if any variables were processed
        if num_variables_in_group > 0
            _apply_spike_threshold_and_remove_for_group!(var_group, high_frequency_data, 
                                                        group_combined_deviations, group_idx, logger, timestamps)
        else
            @warn "No variables were processed in group $group_idx: '$(var_group.name)'"
        end
    end

    if despiking.couple_sonic && !isempty(sonic_vars)
        _couple_sonic_removals!(high_frequency_data, sonic_vars, valid_before, logger)
    end
end

"""
    _couple_sonic_removals!(high_frequency_data, sonic_vars, valid_before, logger) -> Int

Set all `sonic_vars` to NaN at samples where despiking removed at least one of them
(valid before despiking, NaN after). Returns the number of samples affected.
"""
function _couple_sonic_removals!(high_frequency_data::DimArray, sonic_vars, valid_before, logger)
    n = length(dims(high_frequency_data, Ti))
    removed = falses(n)
    for (v, was_valid) in zip(sonic_vars, valid_before)
        removed .|= was_valid .& isnan.(collect(high_frequency_data[Var=At(v)]))
    end
    n_removed = count(removed)
    n_removed == 0 && return 0
    nan = convert(eltype(high_frequency_data), NaN)
    n_extra = 0
    for v in sonic_vars
        col = @view high_frequency_data[Var=At(v)]
        n_extra += count(removed .& .!isnan.(col))
        col[removed] .= nan
    end
    println("    Sonic coupling: $n_removed samples with a removed sonic variable, $n_extra extra values set to NaN")
    logger === nothing || log_event!(logger, :despiking, :sonic_coupling;
                                     samples=n_removed, extra_values=n_extra)
    return n_removed
end

"""
    calculate_rolling_median(data::AbstractVector, window_size::Int) -> Vector

Calculate rolling median with center alignment, handling NaN values.
Optimized version using partialsort! for O(n) median computation instead of full sort.
"""
function calculate_rolling_median(data::AbstractVector, window_size::Int)
    n = length(data)
    result = similar(data)
    half_window = window_size ÷ 2
    
    # Pre-allocate buffer for valid data to avoid repeated allocations
    max_possible_window = min(2 * half_window + 1, n)
    valid_buffer = Vector{eltype(data)}(undef, max_possible_window)
    
    @inbounds for i in 1:n
        # Define window bounds with center alignment
        start_idx = max(1, i - half_window)
        end_idx = min(n, i + half_window)
        
        # Extract valid (non-NaN) data directly into pre-allocated buffer
        valid_count = 0
        for j in start_idx:end_idx
            val = data[j]
            if !isnan(val)
                valid_count += 1
                valid_buffer[valid_count] = val
            end
        end
        
        if valid_count > 0
            # Use partialsort! for O(n) median instead of O(n log n) full sort
            result[i] = _fast_median!(valid_buffer, valid_count)
        else
            result[i] = NaN
        end
    end
    
    return result
end

"""
    _fast_median!(buffer::Vector, count::Int) -> eltype(buffer)

Compute median of first `count` elements in buffer using partialsort!.
This is O(n) on average vs O(n log n) for full sort.
Modifies buffer in-place.
"""
@inline function _fast_median!(buffer::Vector{T}, count::Int) where T
    if count == 1
        return @inbounds buffer[1]
    elseif count == 2
        return @inbounds (buffer[1] + buffer[2]) / 2
    end
    
    mid = (count + 1) ÷ 2
    if isodd(count)
        # For odd count, partialsort! to find the middle element
        return @inbounds partialsort!(@view(buffer[1:count]), mid)
    else
        # For even count, need the two middle elements
        # partialsort! with a range is efficient for this
        partialsort!(@view(buffer[1:count]), mid:mid+1)
        return @inbounds (buffer[mid] + buffer[mid+1]) / 2
    end
end

"""
    calculate_pattern_deviation(absolute_deviations::AbstractVector) -> Vector

Calculate pattern-adjusted deviations using neighboring points:
pattern_deviation = |deviation| - 0.5 * (|neighbor_left| + |neighbor_right|)
Optimized version with vectorized operations and efficient boundary handling.
"""
function calculate_pattern_deviation(absolute_deviations::AbstractVector)
    n = length(absolute_deviations)
    if n < 3
        # For very small arrays, use simple approach
        return abs.(absolute_deviations) .- 0.5 .* abs.(absolute_deviations)
    end
    
    pattern_adjusted_deviations = similar(absolute_deviations)
    
    # Handle boundaries and interior with @inbounds for performance
    @inbounds begin
        # First element: only right neighbor
        pattern_adjusted_deviations[1] = abs(absolute_deviations[1]) - 0.5 * abs(absolute_deviations[2])
        
        # Last element: only left neighbor  
        pattern_adjusted_deviations[n] = abs(absolute_deviations[n]) - 0.5 * abs(absolute_deviations[n-1])
        
        # Interior points: compute abs inline to avoid allocation
        @simd for i in 2:(n-1)
            abs_i = abs(absolute_deviations[i])
            neighbor_average = 0.5 * (abs(absolute_deviations[i-1]) + abs(absolute_deviations[i+1]))
            pattern_adjusted_deviations[i] = abs_i - neighbor_average
        end
    end
    
    return pattern_adjusted_deviations
end

# Helper functions for improved readability

"""
    _calculate_window_size(high_frequency_data::DimArray, window_minutes::Real) -> Int

Calculate window size in data points from time dimension and desired window duration.
"""
function _calculate_window_size(high_frequency_data::DimArray, window_minutes)
    time_dimension = dims(high_frequency_data, Ti)
    
    if length(time_dimension) < 2
        @warn "Insufficient data points for despiking"
        return 3  # Minimum window size
    end
    
    # Calculate sampling frequency (assuming regular intervals)
    time_step = time_dimension[2] - time_dimension[1]
    time_step_seconds = Dates.value(Millisecond(time_step)) / 1000.0  # Convert Millisecond to seconds (be sure to use Millisecond)
    sampling_frequency_hz = 1.0 / time_step_seconds
    
    # Calculate window size in data points
    window_size_points = round(Int, window_minutes * 60 * sampling_frequency_hz)
    window_size_points = max(window_size_points, 3)  # Ensure minimum window size
    
    # Adjust window size if data is too small (use at most 1/3 of data length)
    total_data_points = length(time_dimension)
    if window_size_points > total_data_points ÷ 3
        window_size_points = max(total_data_points ÷ 3, 3)
        @warn "Adjusted window size to $window_size_points points due to limited data ($total_data_points points)"
    end
    
    return window_size_points
end

"""
    _calculate_mad_normalized_deviations(variable_data::AbstractVector, window_size_points::Int) -> Vector

Calculate MAD-normalized deviations for a single variable using the Sigmund et al. algorithm.
Optimized version that minimizes allocations and reuses intermediate arrays.
"""
function _calculate_mad_normalized_deviations(variable_data::AbstractVector, window_size_points::Int; mad_floor::Real=0.0)
    n = length(variable_data)
    
    # Calculate rolling median
    rolling_median_values = calculate_rolling_median(variable_data, window_size_points)
    
    # Calculate absolute deviations from rolling median (reuse this array)
    absolute_deviations = similar(variable_data)
    @inbounds @simd for i in 1:n
        absolute_deviations[i] = abs(variable_data[i] - rolling_median_values[i])
    end
    
    # Calculate MAD (Median Absolute Deviation)
    mad_values = calculate_rolling_median(absolute_deviations, window_size_points)
    
    # Apply pattern detection using neighboring points (modifies absolute_deviations in-place)
    pattern_adjusted_deviations = calculate_pattern_deviation(absolute_deviations)
    
    # Normalize by MAD with small epsilon to avoid division by zero (reuse pattern_adjusted_deviations)
    floor_val = mad_floor > 0 ? mad_floor : 0.0
    @inbounds @simd for i in 1:n
        mad_val = mad_values[i]
        if floor_val > 0
            mad_val = max(mad_val, floor_val)
        end
        if mad_val == 0
            mad_val = 1e-10
        end
        pattern_adjusted_deviations[i] = pattern_adjusted_deviations[i] / mad_val
    end
    
    return pattern_adjusted_deviations
end

"""
    _apply_spike_threshold_and_remove_for_group!(variable_group::VariableGroup, high_frequency_data::DimArray, 
                                               combined_normalized_deviations::Vector, group_index::Int)

Apply spike detection threshold for a specific variable group and remove detected spikes by setting them to NaN.
Uses the group's own spike threshold. Optimized version with pre-computed threshold and efficient indexing.
"""
function _apply_spike_threshold_and_remove_for_group!(variable_group::VariableGroup, high_frequency_data::DimArray, 
                                                     combined_normalized_deviations::Vector, group_index::Int,
                                                     logger, timestamps)
    # Pre-compute normalized threshold (constant for this group)
    normalized_threshold = variable_group.spike_threshold / 0.6745
    
    # Find spike indices efficiently
    n_time_points = length(combined_normalized_deviations)
    spikes_detected_count = 0
    spike_indices = logger === nothing ? nothing : Int[]
    
    # Pre-cache variable views to avoid repeated dimension lookups in the hot loop
    var_dims = dims(high_frequency_data, Var)
    variable_views = [@view high_frequency_data[Var=At(vname)] for vname in variable_group.variables if vname ∈ var_dims]
    
    # Count spikes and apply threshold in single pass
    @inbounds for time_index in 1:n_time_points
        if abs(combined_normalized_deviations[time_index]) >= normalized_threshold
            spikes_detected_count += 1
            logger === nothing || push!(spike_indices, time_index)
            
            # Set spikes to NaN for all variables in this group immediately
            for var_view in variable_views
                var_view[time_index] = NaN
            end
        end
    end
    
    # Report results
    if spikes_detected_count > 0
        println("    Group $group_index ('$(variable_group.name)'): Detected $spikes_detected_count spikes, setting to NaN")
        if logger !== nothing && spike_indices !== nothing && !isempty(spike_indices)
            log_index_runs!(logger, :despiking, :spike, nothing, timestamps, spike_indices;
                            include_run_length=true, group=variable_group.name,
                            threshold=variable_group.spike_threshold,
                            variables=variable_group.variables)
        end
    else
        println("    Group $group_index ('$(variable_group.name)'): No spikes detected")
    end
end

"""
    _remove_ts_median_outliers!(high_frequency_data, window_size_points, limit, logger, timestamps) -> Int

Set Ux, Uy, Uz and Ts to NaN wherever `|Ts - rolling median(Ts)| > limit` (K), using a
centered rolling median over `window_size_points`. Returns the number of removed samples.
"""
function _remove_ts_median_outliers!(high_frequency_data::DimArray, window_size_points::Int,
                                     limit::Real, logger, timestamps)
    :Ts in dims(high_frequency_data, Var) || return 0
    ts = @view high_frequency_data[Var=At(:Ts)]
    rolling_median_values = calculate_rolling_median(ts, window_size_points)
    outliers = findall(i -> !isnan(ts[i]) && !isnan(rolling_median_values[i]) &&
                            abs(ts[i] - rolling_median_values[i]) > limit, eachindex(ts))
    isempty(outliers) && return 0

    nan = convert(eltype(high_frequency_data), NaN)
    for v in (:Ux, :Uy, :Uz, :Ts)
        v in dims(high_frequency_data, Var) || continue
        col = @view high_frequency_data[Var=At(v)]
        col[outliers] .= nan
    end
    println("    Ts median check: removed $(length(outliers)) samples with |Ts - median| > $limit K")
    if logger !== nothing && timestamps !== nothing
        log_index_runs!(logger, :despiking, :ts_median_deviation, :Ts, timestamps, outliers;
                        include_run_length=true, limit=limit, affected_variables=[:Ux, :Uy, :Uz, :Ts])
    end
    return length(outliers)
end
