export H2OCalibration
export correct_gas_analyzer!
export get_calibration_coefficients

using Polynomials
using Statistics



"""
    get_calibration_coefficients(sensor)

Extract H2O calibration coefficients from sensor if available.
Returns nothing if sensor doesn't have calibration coefficients.
"""
function get_calibration_coefficients(sensor::AbstractSensor)
    return hasfield(typeof(sensor), :calibration_coefficients) ?
           sensor.calibration_coefficients : nothing
end

"""
    H2OCalibration

Gas analyzer correction for H2O measurements using calibration coefficients.
Implements bias correction based on relative humidity reference measurements.
Calibration coefficients are extracted from the sensor during correction.

Fields:
- `h2o_variable`: H2O variable name, must be in high frequency data (default: :H2O)
- `pressure_var`: Pressure variable name, must be in high frequency data (default: :P)
- `temp_var`: Temperature variable name in slow data (default: :TA)
- `rh_var`: Relative humidity variable name in slow data (default: :RH)
"""
@kwdef struct H2OCalibration <: AbstractGasAnalyzer
    h2o_variable::Symbol = :H2O
    pressure_var::Symbol = :P
    temp_var::Symbol = :TA
    rh_var::Symbol = :RH
end

"""
    compute_h2o_concentration(RH, TA)

Compute H2O concentration from relative humidity (RH) and temperature (TA).

Parameters:
- RH: Relative humidity in percentage
- TA: Temperature in degrees Celsius

Returns:
- H2O concentration in mmol m⁻³
"""
function compute_h2o_concentration(relative_humidity, air_temperature)
    saturation_vapor_pressure = 611.2 * exp(17.67 * air_temperature / (air_temperature + 243.5)) * (relative_humidity / 100)  # Pa
    h2o_concentration = (1000 * saturation_vapor_pressure / (8.314 * (air_temperature + 273.15)))  # mmol m⁻³
    return h2o_concentration
end

"""
    solve_polynomial_absorptance(y, coeffs)

Solve the cubic polynomial to find absorptance from normalized y value.
Polynomial: A*a + B*a² + C*a³ + y = 0

Parameters:
- y: Normalized y value (negated)
- coeffs: H2OCalibrationCoefficients struct with polynomial coefficients A, B, C

Returns:
- Real root of the polynomial (absorptance)
"""
function solve_polynomial_absorptance(normalized_y_value, calibration_coefficients)
    if isnan(normalized_y_value)
        return NaN
    end

    # Create polynomial: C*a³ + B*a² + A*a + y = 0
    polynomial = Polynomial([normalized_y_value, calibration_coefficients.A, calibration_coefficients.B, calibration_coefficients.C])
    polynomial_roots = roots(polynomial)

    # Find the real root (should be only one for physical solution)
    real_roots = [root for root in polynomial_roots if abs(imag(root)) < 1e-10]

    if isempty(real_roots)
        return NaN
    end

    return real(real_roots[1])
end

"""
    _find_nearest_lf_index(hf_time, lf_times) -> Vector{Int}

For each high-frequency timestamp, find the index of the nearest (by absolute
time difference) low-frequency timestamp. Returns a vector of LF indices, one
per HF point.
"""
function _find_nearest_lf_index(hf_times, lf_times)
    n_hf = length(hf_times)
    n_lf = length(lf_times)
    indices = Vector{Int}(undef, n_hf)

    # Use a running pointer since both are sorted
    lf_ptr = 1
    for i in 1:n_hf
        t = hf_times[i]
        # Advance pointer while the next LF point is closer
        while lf_ptr < n_lf &&
              abs(Dates.value(lf_times[lf_ptr + 1] - t)) <= abs(Dates.value(lf_times[lf_ptr] - t))
            lf_ptr += 1
        end
        indices[i] = lf_ptr
    end

    return indices
end

"""
    resample_to_low_frequency(high_frequency_data, h2o_var, pressure_var, lf_times, hf_times)

Resample high-frequency H2O and pressure data to low-frequency timestamps using
timestamp-based matching. For each LF timestamp, averages HF points whose
nearest LF neighbour is that timestamp.

Returns:
- low_frequency_h2o_averages: Low-frequency averaged H2O data
- low_frequency_pressure_averages: Low-frequency averaged pressure data (in Pa)
"""
function resample_to_low_frequency(high_frequency_data, h2o_variable_name, pressure_variable_name,
                                   lf_times, hf_times, hf_to_lf_indices)
    n_lf = length(lf_times)
    T = eltype(high_frequency_data)
    low_frequency_h2o_averages = fill(T(NaN), n_lf)
    low_frequency_pressure_averages = fill(T(NaN), n_lf)

    h2o_col = collect(high_frequency_data[Var=At(h2o_variable_name)])
    pres_col = collect(high_frequency_data[Var=At(pressure_variable_name)])

    # Accumulate sums per LF bin
    h2o_sum  = zeros(T, n_lf)
    pres_sum = zeros(T, n_lf)
    counts   = zeros(Int, n_lf)

    for i in eachindex(hf_to_lf_indices)
        lf_idx = hf_to_lf_indices[i]
        h = h2o_col[i]
        p = pres_col[i]
        if !isnan(h) && !isnan(p)
            h2o_sum[lf_idx]  += h
            pres_sum[lf_idx] += p
            counts[lf_idx]   += 1
        end
    end

    for j in 1:n_lf
        if counts[j] > 0
            low_frequency_h2o_averages[j]      = h2o_sum[j] / counts[j]
            low_frequency_pressure_averages[j]  = pres_sum[j] / counts[j]
        end
    end

    # Convert pressure to Pa
    low_frequency_pressure_averages .*= 1000

    return low_frequency_h2o_averages, low_frequency_pressure_averages
end

"""
    calculate_reference_absorptances(low_frequency_h2o_averages, low_frequency_pressure_averages, relative_humidity_data, temperature_data,
                                          calibration_coefficients)

Calculate reference absorptances from LI measurements and RH/temperature data.

Returns:
- li_a_raw_lf: LI raw absorptances at low frequency
- rh_a_raw_lf: RH reference raw absorptances at low frequency
"""
function calculate_reference_absorptances(h2o_low_freq_averages, pressure_low_freq_averages, relative_humidity_data, temperature_data,
                                          calibration_coefficients)
    # Compute reference H2O concentration from RH and temperature
    reference_h2o_concentrations = compute_h2o_concentration.(relative_humidity_data, temperature_data)

    # Calculate normalized y values for polynomial
    licor_normalized_y_values = h2o_low_freq_averages ./ pressure_low_freq_averages .* 1000  # mmol m⁻³ kPa⁻¹
    reference_normalized_y_values = reference_h2o_concentrations ./ pressure_low_freq_averages .* 1000  # mmol m⁻³ kPa⁻¹

    # Calculate absorptances using polynomial
    licor_absorptances = [solve_polynomial_absorptance(-y_value, calibration_coefficients) for y_value in licor_normalized_y_values]
    reference_absorptances = [solve_polynomial_absorptance(-y_value, calibration_coefficients) for y_value in reference_normalized_y_values]

    # Calculate raw absorptances
    licor_raw_absorptances_low_freq = licor_absorptances .* pressure_low_freq_averages ./ 1000 ./ calibration_coefficients.H20_Span
    reference_raw_absorptances_low_freq = reference_absorptances .* pressure_low_freq_averages ./ 1000 ./ calibration_coefficients.H20_Span

    return licor_raw_absorptances_low_freq, reference_raw_absorptances_low_freq
end

"""
    interpolate_to_high_frequency(li_a_raw_lf, rh_a_raw_lf, hf_to_lf_indices)

Interpolate low-frequency absorptances to high-frequency using the precomputed
nearest-neighbour index mapping.

Returns:
- li_a_raw_hf: LI raw absorptances at high frequency
- rh_a_raw_hf: RH reference raw absorptances at high frequency
"""
function interpolate_to_high_frequency(licor_raw_absorptances_low_freq, reference_raw_absorptances_low_freq, hf_to_lf_indices)
    n_hf = length(hf_to_lf_indices)
    licor_raw_absorptances_high_freq = Vector{eltype(licor_raw_absorptances_low_freq)}(undef, n_hf)
    reference_raw_absorptances_high_freq = Vector{eltype(licor_raw_absorptances_low_freq)}(undef, n_hf)

    for i in 1:n_hf
        lf_idx = hf_to_lf_indices[i]
        licor_raw_absorptances_high_freq[i] = licor_raw_absorptances_low_freq[lf_idx]
        reference_raw_absorptances_high_freq[i] = reference_raw_absorptances_low_freq[lf_idx]
    end

    return licor_raw_absorptances_high_freq, reference_raw_absorptances_high_freq
end

"""
    apply_bias_correction(h2o_data, pres_data, li_a_raw_hf, rh_a_raw_hf, coeffs)

Apply bias correction to high-frequency H2O measurements.

Returns:
- h2o_corrected: Bias-corrected H2O concentrations
"""
function apply_bias_correction(h2o_measurements, pressure_measurements, licor_raw_absorptances_high_freq, reference_raw_absorptances_high_freq, calibration_coefficients)
    # Calculate high-frequency normalized y and absorptance
    licor_normalized_y_high_freq = h2o_measurements ./ pressure_measurements .* 1000
    licor_absorptances_high_freq = [solve_polynomial_absorptance(-y_value, calibration_coefficients) for y_value in licor_normalized_y_high_freq]
    licor_raw_absorptances_calculated = licor_absorptances_high_freq .* pressure_measurements ./ 1000 ./ calibration_coefficients.H20_Span

    # Apply bias correction formula
    corrected_absorptances_high_freq = ((1 .- reference_raw_absorptances_high_freq) .* licor_raw_absorptances_calculated .- licor_raw_absorptances_high_freq .+ reference_raw_absorptances_high_freq) ./
                                      (1 .- licor_raw_absorptances_high_freq)

    # Convert back to concentration
    normalized_corrected_absorptances = corrected_absorptances_high_freq ./ pressure_measurements .* 1000 .* calibration_coefficients.H20_Span
    normalized_corrected_y_values = calibration_coefficients.A .* normalized_corrected_absorptances .+ calibration_coefficients.B .* normalized_corrected_absorptances .^ 2 .+
                                   calibration_coefficients.C .* normalized_corrected_absorptances .^ 3
    h2o_corrected_concentrations = normalized_corrected_y_values .* pressure_measurements ./ 1000  # mmol/m³

    # Round to 1 decimal place
    return round.(h2o_corrected_concentrations, digits=1)
end

"""
    correct_gas_analyzer!(gas_analyzer::H2OCalibration, high_frequency_data, low_frequency_data, sensor; kwargs...)

Apply H2O calibration correction to gas analyzer measurements.

This function:
1. Extracts calibration coefficients from gas_analyzer or sensor
2. Builds a timestamp-based mapping between high-frequency and low-frequency data
3. Resamples high-frequency data to low-frequency intervals
4. Computes reference H2O from RH and temperature in slow data
5. Calculates absorptances using calibration polynomial
6. Applies bias correction to high-frequency H2O measurements
7. Replaces original H2O variable with corrected data in-place

Parameters:
- gas_analyzer: H2OCalibration configuration
- high_frequency_data: High-frequency measurements (DimArray)
- low_frequency_data: Low-frequency/slow measurements (DimArray)
- sensor: Sensor instance (used to extract calibration coefficients if not provided in gas_analyzer)
"""
function correct_gas_analyzer!(gas_analyzer::H2OCalibration, high_frequency_data::DimArray,
                               low_frequency_data::DimArray, sensor; kwargs...)
    # Get calibration coefficients from sensor
    calibration_coefficients = get_calibration_coefficients(sensor)
    if calibration_coefficients === nothing
        @warn "No calibration coefficients found in sensor $(typeof(sensor)). H2O calibration will be skipped. Only LICOR sensors with calibration coefficients support H2O correction."
        return nothing
    end

    println("Using calibration coefficients from sensor: $(typeof(sensor))")

    # Get time dimensions
    hf_times = collect(dims(high_frequency_data, Ti))
    lf_times = collect(dims(low_frequency_data, Ti))

    # Check temporal overlap
    hf_start, hf_end = extrema(hf_times)
    lf_start, lf_end = extrema(lf_times)
    if lf_end < hf_start || lf_start > hf_end
        @warn "No temporal overlap between high-frequency ($(hf_start)–$(hf_end)) and low-frequency ($(lf_start)–$(lf_end)) data. H2O calibration skipped."
        return nothing
    end

    # Build timestamp-based mapping: for each HF point, find nearest LF point
    hf_to_lf_indices = _find_nearest_lf_index(hf_times, lf_times)

    # Extract required variables
    temperature_data = collect(low_frequency_data[Var(At(gas_analyzer.temp_var))])
    relative_humidity_data = collect(low_frequency_data[Var(At(gas_analyzer.rh_var))]) .* 100  # Convert to percentage
    h2o_variable = gas_analyzer.h2o_variable
    pressure_variable = gas_analyzer.pressure_var

    # Step 1: Resample high-frequency data to low-frequency (timestamp-aligned)
    h2o_low_freq_averages, pressure_low_freq_averages = resample_to_low_frequency(
        high_frequency_data, h2o_variable, pressure_variable,
        lf_times, hf_times, hf_to_lf_indices)

    # Step 2: Calculate reference absorptances
    licor_raw_absorptances_low_freq, reference_raw_absorptances_low_freq = calculate_reference_absorptances(
        h2o_low_freq_averages, pressure_low_freq_averages,
        relative_humidity_data, temperature_data, calibration_coefficients)

    # Step 3: Interpolate absorptances to high frequency (timestamp-aligned)
    licor_raw_absorptances_high_freq, reference_raw_absorptances_high_freq = interpolate_to_high_frequency(
        licor_raw_absorptances_low_freq, reference_raw_absorptances_low_freq, hf_to_lf_indices)

    # Step 4: Apply bias correction to high-frequency data
    h2o_measurements = collect(high_frequency_data[Var(At(gas_analyzer.h2o_variable))])
    pressure_measurements = collect(high_frequency_data[Var(At(gas_analyzer.pressure_var))]) .* 1000  # Convert to Pa

    h2o_corrected_concentrations = apply_bias_correction(h2o_measurements, pressure_measurements, licor_raw_absorptances_high_freq, reference_raw_absorptances_high_freq,
                                                        calibration_coefficients)

    # Step 5: Replace the original H2O variable with corrected data in-place
    high_frequency_data[Var(At(gas_analyzer.h2o_variable))] = h2o_corrected_concentrations

    println("H2O calibration correction applied. Variable $h2o_variable has been corrected in-place.")

    return nothing
end
