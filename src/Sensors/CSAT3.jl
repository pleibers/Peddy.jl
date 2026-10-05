"""
     CSAT3(; diag_sonic=4096, diag_margin=0)

 Campbell Scientific CSAT3 sonic anemometer.

 The `diag_sonic` threshold controls which diagnostic values are considered invalid.
 During `check_diagnostics!`, records exceeding this threshold are set to `NaN` for
 wind components and sonic temperature, together with `diag_margin` records on
 either side of each flagged record.
 """
@kwdef struct CSAT3 <: AbstractSensor
    diag_sonic::Int = 4096
    diag_margin::Int = 0
end
needs_data_cols(sensor::CSAT3) = (:diag_sonic, :Ux, :Uy, :Uz, :Ts)
has_variables(sensor::CSAT3) = (:diag_sonic, :Ux, :Uy, :Uz, :Ts)

function check_diagnostics!(sensor::CSAT3, data::DimArray; kwargs...)
    # Now sonic and csat diagnostics both discard all Ux, Uy, Uz, Ts
    diag_sonic_col = view(data, Var(At(:diag_sonic)))
    logger = get(kwargs, :logger, nothing)
    flagged = findall(>(sensor.diag_sonic), diag_sonic_col)
    n_margin = discard_sonic_with_margin!(data, flagged, sensor.diag_margin)
    if logger !== nothing && !isempty(flagged)
        times = collect(dims(data, Ti))
        log_index_runs!(logger, :quality_control, :diagnostic_flag, :diag_sonic, times, flagged;
                        include_run_length=true, threshold=sensor.diag_sonic, affected_variables=[:Ux, :Uy, :Uz, :Ts])
        n_margin > 0 && log_event!(logger, :quality_control, :diagnostic_margin;
                                   variable=:diag_sonic, margin=sensor.diag_margin, extra_samples=n_margin)
    end
end