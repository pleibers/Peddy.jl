# If it is very confusing we can get rid of has_variables in favor of not using diag stuff
export CSAT3, CSAT3B, LICOR, IRGASON
export default_calibration_coefficients

"""
    discard_sonic_with_margin!(data, flagged, margin) -> Int

Set Ux, Uy, Uz and Ts to NaN at the `flagged` indices and at `margin` samples on
either side of each. During a sonic fault the anemometer alternates between flagged
and seemingly clean samples, and those neighbours are often wrong too.
Returns the number of extra samples removed by the margin.
"""
function discard_sonic_with_margin!(data::DimArray, flagged::AbstractVector{Int}, margin::Int)
    isempty(flagged) && return 0
    n = length(dims(data, Ti))
    bad = falses(n)
    for i in flagged
        bad[max(1, i - margin):min(n, i + margin)] .= true
    end
    nan = convert(eltype(data), NaN)
    for v in (:Ux, :Uy, :Uz, :Ts)
        col = view(data, Var(At(v)))
        col[bad] .= nan
    end
    return count(bad) - length(unique(flagged))
end

include("CSAT3.jl")
include("CSAT3B.jl")
include("LICOR.jl")
include("IRGASON.jl")
