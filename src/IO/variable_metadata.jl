"""
    VariableMetadata(; standard_name, unit="", long_name="", description="")

Lightweight container for per-variable metadata used by output backends.

Fields:
- `standard_name::String`: CF-style standard name or canonical name.
- `unit::String`: Physical unit (free text, e.g. "m s^-1").
- `long_name::String`: Human-readable name.
- `description`: Free-form description.
"""
@kwdef struct VariableMetadata
    standard_name::String
    unit::String = ""
    long_name::String = ""
    description = ""
end

const DEFAULT_VARIABLE_METADATA = Dict{Symbol,VariableMetadata}(
    :timestamp => VariableMetadata(
        standard_name = "timestamp",
        unit = "",
        long_name = "Timestamp",
        description = "Timestamp in the ISO format"
    ),
    :Ux => VariableMetadata(
        standard_name = "eastward_wind",
        unit = "m s^-1",
        long_name = "u wind component",
        description = "East-west wind component"
    ),
    :Uy => VariableMetadata(
        standard_name = "northward_wind",
        unit = "m s^-1",
        long_name = "v wind component",
        description = "North-south wind component"
    ),
    :Uz => VariableMetadata(
        standard_name = "upward_air_velocity",
        unit = "m s^-1",
        long_name = "w wind component",
        description = "Vertical wind component"
    ),
    :Ts => VariableMetadata(
        standard_name = "sonic_temperature",
        unit = "oC",
        long_name = "sonic temperature",
        description = "Sonic temperature"
    ),
    :CO2 => VariableMetadata(
        standard_name = "carbon_dioxide",
        unit = "mmol/m3",
        long_name = "carbon dioxide",
        description = "Carbon dioxide"
    ),
    :H2O => VariableMetadata(
        standard_name = "water_vapor",
        unit = "mmol/m3",
        long_name = "water vapor",
        description = "Water vapor"
    ),
    :P => VariableMetadata(
        standard_name = "air_pressure",
        unit = "kPa",
        long_name = "air pressure",
        description = "Air pressure"
    ),
    :RH => VariableMetadata(
        standard_name = "relative_humidity",
        unit = "%",
        long_name = "relative humidity",
        description = "Relative humidity"
    )
)

"""
    metadata_for(name::Union{Symbol,AbstractString}) -> VariableMetadata

Return metadata for a variable. If the variable is not present in the
`DEFAULT_VARIABLE_METADATA` registry, a generic `VariableMetadata` is returned
using the variable name for both `standard_name` and `long_name` and empty
unit/description.
"""
function metadata_for(name::Union{Symbol,AbstractString})
    sym = name isa Symbol ? name : Symbol(name)
    if haskey(DEFAULT_VARIABLE_METADATA, sym)
        return DEFAULT_VARIABLE_METADATA[sym]
    else
        s = String(name)
        return VariableMetadata(
            standard_name = s,
            unit = "",
            long_name = s,
            description = "",
        )
    end
end

"""
    get_default_metadata() -> Dict{Symbol,VariableMetadata}

Return the default metadata dictionary used by output backends.
"""
get_default_metadata() = DEFAULT_VARIABLE_METADATA