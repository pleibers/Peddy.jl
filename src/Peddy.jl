module Peddy

using Dates
using Reexport
@reexport using DimensionalData
using DimensionalData: @dim

abstract type PipelineStep end

"""
    AbstractSensor

Abstract supertype for sensors. Concrete sensor types (e.g. `CSAT3`, `IRGASON`, `LICOR`)
define how diagnostics are interpreted and which variables are required/present.
"""
abstract type AbstractSensor end

"""
    check_diagnostics!(sensor::AbstractSensor, data::DimArray; kwargs...)

Sensor-specific diagnostics check. Implementations may set invalid records to `NaN`
and optionally log affected indices via the processing logger.
"""
function check_diagnostics! end

"""
    AbstractQC

Abstract supertype for quality control (QC) pipeline steps.
"""
abstract type AbstractQC <: PipelineStep end

"""
    quality_control!(qc::AbstractQC, high_frequency_data, low_frequency_data, sensor; kwargs...)

Apply quality control to the data in-place.
"""
function quality_control! end

"""
    AbstractDespiking

Abstract supertype for despiking pipeline steps.
"""
abstract type AbstractDespiking <: PipelineStep end

"""
    despike!(desp::AbstractDespiking, high_frequency_data, low_frequency_data; kwargs...)

Apply despiking to the data in-place.
"""
function despike! end

"""
    AbstractGapFilling

Abstract supertype for gap-filling pipeline steps.
"""
abstract type AbstractGapFilling <: PipelineStep end

"""
    fill_gaps!(gap::AbstractGapFilling, high_frequency_data, low_frequency_data; kwargs...)

Fill small gaps (represented as `NaN`) in-place.
"""
function fill_gaps! end

"""
    AbstractMakeContinuous

Abstract supertype for steps that enforce a continuous time axis.
"""
abstract type AbstractMakeContinuous <: PipelineStep end

"""
    make_continuous!(step::AbstractMakeContinuous, high_frequency_data, low_frequency_data; kwargs...)

Insert missing timestamps and fill inserted rows with `NaN`.
"""
function make_continuous! end

"""
    AbstractGasAnalyzer

Abstract supertype for gas analyzer correction steps.
"""
abstract type AbstractGasAnalyzer <: PipelineStep end

"""
    correct_gas_analyzer!(step::AbstractGasAnalyzer, high_frequency_data, low_frequency_data, sensor; kwargs...)

Apply gas analyzer corrections (e.g. H2O correction).
"""
function correct_gas_analyzer! end

"""
    AbstractDoubleRotation

Abstract supertype for coordinate rotation steps.
"""
abstract type AbstractDoubleRotation <: PipelineStep end

"""
    rotate!(step::AbstractDoubleRotation, high_frequency_data, low_frequency_data; kwargs...)

Apply coordinate rotation to wind components in-place.
"""
function rotate! end

"""
    AbstractMRD

Abstract supertype for multi-resolution decomposition (MRD) steps.
"""
abstract type AbstractMRD <: PipelineStep end

"""
    decompose!(step::AbstractMRD, high_frequency_data, low_frequency_data; kwargs...)

Run a multi-resolution decomposition. Implementations store results inside `step`.
"""
function decompose! end

"""
    AbstractOutput

Abstract supertype for output writers.
"""
abstract type AbstractOutput <: PipelineStep end

"""
    write_data(output::AbstractOutput, high_frequency_data, low_frequency_data=nothing; kwargs...)

Write processed data to an output sink (files, memory, etc.).
"""
function write_data end

const OptionalPipelineStep = Union{Nothing,PipelineStep}

@dim Var "Variables"
export Var

# Compute mean while ignoring NaN values. Returns NaN if all values are NaN.
function mean_skipnan(arr)
    s = zero(eltype(arr))
    c = 0
    # @inbounds
    for v in arr
        if !isnan(v)
            s += v
            c += 1
        end
    end
    return c == 0 ? NaN : s / c
end

include("Sensors/sensors.jl")

include("logging.jl")

include("pipeline.jl")

include("IO/IO.jl")

include("QC/QC.jl")
include("h2o_correction.jl")
include("despiking.jl")
include("make_continuous.jl")
include("interpolation.jl")
include("double_rotation.jl")
include("MRD/mrd.jl")
include("MRD/mrd_plotting.jl")

export AbstractSensor, AbstractQC, AbstractDespiking, AbstractGapFilling, AbstractMakeContinuous,
       AbstractGasAnalyzer, AbstractDoubleRotation, AbstractMRD, AbstractOutput
export write_data
export check_diagnostics!, quality_control!, despike!, fill_gaps!, correct_gas_analyzer!,
    rotate!, decompose!, make_continuous!
export AbstractProcessingLogger, ProcessingLogger, NoOpLogger,
       log_event!, record_stage_time!, write_processing_log,
       log_index_runs!, log_mask_runs!, is_logging_enabled,
       log_metadata!, reset!

end
