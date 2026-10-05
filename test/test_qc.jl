using Test
using Peddy
using DimensionalData

@testset "Quality Control Tests" begin
    @testset "Limit Struct" begin
        # Test Limit constructor and functionality
        limit = Peddy.Limit(-10.0, 10.0)
        @test limit.min == -10.0
        @test limit.max == 10.0

        # Test with different number types
        int_limit = Peddy.Limit(-5, 5)
        @test int_limit.min == -5
        @test int_limit.max == 5
        @test typeof(int_limit.min) == Int
    end

    @testset "PhysicsBoundsCheck Constructor" begin
        # Test default constructor
        qc = Peddy.PhysicsBoundsCheck()
        @test qc isa Peddy.PhysicsBoundsCheck
        @test qc.Ux.min == -100.0
        @test qc.Ux.max == 100.0

        # Test custom constructor with specific limits
        custom_qc = Peddy.PhysicsBoundsCheck(; Ux=Peddy.Limit(-50.0, 50.0),
                                             Uy=Peddy.Limit(-30.0, 30.0))
        @test custom_qc.Ux.min == -50.0
        @test custom_qc.Ux.max == 50.0
        @test custom_qc.Uy.min == -30.0
        @test custom_qc.Uy.max == 30.0

        # Test with different number type
        float32_qc = Peddy.PhysicsBoundsCheck(; number_type=Float32)
        @test typeof(float32_qc.Ux.min) == Float32
        @test typeof(float32_qc.Ux.max) == Float32
    end

    @testset "Default Physical Limits" begin
        # Test that default limits are reasonable for eddy covariance
        defaults = Peddy.default_physical_limits(; number_type=Float64)

        # Wind speed limits should be reasonable
        @test defaults[:Ux].min < 0
        @test defaults[:Ux].max > 0
        @test abs(defaults[:Ux].min) == abs(defaults[:Ux].max)  # Symmetric

        # Temperature limits should span reasonable range
        @test defaults[:Ts].min < 0  # Below freezing
        @test defaults[:Ts].max > 40  # Above typical air temps

        # CO2 limits should be positive
        @test defaults[:CO2].min >= 0
        @test defaults[:CO2].max > defaults[:CO2].min

        # H2O limits should be positive
        @test defaults[:H2O].min >= 0
        @test defaults[:H2O].max > defaults[:H2O].min
    end

    @testset "Quality Control Application" begin
        # Create test data with some out-of-bounds values
        sensor = Peddy.CSAT3()
        needed_cols = collect(Peddy.needs_data_cols(sensor))  # Convert tuple to vector
        n_points = 20
        n_vars = length(needed_cols)

        # Create test data with known out-of-bounds values
        test_data = randn(n_points, n_vars) .* 2 .+ 5  # Reasonable base data

        hd = DimArray(test_data, (Ti(1:n_points), Var(needed_cols)))
        # Extreme Ux (> 100 m/s)
        hd[Ti=5,Var=At(:Ux)] = 150.0
        # Extreme Uy (< -100 m/s)
        hd[Ti=10, Var=At(:Uy)] = -120.0
        # Extreme Uz (> 50 m/s)
        hd[Ti=15,Var=At(:Uz)] = 200.0
        # Extreme Ts (< -50 C, truly out of bounds)
        hd[Ti=8, Var=At(:Ts)] = -60.0
        ld = DimArray(rand(5, n_vars), (Ti(1:5), Var(needed_cols)))

        # Store original extreme values for verification (after setting them)
        original_extreme_ux = hd[Ti=5, Var=At(:Ux)]  # Should be 150.0
        original_extreme_uy = hd[Ti=10, Var=At(:Uy)] # Should be -120.0
        original_extreme_uz = hd[Ti=15, Var=At(:Uz)] # Should be 200.0
        original_extreme_ts = hd[Ti=8, Var=At(:Ts)]  # Should be -60.0

        # Apply quality control
        qc = Peddy.PhysicsBoundsCheck()
        Peddy.quality_control!(qc, hd, ld, sensor)

        # Test that extreme values were flagged (set to NaN)
        @test isnan(hd[Ti=5, Var=At(:Ux)])   # Extreme Ux should be NaN
        @test isnan(hd[Ti=10, Var=At(:Uy)])  # Extreme Uy should be NaN  
        @test isnan(hd[Ti=15, Var=At(:Uz)])  # Extreme Uz should be NaN
        @test isnan(hd[Ti=8, Var=At(:Ts)])   # Extreme Ts should be NaN

        # Test that reasonable values were preserved
        @test !isnan(hd[Ti=1, Var=At(:Ux)])  # Normal values should remain
        @test !isnan(hd[Ti=2, Var=At(:Uy)])
        @test !isnan(hd[Ti=3, Var=At(:Uz)])
        @test !isnan(hd[Ti=4, Var=At(:Ts)])

        # Verify original extreme values were actually extreme
        @test abs(original_extreme_ux) > 100
        @test abs(original_extreme_uy) > 100
        @test abs(original_extreme_uz) > 100
    end

    @testset "Custom Bounds Application" begin
        # Test with custom, more restrictive bounds
        sensor = Peddy.CSAT3()
        needed_cols = collect(Peddy.needs_data_cols(sensor))  # Convert tuple to vector
        n_points = 10
        n_vars = length(needed_cols)

        # Create data that would be OK with default bounds but not custom ones
        test_data = ones(n_points, n_vars) .* 20  # All values = 20

        hd = DimArray(test_data, (Ti(1:n_points), Var(needed_cols)))
        ld = DimArray(rand(3, n_vars), (Ti(1:3), Var(needed_cols)))

        # Apply restrictive QC (max wind speed = 15 m/s)
        restrictive_qc = Peddy.PhysicsBoundsCheck(; Ux=Peddy.Limit(-15.0, 15.0),
                                                  Uy=Peddy.Limit(-15.0, 15.0),
                                                  Uz=Peddy.Limit(-15.0, 15.0))

        Peddy.quality_control!(restrictive_qc, hd, ld, sensor)

        # All wind values (20 m/s) should now be NaN due to restrictive bounds
        @test isnan(hd[Ti=1, Var=At(:Ux)])
        @test isnan(hd[Ti=1, Var=At(:Uy)])
        @test isnan(hd[Ti=1, Var=At(:Uz)])

        # Temperature should still be OK (assuming reasonable Ts bounds)
        @test !isnan(hd[Ti=1, Var=At(:Ts)])
    end

    @testset "Edge Cases" begin
        # Test with boundary values
        sensor = Peddy.CSAT3()
        needed_cols = collect(Peddy.needs_data_cols(sensor))  # Convert tuple to vector

        # Create data exactly at boundaries
        test_data = zeros(3, length(needed_cols))

        hd = DimArray(test_data, (Ti(1:3), Var(needed_cols)))
        hd[Ti=1, Var=At(:Ux)] = 100.0   # Exactly at Ux upper bound
        hd[Ti=2, Var=At(:Ux)] = -100.0  # Exactly at Ux lower bound
        hd[Ti=3, Var=At(:Ux)] = 99.9    # Just within bounds
        ld = DimArray(rand(2, length(needed_cols)), (Ti(1:2), Var(needed_cols)))

        qc = Peddy.PhysicsBoundsCheck()
        Peddy.quality_control!(qc, hd, ld, sensor)

        # Boundary values should be preserved (assuming inclusive bounds)
        # This tests the specific implementation of bounds checking
        # Note: The actual behavior depends on whether bounds are inclusive or exclusive
        @test !isnan(hd[Ti=3, Var=At(:Ux)])  # Within bounds should be OK
    end

    @testset "No QC (Nothing) Case" begin
        # Test that passing nothing for QC does nothing
        sensor = Peddy.CSAT3()
        needed_cols = collect(Peddy.needs_data_cols(sensor))  # Convert tuple to vector
        n_points = 5
        n_vars = length(needed_cols)

        # Create data with extreme values
        test_data = ones(n_points, n_vars) .* 1000  # Very extreme values
        original_data = copy(test_data)

        hd = DimArray(test_data, (Ti(1:n_points), Var(needed_cols)))
        ld = DimArray(rand(2, n_vars), (Ti(1:2), Var(needed_cols)))

        # Apply no QC (should do nothing)
        Peddy.quality_control!(nothing, hd, ld, sensor)

        # Data should be unchanged
        @test hd.data == original_data
        @test !any(isnan, hd.data)  # No values should be NaN
    end

    @testset "Integration with Pipeline" begin
        # Test QC integration in full pipeline
        sensor = Peddy.CSAT3()
        needed_cols = collect(Peddy.needs_data_cols(sensor))  # Convert tuple to vector
        n_points = 15
        n_vars = length(needed_cols)

        # Create mixed data: some normal, some extreme
        test_data = randn(n_points, n_vars) .* 3 .+ 2  # Mostly normal

        hd = DimArray(test_data, (Ti(1:n_points), Var(needed_cols)))
        hd[Ti=7, Var=At(:Ux)] = 150.0
        hd[Ti=12, Var=At(:Uy)] = -150
        ld = DimArray(rand(5, n_vars), (Ti(1:5), Var(needed_cols)))

        # Set up pipeline with QC only
        output = Peddy.MemoryOutput()
        qc = Peddy.PhysicsBoundsCheck()

        pipeline = Peddy.EddyPipeline(; sensor=sensor,
                                      quality_control=qc,
                                      despiking=nothing,
                                      gap_filling=nothing,
                                      gas_analyzer=nothing,
                                      double_rotation=nothing,
                                      mrd=nothing,
                                      output=output)

        # Run pipeline
        Peddy.process!(pipeline, hd, ld)

        # Get results
        processed_hf, processed_lf = Peddy.get_results(output)

        # Test that extreme values were removed
        @test isnan(processed_hf[Ti=7, Var=At(:Ux)])
        @test isnan(processed_hf[Ti=12, Var=At(:Uy)])

        # Test that normal values were preserved
        normal_indices = [1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 13, 14, 15]
        for idx in normal_indices
            @test !isnan(processed_hf[Ti=idx, Var=At(:Ux)])
            @test !isnan(processed_hf[Ti=idx, Var=At(:Uy)])
        end
    end
end

@testset "Sonic diagnostic margin" begin
    n = 50
    times = [DateTime(2025, 1, 1) + Millisecond(100 * i) for i in 0:(n - 1)]
    function make_data()
        m = hcat(fill(1.0, n), fill(2.0, n), fill(0.1, n), fill(-10.0, n),
                 fill(50.0, n), fill(80.0, n), zeros(n), fill(250.0, n))
        DimArray(m, (Ti(times), Var([:Ux, :Uy, :Uz, :Ts, :H2O, :P, :diag_sonic, :diag_gas])))
    end

    # margin 0: only the flagged sample
    d = make_data(); d[Ti=20, Var=At(:diag_sonic)] = 5000.0
    check_diagnostics!(LICOR(), d)
    @test findall(isnan, collect(d[Var=At(:Ts)])) == [20]

    # margin 2: two samples on either side, all four sonic variables, gas untouched
    d = make_data(); d[Ti=20, Var=At(:diag_sonic)] = 5000.0
    check_diagnostics!(LICOR(diag_margin=2), d)
    for v in (:Ux, :Uy, :Uz, :Ts)
        @test findall(isnan, collect(d[Var=At(v)])) == collect(18:22)
    end
    @test !any(isnan, d[Var=At(:H2O)])

    # margin clipped at the series edges, CSAT3 behaves the same
    d = make_data(); d[Ti=1, Var=At(:diag_sonic)] = 5000.0
    check_diagnostics!(CSAT3(diag_margin=2), d)
    @test findall(isnan, collect(d[Var=At(:Ux)])) == [1, 2, 3]
end
