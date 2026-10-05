using Test
using Peddy
using DimensionalData
using Statistics
using Dates

@testset "Despiking Tests" begin
    @testset "VariableGroup Construction" begin
        # Test VariableGroup construction
        wind_group = Peddy.VariableGroup("Wind Components", [:Ux, :Uy, :Uz])
        @test wind_group.name == "Wind Components"
        @test wind_group.variables == [:Ux, :Uy, :Uz]
        @test wind_group.spike_threshold == 6.0  # Default threshold
        
        # Test VariableGroup with custom threshold
        temp_group = Peddy.VariableGroup("Temperature", [:Ts], spike_threshold=5.0)
        @test temp_group.name == "Temperature"
        @test temp_group.variables == [:Ts]
        @test temp_group.spike_threshold == 5.0
    end
    
    @testset "SimpleSigmundDespiking Construction" begin
        # Test default construction
        despiking = SimpleSigmundDespiking()
        @test despiking.window_minutes == 5.0
        @test length(despiking.variable_groups) == 1
        @test despiking.variable_groups[1].name == "Default Sonic"
        @test despiking.variable_groups[1].variables == [:Ux, :Uy, :Uz, :Ts]
        @test despiking.variable_groups[1].spike_threshold == 6.0
        
        # Test custom construction with variable groups
        wind_group = Peddy.VariableGroup("Wind", [:Ux, :Uy, :Uz], spike_threshold=6.0)
        temp_group = Peddy.VariableGroup("Temperature", [:Ts], spike_threshold=5.0)
        
        custom_despiking = SimpleSigmundDespiking(
            window_minutes=3.0,
            variable_groups=[wind_group, temp_group]
        )
        @test custom_despiking.window_minutes == 3.0
        @test length(custom_despiking.variable_groups) == 2
        @test custom_despiking.variable_groups[1].name == "Wind"
        @test custom_despiking.variable_groups[2].name == "Temperature"
        
        # Test that it implements AbstractDespiking
        @test despiking isa Peddy.AbstractDespiking
    end
    
    @testset "Fast Median Internal" begin
        # Test _fast_median! internal helper
        # Case 1: single element
        buf1 = [10.0]
        @test Peddy._fast_median!(buf1, 1) == 10.0
        
        # Case 2: two elements
        buf2 = [10.0, 20.0]
        @test Peddy._fast_median!(buf2, 2) == 15.0
        
        # Case 3: odd elements
        buf3 = [30.0, 10.0, 20.0, 40.0, 50.0]
        @test Peddy._fast_median!(buf3, 5) == 30.0
        
        # Case 4: even elements
        buf4 = [30.0, 10.0, 20.0, 40.0]
        @test Peddy._fast_median!(buf4, 4) == 25.0
        
        # Case 5: partial buffer
        buf5 = [3.0, 1.0, 2.0, 100.0, 100.0]
        @test Peddy._fast_median!(buf5, 3) == 2.0
    end

    @testset "Window Size Calculation" begin
        # Create dummy data with 10Hz sampling (2000 points)
        times = [DateTime(2023, 1, 1) + Millisecond(i * 100) for i in 0:1999]
        data = DimArray(randn(2000, 1), (Ti(times), Var([:Ux])))
        
        # 1 minute window at 10Hz = 600 points
        # 600 < 2000/3 (666.6), so it shouldn't be adjusted
        @test Peddy._calculate_window_size(data, 1.0) == 600
        
        # Limited data case (1/3 of total)
        # 20 minute window would be 12000 points, but total is 2000. 
        # Should be capped at 2000 ÷ 3 = 666
        @test_logs (:warn, r"Adjusted window size") Peddy._calculate_window_size(data, 20.0) == 666
    end
    
    @testset "Pattern Deviation Calculation" begin
        # Test pattern deviation calculation
        df_di = [1.0, 2.0, 10.0, 3.0, 4.0]  # Spike at position 3
        df_hat = Peddy.calculate_pattern_deviation(df_di)
        
        # The spike should have higher pattern deviation
        @test df_hat[3] > df_hat[1]
        @test df_hat[3] > df_hat[2]
        @test df_hat[3] > df_hat[4]
        @test df_hat[3] > df_hat[5]
        
        # Test boundary conditions
        @test !isnan(df_hat[1])  # First element
        @test !isnan(df_hat[end])  # Last element
        
        # Test with uniform deviations (no spikes) - should be zero deviations
        uniform_deviations = [0.0, 0.0, 0.0, 0.0, 0.0]  # No deviations from median
        uniform_hat = Peddy.calculate_pattern_deviation(uniform_deviations)
        @test all(x -> abs(x) < 1e-10, uniform_hat)  # Should be near zero
    end
    
    @testset "Spike Detection with Synthetic Data" begin
        # Create synthetic data with known spikes
        n_points = 1000
        freq_hz = 10.0  # 10 Hz sampling
        dt_ms = round(Int, 1000 / freq_hz)  # Millisecond
        
        # Create time dimension
        start_time = DateTime(2023, 1, 1, 12, 0, 0)
        times = [start_time + Millisecond(i * dt_ms) for i in 0:(n_points-1)]
        
        # Create clean synthetic data
        clean_ux = 2.0 .+ 0.5 .* sin.(2π .* (1:n_points) ./ (60 * freq_hz))  # 1-minute oscillation
        clean_uy = 1.0 .+ 0.3 .* cos.(2π .* (1:n_points) ./ (30 * freq_hz))  # 30-second oscillation
        clean_uz = 0.1 .+ 0.1 .* randn(n_points)  # Small random variations
        clean_ts = 20.0 .+ 2.0 .* sin.(2π .* (1:n_points) ./ (120 * freq_hz))  # 2-minute temperature cycle
        
        # Add spikes at known locations
        spike_indices = [100, 300, 500, 700]
        spike_ux = copy(clean_ux)
        spike_uy = copy(clean_uy)
        spike_uz = copy(clean_uz)
        spike_ts = copy(clean_ts)
        
        # Add large spikes
        for idx in spike_indices
            spike_ux[idx] += 10.0  # Large positive spike
            spike_uy[idx] -= 8.0   # Large negative spike
            spike_uz[idx] += 5.0   # Medium spike
            spike_ts[idx] += 15.0  # Temperature spike
        end
        
        # Create DimArray in proper format
        test_data = hcat(spike_ux, spike_uy, spike_uz, spike_ts)
        high_freq_data = DimArray(
            test_data,
            (Ti(times), Var([:Ux, :Uy, :Uz, :Ts]))
        )
        
        # Apply despiking with custom variable group
        sonic_group = Peddy.VariableGroup("Sonic Variables", [:Ux, :Uy, :Uz, :Ts], spike_threshold=4.0)
        despiking = SimpleSigmundDespiking(window_minutes=2.0, variable_groups=[sonic_group])
        
        # Count original NaN values
        original_nan_count = sum(isnan.(spike_ux)) + sum(isnan.(spike_uy)) + 
                           sum(isnan.(spike_uz)) + sum(isnan.(spike_ts))
        
        # Apply despiking (this should work without low_frequency_data for this test)
        low_freq_data = nothing
        despike!(despiking, high_freq_data, low_freq_data)
        
        # Count NaN values after despiking
        final_nan_count = sum(isnan.(high_freq_data[Var=At(:Ux)])) + sum(isnan.(high_freq_data[Var=At(:Uy)])) + 
                         sum(isnan.(high_freq_data[Var=At(:Uz)])) + sum(isnan.(high_freq_data[Var=At(:Ts)]))
        
        # Should have detected and removed some spikes
        @test final_nan_count > original_nan_count
        
        # Check that at least some of the known spike locations were detected
        detected_spikes = 0
        for idx in spike_indices
            if isnan(high_freq_data[Var=At(:Ux)][idx]) || isnan(high_freq_data[Var=At(:Uy)][idx]) || 
               isnan(high_freq_data[Var=At(:Uz)][idx]) || isnan(high_freq_data[Var=At(:Ts)][idx])
                detected_spikes += 1
            end
        end
        
        @test detected_spikes > 0  # Should detect at least some spikes
        println("Detected $detected_spikes out of $(length(spike_indices)) artificial spikes")
    end
    
    @testset "H2O Variable Processing" begin
        # Create test data with H2O variable
        n_points = 500
        freq_hz = 10.0
        dt_ms = round(Int, 1000 / freq_hz)
        
        start_time = DateTime(2023, 1, 1, 12, 0, 0)
        times = [start_time + Millisecond(i * dt_ms) for i in 0:(n_points-1)]
        
        # Create synthetic data
        ux_data = 2.0 .+ 0.5 .* randn(n_points)
        uy_data = 1.0 .+ 0.3 .* randn(n_points)
        uz_data = 0.1 .+ 0.1 .* randn(n_points)
        ts_data = 20.0 .+ 2.0 .* randn(n_points)
        h2o_data = 15.0 .+ 3.0 .* randn(n_points)
        
        # Add spikes to H2O data
        h2o_data[100] += 50.0  # Large H2O spike
        h2o_data[200] -= 30.0  # Large negative H2O spike
        
        test_data = hcat(ux_data, uy_data, uz_data, ts_data, h2o_data)
        high_freq_data = DimArray(
            test_data,
            (Ti(times), Var([:Ux, :Uy, :Uz, :Ts, :LI_H2Om_corr]))
        )
        
        # Apply despiking with separate groups for sonic and H2O
        sonic_group = Peddy.VariableGroup("Sonic", [:Ux, :Uy, :Uz, :Ts], spike_threshold=3.0)
        h2o_group = Peddy.VariableGroup("H2O", [:LI_H2Om_corr], spike_threshold=3.0)
        despiking = SimpleSigmundDespiking(window_minutes=1.0, variable_groups=[sonic_group, h2o_group])
        low_freq_data = nothing
        
        original_h2o_nan = sum(isnan.(h2o_data))
        despike!(despiking, high_freq_data, low_freq_data)
        final_h2o_nan = sum(isnan.(high_freq_data[Var=At(:LI_H2Om_corr)]))
        
        # Should have detected H2O spikes
        @test final_h2o_nan > original_h2o_nan
        println("H2O spikes detected: $(final_h2o_nan - original_h2o_nan)")
    end
    
    @testset "Edge Cases and Error Handling" begin
        # Test with insufficient data points
        n_points = 5
        times = [DateTime(2023, 1, 1, 12, 0, 0) + Second(i) for i in 1:n_points]
        
        test_data = hcat([1.0, 2.0, 3.0, 4.0, 5.0], [1.0, 2.0, 3.0, 4.0, 5.0])
        small_data = DimArray(
            test_data,
            (Ti(times), Var([:Ux, :Uy]))
        )
        
        despiking = SimpleSigmundDespiking(window_minutes=10.0)  # Very large window
        low_freq_data = nothing
        
        # Should handle gracefully (may generate warnings about window size)
        @test_logs (:warn, r"Adjusted window size") despike!(despiking, small_data, low_freq_data)
        
        # Test with missing variables
        test_data = reshape([1.0, 2.0, 3.0, 4.0, 5.0], 5, 1)
        partial_data = DimArray(
            test_data,
            (Ti(times), Var([:Ux]))
        )
        
        @test_logs (:warn, r"Adjusted window size") despike!(despiking, partial_data, low_freq_data)
        
        # Test with all NaN data
        test_data = hcat(fill(NaN, 100), fill(NaN, 100))
        nan_data = DimArray(
            test_data,
            (Ti([DateTime(2023, 1, 1, 12, 0, 0) + Second(i) for i in 1:100]), Var([:Ux, :Uy]))
        )
        
        @test_logs (:warn, r"Adjusted window size") despike!(despiking, nan_data, low_freq_data)
    end
    
    @testset "MAD Floor Functionality" begin
        n = 100
        times = [DateTime(2023, 1, 1) + Second(i) for i in 1:n]
        # Data with zero variance (all 1.0)
        data_raw = fill(1.0, n)
        # Add one spike
        data_raw[50] = 10.0
        
        da = DimArray(reshape(copy(data_raw), n, 1), (Ti(times), Var([:Ux])))
        
        # Without MAD floor, zero variance might lead to issues or very sensitive detection
        # With MAD floor = 1.0, the threshold will be much higher
        despiking = SimpleSigmundDespiking(window_minutes=1.0, use_mad_floor=true, 
                                          mad_floor=Dict(:Ux => 1.0))
        
        despike!(despiking, da, nothing)
        
        # The spike at 50 is 9.0 deviation. 
        # Normalized threshold is 6.0 / 0.6745 ≈ 8.89
        # If MAD floor is 1.0, deviation 9.0 / 1.0 = 9.0 >= 8.89 -> Spike detected
        @test isnan(da[Ti=At(times[50]), Var=At(:Ux)])
        
        # Test with even higher floor that prevents detection
        # IMPORTANT: Use a FRESH copy of data_raw
        da2 = DimArray(reshape(copy(data_raw), n, 1), (Ti(times), Var([:Ux])))
        despiking2 = SimpleSigmundDespiking(window_minutes=1.0, use_mad_floor=true, 
                                           mad_floor=Dict(:Ux => 2.0))
        despike!(despiking2, da2, nothing)
        # 9.0 / 2.0 = 4.5 < 8.89 -> No spike detected
        @test !isnan(da2[Ti=At(times[50]), Var=At(:Ux)])
    end

    @testset "SimpleSigmundDespiking with Float32" begin
        # Test constructor with Float32
        despiking_f32 = SimpleSigmundDespiking(number_type=Float32)
        @test despiking_f32 isa SimpleSigmundDespiking{Float32}
        @test despiking_f32.window_minutes isa Float32
        
        # Test despiking with Float32 data
        n = 100
        times = [DateTime(2023, 1, 1) + Second(i) for i in 1:n]
        data = rand(Float32, n, 1)
        data[50] = 100.0f0 # Huge spike
        
        da = DimArray(data, (Ti(times), Var([:Ux])))
        despike!(despiking_f32, da, nothing)
        
        @test isnan(da[Ti=At(times[50]), Var=At(:Ux)])
        @test eltype(da) == Float32
    end

    @testset "Logger Integration" begin
        n = 100
        times = [DateTime(2023, 1, 1) + Second(i) for i in 1:n]
        data = randn(n)
        data[20] = 50.0 # Spike
        
        da = DimArray(reshape(data, n, 1), (Ti(times), Var([:Ux])))
        
        logger = Peddy.ProcessingLogger()
        despiking = SimpleSigmundDespiking(window_minutes=1.0)
        
        despike!(despiking, da, nothing; logger=logger)
        
        @test length(logger.entries) > 0
        @test any(e -> e.stage == :despiking && e.category == :spike, logger.entries)
        
        # Check that it logged the correct group
        spike_entries = filter(e -> e.category == :spike, logger.entries)
        @test !isempty(spike_entries)
        @test spike_entries[1].details[:group] == "Default Sonic"
    end

    @testset "Small Data Pattern Deviation" begin
        # Test n < 3 case
        d1 = [1.0]
        @test Peddy.calculate_pattern_deviation(d1) == [0.5]
        d2 = [1.0, 2.0]
        @test Peddy.calculate_pattern_deviation(d2) == [0.5, 1.0]
    end

    @testset "Rolling Median All NaNs" begin
        data = [NaN, NaN, NaN]
        res = Peddy.calculate_rolling_median(data, 3)
        @test all(isnan, res)
    end
end

@testset "Ts rolling-median check" begin
    n = 3000  # 5 min at 10 Hz
    times = [DateTime(2025, 4, 25) + Millisecond(100 * i) for i in 0:(n - 1)]
    function make_data()
        m = hcat(fill(5.0, n), fill(1.0, n), 0.1 .* sin.(1:n), -20.0 .+ 0.01 .* sin.(1:n))
        DimArray(m, (Ti(times), Var([:Ux, :Uy, :Uz, :Ts])))
    end
    window = Peddy._calculate_window_size(make_data(), 5.0)

    d = make_data(); d[Ti=1000, Var=At(:Ts)] += 12.0; d[Ti=2000, Var=At(:Ts)] += 8.0
    nremoved = Peddy._remove_ts_median_outliers!(d, window, 10.0, nothing, nothing)
    @test nremoved == 1
    for v in (:Ux, :Uy, :Uz, :Ts)
        @test isnan(d[Ti=1000, Var=At(v)])
        @test !isnan(d[Ti=2000, Var=At(v)])
    end

    # option off by default: despike! does not run the check
    @test SimpleSigmundDespiking().max_ts_deviation === nothing
    @test SimpleSigmundDespiking(max_ts_deviation=10).max_ts_deviation == 10.0
end
