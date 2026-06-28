function plot_observed_lsq_fft_fit(star, inc, R_m, W, n_freq, n_ζ, v_z_borders, Δv_z, obs_v_z, obs_profile)
    fig = Figure()

    ζs = [2i_ζ/n_ζ - 1/n_ζ for i_ζ = 1:(n_ζ÷2)]
    v_z_start, v_z_end = v_z_borders
    v_z_step = (v_z_end - v_z_start)/n_freq
    v_zs = [v_z_start + v_z_step*i_v_z - v_z_step/2 for i_v_z = 1:n_freq]

    obs_int = Spline1D(obs_v_z*1e5, obs_profile, k = 1)
    observed_profile = obs_int.(v_zs)
    
    ker_ax = Axis(fig[1,1][1,1], title = "Model kernel")
    # ref_ker_ax = Axis(fig[1,2][1,1], title = "Reference kernel")
    
    model_pars = Observable(Float64[inc, R_m, W])
    hotspot_val = Observable(1.0)
    
    inc_sl = Slider(fig[1,2][1,2], range = 5:85, startvalue = inc)
    inc_sl_txt = Label(fig[1,2][1,1], text = @lift @sprintf("i = %.0f", getindex($model_pars, 1)))

    R_m_sl = Slider(fig[1,2][2,2], range = 2:0.1:7, startvalue = float(R_m))
    R_m_sl_txt = Label(fig[1,2][2,1], text = @lift @sprintf("R_m = %.1f", getindex($model_pars, 2)))

    W_sl = Slider(fig[1,2][3,2], range = 0.01:0.01:2, startvalue = float(W))
    W_sl_txt = Label(fig[1,2][3,1], text = @lift @sprintf("W = %.2f", getindex($model_pars, 3)))

    hotspot_sl = Slider(fig[1,2][4,2], range = 0:0.01:2, startvalue = 1.0)
    W_sl_txt = Label(fig[1,2][4,1], text = @lift @sprintf("H = %.2f", $hotspot_val))

    on(inc_sl.value) do sl_inc
        model_pars[] = [sl_inc, model_pars[][2], model_pars[][3]]
    end

    on(R_m_sl.value) do sl_R_m
        model_pars[] = [model_pars[][1], sl_R_m, model_pars[][3]]
    end

    on(W_sl.value) do sl_W
        model_pars[] = [model_pars[][1], model_pars[][2], sl_W]
    end

    on(hotspot_sl.value) do sl_H
        hotspot_val[] = 10^sl_H
    end

    ker_btn = Button(fig[1,2][5,1], label = "Model", tellwidth = false)
    fit_btn = Button(fig[1,2][5,2], label = "Fit", tellwidth = false)
    # noise_btn = Button(fig[5,2][1,4], label = "Resample", tellwidth = false)

    kernel_matrix = Observable(kernel_ζ_half(calc_emission_kernel_matrix(star, TTauUtils.GeometryAndOrientations.DipoleGeometry(R_m - W/2, R_m + W/2), 
                            Orientation(inc,0,0), Δv_z, n_freq, n_ζ, v_z_borders; n_Rm = 40, n_vz = 40, kernel_smooth = 0)))

    absorption_data = zeros(n_freq)

    geometry = TTauUtils.GeometryAndOrientations.DipoleGeometry(R_m - W/2, R_m + W/2)
    orientation = Orientation(inc)
    absorption_data = TTauUtils.Profiles.calc_absorption_profile(star, geometry, orientation, v_zs, Δv_z, 0.01, 0.01, hotspot_val[])
    
    # for i_freq = 1:n_freq
    #     print("\e[2K\e[1G abs calc: $i_freq")
    #     v_z_freq = v_z_start + v_z_step*i_freq - v_z_step/2
    #     absorption_data[i_freq] = TTauUtils.Profiles.calc_absorption_advanced(star, TTauUtils.GeometryAndOrientations.DipoleGeometry(R_m - W/2, R_m + W/2),
    #                                  Orientation(inc,0,0), v_z_freq, Δv_z; n_ζ = 4n_ζ, n_Rm = 20, n_vz = 20)
    # end

    absorption = Observable(absorption_data)

    on(ker_btn.clicks) do n
        inc, R_m, W = model_pars[]
        @time kernel_matrix[] = kernel_ζ_half(calc_emission_kernel_matrix(star, TTauUtils.GeometryAndOrientations.DipoleGeometry(R_m - W/2, R_m + W/2), 
                        Orientation(inc,0,0), Δv_z, n_freq, n_ζ, v_z_borders; n_Rm = 40, n_vz = 40, kernel_smooth = 0))

        absorption_data = zeros(n_freq)

        geometry = TTauUtils.GeometryAndOrientations.DipoleGeometry(R_m - W/2, R_m + W/2)
        orientation = Orientation(inc)
        @time absorption_data = TTauUtils.Profiles.calc_absorption_profile(star, geometry, orientation, v_zs, Δv_z, 0.01, 0.01, hotspot_val[])

        # for i_freq = 1:n_freq
        #     print("\e[2K\e[1G abs calc: $i_freq")
        #     v_z_freq = v_z_start + v_z_step*i_freq - v_z_step/2
        #     absorption_data[i_freq] = TTauUtils.Profiles.calc_absorption_advanced(star, TTauUtils.GeometryAndOrientations.DipoleGeometry(R_m - W/2, R_m + W/2),
        #                              Orientation(inc,0,0), v_z_freq, Δv_z; n_ζ = 4n_ζ, n_Rm = 20, n_vz = 20)
        # end

        absorption[] = absorption_data
    end

    hm_ker = heatmap!(ker_ax, v_zs, ζs,  (@lift log10.($kernel_matrix)), colorrange = (-1,2.5))
    Colorbar(fig[1,1][1,2], hm_ker)

    source_function = x -> 10 ^ Spline1D([0,0.2,0.3,0.6,0.7,1], [-1.5,-0.6,-0.2,-0.07,-0.03,0.0] .- 1.0, k = 1)(x)
    source_function_initial = Observable(source_function.(ζs))

    # emission_ref = @lift ($kernel_matrix_ref/n_ζ) * $source_function_initial
    # # println(size(to_value(emission_ref)))

    source_ax = Axis(fig[2,1], title = "Source function")
    profile_ax = Axis(fig[2,2], title = "Profiles")

    # on(noise_btn.clicks) do n
    #     noise_observation[] = 1/signal_to_noise[]*(2*randn(n_freq) .- 1)
    # end

    # profile_ref = @lift $emission_ref - $absorption_ref .+ 1 + $noise_observation

    # lines!(source_ax, ζs, @lift(log10.($source_function_initial)))
    lines!(profile_ax, v_zs, observed_profile)
    lines!(profile_ax, v_zs, @lift(1 .- $absorption), linestyle = :dash)

    source_width = lift(source_ax.scene.viewport) do viewport
        widths(viewport)[1]
    end

    emission_assumed = Observable(observed_profile + absorption[])

    on(fit_btn.clicks) do n
        emission_assumed[] = observed_profile + absorption[]
    end

    # emission_assumed = lift(profile_ref, absorption) do profile_ref, absorption
    #     emission_assumed[] = profile_ref + absorption[]
    # end

    pow = Observable(2.0)
    damp = Observable(1.0)

    damp_sl = Slider(fig[3,1][1,2], range = -1:0.005:3, startvalue = log10.(damp[]), tellheight = false)
    damp_sl_txt = Label(fig[3,1][1,1], text = @lift( @sprintf("damp = %.2f", $damp)), tellheight = false)

    pow_sl = Slider(fig[3,1][2,2], range = 0:0.1:5, startvalue = pow[], tellheight = false)
    pow_sl_txt = Label(fig[3,1][2,1], text = @lift( @sprintf("pow = %.2f", $pow)), tellheight = false)

    on(damp_sl.value) do damp_val
        damp[] = 10^damp_val
    end

    on(pow_sl.value) do pow_val
        pow[] = pow_val
    end

    center = 0*1e5

    restored_source_points = lift(emission_assumed) do emission_assumed
        @time restored_source_data = fft_regularization_least_squares_fit_source_emission_profile_no_plot(emission_assumed, kernel_matrix[], v_z_borders,
                                                                         pow[], damp[], center)
        points = zeros(2, n_ζ÷2)
        points[1,:] = ζs
        points[2,:] = restored_source_data
        points
    end

    source_function_vals = lift(restored_source_points) do source_points
        10 .^ source_points[2,:]
        # source_ζs = source_points[1,:]
        # sources = source_points[2,:]
        # source_function = x -> 10 ^ Spline1D(source_ζs, sources, k = 2)(x)
        # source_function_vals = source_function.(ζs)
    end

    emission_mod = lift(source_function_vals) do source_function_vals
        2*kernel_matrix[] / n_ζ * source_function_vals .+ 1.0
    end

    profile_mod = lift(emission_mod) do emission_mod
        emission_mod - absorption[]
    end

    δ = lift(profile_mod) do profile_mod
        √( 1/(n_freq-1)*sum((profile_mod[abs.(v_zs) .> center] - observed_profile[abs.(v_zs) .> center]) .^ 2))
    end

    on(δ) do δ
        profile_ax.title = "Profiles, δ = $δ"
    end

    lines!(source_ax, ζs, @lift(log10.($source_function_vals)))
    lines!(profile_ax, v_zs, profile_mod)
    lines!(profile_ax, v_zs, emission_mod, linestyle = :dash)

    fft_ax = Axis(fig[3,2])

    fft_source = lift(source_function_vals) do source_function 
        fft_data = vcat(source_function, reverse(source_function))
        fft_data = vcat(fft_data, 2fft_data[1] .- fft_data)
        (abs.(fftshift(fft(fft_data))))
    end

    ζ_step = 2/n_ζ
    freqs = fftshift(fftfreq(2*n_ζ, 1/ζ_step))
    max_freq = maximum(freqs)
    abs_freqs_normed = lift(damp, pow)  do damp, pow
        abs.(freqs/(max_freq)) .^ pow / damp
    end

    lines!(fft_ax, freqs[2:2:end], @lift log10.($fft_source)[2:2:end])
    lines!(fft_ax, freqs[2:2:end], @lift log10.($abs_freqs_normed)[2:2:end])
    # scatter!(source_ax, @lift(getindex($restored_source_points, 1, :)), @lift(getindex($restored_source_points, 2, :)))


    fig
end