using GLMakie

include("../include.jl")

function plot_observed_lsq_fft_fit(star, inc, R_m, W, n_freq, n_ζ, v_z_borders, Δv_z, obs_v_z, obs_profile)
    fig = Figure()

    ζs = [4i_ζ/n_ζ - 2/n_ζ for i_ζ = 1:(n_ζ÷4)]
    v_z_start, v_z_end = v_z_borders
    v_z_step = (v_z_end - v_z_start)/n_freq
    v_zs = [v_z_start + v_z_step*i_v_z - v_z_step/2 for i_v_z = 1:n_freq]

    obs_int = LinearInterpolation(obs_profile, obs_v_z)
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

    kernel_matrix = Observable(kernel_ζ_n_half(calc_emission_kernel_matrix(star, DipoleGeometry(R_m - W/2, R_m + W/2), 
                            Orientation(inc,0,0), Δv_z, n_freq, n_ζ, v_z_borders; n_Rm = 40, n_vz = 40),2))

    absorption_data = zeros(n_freq)

    geometry = DipoleGeometry(R_m - W/2, R_m + W/2)
    orientation = Orientation(inc)
    absorption_data = calc_absorption_profile(star, geometry, orientation, v_zs, Δv_z, 0.01, 0.01, hotspot_val[])
    
    # for i_freq = 1:n_freq
    #     print("\e[2K\e[1G abs calc: $i_freq")
    #     v_z_freq = v_z_start + v_z_step*i_freq - v_z_step/2
    #     absorption_data[i_freq] = TTauUtils.Profiles.calc_absorption_advanced(star, TTauUtils.GeometryAndOrientations.DipoleGeometry(R_m - W/2, R_m + W/2),
    #                                  Orientation(inc,0,0), v_z_freq, Δv_z; n_ζ = 4n_ζ, n_Rm = 20, n_vz = 20)
    # end

    absorption = Observable(absorption_data)

    on(ker_btn.clicks) do n
        inc, R_m, W = model_pars[]
        @time kernel_matrix[] = kernel_ζ_n_half(calc_emission_kernel_matrix(star, DipoleGeometry(R_m - W/2, R_m + W/2), 
                        Orientation(inc,0,0), Δv_z, n_freq, n_ζ, v_z_borders; n_Rm = 40, n_vz = 40),2)

        absorption_data = zeros(n_freq)

        geometry = DipoleGeometry(R_m - W/2, R_m + W/2)
        orientation = Orientation(inc)
        @time absorption_data = calc_absorption_profile(star, geometry, orientation, v_zs, Δv_z, 0.01, 0.01, hotspot_val[])

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

    source_function = x -> 10 ^ AkimaInterpolation([-1.5,-0.6,-0.2,-0.07,-0.03,0.0] .- 1.0, [0,0.2,0.3,0.6,0.7,1])(x)
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
        @time restored_source_data = fit_emission_profile(emission_assumed, kernel_matrix[], v_z_borders,
                                                                         pow[], damp[], center)
        points = zeros(2, n_ζ÷4)
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
        4*kernel_matrix[] / n_ζ * source_function_vals .+ 1.0
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

    ζ_step = 4/n_ζ
    freqs = fftshift(fftfreq(n_ζ, 1/ζ_step))
    max_freq = maximum(freqs)
    abs_freqs_normed = lift(damp, pow)  do damp, pow
        abs.(freqs/(max_freq)) .^ pow / damp
    end

    lines!(fft_ax, freqs[2:2:end], @lift log10.($fft_source)[2:2:end])
    lines!(fft_ax, freqs[2:2:end], @lift log10.($abs_freqs_normed)[2:2:end])
    # scatter!(source_ax, @lift(getindex($restored_source_points, 1, :)), @lift(getindex($restored_source_points, 2, :)))


    fig
end

function plot_observed_saved_kernel(star, obs_v_z, obs_profile, n_ζ; kernels_dir = "kernels")
    fig = Figure()

    incs, R_ins, Ws, Δv_zs, Hs = load_kernel_models_grid(; kernels_dir)
    v_z_borders, n_freq, n_ζ_calc = load_kernel_grid(; kernels_dir)

    n_inc, n_Rin, n_W, n_Δvz, n_H = length.([incs, R_ins, Ws, Δv_zs, Hs])

    n_ζ_half = round(Int, log2(n_ζ_calc/n_ζ))

    ζs = [i_ζ/n_ζ - 1/2n_ζ for i_ζ = 1:n_ζ]
    v_z_start, v_z_end = v_z_borders
    v_z_step = (v_z_end - v_z_start)/n_freq
    v_zs = [v_z_start + v_z_step*i_v_z - v_z_step/2 for i_v_z = 1:n_freq]

    obs_int = LinearInterpolation(obs_profile, obs_v_z, extrapolation = ExtrapolationType.Constant)
    observed_profile = obs_int.(v_zs)

    i_inc = n_inc÷2; i_Rin = n_Rin÷2; i_W = n_W÷2; i_Δvz = n_Δvz÷2; i_H = n_H÷2

    model_index = Observable([i_inc, i_Rin, i_W, i_Δvz, i_H])

    kernel_matrix = lift(model_index) do model_index
        kernel_ζ_half(kernel_ζ_half(load_kernel_from_index(model_index[1:end-1]...; kernels_dir)))
    end

    absorption = lift(model_index) do model_index
        load_absorption_from_index(model_index...; kernels_dir)
    end
    
    ker_ax = Axis(fig[1,1][1,1], title = "Model kernel")
    # ref_ker_ax = Axis(fig[1,2][1,1], title = "Reference kernel")
    
    mod_inc_sl = Slider(fig[1,2][1,2], range = 1:n_inc, startvalue = i_inc, update_while_dragging = true)
    mod_inc_sl_txt = Label(fig[1,2][1,1], text = @lift @sprintf("i = %.1f", incs[getindex($model_index, 1)]))

    mod_R_in_sl = Slider(fig[1,2][2,2], range = 1:n_Rin, startvalue = i_Rin, update_while_dragging = true)
    mod_R_in_sl_txt = Label(fig[1,2][2,1], text = @lift @sprintf("R_m = %.1f", R_ins[getindex($model_index, 2)]))

    mod_W_sl = Slider(fig[1,2][3,2], range = 1:n_W, startvalue = i_W, update_while_dragging = true)
    mod_W_sl_txt = Label(fig[1,2][3,1], text = @lift @sprintf("W = %.1f", Ws[getindex($model_index, 3)]))

    mod_Δv_z_sl = Slider(fig[1,2][4,2], range = 1:n_Δvz, startvalue = i_Δvz, update_while_dragging = true)
    mod_Δv_z_sl_txt = Label(fig[1,2][4,1], text = @lift @sprintf("Δv_z = %.1f", Δv_zs[getindex($model_index, 4)]))

    mod_H_sl = Slider(fig[1,2][5,2], range = 1:n_H, startvalue = i_H, update_while_dragging = true)
    mod_H_sl_txt = Label(fig[1,2][5,1], text = @lift @sprintf("H = %.1f", Hs[getindex($model_index, 5)]))

    on(mod_inc_sl.value) do mod_i_inc
        model_index[][1] = mod_i_inc
        notify(model_index)
    end

    on(mod_R_in_sl.value) do mod_i_Rin
        model_index[][2] = mod_i_Rin
        notify(model_index)
    end

    on(mod_W_sl.value) do mod_i_W
        model_index[][3] = mod_i_W
        notify(model_index)
    end

    on(mod_Δv_z_sl.value) do mod_i_Δvz
        model_index[][4] = mod_i_Δvz
        notify(model_index)
    end

    on(mod_H_sl.value) do mod_i_H
        model_index[][5] = mod_i_H
        notify(model_index)
    end

    # noise_btn = Button(fig[5,2][1,4], label = "Resample", tellwidth = false)

    kernel_matrix = Observable(kernel_ζ_n_half(load_kernel_from_index(i_inc, i_Rin, i_W, i_Δvz; kernels_dir),n_ζ_half))

    absorption_data = load_absorption_from_index(i_inc, i_Rin, i_W, i_Δvz, i_H; kernels_dir)

    absorption = Observable(absorption_data)

    on(model_index) do model_index
        i_inc, i_Rin, i_W, i_Δvz, i_H = model_index
        kernel_matrix[] = kernel_ζ_n_half(load_kernel_from_index(i_inc, i_Rin, i_W, i_Δvz; kernels_dir), n_ζ_half)

        absorption[] = load_absorption_from_index(i_inc, i_Rin, i_W, i_Δvz, i_H; kernels_dir)
    end

    hm_ker = heatmap!(ker_ax, v_zs, ζs,  (@lift log10.($kernel_matrix)), colorrange = (-1,2.5))
    Colorbar(fig[1,1][1,2], hm_ker)

    source_function = x -> 10 ^ AkimaInterpolation([-1.5,-0.6,-0.2,-0.07,-0.03,0.0] .- 1.0, [0,0.2,0.3,0.6,0.7,1])(x)
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

    on(absorption) do absorption
        emission_assumed[] = observed_profile + absorption
    end

    # emission_assumed = lift(profile_ref, absorption) do profile_ref, absorption
    #     emission_assumed[] = profile_ref + absorption[]
    # end

    pow = Observable(3.0)
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

    restored_source_points = lift(emission_assumed, pow, damp) do emission_assumed, pow, damp
        restored_source_data = fit_emission_profile(emission_assumed, kernel_matrix[], v_z_borders,
                                                                         pow, damp, center)
        points = zeros(2, n_ζ)
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
        kernel_matrix[] / n_ζ * source_function_vals .+ 1.0
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

    ζ_step = 1/n_ζ
    freqs = fftshift(fftfreq(4n_ζ, 1/ζ_step))
    max_freq = maximum(freqs)
    abs_freqs_normed = lift(damp, pow)  do damp, pow
        abs.(freqs/(max_freq)) .^ pow / damp
    end

    lines!(fft_ax, freqs[2:2:end], @lift log10.($fft_source)[2:2:end])
    lines!(fft_ax, freqs[2:2:end], @lift log10.($abs_freqs_normed)[2:2:end])
    # scatter!(source_ax, @lift(getindex($restored_source_points, 1, :)), @lift(getindex($restored_source_points, 2, :)))


    fig
end

source_known = open("../hart_to_fit/mag_3-4_15_32_S.bin", "r") do io
    read_array(io, Float64)
end

io = open("../hart_to_fit/mag_3-4_15_32_P.bin", "r") 
v_zs = reverse(read_array(io, Float64))
profile = reverse(read_array(io, Float64)) .+ 0.00*randn(length(v_zs))
profile[1] = 1.0; profile[end] = 1.0
close(io)

star = Star("CTTS", 2, 0.8, 4000, 15)

inc = 60

n=128
fig = plot_observed_saved_kernel(star, v_zs, profile, 32; kernels_dir="../kernels_bin")
# fig = plot_observed_lsq_fft_fit(star, 15, 3.5, 1.0, n, n, (-3.5e7, 3.5e7), 2e6, v_zs, profile)


R_ms = [3,3.5,4]
n_S = size(source_known)[2]
for i_Rm = 1:3
    R_m = R_ms[i_Rm]
    rs = collect(range(1,R_m,n_S))
    θs = asin.(sqrt.(rs/R_m))
    ζs = zeros(n_S)
    for i_ζ = 1:n_S
        r_m, ζ = spheretogridproper(DipoleGeometry(3,4), rs[i_ζ], θs[i_ζ])
        ζs[i_ζ] = ζ
    end
    lines!(fig.content[13], ζs, log10.(source_known[i_Rm,:]))
end

scene = display(fig)
wait(scene)