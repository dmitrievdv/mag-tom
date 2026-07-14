function fft_model_jacob!(J, lg_source, v_zs, kernel_matrix, n_ζ, n, planed_fft, center_window, abs_freqs_normed)
    J .= 0.0

    for i = 1:n
        if abs(v_zs[i]) < center_window
            J[i,:] .= 0.0
            continue
        end
        for j = 1:n_ζ
            J[i,j] = kernel_matrix[i,j]*log(10)*10^lg_source[j]/n_ζ
        end
    end

    fft_data = zeros(4n_ζ)
    fft_data[1:n_ζ] = lg_source
    fft_data[n_ζ+1:2n_ζ] = reverse(lg_source)
    fft_data[2n_ζ+1:4n_ζ] = 2fft_data[1] .- fft_data[1:2n_ζ] 

    fft_lg_source = planed_fft * fft_data
    fft_lg_source_abs = abs.(fft_lg_source)

    for i = 1:4n_ζ
        if fft_lg_source_abs[i] ≈ 0.0 
            continue
        end
        for j_1 = 2n_ζ+1:4n_ζ
            ang = 2π*(i-1)/4n_ζ*(j_1-1)
            J[n+i, 1] += 2(real(fft_lg_source[i])*cos(ang) - imag(fft_lg_source[i])*sin(ang))
        end
        for j = 1:n_ζ
            for j_2 in (j, 2n_ζ+1-j, 2n_ζ+j, 4n_ζ+1-j)
                ang = 2π*(i-1)/4n_ζ*(j_2-1)
                sign = (j_2 > 2n_ζ) ? -1 : 1
                J[n+i,j] += sign*(real(fft_lg_source[i])*cos(ang) - 
                                imag(fft_lg_source[i])*sin(ang))
            end
        end
        J[n+i,:] *= abs_freqs_normed[i]/fft_lg_source_abs[i]
    end
end

function fft_model!(out, lg_source, v_zs, kernel_matrix, n_ζ, n, planed_fft, center_window, emission_profile, abs_freqs_normed)
    source = 10 .^ lg_source
    model_profile = (kernel_matrix/n_ζ) * source .+ 1.0

    fft_data = vcat(lg_source, reverse(lg_source))
    fft_data = vcat(fft_data, 2fft_data[1] .- fft_data)

    fft_lg_source = planed_fft * fft_data
    fft_lg_source_abs = abs.(fft_lg_source)

    for i_v_z = 1:n
        if abs(v_zs[i_v_z]) < center_window
            out[i_v_z] = 0.0
        end
        out[i_v_z] = model_profile[i_v_z] - emission_profile[i_v_z]
    end
        
    for i_fft = 1:4n_ζ
        out[i_fft + n] = fft_lg_source_abs[i_fft] * abs_freqs_normed[i_fft]
    end
end

function fit_emission_profile(emission_profile, kernel_matrix, v_z_borders, pow = 1, damp = 20, center_window = 0)
    n, n_ζ = size(kernel_matrix)
    v_z_start, v_z_end = v_z_borders
    v_z_step = (v_z_end - v_z_start)/n
    v_zs = [v_z_start + v_z_step*i_v_z - v_z_step/2 for i_v_z = 1:n]
    ζs = [i_ζ/n_ζ - 1/2n_ζ for i_ζ = 1:n_ζ]
    ζ_step = 1/n_ζ

    lg_source_0 = zeros(n_ζ)

    lg_source_for_fft = vcat(lg_source_0, reverse(lg_source_0))
    lg_source_for_fft = vcat(lg_source_for_fft, -lg_source_for_fft .+ 2lg_source_for_fft[1])

    planed_fft = plan_fft(lg_source_for_fft)
    freqs = fftfreq(4*n_ζ, 1/ζ_step)
    max_freq = maximum(freqs)

    abs_freqs_normed = abs.(freqs/(max_freq)) .^ pow / damp

    

    prob = LeastSquaresProblem(x = lg_source_0, 
                              f! = (out, x) -> fft_model!(out, x, v_zs, kernel_matrix, n_ζ, n, planed_fft, center_window, emission_profile, abs_freqs_normed), 
                              g! = (J, x) -> fft_model_jacob!(J, x, v_zs, kernel_matrix, n_ζ, n, planed_fft, center_window, abs_freqs_normed),
                              output_length = n+4n_ζ)

    res = optimize!(prob, LevenbergMarquardt())
    lg_source_function_restored = res.minimizer
    return lg_source_function_restored
end

function fit_known(profiles_v_zs, profiles; kernels_dir = "kernels")
    incs, R_ins, Ws, Δv_zs, Hs = load_kernel_models_grid(; kernels_dir)
    v_z_borders, n_freq, n_ζ = load_kernel_grid(; kernels_dir)

    target_n_freq = 256
    target_n_ζ = 64

    n_ζ_half = round(Int, log2(n_ζ/target_n_ζ))
    n_freq_half = round(Int, log2(n_freq/target_n_freq))

    n_inc, n_Rin, n_W, n_Δvz, n_H = length.([incs, R_ins, Ws, Δv_zs, Hs])
    n_lines, n_prof = size(profiles)

    residuals = zeros(n_lines, n_inc, n_Rin, n_W, n_Δvz, n_H)
    profiles_model_grid = fill(1.0, (n_lines, n_freq))
     

    for i_line = 1:n_lines
        v_z_start, v_z_end = v_z_borders
        v_z_step = (v_z_end - v_z_start)/n_freq
        v_zs = [v_z_start + v_z_step*i_v_z - v_z_step/2 for i_v_z = 1:n_freq]

        profile_int = LinearInterpolation(profiles[i_line, :], profiles_v_zs[i_line, :])
       

        for i_vz in 1:n_freq
            v_z = v_zs[i_vz]
            if minimum(profiles_v_zs[i_line, :]) < v_z < maximum(profiles_v_zs[i_line, :])
                profiles_model_grid[i_line, i_vz] = profile_int(v_z)
            end
        end
        
    end

    for i_inc = 1:n_inc, i_Rin = 1:n_Rin, i_W = 1:n_W, i_Δvz = 1:n_Δvz
        kernel_matrix = kernel_ζ_n_half(load_kernel_from_index(i_inc, i_Rin, i_W, i_Δvz; kernels_dir), n_ζ_half)
        # println(size(kernel_matrix))
        for i_H = 1:n_H
            # print("")
            absorption = load_absorption_from_index(i_inc, i_Rin, i_W, i_Δvz, i_H; kernels_dir)
            Threads.@threads for i_line = 1:n_lines
                emission_profile = profiles_model_grid[i_line, :] + absorption
                lg_restored_source = fit_emission_profile(emission_profile, kernel_matrix, v_z_borders, 2.0, 1.0, 0.0)
                restored_profile = (kernel_matrix / target_n_ζ) * (10 .^ lg_restored_source) - absorption .+ 1.0
                δ = √((sum( (profiles_model_grid[i_line, :] - restored_profile) .^ 2 )) / n_freq)
            
                residuals[i_line, i_inc, i_Rin, i_W, i_Δvz, i_H] = δ
            end
        end
        print("\e[2K\e[1Gi = $(incs[i_inc]) R_in = $(R_ins[i_Rin]) W = $(Ws[i_W]) Δv_z = $(Δv_zs[i_Δvz])")
    end
        
    print("\n")
    return residuals
end