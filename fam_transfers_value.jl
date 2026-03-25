###############################################################################################
# School
###############################################################################################

function VSj(vjp1, model, j)
    (; beta, gamma, shock_resources) = model

    @threads for idx in eachindex(tasks_idx)
        (R, t, e, i_a, i_z) = tasks_idx[idx]
        
        for shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2

            resources = shock_resources[j, R, t, e, i_a, i_z, shock_in, shock_out, past_in, past_out]
        
            result = optimize(ap1 -> - (u(resources - ap1, gamma) + beta *  sj* EVS_jp1(model, vjp1, R, t, ap1, i_z)),
                            0.0, resources,  Brent(); rel_tol=1e-4, abs_tol=1e-4)

            VSj[R, t, e, i_a, i_z, shock_in, shock_out, past_in, past_out] = -result.minimum
            SPFj[R, t, e, i_a, i_z, shock_in, shock_out, past_in, past_out] = result.minimizer
            
        end
    end

    return VSj, SPFj
end

function EVS_jp1(model, vjp1, R, t, ap1, i_z)
    (; prob_shocks ) = model

    expected_value = 0.0
    for shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2
        expected_value += prob_shocks[jp1, R, m, n, e, i_z, shock_in, shock_out, past_in, past_out] * vjp1(R, t, ap1, i_z, shock_in, shock_out, past_in, past_out)
    end

    return expected_value
end

###############################################################################################
# Working
###############################################################################################

function Vj(vjp1, model, j)
    (; beta, gamma, shock_resources) = model

    # Create interpolation object
    v_itp = [LinearInterpolation((a_grid, z_grid[R]), vjp1[R, m, n, e, :, :, shock_in, shock_out, past_in, past_out], extrapolation_bc=Interpolations.Flat()) 
                                for R in Race, m in marital_status, n in fam_size, e in ed_type, shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2];

    @threads for idx in eachindex(tasks_idx)
        (R, m, n, e, i_a, i_z) = tasks_idx[idx]
        vjp1_itp = v_itp[R, m, n, e, :, :, :, :, :, :]

        for shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2

            resources = shock_resources[j, R, m, n, e, i_a, i_z, shock_in, shock_out, past_in, past_out]
        
            if j < 23
                result = optimize(ap1 -> - (u(resources - ap1, gamma) + beta * EV_jp1(model, vjp1_itp, j, R, m, n, e, past_in, past_out, ap1, i_z)),
                            0.0, resources,  Brent(); rel_tol=1e-4, abs_tol=1e-4)
            else
                result = optimize(ap1 -> - (u(resources - ap1, gamma) + beta * EW_jp1(model, vjp1_itp, j, R, m, n, e, past_in, past_out, ap1, i_z)),
                            0.0, resources,  Brent(); rel_tol=1e-4, abs_tol=1e-4)
            end

            Vj[R, m, n, e, i_a, i_z, shock_in, shock_out, past_in, past_out] = -result.minimum
            PFj[R, m, n, e, i_a, i_z, shock_in, shock_out, past_in, past_out] = result.minimizer
            
        end
    end

    return Vj, PFj
end

function EV_jp1(model, vjp1, j, R, m, n, e, past_in, past_out, ap1, i_z)
    (; Pimat, prob_shocks, zpnts) = model
    jp1 = j + 1

    expected_value = 0.0
    for i_zp1 in 1:zpnts
        pi_z =Pimat[i_z, i_zp1]
        for shock_in in 1:2, shock_out in 1:2, past_in_next in 1:2, past_out_next in 1:2
            if past_in == 2 && past_in_next == 1
                continue  
            end
            if past_out == 2 && past_out_next == 1
                continue 
            end
            expected_value += pi_z * prob_shocks[jp1, R, m, n, e, i_zp1, shock_in, shock_out, past_in_next, past_out_next] * vjp1[shock_in, shock_out, past_in_next, past_out_next](ap1, i_zp1)
        end
    end

    return expected_value
end

###############################################################################################
# Retirement
###############################################################################################

function Wj(wjp1, model, j)
    (; survival_risk, beta, gamma , shock_resources) = model

    # Create interpolation object
    if wjp1 != nothing
        w_itp = [LinearInterpolation((a_grid, z_grid[R]), wjp1[R, m, n, e, :, :, shock_in, shock_out, past_in, past_out], extrapolation_bc=Interpolations.Flat())
                for R in Race, m in marital_status, n in fam_size, e in ed_type, shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2]
    else
        w_itp = nothing
    end 

    @threads for idx in eachindex(tasks_idx)
        (R, m, n, e, i_a, i_z) = tasks_idx[idx]
        sj = survival_risk[j, R]
        wjp1_itp = w_itp[R, m, n, e, :, :, :, :, :, :]

        for shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2

            resources = shock_resources[j, R, m, n, e, i_a, i_z, shock_in, shock_out, past_in, past_out]

            if j == 34
                # Last period of life, no future value
                Wj[R, m, n, e, i_a, i_z, shock_in, shock_out, past_in, past_out] = u(resources, gamma)
                WPFj[R, m, n, e, i_a, i_z, shock_in, shock_out, past_in, past_out] = 0.0
            else
                result = optimize(ap1 -> - (u(resources - ap1, gamma) + beta *  sj* EW_jp1(model, wjp1_itp, j, R, m, n, e, past_in, past_out,ap1, i_z)),
                            0.0, resources,  Brent(); rel_tol=1e-4, abs_tol=1e-4)
                Wj[R, m, n, e, i_a, i_z, shock_in, shock_out, past_in, past_out] = -result.minimum
                WPFj[R, m, n, e, i_a, i_z, shock_in, shock_out, past_in, past_out] = result.minimizer
            end
        end
    end

    return Wj, WPFj
end

function EW_jp1(model, wjp1_itp, j, R, m, n, e, past_in, past_out, ap1, i_z)
    (; survival_risk, prob_shocks ) = model
    jp1 = j + 1
    sjp1 = survival_risk[jp1, R]

    expected_value = 0.0
    for shock_in in 1:2, shock_out in 1:2, past_in_next in 1:2, past_out_next in 1:2
        if past_in == 2 && past_in_next == 1
            continue 
        end
        if past_out == 2 && past_out_next == 1
            continue  
        end
        expected_value += prob_shocks[jp1, R, m, n, e, i_z, shock_in, shock_out, past_in_next, past_out_next] * wjp1_itp[shock_in, shock_out, past_in_next, past_out_next](ap1, i_z)
    end

    return expected_value * sjp1
end
