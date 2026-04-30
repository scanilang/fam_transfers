# Value function for those that chose not to attend college
# Positive assets only
# e = 0 no college

###############################################################################################
# Working
###############################################################################################

function Vnc1j_solve(vnc1jp1, vnc2jp1, Vj_nc1, PFj_nc1, model, j)
    (; beta, gamma,a_grid_nocollege, tasks_idx_nc1, shock_resources_nc, z_grid, fam_shock_period, Race, marital_status, fam_size, fam_type) = model

    fill!(Vj_nc1, 0f0)
    fill!(PFj_nc1, 0f0)

    # When creating interpolation objects for Vj+1:
    

    if j < fam_shock_period
        vnc2_itp = nothing
        vnc1jp1_itp = [LinearInterpolation((a_grid_nocollege, z_grid[R]), vnc1jp1[R, t, :, :, shock_in, shock_out, past_in, past_out], extrapolation_bc=Flat()) 
           for R in Race, t in fam_type, shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2]
    else
        vnc2_itp = [LinearInterpolation((a_grid_nocollege, z_grid[R]), vnc2jp1[R, m, n, t, :, :, shock_in, shock_out, past_in, past_out], extrapolation_bc=Flat()) 
               for R in Race, m in marital_status, n in fam_size, t in fam_type, shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2]
        vnc1jp1_itp = nothing
    end

    @threads for idx in eachindex(tasks_idx_nc1)
        (R, t, i_a, i_z) = tasks_idx_nc1[idx]

        if j < fam_shock_period
            vnc1_itp = vnc1jp1_itp[R, t, :, :, :, :]

        else
            vnc1_itp = nothing
        end

        for shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2

            # hardcode m=1, n=1 since single and no family
            net_resources = shock_resources_nc[j, R, 1, 1, t, i_a, i_z, shock_in, shock_out, past_in, past_out]

            if j < fam_shock_period
                result = optimize(ap1 -> -(u(net_resources - ap1, gamma) + 
                     beta * EVnc1_jp1(model, vnc1_itp, j, R, t, ap1, i_z, shock_in, shock_out, past_in, past_out)),
                     0.0, net_resources,Brent(); rel_tol=1e-4, abs_tol=1e-4)
            else
                result = optimize(ap1 -> -(u(net_resources - ap1, gamma) + 
                     beta * EVnc_family_jp1(model, vnc2_itp, j, R, t, ap1, i_z, shock_in, shock_out, past_in, past_out)),
                     0.0, net_resources,Brent(); rel_tol=1e-4, abs_tol=1e-4)
 
            end
        
            Vj_nc1[R, t, i_a, i_z, shock_in, shock_out, past_in, past_out] = -result.minimum
            PFj_nc1[R, t, i_a, i_z, shock_in, shock_out, past_in, past_out] = result.minimizer
        end
    end
    
    return Vj_nc1, PFj_nc1
end

function EVnc1_jp1(model, vjp1, j, R, t, ap1, i_z, shock_in, shock_out, past_in, past_out)
    (; Pimat, zpnts, y_values, ra_w, ra_b) = model
    jp1 = j + 1
    a_income = R == 1 ? ap1 * ra_w : ap1 * ra_b

    expected_value = 0.0
    for i_zp1 in 1:zpnts
        pi_z =Pimat[R][i_z, i_zp1]
        y = y_values[R, jp1, 1, 1, i_zp1]
        for shock_in_next in 1:2, shock_out_next in 1:2, past_in_next in 1:2, past_out_next in 1:2
            # update past in and past out flags based on past flags and past shocks
            if past_in == 2 && past_in_next == 1 || shock_in == 2 && past_in_next == 1 
                continue  
            end
            if past_out == 2 && past_out_next == 1 || shock_out == 2 && past_out_next == 1 
                continue 
            end

            # probability of shock in and shock out next period
            shock_out_prob= shocks_out_prob(R,1,1,jp1,y, a_income, 1, t, past_in-1, past_out-1)
            shock_in_prob = shocks_in_prob(R,1,1,jp1,y, a_income, 1, t, past_in-1, past_out-1)
            if shock_in_next == 2
                shock_in_next_prob = shock_in_prob
            else
                shock_in_next_prob = 1 - shock_in_prob
            end
            if shock_out_next == 2
                shock_out_next_prob = shock_out_prob
            else
                shock_out_next_prob = 1 - shock_out_prob
            end

            # expected value contribution from next period state
            expected_value += pi_z * shock_out_next_prob *  shock_in_next_prob* vjp1[shock_in_next, shock_out_next, past_in_next, past_out_next](ap1, i_zp1)
        end
    end

    return expected_value
end


# Expected value with family transition
function EVnc_family_jp1(model, vnc2_itp, j, R, t, ap1, i_z, shock_in, shock_out, past_in, past_out)
    (; Pimat, zpnts, y_values, ra_w, ra_b) = model
    jp1 = j + 1
    a_income = R == 1 ? ap1 * ra_w : ap1 * ra_b

    outcomes = family_shock_probs[(R, e_college(1))]

    expected_value = 0.0
    for (m_next, n_next, t_next, prob_fam) in outcomes
        if t_next !== t
            continue
        end

        for i_zp1 in 1:zpnts
            pi_z =Pimat[R][i_z, i_zp1]
            y = y_values[R, jp1, m_next, 1, i_zp1]
            for shock_in_next in 1:2, shock_out_next in 1:2, past_in_next in 1:2, past_out_next in 1:2
                # update past in and past out flags based on past flags and past shocks
                if past_in == 2 && past_in_next == 1 || shock_in == 2 && past_in_next == 1 || past_out == 2 && past_out_next == 1 || shock_out == 2 && past_out_next == 1
                    continue  
                end

                # probability of shock in and shock out next period
                shock_out_prob= shocks_out_prob(R,n_next,m_next,jp1,y, a_income, 1, t_next, past_in-1, past_out-1)
                shock_in_prob = shocks_in_prob(R,n_next,m_next,jp1,y, a_income, 1, t_next, past_in-1, past_out-1)
                if shock_in_next == 2
                    shock_in_next_prob = shock_in_prob
                else
                    shock_in_next_prob = 1 - shock_in_prob
                end
                if shock_out_next == 2
                    shock_out_next_prob = shock_out_prob    
                else
                    shock_out_next_prob = 1 - shock_out_prob
                end

                # expected value contribution from next period state
                expected_value += prob_fam * pi_z * shock_out_next_prob *  shock_in_next_prob* vnc2_itp[R, m_next, n_next, t_next, shock_in_next, shock_out_next, past_in_next, past_out_next](ap1, i_zp1)
       
            end
        end
    end

    return expected_value
end

function Vnc2j_solve(vnc2jp1, wncjp1, Vj_nc2, PFj_nc2, model, j)
    (; beta, gamma,a_grid_nocollege, tasks_idx_nc2, shock_resources_nc, z_grid, Race, marital_status, fam_size, fam_type, working_years) = model

    fill!(Vj_nc2, 0f0)
    fill!(PFj_nc2, 0f0)

    # When creating interpolation objects for Vj+1:
    
    if j < working_years
        wnc_itp = nothing
        vnc2_itp = [LinearInterpolation((a_grid_nocollege, z_grid[R]), vnc2jp1[R, m, n, t, :, :, shock_in, shock_out, past_in, past_out], extrapolation_bc=Flat()) 
           for R in Race, m in marital_status, n in fam_size, t in fam_type, shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2]
    else
        wnc_itp = [LinearInterpolation((a_grid_nocollege, z_grid[R]), wncjp1[R, m, n, t, :, :, shock_in, shock_out, past_in, past_out], extrapolation_bc=Flat()) 
               for R in Race, m in marital_status, n in 1:2, t in fam_type, shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2]
        vnc2_itp = nothing
    end

    @threads for idx in eachindex(tasks_idx_nc2)
        (R, m, n, t, i_a, i_z) = tasks_idx_nc2[idx]

        if j < working_years
            wncjp1_itp = nothing
            vnc2jp1_itp = vnc2_itp[R, m, n, t, :, :, :, :]

        else    
            wncjp1_itp = wnc_itp[R, m, m, t, :, :, :, :]
            vnc2jp1_itp = nothing
        end

        for shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2

            net_resources = shock_resources_nc[j, R, m, n,t, i_a, i_z, shock_in, shock_out, past_in, past_out]

            if j < working_years
                result = optimize(ap1 -> -(u(net_resources - ap1, gamma) + 
                     beta * EVnc2_jp1(model, vnc2jp1_itp, j, R, m, n, t, ap1, i_z, shock_in, shock_out, past_in, past_out)),
                     0.0, net_resources,Brent(); rel_tol=1e-4, abs_tol=1e-4)
            else
                result = optimize(ap1 -> -(u(net_resources - ap1, gamma) + 
                     beta * EWnc_jp1(model, wncjp1_itp, j, R, m, m, t, ap1, i_z, shock_in, shock_out, past_in, past_out)),
                     0.0, net_resources, Brent(); rel_tol=1e-4, abs_tol=1e-4)
            end
        
            Vj_nc2[R, m, n, t, i_a, i_z, shock_in, shock_out, past_in, past_out] = -result.minimum
            PFj_nc2[R, m, n, t, i_a, i_z, shock_in, shock_out, past_in, past_out] = result.minimizer
        end
    end
    
    return Vj_nc2, PFj_nc2
end

# Expected family with no family transition
function EVnc2_jp1(model, vjp1, j, R, m, n, t,ap1, i_z, shock_in, shock_out, past_in, past_out)
    (; Pimat, zpnts, y_values, ra_w, ra_b) = model
    jp1 = j + 1
    a_income = R == 1 ? ap1 * ra_w : ap1 * ra_b

    expected_value = 0.0
    for i_zp1 in 1:zpnts
        pi_z =Pimat[R][i_z, i_zp1]
        y = y_values[R, jp1, m, 1, i_zp1]
        for shock_in_next in 1:2, shock_out_next in 1:2, past_in_next in 1:2, past_out_next in 1:2
            # update past in and past out flags based on past flags and past shocks
            if past_in == 2 && past_in_next == 1 || shock_in == 2 && past_in_next == 1 
                continue  
            end
            if past_out == 2 && past_out_next == 1 || shock_out == 2 && past_out_next == 1 
                continue 
            end

            # probability of shock in and shock out next period
            shock_out_prob= shocks_out_prob(R,n,m,jp1,y, a_income, 1, t, past_in-1, past_out-1)
            shock_in_prob = shocks_in_prob(R,n,m,jp1,y, a_income, 1, t, past_in-1, past_out-1)
            if shock_in_next == 2
                shock_in_next_prob = shock_in_prob
            else
                shock_in_next_prob = 1 - shock_in_prob
            end
            if shock_out_next == 2
                shock_out_next_prob = shock_out_prob
            else
                shock_out_next_prob = 1 - shock_out_prob
            end

            # expected value contribution from next period state
            expected_value += pi_z * shock_out_next_prob *  shock_in_next_prob* vjp1[shock_in_next, shock_out_next, past_in_next, past_out_next](ap1, i_zp1)
        end
    end

    return expected_value
end



###############################################################################################
# Retirement
###############################################################################################

function Wncj(wncjp1, Wj_nc, WPFj_nc, model, j)
    (; survival_risk, beta, gamma , shock_resources_nc, a_grid_nocollege, z_grid, tasks_idx_nc2, jpnts, Race, marital_status, fam_type) = model

    fill!(Wj_nc, 0f0)
    fill!(WPFj_nc, 0f0)

    # Create interpolation object
    if j < jpnts
        wnc_itp = [LinearInterpolation((a_grid_nocollege, z_grid[R]), wncjp1[R, m, n, t,  :, :, shock_in, shock_out, past_in, past_out], extrapolation_bc=Interpolations.Flat())
                for R in Race, m in marital_status, n in 1:2, t in fam_type, shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2]
    else
        wnc_itp = nothing
    end

    @threads for idx in eachindex(tasks_idx_nc2)
        (R, m, n, t, i_a, i_z) = tasks_idx_nc2[idx]
        if n > 2 || (m== 2 && n == 1)
            continue  # Skip family sizes that don't exist in retirement phase
        end

        sj = survival_risk[j - 42, R]
        
        if j < jpnts
            wncjp1_itp = wnc_itp[R, m, m, t, :, :, :, :]
        else
            wncjp1_itp = nothing
        end

        for shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2

            net_resources = shock_resources_nc[j, R, m, m, t, i_a, i_z, shock_in, shock_out, past_in, past_out]

            if j == jpnts
                # Last period of life, no future value
                Wj_nc[R, m, m, t, i_a, i_z, shock_in, shock_out, past_in, past_out] = u(net_resources, gamma)
                WPFj_nc[R, m, m, t, i_a, i_z, shock_in, shock_out, past_in, past_out] = 0.0
            else
                result = optimize(ap1 -> - (u(net_resources - ap1, gamma) + beta *  sj* EWnc_jp1(model, wncjp1_itp, j, R, m, m,  t, ap1, i_z, shock_in, shock_out, past_in, past_out)),
                            0.0, net_resources, Brent(); rel_tol=1e-4, abs_tol=1e-4)
                Wj_nc[R, m, m, t, i_a, i_z, shock_in, shock_out, past_in, past_out] = -result.minimum
                WPFj_nc[R, m, m, t, i_a, i_z, shock_in, shock_out,past_in,past_out] = result.minimizer
            end
        end
    end

    return Wj_nc, WPFj_nc
end

function EWnc_jp1(model, wncjp1_itp, j, R, m, n, t, ap1, i_z, shock_in, shock_out, past_in, past_out)

    (; y_values, ra_w, ra_b) = model
    jp1 = j + 1
    a_income = R == 1 ? ap1 * ra_w : ap1 * ra_b

    expected_value = 0.0
    y = y_values[R, j, m, 1, i_z]

    for shock_in_next in 1:2, shock_out_next in 1:2, past_in_next in 1:2, past_out_next in 1:2
        if past_in == 2 && past_in_next == 1 || shock_in == 2 && past_in_next == 1 || past_out == 2 && past_out_next == 1 || shock_out == 2 && past_out_next == 1
            continue  
        end

        # probability of shock in and shock out next period
        shock_out_prob= shocks_out_prob(R,n,m,jp1,y, a_income, 1, t, past_in-1, past_out-1  )
        shock_in_prob = shocks_in_prob(R,n,m,jp1,y, a_income, 1, t, past_in-1, past_out-1)
        
        if shock_in_next == 2
            shock_in_next_prob = shock_in_prob
        else
            shock_in_next_prob = 1 - shock_in_prob
        end

        if shock_out_next == 2
            shock_out_next_prob = shock_out_prob    
        else
            shock_out_next_prob = 1 - shock_out_prob
        end

        # Expected value contribution from next period state
        expected_value += shock_out_next_prob * shock_in_next_prob* wncjp1_itp[shock_in_next, shock_out_next, past_in_next, past_out_next](ap1, i_z)
        
    end

    return expected_value
end
