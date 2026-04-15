# Value function for those that chose not to attend college
# Positive assets only
# e = 0 no college

###############################################################################################
# Working
###############################################################################################

function Vncj_solve(vncjp1, wncjp1, Vj_nc, PFj_nc, model, j)
    (; beta, gamma,a_grid_nocollege, tasks_idx_nc, shock_resources_nc, fam_shock_period) = model

    fill!(Vj_nc, 0f0)
    fill!(PFj_nc, 0f0)

    # When creating interpolation objects for Vj+1:
    vnc_itp = [LinearInterpolation((a_grid_nocollege, z_grid[R]), vncjp1[R, m, n, t, :, :, shock_in, shock_out, past_in, past_out], extrapolation_bc=Flat()) 
           for R in Race, m in marital_status, n in fam_size, t in fam_type, shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2]

    if j < working_years
        wnc_itp = nothing
    else
        wnc_itp = LinearInterpolation((a_grid_nocollege, z_grid[R]), wncjp1[R, m, n, t, :, :, :, :], extrapolation_bc=Flat())
    end

    @threads for idx in eachindex(tasks_idx_nc)
        (R, m, n, t, i_a, i_z) = tasks_idx_nc[idx]

        vncjp1_itp = vnc_itp[R, m, n, t, :, :, :, :]
        if j < working_years
            wncjp1_itp = nothing
        else    
            wncjp1_itp = wnc_itp[R, m, m, t, :, :, :, :]
        end

        for shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2

            net_resources = shock_resources_nc[j, R, m, n,t, 1, i_a, i_z, shock_in, shock_out, past_in, past_out]

            if j == fam_shock_period
                result = optimize(ap1 -> -(u(net_resources - ap1, gamma) + 
                     beta * EVnc_family_jp1(model, vnc_itp, j, R, ap1, i_z, shock_in, shock_out, past_in, past_out)),
                     0.0, net_resources,Brent(); rel_tol=1e-4, abs_tol=1e-4)

            elseif j < working_years
                result = optimize(ap1 -> -(u(net_resources - ap1, gamma) + 
                     beta * EVnc_jp1(model, vncjp1_itp, j, R, m, n, ap1, i_z, shock_in, shock_out, past_in, past_out)),
                     0.0, net_resources,Brent(); rel_tol=1e-4, abs_tol=1e-4)
            else
                result = optimize(ap1 -> -(u(net_resources - ap1, gamma) + 
                     beta * EWnc_jp1(model, wncjp1_itp, j, R, m, n, t, ap1, i_z, shock_in, shock_out, past_in, past_out)),
                     0.0, net_resources, Brent(); rel_tol=1e-4, abs_tol=1e-4)
            end
        
            Vj_nc[R, m, n, t, i_a, i_z, shock_in, shock_out, past_in, past_out] = -result.minimum
            PFj_nc[R, m, n, t, i_a, i_z, shock_in, shock_out, past_in, past_out] = result.minimizer
        end
    end
    
    return Vj_nc, PFj_nc
end

# Expected family with no family transition
function EVnc_jp1(model, vjp1, j, R, m, n, ap1, i_z, shock_in, shock_out, past_in, past_out)
    (; Pimat, zpnts, y_values) = model
    jp1 = j + 1
    a_income = R == 1 ? ap1 * ra_w : ap1 * ra_b

    expected_value = 0.0
    for i_zp1 in 1:zpnts
        pi_z =Pimat[i_z, i_zp1]
        y = y_values[R, j, m, 1, i_z]
        for shock_in_next in 1:2, shock_out_next in 1:2, past_in_next in 1:2, past_out_next in 1:2
            # update past in and past out flags based on past flags and past shocks
            if past_in == 2 && past_in_next == 1 || shock_in == 2 && past_in_next == 1 
                continue  
            end
            if past_out == 2 && past_out_next == 1 || shock_out == 2 && past_out_next == 1 
                continue 
            end

            # probability of shock in and shock out next period
            shock_out_prob= shocks_out_prob(R,n,m,jp1,y, a_income, 1, t, past_in, past_out)
            shock_in_prob = shocks_in_prob(R,n,m,jp1,y, a_income, 1, t, past_in, past_out)
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
function EVnc_family_jp1(model, vnc_itp, j, R,  ap1, i_z, shock_in, shock_out, past_in, past_out)
    (; Pimat, zpnts, y_values) = model
    jp1 = j + 1
    a_income = R == 1 ? ap1 * ra_w : ap1 * ra_b

    outcomes = family_shock_probs[(R, e_college(1))]

    expected_value = 0.0
    for (m_next, n_next, t_next, prob_fam) in outcomes

        for i_zp1 in 1:zpnts
            pi_z =Pimat[i_z, i_zp1]
            y = y_values[R, jp1, m_next, 1, i_z]
            for shock_in_next in 1:2, shock_out_next in 1:2, past_in_next in 1:2, past_out_next in 1:2
                # update past in and past out flags based on past flags and past shocks
                if past_in == 2 && past_in_next == 1 || shock_in == 2 && past_in_next == 1 || past_out == 2 && past_out_next == 1 || shock_out == 2 && past_out_next == 1
                    continue  
                end

                # probability of shock in and shock out next period
                shock_out_prob= shocks_out_prob(R,n_next,m_next,jp1,y, a_income, 1, t_next, past_in, past_out)
                shock_in_prob = shocks_in_prob(R,n_next,m_next,jp1,y, a_income, 1, t_next, past_in, past_out)
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
                expected_value += prob_fam * pi_z * shock_out_next_prob *  shock_in_next_prob* vnc_itp[R, m_next, n_next, t_next, shock_in_next, shock_out_next, past_in_next, past_out_next](ap1, i_zp1)
       
            end
        end
    end

    return expected_value
end


###############################################################################################
# Retirement
###############################################################################################

function Wncj(wncjp1, Wj_nc, WPFj_nc, model, j)
    (; survival_risk, beta, gamma , shock_resources_nc, a_grid_nocollege, tasks_idx_nc, jpnts) = model

    fill!(Wj_nc, 0f0)
    fill!(WPFj_nc, 0f0)

    # Create interpolation object
    wnc_itp = [LinearInterpolation((a_grid_nocollege, z_grid[R]), wncjp1[R, m, n, t,  :, :, shock_in, shock_out, past_in, past_out], extrapolation_bc=Interpolations.Flat())
                for R in Race, m in marital_status, n in fam_size, t in fam_type, shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2]

    @threads for idx in eachindex(tasks_idx_nc)
        (R, m, n, t, i_a, i_z) = tasks_idx_nc[idx]
        sj = survival_risk[j, R]
        wncjp1_itp = wnc_itp[R, m, n, t, :, :, :, :, :]

        for shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2

            net_resources = shock_resources_nc[j, R, m, n, 1, i_a, i_z, shock_in, shock_out, past_in, past_out]

            if j == jpnts
                # Last period of life, no future value
                Wj_nc[R, m, n, t, i_a, i_z, shock_in, shock_out, past_in, past_out] = u(net_resources, gamma)
                WPFj_nc[R, m, n, t, i_a, i_z, shock_in, shock_out, past_in, past_out] = 0.0
            else
                result = optimize(ap1 -> - (u(net_resources - ap1, gamma) + beta *  sj* EWnc_jp1(model, wncjp1_itp, j, R, m, n,  t, ap1, i_z, shock_in, shock_out, past_in, past_out)),
                            0.0, net_resources,  Brent(); rel_tol=1e-4, abs_tol=1e-4)
                Wj_nc[R, m, n, t, i_a, i_z, shock_in, shock_out, past_in, past_out] = -result.minimum
                WPFj_nc[R, m, n, t, i_a, i_z, shock_in, shock_out,past_in,past_out] = result.minimizer
            end
        end
    end

    return Wj_nc, WPFj_nc
end

function EWnc_jp1(model, wncjp1_itp, j, R, m, n, t, ap1, i_z, shock_in, shock_out, past_in, past_out)

    (; y_values) = model
    jp1 = j + 1
    a_income = R == 1 ? ap1 * ra_w : ap1 * ra_b

    expected_value = 0.0
    y = y_values[R, j, m, 1, i_z]

    for shock_in_next in 1:2, shock_out_next in 1:2, past_in_next in 1:2, past_out_next in 1:2
        if past_in == 2 && past_in_next == 1 || shock_in == 2 && past_in_next == 1 || past_out == 2 && past_out_next == 1 || shock_out == 2 && past_out_next == 1
            continue  
        end

        # probability of shock in and shock out next period
        shock_out_prob= shocks_out_prob(R,n,m,jp1,y, a_income, 1, t, past_in, past_out)
        shock_in_prob = shocks_in_prob(R,n,m,jp1,y, a_income, 1, t, past_in, past_out)
        
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
