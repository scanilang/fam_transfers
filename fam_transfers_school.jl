###############################################################################################
# School
###############################################################################################

function VSj_first_period(vsjp1, model)
    (; beta, gamma, r, r_loan, school_a_grid, d_limit,
       tuition_2yr, tuition_4yr) = model

    vsjp1_itp = [LinearInterpolation((a_grid), vsjp1[R, m, n, t, e, :, shock_in, shock_out, past_in, past_out], extrapolation_bc=Interpolations.Flat()) 
                                for R in Race, m in marital_status, n in fam_size, t in fam_type, e in ed_type, shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2]

    @threads for idx in eachindex(school_tasks_idx)
        (R, t, e, i_a) = school_tasks_idx[idx]
        
        a = school_a_grid[i_a]
        tuition = e == 2 ? tuition_2yr : tuition_4yr
        
        # Total expected parental transfer (lump sum) (based on parents characteristics and student's degree choice)
        prob_help = edu_transfer_prob(R, n, m, e, y, a_income, e, t, degree_choice)
        amt_help  = edu_transfer_amount(R, n, m, e, y, a_income, e, t, degree_choice)
        
        for shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2, prob_help in 1:2

            edu_transfer = prob_help == 2 ? amt_help : 0.0
            intervivos_transfers = net_transfers[1, R, m, n, e, i_a, 1, shock_in, shock_out, past_in, past_out] # Assuming z_idx = 1 for first period?
            
            if a >= 0
                resources = a * (1 + r) + edu_transfer - tuition/e + intervivos_transfers
            else
                resources = a * (1 + r_loan) + edu_transfer - tuition/e + intervivos_transfers
            end
            
            borrow_floor = -d_limit[1, R, e]
            ub = resources - 0.001
            lb = max(borrow_floor, ub - 1e6)
            
            if ub <= lb
                VSj_arr[R, t, e, i_a, shock_in] = -1e10
                SPFj_arr[R, t, e, i_a, shock_in] = lb
            else
                result = optimize(
                    ap1 -> -(u(max(resources - ap1, 0.001), gamma) +
                             beta * EVS_jp1((model, vsjp1_itp, R, t, e, ap1, i_z))),
                    lb, ub, Brent(); rel_tol=1e-4, abs_tol=1e-4)
                
                VSj_arr[R, t, e, i_a, shock_in] = -result.minimum
                SPFj_arr[R, t, e, i_a, shock_in] = result.minimizer
            end
        end
    end
    return VSj_arr, SPFj_arr
end

function VSj_enrolled(vsjp1_itp, vjp1, model, j)
    (; beta, gamma, r, r_loan, a_grid_college, d_limit,
       tuition_2yr, tuition_4yr) = model

    @threads for idx in eachindex(school_tasks_idx)
        (R, t, e, i_a) = school_tasks_idx[idx]
        
        a = a_grid_college[i_a]
        tuition = e == 2 ? tuition_2yr : tuition_4yr
        
        if a >= 0
            resources = a * (1 + r) - tuition
        else
            resources = a * (1 + r_loan) - tuition
        end
        
        borrow_floor = -d_limit[j, R, e]
        ub = resources - 0.001
        lb = max(borrow_floor, ub - 1e6)
        
        if ub <= lb
            VSj_arr[R, t, e, i_a] = -1e10
            SPFj_arr[R, t, e, i_a] = lb
        else
            if j == 2 & e == 2
               result = optimize(
                    ap1 -> -(u(max(resources - ap1, 0.001), gamma) +
                         beta * EV_jp1(model, vjp1, R, t, e, ap1)),
                    lb, ub, Brent(); rel_tol=1e-4, abs_tol=1e-4)
            else 
                result = optimize(
                    ap1 -> -(u(max(resources - ap1, 0.001), gamma) +
                         beta * EVS_jp1(model, vsjp1_itp, R, t, e, ap1)),
                    lb, ub, Brent(); rel_tol=1e-4, abs_tol=1e-4)
            end
            
            VSj_arr[R, t, e, i_a] = -result.minimum
            SPFj_arr[R, t, e, i_a] = result.minimizer
        end
    end
    return VSj_arr, SPFj_arr
end


function EVS_jp1(model, vsjp1_itp, R, t, e, ap1)
    (; prob_shocks ) = model

    expected_value = 0.0
    for shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2
        expected_value += prob_shocks[jp1, R, m, n, e, i_z, shock_in, shock_out, past_in, past_out] * vsjp1_itp[R, m, n, t, e, :, :, shock_in, shock_out, past_in, past_out](ap1)
    end

    return expected_value
end
# Value function for those that chose not to attend college
# Positive assets only
# e = 0 no college

###############################################################################################
# Working
###############################################################################################

function Vcj_solve(vcjp1, wcjp1, Vj_c, PFj_c, model, j)
    (; beta, gamma,a_grid_college, tasks_idx_c) = model

    fill!(Vj_c, 0f0)
    fill!(PFj_c, 0f0)

    # When creating interpolation objects for Vj+1:
    vc_itp = LinearInterpolation((a_grid_college, z_grid[R]), vcjp1[R, m, n, t, :, :, :, :], extrapolation_bc=Flat())

    if j < working_years
        wc_itp = nothing
    else
        wc_itp = LinearInterpolation((a_grid_college, z_grid[R]), wcjp1[R, m, n, t, :, :, :, :], extrapolation_bc=Flat())
    end

    # Lower bound: no college can't borrow
    lb = 0.0

    @threads for idx in eachindex(tasks_idx_c)
        (R, m, n, i_a, i_z) = tasks_idx_c[idx]

        vcjp1_itp = vc_itp[R, m, n, t, :, :, :, :]
        if j < working_years
            wcjp1_itp = nothing
        else    
            wcjp1_itp = wc_itp[R, m, m, t, :, :, :, :]
        end

        for shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2

            net_resources = shock_resources[j, R, m, n,t, 1, i_a, i_z, shock_in, shock_out, past_in, past_out]

            if j == fam_shock_age
                result = optimize(ap1 -> -(u(net_resources - ap1, gamma) + 
                     beta * EVc_family_jp1(model, vc_itp, j, R, m, n, ap1, i_z, shock_in, shock_out, past_in, past_out)),
                     lb, net_resources,Brent(); rel_tol=1e-4, abs_tol=1e-4)

            elseif j < working_years
                result = optimize(ap1 -> -(u(net_resources - ap1, gamma) + 
                     beta * EVc_jp1(model, vcjp1_itp, j, R, m, n, ap1, i_z, shock_in, shock_out, past_in, past_out)),
                     lb, net_resources,Brent(); rel_tol=1e-4, abs_tol=1e-4)
            else
                result = optimize(ap1 -> -(u(net_resources - ap1, gamma) + 
                     beta * EWc_jp1(model, wcjp1_itp, j, R, m, m, ap1, i_z, shock_in, shock_out, past_in, past_out)),
                     lb, net_resources, Brent(); rel_tol=1e-4, abs_tol=1e-4)
            end
        
            Vj_c[R, m, n, t, i_a, i_z, shock_in, shock_out, past_in, past_out] = -result.minimum
            PFj_c[R, m, n, t, i_a, i_z, shock_in, shock_out, past_in, past_out] = result.minimizer
        end
    end
    
    return Vj_c, PFj_c
end

# Expected family with no family transition
function EVc_jp1(model, vjp1, j, R, m, n, ap1, i_z, shock_in, shock_out, past_in, past_out)
    (; Pimat, zpnts, y_values) = model
    jp1 = j + 1
    a_income = R == 1 ? ap1 * ra_w : ap1 * ra_b

    expected_value = 0.0
    for i_zp1 in 1:zpnts
        pi_z =Pimat[i_z, i_zp1]
        y = y_values[R, j, m, t, e, i_z]
        for shock_in_next in 1:2, shock_out_next in 1:2, past_in_next in 1:2, past_out_next in 1:2
            # update past in and past out flags based on past flags and past shocks
            if past_in == 2 && past_in_next == 1 || shock_in == 2 && past_in_next == 1 
                continue  
            end
            if past_out == 2 && past_out_next == 1 || shock_out == 2 && past_out_next == 1 
                continue 
            end

            # probability of shock in and shock out next period
            shock_out= shocks_out_prob(r,n,m,jp1,y, a_income, e, t, past_in, past_out)
            shock_in = shocks_in_prob(r,n,m,jp1,y, a_income, e, t, past_in, past_out)
            if shock_in_next == 2
                shock_in_next_val = shock_in
            else
                shock_in_next_val = 1 - shock_in
            end
            if shock_out_next == 2
                shock_out_next_val = shock_out
            else
                shock_out_next_val = 1 - shock_out
            end

            # expected value contribution from next period state
            expected_value += pi_z * shock_out_next_val *  shock_out_next_val* vjp1[shock_in, shock_out, past_in_next, past_out_next](ap1, i_zp1)
        end
    end

    return expected_value
end

# Expected value with family transition
function EVc_family_jp1(model, vc_itp, j, R, m, n, ap1, i_z, shock_in, shock_out, past_in, past_out)
    (; Pimat, zpnts, y_values) = model
    jp1 = j + 1
    a_income = R == 1 ? ap1 * ra_w : ap1 * ra_b

    expected_value = 0.0
    for n_next in 1:6, t_next in fam_type, m_next in marital_status
        if n_next == 1 && m_next == 2 || n_next > 1 && m_next == 1
            continue
        end

        prob_family_transition = family_transition_prob(n, t, m, n_next, t_next, m_next)

        for i_zp1 in 1:zpnts
            pi_z =Pimat[i_z, i_zp1]
            y = y_values[R, j, m, t, e, i_z]
            for shock_in_next in 1:2, shock_out_next in 1:2, past_in_next in 1:2, past_out_next in 1:2
                # update past in and past out flags based on past flags and past shocks
                if past_in == 2 && past_in_next == 1 || shock_in == 2 && past_in_next == 1 || past_out == 2 && past_out_next == 1 || shock_out == 2 && past_out_next == 1
                    continue  
                end

                # probability of shock in and shock out next period
                shock_out= shocks_out_prob(r,n_next,m_next,jp1,y, a_income, e, t_next, past_in, past_out)
                shock_in = shocks_in_prob(r,n_next,m_next,jp1,y, a_income, e, t_next, past_in, past_out)
                if shock_in_next == 2
                    shock_in_next_prob = shock_in
                else
                    shock_in_next_prob = 1 - shock_in
                end
                if shock_out_next == 2
                    shock_out_next_prob = shock_out
                else
                    shock_out_next_prob = 1 - shock_out
                end

                # expected value contribution from next period state
                expected_value += prob_family_transition * pi_z * shock_out_next_prob *  shock_in_next_prob* vc_itp[R, m_next, n_next, t_next, shock_in, shock_out, past_in_next, past_out_next](ap1, i_zp1)
       
            end
        end
    end

    return expected_value
end


###############################################################################################
# Retirement
###############################################################################################

function Wcj(wcjp1, Wj_c, WPFj_c, model, j)
    (; survival_risk, beta, gamma , shock_resources, a_grid_college, tasks_idx_c, jpnts) = model

    fill!(Wj_c, 0f0)
    fill!(WPFj_c, 0f0)

    # Create interpolation object
    wc_itp = [LinearInterpolation((a_grid_college, z_grid[R]), wcjp1[R, m, n, t,  :, :, shock_in, shock_out, past_in, past_out], extrapolation_bc=Interpolations.Flat())
                for R in Race, m in marital_status, n in fam_size, t in fam_type, shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2]

    @threads for idx in eachindex(tasks_idx_c)
        (R, m, n, i_a, i_z) = tasks_idx_c[idx]
        sj = survival_risk[j, R]
        wcjp1_itp = wc_itp[R, m, n, t, :, :, :, :, :]

        for shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2

            net_resources = shock_resources[j, R, m, n, 1, i_a, i_z, shock_in, shock_out, past_in, past_out]

            if j == jpnts
                # Last period of life, no future value
                Wj_c[R, m, n, t, i_a, i_z, shock_in, shock_out, past_in, past_out] = u(net_resources, gamma)
                WPFj_c[R, m, n, t, i_a, i_z, shock_in, shock_out, past_in, past_out] = 0.0
            else
                result = optimize(ap1 -> - (u(net_resources - ap1, gamma) + beta *  sj* EWc_jp1(model, wcjp1_itp, j, R, m, n, ap1, i_z, shock_in, shock_out, past_in, past_out)),
                            0.0, net_resources,  Brent(); rel_tol=1e-4, abs_tol=1e-4)
                Wj_c[R, m, n, t, i_a, i_z, shock_in, shock_out, past_in, past_out] = -result.minimum
                WPFj_c[R, m, n, t, i_a, i_z, shock_in, shock_out,past_in,past_out] = result.minimizer
            end
        end
    end

    return Wj_c, WPFj_c
end

function EWc_jp1(model, wcjp1_itp, j, R, m, n, ap1, i_z, shock_in, shock_out, past_in, past_out)

    (; y_values) = model
    jp1 = j + 1
    a_income = R == 1 ? ap1 * ra_w : ap1 * ra_b

    expected_value = 0.0
    y = y_values[R, j, m, t, e, i_z]

    for shock_in_next in 1:2, shock_out_next in 1:2, past_in_next in 1:2, past_out_next in 1:2
        if past_in == 2 && past_in_next == 1 || shock_in == 2 && past_in_next == 1 || past_out == 2 && past_out_next == 1 || shock_out == 2 && past_out_next == 1
            continue  
        end

        # probability of shock in and shock out next period
        shock_out= shocks_out_prob(r,n,m,jp1,y, a_income, e, t, past_in, past_out)
        shock_in = shocks_in_prob(r,n,m,jp1,y, a_income, e, t, past_in, past_out)
        
        if shock_in_next == 2
            shock_in_next_prob = shock_in
        else
            shock_in_next_prob = 1 - shock_in
        end

        if shock_out_next == 2
            shock_out_next_prob = shock_out
        else
            shock_out_next_prob = 1 - shock_out
        end

        # Expected value contribution from next period state
        expected_value += shock_out_next_prob * shock_in_next_prob* wcjp1_itp[shock_in, shock_out, past_in_next, past_out_next](ap1, i_z)
        
    end

    return expected_value
end
