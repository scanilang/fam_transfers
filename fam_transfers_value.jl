###############################################################################################
# School
###############################################################################################

function VSj(vjp1, model, j)
    (; beta, gamma, r, r_loan, d_limit, school_a_grid, z_grid,
       tuition_2yr, tuition_4yr, Race, fam_type, ed_type) = model

    @threads for idx in eachindex(school_tasks_idx)
        (R, t, e, i_a, i_z) = school_tasks_idx[idx]
        
        # e=2 is 2yr, e=3 is 4yr (e=1 is no college, shouldn't be here)
        tuition = e == 2 ? tuition_2yr : tuition_4yr
        
        a = school_a_grid[i_a]
        
        # Resources: asset return + parental transfer - tuition
        parental_transfer = expected_edu_transfer(R, n, m, j, y, a_income, e, t, e - 1)
        # e-1 maps ed_type to degree_choice: e=2→1(2yr), e=3→2(4yr)
        
        if a >= 0
            resources = a * (1 + r) + parental_transfer - tuition
        else
            resources = a * (1 + r_loan) + parental_transfer - tuition
        end
        
        # Natural borrowing limit at this age
        borrow_floor = -d_limit[j, R, e]
        
        # Ensure feasible optimization bounds
        ub = resources - 0.001
        lb = max(borrow_floor, ub - 1e6)  # lb must be < ub
        
        if ub <= lb
            # Can't afford anything — assign very low value
            VSj[R, t, e, i_a, i_z] = -1e10
            SPFj[R, t, e, i_a, i_z] = lb
        else
            result = optimize(
                ap1 -> -(u(max(resources - ap1, 0.001), gamma) + 
                         beta * EVS_jp1(model, vjp1, R, t, e, ap1, i_z)),
                lb, ub, Brent(); rel_tol=1e-4, abs_tol=1e-4)
            
            VSj[R, t, e, i_a, i_z] = -result.minimum
            SPFj[R, t, e, i_a, i_z] = result.minimizer
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
    v_itp = [LinearInterpolation((a_grid, z_grid[R]), vjp1[R, m, n, t, e, :, :, shock_in, shock_out, past_in, past_out], extrapolation_bc=Interpolations.Flat()) 
                                for R in Race, m in marital_status, n in fam_size, t in fam_type, e in ed_type, shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2];

    @threads for idx in eachindex(tasks_idx)
        (R, m, n, e, i_a, i_z) = tasks_idx[idx]
        vjp1_itp = v_itp[R, m, n, t, e, :, :, :, :, :]
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
    (; Pimat, zpnts) = model
    jp1 = j + 1
    a_income = R == 1 ? ap1 * ra_w : ap1 * ra_b

    expected_value = 0.0
    for i_zp1 in 1:zpnts
        pi_z =Pimat[i_z, i_zp1]
        y = g(R, j, r, n) * z_grid[R][i_z]
        for shock_in_next in 1:2, shock_out_next in 1:2, past_in_next in 1:2, past_out_next in 1:2
            if past_in == 2 && past_in_next == 1
                continue  
            end
            if past_out == 2 && past_out_next == 1
                continue 
            end
            shock_out= shocks_out_prob(β_white_probit_out, β_black_probit_out, n, m, jp1, y_next, a_income, e, t, past_in_next, past_out) 
            shock_in = shocks_in_prob(β_white_probit_in, β_black_probit_in, n, m, j, y_next, a_income, e, t, past_in_next, past_out)
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
            expected_value += pi_z * shock_out_next_val *  shock_out_next_val* vjp1[shock_in, shock_out, past_in_next, past_out_next](ap1, i_zp1)
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
        w_itp = [LinearInterpolation((a_grid, z_grid[R]), wjp1[R, m, n, t, e, :, :, shock_in, shock_out, past_in, past_out], extrapolation_bc=Interpolations.Flat())
                for R in Race, m in marital_status, n in fam_size, t in fam_type, e in ed_type, shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2]
    else
        w_itp = nothing
    end 

    @threads for idx in eachindex(tasks_idx)
        (R, m, n, e, i_a, i_z) = tasks_idx[idx]
        sj = survival_risk[j, R]
        wjp1_itp = w_itp[R, m, n, t, e, :, :, :, :, :]
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
    (; survival_risk ) = model
    jp1 = j + 1
    a_income = R == 1 ? ap1 * ra_w : ap1 * ra_b

    expected_value = 0.0
    y_next = g(R, j, e, n) * z_grid[R][i_z]
    for shock_in_next in 1:2, shock_out_next in 1:2, past_in_next in 1:2, past_out_next in 1:2
            if past_in == 2 && past_in_next == 1
                continue  
            end
            if past_out == 2 && past_out_next == 1
                continue 
            end
            shock_out= shocks_out_prob(β_white_probit_out, β_black_probit_out, n, m, jp1, y_next, a_income, e, t, past_in_next, past_out) 
            shock_in = shocks_in_prob(β_white_probit_in, β_black_probit_in, n, m, j, y_next, a_income, e, t, past_in_next, past_out)
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
            expected_value += shock_out_next_val * shock_out_next_val* wjp1_itp[shock_in, shock_out, past_in_next, past_out_next](ap1, i_z)
        end

    return expected_value
end
