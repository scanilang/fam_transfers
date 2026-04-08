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
    (; beta, gamma, r, r_loan, school_a_grid, d_limit,
       tuition_2yr, tuition_4yr) = model

    @threads for idx in eachindex(school_tasks_idx)
        (R, t, e, i_a) = school_tasks_idx[idx]
        
        a = school_a_grid[i_a]
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

###############################################################################################
# Working
###############################################################################################

function Vj(vjp1, model, j)
    (; beta, gamma, r, r_loan, a_grid_nocollege, a_grid_college,
       shock_resources) = model

       # When creating interpolation objects for Vj+1:
    vjp1_itp = if e == 1
        LinearInterpolation((a_grid_nocollege, z_grid[R]), 
            vjp1[R, m, n, t, e, :, :], extrapolation_bc=Flat())
    else
        LinearInterpolation((a_grid_college, z_grid[R]), 
        vjp1[R, m, n, t, e, :, :], extrapolation_bc=Flat())
    end

    @threads for idx in eachindex(tasks_idx)
        (R, m, n, e, i_a, i_z) = tasks_idx[idx]
        
        # Select grid based on education
        a_grid_e = e == 1 ? a_grid_nocollege : a_grid_college
        a = a_grid_e[i_a]
        
        # Interest rate depends on sign of assets
        if a >= 0
            a_income = R == 1 ? a * ra_w : a * ra_b
        else
            a_income = a * (1 + r_loan)  # debt accrues at loan rate
        end
        
        y = g(R, j, e, m) * z_grid[R][i_z]
        resources = y + a_income - tax_y(y, m)
        
        # Lower bound: no college can't borrow, college can carry debt
        if e == 1
            lb = 0.0
        else
            lb = -d_limit[j, R, e]  # natural borrowing limit
        end
        
        result = optimize(
            ap1 -> -(u(max(resources - ap1, 0.001), gamma) + 
                     beta * EV_jp1(model, vjp1_itp, j, R, m, n, e, ap1, i_z)),
            lb, resources - 0.001,
            Brent(); rel_tol=1e-4, abs_tol=1e-4)
        
        Vj_arr[R, m, n, t, e, i_a, i_z] = -result.minimum
        PFj_arr[R, m, n, t, e, i_a, i_z] = result.minimizer
    end
    
    return Vj_arr, PFj_arr
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
                result = optimize(ap1 -> - (u(resources - ap1, gamma) + beta *  sj* EW_jp1(wjp1_itp, j, R, m, n, e, past_in, past_out,ap1, i_z)),
                            0.0, resources,  Brent(); rel_tol=1e-4, abs_tol=1e-4)
                Wj[R, m, n, e, i_a, i_z, shock_in, shock_out, past_in, past_out] = -result.minimum
                WPFj[R, m, n, e, i_a, i_z, shock_in, shock_out, past_in, past_out] = result.minimizer
            end
        end
    end

    return Wj, WPFj
end

function EW_jp1(wjp1_itp, j, R, m, n, e, past_in, past_out, ap1, i_z)

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
            expected_value += shock_out_next_val * shock_out_next_val* wjp1_itp[shock_in, shock_out, past_in_next, past_out_next](ap1, i_z)
        end

    return expected_value
end
