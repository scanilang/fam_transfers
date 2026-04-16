###############################################################################################
# School
###############################################################################################
#prob_help = edu_transfer_prob(R, n, m, e, y, a_income, e, t, degree_choice)

function VSj_first_period(vsjp1, Vsj_1, PFsj_1, model)
    (; beta, gamma, ra_w, ra_b, tasks_idx_c1, a_grid_nocollege, d_limit,
       tuition_2yr, tuition_4yr, y_values) = model

    fill!(Vsj_1, 0f0)
    fill!(PFsj_1, 0f0)

    vsjp1_itp = [LinearInterpolation((a_grid_nocollege), vsjp1[R, t, degree, :], extrapolation_bc=Interpolations.Flat()) 
                                for R in Race, t in fam_type, degree in 1:2]

    @threads for idx in eachindex(tasks_idx_c1)
       (R, m, n, t, e, i_a, i_z) = tasks_idx_c1[idx]
        
        # Skip negative assets in first period since they can't borrow yet    
        a = a_grid_nocollege[i_a]

        r = R == 1 ? ra_w : ra_b
        for edu_help in 1:2, degree in 1:2

            degree_choice  = [2, 4][degree]
            tuition = degree_choice == 2 ? tuition_2yr/2 : tuition_4yr/4
``
            if edu_help == 2
                # Parental transfer (lump sum) (based on parents characteristics and student's degree choice)
                a_income = a* r
                y = y_values[R, 43, m, e, i_z] # Parents age when child is 18
                edu_transfer  = edu_transfer_amount(R, y, a_income, e, degree_choice)
            else
                edu_transfer = 0.0
            end
            
            resources = a * (1 + r) + edu_transfer - tuition 
            
            borrow_floor = -d_limit[1, R, degree]
            ub = resources
            lb = max(borrow_floor, ub - 1e6)
            
            if ub <= lb
                Vsj_1[R, m, e, i_a, i_z, edu_help, degree] = -1e10
                PFsj_1[R, m, e, i_a, i_z, edu_help, degree] = lb
            else
                result = optimize(
                    ap1 -> -(u(max(resources - ap1, 0.001), gamma) +
                             beta * vsjp1_itp[R, t, degree](ap1)),
                    lb, ub, Brent(); rel_tol=1e-4, abs_tol=1e-4)
                
                Vsj_1[R, m, e, i_a, i_z, edu_help, degree] = -result.minimum
                PFsj_1[R, m, e, i_a, i_z, edu_help, degree] = result.minimizer
            end
        end
    end
    return Vsj_1, PFsj_1
end

function VSj_enrolled(vsjp1, vjp1, model, j)
    (; beta, gamma, r, r_loan, a_grid_college, d_limit,
       tuition_2yr, tuition_4yr, tasks_idx_c) = model

    fill!(Vsj, 0f0)
    fill!(PFsj, 0f0)

    vsjp1_itp = [LinearInterpolation((a_grid_college), vsjp1[R, t, degree, :], extrapolation_bc=Interpolations.Flat()) 
                                for R in Race, t in fam_type, degree in 1:2]

    if j < working_years
        vjp1_itp = nothing
    else
        vjp1_itp = [LinearInterpolation((a_grid_college, z_grid[R]), vjp1[R, m, n, t, e, :, :, shock_in, shock_out, past_in, past_out], extrapolation_bc=Flat()) 
               for R in Race, m in marital_status, n in fam_size, t in fam_type, e in 1:2, shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2]
    end

    @threads for idx in eachindex(tasks_idx_c)
       (R, t, e, i_a) = tasks_idx_c[idx]
       degree = e - 1
        
        a = a_grid_college[i_a]
        tuition = e == 2 ? tuition_2yr/e : tuition_4yr/e
        r = R == 1 ? ra_w : ra_b

        if a >= 0
            resources = a * (1 + r) - tuition
        else
            resources = a * (1 + r_loan) - tuition
        end

        borrow_floor = -d_limit[j, R, degree]
        ub = resources - 0.001
        lb = max(borrow_floor, ub - 1e6)
        
        if ub <= lb
            Vsj[R, t, degree, i_a] = -1e10
            PFsj[R, t, degree, i_a] = lb
        else
            if j == 2 && e == 1 || j == 4 && e == 2
               # Transition to working phase
                result = optimize(
                            ap1 -> -(u(max(resources - ap1, 0.001), gamma) +
                                     beta * EV_graduation(model, vjp1_itp, R, t, degree, ap1)),
                            lb, ub, Brent(); rel_tol=1e-4, abs_tol=1e-4)
            else 
                result = optimize(
                    ap1 -> -(u(max(resources - ap1, 0.001), gamma) +
                         beta * vsjp1_itp[R, t, degree](ap1)),
                    lb, ub, Brent(); rel_tol=1e-4, abs_tol=1e-4)
            end
            
            Vsj[R, t, degree, i_a] = -result.minimum
            PFsj[R, t, degree, i_a] = result.minimizer
        end
    end
    return Vsj, PFsj
end

# Value function for those that chose not to attend college
# Positive assets only
# e = 0 no college
function EV_graduation(model, vjp1_itp, R, t, degree, ap1)
    (; ra_w, ra_b, zpnts, z_grid, y_values,
       shocks_in_prob, shocks_out_prob) = model
    
    # At graduation: single (m=1), no kids (n=1), no transfer history
    m = 1
    n = 1
    past_in = 0
    past_out = 0
    
    a_income = ap1 >= 0 ? (R == 1 ? ap1 * ra_w : ap1 * ra_b) : 0.0
    
    # At graduation, which working j? 2yr → j=3 (age 20), 4yr → j=5 (age 22)
    j_grad_work = degree == 1 ? 3 : 5
    
    expected_value = 0.0
    
    # Integrate over initial z (uniform — no prior work history)
    for i_z in 1:zpnts
        pi_z = 1.0 / zpnts
        z_val = z_grid[R][i_z]
        y = y_values[R, j_grad_work, m, degree, i_z]
        
        # Integrate over first-period transfer shocks
        prob_in  = shocks_in_prob(R, n, m, j_grad_work, y, a_income, degree, t, past_in, past_out)
        prob_out = shocks_out_prob(R, n, m, j_grad_work, y, a_income, degree, t, past_in, past_out)
        
        for shock_in in 1:2, shock_out in 1:2
            p_in  = shock_in  == 2 ? prob_in  : 1 - prob_in
            p_out = shock_out == 2 ? prob_out : 1 - prob_out
            
            expected_value += pi_z * p_in * p_out *
                vjp1_itp[R, m, n, t, degree, shock_in, shock_out, past_in, past_out](ap1, z_val)
        end
    end
    
    return expected_value
end

###############################################################################################
# Working
###############################################################################################

function Vcj_solve(vcjp1, wcjp1, Vj_c, PFj_c, model, j)
    (; beta, gamma,a_grid_college, tasks_idx_c, shock_resources_c, d_limit, fam_shock_period) = model

    fill!(Vj_c, 0f0)
    fill!(PFj_c, 0f0)

    # When creating interpolation objects for Vj+1:
    vc_itp = [LinearInterpolation((a_grid_college, z_grid[R]), vcjp1[R, m, n, t, e, :, :, shock_in, shock_out, past_in, past_out], extrapolation_bc=Flat()) 
           for R in Race, m in marital_status, n in fam_size, t in fam_type, e in 1:2, shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2]

    if j < working_years
        wc_itp = nothing
    else
        wc_itp = [LinearInterpolation((a_grid_college, z_grid[R]), wcjp1[R, m, n, t, e, :, :, shock_in, shock_out, past_in, past_out], extrapolation_bc=Flat()) 
               for R in Race, m in marital_status, n in fam_size, t in fam_type, e in 1:2, shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2]
    end



    @threads for idx in eachindex(tasks_idx_c)
        (R, m, n, t, e, i_a, i_z) = tasks_idx_c[idx]

        vcjp1_itp = vc_itp[R, m, n, t, e-1, :, :, :, :]
        if j < working_years
            wcjp1_itp = nothing
        else    
            wcjp1_itp = wc_itp[R, m, n, t, e-1, :, :, :, :]
        end

        # Lower bound: natural borrowing limit
        lb = -d_limit[j, R, e-1]

        for shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2

            net_resources = shock_resources_c[j, R, m, n,t, e-1, i_a, i_z, shock_in, shock_out, past_in, past_out]

            if j == fam_shock_period
                result = optimize(ap1 -> -(u(net_resources - ap1, gamma) + 
                     beta * EVc_family_jp1(model, vc_itp, j, R, e, ap1, i_z, shock_in, shock_out, past_in, past_out)),
                     lb, net_resources,Brent(); rel_tol=1e-4, abs_tol=1e-4)

            elseif j < working_years
                result = optimize(ap1 -> -(u(net_resources - ap1, gamma) + 
                     beta * EVc_jp1(model, vcjp1_itp, j, R, m, n, e, t, ap1, i_z, shock_in, shock_out, past_in, past_out)),
                     lb, net_resources,Brent(); rel_tol=1e-4, abs_tol=1e-4)
            else
                result = optimize(ap1 -> -(u(net_resources - ap1, gamma) + 
                     beta * EWc_jp1(model, wcjp1_itp, j, R, m, m, e, t, ap1, i_z, shock_in, shock_out, past_in, past_out)),
                     lb, net_resources, Brent(); rel_tol=1e-4, abs_tol=1e-4)
            end
        
            Vj_c[R, m, n, t, e-1, i_a, i_z, shock_in, shock_out, past_in, past_out] = -result.minimum
            PFj_c[R, m, n, t, e-1, i_a, i_z, shock_in, shock_out, past_in, past_out] = result.minimizer
        end
    end
    
    return Vj_c, PFj_c
end

# Expected family with no family transition
function EVc_jp1(model, vjp1, j, R, m, n, t, e, ap1, i_z, shock_in, shock_out, past_in, past_out)
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
            shock_out_prob= shocks_out_prob(R,n,m,jp1,y, a_income, e, t, past_in, past_out)
            shock_in_prob = shocks_in_prob(R,n,m,jp1,y, a_income, e, t, past_in, past_out)
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
function EVc_family_jp1(model, vc_itp, j, R, e, ap1, i_z, shock_in, shock_out, past_in, past_out)
    (; Pimat, zpnts, y_values) = model
    jp1 = j + 1
    a_income = R == 1 ? ap1 * ra_w : ap1 * ra_b

    outcomes = family_shock_probs[(R, e_college(e))]
    e_idx = e == 2 ? 1 : 2
    expected_value = 0.0
    for (m_next, n_next, t_next, prob_fam) in outcomes

        for i_zp1 in 1:zpnts
            pi_z =Pimat[i_z, i_zp1]
            y = y_values[R, jp1, m_next, e_idx, i_zp1]
            for shock_in_next in 1:2, shock_out_next in 1:2, past_in_next in 1:2, past_out_next in 1:2
                # update past in and past out flags based on past flags and past shocks
                if past_in == 2 && past_in_next == 1 || shock_in == 2 && past_in_next == 1 || past_out == 2 && past_out_next == 1 || shock_out == 2 && past_out_next == 1
                    continue  
                end

                # probability of shock in and shock out next period
                shock_out_prob= shocks_out_prob(R,n_next,m_next,jp1,y, a_income, e, t_next, past_in, past_out)
                shock_in_prob = shocks_in_prob(R,n_next,m_next,jp1,y, a_income, e, t_next, past_in, past_out)
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
                expected_value += prob_fam * pi_z * shock_out_next_prob *  shock_in_next_prob* vc_itp[R, m_next, n_next, t_next, e, shock_in_next, shock_out_next, past_in_next, past_out_next](ap1, i_zp1)
    
            end
        end
    end

    return expected_value
end


###############################################################################################
# Retirement
###############################################################################################

function Wcj(wcjp1, Wj_c, WPFj_c, model, j)
    (; survival_risk, beta, gamma , shock_resources_c, a_grid_college, tasks_idx_c, jpnts, d_limit) = model

    fill!(Wj_c, 0f0)
    fill!(WPFj_c, 0f0)

    # Create interpolation object
    wc_itp = [LinearInterpolation((a_grid_college, z_grid[R]), wcjp1[R, m, n, t, e,   :, :, shock_in, shock_out, past_in, past_out], extrapolation_bc=Interpolations.Flat())
                for R in Race, m in marital_status, n in fam_size, t in fam_type, e in 1:2, shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2]


    @threads for idx in eachindex(tasks_idx_c)
        (R, m, n, t, e, i_a, i_z) = tasks_idx_c[idx]
        sj = survival_risk[j, R]
        wcjp1_itp = wc_itp[R, m, n, t,e-1, :, :, :, :]
        
        # Lower bound: natural borrowing limit
        lb = -d_limit[j, R, e-1]

        for shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2

            net_resources = shock_resources_c[j, R, m, n, e-1, i_a, i_z, shock_in, shock_out, past_in, past_out]

            if j == jpnts
                # Last period of life, no future value
                Wj_c[R, m, n, t, i_a, i_z, shock_in, shock_out, past_in, past_out] = u(net_resources, gamma)
                WPFj_c[R, m, n, t, i_a, i_z, shock_in, shock_out, past_in, past_out] = 0.0
            else
                result = optimize(ap1 -> - (u(net_resources - ap1, gamma) + beta *  sj* EWc_jp1(model, wcjp1_itp, j, R, m, n, e, t, ap1, i_z, shock_in, shock_out, past_in, past_out)),
                            lb, net_resources,  Brent(); rel_tol=1e-4, abs_tol=1e-4)
                Wj_c[R, m, n, t, i_a, i_z, shock_in, shock_out, past_in, past_out] = -result.minimum
                WPFj_c[R, m, n, t, i_a, i_z, shock_in, shock_out,past_in,past_out] = result.minimizer
            end
        end
    end

    return Wj_c, WPFj_c
end

function EWc_jp1(model, wcjp1_itp, j, R, m, n, t, e,ap1, i_z, shock_in, shock_out, past_in, past_out)

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
        shock_out_prob= shocks_out_prob(R,n,m,jp1,y, a_income, e, t, past_in, past_out)
        shock_in_prob = shocks_in_prob(R,n,m,jp1,y, a_income, e, t, past_in, past_out)
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
        expected_value += shock_out_next_prob *  shock_in_next_prob* wcjp1_itp[shock_in_next, shock_out_next, past_in_next, past_out_next](ap1, i_z)

    end

    return expected_value
end
