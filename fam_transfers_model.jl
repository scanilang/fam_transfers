function model_create(;
    gamma = 2.0,
    beta = 0.96,
    r = 0.04,
    rb = 0.05,
    ra_w = 6.09,
    ra_b = 2.91,
    tuition_2yr = 5000.0, #yearly tuition cost for 2-year degree
    tuition_4yr = 15000.0, #yearly tuition cost for 4-year degree
    r_loan = 0.06,
    survival_risk = survival_risk_full,
    Race = [1, 2],
    marital_status = [1, 2],
    fam_type = [1, 2, 3],
    fam_size = [1, 2, 3, 4, 5],
    ed_type = [1, 2, 3],
    jpnts = 85 - 18,
    working_years = 60 - 18,
    fam_shock_period = 25 - 18,
    zpnts = 7,
    ρ = [0.807, 0.772],
    eps_std = [0.436, 0.535],
    apnts = 13,
    a_min = 7,
    a_max = 14,
    tax_a = 0.2,
    family_shock_probs = family_shock_probs
    )

    # Income Process
    MC_1 = rouwenhorst(zpnts, ρ[1], eps_std[1], 0)
    MC_2 = rouwenhorst(zpnts, ρ[2], eps_std[2], 0)
    Pimat = [MC_1.p, MC_2.p]
    z_grid = [exp.(MC_1.state_values), exp.(MC_2.state_values)]

    # Asset grids — DEFINE FIRST
    a_grid = [0.0; 500.0; polyexpandgrid(apnts-2, 2000.0, exp(a_max), 2.75)]
    a_grid_nocollege = a_grid  # non-negative only

    # Borrowing limit FIRST (needed for school_a_grid)
    phi = 0.22
    d_limit = zeros(jpnts, 2, 2)
    for R in Race, e in [2, 3]
        e_idx = e - 1
        d_limit[:, R, e_idx] = compute_borrowing_limit(z_grid, r_loan, jpnts, R, e, phi)
    end

    max_debt = maximum(d_limit)
    d_points = [-max_debt, -max_debt*0.75, -max_debt*0.5, -max_debt*0.25, -5000.0, -1000.0, 0.0]
    a_grid_college = [d_points; a_grid[2:end]]
    school_a_grid = a_grid_college  # same grid for school and college working

    # Grid sizes — NOW define these
    apnts_nc = length(a_grid_nocollege)
    apnts_c  = length(a_grid_college)

    # Task indices — NOW define these
    tasks_idx_nc2 = Vector{NTuple{6, Int64}}()
    for R in Race, m in marital_status, n in fam_size, t in fam_type, i_a in 1:apnts_nc, i_z in 1:zpnts
        push!(tasks_idx_nc2, (R, m, n, t, i_a, i_z))
    end

    tasks_idx_nc1 = Vector{NTuple{6, Int64}}()
    for R in Race, t in fam_type, i_a in 1:apnts_nc, i_z in 1:zpnts
        push!(tasks_idx_nc1, (R, t, i_a, i_z))
    end

    tasks_idx_c2 = Vector{NTuple{7, Int64}}()
    for R in Race, m in marital_status, n in fam_size, t in fam_type, degree in 1:2, i_a in 1:apnts_c, i_z in 1:zpnts
        push!(tasks_idx_c2, (R, m, n, t, degree, i_a, i_z))
    end

    tasks_idx_c1 = Vector{NTuple{7, Int64}}()
    for R in Race, t in fam_type, degree in 1:2, i_a in 1:apnts_c, i_z in 1:zpnts
        push!(tasks_idx_c1, (R, t, degree, i_a, i_z))
    end

    tasks_idx_s1 = Vector{NTuple{7, Int64}}()
    for R in Race, m in marital_status, t in fam_type, e in 1:3, i_a_p in 1:apnts_nc, i_a in 1:apnts_nc, i_z in 1:zpnts
        push!(tasks_idx_s1, (R, m, t, e, i_a_p, i_a, i_z))
    end

    # School task index — simple (R, t, i_a) only
    tasks_idx_s2 = Vector{NTuple{4, Int64}}()
    for R in Race, t in fam_type, degree in 1:2, i_a in 1:apnts_c
        push!(tasks_idx_s2, (R, t, degree, i_a))
    end

    # Precompute y_values
    y_values = zeros(Float64, 2, jpnts, 2, 3, zpnts)
    for R in Race, j in 1:jpnts, m in marital_status, e in ed_type, i_z in 1:zpnts
        if j <= working_years
            y_values[R, j, m, e, i_z] = g(R, j, e, m) * z_grid[R][i_z]
        else
            y_values[R, j, m, e, i_z] = g(R, working_years, e, m) * z_grid[R][i_z] * 0.4
        end
    end

    # Precomputation for no-college and college paths
    shock_resources_nc = zeros(Float64, jpnts, 2, 2, 5, 3, apnts_nc, zpnts, 2, 2, 2, 2)
    net_transfers_nc   = zeros(Float64, jpnts, 2, 2, 5, 3, apnts_nc, zpnts, 2, 2, 2, 2)
    prob_shocks_nc     = zeros(Float64, jpnts, 2, 2, 5, 3, apnts_nc, zpnts, 2, 2, 2, 2)

    for j in 1:jpnts, R in Race, m in marital_status, n in fam_size, t in fam_type
        for i_a in 1:apnts_nc
            a        = a_grid_nocollege[i_a]
            r        = R == 1 ? ra_w : ra_b
            a_income = a * r
            a_next   = a * (1 + r * (1 - tax_a))

            for i_z in 1:zpnts
                y     = y_values[R, j, m, 1, i_z]
                y_tax = j <= working_years ? tax_y(y, m) : 0.0

                for shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2
                    if j <= working_years
                        shock_in_amount  = transfers_in_amount(R, n, m, j, y, a_income, 1, t)
                        shock_out_amount = transfers_out_amount(R, n, m, j, y, a_income, 1, t)
                        prob_in  = shocks_in_prob(R, n, m, j, y, a_income, 1, t, past_in, past_out)
                        prob_out = shocks_out_prob(R, n, m, j, y, a_income, 1, t, past_in, past_out)
                    else
                        shock_in_amount  = 0.0
                        shock_out_amount = 0.0
                        prob_in          = 0.0
                        prob_out         = 0.0
                    end

                    shock_in_prob  = shock_in  == 2 ? prob_in  : 1 - prob_in
                    shock_out_prob = shock_out == 2 ? prob_out : 1 - prob_out
                    transfer_in    = (shock_in  - 1) * shock_in_amount

                    # Apply cap to transfer_out based on realized transfer_in
                    shock_out_effective = effective_transfer_out(m,n, shock_out_amount, y, y_tax, a_income, transfer_in)   
                    transfer_out   = (shock_out - 1) * shock_out_effective

                    shock_resources_nc[j, R, m, n, t, i_a, i_z, shock_in, shock_out, past_in, past_out] =
                        y - y_tax + transfer_in - transfer_out + a_next
                    net_transfers_nc[j, R, m, n, t, i_a, i_z, shock_in, shock_out, past_in, past_out] =
                        transfer_in - transfer_out
                    prob_shocks_nc[j, R, m, n, t, i_a, i_z, shock_in, shock_out, past_in, past_out] =
                        shock_in_prob * shock_out_prob
                end
            end
        end
    end

    shock_resources_c = zeros(Float64, jpnts, 2, 2, 5, 3, 2, apnts_c, zpnts, 2, 2, 2, 2)
    net_transfers_c   = zeros(Float64, jpnts, 2, 2, 5, 3, 2, apnts_c, zpnts, 2, 2, 2, 2)
    prob_shocks_c     = zeros(Float64, jpnts, 2, 2, 5, 3, 2, apnts_c, zpnts, 2, 2, 2, 2)

    for j in 1:jpnts, R in Race, m in marital_status, n in fam_size, t in fam_type, e in 2:3
        e_idx = e - 1
        for i_a in 1:apnts_c
            a = a_grid_college[i_a]
            r = R == 1 ? ra_w : ra_b
            if a >= 0
                a_income = a * r
                a_next   = a * (1 + r * (1 - tax_a))
            else
                a_income = 0.0
                a_next   = a * (1 + r_loan)
            end

            for i_z in 1:zpnts
                y     = y_values[R, j, m, e, i_z]
                y_tax = j <= working_years ? tax_y(y, m) : 0.0

                for shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2
                    if j <= working_years
                        shock_in_amount  = transfers_in_amount(R, n, m, j, y, a_income, e, t)
                        shock_out_amount = transfers_out_amount(R, n, m, j, y, a_income, e, t)
                        prob_in  = shocks_in_prob(R, n, m, j, y, a_income, e, t, past_in, past_out)
                        prob_out = shocks_out_prob(R, n, m, j, y, a_income, e, t, past_in, past_out)
                    else
                        shock_in_amount  = 0.0
                        shock_out_amount = 0.0
                        prob_in          = 0.0
                        prob_out         = 0.0
                    end

                    shock_in_prob  = shock_in  == 2 ? prob_in  : 1 - prob_in
                    shock_out_prob = shock_out == 2 ? prob_out : 1 - prob_out
                    transfer_in    = (shock_in  - 1) * shock_in_amount

                    # Apply cap to transfer_out based on realized transfer_in
                    c_floor = m == 1 ? 5400.0 : 10800.0
                    shock_out_effective = effective_transfer_out(m,n, shock_out_amount, y, y_tax, a_income, transfer_in)   
                    transfer_out   = (shock_out - 1) * shock_out_effective

                    shock_resources_c[j, R, m, n, t, e_idx, i_a, i_z, shock_in, shock_out, past_in, past_out] =
                        y - y_tax + transfer_in - transfer_out + a_next
                    net_transfers_c[j, R, m, n, t, e_idx, i_a, i_z, shock_in, shock_out, past_in, past_out] =
                        transfer_in - transfer_out
                    prob_shocks_c[j, R, m, n, t, e_idx, i_a, i_z, shock_in, shock_out, past_in, past_out] =
                        shock_in_prob * shock_out_prob
                end
            end
        end
    end

    return (; tuition_2yr, tuition_4yr, r_loan,
              a_grid_nocollege, a_grid_college, school_a_grid, apnts_nc, apnts_c,
              working_years, jpnts, fam_shock_period,
              family_shock_probs, fam_type, r, rb, ra_w, ra_b, gamma, beta, tax_a,
              survival_risk, Pimat, z_grid, a_grid, d_limit,
              tasks_idx_nc1, tasks_idx_nc2, tasks_idx_s2, tasks_idx_c1, tasks_idx_c2, tasks_idx_s1,
              y_values, shock_resources_nc, net_transfers_nc, prob_shocks_nc,
              shock_resources_c, net_transfers_c, prob_shocks_c,
              Race, marital_status, fam_size, ed_type, zpnts)
end