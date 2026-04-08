function model_create(;
    gamma = 2.0,
    beta = 0.96,
    r = 0.04,
    rb = 0.05,
    ra_w = 6.09,
    ra_b = 2.91,

    survival_risk = survival_risk_full,

    # Demographic variables
    Race = [1, 2],
    marital_status = [1, 2],
    fam_size = [1, 2, 3, 4, 5, 6],
    ed_type = [1, 2, 3],
    jpnts = 85 - 18,
    working_years = 60 - 18,

    # Income variables
    zpnts = 7,  # Number of income grid points
    ρ = [0.807, 0.772],   # [white, black]
    eps_std = [0.436, 0.535],

    # Asset variables
    apnts = 13,
    a_min = 7,
    a_max = 14,

    # Government
    tax_a = 0.2, # capital income tax


    )
    
    # 2 year interest rates
    #rb = (1 + r)^2 - 1
    #ra_w = (1 + ra_w)^2 - 1
    #ra_b = (1 + ra_b)^2 - 1

    # Income Process
    MC_1 = rouwenhorst(zpnts, ρ[1] , eps_std[1], 0) 
    MC_2 =  rouwenhorst(zpnts, ρ[2] , eps_std[2], 0) 
    Pimat = [MC_1.p, MC_2.p]
    z_grid = [exp.(MC_1.state_values), exp.(MC_2.state_values)]

    # Asset grid
    a_grid = [0.0; 500.0; polyexpandgrid(apnts-2, 2000.0, exp(a_max), 2.75)]

    # State space for parallelization
    tasks_idx = Vector{NTuple{6, Int64}}()
    for R in Race, m in marital_status, n in fam_size, e in ed_type, i_a in 1:apnts, i_z in 1:zpnts
        push!(tasks_idx, (R, m, n, e, i_a, i_z))
    end

    # Borrowing limit
    d_limit = zeros(jpnts, length(Race), length(ed_type), length(marital_status), 2)  

    for R in Race, e in ed_type, m in marital_status, dc in 1:2
        d_limit[:, R, e, m, dc] = compute_natural_borrowing_limit(model, R, e, m, dc)
    end

    max_debt = maximum(d_limit)  # largest possible debt across all types
    d_points = [-max_debt, -max_debt*0.75, -max_debt*0.5, -max_debt*0.25, -5000.0, -1000.0, 0.0]
    school_a_grid = [d_points; a_grid[2:end]]

    # precompute resources after shocks and probability of shocks
    shock_resources = zeros(Float64, working_years, 2, 2, 6, 4, 3,  apnts, zpnts, 2, 2, 2, 2) # (j, R, m, n, t, e, i_a, i_z, shock_in, shock_out, past_in, past_out)
    net_transfers = zeros(Float64, working_years, 2, 2, 6, 4, 3,  apnts, zpnts, 2, 2, 2, 2) # (j, R, m, n, t, e, i_a, i_z, shock_in, shock_out, past_in, past_out)
    prob_shocks = zeros(Float64, working_years, 2, 2, 4, 6, 3, apnts, zpnts, 2, 2, 2, 2) # (j, R, m, n,t,  e, i_a, i_z, shock_in, shock_out, past_in, past_out)

    for j in 1:working_years, R in Race, m in marital_status, n in fam_size, e in ed_type, t in fam_type
        for i_a in 1:apnts
            a = a_grid[i_a]
            a_income = R == 1 ? a * ra_w : a * ra_b
            for i_z in 1:zpnts
                for shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2
                    if j < working_years
                        y = g(R, j, e, m) * z_grid[R][i_z]
                        y_tax = tax_y(y, m)
                    else
                        y = g(R, 22, e, m) * z_grid[R][i_z] * 0.4 # assume income drops to 40% of last working year in retirement
                        y_tax = 0.0
                    end

                    shock_in_amount = transfers_in_amount(r, n, m, j, y, a_income, e, t)
                    shock_out_amount = transfers_out_amount(r,n,m,j,y, a_income, e, t)
                    prob_in = shocks_in_prob(r, n, m, j, y, a_income, e, t, past_in, past_out)
                    prob_out = shocks_out_prob(r, n, m, j, y, a_income, e, t, past_in, past_out)

                    if R == 1
                        if a > 0
                            a_next = a * (1 + ra_w * (1 - tax_a))
                        else
                            a_next = a * (1+rb)
                        end
                    else
                        if a > 0
                            a_next = a * (1 + ra_b * (1 - tax_a))
                        else
                            a_next = a * (1+rb)
                        end
                    end
                    if shock_in == 2
                        shock_in_prob = prob_in
                    else
                        shock_in_prob = 1 - prob_in
                    end
                    if shock_out == 2
                        shock_out_prob = prob_out
                    else
                        shock_out_prob = 1 - prob_out
                    end

                    
                    shock_resources[j, R, m, n,t, e, i_a, i_z, shock_in, shock_out, past_in, past_out] = y - y_tax + shock_in_amount + shock_out_amount + a_next
                    net_transfers[j, R, m, n,t, e, i_a, i_z, shock_in, shock_out, past_in, past_out] = shock_in_amount + shock_out_amount
                    prob_shocks[j, R, m, n, t, e, i_a, i_z, shock_in, shock_out, past_in, past_out] = shock_in_prob * shock_out_prob
                end
            end
        end
    end     


    return (; r, rb, ra_w, ra_b, gamma, beta, tax_a, survival_risk, Pimat, z_grid, a_grid, school_a_grid, d_limit, tasks_idx, shock_resources, net_transfers, prob_shocks)
end
