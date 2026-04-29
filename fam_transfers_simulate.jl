
function education_decision(model, Vsj_1, vjp1, R, m, n, t, e, i_a, i_z)
    (; a_grid_nocollege, y_values) = model
    y = y_values[R, 43, m, e, i_z]
    a_income = R == 1 ? a_grid_nocollege[i_a] * ra_w : a_grid_nocollege[i_a] * ra_b

    # No college: enter working phase immediately at j=1 with e=1
    # Need to integrate over initial z and transfer shocks for no-college
    v_nocollege = EVnc_jp1(model, vjp1, 0, R, 0, 1, a_grid_nocollege[i_a], i_z, 1, 1, 1, 1)

    # 2yr: integrate over education transfer shock
    prob_help_2yr = edu_transfer_prob(R, m,n, y, a_income, e, t, 2)  # degree_choice=2
    v_2yr = prob_help_2yr       * Vsj_1[R, m, e, i_a, i_z, 2, 1] +  # help
            (1 - prob_help_2yr) * Vsj_1[R, m, e, i_a, i_z, 1, 1]     # no help

    # 4yr: integrate over education transfer shock
    prob_help_4yr = edu_transfer_prob(R, m,n, y, a_income, e, t, 4)  # degree_choice=4
    v_4yr = prob_help_4yr       * Vsj_1[R, m, e, i_a, i_z, 2, 2] +  # help
            (1 - prob_help_4yr) * Vsj_1[R, m, e, i_a, i_z, 1, 2]     # no help

    return argmax((v_nocollege, v_2yr, v_4yr))
end

function family_shock(model, R, e, t_parent, rng)
    # t_parent: parental family type (1=both_low, 2=both_mid, 3=both_high, 4=single)
    # Returns: (m_next, n_next, t_next) for the agent at age 25
    
    (; fam_transition_matrix) = model
    # fam_transition_matrix[R, e, t_parent] is a vector of probabilities over 
    # all valid (m, n, t) combinations

    outcomes = [(m, n, t) for m in 1:2, n in 1:6, t in 1:4
                if !(m == 1 && n > 1) && !(m == 2 && n < 1)]
    
    probs = fam_transition_matrix[R, e, t_parent]  # length = number of valid outcomes
    u = rand(rng)
    cumprob = 0.0
    for (i, p) in enumerate(probs)
        cumprob += p
        if u <= cumprob
            return outcomes[i]
        end
    end
    return outcomes[end]  # fallback
end

function famtransfer_simulate(model, n_agents, n_periods, random_seed)
    
    # Model Parameters
    (; Race, marital_status, fam_size, ed_type, fam_type, a_grid, z_grid, ra_w, ra_b, tax_y, tax_a, transfers_in_amount, transfers_out_amount, shocks_in_prob, shocks_out_prob) = model
    
    rng = MersenneTwister(random_seed)

    # ---------------------------
    # Set Up
    # ---------------------------

    # Simulation result placeholders
    total_periods = n_periods + 1
    assets = zeros(n_agents, total_periods + 1)
    debt = zeros(n_agents, total_periods + 1)
    income = zeros(n_agents, total_periods)
    age = zeros(Int, n_agents, total_periods)
    z_mat = zeros(n_agents, total_periods)
    z_idx_mat = zeros(Int, n_agents, total_periods)
    past_in_flag = zeros(Int, n_agents, total_periods)
    past_out_flag = zeros(Int, n_agents, total_periods)
    shock_in_flag = zeros(Int, n_agents, total_periods)
    shock_out_flag = zeros(Int, n_agents, total_periods)
    transfer_in_mat = zeros(n_agents, total_periods)
    transfer_out_mat = zeros(n_agents, total_periods)
    Race_mat = zeros(Int, n_agents)
    fam_type_mat = zeros(n_agents, total_periods)
    e_mat = zeros(Int, n_agents)

    # Precompute draws
    survival_draws = rand(rng, n_agents, n_periods)
    shock_in_draws = rand(rng, n_agents, n_periods)
    shock_out_draws = rand(rng, n_agents, n_periods)
    z_draws = rand(rng, n_agents, n_periods)

    # ---------------------------
    # Initial conditions
    # ---------------------------

    # ---------------------------
    # Education decision at j=0
    # ---------------------------
    for i in 1:n_agents
        R = Race_mat[i]
        m = marital_status[i]
        t = fam_type_mat[i, 1]
        e = ed_type[i]
        i_a = searchsortedfirst(a_grid_nocollege, assets[i, 1])
        e_mat[i] = education_decision(model, Vsj_1, vjp1, R, m, t, e, i_a)
    end
    # ---------------------------
    # Simulation loop
    # ---------------------------
    for period in 2:total_periods
        prev_period = period - 1
        next_period = period + 1
        j = period - 1

        age_j = floor(Int, j + 17) # age starts at 18 in period 1
        @. age[:,period] = ifelse(age[:,period] != -1, age_j, -1)

        # ---------------------------
        # Load interpolations ONCE for this j
        # ---------------------------            

        if j <= jpnts
            vj, pfj = get_itp_function(model, j, savepath);
        end

        @threads for i in 1:n_agents
            # Skip dead agent
            if age[i,period] == -1
                continue 
            end

            # ---------------------------
            # Family Shock
            # ---------------------------
            if age[i,period] == 25
                m, n, t = family_shock() 
            end

            # ---------------------------
            # Current state variables used for solving 
            # ---------------------------
            R = Race_mat[i]
            t = fam_type_mat[i]  
            m = marital_status[i]
            n = fam_size[i]
            e = ed_type[i]
            a = assets[i,period]
            z = z_mat[i,period]

            # ---------------------------
            # Update z, income, and post-tax resources
            # ---------------------------
            if j >= 1 && j < 121
                # working and quarterly income
                idx = z_idx_mat[i, prev_period]
                next_state = min(searchsortedfirst(z_pimat[R][idx,:], z_draws[i, period]), zpnts)
                z_idx_mat[i, period] = next_state
                z_mat[i, period] = z_grid[R][next_state] 
                income[i,period] = g(r,j,e, m) * z_mat[i,period]
            elseif j == 121
                # retirement and quarterly income
                z_mat[i,period] = z_mat[i,prev_period]
                income[i,period] = income[i,prev_period] * 0.4
            elseif j > 121
                z_mat[i,period] = z_mat[i,prev_period]
                income[i,period] = income[i,prev_period] 
            end


            # ---------------------------
            # Working age decisions (j ≤ 140)
            # ---------------------------
            if j <= jpnts
                # --- Assign value function interpolators for agent ---
                vj_itp = vj[R,m, n,t,e,:,:,:,:,:,:];
                pfj_itp = pfj[R,m, n,t,e,:,:,:,:,:,:];
                # --- Solve for optimal decision ---
                ap1 = pfj_itp(a, z_mat[i,period])
                # --- Update assets and debt based on optimal decision ---
                if ap1 >= 0
                    assets[i, period] = ap1
                else
                    assets[i, period] = 0.0
                end
            end

            # ---------------------------
            # Update survival
            # ---------------------------
            if survival_draws[i,prev_period] <= (1 - survival_risk_full[j,R]) && next_period <= total_periods
                age[i, next_period] = -1
            end
        end
    end

    # ---------------------------
    # Post-process results
    # ---------------------------

    # Pack results
    simulation_results = (
        Race_mat = Race_mat, 
        assets = assets, 
        fam_type_mat = fam_type_mat, 
        age = age, 
        income = income, 
        z_mat = z_mat,
        z_idx_mat = z_idx_mat,
        shock_in_flag = shock_in_flag,
        shock_out_flag = shock_out_flag,
        transfer_in_mat = transfer_in_mat,
        transfer_out_mat = transfer_out_mat,
        past_in_flag = past_in_flag,
        past_out_flag = past_out_flag,

    )
    return simulation_results
end

function get_itp_function(model, j, savepath)

    (; Race, fam_type, z_grid, marital_status, fam_size, ed_type, a_grid, d_grid) = model

    data = open(deserialize, joinpath(savepath, "File_7_age_$(j).jls"));
    Vj = data.Vj;
    PFj = data.PFj;

    vj_itp = [LinearInterpolation((a_grid, z_grid[R]), Vj[R, m, n, t, e, :, :, shock_in, shock_out, past_in, past_out], extrapolation_bc=Interpolations.Flat()) 
                                for R in Race, m in marital_status, n in fam_size, t in fam_type, e in ed_type, shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2];

    pfj_itp = [LinearInterpolation((a_grid, z_grid[R]), PFj[R, m, n, t, e, :, :, shock_in, shock_out, past_in, past_out], extrapolation_bc=Interpolations.Flat()) 
                                for R in Race, m in marital_status, n in fam_size, t in fam_type, e in ed_type, shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2];

    return vj_itp, pfj_itp
end

