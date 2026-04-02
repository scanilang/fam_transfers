

function famtransfer_simulate(model, n_agents, seed)
    
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
    fam_type_mat = zeros(Int, n_agents)

    # Precompute draws
    survival_draws = rand(rng, n_agents, n_periods)
    shock_in_draws = rand(rng, n_agents, n_periods)
    shock_out_draws = rand(rng, n_agents, n_periods)
    z_draws = rand(rng, n_agents, n_periods)

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
            # Current state variables used for solving 
            # ---------------------------
            R = Race_mat[i]
            t = fam_type_mat[i]  
            m = marital_status[i]
            n = fam_size[i]
            e = ed_type[i]
            a = assets[i,period]
            d = debt[i,period]
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
                    debt[i, period] = 0.0
                else
                    assets[i, period] = 0.0
                    debt[i, period] = -ap1      
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
        debt = debt, 
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

    vj_itp = [LinearInterpolation((a_grid, z_grid[R]), Vj[R, m, n, e, :, :, shock_in, shock_out, past_in, past_out], extrapolation_bc=Interpolations.Flat()) 
                                for R in Race, m in marital_status, n in fam_size, t in fam_type, e in ed_type, shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2];

    pfj_itp = [LinearInterpolation((a_grid, z_grid[R]), PFj[R, m, n, e, :, :, shock_in, shock_out, past_in, past_out], extrapolation_bc=Interpolations.Flat()) 
                                for R in Race, m in marital_status, n in fam_size, t in fam_type, e in ed_type, shock_in in 1:2, shock_out in 1:2, past_in in 1:2, past_out in 1:2];

    return vj_itp, pfj_itp
end

function education_decision(model, R, t)
    v_1 = vj_itp[R, 1, n, t, 1, shock_in, shock_out, past_in, past_out](a, z) # education choice 1
    v_2 = vj_itp[R, 1, n, t, 2, shock_in, shock_out, past_in, past_out](a, z) # education choice 2
    v_3 = vj_itp[R, 1, n, t, 3, shock_in, shock_out, past_in, past_out](a, z) # education choice 3
    return argmax((v_1, v_2, v_3))
    
end