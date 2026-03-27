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
    Race_mat = zeros(Int, n_agents)
    fam_type_mat = zeros(Int, n_agents)

    # Precompute draws
    (; expenses, z_draws, z_pimat, b_draws, w7_draws, w13_draws, df_draws, survival_draws, wage_garnishment_draws) = precompute_draws(n_agents, n_periods, Race_mat, Race, fam_type_mat, eps_std , kpnts, Pimat, k_probs_quarterly, rng)




end



function precompute_draws(n_agents, n_periods, Race_mat, Race, fam_type_mat, Pimat, k_probs_quarterly, rng)

    # probability drawn precompute
    survival_draws = rand(rng, n_agents, n_periods)
    shock_in_draws = rand(rng, n_agents, n_periods)
    shock_out_draws = rand(rng, n_agents, n_periods)

    # initial z draw
    z_draws = rand(rng, n_agents, n_periods)
    z_pimat= [cumsum(Pimat[R], dims = 2) for R in Race]

    # precompute shock drawn based on probabilities

    prob_shocks[j, R, m, n, t, e, i_a, i_z, shock_in, shock_out, past_in, past_out]
    
    return (; survival_draws, shock_in_draws, shock_out_draws, z_draws, z_pimat, expenses)
end
