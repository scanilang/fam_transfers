function solve_model(model)
    (; apnts_nc, apnts_c, zpnts, working_years, jpnts, fam_shock_period) = model

    # Dimensions
    # j ranges: retirement 44-68 (25 periods), working 1-43 (43 periods)
    # fam_shock_period = 7 (age 25 = j=8, so last pre-family period is j=7)
    # n_retirement = jpnts - working_years = 68 - 43 = 25

    n_retirement = jpnts - working_years

    # -----------------------------------------------------------------------
    # NO COLLEGE PATH
    # -----------------------------------------------------------------------

    # --- Retirement: W_nc (ages 61-85, j=44 to j=68, indexed 1 to 25) ---
    # State: (R, m, n, t, i_a, i_z, shocks)
    Wj_nc  = zeros(Float32, n_retirement, 2, 2, 5, 3, apnts_nc, zpnts, 2, 2, 2, 2)
    WPFj_nc = copy(Wj_nc)

    # Terminal period j=68 (index n_retirement)
    Wj_nc[n_retirement, :,:,:,:,:,:,:,:,:,:], WPFj_nc[n_retirement, :,:,:,:,:,:,:,:,:,:] = Wncj(nothing, model, jpnts)

    for j in (jpnts-1):-1:(working_years+1)
        idx = j - working_years
        W_jp1 = @view Wj_nc[idx+1, :,:,:,:,:,:,:,:,:,:]
        Wj_nc[idx, :,:,:,:,:,:,:,:,:,:], WPFj_nc[idx, :,:,:,:,:,:,:,:,:,:] =
            Wncj(W_jp1, model, j)
    end

    # --- Post-family-shock working: Vnc2 (j = fam_shock_period+1 to working_years) ---
    # = j=8 to j=43 (36 periods)
    # State: (R, m, n, t, i_a, i_z, shocks)
    n_nc2 = working_years - fam_shock_period   # = 43 - 7 = 36
    Vj_nc2  = zeros(Float32, n_nc2, 2, 2, 5, 3, apnts_nc, zpnts, 2, 2, 2, 2)
    PFj_nc2 = copy(Vj_nc2)

    # Last working period j=43 (idx n_nc2)
    W_ret_first = @view Wj_nc[1, :,:,:,:,:,:,:,:,:,:]  # j=44
    Vj_nc2[n_nc2, :,:,:,:,:,:,:,:,:,:], PFj_nc2[n_nc2, :,:,:,:,:,:,:,:,:,:] =
        Vnc2_solve(nothing, W_ret_first, model, working_years)

    for j in (working_years-1):-1:(fam_shock_period+1)
        idx = j - fam_shock_period
        V_jp1 = @view Vj_nc2[idx+1, :,:,:,:,:,:,:,:,:,:]
        Vj_nc2[idx, :,:,:,:,:,:,:,:,:,:], PFj_nc2[idx, :,:,:,:,:,:,:,:,:,:] =
            Vnc2_solve(V_jp1, nothing, model, j)
    end

    # --- Pre-family-shock working: Vnc1 (j=1 to j=fam_shock_period) ---
    # State: (R, t, i_a, i_z, shocks) — m=1, n=1 fixed
    n_nc1 = fam_shock_period  # = 7
    Vj_nc1  = zeros(Float32, n_nc1, 2, 3, apnts_nc, zpnts, 2, 2, 2, 2)
    PFj_nc1 = copy(Vj_nc1)

    # Last pre-family period j=7 (idx n_nc1): continuation is Vnc2[1] (j=8) with family shock
    Vnc2_first = @view Vj_nc2[1, :,:,:,:,:,:,:,:,:,:]
    Vj_nc1[n_nc1, :,:,:,:,:,:,:,:,:], PFj_nc1[n_nc1, :,:,:,:,:,:,:,:,:] =
        Vnc1j_solve(nothing, Vnc2_first, model, fam_shock_period)

    for j in (fam_shock_period-1):-1:1
        V_jp1 = @view Vj_nc1[j+1, :,:,:,:,:,:,:,:,:]
        Vj_nc1[j, :,:,:,:,:,:,:,:,:], PFj_nc1[j, :,:,:,:,:,:,:,:,:] =
            Vnc1j_solve(V_jp1, nothing, model, j)
    end

    # -----------------------------------------------------------------------
    # COLLEGE PATH
    # -----------------------------------------------------------------------

    # --- Retirement: W_c ---
    # State: (R, m, n, t, e_idx, i_a, i_z, shocks)
    Wj_c  = zeros(Float32, n_retirement, 2, 2, 5, 3, 2, apnts_c, zpnts, 2, 2, 2, 2)
    WPFj_c = copy(Wj_c)

    Wj_c[n_retirement, :,:,:,:,:,:,:,:,:,:,:], WPFj_c[n_retirement, :,:,:,:,:,:,:,:,:,:,:] =
        Wcj(nothing, model, jpnts)

    for j in (jpnts-1):-1:(working_years+1)
        idx = j - working_years
        W_jp1 = @view Wj_c[idx+1, :,:,:,:,:,:,:,:,:,:,:]
        Wj_c[idx, :,:,:,:,:,:,:,:,:,:,:], WPFj_c[idx, :,:,:,:,:,:,:,:,:,:,:] =
            Wcj(W_jp1, model, j)
    end

    # --- Post-family-shock college working: Vc2 (j=8 to j=43) ---
    # State: (R, m, n, t, e_idx, i_a, i_z, shocks)
    n_c2 = working_years - fam_shock_period  # = 36
    Vj_c2  = zeros(Float32, n_c2, 2, 2, 5, 3, 2, apnts_c, zpnts, 2, 2, 2, 2)
    PFj_c2 = copy(Vj_c2)

    W_ret_first_c = @view Wj_c[1, :,:,:,:,:,:,:,:,:,:,:]
    Vj_c2[n_c2, :,:,:,:,:,:,:,:,:,:,:], PFj_c2[n_c2, :,:,:,:,:,:,:,:,:,:,:] =
        Vc2j_solve(nothing, W_ret_first_c, model, working_years)

    for j in (working_years-1):-1:(fam_shock_period+1)
        idx = j - fam_shock_period
        V_jp1 = @view Vj_c2[idx+1, :,:,:,:,:,:,:,:,:,:,:]
        Vj_c2[idx, :,:,:,:,:,:,:,:,:,:,:], PFj_c2[idx, :,:,:,:,:,:,:,:,:,:,:] =
            Vc2j_solve(V_jp1, nothing, model, j)
    end

    # --- Pre-family-shock college working: Vc1 (j=3 or j=5 through j=7) ---
    # 2yr grads enter at j=3, 4yr grads enter at j=5
    # Solve for all j=1 to j=7 to simplify (unreachable slots stay zero)
    # State: (R, t, degree, i_a, i_z, shocks)
    n_c1 = fam_shock_period  # = 7
    Vj_c1  = zeros(Float32, n_c1, 2, 3, 2, apnts_c, zpnts, 2, 2, 2, 2)
    PFj_c1 = copy(Vj_c1)

    # j=7: continuation is Vc2[1] (j=8) with family shock
    Vc2_first = @view Vj_c2[1, :,:,:,:,:,:,:,:,:,:,:]
    Vj_c1[n_c1, :,:,:,:,:,:,:,:,:,:], PFj_c1[n_c1, :,:,:,:,:,:,:,:,:,:] =
        Vc1j_solve(nothing, Vc2_first, model, fam_shock_period)

    for j in (fam_shock_period-1):-1:1
        V_jp1 = @view Vj_c1[j+1, :,:,:,:,:,:,:,:,:,:]
        Vj_c1[j, :,:,:,:,:,:,:,:,:,:], PFj_c1[j, :,:,:,:,:,:,:,:,:,:] =
            Vc1j_solve(V_jp1, nothing, model, j)
    end

    # -----------------------------------------------------------------------
    # SCHOOL PATH
    # -----------------------------------------------------------------------
    # 2yr: school j=1,2 → graduate, enter Vj_c1[3] at j=3
    # 4yr: school j=1,2,3,4 → graduate, enter Vj_c1[5] at j=5

    # Vsj dims: (j, R, t_own, degree, i_a_school)
    n_school_periods = 4
    Vsj  = zeros(Float32, n_school_periods, 2, 3, 2, apnts_c)
    PFsj = copy(Vsj)

    # j=4: last period for 4yr graduate. Continuation is Vj_c1[5] (first working period for 4yr)
    Vc1_j5 = @view Vj_c1[5, :,:,:,:,:,:,:,:,:,:]
    Vsj[4, :,:,:,:], PFsj[4, :,:,:,:] = VSj_enrolled(
        zeros(Float32, 2, 3, 2, apnts_c),   # vsjp1 not used since 4yr graduates here
        Vc1_j5,
        Vsj[4, :,:,:,:], PFsj[4, :,:,:,:],
        model, 4)

    # j=3: 4yr still enrolled, no graduation
    Vsj[3, :,:,:,:], PFsj[3, :,:,:,:] = VSj_enrolled(
        @view(Vsj[4, :,:,:,:]),
        nothing,
        Vsj[3, :,:,:,:], PFsj[3, :,:,:,:],
        model, 3)

    # j=2: 2yr graduates (enters Vj_c1[3]); 4yr still enrolled
    Vc1_j3 = @view Vj_c1[3, :,:,:,:,:,:,:,:,:,:]
    Vsj[2, :,:,:,:], PFsj[2, :,:,:,:] = VSj_enrolled(
        @view(Vsj[3, :,:,:,:]),
        Vc1_j3,
        Vsj[2, :,:,:,:], PFsj[2, :,:,:,:],
        model, 2)

    # j=1: first period with education transfer realization
    # Vsj_1 dims: (R, m_p, t_parent, e_p, i_a_p, i_a, i_z_p, edu_shock, degree)
    # Matches tasks_idx_s1 in setup (apart from adding edu_shock, degree)
    Vsj_1  = zeros(Float32, 2, 2, 3, 3, apnts_nc, apnts_nc, zpnts, 2, 2)
    PFsj_1 = copy(Vsj_1)

    Vsj_1, PFsj_1 = VSj_first_period(
        @view(Vsj[2, :,:,:,:]),
        Vsj_1, PFsj_1, model)

    # -----------------------------------------------------------------------
    # Return
    # -----------------------------------------------------------------------
    return (
        # No college
        Vj_nc1  = Vj_nc1,  PFj_nc1 = PFj_nc1,
        Vj_nc2  = Vj_nc2,  PFj_nc2 = PFj_nc2,
        Wj_nc   = Wj_nc,   WPFj_nc = WPFj_nc,
        # College working
        Vj_c1   = Vj_c1,   PFj_c1  = PFj_c1,
        Vj_c2   = Vj_c2,   PFj_c2  = PFj_c2,
        Wj_c    = Wj_c,    WPFj_c  = WPFj_c,
        # School
        Vsj     = Vsj,     PFsj    = PFsj,
        Vsj_1   = Vsj_1,   PFsj_1  = PFsj_1,
    )
end