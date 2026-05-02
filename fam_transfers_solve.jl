function solve_model(model, savepath)
    (; apnts_nc, apnts_c, zpnts, working_years,n_retirement, jpnts, fam_shock_period) = model

    # -----------------------------------------------------------------------
    # NO COLLEGE PATH
    # -----------------------------------------------------------------------
    if !isfile(joinpath(savepath, "NoSchool.jls")) 
        
        # --- Retirement W_nc ---
        println("Solving No College Retirement...")
        # Dims: (j, R, m, n, t, i_a, i_z, shock_in, shock_out, past_in, past_out) = 10 non-j
        Wj_nc  = zeros(Float32, 2, 2, 5, 3, apnts_nc, zpnts, 2, 2, 2, 2)
        WPFj_nc = copy(Wj_nc)

        W_nc  = zeros(Float32, n_retirement, 2, 2, 5, 3, apnts_nc, zpnts, 2, 2, 2, 2)
        WPF_nc = copy(W_nc)

        W_nc[n_retirement,  :,:,:,:,:,:,:,:,:,:], WPF_nc[n_retirement, :,:,:,:,:,:,:,:,:,:] = Wncj(nothing, Wj_nc, WPFj_nc, model, jpnts)

        for j in (jpnts-1):-1:(working_years+1)
            idx = j - working_years
            W_jp1 = @view W_nc[idx+1, :,:,:,:,:,:,:,:,:,:]
            W_nc[idx,  :,:,:,:,:,:,:,:,:,:], 
            WPF_nc[idx, :,:,:,:,:,:,:,:,:,:] = Wncj(W_jp1, Wj_nc, WPFj_nc, model, j)
        end

        # --- Post-family-shock working Vnc2 ---
        println("Solving No College Post-Family-Shock Working...")
        # Dims: (j, R, m, n, t, i_a, i_z, shocks) = 10 non-j
        n_nc2 = working_years - fam_shock_period  # = 36
        V_nc2  = zeros(Float32, n_nc2, 2, 2, 5, 3, apnts_nc, zpnts, 2, 2, 2, 2)
        PF_nc2 = copy(V_nc2)

        Vj_nc2  = zeros(Float32, 2, 2, 5, 3, apnts_nc, zpnts, 2, 2, 2, 2)
        PFj_nc2 = copy(Vj_nc2)

        W_ret_first = @view W_nc[1, :,:,:,:,:,:,:,:,:,:]
        V_nc2[n_nc2, :,:,:,:,:,:,:,:,:,:], PF_nc2[n_nc2, :,:,:,:,:,:,:,:,:,:] = Vnc2j_solve(nothing, W_ret_first, Vj_nc2, PFj_nc2, model, working_years)

        for j in (working_years-1):-1:(fam_shock_period+1)
            idx = j - fam_shock_period
            V_jp1 = @view V_nc2[idx+1, :,:,:,:,:,:,:,:,:,:]
            V_nc2[idx, :,:,:,:,:,:,:,:,:,:], PF_nc2[idx, :,:,:,:,:,:,:,:,:,:] = Vnc2j_solve(V_jp1, nothing, Vj_nc2, PFj_nc2, model, j)
        end

        # --- Pre-family-shock working Vnc1 ---
        println("Solving No College Pre-Family-Shock Working...")
        # Dims: (j, R, t, i_a, i_z, shocks) = 8 non-j
        n_nc1 = fam_shock_period  # = 7
        Vj_nc1  = zeros(Float32, 2, 3, apnts_nc, zpnts, 2, 2, 2, 2)
        PFj_nc1 = copy(Vj_nc1)

        V_nc1  = zeros(Float32, n_nc1, 2, 3, apnts_nc, zpnts, 2, 2, 2, 2)
        PF_nc1 = copy(V_nc1)

        Vnc2_first = @view V_nc2[1, :,:,:,:,:,:,:,:,:,:]
        V_nc1[n_nc1, :,:,:,:,:,:,:,:], PF_nc1[n_nc1, :,:,:,:,:,:,:,:] = Vnc1j_solve(nothing, Vnc2_first, Vj_nc1, PFj_nc1, model, fam_shock_period)

        for j in (fam_shock_period-1):-1:1
            V_jp1 = @view V_nc1[j+1, :,:,:,:,:,:,:,:]
            V_nc1[j, :,:,:,:,:,:,:,:], PF_nc1[j, :,:,:,:,:,:,:,:] = Vnc1j_solve(V_jp1, nothing, Vj_nc1, PFj_nc1, model, j)
        end
        
        open(joinpath(savepath, "NoSchool.jls"), "w") do io
                serialize(io, (V_nc1 = V_nc1, PF_nc1 = PF_nc1,
            V_nc2 = V_nc2, PF_nc2 = PF_nc2,
            W_nc  = W_nc,  WPF_nc = WPF_nc))
            end
    end

    # -----------------------------------------------------------------------
    # COLLEGE PATH
    # -----------------------------------------------------------------------

    # --- Retirement W_c ---
    println("Solving College Retirement...")
    # Dims: (j, R, m, n, t, e_idx, i_a, i_z, shocks) = 11 non-j
    Wj_c  = zeros(Float32, 2, 2, 5, 3, 2, apnts_nc, zpnts, 2, 2, 2, 2)
    WPFj_c = copy(Wj_c)

    W_c  = zeros(Float32, n_retirement,2, 2, 5, 3, 2, apnts_nc, zpnts, 2, 2, 2, 2)
    WPF_c = copy(W_c)

    W_c[n_retirement,  :,:,:,:,:,:,:,:,:,:,:], WPF_c[n_retirement, :,:,:,:,:,:,:,:,:,:,:] = Wcj(nothing, Wj_c, WPFj_c, model, jpnts)

    for j in (jpnts-1):-1:(working_years+1)
        idx = j - working_years
        W_jp1 = @view W_c[idx+1, :,:,:,:,:,:,:,:,:,:,:]
        W_c[idx,  :,:,:,:,:,:,:,:,:,:,:], WPF_c[idx, :,:,:,:,:,:,:,:,:,:,:] = Wcj(W_jp1, Wj_c, WPFj_c, model, j)
    end

    # --- Post-family-shock college working Vc2 ---
    println("Solving College Post-Family-Shock Working...")
    # Dims: (j, R, m, n, t, e_idx, i_a, i_z, shocks) = 11 non-j
    n_c2 = working_years - fam_shock_period  # = 36
    Vj_c2  = zeros(Float32, 2, 2, 5, 3, 2, apnts_c, zpnts, 2, 2, 2, 2)
    PFj_c2 = copy(Vj_c2)

    V_c2  = zeros(Float32,n_c2, 2, 2, 5, 3, 2, apnts_c, zpnts, 2, 2, 2, 2)
    PF_c2 = copy(V_c2)

    W_ret_first_c = @view W_c[1, :,:,:,:,:,:,:,:,:,:,:]
    V_c2[n_c2, :,:,:,:,:,:,:,:,:,:,:], PF_c2[n_c2, :,:,:,:,:,:,:,:,:,:,:] = Vc2j_solve(nothing, W_ret_first_c, Vj_c2, PFj_c2, model, working_years)

    for j in (working_years-1):-1:(fam_shock_period+1)
        idx = j - fam_shock_period
        V_jp1 = @view V_c2[idx+1, :,:,:,:,:,:,:,:,:,:,:]
        V_c2[idx, :,:,:,:,:,:,:,:,:,:,:], PF_c2[idx, :,:,:,:,:,:,:,:,:,:,:] = Vc2j_solve(V_jp1, nothing, Vj_c2, PFj_c2, model, j)
    end

    # --- Pre-family-shock college working Vc1 ---
    println("Solving College Pre-Family-Shock Working...")
    # Dims: (j, R, t, degree, i_a, i_z, shocks) = 9 non-j
    n_c1 = fam_shock_period - 2  # = 7
    Vj_c1  = zeros(Float32, 2, 3, 2, apnts_c, zpnts, 2, 2, 2, 2)
    PFj_c1 = copy(Vj_c1)

    V_c1  = zeros(Float32, n_c1, 2, 3, 2, apnts_c, zpnts, 2, 2, 2, 2)
    PF_c1 = copy(V_c1)

    Vc2_first = @view V_c2[1, :,:,:,:,:,:,:,:,:,:,:]
    V_c1[n_c1, :,:,:,:,:,:,:,:,:], PF_c1[n_c1, :,:,:,:,:,:,:,:,:] = Vc1j_solve(nothing, Vc2_first, Vj_c1, PFj_c1, model, fam_shock_period)

    for j in (fam_shock_period-3):-1:1
        V_jp1 = @view V_c1[j+1, :,:,:,:,:,:,:,:,:]
        V_c1[j, :,:,:,:,:,:,:,:,:], PF_c1[j, :,:,:,:,:,:,:,:,:] = Vc1j_solve(V_jp1, nothing, Vj_c1, PFj_c1, model, j)
    end

    # -----------------------------------------------------------------------
    # SCHOOL PATH
    # -----------------------------------------------------------------------
    println("Solving School Path...")
    # Vsj dims: (j, R, t_own, degree, i_a_school) = 4 non-j
    n_school_periods = 3
    Vsj  = zeros(Float32, 2, 3, 2, apnts_c)
    PFsj = copy(Vsj)

    Vs  = zeros(Float32, n_school_periods, 2, 3, 2, apnts_c)
    PFs = copy(Vs)

    # j=4: last period for 4yr — graduate into Vc1[5]
    Vc1_j3 = @view V_c1[3, :,:,:,:,:,:,:,:,:]
    Vs[3, :,:,:,:], PFs[3, :,:,:,:] = VSj_enrolled(nothing, Vc1_j3, Vsj, PFsj, model, 4)

    # j=3: 4yr still enrolled
    Vs[2, :,:,:,:], PFs[2, :,:,:,:] = VSj_enrolled(@view(Vs[3, :,:,:,:]),nothing, Vsj, PFsj,model, 3)

    # j=2: 2yr graduates into Vc1[3]; 4yr still enrolled
    Vc1_j1 = @view Vj_c1[1, :,:,:,:,:,:,:,:,:]
    Vs[1, :,:,:,:], PFs[1, :,:,:,:] = VSj_enrolled(@view(Vs[2, :,:,:,:]),Vc1_j1,Vsj, PFsj,model, 2)

    # j=1: first period
    # Vsj_1 dims: (R, m_p, t_parent, e_p, i_a_p, i_a, i_z_p, edu_shock, degree)
    Vs_1  = zeros(Float32, 2, 2, 3, 3, apnts_nc, apnts_nc, zpnts, 2, 2)
    PFs_1 = copy(Vs_1)

    Vs_1, PFs_1 = VSj_first_period(@view(Vs[1, :,:,:,:]), Vs_1, PFs_1, model)

    open(joinpath(savepath, "School.jls"), "w") do io
            serialize(io, (V_c1  = V_c1,  PF_c1  = PF_c1,
        V_c2  = V_c2,  PF_c2  = PF_c2,
        W_c   = W_c,   WPF_c  = WPF_c,
        Vs    = Vs,    PFs    = PFs,
        Vs_1  = Vs_1,  PFs_1  = PFs_1))#
        end

    return (
        V_nc1 = V_nc1, PF_nc1 = PF_nc1,
        V_nc2 = V_nc2, PF_nc2 = PF_nc2,
        W_nc  = W_nc,  WPF_nc = WPF_nc,
        V_c1  = V_c1,  PF_c1  = PF_c1,
        V_c2  = V_c2,  PF_c2  = PF_c2,
        W_c   = W_c,   WPF_c  = WPF_c,
        Vs    = Vs,    PFs    = PFs,
        Vs_1  = Vs_1,  PFs_1  = PFs_1
    )
end