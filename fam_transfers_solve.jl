function solve_model(model)
    (; apnts_nc, apnts_c, zpnts) = model

    # -----------------------------------------------------------------------
    # Dimensions:
    # R: 2 (white, black)
    # m: 2 (single, married)  
    # n: 5 (family size 1-5)
    # t: 3 (low, mid, high)
    # e: 3 (no college=1, 2yr=2, 4yr=3) — only matters post-graduation
    # shock_in, shock_out, past_in, past_out: 2 each
    # -----------------------------------------------------------------------

    # -----------------------------------------------------------------------
    # NO COLLEGE PATH
    # -----------------------------------------------------------------------

    # --- Retirement (ages 61-84, j=44 to j=67, 24 periods) ---
    Wj_nc  = zeros(Float32, 24, 2, 2, 5, 3, apnts_nc, zpnts, 2, 2, 2, 2)
    WPFj_nc = copy(Wj_nc)

    # Terminal period j=67 (age 84)
    Wj_nc[24,:,:,:,:,:,:,:,:,:,:], WPFj_nc[24,:,:,:,:,:,:,:,:,:,:] =
        Wncj(nothing, model, 67)

    # Backward through retirement
    for j in 66:-1:44
        idx = j - 43   # maps j=44->1, j=66->23, terminal already at 24
        W_jp1 = @view Wj_nc[idx+1, :,:,:,:,:,:,:,:,:,:]
        Wj_nc[idx,:,:,:,:,:,:,:,:,:,:], WPFj_nc[idx,:,:,:,:,:,:,:,:,:,:] =
            Wncj(W_jp1, model, j)
    end

    # --- Working no college (ages 18-60, j=1 to j=43) ---
    Vj_nc  = zeros(Float32, 43, 2, 2, 5, 3, apnts_nc, zpnts, 2, 2, 2, 2)
    PFj_nc = copy(Vj_nc)

    # Last working period j=43 (age 60) — continuation is first retirement period
    W_ret_first = @view Wj_nc[1, :,:,:,:,:,:,:,:,:,:]   # j=44, idx=1
    Vj_nc[43,:,:,:,:,:,:,:,:,:,:], PFj_nc[43,:,:,:,:,:,:,:,:,:,:] =
        Vncj_solve(W_ret_first, nothing, model, 43)

    # Backward through working life
    for j in 42:-1:1
        V_jp1 = @view Vj_nc[j+1, :,:,:,:,:,:,:,:,:,:]
        Wj_nc_first = @view Wj_nc[1, :,:,:,:,:,:,:,:,:,:]   # passed but not used if j < 43
        Vj_nc[j,:,:,:,:,:,:,:,:,:,:], PFj_nc[j,:,:,:,:,:,:,:,:,:,:] =
            Vncj_solve(V_jp1, Wj_nc_first, model, j)
    end

    # -----------------------------------------------------------------------
    # COLLEGE PATH
    # -----------------------------------------------------------------------

    # --- Retirement (same grid structure, college asset grid) ---
    Wj_c  = zeros(Float32, 24, 2, 2, 5, 3, 2, apnts_c, zpnts, 2, 2, 2, 2)
    WPFj_c = copy(Wj_c)

    # e dimension: 2 (2yr=1, 4yr=2) — only two college types
    Wj_c[24,:,:,:,:,:,:,:,:,:,:,:], WPFj_c[24,:,:,:,:,:,:,:,:,:,:,:] =
        Wcj(nothing, model, 67)

    for j in 66:-1:44
        idx = j - 43
        W_jp1 = @view Wj_c[idx+1, :,:,:,:,:,:,:,:,:,:,:]
        Wj_c[idx,:,:,:,:,:,:,:,:,:,:,:], WPFj_c[idx,:,:,:,:,:,:,:,:,:,:,:] =
            Wcj(W_jp1, model, j)
    end

    # --- Working college (ages 22-60 for 4yr, 20-60 for 2yr, j=3 or j=5 to j=43) ---
    # Solve for all j=1 to j=43 — school functions will handle which j is valid
    Vj_c  = zeros(Float32, 43, 2, 2, 5, 3, 2, apnts_c, zpnts, 2, 2, 2, 2)
    PFj_c = copy(Vj_c)

    # Last working period
    W_ret_first_c = @view Wj_c[1, :,:,:,:,:,:,:,:,:,:,:]
    Vj_c[43,:,:,:,:,:,:,:,:,:,:,:], PFj_c[43,:,:,:,:,:,:,:,:,:,:,:] =
        Vcj_solve(W_ret_first_c, nothing, model, 43)

    for j in 42:-1:1
        V_jp1 = @view Vj_c[j+1, :,:,:,:,:,:,:,:,:,:,:]
        Vj_c[j,:,:,:,:,:,:,:,:,:,:,:], PFj_c[j,:,:,:,:,:,:,:,:,:,:,:] =
            Vcj_solve(V_jp1, W_ret_first_c, model, j)
    end

    # -----------------------------------------------------------------------
    # SCHOOL PATH
    # j=1: age 18 — first period (education decision made here)
    # j=2: age 19 — 2yr enrolled OR 2yr just graduated (enters Vj_c at j=2)
    # j=3: age 20 — 4yr enrolled year 2
    # j=4: age 21 — 4yr enrolled year 3  
    # j=5: age 22 — 4yr enrolled year 4 / graduation (enters Vj_c at j=5)
    # -----------------------------------------------------------------------

    # School state space: (R, t_parent, degree, i_a)
    # t_parent: 3 (low, mid, high) — parental background during school
    # degree: 2 (2yr=1, 4yr=2)
    # No (m, n, shock) dims during school — agent not yet household head

    # Enrolled periods — solve backward from graduation
    # 2yr: one enrolled period j=2, graduates into Vj_c[2,:] at j=3
    # 4yr: three enrolled periods j=2,3,4, graduates into Vj_c[5,:] at j=5

    # Arrays: (R, t_parent, degree, i_a)
    n_school_periods = 4   # j=1 to j=4 (age 18-21), j=5 is graduation for 4yr
    Vsj  = zeros(Float32, n_school_periods, 2, 3, 2, apnts_c)
    PFsj = copy(Vsj)

    # j=4 (age 21): 4yr final enrolled year — continuation is Vj_c[5] (working at 22)
    # j=2 (age 19): 2yr final enrolled year — continuation is Vj_c[3] (working at 20)
    # Graduation transitions: 2yr exits at j=3, 4yr exits at j=5

    # Solve last enrolled period for 4yr (j=4, graduating into j=5 working)
    V_grad_4yr = @view Vj_c[5, :,:,:,:,:,:,:,:,:,:,:]   # working at age 22
    Vsj[4,:,:,:,:], PFsj[4,:,:,:,:] = VSj_enrolled(V_grad_4yr, nothing, model, 4)

    # j=3 for 4yr (still enrolled), j=2 is graduation year for 2yr
    V_grad_2yr = @view Vj_c[3, :,:,:,:,:,:,:,:,:,:,:]   # working at age 20
    Vsj[3,:,:,:,:], PFsj[3,:,:,:,:] = VSj_enrolled(V_grad_2yr, @view(Vsj[4,:,:,:,:]), model, 3)

    # j=2: 2yr graduation year, 4yr second enrolled year
    Vsj[2,:,:,:,:], PFsj[2,:,:,:,:] = VSj_enrolled(V_grad_2yr, @view(Vsj[3,:,:,:,:]), model, 2)

    # j=1: first period — education decision
    # Vsj[2] is the value of being enrolled next period
    # Vj_nc[1] is the value of no college at j=1
    Vsj_1  = zeros(Float32, 2, 2, 5, 3, 3, apnts_nc, zpnts, 2, 2)
    PFsj_1 = copy(Vsj_1)

    Vsj_1, PFsj_1 = VSj_first_period(
        @view(Vsj[2,:,:,:,:]),    # continuation if enrolled
        Vsj_1, PFsj_1, model
    )

    # -----------------------------------------------------------------------
    # Return everything needed for simulation and education decision
    # -----------------------------------------------------------------------
    return (
        # No college
        Vj_nc  = Vj_nc,
        PFj_nc = PFj_nc,
        Wj_nc  = Wj_nc,
        WPFj_nc = WPFj_nc,
        # College working
        Vj_c   = Vj_c,
        PFj_c  = PFj_c,
        Wj_c   = Wj_c,
        WPFj_c = WPFj_c,
        # School
        Vsj    = Vsj,
        PFsj   = PFsj,
        Vsj_1  = Vsj_1,
        PFsj_1 = PFsj_1,
    )
end