function solve_model(model)
    (; Race, marital_status, fam_size, ed_type, apnts, zpnts, tasks_idx) = model

    # Solve Retirement
    Wj = zeros(Float64, 11, 2, 2, 6, 3, apnts, zpnts, 2, 2, 2, 2) # (j, R, m, n, e, i_a, i_z, shock_in, shock_out, past_in, past_out)
    PWj = copy(Wj)
    Wj[11, :, :, :, :, :, :, :, :, :, :], PWj[11, :, :, :, :, :, :, :, :, :, :]  = Wj_solve(nothing, model, 34)

    for j in 33:-1:23
        W_jp1 = @view W2[j+1, :, :, :, :, :, :, :, :, :, :];
        Wj[j, :, :, :, :, :, :, :, :, :, :], PWj[j, :, :, :, :, :, :, :, :, :, :] = Wj(W_jp1, model, j)
    end

    # Solve Working
    Vj = zeros(Float64, 11, 2, 2, 6, 3, apnts, zpnts, 2, 2, 2, 2)
    PFj = copy(Vj)
    Vj[23, :, :, :, :, :, :, :, :, :, :], PFj[23, :, :, :, :, :, :, :, :, :, :] = Vj_solve(Wj[23, :, :, :, :, :, :, :, :, :, :], model, 23)

    for j in 22:-1:1
        V_jp1 = @view Vj[j+1, :, :, :, :, :, :, :, :, :, :];
        Vj[j, :, :, :, :, :, :, :, :, :, :], PFj[j, :, :, :, :, :, :, :, :, :, :] = Vj(V_jp1, model, j)
    end

    return (; Vj, PFj, Wj, PWj)
end