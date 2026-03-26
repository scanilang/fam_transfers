###############################################################################################
# Working
###############################################################################################
# OECD scale
function oecd(m,n)
    # oecd equivalence scale (1 + 0.7(Adults - 1) + 0.5 kids)
    # m = `1 is married
    # n family size`
    return (1 + 0.7(m) + 0.5(n-(m + 1)))
end

# CRRA utility function
function u(c,gamma)
    Val = ((c)^(1-gamma)) / (1-gamma)
    #if Val == -Inf
    #    Val = -100000
    #end
    return Val
end

# CRRA utility function scaled
function u_hh(c,gamma, m,n)
    Val = ((c/oecd(m,n))^(1-gamma)) / (1-gamma)
    return n * Val
end

###############################################################################################
# Income
###############################################################################################


function g(r,j,e, n)
    mid_age = (22.5 + 4*j)

    if r == 2
        g = 4.469 + 0.4817*mid_age - 0.01037*mid_age^2 + 0.00007375*mid_age^3 + 0.06854*n -0.08079
    else
        g = 9.509 +  0.1356*mid_age - 0.001852*(mid_age^2)  + 0.000005373*(mid_age^3) + 0.07886*n + 0.04664
    end
    return g
end

function tax_y(y, m)
    div_y = y/2
    
    if m == 1
        lambda = 0.08
        tau_y = 0.115
    else 
        lambda = 0.07
        tau_y = 0.07
    end
    peryear = div_y - (1-lambda)*div_y^(1-tau_y)
    return peryear * 2
end

###############################################################################################
# Transfer Functions
###############################################################################################
using DataFrames
if pwd() == "/Users/scanilang/Documents/econ/umn/family_transfers/2026"
    transfer_probit = CSV.read("/Users/scanilang/Documents/econ/umn/family_transfers/data/transfer_probit_results.csv", DataFrame)
    transfer_amount = CSV.read("/Users/scanilang/Documents/econ/umn/family_transfers/data/transfer_amount_results.csv", DataFrame)
else
    transfer_probit = CSV.read("/users/4/canil007/bankruptcy/family_transfers/Data/transfer_probit_results.csv", DataFrame)
    transfer_amount = CSV.read("/users/4/canil007/bankruptcy/family_transfers/Data/transfer_amount_results.csv", DataFrame)
end

β_white_probit_in   = Tuple(transfer_probit.white_probit_in)
β_black_probit_in   = Tuple(transfer_probit.black_probit_in)
β_white_probit_out   = Tuple(transfer_probit.white_probit_out)
β_black_probit_out   = Tuple(transfer_probit.black_probit_out)
β_white_transfer_in  = Tuple(transfer_amount.white_transfer_in)
β_black_transfer_in = Tuple(transfer_amount.black_transfer_in)
β_black_transfer_out  = Tuple(transfer_amount.black_transfer_out)
β_white_transfer_out     = Tuple(transfer_amount.white_transfer_out)


function shocks_out_prob(β_white_probit_out, β_black_probit_out, n, m, j, y, a_income, e, t, past_in, past_out)
    age = j + 17
    e_1 = e == 1 ? 1 : 0
    e_0 = e == 0 ? 1 : 0
    f_1 = t == 1 ? 1 : 0
    f_2 = t == 2 ? 1 : 0
    f_3 = t == 3 ? 1 : 0

    if r == 2 
        val = β_black_probit_out[1] + β_black_probit_out[2]*log(y) + β_black_probit_out[3]*log(a_income +1 ) + β_black_probit_out[4]*age +  
        β_black_probit_out[5]*age^2 + β_black_probit_out[6]*n + β_black_probit_out[7]*e_0 + β_black_probit_out[8]*e_1 + β_black_probit_out[9]*m
        β_black_probit_out[10]*f_1 + β_black_probit_out[11]*f_2 + β_black_probit_out[12]*f_3 + β_black_probit_out[13]*past_in + β_black_probit_out[14]*past_out + β_white_probit_out[20]
    else
        val = β_white_probit_out[1] + β_white_probit_out[2]*log(y) + β_white_probit_out[3]*log(a_income +1 ) + β_white_probit_out[4]*age + 
        β_white_probit_out[5]*age^2 + β_white_probit_out[6]*n + β_white_probit_out[7]*e_0 + β_white_probit_out[8]*e_1 + β_white_probit_out[9]*m
        β_white_probit_out[10]*f_1 + β_white_probit_out[11]*f_2 + β_white_probit_out[12]*f_3 + β_white_probit_out[13]*past_in + β_white_probit_out[14]*past_out + β_white_probit_out[20]
    end

    return cdf(Normal(), val)
end

function transfers_out_amount(β_white_transfer_out, β_black_transfer_out,r,n,m, j, y, a_income, e, t)
    age = j + 17
    e_1 = e == 1 ? 1 : 0
    e_0 = e == 0 ? 1 : 0
    f_1 = t == 1 ? 1 : 0
    f_2 = t == 2 ? 1 : 0
    f_3 = t == 3 ? 1 : 0

    if r == 2 
        val = β_black_transfer_out[1] + β_black_transfer_out[2]*log(y) + β_black_transfer_out[3]*log(a_income +1 ) + β_black_transfer_out[4]*age + 
        β_black_transfer_out[5]*age^2 + β_black_transfer_out[6]*n + β_black_transfer_out[7]*e_0 + β_black_transfer_out[8]*e_1 + β_black_transfer_out[9]*m
        β_black_transfer_out[10]*f_1 + β_black_transfer_out[11]*f_2 + β_black_transfer_out[12]*f_3 + β_black_transfer_out[14]
    else
        val = β_white_transfer_out[1] + β_white_transfer_out[2]*log(y) + β_white_transfer_out[3]*log(a_income +1 ) + β_white_transfer_out[4]*age + 
        β_white_transfer_out[5]*age^2 + β_white_transfer_out[6]*n + β_white_transfer_out[7]*e_0 + β_white_transfer_out[8]*e_1 + β_white_transfer_out[9]*m
        β_white_transfer_out[10]*f_1 + β_white_transfer_out[11]*f_2 + β_white_transfer_out[12]*f_3 + β_white_transfer_out[14]
    end
    return exp(val)
end

function shocks_in_prob(β_white_probit_in, β_black_probit_in,r,n,m,j,y, a_income, e, t, past_in, past_out)
    age = j + 17
    e_1 = e == 1 ? 1 : 0
    e_0 = e == 0 ? 1 : 0
    f_1 = t == 1 ? 1 : 0
    f_2 = t == 2 ? 1 : 0
    f_3 = t == 3 ? 1 : 0

    if r == 2 
        val = β_black_probit_in[1] + β_black_probit_in[2]*log(y) + β_black_probit_in[3]*log(a_income) + β_black_probit_in[4]*age +  
        β_black_probit_in[5]*age^2 + β_black_probit_in[6]*n + β_black_probit_in[7]*e_0 + β_black_probit_in[8]*e_1 + β_black_probit_in[9]*m
        β_black_probit_in[10]*f_1 + β_black_probit_in[11]*f_2 + β_black_probit_in[12]*f_3 + β_black_probit_in[13]*past_in + β_black_probit_in[14]*past_out + β_white_probit_in[20]
    else
        val = β_white_probit_in[1] + β_white_probit_in[2]*log(y) + β_white_probit_in[3]*log(a_income) + β_white_probit_in[4]*age + 
        β_white_probit_in[5]*age^2 + β_white_probit_in[6]*n + β_white_probit_in[7]*e_0 + β_white_probit_in[8]*e_1 + β_white_probit_in[9]*m
        β_white_probit_in[10]*f_1 + β_white_probit_in[11]*f_2 + β_white_probit_in[12]*f_3 + β_white_probit_in[13]*past_in + β_white_probit_in[14]*past_out + β_white_probit_in[20]
    end

    return cdf(Normal(), val)
end

function transfers_in_amount(β_white_transfer_in, β_black_transfer_in,r,n,m,j,y, a_income, e, t)
    age = j + 17
    e_1 = e == 1 ? 1 : 0
    e_0 = e == 0 ? 1 : 0
    f_1 = t == 1 ? 1 : 0
    f_2 = t == 2 ? 1 : 0
    f_3 = t == 3 ? 1 : 0

    if r == 2 
        val = β_black_transfer_in[1] + β_black_transfer_in[2]*log(y) + β_black_transfer_in[3]*log(a_income +1 ) + β_black_transfer_in[4]*age + 
        β_black_transfer_in[5]*age^2 + β_black_transfer_in[6]*n + β_black_transfer_in[7]*e_0 + β_black_transfer_in[8]*e_1 + β_black_transfer_in[9]*m
        β_black_transfer_in[10]*f_1 + β_black_transfer_in[11]*f_2 + β_black_transfer_in[12]*f_3 + β_black_transfer_in[14]
    else
        val = β_white_transfer_in[1] + β_white_transfer_in[2]*log(y) + β_white_transfer_in[3]*log(a_income +1 ) + β_white_transfer_in[4]*age + 
        β_white_transfer_in[5]*age^2 + β_white_transfer_in[6]*n + β_white_transfer_in[7]*e_0 + β_white_transfer_in[8]*e_1 + β_white_transfer_in[9]*m
        β_white_transfer_in[10]*f_1 + β_white_transfer_in[11]*f_2 + β_white_transfer_in[12]*f_3 + β_white_transfer_in[14]
    end
    return exp(val)
end

###############################################################################################
# Survival Probabilities
###############################################################################################
# survival probabilities
if pwd() == "/Users/scanilang/Documents/econ/umn/bankruptcy/model/backwards induction"
    survival_risk_df = CSV.read("/Users/scanilang/Documents/econ/umn/bankruptcy/model/Death_rate_by_age_race.csv", DataFrame; limit = 3)
else
    survival_risk_df = CSV.read("/users/4/canil007/bankruptcy/bankruptcy/Data/Death_rate_by_age_race.csv", DataFrame; limit = 3)
end

survival_risk_df[(1:3),(2:3)] = survival_risk_df[(1:3),(2:3)] ./ 100000 # rate per 100,000
survival_risk_df[(1:3),(2:3)] = 1 .- (1 .-survival_risk_df[(1:3),(2:3)]).^(1/4) # convert to quarterly death rates
survival_risk = Matrix{Float64}(survival_risk_df[:, 2:3]) 


survival_risk_full = ones(11, 3)
survival_risk_full[121:140, 2] .= survival_risk[1, 1]
survival_risk_full[121:140, 3] .= survival_risk[1, 2]
survival_risk_full[141:180, 2] .= survival_risk[2, 1]
survival_risk_full[141:180, 3] .= survival_risk[2, 2]
survival_risk_full[181:200, 2] .= survival_risk[3, 1]
survival_risk_full[181:200, 3] .= survival_risk[3, 2]

survival_risk_full[121:200, 2:3] = 1 .- survival_risk_full[121:200, 2:3]

survival_risk_full = survival_risk_full[:, 2:3]