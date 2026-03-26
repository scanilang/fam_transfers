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

function g(r,j,n)
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


if pwd() == "/Users/scanilang/Documents/econ/umn/family_transfers/2026"
    transfer_probit = CSV.read("/Users/scanilang/Documents/econ/umn/family_transfers/data/transfer_probit.csv", DataFrame)
    transfer_amount = CSV.read("/Users/scanilang/Documents/econ/umn/family_transfers/data/transfer_amount.csv", DataFrame)
else
    transfer_probit = CSV.read("/users/4/canil007/bankruptcy/family_transfers/Data/transfer_probit.csv", DataFrame)
    transfer_amount = CSV.read("/users/4/canil007/bankruptcy/family_transfers/Data/transfer_amount.csv", DataFrame)
end

β_white_probit_in   = Tuple(transfer_probit.white_probit_in)
β_black_probit_in   = Tuple(transfer_probit.black_probit_in)
β_white_probit_out   = Tuple(transfer_probit.white_probit_out)
β_black_probit_out   = Tuple(transfer_probit.black_probit_out)
β_white_transfer_in  = Tuple(transfer_amount.white_transfer_in)
β_black_transfer_in = Tuple(transfer_amount.black_transfer_in)
β_black_transfer_out  = Tuple(transfer_amount.black_transfer_out)
β_white_transfer_out     = Tuple(transfer_amount.white_transfer_out)

struct ModelCoefs
    β_white_probit_in::NTuple{19,Float64}
    β_black_probit_in::NTuple{19,Float64}
    β_white_probit_out::NTuple{19,Float64}
    β_black_probit_out::NTuple{19,Float64}
    β_white_transfer_in::NTuple{13,Float64}
    β_black_transfer_in::NTuple{13,Float64}
    β_black_transfer_out::NTuple{13,Float64}
    β_white_transfer_out::NTuple{13,Float64}
end

coefs = ModelCoefs(
    β_white_probit_in,
    β_black_probit_in,
    β_white_probit_out,
    β_black_probit_out,
    β_white_transfer_in,
    β_black_transfer_in,
    β_black_transfer_out,
    β_white_transfer_out
)


function shocks_out_prob(r,n,m, j, y, past_in, past_out)
    mid_age = 22.5 + 4*j
    if r == 2 
        val = - 2.25826 - 0.11377 * n - 0.04575 *m + (0.01669455 *mid_age)  - (0.0005308654 *mid_age^2) - past_in*0.50769 + past_out* 1.45337 +log(y) *0.32320 -0.53053 - 1.87701 # transfer min, income max 180000
    else
        val = -4.463962 - 0.09363 *n + 0.18690*m + (0.0397062151 *mid_age)  - (0.0003243744 *mid_age^2) -past_in*0.08483 + past_out*1.73401 + log(y) * 0.24322 - 1.71271 + 1.04578 # # transfer min, income max 180000

    end

    return cdf(Normal(), val)
end

function transfers_out_amount(r,n,m, j, y)
    mid_age = 22.5 + 4*j

    if r == 2
        return exp(-0.9451237 - 0.12296*n - 0.23850*m + (0.0419915484 *mid_age) - (0.0003287593 *mid_age^2) +log(y) * 0.62658 + 0.19965 + 0.17305) # transfer of at least 200, income max 180000
    else
        return exp(-3.00629 -0.134370 *n - 0.054113*m + (0.0441335028 *mid_age)  - (0.0003151821 *mid_age^2) +log(y) * 0.788086 + 0.714118 + 0.212124 ) # transfer of at least 200, income max 180000
    end

end

function shocks_in_prob(r,n,m,j,y, past_in, past_out)
    mid_age = 22.5 + 4*j

    if r == 2
      val = 1.458929 - 0.09941*n - 0.15475*m  - (0.1211193934 *mid_age)  + (0.0007720176 *mid_age^2) +past_in* 1.46933 - past_out*0.06697 - log(y) * 0.18051 + 2.96273 - 0.27434 # transfer of at least 200
    else
      val = 3.600208 -0.04262*n + 0.06188*m - (0.1096650013 *mid_age) + (0.0006706238 *mid_age^2) + past_in*1.58966 - past_out*0.32242 - log(y) * 0.30893 + 1.82198 + 0.07752 # transfer of at least 200, income max 180000
    end
    return cdf(Normal(), val)
end

function transfers_in_amount(r,n,m,j,y)
    mid_age = 22.5 + 4*j

    if r == 2
        return exp(3.128561 - 0.06088*n - 0.17337*m + (0.098899763 *mid_age)  - (0.001003151*mid_age^2) + log(y) * 0.11010 + 0.89794 -0.13826) # transfers of at least 200
    else
        return exp(3.041345 - 0.07338*n + 0.11436 *m + (0.0518311507 *mid_age)  - (0.0005127471  *mid_age^2) + log(y) *0.30850 +0.31655 + 0.44215) # transfer of at least 200, income max 180000
    end
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