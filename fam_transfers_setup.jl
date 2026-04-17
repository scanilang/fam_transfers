using DataFrames

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
if pwd() == "/Users/scanilang/Documents/econ/umn/family_transfers/2026"
    income_results = CSV.read("/Users/scanilang/Documents/econ/umn/family_transfers/data/income_results.csv", DataFrame)
else

    income_results = CSV.read("/users/4/canil007/bankruptcy/family_transfers/Data/income_results.csv", DataFrame)
end

β_black_income   = Tuple(income_results.black_income_process)
β_white_income   = Tuple(income_results.white_income_process)

function g( r,j,e, m)
    age = j + 17
    e_2 = e == 2 ? 1 : 0
    e_0 = e == 1 ? 1 : 0
    m = m - 1  # married=2 becomes 1, single=1 becomes 0

    if r == 2
        g = β_black_income[1] + β_black_income[2]*age + β_black_income[3]*age^2 + β_black_income[4]*age^3 + β_black_income[5]*e_0 + β_black_income[6]*e_2 + β_black_income[7]*m + β_black_income[26]
    else
        g = β_white_income[1] + β_white_income[2]*age + β_white_income[3]*age^2 + β_white_income[4]*age^3 + β_white_income[5]*e_0 + β_white_income[6]*e_2 + β_white_income[7]*m + β_white_income[26]
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
# Get Regression Coefficients
###############################################################################################

if pwd() == "/Users/scanilang/Documents/econ/umn/family_transfers/2026"
    transfer_probit = CSV.read("/Users/scanilang/Documents/econ/umn/family_transfers/data/transfer_probit_results.csv", DataFrame)
    transfer_amount = CSV.read("/Users/scanilang/Documents/econ/umn/family_transfers/data/transfer_amount_results.csv", DataFrame)
    edu_transfer_probit = CSV.read("/Users/scanilang/Documents/econ/umn/family_transfers/data/edu_transfer_probit_results.csv", DataFrame)
    edu_transfer_amount = CSV.read("/Users/scanilang/Documents/econ/umn/family_transfers/data/edu_transfer_amount_results.csv", DataFrame)
else
    transfer_probit = CSV.read("/users/4/canil007/bankruptcy/family_transfers/Data/transfer_probit_results.csv", DataFrame)
    transfer_amount = CSV.read("/users/4/canil007/bankruptcy/family_transfers/Data/transfer_amount_results.csv", DataFrame)
    edu_transfer_probit = CSV.read("/users/4/canil007/bankruptcy/family_transfers/Data/edu_transfer_probit_results.csv", DataFrame)
    edu_transfer_amount = CSV.read("/users/4/canil007/bankruptcy/family_transfers/Data/edu_transfer_amount_results.csv", DataFrame)
end

###############################################################################################
# Education Transfer Functions
###############################################################################################

β_edu_probit = Tuple(edu_transfer_probit.pooled_probit_edu)
β_edu_amount = Tuple(edu_transfer_amount.pooled_transfer_edu)

function edu_transfer_prob(R, m, n, y, a_income, e, t, degree_choice)
    race_white = R == 1 ? 1 : 0
    is_4yr     = degree_choice == 2 ? 0 : 1
    e_0 = e == 1 ? 1 : 0   # No College
    e_2 = e == 2 ? 1 : 0   # Some College
    m_single = m == 1 ? 1 : 0
    f_1 = t == 1 ? 1 : 0   # both_low
    f_2 = t == 2 ? 1 : 0   # both_mid

    val = β_edu_probit[1] +                    # intercept
          β_edu_probit[2] * log(y + 1) +        # log_nonasset_income
          β_edu_probit[3] * log(a_income + 1) + # log_asset_income
          β_edu_probit[4] * is_4yr +            # degree_type_final4yr
          β_edu_probit[5] * n +                 # Family_Unit_Size
          β_edu_probit[6] * e_0 +               # Head_CollegeNo College
          β_edu_probit[7] * e_2 +               # Head_CollegeSome College
          β_edu_probit[8] * m_single +           # Marital_StatusSingle
          β_edu_probit[9] * f_1 +               # family_typeboth_low
          β_edu_probit[10] * f_2 +              # family_typeboth_mid
          β_edu_probit[11] * race_white          # Race_HeadWhite
          # [13] and [14] enroll_era — omitted, set to reference category

    return cdf(Normal(), val)
end

# OLS: E(log amount | parents help, attending college)
# Coefficient order from R:
# [1] (Intercept)              [2] log_nonasset_income     [3] log_asset_income
# [4] Head_CollegeNo College   [5] Head_CollegeSome College [6] degree_type_final4yr
# [7] enroll_era2009-2012      [8] enroll_erapre-2002       [9] Race_HeadWhite

function edu_transfer(R, y, a_income, e, degree_choice)
    race_white = R == 1 ? 1 : 0
    is_4yr     = degree_choice == 2 ? 0 : 1
    e_0 = e == 1 ? 1 : 0   # No College
    e_2 = e == 2 ? 1 : 0   # Some College

    val = β_edu_amount[1] +                    # intercept
          β_edu_amount[2] * log(y + 1) +        # log_nonasset_income
          β_edu_amount[3] * log(a_income + 1) + # log_asset_income
          β_edu_amount[4] * e_0 +               # Head_CollegeNo College
          β_edu_amount[5] * e_2 +               # Head_CollegeSome College
          β_edu_amount[6] * is_4yr +            # degree_type_final4yr
          β_edu_amount[7] * race_white           # Race_HeadWhite
          # [7] and [8] enroll_era — omitted, set to reference category
          # NOTE: shift indices if you keep enroll_era in export

    return exp(val)
end

# Expected education transfer: P(help) × E(amount | help)
function expected_edu_transfer(r, n, m, y, a_income, e, t, degree_choice)
    prob   = edu_transfer_prob(r, n, m, y, a_income, e, t, degree_choice)
    amount = edu_transfer(r, y, a_income, e, degree_choice)
    return prob * amount
end

###############################################################################################
# Intervivos Transfer Functions
###############################################################################################

β_white_probit_in   = Tuple(transfer_probit.white_probit_in)
β_black_probit_in   = Tuple(transfer_probit.black_probit_in)
β_white_probit_out   = Tuple(transfer_probit.white_probit_out)
β_black_probit_out   = Tuple(transfer_probit.black_probit_out)
β_white_transfer_in  = Tuple(transfer_amount.white_transfer_in)
β_black_transfer_in = Tuple(transfer_amount.black_transfer_in)
β_black_transfer_out  = Tuple(transfer_amount.black_transfer_out)
β_white_transfer_out     = Tuple(transfer_amount.white_transfer_out)

β_white_edu_probit_in   = Tuple(edu_transfer_probit.white_edu_probit_in)
β_black_edu_probit_in   = Tuple(edu_transfer_probit.black_edu_probit_in)
β_white_edu_transfer_in  = Tuple(edu_transfer_amount.white_edu_transfer_in)
β_black_edu_transfer_in = Tuple(edu_transfer_amount.black_edu_transfer_in)
##################### Transfer out
function shocks_out_prob(r, n, m, j, y, a_income, e, t, past_in, past_out)
    age = j + 17
    e_2 = e == 2 ? 1 : 0
    e_0 = e == 1 ? 1 : 0
    f_1 = t == 1 ? 1 : 0
    f_2 = t == 2 ? 1 : 0
    m = m - 1  # married=2 becomes 1, single=1 becomes 0
    if r == 2 
        val = β_black_probit_out[1] + β_black_probit_out[2]*log(y) + β_black_probit_out[3]*log(a_income +1 ) + β_black_probit_out[4]*age +  
        β_black_probit_out[5]*age^2 + β_black_probit_out[6]*n + β_black_probit_out[7]*e_0 + β_black_probit_out[8]*e_2 + β_black_probit_out[9]*m + 
        β_black_probit_out[10]*f_1 + β_black_probit_out[11]*f_2  + β_black_probit_out[12]*past_in + β_black_probit_out[13]*past_out + β_white_probit_out[19]
    else
        val = β_white_probit_out[1] + β_white_probit_out[2]*log(y) + β_white_probit_out[3]*log(a_income +1 ) + β_white_probit_out[4]*age + 
        β_white_probit_out[5]*age^2 + β_white_probit_out[6]*n + β_white_probit_out[7]*e_0 + β_white_probit_out[8]*e_2 + β_white_probit_out[9]*m + 
        β_white_probit_out[10]*f_1 + β_white_probit_out[11]*f_2  + β_white_probit_out[12]*past_in + β_white_probit_out[13]*past_out + β_white_probit_out[19]
    end

    return cdf(Normal(), val)
end

function transfers_out_amount(r,n,m, j, y, a_income, e, t)
    age = j + 17
    e_2 = e == 2 ? 1 : 0
    e_0 = e == 1 ? 1 : 0
    f_1 = t == 1 ? 1 : 0
    f_2 = t == 2 ? 1 : 0
    m = m - 1  # married=2 becomes 1, single=1 becomes 0

    if r == 2 
        val = β_black_transfer_out[1] + β_black_transfer_out[2]*log(y) + β_black_transfer_out[3]*log(a_income +1 ) + β_black_transfer_out[4]*age + 
        β_black_transfer_out[5]*age^2 + β_black_transfer_out[6]*n + β_black_transfer_out[7]*e_0 + β_black_transfer_out[8]*e_2 + β_black_transfer_out[9]*m + 
        β_black_transfer_out[10]*f_1 + β_black_transfer_out[11]*f_2  + β_black_transfer_out[13]
    else
        val = β_white_transfer_out[1] + β_white_transfer_out[2]*log(y) + β_white_transfer_out[3]*log(a_income +1 ) + β_white_transfer_out[4]*age + 
        β_white_transfer_out[5]*age^2 + β_white_transfer_out[6]*n + β_white_transfer_out[7]*e_0 + β_white_transfer_out[8]*e_2 + β_white_transfer_out[9]*m + 
        β_white_transfer_out[10]*f_1 + β_white_transfer_out[11]*f_2 + β_white_transfer_out[13]
    end
    return exp(val)
end

##################### Transfer in
function shocks_in_prob(r,n,m,j,y, a_income, e, t, past_in, past_out)
    age = j + 17
    e_2 = e == 2 ? 1 : 0
    e_0 = e == 1 ? 1 : 0
    f_1 = t == 1 ? 1 : 0
    f_2 = t == 2 ? 1 : 0
    m = m - 1  # married=2 becomes 1, single=1 becomes 0

    if r == 2 
        val = β_black_probit_in[1] + β_black_probit_in[2]*log(y) + β_black_probit_in[3]*log(a_income) + β_black_probit_in[4]*age +  
        β_black_probit_in[5]*age^2 + β_black_probit_in[6]*n + β_black_probit_in[7]*e_0 + β_black_probit_in[8]*e_2 + β_black_probit_in[9]*m + 
        β_black_probit_in[10]*f_1 + β_black_probit_in[11]*f_2  + β_black_probit_in[12]*past_in + β_black_probit_in[13]*past_out + β_white_probit_in[19]
    else
        val = β_white_probit_in[1] + β_white_probit_in[2]*log(y) + β_white_probit_in[3]*log(a_income) + β_white_probit_in[4]*age + 
        β_white_probit_in[5]*age^2 + β_white_probit_in[6]*n + β_white_probit_in[7]*e_0 + β_white_probit_in[8]*e_2 + β_white_probit_in[9]*m + 
        β_white_probit_in[10]*f_1 + β_white_probit_in[11]*f_2 + β_white_probit_in[12]*past_in + β_white_probit_in[13]*past_out + β_white_probit_in[19]
    end

    return cdf(Normal(), val)
end

function transfers_in_amount(r,n,m,j,y, a_income, e, t)
    age = j + 17
    e_2 = e == 2 ? 1 : 0
    e_0 = e == 1 ? 1 : 0
    f_1 = t == 1 ? 1 : 0
    f_2 = t == 2 ? 1 : 0
    m = m - 1  # married=2 becomes 1, single=1 becomes 0

    if r == 2 
        val = β_black_transfer_in[1] + β_black_transfer_in[2]*log(y) + β_black_transfer_in[3]*log(a_income +1 ) + β_black_transfer_in[4]*age + 
        β_black_transfer_in[5]*age^2 + β_black_transfer_in[6]*n + β_black_transfer_in[7]*e_0 + β_black_transfer_in[8]*e_2 + β_black_transfer_in[9]*m + 
        β_black_transfer_in[10]*f_1 + β_black_transfer_in[11]*f_2  + β_black_transfer_in[13]
    else
        val = β_white_transfer_in[1] + β_white_transfer_in[2]*log(y) + β_white_transfer_in[3]*log(a_income +1 ) + β_white_transfer_in[4]*age + 
        β_white_transfer_in[5]*age^2 + β_white_transfer_in[6]*n + β_white_transfer_in[7]*e_0 + β_white_transfer_in[8]*e_2 + β_white_transfer_in[9]*m + 
        β_white_transfer_in[10]*f_1 + β_white_transfer_in[11]*f_2  + β_white_transfer_in[13]
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


###############################################################################################
# Natural Borrowing Limit
###############################################################################################
function compute_borrowing_limit(z_grid, r_loan, jpnts, R, e, phi)
    z_min  = z_grid[R][1]          # lowest productivity draw z̄
    j_grad = e == 2 ? 3 : 5        # 2yr graduates j=3, 4yr graduates j=5

    # chi[j] = PV of minimum future earnings from j onward
    chi = zeros(jpnts + 1)          # chi[jpnts+1] = 0 (terminal)
    for j in jpnts:-1:1
        y_min  = j >= j_grad ? g(R, j, e, 1) * z_min : 0.0
        chi[j] = y_min + chi[j + 1] / (1 + r_loan)
    end

    # Borrowing limit: phi * PV of next period's earnings stream
    d_limit = zeros(jpnts)
    for j in 1:jpnts
        d_limit[j] = phi * chi[j + 1] / (1 + r_loan)
    end

    return d_limit
end

"""
function compute_natural_borrowing_limit(z_grid, r_loan, jpnts, R, e)
    z_min   = z_grid[R][1]
    c_floor = 0.001
    j_grad  = e == 2 ? 3 : 5    # 2yr graduates j=3, 4yr graduates j=5

    d_limit = zeros(jpnts)
    for j in (jpnts - 1):-1:1
        if j >= j_grad
            y_min      = g(R, j, e, 1) * z_min   # worst case: single (m=1), min z
            d_limit[j] = (y_min - c_floor + d_limit[j + 1]) / (1 + r_loan)
        else
            d_limit[j] = (-c_floor + d_limit[j + 1]) / (1 + r_loan)
        end
    end
    return d_limit
end
"""

###############################################################################################
# Family Shock Proabilities
###############################################################################################
# Load shock table
fam_shock_df = CSV.read("data/family_shock_table.csv", DataFrame)

# Build lookup: (R, e_college) -> vector of (m, n, t_idx, prob)
t_map = Dict("low" => 1, "mid" => 2, "high" => 3)

# e_college: collapse to 0 vs 1 to match shock table
e_college(e) = e < 3 ? 0 : 1

family_shock_probs = Dict{Tuple{Int,Int}, Vector{Tuple{Int,Int,Int,Float64}}}()

for row in eachrow(fam_shock_df)
    key = (row.R, row.e)
    outcome = (row.m, row.n, t_map[row.t], Float64(row.prob))
    push!(get!(family_shock_probs, key, Tuple{Int,Int,Int,Float64}[]), outcome)
end

function compute_t_own(y_p, a_inc_p)
    # Thresholds from PSID data — match how family_type was coded in 1.1
    total_p = y_p + a_inc_p
    if total_p < y_low_threshold
        return 1  # both_low (poor)
    elseif total_p > y_high_threshold
        return 3  # both_high (well off)
    else
        return 2  # both_mid
    end
end