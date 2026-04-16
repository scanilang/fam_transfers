using Interpolations
import Optim: optimize
import QuantEcon: rouwenhorst
import Roots: find_zero, Brent
import Distributions: cdf, Normal
using Random
import CSV
import DataFrames

include("fam_transfers_simulate.jl")
include("fam_transfers_setup.jl")
include("fam_transfers_noschool.jl")
include("fam_transfers_school.jl")
include("fam_transfers_solve.jl")

model = model_create();
solution = solve_model(model)
