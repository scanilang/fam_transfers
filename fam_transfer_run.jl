using Base.Threads
using Interpolations
import Optim: optimize, Brent
import QuantEcon: rouwenhorst
import Roots: find_zero
import Distributions: cdf, Normal
using Random
import CSV
using DataFrames
using Serialization   # for save/load
using Dates           # for timestamp in savepath

include("fam_transfers_setup.jl")
include("fam_transfers_model.jl")
include("fam_transfers_noschool.jl")
include("fam_transfers_school.jl")
include("fam_transfers_solve.jl")
include("fam_transfers_simulate.jl")

model = model_create();
savepath = "/Users/scanilang/Documents/econ/umn/family_transfers/value_policy/"
solution = solve_model(model, savepath)