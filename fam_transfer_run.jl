using Interpolations
import Optim: optimize
import QuantEcon: rouwenhorst
import Roots: find_zero, Brent
import Distributions: cdf, Normal
using Random

include("fam_transfers_simulate.jl")
include("fam_transfers_setup.jl")
include("fam_transfers_2026.jl")