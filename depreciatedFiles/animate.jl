# For file IO
using DelimitedFiles
include("helpers/fileHandling.jl")

# For animation
using Plots
using ColorSchemes
using Printf
include("helpers/animations.jl")

animateVorticity("output/testing", "test.gif", (-5, 5), 10)