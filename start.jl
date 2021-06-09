include("main.jl")

h = 1
iters = 1e3
file_IC = "vorticityFreqIC.csv"
type_IC = "vorticityFreq"	# One of {"velocity", "vorticityFreq"}
viz = "all"   				# One of {"velocity", "vorticity", "vorticityFreq", "all"}
lims_vorticity = (-0.1, 0.1)
lims_vorticityFreq = (-1, 1)

main(h, iters, file_IC, type_IC, viz, lims_vorticity, lims_vorticityFreq)