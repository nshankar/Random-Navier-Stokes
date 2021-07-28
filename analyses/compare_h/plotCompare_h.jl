using Plots
using ColorSchemes
#include("../../helpers/animations.jl")
# For file IO
using DelimitedFiles
include("../../helpers/fileHandling.jl")


function plotCompare_h()
	path = "/Users/nikhil/Documents/Duke/Research/Summer 2021/Random-Navier-Stokes/output/"
	folders = ["h=5e-5_1000cycles/", "h=1e-4_200cycles/", "h=4e-4_50cycles/", "h=8e-4_25cycles/", "h=1e-3_20cycles/"]
	repeats = [1, 2, 8, 16, 20]
	N = 61
	cycles = 400

	qMultipleSeries = zeros(N, N, cycles, length(folders))

	# Special case for j = 1
	folder = folders[1]
	h, N, ncycles, qSeries = getVorticityData(path*folder)
	qMultipleSeries[:,:,:,1] = qSeries[:,:,1:400]

	for j = 2:length(folders)
		folder = folders[j]
		h, N, ncycles, qSeries = getVorticityData(path*folder)
		qMultipleSeries[:,:,:,j] = repeat(qSeries[:,:,1:end-1], inner=(1,1, repeats[j]))
	end

	grid = LinRange(0, 2*pi, N)

	anim = @animate for i = 1:cycles
			p1 = heatmap(grid, grid, qMultipleSeries[:,:,i,1],
				c = :delta,	aspect_ratio=1, clims=(-5,5), 
				show=false, title="h=5e-5", legend=false)
		p2 = heatmap(grid, grid, qMultipleSeries[:,:,i,2],
				c = :delta,	aspect_ratio=1, clims=(-5,5), 
				show=false, title="h=1e-4", legend=false)
		p3 = heatmap(grid, grid, qMultipleSeries[:,:,i,3],
				c = :delta,	aspect_ratio=1, clims=(-5,5), 
				show=false, title="h=4e-4", legend=false)
		p4 = heatmap(grid, grid, qMultipleSeries[:,:,i,4],
				c = :delta,	aspect_ratio=1, clims=(-5,5), 
				show=false, title="h=8e-4", legend=false)

		p5 = heatmap(grid, grid, qMultipleSeries[:,:,i,5],
				c = :delta,	aspect_ratio=1, clims=(-5,5), 
				show=false, title="h=1e-3", legend=false)

		title = plot(title="Compare Evolutions",framestyle=nothing,showaxis=false,xticks=false,yticks=false,margin=0Plots.px)
		plot(title, p1, p2, p3, p4, p5, layout=@layout([A{0.05h}; [B C D]; [E F]]), size=(1350, 1100))
	end
	gif(anim, "analyses/compare_h/compare_h.gif", fps = cycles/10)

end

# Important for plotting behavior
ENV["GKSwstype"]="nul"
plotCompare_h()