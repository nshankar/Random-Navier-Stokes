using Plots
using FFTW
using DelimitedFiles
using ColorSchemes

function animateVorticity(filein, fileout, lims)
	# read in data as Nx?xT matrix
	# irfft first 2 coords
	# animate real data
	f = open(filein, "r")
	data, header = readdlm(f, ',', ComplexF64, header=true)
	
	# Think about if there is a better way to do this
	h,N = split(header[1])
	h = parse(Float64,h)
	N = Int(parse(Float64,N))
	indexer = Int((N-1)/2) + 1

	close(f)

	grid = LinRange(0, 2*pi, N)
	anim = @animate for i = 1:Int(size(data)[1]/indexer)
		qhat = @view data[(i-1)*indexer + 1:i*indexer, :]
		heatmap(grid, grid, irfft(qhat, N)',
				c = :delta,	aspect_ratio=1, clims=lims, 
				show=false, title="Vorticity Field")
	end
	gif(anim, fileout, fps = 10)

end

# Important for plotting behavior
ENV["GKSwstype"]="nul"
animateVorticity("output/output5", "test5.gif", (-5, 5))