using PyPlot

filename = "data_10.dat"

lines = open(filename, "r") do fp
	readlines(fp)
end	

d1 = Array{Float64, 1}()
d2 = Array{Float64, 1}()
for line in lines
  d = parse.(Float64, split(line))
  push!(d1, d[1])
  push!(d2, d[2])
end

figure(figsize=(8,6))
plot(d1, d2, "-")
show()
