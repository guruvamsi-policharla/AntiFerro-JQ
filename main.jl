#Includes
using Distributed
#addprocs(Sys.CPU_THREADS)
addprocs(8)
println(nprocs())

#@everywhere include("/home/guru/repos/AntiFerro-Lattices/skyrm_aux.jl")
#@everywhere include("/home/guru/repos/AntiFerro-Lattices/error_aux.jl")
#@everywhere include("/home/guru/repos/AntiFerro-Lattices/energy_aux.jl")
#@everywhere include("/home/guru/repos/AntiFerro-Lattices/lat_aux.jl")

@everywhere include("skyrm_aux.jl")
@everywhere include("error_aux.jl")
@everywhere include("energy_aux.jl")
@everywhere include("lat_aux.jl")

@everywhere using SharedArrays
@everywhere using Distributions
@everywhere using StatsBase
@everywhere using LinearAlgebra
@everywhere using JLD2
using Dates
#driver
Tmin = 0.1
Tchange = 0.2
Tmax = 2
N = 4
Temperature = Tmin:Tchange:Tmax
J_space = 0:0.2:2

E_temp = SharedArray{Float64,6}(length(Temperature),length(J_space),4,3,2,nprocs()-1)
mag_temp = SharedArray{Float64,6}(length(Temperature),length(J_space),4,3,2,nprocs()-1)
skyrm_temp = SharedArray{Float64,6}(length(Temperature),length(J_space),4,3,2,nprocs()-1)
magbind_temp = SharedArray{Float64,5}(length(Temperature),length(J_space),4,2,nprocs()-1)
skyrmbind_temp = SharedArray{Float64,5}(length(Temperature),length(J_space),4,2,nprocs()-1)

proc_complete = SharedArray{Int,1}(nprocs())

for i in 1:nprocs()
    proc_complete[i] = 0
end
@distributed for i in 2:nprocs()
	E_temp[:,:,:,:,:,i-1],
    skyrm_temp[:,:,:,:,:,i-1],
	mag_temp[:,:,:,:,:,i-1],
	magbind_temp[:,:,:,:,i-1],
	skyrmbind_temp[:,:,:,:,i-1] = fetch(@spawnat i montecarlo(Temperature,N,J_space))
    proc_complete[i] = 1
end

proc_complete[1] = 1

for i in 1:5000
    if(mean(proc_complete) == 1)
        println(proc_complete)
        @save "JQdata"*string(N)*"x"*string(N)*"fullresbind"*string(Dates.now())*".jld2" E_temp skyrm_temp mag_temp magbind_temp skyrmbind_temp Temperature N J_space
	break
    end
    println(proc_complete)
    sleep(35)
end

println("ending")
