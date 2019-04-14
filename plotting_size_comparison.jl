using PyPlot
using Statistics
using JLD2
using PlotUtils
using Images
palsize = 30
cm = cgrad(:plasma);
space = 1
cm_temp = [cm[i] for i in 1:1:200]
palette = Array{Float64,2}(undef,palsize,3)
for i in 1:palsize
    palette[i,:] = [red(cm_temp[i]),green(cm_temp[i]),blue(cm_temp[i])]
end

nanmean(x) = mean(filter(!isnan,x))
nanmean(x,y) = mapslices(nanmean,x,y)

function process_data(f)

    E_temp = f["E_temp"].s
    mag_temp = f["mag_temp"].s
    skyrm_temp = f["skyrm_temp"].s
    magbind_temp = f["magbind_temp"].s
    skyrmbind_temp = f["skyrmbind_temp"].s
    N = f["N"]
    Temperature = f["Temperature"]
    J_space = f["J_space"]

    skyrm = Array{Float64,5}(undef,length(Temperature),length(J_space),4,3,2)
    mag = Array{Float64,5}(undef,length(Temperature),length(J_space),4,3,2)
    E = Array{Float64,5}(undef,length(Temperature),length(J_space),4,3,2)

    magbind = Array{Float64,4}(undef,length(Temperature),length(J_space),4,2)
    skyrmbind = Array{Float64,4}(undef,length(Temperature),length(J_space),4,2)

    E[:,:,:,:,1] = reshape(mean(E_temp[:,:,:,:,1,:],dims=5),(size(E_temp,1),size(E_temp,2),4,3))
    E[:,:,:,:,2] = reshape(sqrt.(mean(E_temp[:,:,:,:,2,:].^2,dims=5)),(size(E_temp,1),size(E_temp,2),4,3))

    mag[:,:,:,:,1] = reshape(mean(mag_temp[:,:,:,:,1,:],dims=5),(size(mag_temp,1),size(mag_temp,2),4,3))
    mag[:,:,:,:,2] = reshape(sqrt.(mean(mag_temp[:,:,:,:,2,:].^2,dims=5)),(size(mag_temp,1),size(mag_temp,2),4,3))

    skyrm[:,:,:,:,1] = reshape(mean(skyrm_temp[:,:,:,:,1,:],dims=5),(size(skyrm_temp,1),size(skyrm_temp,2),4,3))
    skyrm[:,:,:,:,2] = reshape(sqrt.(mean(skyrm_temp[:,:,:,:,2,:].^2,dims=5)),(size(skyrm_temp,1),size(skyrm_temp,2),4,3))

    magbind[:,:,:,1] = reshape(mean(magbind_temp[:,:,:,1,:],dims=4),(size(magbind_temp,1),size(magbind_temp,2),4))
    magbind[:,:,:,2] = reshape(sqrt.(mean(magbind_temp[:,:,:,2,:].^2,dims=4)),(size(magbind_temp,1),size(magbind_temp,2),4))

    skyrmbind[:,:,:,1] = reshape(meanfinite(skyrmbind_temp[:,:,:,1,:],4),(size(skyrmbind_temp,1),size(skyrmbind_temp,2),4))
    skyrmbind[:,:,:,2] = reshape(sqrt.(meanfinite(skyrmbind_temp[:,:,:,2,:].^2,4)),(size(skyrmbind_temp,1),size(skyrmbind_temp,2),4))

    return N, Temperature, J_space, E, mag, skyrm, magbind, skyrmbind
end

function plot_single(N, Temperature, J_space, E, mag, skyrm, magbind, skyrmbind)
    #skyrmj1j2
    for jj in 1:3
        figure()
        for ii in 1:4
            subplot(2,2,ii)
            for i in 1:1:length(Temperature)
                errorbar(J_space,skyrm[i,:,ii,jj,1],skyrm[i,:,ii,jj,2],fmt="o",linestyle="-",color=palette[mod(3*i-3,palsize)+1,:])
            end
            if ii == 1
                title("Skyrmion 00 - "*string(N)*"x"*string(N),fontsize = 17)
            elseif ii == 2
                title(L"Skyrmion $(0\pi+\pi 0)/2$ - "*string(N)*"x"*string(N),fontsize = 17)
                legend("T = ".*string.(Temperature[1:1:end]),bbox_to_anchor=[1.05,1],loc=2,ncol = 1)
            elseif ii == 3
                title(L"Skyrmion $\pi 0$ - "*string(N)*"x"*string(N),fontsize = 17)
            elseif ii == 4
                title(L"Skyrmion $\pi \pi$ - "*string(N)*"x"*string(N),fontsize = 17)
            end
            #legend("T = ".*string.(Temperature[1:1:end]),bbox_to_anchor=[1.05,1],loc=2,ncol = 1)
            if ii>2
                xlabel(L"$J_2/J_1$",fontsize = 14)
            end
            if mod(ii,2)==1
                if jj == 1
                    ylabel(L"|skyrm|",fontsize = 14)
                elseif jj == 2
                    ylabel(L"$\langle skyrm^2 \rangle$",fontsize = 14)
                elseif jj == 3
                    ylabel(L"skyrm^4",fontsize = 14)
                end
            end
            axvline(x=0.5,linestyle="-.",color="r")
            grid("on")
        end
    end

    #=
    #magj1j2
    for jj in 1:3
        figure()
        for ii in 1:4
            subplot(2,2,ii)
            for i in 1:1:length(Temperature)
                errorbar(J_space,mag[i,:,ii,jj,1],mag[i,:,ii,jj,2],fmt="o",linestyle="-",color=palette[mod(3*i-3,palsize)+1,:])
            end
            if ii == 1
                title("Magnetisation 00 - "*string(N)*"x"*string(N),fontsize = 17)
            elseif ii == 2
                title(L"Magnetisation $(0\pi+ \pi 0)/2$ - "*string(N)*"x"*string(N),fontsize = 17)
                legend("T = ".*string.(Temperature[1:1:end]),bbox_to_anchor=[1.05,1],loc=2,ncol = 1)
            elseif ii == 3
                title(L"Magnetisation $\pi 0$ - "*string(N)*"x"*string(N),fontsize = 17)
            elseif ii == 4
                title(L"Magnetisation $\pi \pi$ - "*string(N)*"x"*string(N),fontsize = 17)
            end
            if ii>2
                xlabel(L"$J_2/J_1$",fontsize = 14)
            end
            if mod(ii,2)==1
                if jj == 1
                    ylabel(L"|mag|",fontsize = 14)
                elseif jj == 2
                    ylabel(L"$ \langle mag^2 \rangle$",fontsize = 14)
                elseif jj == 3
                    ylabel(L"$ mag^4$",fontsize = 14)
                end
            end
            axvline(x=0.5,linestyle="-.",color="r")
            grid("on")
        end
    end

    #Ej1j2
    for jj in 1:3
        figure()
        for ii in 1:4
            subplot(2,2,ii)
            for i in 1:1:length(Temperature)
                errorbar(J_space,E[i,:,ii,jj,1],E[i,:,ii,jj,2],fmt="o",linestyle="-",color=palette[mod(3*i-3,palsize)+1,:])
            end
            if ii == 1
                title("Energy 00 - "*string(N)*"x"*string(N),fontsize = 17)
            elseif ii == 2
                title(L"Energy $(0\pi+ \pi 0)/2$ - "*string(N)*"x"*string(N),fontsize = 17)
                legend("T = ".*string.(Temperature[1:1:end]),bbox_to_anchor=[1.05,1],loc=2,ncol = 1)
            elseif ii == 3
                title(L"Energy $\pi 0$ - "*string(N)*"x"*string(N),fontsize = 17)
            elseif ii == 4
                title(L"Energy $\pi \pi$ - "*string(N)*"x"*string(N),fontsize = 17)
            end
            if ii>2
                xlabel(L"$J_2/J_1$",fontsize = 14)
            end
            if mod(ii,2)==1
                if jj == 1
                    ylabel(L"E",fontsize = 14)
                elseif jj == 2
                    ylabel(L"$ \langle E^2 \rangle$",fontsize = 14)
                elseif jj == 3
                    ylabel(L"$ \langle E^4 \rangle$",fontsize = 14)
                end
            end
            axvline(x=0.5,linestyle="-.",color="r")
            grid("on")
        end
    end
    =#
end

function plot_list(N_list, Temperature_list, J_space_list, E_list, mag_list, skyrm_list, magbind_list, skyrmbind_list)

    #skyrmj1j2
    for jj in 1:3
        figure()
        for ii in 1:4
            subplot(2,2,ii)
            for j in 1:length(N_list)
                N = N_list[j]
                Temperature = Temperature_list[j]
                J_space = J_space_list[j]
                E = E_list[j]
                mag = mag_list[j]
                skyrm = skyrm_list[j]
                magbind = magbind_list[j]
                skyrmbind = skyrmbind_list[j]
                for i in 1:1:1
                    errorbar(J_space,skyrm[i,:,ii,jj,1],skyrm[i,:,ii,jj,2],fmt="o",linestyle="-")
                end
            end
            if ii == 1
                title("Skyrmion 00",fontsize = 17)
            elseif ii == 2
                title(L"Skyrmion $(0\pi+\pi 0)/2$",fontsize = 17)
                legend("N = ".*string.(N_list[1:1:end]),bbox_to_anchor=[1.05,1],loc=2,ncol = 1)
            elseif ii == 3
                title(L"Skyrmion $\pi 0$",fontsize = 17)
            elseif ii == 4
                title(L"Skyrmion $\pi \pi$",fontsize = 17)
            end
            #legend("T = ".*string.(Temperature[1:1:end]),bbox_to_anchor=[1.05,1],loc=2,ncol = 1)
            if ii>2
                xlabel(L"$J_2/J_1$",fontsize = 14)
            end
            if mod(ii,2)==1
                if jj == 1
                    ylabel(L"|skyrm|",fontsize = 14)
                elseif jj == 2
                    ylabel(L"$\langle skyrm^2 \rangle$",fontsize = 14)
                elseif jj == 3
                    ylabel(L"skyrm^4",fontsize = 14)
                end
            end
            axvline(x=0.5,linestyle="-.",color="r")
            grid("on")
        end
    end


    #magj1j2
    for jj in 1:3
        figure()
        for ii in 1:4
            subplot(2,2,ii)
            for j in 1:length(N_list)
                N = N_list[j]
                Temperature = Temperature_list[j]
                J_space = J_space_list[j]
                E = E_list[j]
                mag = mag_list[j]
                skyrm = skyrm_list[j]
                magbind = magbind_list[j]
                skyrmbind = skyrmbind_list[j]
                for i in 1:1:1
                    errorbar(J_space,mag[i,:,ii,jj,1],mag[i,:,ii,jj,2],fmt="o",linestyle="-")
                end
            end
            if ii == 1
                title("Magnetisation 00",fontsize = 17)
            elseif ii == 2
                title(L"Magnetisation $(0\pi+ \pi 0)/2$",fontsize = 17)
                legend("N = ".*string.(N_list[1:1:end]),bbox_to_anchor=[1.05,1],loc=2,ncol = 1)
            elseif ii == 3
                title(L"Magnetisation $\pi 0$",fontsize = 17)
            elseif ii == 4
                title(L"Magnetisation $\pi \pi$",fontsize = 17)
            end
            if ii>2
                xlabel(L"$J_2/J_1$",fontsize = 14)
            end
            if mod(ii,2)==1
                if jj == 1
                    ylabel(L"|mag|",fontsize = 14)
                elseif jj == 2
                    ylabel(L"$ \langle mag^2 \rangle$",fontsize = 14)
                elseif jj == 3
                    ylabel(L"$ mag^4$",fontsize = 14)
                end
            end
            axvline(x=0.5,linestyle="-.",color="r")
            grid("on")
        end
    end


    #Ej1j2
    for jj in 1:3
        figure()
        for ii in 1:4
            subplot(2,2,ii)
            for j in 1:length(N_list)
                N = N_list[j]
                Temperature = Temperature_list[j]
                J_space = J_space_list[j]
                E = E_list[j]
                mag = mag_list[j]
                skyrm = skyrm_list[j]
                magbind = magbind_list[j]
                skyrmbind = skyrmbind_list[j]
                for i in 1:1:1
                    errorbar(J_space,E[i,:,ii,jj,1],E[i,:,ii,jj,2],fmt="o",linestyle="-")
                end
            end

            if ii == 1
                title("Energy 00 - "*string(N)*"x"*string(N),fontsize = 17)
            elseif ii == 2
                title(L"Energy $(0\pi+ \pi 0)/2$ - "*string(N)*"x"*string(N),fontsize = 17)
                legend("N = ".*string.(N_list[1:1:end]),bbox_to_anchor=[1.05,1],loc=2,ncol = 1)
            elseif ii == 3
                title(L"Energy $\pi 0$ - "*string(N)*"x"*string(N),fontsize = 17)
            elseif ii == 4
                title(L"Energy $\pi \pi$ - "*string(N)*"x"*string(N),fontsize = 17)
            end
            if ii>2
                xlabel(L"$J_2/J_1$",fontsize = 14)
            end
            if mod(ii,2)==1
                if jj == 1
                    ylabel(L"E",fontsize = 14)
                elseif jj == 2
                    ylabel(L"$ \langle E^2 \rangle$",fontsize = 14)
                elseif jj == 3
                    ylabel(L"$ \langle E^4 \rangle$",fontsize = 14)
                end
            end
            axvline(x=0.5,linestyle="-.",color="r")
            grid("on")
        end
    end


    #SKYRMBIND
    figure()
    for ii in 1:4
        subplot(2,2,ii)
        for j in 1:length(N_list)
            N = N_list[j]
            Temperature = Temperature_list[j]
            J_space = J_space_list[j]
            E = E_list[j]
            mag = mag_list[j]
            skyrm = skyrm_list[j]
            magbind = magbind_list[j]
            skyrmbind = skyrmbind_list[j]
            for i in 1:1:1
                errorbar(J_space,skyrmbind[i,:,ii,1],skyrmbind[i,:,ii,2],fmt="o",linestyle="-")
            end
        end

        if ii == 1
            title("Skyrmion Binder 00",fontsize = 17)
        elseif ii == 2
            title(L"Skyrmion Binder $(0\pi+\pi 0)/2$",fontsize = 17)
            legend("N = ".*string.(N_list[1:1:end]),bbox_to_anchor=[1.05,1],loc=2,ncol = 1)
        elseif ii == 3
            title(L"Skyrmion Binder $\pi 0$",fontsize = 17)
        elseif ii == 4
            title(L"Skyrmion Binder $\pi \pi$",fontsize = 17)
        end
        if ii>2
            xlabel(L"$J_1/J_2$",fontsize = 14)
        end
        if mod(ii,2)==1
                ylabel(L"$\langle skyrm^2 \rangle ^2 / \langle skyrm^4 \rangle$",fontsize = 14)
        end
        axvline(x=0.5,linestyle="-.",color="r")
        grid("on")
    end


    #MAGBINDER
    figure()
    for ii in 1:4
        subplot(2,2,ii)
        for j in 1:length(N_list)
            N = N_list[j]
            Temperature = Temperature_list[j]
            J_space = J_space_list[j]
            E = E_list[j]
            mag = mag_list[j]
            skyrm = skyrm_list[j]
            magbind = magbind_list[j]
            skyrmbind = skyrmbind_list[j]
            for i in 1:1:1
                errorbar(J_space,magbind[i,:,ii,1],magbind[i,:,ii,2],fmt="o",linestyle="-")
            end
        end
        if ii == 1
            title("Magnetisation Binder 00",fontsize = 17)
        elseif ii == 2
            title(L"Magnetisation Binder $(0\pi+ \pi 0)/2$",fontsize = 17)
            legend("N = ".*string.(N_list[1:1:end]),bbox_to_anchor=[1.05,1],loc=2,ncol = 1)
        elseif ii == 3
            title(L"Magnetisation Binder $\pi 0$",fontsize = 17)
        elseif ii == 4
            title(L"Magnetisation Binder $\pi \pi$",fontsize = 17)
        end
        if ii>2
            xlabel(L"$J_1/J_2$",fontsize = 14)
        end
        if mod(ii,2)==1
                ylabel(L"$\langle mag^2 \rangle ^2\langle mag^4 \rangle/$",fontsize = 14)
        end
        axvline(x=0.5,linestyle="-.",color="r")
        grid("on")
    end


    #MAGBINDERXT
    figure()
    for ii in 1:4
        subplot(2,2,ii)
        for j in 1:length(N_list)
            N = N_list[j]
            Temperature = Temperature_list[j]
            J_space = J_space_list[j]
            E = E_list[j]
            mag = mag_list[j]
            skyrm = skyrm_list[j]
            magbind = magbind_list[j]
            skyrmbind = skyrmbind_list[j]
            for i in 15:15:15
                errorbar(Temperature,magbind[:,i,ii,1],magbind[:,i,ii,2],fmt="o",linestyle="-")
            end
        end
        if ii == 1
            title("Magnetisation Binder 00 J1/J2 = 1",fontsize = 17)
        elseif ii == 2
            title(L"Magnetisation Binder $(0\pi+ \pi 0)/2$",fontsize = 17)
            legend("N = ".*string.(N_list[1:1:end]),bbox_to_anchor=[1.05,1],loc=2,ncol = 1)
        elseif ii == 3
            title(L"Magnetisation Binder $\pi 0$",fontsize = 17)
        elseif ii == 4
            title(L"Magnetisation Binder $\pi \pi$",fontsize = 17)
        end
        if ii>2
            xlabel("T",fontsize = 14)
        end
        if mod(ii,2)==1
                ylabel(L"$\langle mag^2 \rangle ^2 \langle mag^4 \rangle/$",fontsize = 14)
        end
        axvline(x=0.5,linestyle="-.",color="r")
        grid("on")
    end
end


function plot_arr(N_arr,Temperature,J_space,E_arr,mag_arr,skyrm_arr,magbind_arr,skyrmbind_arr)
    #mag_sizes
    for jj in 1:3
        figure()
        for ii in 1:4
            subplot(2,2,ii)
            for j in 1:length(Temperature)
                for i in 15:15:15
                    errorbar(inv.(N_arr),mag_arr[j,i,ii,jj,1,:],mag_arr[j,i,ii,jj,2,:],fmt="o",linestyle="-")
                end
            end
            if ii == 1
                title("Magnetisation 00 J2/J1 = 1",fontsize = 17)
            elseif ii == 2
                title(L"Magnetisation $(0\pi+ \pi 0)/2$",fontsize = 17)
                legend("T = ".*string.(Temperature[1:1:end]),bbox_to_anchor=[1.05,1],loc=2,ncol = 1)
            elseif ii == 3
                title(L"Magnetisation $\pi 0$",fontsize = 17)
            elseif ii == 4
                title(L"Magnetisation $\pi \pi$",fontsize = 17)
            end
            if ii>2
                xlabel(L"$J_2/J_1$",fontsize = 14)
            end
            if mod(ii,2)==1
                if jj == 1
                    ylabel(L"|mag|",fontsize = 14)
                elseif jj == 2
                    ylabel(L"$ \langle mag^2 \rangle$",fontsize = 14)
                elseif jj == 3
                    ylabel(L"$ mag^4$",fontsize = 14)
                end
            end
            grid("on")
        end
    end

end

################################################################################
#=
f=jldopen("/home/vamsi/Github/AntiFerro-Lattices/Data/4x4/data4x4fullresbind2019-03-13T05:26:19.273.jld2","r")
N, Temperature, J_space, E, mag, skyrm, magbind, skyrmbind = process_data(f)

N_list = [N]
Temperature_list = [Temperature]
J_space_list = [J_space]
E_list = [E]
mag_list = [mag]
skyrm_list = [skyrm]
magbind_list = [magbind]
skyrmbind_list = [skyrmbind]

f=jldopen("/home/vamsi/Github/AntiFerro-Lattices/Data/16x16/data16x16fullresbind2019-03-15T10:57:26.978.jld2","r")
N, Temperature, J_space, E, mag, skyrm, magbind, skyrmbind = process_data(f)

N_list = push!(N_list,N)
Temperature_list = push!(Temperature_list,Temperature)
J_space_list = push!(J_space_list,J_space)
E_list = push!(E_list,E)
mag_list = push!(mag_list,mag)
skyrm_list = push!(skyrm_list,skyrm)
magbind_list = push!(magbind_list,magbind)
skyrmbind_list = push!(skyrmbind_list,skyrmbind)

f=jldopen("/home/vamsi/Github/AntiFerro-Lattices/Data/32x32/data32x32fullresbind2019-03-16T01:46:14.997.jld2","r")

N, Temperature, J_space, E, mag, skyrm, magbind, skyrmbind = process_data(f)

N_list = push!(N_list,N)
Temperature_list = push!(Temperature_list,Temperature)
J_space_list = push!(J_space_list,J_space)
E_list = push!(E_list,E)
mag_list = push!(mag_list,mag)
skyrm_list = push!(skyrm_list,skyrm)
magbind_list = push!(magbind_list,magbind)
skyrmbind_list = push!(skyrmbind_list,skyrmbind)
=#


################################################################################



f=jldopen("/home/vamsi/Github/AntiFerro-Lattices/Data/4x4/data4x4fullresbind2019-03-13T05:26:53.877.jld2","r")
N, Temperature, J_space, E, mag, skyrm, magbind, skyrmbind = process_data(f)

N_list = [N]
Temperature_list = [Temperature]
J_space_list = [J_space]
E_list = [E]
mag_list = [mag]
skyrm_list = [skyrm]
magbind_list = [magbind]
skyrmbind_list = [skyrmbind]

N_arr = N

E_arr = E
mag_arr = mag
skyrm_arr = skyrm
magbind_arr = magbind
skyrmbind_arr = skyrmbind

f=jldopen("/home/vamsi/Github/AntiFerro-Lattices/Data/16x16/data16x16fullresbind2019-03-15T15:12:44.616.jld2","r")
N, Temperature, J_space, E, mag, skyrm, magbind, skyrmbind = process_data(f)

N_list = push!(N_list,N)
Temperature_list = push!(Temperature_list,Temperature)
J_space_list = push!(J_space_list,J_space)
E_list = push!(E_list,E)
mag_list = push!(mag_list,mag)
skyrm_list = push!(skyrm_list,skyrm)
magbind_list = push!(magbind_list,magbind)
skyrmbind_list = push!(skyrmbind_list,skyrmbind)


N_arr = cat(N_arr,N,dims = ndims(N)+1)
E_arr = cat(E_arr,E,dims = ndims(E)+1)
mag_arr = cat(mag_arr,mag,dims = ndims(mag)+1)
skyrm_arr = cat(skyrm_arr,skyrm,dims = ndims(skyrm)+1)
magbind_arr = cat(magbind_arr,magbind,dims = ndims(magbind)+1)
skyrmbind_arr = cat(skyrmbind_arr,skyrmbind,dims = ndims(skyrmbind)+1)

f=jldopen("/home/vamsi/Github/AntiFerro-Lattices/Data/32x32/data32x32fullresbind2019-03-16T02:11:03.114.jld2","r")
N, Temperature, J_space, E, mag, skyrm, magbind, skyrmbind = process_data(f)

N_list = push!(N_list,N)
Temperature_list = push!(Temperature_list,Temperature)
J_space_list = push!(J_space_list,J_space)
E_list = push!(E_list,E)
mag_list = push!(mag_list,mag)
skyrm_list = push!(skyrm_list,skyrm)
magbind_list = push!(magbind_list,magbind)
skyrmbind_list = push!(skyrmbind_list,skyrmbind)


N_arr = cat(N_arr,N,dims = ndims(N)+1)
E_arr = cat(E_arr,E,dims = ndims(E)+1)
mag_arr = cat(mag_arr,mag,dims = ndims(mag)+1)
skyrm_arr = cat(skyrm_arr,skyrm,dims = ndims(skyrm)+1)
magbind_arr = cat(magbind_arr,magbind,dims = ndims(magbind)+1)
skyrmbind_arr = cat(skyrmbind_arr,skyrmbind,dims = ndims(skyrmbind)+1)


#plot_list(N_list, Temperature_list, J_space_list, E_list, mag_list, skyrm_list, magbind_list, skyrmbind_list)

plot_arr(N_arr,Temperature,J_space,E_arr,mag_arr,skyrm_arr,magbind_arr,skyrmbind_arr)
