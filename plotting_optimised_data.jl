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

f=jldopen("/home/vamsi/Github/AntiFerro-JQ port/Data/JQ0_2data4x4fullresbind2019-04-14T19:41:07.676.jld2","r")
E_temp = f["E_temp"].s
mag_temp = f["mag_temp"].s
skyrm_temp = f["skyrm_temp"].s
magbind_temp = f["magbind_temp"].s
skyrmbind_temp = f["skyrmbind_temp"].s
N = f["N"]
Temperature = f["Temperature"]
J_space = f["J_space"]

jstart = 1
jend = length(J_space)

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

#=
#magtemp
for jj in 1:3
figure()
    for ii in 1:4
        #figure()
        subplot(2,2,ii)
        for i in jstart:jend #Jspace
            if ii == 2
                errorbar(Temperature,(mag[:,i,ii,jj]+mag[:,i,ii+1,jj])/2,yerr = (mag_err[:,i,ii,jj]+mag_err[:,i,ii+1,jj])/2,fmt="o",linestyle="-",color=palette[mod(i,palsize)+1,:])
            else
                errorbar(Temperature,mag[:,i,ii,jj],yerr = mag_err[:,i,ii,jj],fmt="o",linestyle="-",color=palette[mod(i,palsize)+1,:])
            end
        end
        if ii == 1
            title("Magnetisation 00 - "*string(N)*"x"*string(N),fontsize = 17)
        elseif ii == 2
            title("Magnetisation (0pi+pi0)/2 - "*string(N)*"x"*string(N),fontsize = 17)
            legend("J2/J1 = ".*string.(J_space[jstart:jend]),bbox_to_anchor=[1.05,1],loc=2,ncol = 1)
        elseif ii == 3
            title("Magnetisation pi0 - "*string(N)*"x"*string(N),fontsize = 17)
        elseif ii == 4
            title("Magnetisation pipi - "*string(N)*"x"*string(N),fontsize = 17)
        end
        #legend("J1/J2 = ".*string.(J_space[jstart:jend]),bbox_to_anchor=[1.05,1],loc=2,ncol = 1)
        if ii>2
            xlabel("Temperature",fontsize = 14)
        end
        if mod(ii,2)==1
            if jj == 1
                ylabel(L"|mag|",fontsize = 14)
            elseif jj == 2
                ylabel(L"$ mag^2$",fontsize = 14)
            elseif jj == 3
                ylabel(L"$ mag^4$",fontsize = 14)
            end
        end
        grid("on")
    end
end

#skyrmtemp
for jj in 1:3
    figure()
    for ii in 1:4
        #figure()
        subplot(2,2,ii)
        for i in jstart:jend
            #errorbar(Temperature,mean(data["mag"*string(i)],2),mean(data["mag_err"*string(i)],2))
            #if ii == 2
            #    errorbar(Temperature,(skyrm[:,i,ii,jj]+skyrm[:,i,ii+1,jj])/2,yerr = (skyrm_err[:,i,ii,jj]+skyrm_err[:,i,ii+1,jj])/2,fmt="o",linestyle="-",color=palette[mod(i,palsize)+1,:])
            #else
                errorbar(Temperature,skyrm[:,i,ii,jj],yerr = skyrm_err[:,i,ii,jj],fmt="o",linestyle="-",color=palette[mod(i,palsize)+1,:])
            #end
        end
        if ii == 1
            title("Skyrmion 00 - "*string(N)*"x"*string(N))
        elseif ii == 2
            title("Skyrmion (0pi+pi0)/2 - "*string(N)*"x"*string(N))
            legend("J1/J2 = ".*string.(J_space[jstart:jend]),bbox_to_anchor=[1.05,1],loc=2,ncol = 1)
        elseif ii == 3
            title("Skyrmion pi0 - "*string(N)*"x"*string(N))
        elseif ii == 4
            title("Skyrmion pipi - "*string(N)*"x"*string(N))
        end
        if ii>2
            xlabel("Temperature")
        end
        if mod(ii,2)==1
            if jj == 1
                ylabel("abs(skyrm)")
            elseif jj == 2
                ylabel("skyrm.^2")
            elseif jj == 3
                ylabel("skyrm.^4")
            end
        end
        grid("on")
    end
end
=#

#=
#magj1j2
for jj in 1:3
    figure()
    for ii in 1:4
        subplot(2,2,ii)
        for i in 1:1:length(Temperature)
            errorbar(J_space[jstart:jend],mag[i,jstart:jend,ii,jj,1],mag[i,jstart:jend,ii,jj,2],fmt="o",linestyle="-",color=palette[mod(3*i-3,palsize)+1,:])
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
            xlabel(L"$Q/J$",fontsize = 14)
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
        #axvline(x=0.5,linestyle="-.",color="r")
        grid("on")
    end
end

#skyrmj1j2
for jj in 1:3
    figure()
    for ii in 1:4
        subplot(2,2,ii)
        for i in 1:1:length(Temperature)
            errorbar(J_space[jstart:jend],skyrm[i,jstart:jend,ii,jj,1],skyrm[i,jstart:jend,ii,jj,2],fmt="o",linestyle="-",color=palette[mod(3*i-3,palsize)+1,:])
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
            xlabel(L"$Q/J$",fontsize = 14)
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
        #axvline(x=0.5,linestyle="-.",color="r")
        grid("on")
    end
end
=#

#=
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
            xlabel(L"$Q/J$",fontsize = 14)
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


#Susceptibilityj1j2

#=
figure()
for ii in 1:4
    jj = 2 #M^2 terms
    subplot(2,2,ii)
    for i in 1:1:length(Temperature)
        errorbar(J_space,mag[i,:,ii,jj,1].*(N^4/Temperature[i]),mag[i,:,ii,jj,2].*(N^4/Temperature[i]),fmt="o",linestyle="-",color=palette[mod(3*i-3,palsize)+1,:])
    end
    if ii == 1
        title("Susceptibility 00 - "*string(N)*"x"*string(N),fontsize = 17)
    elseif ii == 2
        title(L"Susceptibility $(0\pi+ \pi 0)/2$ - "*string(N)*"x"*string(N),fontsize = 17)
        legend("T = ".*string.(Temperature[1:1:end]),bbox_to_anchor=[1.05,1],loc=2,ncol = 1)
    elseif ii == 3
        title(L"Susceptibility $\pi 0$ - "*string(N)*"x"*string(N),fontsize = 17)
    elseif ii == 4
        title(L"Susceptibility $\pi \pi$ - "*string(N)*"x"*string(N),fontsize = 17)
    end
    if ii>2
        xlabel(L"$Q/J$",fontsize = 14)
    end
    if mod(ii,2)==1
        if jj == 1
            ylabel(L"|mag|",fontsize = 14)
        elseif jj == 2
            ylabel(L"$ \langle \chi \rangle$",fontsize = 14)
        elseif jj == 3
            ylabel(L"$ mag^4$",fontsize = 14)
        end
    end
    #axvline(x=0.5,linestyle="-.",color="r")
    grid("on")
end
=#


#SkyrmSusceptibilityj1j2

figure()
for ii in 1:4
    jj = 2 #M^2 terms
    subplot(2,2,ii)
    for i in 1:1:length(Temperature)
        errorbar(J_space,skyrm[i,:,ii,jj,1].*(N^4/Temperature[i]),skyrm[i,:,ii,jj,2].*(N^4/Temperature[i]),fmt="o",linestyle="-",color=palette[mod(3*i-3,palsize)+1,:])
    end
    if ii == 1
        title("Susceptibility 00 - "*string(N)*"x"*string(N),fontsize = 17)
    elseif ii == 2
        title(L"Susceptibility $(0\pi+ \pi 0)/2$ - "*string(N)*"x"*string(N),fontsize = 17)
        legend("T = ".*string.(Temperature[1:1:end]),bbox_to_anchor=[1.05,1],loc=2,ncol = 1)
    elseif ii == 3
        title(L"Susceptibility $\pi 0$ - "*string(N)*"x"*string(N),fontsize = 17)
    elseif ii == 4
        title(L"Susceptibility $\pi \pi$ - "*string(N)*"x"*string(N),fontsize = 17)
    end
    if ii>2
        xlabel(L"$Q/J$",fontsize = 14)
    end
    if mod(ii,2)==1
        if jj == 1
            ylabel(L"|mag|",fontsize = 14)
        elseif jj == 2
            ylabel(L"$ \langle \chi_{skyrm} \rangle$",fontsize = 14)
        elseif jj == 3
            ylabel(L"$ mag^4$",fontsize = 14)
        end
    end
    #axvline(x=0.5,linestyle="-.",color="r")
    grid("on")
end


#SKYRMBIND

#=
figure()
for ii in 1:4
    subplot(2,2,ii)
    for i in 1:1:length(Temperature)
        errorbar(J_space[jstart:jend],skyrmbind[i,jstart:jend,ii,1],skyrmbind[i,jstart:jend,ii,2],fmt="o",linestyle="-",color=palette[mod(3*i-3,palsize)+1,:])
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
    if ii>2
        xlabel(L"$Q/J$",fontsize = 14)
    end
    if mod(ii,2)==1
            ylabel(L"$\langle skyrm^2 \rangle ^2 / \langle skyrm^4 \rangle$",fontsize = 14)
    end
    #axvline(x=0.5,linestyle="-.",color="r")
    grid("on")
end


#MAGBINDER


figure()
for ii in 1:4
    subplot(2,2,ii)
    for i in 1:1:length(Temperature)
        errorbar(J_space[jstart:jend],magbind[i,jstart:jend,ii,1],magbind[i,jstart:jend,ii,2],fmt="o",linestyle="-",color=palette[mod(3*i-3,palsize)+1,:])
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
        xlabel(L"$Q/J$",fontsize = 14)
    end
    if mod(ii,2)==1
            ylabel(L"$\langle mag^2 \rangle ^2 / \langle mag^4 \rangle$",fontsize = 14)
    end
    #axvline(x=0.5,linestyle="-.",color="r")
    grid("on")
end
=#
