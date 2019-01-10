function sample_gauss(v)
    #Input is the central vector around which we flip
    var = 0.2
    v_new = [rand(Normal(v[1],var)),rand(Normal(v[2],var)),rand(Normal(v[3],var))]
    v_new = v_new/sqrt(dot(v_new,v_new))
    return convert(Vector,v_new)
end

function sample_uni()
    x = [rand(Uniform(-1,1)), rand(Uniform(-1,1)), rand(Uniform(-1,1))]
    x = x/sqrt(dot(x,x))
    return convert(Vector,x)
end

function initialise(M::Int, N::Int)
    """ Initialising the lattice with random values """
    #lat = Array{Float64, 3}(N, N)
    lat = Array{Vector{Float64},2}(undef,M, N);
    for i = 1:M
        for j = 1:N
            lat[i,j] = sample_uni()
        end
    end
    return lat
end

function initialise_gauss(M::Int, N::Int)
    """ Initialising the lattice with random values """
    #lat = Array{Float64, 3}(N, N)
    lat = Array{Vector{Float64},2}(M, N);
    for i = 1:M
        for j = 1:N
            lat[i,j] = sample_gauss([1,0,0])
        end
    end
    return lat
end

function four_trans(lat)
    x = zeros(M,N)
    y = zeros(M,N)
    z = zeros(M,N)
    for i in 1:M
      for j in 1:N
        x[i,j] = lat[i,j][1]
        y[i,j] = lat[i,j][2]
        z[i,j] = lat[i,j][3]
      end
    end

    ftx = abs.(fft(x))
    fty = abs.(fft(y))
    ftz = abs.(fft(z))

    ft = abs.(sqrt.(ftx.*ftx + fty.*fty + ftz.*ftz))
return ft
end

function lat_transform(lat,latindex)
    N = size(lat,1)
    if latindex == 2
        for ii in 1:N
            for jj in 1:N
                lat[ii,jj] = (-1)^(jj).*lat[ii,jj]
            end
        end
    elseif latindex == 3
        for ii in 1:N
            for jj in 1:N
                lat[ii,jj] = (-1)^(ii).*lat[ii,jj]
            end
        end
    elseif latindex == 4
        for ii in 1:N
            for jj in 1:N
                lat[ii,jj] = (-1)^(ii+jj).*lat[ii,jj]
            end
        end
    end
end

function transient_results(lat, transient::Int, J, T)
    """Takes lat as input and removes initial transients by running for transient number of steps"""
    M = size(lat,1)
    N = size(lat,2)
    for i = 1:transient
        for j = 1:M*N
                x = rand(1:M)
                y = rand(1:N)
                test_flip(x,y,J,lat,T)
        end
    end
end

function montecarlo(Temperature,N,J_space)
    mcs = 25
    M = N

    normalisation=(1.0/float(M*N))

    JM_vec = zeros(length(Temperature),length(J_space),4,3)
    JM_vec_err = zeros(length(Temperature),length(J_space),4,3)
    Jskyrm_vec = zeros(length(Temperature),length(J_space),4,3)
    Jskyrm_vec_err = zeros(length(Temperature),length(J_space),4,3)
    Jmagbind_vec = zeros(length(Temperature),length(J_space),4)
    Jskyrmbind_vec = zeros(length(Temperature),length(J_space),4)
    Jmagbind_vec_err = zeros(length(Temperature),length(J_space),4)
    Jskyrmbind_vec_err = zeros(length(Temperature),length(J_space),4)
#################################################################
    M_vec = zeros(length(Temperature),2,3,4)
    M_jack = zeros(mcs,3,4)

    skyrm_vec = zeros(length(Temperature),2,3,4)
    skyrm_jack = zeros(mcs,3,4)

    magbind_vec = zeros(length(Temperature),2,4)
    skyrmbind_vec = zeros(length(Temperature),2,4)

    #autocor_vec = 0
    Jcount = 1
    for J in J_space
        lat = initialise(M,N)
        Tcount = 1
        for T in Temperature
            transient_results(lat,3000,J,T)
            E = total_energy(J,lat)
            for i in 1:mcs
                for j in 1:M*N
                    x = rand(1:M)
                    y = rand(1:N)
                    E_0 = energy_pos(x,y,J,lat)
                    if(test_flip(x,y,J,lat,T))
                        E = E + energy_pos(x,y,J,lat) - E_0
                    end
                end

                for latindex in 1:4
                    if latindex == 1
                        skyrm_num = skyrmion_number(lat)
                        Mag = total_mag(lat)
                    elseif latindex == 2
                        lat_transform(lat,latindex)
                        skyrm_num = skyrmion_number(lat)
                        Mag = total_mag(lat)
                        lat_transform(lat,latindex)
                    elseif latindex == 3
                        lat_transform(lat,latindex)
                        skyrm_num = skyrmion_number(lat)
                        Mag = total_mag(lat)
                        lat_transform(lat,latindex)
                    elseif latindex == 4
                        lat_transform(lat,latindex)
                        skyrm_num = skyrmion_number(lat)
                        Mag = total_mag(lat)
                        lat_transform(lat,latindex)
                    end

                    skyrm_jack[i,1,latindex] = abs(skyrm_num*normalisation)
                    skyrm_jack[i,2,latindex] = (skyrm_num*normalisation).^2
                    skyrm_jack[i,3,latindex] = (skyrm_num*normalisation).^4

                    M_jack[i,1,latindex] = (norm(Mag)*normalisation)
                    M_jack[i,2,latindex] = (norm(Mag)*normalisation).^2
                    M_jack[i,3,latindex] = (norm(Mag)*normalisation).^4

                end
            end

            for jj in 1:4
                for ii in 1:3
                    if jj == 2
                        skyrm_vec[Tcount,1,ii,jj], skyrm_vec[Tcount,2,ii,jj] = jackknife((skyrm_jack[:,ii,jj]+skyrm_jack[:,ii,jj+1])./2)
                        M_vec[Tcount,1,ii,jj], M_vec[Tcount,2,ii,jj] = jackknife((M_jack[:,ii,jj]+M_jack[:,ii,jj+1])./2)
                    else
                        skyrm_vec[Tcount,1,ii,jj], skyrm_vec[Tcount,2,ii,jj] = jackknife(skyrm_jack[:,ii,jj])
                        M_vec[Tcount,1,ii,jj], M_vec[Tcount,2,ii,jj] = jackknife(M_jack[:,ii,jj])
                    end
                end
                if jj == 2
                    magbind_vec[Tcount,1,jj],magbind_vec[Tcount,2,jj] = bindjack((M_jack[:,3,jj]+M_jack[:,3,jj+1])./2,(M_jack[:,2,jj]+M_jack[:,2,jj+1])./2)
                    skyrmbind_vec[Tcount,1,jj],skyrmbind_vec[Tcount,2,jj] = bindjack((skyrm_jack[:,3,jj]+skyrm_jack[:,3,jj+1])./2,(skyrm_jack[:,2,jj]+skyrm_jack[:,2,jj+1])./2)
                else
                    magbind_vec[Tcount,1,jj],magbind_vec[Tcount,2,jj] = bindjack(M_jack[:,3,jj],M_jack[:,2,jj])
                    skyrmbind_vec[Tcount,1,jj],skyrmbind_vec[Tcount,2,jj] = bindjack(skyrm_jack[:,3,jj],skyrm_jack[:,2,jj])
                end
	    end
            Tcount = Tcount + 1

            #if T == Tmin
            #    autocor_vec = autocor(skyrm_jack)
            #end
        end

        for jj in 1:4
            for ii in 1:3
                Jskyrm_vec[:,Jcount,jj,ii] = skyrm_vec[:,1,ii,jj]
                Jskyrm_vec_err[:,Jcount,jj,ii] = skyrm_vec[:,2,ii,jj]

                JM_vec[:,Jcount,jj,ii] = M_vec[:,1,ii,jj]
                JM_vec_err[:,Jcount,jj,ii] = M_vec[:,2,ii,jj]
            end
            Jmagbind_vec[:,Jcount,jj] = magbind_vec[:,1,jj]
            Jmagbind_vec_err[:,Jcount,jj] = magbind_vec[:,2,jj]

            Jskyrmbind_vec[:,Jcount,jj] = skyrmbind_vec[:,1,jj]
            Jskyrmbind_vec_err[:,Jcount,jj] = skyrmbind_vec[:,2,jj]
        end

        Jcount = Jcount + 1
        println("J_space:",J)
    end

    return Jskyrm_vec,Jskyrm_vec_err,JM_vec,JM_vec_err,Jmagbind_vec,Jmagbind_vec_err,Jskyrmbind_vec,Jskyrmbind_vec_err

end
