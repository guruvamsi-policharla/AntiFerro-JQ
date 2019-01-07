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
    lat = Array{Vector{Float64},2}(undef,M, N);
    for i = 1:M
        for j = 1:N
            lat[i,j] = sample_uni()
        end
    end
    return lat
end

function myplus(a,b,N)
    c = mod(a+b,N)
    if c == 0
        c = N
    end
    return c
end

function energy_pos(x, y, Q, lat, a = [0,0,0])
    M = size(lat,1)
    N = size(lat,2)

    up = lat[myplus(x,-1,N),y]
    down = lat[myplus(x,1,N),y]
    left = lat[x,myplus(y,-1,N)]
    right = lat[x,myplus(y,1,N)]
    urc = lat[myplus(x,-1,N),myplus(y,1,N)]
    ulc = lat[myplus(x,-1,N),myplus(y,-1,N)]
    lrc = lat[myplus(x,1,N),myplus(y,1,N)]
    llc = lat[myplus(x,1,N),myplus(y,-1,N)]

    energy = -1*(1-dot(lat[x,y],(right+down)))
    energy = energy - Q*( (1-dot(lat[x,y],down))*(1-dot(right,lrc)) + (1-dot(lat[x,y],right))*(1-dot(down,lrc)) )

    return energy
end

function test_flip(x, y, Q, lat, T)
""" Checks whether energy allows for a flip or not """
    N = size(lat,1)
    de = -energy_pos(x,y,Q,lat)-energy_pos(myplus(x,-1,N),y,Q,lat)-energy_pos(myplus(x,-1,N),myplus(y,-1,N),Q,lat)-energy_pos(x,myplus(y,-1,N),Q,lat)
    a0 = lat[x,y]
    a = sample_gauss(lat[x,y])
    lat[x,y] = a
    de = de + energy_pos(x,y,Q,lat)+energy_pos(myplus(x,-1,N),y,Q,lat)+energy_pos(myplus(x,-1,N),myplus(y,-1,N),Q,lat)+energy_pos(x,myplus(y,-1,N),Q,lat)
    if(de<0)
        return true
    elseif(rand()<exp(-de/T))
        return true
    else
        lat[x,y] = a0
        return false
    end
end

function transient_results(lat, transient::Int, Q, T)
    """Takes lat as input and removes initial transients by running for transient number of steps"""
    M = size(lat,1)
    N = size(lat,2)
    for i = 1:transient
        for j = 1:M*N
                x = rand(1:M)
                y = rand(1:N)
                test_flip(x,y,Q,lat,T)
        end
    end
end

function total_mag(lat)
    #returns total magnetisation vector
    return sum(lat)
end

function total_energy(Q,lat)
    e = 0.0
    for i = 1:size(lat,1)
        for j = 1:size(lat,2)
            e = e + energy_pos(i,j,Q,lat)
        end
    end
    return e/2
end

function jackknife(v)
    s = sum(v)
    n = length(v)
    vec_jack = (s .- v)/(n-1)
    jack_avg = mean(vec_jack)
    jack_err = sqrt((mean(vec_jack.^2) .- jack_avg.^2) * (n-1))
    return jack_avg,jack_err
end

function bindjack(vec4,vec2)
    n = length(vec4)
    s4 = sum(vec4)
    s2 = sum(vec2)
    vec_jack = ((s4 .- vec4)./(s2 .- vec2).^2) .* (n-1)

    jack_avg = mean(vec_jack)
    jack_err = sqrt(abs(mean(vec_jack.^2) .- jack_avg.^2) * (n-1))

    #println("Finished calculating Binder")
    return jack_avg,jack_err
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

function four_trans_skyrm(lat)
    M = size(lat,1)
    N = size(lat,2)
    q = zeros(M,N)
    for i in 1:M
       for j in 1:N
           a = lat[i,j]#centre
           b = lat[i,mod(j,N)+1]#right
           c = lat[mod(i,M)+1,mod(j,N)+1]#rightdown
           d = lat[mod(i,M)+1,j]#down
           q[i,j] = (spher_tri_area(a,b,c) + spher_tri_area(a,c,d))/(4*pi)
       end
    end
    ft = abs.(fft(q))
return ft
end

function trans_skyrm(lat)
    M = size(lat,1)
    N = size(lat,2)
    q = zeros(M,N)
    for i in 1:M
       for j in 1:N
           a = lat[i,j]#centre
           b = lat[i,mod(j,N)+1]#right
           c = lat[mod(i,M)+1,mod(j,N)+1]#rightdown
           d = lat[mod(i,M)+1,j]#down
           q[i,j] = (spher_tri_area(a,b,c) + spher_tri_area(a,c,d))/(4*pi)
       end
    end
return q
end

function skyrmion_number(lat)
    M = size(lat,1)
    N = size(lat,2)
    q = 0
    for i in 1:M
        for j in 1:N
            a = lat[i,j]#centre
            b = lat[i,mod(j,N)+1]#right
            c = lat[mod(i,M)+1,mod(j,N)+1]#rightdown
            d = lat[mod(i,M)+1,j]#down
            q = q + (spher_tri_area(a,b,c) + spher_tri_area(a,c,d))/(4*pi)
        end
    end
    return q
end

function spher_tri_area(a,b,c)
    x = cross(a,b)
    y = cross(b,c)
    z = cross(c,a)
    try
        a1 = acos(dot(x,-z)/norm(x)/norm(z))
        a2 = acos(dot(y,-x)/norm(y)/norm(x))
        a3 = acos(dot(z,-y)/norm(z)/norm(y))
        return (a1 + a2 + a3 - pi)*sign(dot(a,cross(b,c)))
    catch err
        if isa(err,DomainError)
            println(dot(x,-z)/norm(x)/norm(-z),dot(y,-x)/norm(y)/norm(-x),dot(z,-y)/norm(z)/norm(-y))
            a1 = acos(round(dot(x,-z)/norm(x)/norm(z)))
            a2 = acos(round(dot(y,-x)/norm(y)/norm(x)))
            a3 = acos(round(dot(z,-y)/norm(z)/norm(y)))
            return (a1 + a2 + a3 - pi)*sign(dot(a,cross(b,c)))
        end
    end
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

function montecarlo(Temperature,N,Q_space)
    mcs = 40000
    M = N

    normalisation=(1.0/float(M*N))

    JM_vec = zeros(length(Temperature),length(Q_space),4,3)
    JM_vec_err = zeros(length(Temperature),length(Q_space),4,3)
    Jskyrm_vec = zeros(length(Temperature),length(Q_space),4,3)
    Jskyrm_vec_err = zeros(length(Temperature),length(Q_space),4,3)
    Jmagbind_vec = zeros(length(Temperature),length(Q_space),4)
    Jskyrmbind_vec = zeros(length(Temperature),length(Q_space),4)
    Jmagbind_vec_err = zeros(length(Temperature),length(Q_space),4)
    Jskyrmbind_vec_err = zeros(length(Temperature),length(Q_space),4)

    M_vec = zeros(length(Temperature),2,3,4)
    M_jack = zeros(mcs,3,4)

    skyrm_vec = zeros(length(Temperature),2,3,4)
    skyrm_jack = zeros(mcs,3,4)

    magbind_vec = zeros(length(Temperature),2,4)
    skyrmbind_vec = zeros(length(Temperature),2,4)

    #autocor_vec = 0
    Qcount = 1
    for Q in Q_space
        lat = initialise(M,N)
        Tcount = 1
        for T in Temperature
            transient_results(lat,3000,Q,T)
            E = total_energy(Q,lat)
            for i in 1:mcs
                for j in 1:M*N
                    x = rand(1:M)
                    y = rand(1:N)
                    E_0 = energy_pos(x,y,Q,lat)
                    if(test_flip(x,y,Q,lat,T))
                        E = E + energy_pos(x,y,Q,lat) - E_0
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
                #qFT[:,:,Qcount] = qFT[:,:,Qcount] + four_trans_skyrm(lat)
            end

            for jj in 1:4
                for ii in 1:3
                    skyrm_vec[Tcount,1,ii,jj], skyrm_vec[Tcount,2,ii,jj] = jackknife(skyrm_jack[:,ii,jj])
                    M_vec[Tcount,1,ii,jj], M_vec[Tcount,2,ii,jj] = jackknife(M_jack[:,ii,jj])
                end
                magbind_vec[Tcount,1,jj],magbind_vec[Tcount,2,jj] = bindjack(M_jack[:,3,jj],M_jack[:,2,jj])
                skyrmbind_vec[Tcount,1,jj],skyrmbind_vec[Tcount,2,jj] = bindjack(skyrm_jack[:,3,jj],skyrm_jack[:,2,jj])
	    end
            Tcount = Tcount + 1
            println(T)
            #if T == Tmin
            #    autocor_vec = autocor(skyrm_jack)
            #end
        end

        for jj in 1:4
            for ii in 1:3
                Jskyrm_vec[:,Qcount,jj,ii] = skyrm_vec[:,1,ii,jj]
                Jskyrm_vec_err[:,Qcount,jj,ii] = skyrm_vec[:,2,ii,jj]

                JM_vec[:,Qcount,jj,ii] = M_vec[:,1,ii,jj]
                JM_vec_err[:,Qcount,jj,ii] = M_vec[:,2,ii,jj]
            end
            Jmagbind_vec[:,Qcount,jj] = magbind_vec[:,1,jj]
            Jmagbind_vec_err[:,Qcount,jj] = magbind_vec[:,2,jj]

            Jskyrmbind_vec[:,Qcount,jj] = skyrmbind_vec[:,1,jj]
            Jskyrmbind_vec_err[:,Qcount,jj] = skyrmbind_vec[:,2,jj]
        end

        Qcount = Qcount + 1
        println("Q_space:",Q)
    end

    return Jskyrm_vec,Jskyrm_vec_err,JM_vec,JM_vec_err,Jmagbind_vec,Jmagbind_vec_err,Jskyrmbind_vec,Jskyrmbind_vec_err
end
