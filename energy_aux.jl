function total_mag(lat)
    #returns total magnetisation vector
    return sum(lat)
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

    #energy = -2*(1-dot(lat[x,y],(right+down)))#2 is since we will be dividing by 2 later on

    energy = 2*(1-dot(lat[x,y],(right+down)))#Changing sign of J to check for consistency
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
    de = de/2

    if(de<0)
        return true
    elseif(rand()<exp(-de/T))
        return true
    else
        lat[x,y] = a0
        return false
    end
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
