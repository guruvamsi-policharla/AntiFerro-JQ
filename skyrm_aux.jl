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
