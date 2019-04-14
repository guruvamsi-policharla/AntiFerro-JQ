function jackknife(v)
    s = sum(v)
    n = length(v)
    vec_jack = (s .- v)/(n-1)
    jack_avg = mean(vec_jack)
    jack_err = sqrt(abs(mean(vec_jack.^2) .- jack_avg.^2) * (n-1)) #abs since the error is too small
    return jack_avg,jack_err
end

function bindjack(vec4,vec2)
    n = length(vec4)
    s4 = sum(vec4)
    s2 = sum(vec2)
    vec_jack = (((s2 .- vec2).^2) ./(s4 .- vec4)) ./ (n-1)

    jack_avg = mean(vec_jack)
    jack_err = sqrt(abs(mean(vec_jack.^2) .- jack_avg.^2) * (n-1))

    return jack_avg,jack_err
end
