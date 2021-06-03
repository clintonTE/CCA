module F
export f

function f(arg::Vector{Float64})
    return sum(arg)
end

end
