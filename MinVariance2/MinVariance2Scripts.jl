using Revise, Pkg, BenchmarkTools


if  pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end

Pkg.activate(pwd())
using MinVariance2

#need this since we are uing generated functions
if @isdefined(FIRST_RUN)
  revise(MinVariance2)
end

const FIRST_RUN = true

function minvarscript()
  MinVariance2.getparameters!()
end

@time minvarscript()
