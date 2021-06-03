using Revise
using Pkg


if  pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end

Pkg.activate(pwd())
using BayesMNIST

sleep(0.4)
@time BayesMNIST.modelmlp()
