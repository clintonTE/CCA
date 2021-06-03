#using Pkg
#Pkg.activate(pwd())

if  pwd() ∉ LOAD_PATH
    push!(LOAD_PATH,pwd())
end

using Revise, mwe

mwe.func()