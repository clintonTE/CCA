using Revise, Pkg, DataFrames


if  pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end

Pkg.activate(pwd())
using BloombergData

#need this since we are uing generated functions]
if @isdefined(FIRST_RUN)
  revise(BloombergData)
end

const FIRST_RUN = true


function BloombergDataScript()
  BloombergData.bbsession(:geoffpull, livepull=true)
end

@time BloombergDataScript()
