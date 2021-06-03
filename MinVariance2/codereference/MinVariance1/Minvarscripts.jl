using Revise

if  pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end


using Minvar, Random, Revise

#Revise.revise(Minvar)

const WORKERS=2
const CLEAN_WORKERS=false

#set up for parallel operations
initializepar(WORKERS, reset=CLEAN_WORKERS)


Random.seed!(11)
priortype= :uninformative
nburn = 10_000
prefix = priortype == :informative ? "inf_$nburn" : "uninf_$nburn"

#testportfolio(iter=20)
#testcov(iter=2_000)
#testtrisectedportfolio(iter=200, increment=199, displayiter = 0)
#print(mathematica2juliaEq())
#testω(iter=1000)
#testposteriors()
@time E = estimateminvariance(nburn=nburn, nsamples=100_000, parallel=(nburn≠0), recordweights=false,
  saveestimate=true, priortype=priortype, prefix=prefix)
@time makeconvergencegraphs(E, prefix=prefix)

if CLEAN_WORKERS
  cleanup()
end
