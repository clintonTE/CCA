




function cleanup()::Nothing
  pToClean=procs()
  for i in 2:length(pToClean)
    rmprocs(pToClean[i])
  end

  return nothing
end

#adds processes as desired
function initializePar(workers::Int = WORKERS; reset::Bool=false, fixRandom::Bool = true)::Nothing
  #if the reset flag is active we kill all active processes besides the main
  #note this takes a while
  if reset == true
    cleanup()
  end

  #create the processes
  if nworkers() < workers || (nworkers()==1 && workers == 1)
    addprocs(nworkers()==1 ? workers-nworkers()+1 : workers-nworkers())
  end

  @everywhere eval( :(if  pwd() ∉ LOAD_PATH
                push!(LOAD_PATH,pwd())
              end))


  #make sure all processes have access to everything udrf
  #@everywhere eval( :(using FLMethods, FLOptimizer, FLData))
  #@everywhere eval( :(using FLOptimizerPar))
  #@everywhere eval( :(using FLScripts))
  @everywhere eval( :(using Endowment))

  return nothing
end
