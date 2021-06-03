

if pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end

const WORKERS = 8 #or whetever number of parallel routines you want

if nworkers() < WORKERS || (nworkers()==1 && WORKERS == 1)
  addprocs(nworkers()==1?WORKERS-nworkers()+1:WORKERS-nworkers())
end

using F
@everywhere using F

const dim = 10

function mainFunc()::Void
  m::Matrix{Float64} = Matrix{Float64}(dim,dim)  ## a 1000x1000 matrix of floats
  vargs::Vector{Vector{Float64}} = Vector{Vector{Float64}}(dim*dim)
  vans::Vector{Float64} = Vector{Float64}(dim*dim)

  Threads.@threads for i ∈ 1:dim
   for j ∈ 1:dim
       vargs[(i-1)*dim + j] = [i, j]
     end
   end
  vans .= pmap(f, vargs)
  for i ∈ 1:dim
   for j ∈ 1:dim
       m[i,j] = vans[(i-1)*dim + j]
     end
   end

  showall(m)
  return nothing
end

mainFunc()
