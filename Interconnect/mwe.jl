using LinearAlgebra
function mwe()
  local Q::Matrix{Float64}

  #Q = Matrix{Float64}(undef, 800, 800)
  X = rand(200,13)
  f = qr(X)
  @time Q = Matrix(f.Q[:,:])

  X = rand(400,13)
  f = qr(X)
  @time Q = Matrix(f.Q[:,:])

  X = rand(600,13)
  f = qr(X)
  @time Q = Matrix(f.Q[:,:])

  X = rand(800,13)
  f = qr(X)
  @time Q = Matrix(f.Q[:,:])
end

@time mwe()
