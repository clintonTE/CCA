
function combine2tuples(vs::Vector...)
  combined::Vector = vec(collect(Base.product(vs...)))

  return combined
end

function combine2vectors(vs::T...) where T<:Vector
  combinedtuples::Vector = combine2tuples(vs...)

  N::Int = length(combinedtuples)
  M::Int = length(combinedtuples[1])

  #preallocate
  combined::Vector{T} = [T(undef, N) for i ∈ 1:M]

  for c ∈ 1:M, r ∈ 1:N
    combined[c][r] = combinedtuples[r][c]
  end

  return combined
end
