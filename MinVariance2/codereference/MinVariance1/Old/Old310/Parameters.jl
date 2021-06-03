
#holds the main parameter object
mutable struct Parameters
  σ²G::Float64
  ζ²G::Float64
  ζ²P::Float64
  SG::Float64
  SGP::Float64
end

#named constructor for testing purposes
function Parameters(;  σ²G::Float64 = 0.0,  ζ²G::Float64 = 0.0, ζ²P::Float64 = 0.0,
  SG::Float64 = 0.0,  SGP::Float64 = 0.0)

  return Parameters(σ²G, ζ²G, ζ²P, SG, SGP)
end
