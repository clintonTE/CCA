#this structure holds the partitioned portfolio
struct TrisectedPortfolio
  P::Portfolio
  partitions::Vector{Portfolio}
  assignments::Vector{Int}

  wₛ::Vector{Float64} #sums of the weights of the partitions
  S::Matrix{Float64} #covariance matrix of the portfolio
end

#constructor
function TrisectedPortfolio(P::Portfolio, start1=1)::TrisectedPortfolio
  targetPartitionSize::Int = P.U.M ÷ 3
  assignments::Vector{Int} = fill(3, P.U.M)

  #do the initial assignments
  assignments[1:targetPartitionSize] .= 1
  assignments[(targetPartitionSize+1):(2*targetPartitionSize)] .= 2

  #create 3 new portfolios, one for each partition
  #the portfolios are a placeholder until the refresh command is called as part of rotatetrisected
  partitions::Vector{Portfolio} = (i::Int->Portfolio(P, P.wₜ)).(1:3)

  #make the trisected portfolio
  TP::TrisectedPortfolio = TrisectedPortfolio(
    P,
    partitions,
    assignments,
    Vector{Float64}(undef,3),
    Matrix{Float64}(undef,3,3)
  )

  #rotate the portfolio, which will also compute wₛ and S
  rotatetrisected!(TP, start1 - 1)

  return TP
end

#computes the covariance of the partitions and the test portfolio
function Statistics.cov(TG::TrisectedPortfolio, P::Portfolio)::Vector{Float64}
  local SGP::Vector{Float64} = Vector{Float64}(undef, 3)

  for i::Int ∈ 1:3
    SGP[i] = cov(TG.partitions[i], P)
  end

  return SGP
end

#rotate the assignments by a set number
function rotatetrisected!(TP::TrisectedPortfolio, by::Int = 1)::TrisectedPortfolio

  #rotate the assignments
  for i::Int ∈ 1:TP.P.U.M
    TP.assignments[i] = TP.assignments[1 + ((i + by - 1) % (TP.P.U.M))]
  end

  #refresh the remaining contents given the changes
  refreshtrisected!(TP)
  return TP
end

function refreshtrisected!(TP::TrisectedPortfolio)::TrisectedPortfolio
  #get the sum of the weights for each parittion
  TP.wₛ .= (i::Int->sum(TP.P.wₜ .* (i .== TP.assignments))).(1:3)

  #update the portfolio partitions
  (i::Int->Portfolio!(TP.partitions[i], TP.P.wₜ .* (i .== TP.assignments), normalizeto=TP.wₛ[i])).(1:3)

  #update the covariance matrix
  cov!(TP.partitions, TP.S)

  return TP
end

#refreshes the trisected portfolio given a new set of target weights
function refreshtrisected!(TP, wₜ)::TrisectedPortfolio
  Portfolio!(TP.P, wₜ)

  return refreshtrisected!(TP)
end

function getω₁C1(Θ::Parameters, TG::TrisectedPortfolio, SGP::Vector{Float64})
  CSGP::Vector{ComplexF64} = (ComplexF64).(SGP) #this is necessary due to potential roundoff error
  CS::Matrix{ComplexF64} = (ComplexF64).(TG.S)
  Cwₛ::Vector{ComplexF64} = (ComplexF64).(TG.wₛ)

  @inbounds ω₁::Float64 = real(
    (CSGP[2]^2*CS[3,3]*Cwₛ[1] - 2*CS[2,3]*CSGP[2]*CSGP[3]*Cwₛ[1] + CS[2,2]*CSGP[3]^2*Cwₛ[1]
    - CSGP[1]*CSGP[2]*CS[3,3]*Cwₛ[2] +CSGP[1]*CS[2,3]*CSGP[3]*Cwₛ[2] - CSGP[2]*CS[3,3]*Θ.SGP*Cwₛ[1]*Cwₛ[2]
    + CS[2,3]*CSGP[3]*Θ.SGP*Cwₛ[1]*Cwₛ[2] + CSGP[1]*CS[3,3]*Θ.SGP*Cwₛ[2]^2 +CSGP[1]*CS[2,3]*CSGP[2]*Cwₛ[3]
    - CSGP[1]*CS[2,2]*CSGP[3]*Cwₛ[3] + CS[2,3]*CSGP[2]*Θ.SGP*Cwₛ[1]*Cwₛ[3] - CS[2,2]*CSGP[3]*Θ.SGP*Cwₛ[1]*Cwₛ[3]
    -2*CSGP[1]*CS[2,3]*Θ.SGP*Cwₛ[2]*Cwₛ[3] + CSGP[1]*CS[2,2]*Θ.SGP*Cwₛ[3]^2 + CS[1,3]*CSGP[2]*(CSGP[3]*Cwₛ[2]
    - CSGP[2]*Cwₛ[3]) -CS[1,2]*CSGP[3]*(CSGP[3]*Cwₛ[2] - CSGP[2]*Cwₛ[3]) - CS[1,3]*Θ.SGP*Cwₛ[2]*(CSGP[3]*Cwₛ[2]
    - CSGP[2]*Cwₛ[3]) + CS[1,2]*Θ.SGP*Cwₛ[3]*(CSGP[3]*Cwₛ[2] - CSGP[2]*Cwₛ[3]) -(CSGP[3]*Cwₛ[2]
    - CSGP[2]*Cwₛ[3])^2* sqrt((1/(CSGP[3]*Cwₛ[2] - CSGP[2]*Cwₛ[3])^2)*((-CS[1,1])*CSGP[2]^2*CS[3,3]
    + 2*CS[1,1]*CS[2,3]*CSGP[2]*CSGP[3] +CS[1,2]^2*CSGP[3]^2 - CS[1,1]*CS[2,2]*CSGP[3]^2
    - 2*CS[1,2]*CSGP[2]*CS[3,3]*Θ.SGP*Cwₛ[1] + 2*CS[1,2]*CS[2,3]*CSGP[3]*Θ.SGP*Cwₛ[1]
    + CSGP[2]^2*CS[3,3]*Θ.SG*Cwₛ[1]^2 -2*CS[2,3]*CSGP[2]*CSGP[3]*Θ.SG*Cwₛ[1]^2 + CS[2,2]*CSGP[3]^2*Θ.SG*Cwₛ[1]^2
    + CS[2,3]^2*Θ.SGP^2*Cwₛ[1]^2 - CS[2,2]*CS[3,3]*Θ.SGP^2*Cwₛ[1]^2 +2*CS[1,1]*CSGP[2]*CS[3,3]*Θ.SGP*Cwₛ[2]
    - 2*CS[1,1]*CS[2,3]*CSGP[3]*Θ.SGP*Cwₛ[2] - 2*CS[1,2]*CSGP[3]^2*Θ.SG*Cwₛ[1]*Cwₛ[2]
    + 2*CS[1,2]*CS[3,3]*Θ.SGP^2*Cwₛ[1]*Cwₛ[2] +CS[1,1]*CSGP[3]^2*Θ.SG*Cwₛ[2]^2 - CS[1,1]*CS[3,3]*Θ.SGP^2*Cwₛ[2]^2
    + CS[1,3]^2*(CSGP[2] - Θ.SGP*Cwₛ[2])^2 - 2*CS[1,1]*CS[2,3]*CSGP[2]*Θ.SGP*Cwₛ[3] -2*CS[1,2]^2*CSGP[3]*Θ.SGP*Cwₛ[3]
    + 2*CS[1,1]*CS[2,2]*CSGP[3]*Θ.SGP*Cwₛ[3] + 2*CS[1,2]*CSGP[2]*CSGP[3]*Θ.SG*Cwₛ[1]*Cwₛ[3]
    - 2*CS[1,2]*CS[2,3]*Θ.SGP^2*Cwₛ[1]*Cwₛ[3] -2*CS[1,1]*CSGP[2]*CSGP[3]*Θ.SG*Cwₛ[2]*Cwₛ[3]
    + 2*CS[1,1]*CS[2,3]*Θ.SGP^2*Cwₛ[2]*Cwₛ[3] + CS[1,1]*CSGP[2]^2*Θ.SG*Cwₛ[3]^2 + CS[1,2]^2*Θ.SGP^2*Cwₛ[3]^2
    -CS[1,1]*CS[2,2]*Θ.SGP^2*Cwₛ[3]^2 + CSGP[1]^2*(CS[2,3]^2 - CS[2,2]*CS[3,3] + CS[3,3]*Θ.SG*Cwₛ[2]^2
    - 2*CS[2,3]*Θ.SG*Cwₛ[2]*Cwₛ[3] + CS[2,2]*Θ.SG*Cwₛ[3]^2) -2*CSGP[1]*(CS[1,3]*((-CS[2,2])*CSGP[3]
    + CSGP[3]*Θ.SG*Cwₛ[2]^2 + CS[2,3]*(CSGP[2] - Θ.SGP*Cwₛ[2]) + CS[2,2]*Θ.SGP*Cwₛ[3]
    - CSGP[2]*Θ.SG*Cwₛ[2]*Cwₛ[3]) +CS[1,2]*((-CSGP[2])*CS[3,3] + CS[2,3]*CSGP[3] + CS[3,3]*Θ.SGP*Cwₛ[2]
    - CS[2,3]*Θ.SGP*Cwₛ[3] - CSGP[3]*Θ.SG*Cwₛ[2]*Cwₛ[3] + CSGP[2]*Θ.SG*Cwₛ[3]^2) +Cwₛ[1]*(CS[2,3]^2*Θ.SGP
    - CS[2,2]*CS[3,3]*Θ.SGP + CSGP[2]*CS[3,3]*Θ.SG*Cwₛ[2] + CS[2,2]*CSGP[3]*Θ.SG*Cwₛ[3]
    - CS[2,3]*Θ.SG*(CSGP[3]*Cwₛ[2] + CSGP[2]*Cwₛ[3]))) -2*CS[1,3]*(CS[1,2]*(CSGP[2] -
    Θ.SGP*Cwₛ[2])*(CSGP[3] - Θ.SGP*Cwₛ[3]) + Cwₛ[1]*(CS[2,3]*Θ.SGP*(-CSGP[2] + Θ.SGP*Cwₛ[2])
    +CSGP[2]*Θ.SG*((-CSGP[3])*Cwₛ[2] + CSGP[2]*Cwₛ[3]) + CS[2,2]*Θ.SGP*(CSGP[3] - Θ.SGP*Cwₛ[3]))))))/
    (CS[2,2]*(CSGP[3]*Cwₛ[1] - CSGP[1]*Cwₛ[3])^2 +CSGP[2]^2*(CS[3,3]*Cwₛ[1]^2 + Cwₛ[3]*(-2*CS[1,3]*Cwₛ[1]
    + CS[1,1]*Cwₛ[3])) - 2*CSGP[2]*(CSGP[1]*CS[3,3]*Cwₛ[1]*Cwₛ[2] - CS[1,3]*CSGP[3]*Cwₛ[1]*Cwₛ[2]
    -CS[1,2]*CSGP[3]*Cwₛ[1]*Cwₛ[3] - CS[1,3]*CSGP[1]*Cwₛ[2]*Cwₛ[3] + CS[1,1]*CSGP[3]*Cwₛ[2]*Cwₛ[3]
    + CS[1,2]*CSGP[1]*Cwₛ[3]^2 + CS[2,3]*Cwₛ[1]*(CSGP[3]*Cwₛ[1] - CSGP[1]*Cwₛ[3])) +Cwₛ[2]*(CSGP[3]^2*(-2*CS[1,2]*Cwₛ[1]
    + CS[1,1]*Cwₛ[2]) + 2*CSGP[1]*CSGP[3]*(CS[2,3]*Cwₛ[1] - CS[1,3]*Cwₛ[2] + CS[1,2]*Cwₛ[3])
    +CSGP[1]^2*(CS[3,3]*Cwₛ[2] - 2*CS[2,3]*Cwₛ[3])))
    )

  return ω₁
end

function getω₁C2(Θ::Parameters, TG::TrisectedPortfolio, SGP::Vector{Float64})
  CSGP::Vector{ComplexF64} = (ComplexF64).(SGP) #this is necessary due to potential roundoff error
  CS::Matrix{ComplexF64} = (ComplexF64).(TG.S)
  Cwₛ::Vector{ComplexF64} = (ComplexF64).(TG.wₛ)

  @inbounds ω₁::Float64 = real(
    (CSGP[2]^2*CS[3,3]*Cwₛ[1] - 2*CS[2,3]*CSGP[2]*CSGP[3]*Cwₛ[1] + CS[2,2]*CSGP[3]^2*Cwₛ[1]
    - CSGP[1]*CSGP[2]*CS[3,3]*Cwₛ[2] +CSGP[1]*CS[2,3]*CSGP[3]*Cwₛ[2] - CSGP[2]*CS[3,3]*Θ.SGP*Cwₛ[1]*Cwₛ[2]
    + CS[2,3]*CSGP[3]*Θ.SGP*Cwₛ[1]*Cwₛ[2] + CSGP[1]*CS[3,3]*Θ.SGP*Cwₛ[2]^2 +CSGP[1]*CS[2,3]*CSGP[2]*Cwₛ[3]
    - CSGP[1]*CS[2,2]*CSGP[3]*Cwₛ[3] + CS[2,3]*CSGP[2]*Θ.SGP*Cwₛ[1]*Cwₛ[3] - CS[2,2]*CSGP[3]*Θ.SGP*Cwₛ[1]*Cwₛ[3]
    -2*CSGP[1]*CS[2,3]*Θ.SGP*Cwₛ[2]*Cwₛ[3] + CSGP[1]*CS[2,2]*Θ.SGP*Cwₛ[3]^2 + CS[1,3]*CSGP[2]*(CSGP[3]*Cwₛ[2]
    - CSGP[2]*Cwₛ[3]) - CS[1,2]*CSGP[3]*(CSGP[3]*Cwₛ[2] -CSGP[2]*Cwₛ[3]) - CS[1,3]*Θ.SGP*Cwₛ[2]*(CSGP[3]*Cwₛ[2]
    - CSGP[2]*Cwₛ[3]) + CS[1,2]*Θ.SGP*Cwₛ[3]*(CSGP[3]*Cwₛ[2] - CSGP[2]*Cwₛ[3]) +(CSGP[3]*Cwₛ[2]
    - CSGP[2]*Cwₛ[3])^2* sqrt((1/(CSGP[3]*Cwₛ[2] - CSGP[2]*Cwₛ[3])^2)*((-CS[1,1])*CSGP[2]^2*CS[3,3]
    + 2*CS[1,1]*CS[2,3]*CSGP[2]*CSGP[3] +CS[1,2]^2*CSGP[3]^2 - CS[1,1]*CS[2,2]*CSGP[3]^2
    - 2*CS[1,2]*CSGP[2]*CS[3,3]*Θ.SGP*Cwₛ[1] + 2*CS[1,2]*CS[2,3]*CSGP[3]*Θ.SGP*Cwₛ[1]
    + CSGP[2]^2*CS[3,3]*Θ.SG*Cwₛ[1]^2 -2*CS[2,3]*CSGP[2]*CSGP[3]*Θ.SG*Cwₛ[1]^2 + CS[2,2]*CSGP[3]^2*Θ.SG*Cwₛ[1]^2
    + CS[2,3]^2*Θ.SGP^2*Cwₛ[1]^2 - CS[2,2]*CS[3,3]*Θ.SGP^2*Cwₛ[1]^2 +2*CS[1,1]*CSGP[2]*CS[3,3]*Θ.SGP*Cwₛ[2]
    - 2*CS[1,1]*CS[2,3]*CSGP[3]*Θ.SGP*Cwₛ[2] - 2*CS[1,2]*CSGP[3]^2*Θ.SG*Cwₛ[1]*Cwₛ[2]
    + 2*CS[1,2]*CS[3,3]*Θ.SGP^2*Cwₛ[1]*Cwₛ[2] +CS[1,1]*CSGP[3]^2*Θ.SG*Cwₛ[2]^2 - CS[1,1]*CS[3,3]*Θ.SGP^2*Cwₛ[2]^2
    + CS[1,3]^2*(CSGP[2] - Θ.SGP*Cwₛ[2])^2 - 2*CS[1,1]*CS[2,3]*CSGP[2]*Θ.SGP*Cwₛ[3] -2*CS[1,2]^2*CSGP[3]*Θ.SGP*Cwₛ[3]
    + 2*CS[1,1]*CS[2,2]*CSGP[3]*Θ.SGP*Cwₛ[3] + 2*CS[1,2]*CSGP[2]*CSGP[3]*Θ.SG*Cwₛ[1]*Cwₛ[3]
    - 2*CS[1,2]*CS[2,3]*Θ.SGP^2*Cwₛ[1]*Cwₛ[3] -2*CS[1,1]*CSGP[2]*CSGP[3]*Θ.SG*Cwₛ[2]*Cwₛ[3]
    + 2*CS[1,1]*CS[2,3]*Θ.SGP^2*Cwₛ[2]*Cwₛ[3] + CS[1,1]*CSGP[2]^2*Θ.SG*Cwₛ[3]^2 + CS[1,2]^2*Θ.SGP^2*Cwₛ[3]^2
    -CS[1,1]*CS[2,2]*Θ.SGP^2*Cwₛ[3]^2 + CSGP[1]^2*(CS[2,3]^2 - CS[2,2]*CS[3,3] + CS[3,3]*Θ.SG*Cwₛ[2]^2
    - 2*CS[2,3]*Θ.SG*Cwₛ[2]*Cwₛ[3] + CS[2,2]*Θ.SG*Cwₛ[3]^2) -2*CSGP[1]*(CS[1,3]*((-CS[2,2])*CSGP[3]
    + CSGP[3]*Θ.SG*Cwₛ[2]^2 + CS[2,3]*(CSGP[2] - Θ.SGP*Cwₛ[2]) + CS[2,2]*Θ.SGP*Cwₛ[3]
    - CSGP[2]*Θ.SG*Cwₛ[2]*Cwₛ[3]) +CS[1,2]*((-CSGP[2])*CS[3,3] + CS[2,3]*CSGP[3] + CS[3,3]*Θ.SGP*Cwₛ[2]
    - CS[2,3]*Θ.SGP*Cwₛ[3] - CSGP[3]*Θ.SG*Cwₛ[2]*Cwₛ[3] + CSGP[2]*Θ.SG*Cwₛ[3]^2) +Cwₛ[1]*(CS[2,3]^2*Θ.SGP
    - CS[2,2]*CS[3,3]*Θ.SGP + CSGP[2]*CS[3,3]*Θ.SG*Cwₛ[2] + CS[2,2]*CSGP[3]*Θ.SG*Cwₛ[3]
    - CS[2,3]*Θ.SG*(CSGP[3]*Cwₛ[2] + CSGP[2]*Cwₛ[3]))) -2*CS[1,3]*(CS[1,2]*(CSGP[2] -
    Θ.SGP*Cwₛ[2])*(CSGP[3] - Θ.SGP*Cwₛ[3]) + Cwₛ[1]*(CS[2,3]*Θ.SGP*(-CSGP[2] + Θ.SGP*Cwₛ[2])
    +CSGP[2]*Θ.SG*((-CSGP[3])*Cwₛ[2] + CSGP[2]*Cwₛ[3]) + CS[2,2]*Θ.SGP*(CSGP[3] - Θ.SGP*Cwₛ[3]))))))/
    (CS[2,2]*(CSGP[3]*Cwₛ[1] - CSGP[1]*Cwₛ[3])^2 +CSGP[2]^2*(CS[3,3]*Cwₛ[1]^2 + Cwₛ[3]*(-2*CS[1,3]*Cwₛ[1]
    + CS[1,1]*Cwₛ[3])) - 2*CSGP[2]*(CSGP[1]*CS[3,3]*Cwₛ[1]*Cwₛ[2] - CS[1,3]*CSGP[3]*Cwₛ[1]*Cwₛ[2]
    -CS[1,2]*CSGP[3]*Cwₛ[1]*Cwₛ[3] - CS[1,3]*CSGP[1]*Cwₛ[2]*Cwₛ[3] + CS[1,1]*CSGP[3]*Cwₛ[2]*Cwₛ[3]
    + CS[1,2]*CSGP[1]*Cwₛ[3]^2 +CS[2,3]*Cwₛ[1]*(CSGP[3]*Cwₛ[1] - CSGP[1]*Cwₛ[3])) + Cwₛ[2]*(CSGP[3]^2*(-2*CS[1,2]*Cwₛ[1]
    + CS[1,1]*Cwₛ[2]) + 2*CSGP[1]*CSGP[3]*(CS[2,3]*Cwₛ[1] -CS[1,3]*Cwₛ[2] + CS[1,2]*Cwₛ[3])
    + CSGP[1]^2*(CS[3,3]*Cwₛ[2] - 2*CS[2,3]*Cwₛ[3])))
    )

  return ω₁
end

function getω₂(Θ::Parameters, TG::TrisectedPortfolio, SGP::Vector{Float64}, ω₁::Float64)::Float64
  ω₂::Float64 = (
    -((-SGP[3] + SGP[3] * ω₁ * TG.wₛ[1] + Θ.SGP * TG.wₛ[3] - SGP[1] * ω₁ * TG.wₛ[3]) /
    (SGP[3] * TG.wₛ[2] - SGP[2] * TG.wₛ[3]))
  )

  return ω₂
end

function getω₃(Θ::Parameters, TG::TrisectedPortfolio, SGP::Vector{Float64}, ω₁::Float64)::Float64
  ω₃::Float64 = (
    -((-SGP[2] + SGP[2] * ω₁ * TG.wₛ[1] + Θ.SGP * TG.wₛ[2] - SGP[1] * ω₁ * TG.wₛ[2]) /
    (-SGP[3] * TG.wₛ[2] + SGP[2] * TG.wₛ[3]))
  )

  return ω₃
end


#builds the scaling vector ω
#pick the solution with the lowest error
function getωpickbest(Θ::Parameters, TG::TrisectedPortfolio, SGP::Vector{Float64})::Vector{Float64}

  ωcandidate1::Vector{Float64} = Vector{Float64}(undef, 3)
  ωcandidate2::Vector{Float64} = Vector{Float64}(undef, 3)
  ωcandidate1[1]::Float64 = getω₁C1(Θ, TG, SGP)
  ωcandidate2[1]::Float64 = getω₁C2(Θ, TG, SGP)



  #fill the rest of the candidate vectors
  ωcandidate1[2] = getω₂(Θ, TG, SGP, ωcandidate1[1])
  ωcandidate1[3] = getω₃(Θ, TG, SGP, ωcandidate1[1])

  ωcandidate2[2] = getω₂(Θ, TG, SGP, ωcandidate2[1])
  ωcandidate2[3] = getω₃(Θ, TG, SGP, ωcandidate2[1])

  errcandidate1::Float64 = abs(ωcandidate1' * TG.S * ωcandidate1 - Θ.SG)#sum((ωcandidate1 .- 1.0).^2)
  errcandidate2::Float64 = abs(ωcandidate2' * TG.S * ωcandidate2 - Θ.SG)#sum((ωcandidate2 .- 1.0).^2)

  #pick the vector with the lowest error
  ω::Vector{Float64} = errcandidate1 ≤ errcandidate2 ? ωcandidate1 : ωcandidate2

  return ω
end

function getω(Θ::Parameters, TG::TrisectedPortfolio, SGP::Vector{Float64})::Vector{Float64}

  ω::Vector{Float64} = Vector{Float64}(undef, 3)
  ω[1]::Float64 = getω₁C2(Θ, TG, SGP)
  ω[2] = getω₂(Θ, TG, SGP, ω[1])
  ω[3] = getω₃(Θ, TG, SGP, ω[1])

  return ω
end


getω(Θ::Parameters, TG::TrisectedPortfolio, P::Portfolio)::Vector{Float64} =
  getω(Θ, TG, cov(TG, P))

function rescaletrisected!(TG::TrisectedPortfolio, ω::Vector{Float64})::TrisectedPortfolio
  wₜ::Vector{Float64} = TG.P.wₜ

  for i ∈ 1:length(wₜ)
    wₜ[i] *= ω[TG.assignments[i]]
  end

  return refreshtrisected!(TG, wₜ)
end
