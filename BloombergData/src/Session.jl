

function bbsession(bbscripts::Vector{T};
    livepull::Bool = error("livepull must be true or false")) where T<:Val
  S = startbbsession()
  try
    for v ∈ bbscripts
      bbscript(S,v, livepull=livepull)
    end
  finally
    endbbsession(S) #always close the session when done
  end
end

#helper methods to make the above easier to call
bbsession(;kws...) = bbsession(Vector{Val}(); kws...)
bbsession(v::Val;kws...) = bbsession([v]; kws...)
bbsession(bbscripts::Vector{Symbol};kws...) = bbsession([Val{s}() for s ∈ bbscripts]; kws...)
bbsession(s::Symbol;kws...) = bbsession([s]; kws...)


startbbsession() = BLPData.Session(client_mode=BLPData.BLPAPI_CLIENTMODE_SAPI)
endbbsession(S) = BLPData.stop(S)
