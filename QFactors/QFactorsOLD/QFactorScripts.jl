using Revise, Pkg

#if  pwd() ∉ LOAD_PATH
#  push!(LOAD_PATH,pwd())
#end

Pkg.develop(PackageSpec(path=pwd()))

import QFactors


#Pkg.activate(pwd())







#revise(QFactors)

#Pkg.develop(PackageSpec(path=pwd()))

#QFactors.recompressqfactorfiles()
#QFactors.recompressqfactorfile("COMP-A")

@time begin
  QFactors.prepcrsp(refreshcrsp=false)
  #QFactors.prepcomp(refreshcomp=false)
end
