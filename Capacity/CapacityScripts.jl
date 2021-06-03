using Revise, Pkg, DataFrames, Finometrics, BenchmarkTools


if  pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end

Pkg.activate(pwd())
using Capacity

#need this since we are uing generated functions
if @isdefined(FIRST_RUN)
  revise(Capacity)
end

const FIRST_RUN = true
#=throw("1) Determine which spec is better for the tables (see tables0 vs tables)
2) verify the SEs for the AR(1) approach- NOTE- verify the time series are no longer stationary
- maybe put in a check to prevent this from occuring in thef uture? like an INF if the coef is ≥1?
3) Convert the multiple solution simulation graphs into slides
4) Outline the anlytical solutions")=#
#throw("Need to run the funds with refresh
#also, need to do full refresh runs to account for the error in 6/9 month
#momentum")
#throw("Reset the size acceptance threshold to 20% (or smaller)")
#=throw( "- why are MF assets higher than total assets? Doesn't seem right- maybe use the corridor method,
 or try changing up the factors?
   -Try the new domestic only factor
   -Determine if we really need time FE, or an intercept, or nothing. Update the deck accordingly")=#

#throw("now run the new and imrpoved AR(1)")

@warn("
check LARS todo
check timing of shocks to momentum assets

  Why are g and ng producing the same results? repkicate and test
")

@warn("
 - check the mutual fund file and make sure the betas and assets are ok- consider trimming the dates
 via the global setting if needed, may need to aquire new mf file if necessary. Check and
   compare the winsorized results. If it doesn't make a difference, may want to disable
   or shortcircuit the winsorization
 -compute standard errors on G (weighted abs std deviation betas * A / (total assets))
  - then maybe use the delta method to compute the SE of the ratio
 - try the mutual fund assets approach, either just the growth or the bootstrap approach
 - try OOS regressions
 - consider alternate measures of participation
 - check correlation of D/P ratio w/ the lmt measure
 - long run resutls are suspect due to the persistence of the log measure
 - longer horizons- try running the monthyl data on 35/60 years and see what happens
 - short run results seem pretty good- note using G instead of aG makes sense as
 prices will rise with demand
 - the rise in prices is tough to capture because contemperaneous outflows may
 coincide with crisis periods that often have high momentum returns




")
@info("
# Important note: CUDA v2.6.3 runs fastest, v3.x.x runs slow so far.
The easiest way to transition between version is to, in the main environment,
switch to the later version.
The results seem basically equivelent if allowed to run to completion,
  although in one instance I didn't get convergence with the new a version
wrt to the reason for the gap, I tried 1) checking the broadcast
  2) the regression. Implicity, the projection isn't responsible for the
  performance regression either since the gap still exists in the
  level version. Its not limit either.
  It is stable, at leas twith the old version.
  Note the latest run is microscopicalyl different from a bit back,
  but the difference is so tiny w/ correlation 1.0, beta ≈ 1.0. More improtantly,
  the answer is stable as written with a full pass.")
@info "WISHLIST (in rough priority)
REPLICATION and VERIFICATION!!! Mutual funds or momentum?
Mutual funds- some decent correlations. Could try some randomization for placebo test
Ultimately however, if the momentum emasure works well, that is strong evidence that
outweighs even placebo tests
"
#=NOTE OLD @info "I recced this against Lewellen 2015. Returns in Lewellen are a bit lower, and the b/m metric
sd is a bit high. Size and investment (using inv4) are ok. N seems out of control, but who knows
how they are calculating that."=#

#=@ NOTE OLD @info "MCMC- Seems like A₀ is mostly fixed- probably opportunities for performance improvements
  via profiling as well as additional verificaiton that the differences between the rep
  and actual fall with sample size"=#
#@ NOTE OLD @info "a GMM style EIV model will be much more efficient (check rossi slides)"
#note- the uncontrolled version greatly increases power and correlation

function capacityscript()
  #Capacity.gz2zstd("data\\crsp\\crspd35y")
  #Capacity.gz2zstd("data\\crsp\\crspd30y")
  #Capacity.gz2zstd("data\\crsp\\crspd23y")
  #Capacity.gz2zstd("data\\crsp\\crspd20y")
  #Capacity.gz2zstd("data\\crsp\\crspd10y")
  #Capacity.gz2zstd("data\\preliminary\\prelimcwarspm")5
  #Capacity.gz2zstd("data\\crsp\\crspd2y")
  #Capacity.gz2zstd("data\\comp\\compq")
  #Capacity.gz2zstd("data\\ibes\\ibes")
  #Capacity.gz2zstd("data\\hfr\\hfralive")
  #Capacity.gz2zstd("data\\hfr\\hfrdead")


  Capacity.getparameters!() #WARNING- mutates parameter dictionary to match excel parameter file
  #Capacity.formpreliminaryfunds()z
  #Capacity.formpreliminarybeta()


  #Capacity.computemeasure()
  Capacity.analyze()

  #Capacity.capacityio()

  #Capacity.comomentumtable()
  #Capacity.comomentumgraphs()
  #Capacity.fundgraphs()


end

@time capacityscript()
