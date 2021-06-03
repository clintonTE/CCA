if  pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end

if "$(pwd())\\data" ∉ LOAD_PATH
  push!(LOAD_PATH,"$(pwd())\\data")
end


using CSV, DataFrames, CT, Random

const OUTPUT_PATH = pwd() * "\\output" #path where we store the data
const FOOTER_NAME = "footer.tex" #these are just headers for the tex tables
const HEADER_NAME = "header.tex"

function mwe()
  local df::DataFrame
  local b::IOBuffer

  iStream = open("mwe.csv")
  df = CSV.read(iStream)
  println("This is using open: $(df[:NAME])")
  close(iStream)

  iStream = open("mwe.csv")
  b = IOBuffer(read(iStream))
  df = CSV.read(b)
  println("This uses an IOBuffer: $(df[:NAME])")
  close(iStream)

  df = CSV.read("mwe.csv")
  println("This is direct: $(df[:NAME])")
end

#mwe()

function toyRegression()


  Random.seed!(11)
  years::Vector{Int} = collect(2001:2010) #10 years of data
  post::Vector{Int} = [zeros(Int, 5); ones(Int, 5)]
  groups::Vector{Int} = collect(1:100) #100 groups

  yearGroupFE::Vector{Float64} = # a fixed effect for each group and year with a time trend
    rand(length(groups) * length(years)) .*collect(1:(length(groups) * length(years))) ./ 100

  treated::Vector{Int} = [zeros(Int, 60); ones(Int, 40)]

  individualsPerGroup = 10
  individualFE::Vector{Float64} = rand(individualsPerGroup*length(groups))
  treatmentEffect = 2.0




  observationsN::Int = individualsPerGroup*length(groups)*length(years)
  individuals::Vector{Int} = collect(1:(individualsPerGroup*length(groups)))
  individualsN::Int = length(individuals)
  confounder::Vector{Float64} = 1.0 .* ( ones(observationsN) .* collect(1:observationsN) ./ 1000).^2 #time trend confounder
  ε::Vector{Float64} = rand(observationsN)

  irrelevantEffect::Vector{Float64} =  1.0 .* (-1.0 .* rand(observationsN) .* collect(1:observationsN) ./ 2000)
  irrelevant::Vector{Int} = [fill((i -> (i≤4)).(1:individualsPerGroup), length(groups)*length(years))...;]
  irrelevantEffect .*= irrelevant

  #this will hold the information
  df = DataFrame(N = collect(1:observationsN),
    ε = ε,
    outcome=Vector{MFloat64}(undef, observationsN),
    confounder=confounder,
    years=Vector{MSymbol}(undef, observationsN),
    groups=Vector{MSymbol}(undef, observationsN),
    individuals=Vector{MSymbol}(undef, observationsN),
    post=Vector{MInt}(undef, observationsN),
    irrelevant=irrelevant,
    treated=Vector{MInt}(undef, observationsN),
    keep = trues(observationsN)
    )

  ctr::Int = 0 #now construct the simulated sample
  yearGroupCtr::Int = 0
  for y::Int ∈ 1:length(years)
    for g::Int ∈ 1:length(groups)
      yearGroupCtr += 1
      for i::Int ∈ 1:individualsPerGroup
        individual::Int = (g-1)*individualsPerGroup + i #this is the index of the individual
        ctr+=1
        df[ctr, :years] = Symbol(years[y])
        df[ctr, :groups] = Symbol(groups[g])
        df[ctr, :individuals] = Symbol(individuals[individual])
        df[ctr, :post] = post[y]
        df[ctr, :treated] = treated[g]

        #determine if these are the folks we need to exclude in the post sample
        df[ctr, :keep] = !((irrelevant[ctr] == 1) && (post[y] == 1))

        #the actual process
        df[ctr, :outcome] = individualFE[individual] + yearGroupFE[yearGroupCtr] + confounder[ctr] +
          treated[g]*post[y]*treatmentEffect + irrelevantEffect[ctr] + ε[ctr]
      end
    end
  end

  df = df[df[:keep] .== true,:]
  df[:postXtreated]  = df[:treated] .* df[:post] #create interaction terms
  df[:yearXgroup] = ((y::Symbol,g::Symbol)->Symbol(y, "_", g)).(df[:years],df[:groups]) #make an interacted variable

  categorical!(df, :years) #make everything categorical
  categorical!(df, :groups)
  categorical!(df, :individuals)
  categorical!(df, :yearXgroup)



  models::Vector{CTLM} = Vector{CTLM}()
  #XNames::Vector{Vector{Symbol}} = Vector{Vector{Symbol}}() # holds the names of the x variables

  focal::Vector{Symbol} = [:intercept, :post, :postXtreated]
  push!(models, CTLM(df, Meta.parse("post + postXtreated + individuals"), :outcome,
    XNames=focal, YName=:outcome))

  push!(models, CTLM(df[df[:irrelevant] .≠ 1,:], Meta.parse("post + postXtreated + individuals"), :outcome,
    XNames=focal, YName=:outcome))

  focal = [:intercept, :post, :postXtreated, :confounder]
  push!(models, CTLM(df, Meta.parse("post + postXtreated + confounder + individuals"), :outcome,
    XNames=focal, YName=:outcome))

  push!(models, CTLM(df[df[:irrelevant] .≠ 1,:], Meta.parse("post + postXtreated + confounder + individuals"), :outcome,
    XNames=focal, YName=:outcome))

  displayed::Vector{Symbol} = [:post, :postXtreated, :confounder]

  tableInTex::String = tableText::String = texTable(models,
    getModWhiteΣ!,
    displayed,
    titleCaption = "Comparison",
    colNames = [["Just FE", "Dropped", "Just FE + confounder", "Dropped + confounder"]],
    contentRowNames = (string).(displayed),
    #descRowNames = descRowNames,
    #descContent = descContent,
    decimalDigits = 4,
    #stars=true,
    #starStrings = OVERRIDE_STAR_STRINGS,
    #clearMem = USE_AGGRESSIVE_GC,
    caption = "to be written")

  writeTables2File([tableInTex],
      HEADER_NAME, FOOTER_NAME, path=OUTPUT_PATH,
      outName = "comparisonTable.tex")

  #push!(models, CTLM(longDF, "",  :outcome, XNames=XNames[i], YName = YSpecs[i]))
end

@time toyRegression()
