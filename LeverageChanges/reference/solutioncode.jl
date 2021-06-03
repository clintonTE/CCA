

#CSV, Serialize, Dates, Satistics, and Dataframes are essential for working with data.
#GZip is for gzip files. The other packages relate to drawing graphs using Gadfly.
using Revise
using  GZip, CSV, DataFrames, Serialization, Dates,Statistics, Cairo, Fontconfig, Gadfly

#define path names as constants
const DATAPATH = "data"
const OUTPUTPATH = "output"
const CRSPNAME = "CRSPData"
const CRSPEVENTNAME = "CRSPDataEvent"


#make sure we can load things using relative paths
if  pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end


function runSolution(;refresh=true, buildsolution=true)::Nothing

  #load the data
  crsp::DataFrame = loadCRSPData("CRSPData", refresh=refresh)
  crspevent::DataFrame = loadCRSPData("CRSPDataEvent", refresh=refresh)

  #examine the data
#  show(describe(crsp), allrows=true, allcols=true)

  buildsolution && buildStudyPanel(crsp, crspevent)
  makegraph()

  return nothing
end

function loadCRSPData(crspname;datapath=DATAPATH, refresh=true)::DataFrame
  local p::String = "$datapath\\$crspname"
  local df::DataFrame

  #makes a serialization object so we don't have to keep parsing the CSV
  if refresh || (!isfile("$p.jls"))
    df = CSV.read("$p.csv", allowmissing=:all, silencewarnings=true)
    serialize("$p.jls", df)
  end

  #read in the dataframe
  df = deserialize("$p.jls")
  sort!(df, [:PERMNO])
  return df
end

function buildStudyPanel(crsp, crspevent)::Nothing

  local dateFormat::Dates.DateFormat = Dates.DateFormat("yyyymmdd")

  #get the total returns of each security
  crsp[:ret] = ((ret,dlret)->
    ismissing(ret) && ismissing(dlret) ? missing :
      ismissing(ret) ? dlret :
      ismissing(dlret) ? ret :
      (1+ret)*(1+dlret) - 1.
    ).(crsp[:RET],crsp[:DLRET])

  #drop rows with missing data
  crsp = crsp[completecases(crsp[[:ret, :SHRCD, :EXCHCD, :sprtrn]]),:]
  crspevent = crspevent[completecases(crspevent),:]

  #uncomment major exchanges only
  #crsp = crsp[((exchcd,shrcd)->
  #  (exchcd ∈ [1,2,3]) && (shrcd ∈ [10,11])).(crsp[:EXCHCD], crsp[:SHRCD]),:]

  #make actual dates
  crsp[:DATE] = (s::String->Dates.Date(s,dateFormat)).((string).(crsp[:date]))
  crspevent[:DATE] = (s::String->Dates.Date(s,dateFormat)).((string).(crspevent[:EXDT]))

  #form excess returns
  crsp[:excess] = crsp[:ret] .- crsp[:sprtrn]

  sort!(crsp,[:PERMNO,:DATE])
  sort!(crspevent, [:PERMNO, :DATE])

  #make fields to contain the event data
  local dayfields::Vector{Symbol} = (d->d<0 ? Symbol("TM",abs(d)) : Symbol("D",d)).(collect(-20:20))
  for f::Symbol ∈ dayfields
    crspevent[f] = Vector{Union{Float64,Missing}}(missing, size(crspevent,1))
  end

  #We want to be careful with the dates- crsp will not adjust the index
  #return for gaps in the holding period returns! (See R solution)
  crspdates::Vector{Date} = unique(crsp[:,:DATE])
  sort!(crspdates)
  validdateranges::Dict =
    Dict(crspdates[i]=>(crspdates[i-20], crspdates[i+20]) for i ∈ 21:(length(crspdates)-20))

  #we will use this to set the boundaries for the dates
  crspevent[:tokeep]=trues(size(crspevent,1))
  crspevent[:begindate] = Vector{Union{Missing,Date}}(missing, size(crspevent,1))
  crspevent[:enddate] = Vector{Union{Missing,Date}}(missing, size(crspevent,1))

  print("\nStarting # of events: $(size(crspevent,1))\n")

  #build a set of crsp events
  subevents::GroupedDataFrame = groupby(crspevent, :PERMNO)

  #create a lookup file mapping permnos to events
  subeventsindex::Dict = Dict((subevents[i])[1,:PERMNO] => i for i ∈ 1:length(subevents))

  permctr::Int = 0
  for subcrsp::SubDataFrame ∈ groupby(crsp, [:PERMNO])
    permctr += 1

    #this could take a while looking at so many sotcks
    if permctr % 1000 == 0
      println("$permctr permnos complete")
    end

    #make a lookup table mapping dates to crsp rows
    dateindex::Dict = Dict(subcrsp[i,:DATE]=>i for i ∈ 1:(size(subcrsp,1)))
    permno::Int = subcrsp[1,:PERMNO]

    #partition out the events for a single stock
    if haskey(subeventsindex, permno)
      subevent::SubDataFrame = subevents[subeventsindex[permno]]

      #make vectors for fast inner-loop access
      excesses::Vector{SubArray{Union{Float64,Missing},1}} =
        (i::Int->subevent[dayfields[i]]).(1:length(dayfields))
      tokeeps::SubArray{Bool,1} = subevent[:tokeep]
      excess::SubArray{Float64,1} = subcrsp[:excess]

      for i::Int ∈ 1:(size(subevent,1))

        local ind::Int #this willbe used as a row pointer later
        #make sure we are not too close to the ends, and the ex-div date is a trading day
        #strictly this may not be necessary, but it does keep a consistent timeline
        eventdate::Date = subevent[i,:DATE]
        if haskey(dateindex, eventdate) && haskey(validdateranges,eventdate)
          ind = dateindex[eventdate]
        else
          tokeeps[i] = false
        end

        #check the endpoints and make sure we capture the entire range
        #also look up the associated ate range against the valid date ranges
        if ((tokeeps[i]) && (ind > 20) && (ind + 20 ≤ size(subcrsp,1)) &&
          (validdateranges[eventdate] == (subcrsp[ind-20, :DATE], subcrsp[ind+20, :DATE])))

          subevent[i, :begindate] = subcrsp[ind-20, :DATE] #set the date boundaries to ensure no overlap
          subevent[i, :enddate] = subcrsp[ind+20, :DATE]

          #innermost loop- tell julia to fully optimize
          @fastmath @inbounds @simd for j ∈ 1:(length(dayfields))
            excesses[j][i] = excess[ind-21+j]
          end


          #now we check the boundaries
          if i > 1
            nocollision::Bool = (subevent[i-1,:DATE]<subevent[i, :begindate])
            tokeeps[i] =tokeeps[i] && nocollision
            tokeeps[i-1] = tokeeps[i-1] && nocollision
          end

          if i < length(tokeeps)
            nocollision = (subevent[i+1,:DATE]>subevent[i, :enddate])
            tokeeps[i] = tokeeps[i] && nocollision
            tokeeps[i+1] = tokeeps[i+1] && nocollision
          end
        else
          tokeeps[i] = false
        end
      end
    end
  end

  crspevent = crspevent[crspevent[:tokeep],:]
  crspevent = crspevent[completecases(crspevent[dayfields]),:]

  serialize("$DATAPATH\\results.jls", crspevent)

  print("Ending # of events: $(size(crspevent,1))\n")
  return nothing
end

function makegraph()::Nothing
  local dayfields::Vector{Symbol} = (d->d<0 ? Symbol("TM",abs(d)) : Symbol("D",d)).(collect(-20:20))

  #----
  crspevent::DataFrame = deserialize("$DATAPATH\\results.jls")

  meanexcess::DataFrame = DataFrame(day=collect(-20:20),
    excess=vec(Statistics.mean(Matrix{Float64}(crspevent[dayfields]),dims=1)).*10_000)


  p = plot(meanexcess, x=:day, y=:excess, Guide.xlabel("Day"), Guide.ylabel("Excess Returns (bp)"),
    Geom.line, Geom.point, Guide.title("Mean Total Returns Over S&P500 Near Ex-Div Dates"))


  draw(PNG("ExcessReturns.png",9inch, 7inch),p)
  return nothing
end

#to run the full code, set both of these to true
@time runSolution(refresh=false, buildsolution=false)
