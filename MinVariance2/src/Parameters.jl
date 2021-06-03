#the point of this file is to read the parameters as an excel file
#entry point to acquire the parameters
function getparameters!(parameterfile::String=PARAMETER_FILE)

  #this will contain the parameters
  workbook = XLSX.readxlsx("$parameterfile.xlsx")
  (length(PARAM) > 0) && clearcollection(PARAM)

  #iterate across sheets
  for sheet::XLSX.Worksheet ∈ workbook.workbook.sheets
    #special worksheet names that should be ignored
    if startswith("$(sheet.name)", "dep_") || ("$(sheet.name)" ≡ "readme")
      continue
    end

    if startswith("$(sheet.name)", "sub_") #allows for subdictionaries of parameters
      subparam::Dict{Symbol, Any} = Dict{Symbol, Any}()
      subparamname = "$(sheet.name)"[5:end] |> Symbol
      haskey(PARAM, subparamname) && throw(
        "Cannot add subparam dict $subparamname because key already exists in param dictionary")
      PARAM[subparamname] = getparameters(sheet, subparam)
    else
      getparameters(sheet, PARAM)
    end
  end

  close(workbook)

  #keep a log of the parameters for each pass- will want to clean these out periodically
  logparameters()

  #PARAM[] = param

  return nothing
end

#extracts and parses the parameters out of each workbook
function getparameters(sheet::XLSX.Worksheet, param::AbstractDict{Symbol, Any})

  #create a dataframe
  paramdf::DataFrame = XLSX.gettable(sheet) |> DataFrame

  #easier to access the vectors as a dictionary
  paramdfindex = Dict(r.x2=>r.x1 for r ∈ eachrow(paramdf))

  #get the types
  typestrings::Vector{String} = paramdfindex[:type]

  #get the default values
  valuestrings::Vector{String} = paramdfindex[:default]

  #get the labels
  parameterlabels = (Symbol).(paramdfindex[:parameter])

  #now add in the overrides
  for (i, s) ∈ enumerate(paramdfindex[:override])
    if !ismissing(s)
      try
        valuestrings[i] = s
      catch err
        error("Failed to load override. parse failed.
          Label: $(parameterlabels[i]) Type: $(typestrings[i]) Val: $s
          error: $err
          trace:\n$(stacktrace(catch_backtrace()))")
      end
    end
  end



  #now parse the parameters to their correct type
  for (i, s) ∈ enumerate(valuestrings)
    try

      #each parameter needs a unique key
      haskey(param, parameterlabels[i]) && throw(
        "Multiple parameters found with same key.")

      #core parsing logic
      T::Type = parseparam(Type, typestrings[i])
      param[parameterlabels[i]] = parseparam(T, s)

    catch err #throw a meaningful error following failure
      throw("Parameter parse failed.
        Label: $(parameterlabels[i]) Type: $(typestrings[i]) Val: $s
        error: $err")
    end
  end

  return param
end


function logparameters(;logpath = "$(PARAM[:testpath])\\log",
    logname = "cap$(Dates.format(now(),"yyyymmdd_HHMMSS"))")

  cp("$PARAMETER_FILE.xlsx", "$logpath\\MVP_$(logname).xlsx")
  return nothing
end



#parsing scenarios
parseparam(::Type{Nothing}, args...) = nothing
parseparam(::Type{T}, s::String) where T<:Number = something(tryparse(T,s), parseparamalways(T,s))
parseparam(::Type{Type}, s::String) = eval(Meta.parse("""$(s)"""))
parseparam(::Type{String}, s::String) = s
parseparam(::Type{Union{String,Nothing}}, s::String) = s ≡"nothing" ? nothing : s


parseparam(::Type{Date}, i::Int, dateformat::DateFormat) = parseparam(Date, "$i", dateformat)
parseparam(::Type{DateFormat}, s::String) = DateFormat(s)
parseparam(::Type{T}, s::String) where T<:DatePeriod = eval(Meta.parse("""$T($(s))"""))
parseparam(::Type{Date}, s::String, dateformat::DateFormat) = Dates.Date(s, dateformat)

#fallback
parseparam(::Type{T}, s::String) where T = parseparamalways(T,s)
function parseparamalways(::Type{T}, s::String) where T
  x::T = eval(Meta.parse("""($(s))"""))
  return x
end
