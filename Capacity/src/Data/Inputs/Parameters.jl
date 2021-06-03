
#entry point to acquire the parameters
function getparameters!(parameterfile::String=PARAMETER_FILE)

  #this will contain the parameters
  workbook = XLSX.readxlsx(parameterfile)
  (length(PARAM) > 0) && clearcollection(PARAM)

  #iterate across sheets
  for sheet::XLSX.Worksheet ∈ workbook.workbook.sheets
    startswith("$(sheet.name)", "dep_") && continue #allows for depreciaiton of worksheets
    getparameters(sheet, PARAM)
  end

  close(workbook)

  logparameters()

  #PARAM[] = param

  return nothing
end

#extracts and parses the parameters out of each workbook
function getparameters(sheet::XLSX.Worksheet, param::AbstractDict{Symbol, Any})

  #create a dataframe
  paramdf::DataFrame = DataFrame(XLSX.gettable(sheet)...)

  #println(describe(paramdf))

  #easier to access the vectors as a dictionary
  #paramdfindex = Dict(r.x2=>r.x1 for r ∈ eachrow(paramdf))

  #get the types
  typestrings::Vector{String} = paramdf.type

  #get the default values
  valuestrings::Vector{String} = paramdf.default

  #get the labels
  parameterlabels = (Symbol).(paramdf.parameter)

  #now add in the overrides
  for (i, s) ∈ enumerate(paramdf.override)
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
      haskey(param, parameterlabels[i]) && error(
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
end

#try to write out all the parameter info used
#WARNING: constants need to be actively added to this list, so try to avoid them
function formatparameters()
  b = IOBuffer()
  write(b, "Capacity log from time: $(Dates.format(now(),"yyyymmdd_HHMM"))")
  write(b, "\n*****************\n")
  write(b, formatdict(PARAM, title="PARAM dict:"))
  write(b, "\n*****************\n")
  write(b, formatdict(REGRESSION_METHOD, title="REGRESSION_METHOD dict:"))
  write(b, "\n*****************\nConstants: \n")
  write(b, "RANDOM_SEED_VALUE = $RANDOM_SEED_VALUE\n")
  write(b, "PARAMETER_FILE = $(PARAMETER_FILE)\n")
  write(b, "PARALLEL = $(PARALLEL)\n")
  write(b, "LARGE_VAL = $(LARGE_VAL)\n")
  write(b, "LARGE_VAL32 = $(LARGE_VAL32)\n")
  write(b, "SMALL_VAL = $(SMALL_VAL)\n")
  write(b, "SMALL_VAL32 = $(SMALL_VAL32)\n")

  return String(take!(b))
end


function logparameters(;logpath = "$(PARAM[:testpath])\\log",
    logname = "cap-$(Dates.format(now(),"yyyymmdd_HHMM"))")
  paramlog = formatparameters()
  open("$logpath\\$logname.txt", "w+") do f
    write(f, paramlog)
  end

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
