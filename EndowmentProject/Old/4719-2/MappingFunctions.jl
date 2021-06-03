#this file contains the functions in Mapping.jl

#the function naming convention is [filetype]_[year(s)], where all refers to all for the type
#examples:
# all_all, pf_all, pf_2003-2009_2011
# refers to all files, all private foundation files,
# and the 2003-2009 and 2011 private foundation files
# the first function is always all_all, which is the only one with no
# dataframe argument

#holds the 1:1 mapping

const MAPPING = Dict{Symbol, Dict{Symbol, Symbol}}(
  :all_all=>Dict(
    :ein=>:ein,
    :name => :name,
    :fisyr=>:fisyr),

  :pf_all=>Dict(
    :p1intrev => :cashinterest,
    :p1divid => :dividends,
    :p1ginvpf => :grprof,
    :p1othinc => :otherinc,
    :p1totrev => :totrev,
    :p1netinv => :returnamt,
    :p1tcont => :contributions,
    :p1contpd => :distributions,
    :p1texmex => :pfexpcharity,
    :p1totexp => :totexp,
    :p2gvtinv => :treasuries,
    :p2crpstk => :stock,
    :p2crpbnd => :bonds,
    :p2totast => :bookassets,
    :p2tliabl => :bookliabilities,
    ),

  :pf_1997_2015=>Dict(
    :p1invrv => :invgross,
    :p1nadjrv => :pfopinc,
    :p1grents => :rents,
    :p1tinvex => :invexpense,
    :p1tadjex => :pfexpincome,
    :p2cashbv => :cash,
    :p2mortg => :mortgages,
    :p2eyotin => :otherexcrealestate,
    :p3eytfnd => :pfnetassets,
    :p10tasst => :totnoncharassets,
    :taxper => :datefiscal
    ),

  :pf_1998_2015=>Dict(
    :ntee1=>:category
    ),

  :pf_1992_2015=>Dict(
    :p10mir => :pfminreturn,
    :p2tasfmv => :marketassets
    ),

  :pf_2011_2015=>Dict(
    :naics=>:naics
    ),

  :pc_all=>Dict(
    :invinc => :totinv,
    :ntee1 => :category,
    :fundbal => :netassets,
    :totrev => :totrev,
    :secur => :secursales,
    :salesexp => :securexp,
    :exps => :totexp,
    :cont => :contributions,
    :grprof => :grprof,
    :ass_eoy => :bookassets,
    :liab_eoy => :bookliabilities,
    :taxper => :datefiscal
    ),

  :pc_1991_2015=>Dict(
    :progrev => :progrev,
    :fundinc => :specrev,
    :othinc => :othinc
    ),

  :pc_1997_2015=>Dict(
    :salesecn => :securrealized,
    :saleothn => :othrealized,
    :netrent => :netrent
    ),

  :pc_2009_2015=>Dict(
    :naics=>:naics
    ),

  :co_all=>Dict(
    :ntee1 => :category
    ),
  :co_1997_2015=>Dict(
    :taxper => :datefiscal
    ),

  :co_2000_2015=>Dict(
    :invinc => :totinv,
    :fundbal => :netassets,
    :salesexp => :securexp,
    :salesecn => :securrealized,
    :saleothn => :othrealized,
    :totrev2 => :totrev,
    :exps => :totexp,
    :cont => :contributions,
    :progrev => :progrev,
    :grprof => :grprof,
    :fundinc => :specrev,
    :othinc => :othinc,
    :netrent => :netrent,
    :ass_eoy => :bookassets,
    :liab_eoy => :bookliabilities,
    ),

  :co_2011_2015=>Dict(
    :naics=>:naics
    ),
  :co_1989_1999=>Dict(
    :p1naseoy => :netassets,
    :grossrec => :totrev

    )

)

replacemissing(v::Union{Number, Missing} , default::Union{Number, Missing} = 0.) = ismissing(v) ? default : v


function all_all(data::NCCSFile)
  local outDF::DataFrame #the output dataframe
  local fields::Dict #contains the name mapping

  try
    fields = MAPPING[:all_all]
    outDF = data.df[collect(keys(fields))]
    rename!(outDF, fields)

    rename!(outDF, :fisyr=>:fisyrold)
    outDF[:fisyr] = ((s::String)->parse(Int, s)).(outDF[:fisyrold])
    deletecols!(outDF, :fisyrold)

    outDF[:fileyear] = data.meta.fileYear
  catch err
    @error "$err \nMetadata: $(data.meta)"
  end


  #println(outDF[1:3, :name])
  return outDF
end

#public charities
function pc_all(data::NCCSFile, outDF::DataFrame)
  fields::Dict = MAPPING[:pc_all]

  #1:1 mappings
  outDF = hcat(outDF, data.df[collect(keys(fields))])
  rename!(outDF, fields)

  outDF[ :charitytype] = fill(:pc, size(outDF, 1))
  #outDF[ :invrealized] =  outDF[ :totinv] .+ outDF[ :secursales] .-  outDF[ :securexp]
  outDF[ :netincome] =  outDF[ :totrev] .- outDF[:totexp]

  return outDF
end

function pc_1991_2015(data::NCCSFile, outDF::DataFrame)
  fields::Dict = MAPPING[:pc_1991_2015]

  #1:1 mappings
  outDF = hcat(outDF, data.df[collect(keys(fields))])
  rename!(outDF, fields)

  return outDF
end

function pc_1997_2015(data::NCCSFile, outDF::DataFrame)
  fields::Dict = MAPPING[:pc_1997_2015]

  #1:1 mappings
  outDF = hcat(outDF, data.df[collect(keys(fields))])
  rename!(outDF, fields)

  #For data quality reasons measure this two different ways
  outDF[:invrealizeddirect] =  outDF[:securrealized] .+  outDF[:othrealized]  .+
    outDF[:totinv] .+ (replacemissing).(outDF[:netrent])
  outDF[:invrealizedindirect] = outDF[:totrev]  .- outDF[:contributions] .- outDF[:grprof]
  outDF[:invrealizedindirect] .-=  (replacemissing).(outDF[:progrev]) .+
    (replacemissing).(outDF[:specrev]) .+ (replacemissing).(outDF[:othinc])


  return outDF
end

function pc_2009_2015(data::NCCSFile, outDF::DataFrame)
  fields::Dict = MAPPING[:pc_2009_2015]

  #1:1 mappings
  outDF = hcat(outDF, data.df[collect(keys(fields))])
  rename!(outDF, fields)

  return outDF
end

function co_all(data::NCCSFile, outDF::DataFrame)
  #get targeted fields
  fields::Dict = MAPPING[:co_all]

  #1:1 mappings
  outDF = hcat(outDF, data.df[collect(keys(fields))])
  rename!(outDF, fields) #make the appropriate names
  outDF[:charitytype] = fill(:co, size(outDF, 1))

  return outDF
end

function co_1997_2015(data::NCCSFile, outDF::DataFrame)
  #get targeted fields
  fields::Dict = MAPPING[:co_1997_2015]

  #1:1 mappings
  outDF = hcat(outDF, data.df[collect(keys(fields))])
  rename!(outDF, fields) #make the appropriate names

  return outDF
end

function co_2000_2015(data::NCCSFile, outDF::DataFrame)
  #get targeted fields
  fields::Dict = MAPPING[:co_2000_2015]

  #1:1 mappings
  outDF = hcat(outDF, data.df[collect(keys(fields))])
  rename!(outDF, fields) #make the appropriate names

  outDF[ :netincome] =  outDF[ :totrev] .- outDF[:totexp]
  #For data quality reasons measure this two different ways
  outDF[:invrealizeddirect] =  outDF[:securrealized] .+  outDF[:othrealized]  .+
    outDF[:totinv] .+ (replacemissing).(outDF[:netrent])
  outDF[:invrealizedindirect] = outDF[:totrev]  .- outDF[:contributions] .- outDF[:grprof]
  outDF[:invrealizedindirect] .-=  (replacemissing).(outDF[:progrev]) .+
    (replacemissing).(outDF[:specrev]) .+ (replacemissing).(outDF[:othinc])


  return outDF
end

function co_2011_2015(data::NCCSFile, outDF::DataFrame)
  #get targeted fields
  fields::Dict = MAPPING[:co_2011_2015]

  #1:1 mappings
  outDF = hcat(outDF, data.df[collect(keys(fields))])
  rename!(outDF, fields) #make the appropriate names

  return outDF
end

function co_1989_1999(data::NCCSFile, outDF::DataFrame)
  #get targeted fields
  fields::Dict = MAPPING[:co_1989_1999]

  #1:1 mappings
  outDF = hcat(outDF, data.df[collect(keys(fields))])
  rename!(outDF, fields) #make the appropriate names

  return outDF
end

#private foundations
function pf_all(data::NCCSFile, outDF::DataFrame)
  fields::Dict = MAPPING[:pf_all]

  #1:1 mappings
  outDF = hcat(outDF, data.df[collect(keys(fields))])
  rename!(outDF, fields)

  #transformations
  outDF[:returncheckamt] =  outDF[:totrev] .- outDF[:contributions] .- outDF[:grprof]

  outDF[:netassetsbook] = outDF[:bookassets] .- outDF[:bookliabilities]
  outDF[:netincome] =  outDF[ :totrev] .- outDF[:totexp]

  outDF[:charitytype] = fill(:pf, size(outDF, 1))

  return outDF
end
function pf_2011_2015(data::NCCSFile, outDF::DataFrame)
  #get targeted fields
  fields::Dict = MAPPING[:pf_2011_2015]

  #1:1 mappings
  outDF = hcat(outDF, data.df[collect(keys(fields))])
  rename!(outDF, fields) #make the appropriate names


  return outDF
end

#private foundation files 1997-2015
function pf_1998_2015(data::NCCSFile, outDF::DataFrame)
  #get targeted fields
  fields::Dict = MAPPING[:pf_1998_2015]

  #1:1 mappings
  outDF = hcat(outDF, data.df[collect(keys(fields))])
  rename!(outDF, fields) #make the appropriate names


  return outDF
end

#private foundation files 1997-2015
function pf_1997_2015(data::NCCSFile, outDF::DataFrame)
  #get targeted fields
  fields::Dict = MAPPING[:pf_1997_2015]

  #1:1 mappings
  outDF = hcat(outDF, data.df[collect(keys(fields))])
  rename!(outDF, fields) #make the appropriate names

  #other transformations
  outDF[:pfgrossfromnet] = outDF[:returnamt] .- outDF[:invexpense] # =outDF[:invgross]

  return outDF
end

#private foundation files 1992-2015
function pf_1992_2015(data::NCCSFile, outDF::DataFrame)
  #get targeted fields
  fields::Dict = MAPPING[:pf_1992_2015]

  #1:1 mappings
  outDF = hcat(outDF, data.df[collect(keys(fields))])
  rename!(outDF, fields) #make the appropriate names
  outDF[:netassets] = outDF[:marketassets] .- outDF[:bookliabilities]

  return outDF
end
