using Revise
module TaxData

#this is useful if we have mdoules
if pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end

using CSV, DataFrames, Distributions, StatsBase, GZip#, Gadfly#, Cairo, Fontconfig
using  CTCopy, DataSeries, Formatting,  Dates, Serialization, Gadfly

#=NOTE##############
Some companies have multiple PERMNOs. The returns are different, so I am leaving them as is.
=####################
const WORKING_PATH = pwd() * "\\data" #path where we store the data
const DATA_PATH = pwd()
const OUTPUT_PATH = pwd() * "\\output"
const DATA_FILE_NAME = "histtaxes" #data file name
const DATE_FORMAT_STR = "yyyymmdd"

const FOOTER_NAME = "footer.tex" #these are just headers for the tex tables
const HEADER_NAME = "header.tex"
const YAHOO_FILES = ["vwesx", "vfiix"] # LT IG Corp,GNMA

#useful type
const CTSym = Union{Symbol,Nothing}
const PlotContainer = Union{Plot,Gadfly.Compose.Context}

#DEFAULT_OTHER_SERIES = "gfd_2y_gov_yield_avg + gfd_10y_gov_yield_avg + gfd_20y_gov_yield_avg + corp_spread"
#DEFAULT_OTHER_SERIES_SYMS = [:gfd_2y_gov_yield_avg, :gfd_10y_gov_yield_avg, :gfd_20y_gov_yield_avg, :corp_spread]
const DEFAULT_OTHER_SERIES = "gfd_10_20_t_spread_avg + inflation"
const DEFAULT_OTHER_SERIES_SYMS = [:gfd_10_20_t_spread_avg, :inflation]
const DEFAULT_OTHER_SERIES_NAMES = ["10 Minus 20-Year T-Bond", "CPI Inflation"]
const DEFAULT_DEPENDENT = :muni_spread
const DEFAULT_MUNI_FIELD =  :gfd_20y_muni_yield_avg
const OVERRIDE_STAR_STRINGS = ["\\ostar","\\dstar","\\tstar"]

#=const STATE_MARG_FIELD = :mean_tot_marg
const STATE_MAX_FIELD = :ca_us_wage_max_rate
const STATE_GAINS_FIELD = :ca_us_lt_gains_max_rate
const STATE_INT_FIELD =  :ca_us_wage_max_rate=#


#this function reads in the csv files and writes out a julia serial file
#in: nothing, but reads in the datafiles
#out: nothing, but writes the tables
function writeTaxDataToDF()::Nothing
  oStream::IOStream = open("$WORKING_PATH\\$DATA_FILE_NAME.jls", "w")

  #=serialize(oStream, [readtable("$WORKING_PATH\\$DATA_FILE_NAME-m_pre.csv"),
    readtable("$WORKING_PATH\\$DATA_FILE_NAME-q_pre.csv"),
    readtable("$WORKING_PATH\\$DATA_FILE_NAME-y_pre.csv")])=#
  serialize(oStream, [CSV.read("$WORKING_PATH\\$DATA_FILE_NAME-m_pre.csv"),
    CSV.read("$WORKING_PATH\\$DATA_FILE_NAME-q_pre.csv"),
    CSV.read("$WORKING_PATH\\$DATA_FILE_NAME-y_pre.csv")])

  close(oStream)

  return nothing
end

function scaleTaxData!(taxDFm, taxDFq, taxDFy)::NTuple{3,DataFrame}
    #first scale the data
    #need to convert the following fields to decimal
    for df::DataFrame ∈  [taxDFm, taxDFq, taxDFy]

      #assume any types are numeric
      for c ∈ 1:size(df,2)
        if eltype(df[c])==Any || eltype(df[c])==Union{Any,Missing}
          df[c] = Vector{Union{Float64,Missing}}(df[c])
        end
      #  println(eltype(df[c]))
      end

      for s ∈ [:corp_industrial_yield,
        :corp_industrial_yield_avg, :corp_industrial_yield_close, :corp_max_rate, :corp_t_note_spread,
        :corp_t_note_spread_avg, :corp_t_note_spread_close, :cpi_inflation, :cpi_inflation_avg,
        :cpi_inflation_close, :gfd_1y_gov_yield_avg, :gfd_1y_gov_yield_close,
        :gfd_10y_corp_yield_avg, :gfd_10y_corp_yield_close, :gfd_10y_gov_yield_avg,
        :gfd_10y_gov_yield_close, :gfd_20y_gov_yield_avg, :gfd_20y_gov_yield_close,
        :gfd_20y_muni_yield_avg, :gfd_20y_muni_yield_close, :gfd_20y_muni_broad_yield_avg,
        :gfd_20y_muni_broad_yield_close, :gfd_2y_gov_yield_avg,
        :gfd_2y_gov_yield_close, :gfd_lt_corp_yield_avg, :gfd_lt_corp_yield_close,
        :gfd_st_gov_yield_avg, :gfd_st_gov_yield_close, :muni_go_20_yield_avg,
        :muni_go_20_yield_close, :muni_high_rating_yield, :muni_high_rating_yield_avg,
        :muni_high_rating_yield_close, :muni_highest_rating_yield, :muni_highest_rating_yield_avg,
        :muni_highest_rating_yield_close, :t_3mo_yield, :t_3mo_yield_avg, :t_3mo_yield_close,
        :t_note_10y_rate, :t_note_10y_rate_avg, :t_note_10y_rate_close, :t_st_yield,
        :t_st_yield_avg, :t_st_yield_close, :top_div_rate, :top_lt_gains_rate, :us_cash_gifts_marg_rate,
        :us_dividend_marg_rate, :us_interest_marg_rate, :us_interest_marg_rate_10,
        :us_interest_marg_rate_100, :us_interest_marg_rate_1000, :us_interest_marg_rate_20,
        :us_interest_marg_rate_200, :us_interest_marg_rate_40, :us_interest_marg_rate_400,
        :us_interest_marg_rate_5, :us_lt_gain_marg_rate, :us_lt_gains_marg_rate_10,
        :us_lt_gains_marg_rate_100, :us_lt_gains_marg_rate_1000, :us_lt_gains_marg_rate_20,
        :us_lt_gains_marg_rate_200, :us_lt_gains_marg_rate_40, :us_lt_gains_marg_rate_400,
        :us_lt_gains_marg_rate_5, :us_max_rate, :us_mortage_interest_marg_rate,
        :us_pensions_marg_rate, :us_property_tax_marg_rate, :us_propietors_marg_rate,
        :us_qualified_dividends_marg_rate, :us_st_gain_marg_rate, :us_wage_marg_rate,
        :us_wage_marg_rate_10, :us_wage_marg_rate_100, :us_wage_marg_rate_1000,
        :us_wage_marg_rate_20, :us_wage_marg_rate_200, :us_wage_marg_rate_40,
        :us_wage_marg_rate_400, :us_wage_marg_rate_5, :us_st_wage_marg_rate,
        :us_st_interest_marg_rate,	:us_st_dividend_marg_rate, :us_st_qualified_dividends_marg_rate,
        :us_st_st_gain_marg_rate,	:us_st_lt_gain_marg_rate	,:us_st_pensions_marg_rate,
        :us_st_propietors_marg_rate,	:us_st_property_tax_marg_rate,
        :us_st_mortage_interest_marg_rate,	:us_st_cash_gifts_marg_rate,
        :tips_breakeven_inflation_avg, :tips_breakeven_inflation_close]

        namesOfCols::Vector{Symbol} = names(df)

        if (s ∈ namesOfCols)
          checkVal::Float64 = abs(mean(skipmissing(df[s])))
          if checkVal < 0.5 && checkVal >.0001
            println("WARNING: $s may not require scaling but is being scaled (avg: $checkVal)")
          end
          if eltype(df[s]) ∈ [Int, Union{Int64, Missings.Missing}]
            df[s] /= 100.0
          else
            df[s] ./= 100.0 #convert capital gains rate to a decimal
          end
        end
      end

    end

  return taxDFm, taxDFq, taxDFy
end


function fuseDataSeries(series1::Vector{T}, series2::Vector{U};
  viaFunction::Function = pick2)::Vector where {T<:Union{<:Real,Missing}, U<:Union{<:Real,Missing}}

  #pre-allocate the output
  out::Vector = similar(series1)
  #loop through the series and take non-NA values where possible
  for i ∈ 1:length(out)
    if ismissing(series1[i])
      out[i] = series2[i]
    elseif ismissing(series2[i])
      out[i] = series1[i]
    else
      out[i] = viaFunction(series1[i],series2[i]) #resolve conflicts with viaFunction
    end
  end

  return out
end

function fuseDataSeries(df::DataFrame, series1Sym::Symbol, series2Sym::Symbol;
    viaFunction::Function = pick2)

    if eltype(df[series1Sym]) == Any || eltype(df[series1Sym]) == Union{Any,Missing}
      T1::Type = Union{Float64,Missing}
    else
      T1 = eltype(df[series1Sym])
    end

    if eltype(df[series2Sym]) == Any || eltype(df[series2Sym]) == Union{Any,Missing}
      T2::Type = Union{Float64,Missing}
    else
      T2 = eltype(df[series2Sym])
    end

  return fuseDataSeries(Vector{T1}(df[series1Sym]),Vector{T2}(df[series2Sym]),
    viaFunction = viaFunction)::Vector{T1}
end

function processYahooData!(taxDFm, taxDFq, taxDFy)::Nothing

  yahooNames::Vector{String} =  YAHOO_FILES

  for f::String ∈ yahooNames
    yahooDF::DataFrame = getYahooData(WORKING_PATH, f)

    #process the yahoo series into monthly
    taxDFm[ Symbol(f)] = dailyToEOM(taxDFm[:date], yahooDF[:DATE], yahooDF[ :Adj_Close])
    taxDFm[Symbol(f,"_hpr")] = similar(taxDFm[ Symbol(f)])
    taxDFm[2:end,Symbol(f,"_hpr")] = (taxDFm[2:end,Symbol(f)] ./ taxDFm[1:(end-1),Symbol(f)]) .- 1.0

    taxDFq[ Symbol(f)] = dailyToEOQ(taxDFq[:date], yahooDF[:DATE], yahooDF[ :Adj_Close])
    taxDFq[Symbol(f,"_hpr")] = similar(taxDFq[ Symbol(f)])
    taxDFq[2:end,Symbol(f,"_hpr")] = (taxDFq[2:end,Symbol(f)] ./ taxDFq[1:(end-1),Symbol(f)]) .- 1.0

    taxDFy[ Symbol(f)] = dailyToEOY(taxDFy[:year], yahooDF[:DATE], yahooDF[ :Adj_Close])
    taxDFy[Symbol(f,"_hpr")] = similar(taxDFy[ Symbol(f)])
    taxDFy[2:end,Symbol(f,"_hpr")] = (taxDFy[2:end,Symbol(f)] ./ taxDFy[1:(end-1),Symbol(f)]) .- 1.0

  end

  return nothing

end


#this function gets an array of muni-tax values
#in: muni yields, the risk free rate (As a vector), the credit spread
#out: a datavector of muni tax rates
#NOTE: Assume a tax rate of (1-muniYield/(rfr + spread))
function getMuniTaxRate(muniYield::Vector{Union{Float64, Missing}}, riskFree::Vector{Union{Float64, Missing}},
  creditSpread::Vector)::Vector{Union{Float64, Missing}}

  #pre-allocate
  taxRates::Vector{Union{Float64, Missing}} = similar(muniYield)

  #calculate the implied tax rate
  taxRates .= 1.0 .- muniYield ./ (riskFree .+ creditSpread)

  return taxRates
end

#helper function which processes symbol/dataframe inputs
getMuniTaxRate(df::DataFrame, muniYieldSymbol::Symbol, riskFreeSymbol::Symbol,
  creditSpreadSymbol::Symbol)::Vector{Union{Float64, Missing}} =
    getMuniTaxRate(df[muniYieldSymbol], df[riskFreeSymbol], df[creditSpreadSymbol])

#generic functions for merging data if both data points are populated
pick1(x1::T, x2::T)  where T<:Union{Int,Float64} = x1::T
pick2(x1::T, x2::T) where T<:Union{Int,Float64} = x2::T
pickMid(x1::T, x2::T) where T<:Union{Int,Float64} = ((x1 + x2)/T(2))::T

#this funciton exists purely for organizational reasons, and houses the script for
# aggregating the data series relating to muni tax rates and calculating these tax rates
# I chose not to create a generic function for this for code clarity and future flexibiltiy
#IN: Data frames for the monthly, quarterly, and yearly files
#OUT: Modifies the dataframes in place to include the tax rates and intermediate columns
# The data frames are also returned as a 3-tuple
function processMuniData!(taxDFm::DataFrame, taxDFq::DataFrame, taxDFy::DataFrame;
  viaFunction::Function = pickMid)::NTuple{3,DataFrame}

  ###############First tackle the yearly file
  #first build the muni rate series from the fred data
  taxDFy[:zeroSpread] = 0.0
  taxDFy[:fred_muni_yield_close] =
    fuseDataSeries(taxDFy, :muni_high_rating_yield_close, :muni_go_20_yield_close, viaFunction = viaFunction)

  taxDFy[:fred_muni_yield_avg] =
    fuseDataSeries(taxDFy, :muni_high_rating_yield_avg, :muni_go_20_yield_avg, viaFunction = viaFunction)

  #combine the fred corporate yields and corporate spreads with the gfd treasury yields
  taxDFy[:fred_spread_close] = fuseDataSeries(Vector{Union{Float64, Missing}}(
    taxDFy[:corp_industrial_yield_close] .-
    taxDFy[:gfd_10y_gov_yield_close]), taxDFy[:corp_t_note_spread_close], viaFunction=viaFunction)

  taxDFy[:fred_spread_avg] = fuseDataSeries(Vector{Union{Float64, Missing}}(taxDFy[:corp_industrial_yield_avg] .-
    taxDFy[:gfd_10y_gov_yield_avg]), taxDFy[:corp_t_note_spread_avg], viaFunction=viaFunction)

  #we need to make the spreads from the corporate rates
  taxDFy[:gfd_credit_spread_close] =
    taxDFy[:gfd_lt_corp_yield_close] .- taxDFy[:gfd_20y_gov_yield_close]

  #note: the gfd_lt_corp_yield_avg is of duration 20-30 years
  taxDFy[:gfd_credit_spread_avg] =
    taxDFy[:gfd_lt_corp_yield_avg] .- taxDFy[:gfd_20y_gov_yield_avg]

  #define the series that we will create
  muniSymbols::Vector{Symbol} =
    [:fred_muni_yield_close, :fred_muni_yield_avg, :gfd_20y_muni_yield_close,
      :gfd_20y_muni_yield_avg, :gfd_20y_muni_yield_close, :gfd_20y_muni_yield_avg]
  spreadSymbols::Vector{Symbol} =
    #[:fred_spread_close, :fred_spread_avg, :gfd_credit_spread_close,
    [:gfd_credit_spread_close, :gfd_credit_spread_avg,  :gfd_credit_spread_close,
      :gfd_credit_spread_avg,  :gfd_credit_spread_close, :gfd_credit_spread_avg]
  rfrSymbols::Vector{Symbol} =
    [:gfd_10y_gov_yield_close, :gfd_10y_gov_yield_avg, :gfd_10y_gov_yield_close,
      :gfd_10y_gov_yield_avg, :gfd_20y_gov_yield_close, :gfd_20y_gov_yield_avg]
  muniTaxSymbols::Vector{Symbol} =
    [:fred_muni_tax_close, :fred_muni_tax_avg, :gfd_muni_tax_close,
      :gfd_muni_tax_avg, :gfd_20_muni_tax_close, :gfd_20_muni_tax_avg]

  #now calculate the muni implied tax rates
  for i::Int ∈ 1: length(muniTaxSymbols)
    taxDFy[muniTaxSymbols[i]] = getMuniTaxRate(taxDFy, muniSymbols[i], rfrSymbols[i], spreadSymbols[i])
    taxDFy[Symbol(muniTaxSymbols[i],"_0Credit")] = getMuniTaxRate(taxDFy, muniSymbols[i], rfrSymbols[i], :zeroSpread)
  end

  #####the quarterly file is similar to the  yearly file for the muni series
  #first build the muni rate series from the fred data
  taxDFq[:zeroSpread] = 0.0

  taxDFq[:fred_muni_yield_close] =
    fuseDataSeries(taxDFq, :muni_high_rating_yield_close, :muni_go_20_yield_close, viaFunction = viaFunction)

  taxDFq[:fred_muni_yield_avg] =
    fuseDataSeries(taxDFq, :muni_high_rating_yield_avg, :muni_go_20_yield_avg, viaFunction = viaFunction)

    #combine the fred corporate yields and corporate spreads with the gfd treasury yields
  taxDFq[:fred_spread_close] = fuseDataSeries(Vector{Union{Float64, Missing}}(taxDFq[:corp_industrial_yield_close] .-
    taxDFq[:gfd_10y_gov_yield_close]), taxDFq[:corp_t_note_spread_close], viaFunction=viaFunction)

  taxDFq[:fred_spread_avg] = fuseDataSeries(Vector{Union{Float64, Missing}}(taxDFq[:corp_industrial_yield_avg] .-
    taxDFq[:gfd_10y_gov_yield_avg]), taxDFq[:corp_t_note_spread_avg], viaFunction=viaFunction)

  #we need to make the spreads from the corporate rates
  taxDFq[:gfd_credit_spread_close] =
    taxDFq[:gfd_lt_corp_yield_close] .- taxDFq[:gfd_20y_gov_yield_close]

  taxDFq[:gfd_credit_spread_avg] =
    taxDFq[:gfd_lt_corp_yield_avg] .- taxDFq[:gfd_20y_gov_yield_avg]

    #now calculate the muni implied tax rates
    for i::Int ∈ 1: length(muniTaxSymbols)
      taxDFq[muniTaxSymbols[i]] = getMuniTaxRate(taxDFq, muniSymbols[i], rfrSymbols[i], spreadSymbols[i])
      taxDFq[Symbol(muniTaxSymbols[i],"_0Credit")] = getMuniTaxRate(taxDFq, muniSymbols[i], rfrSymbols[i], :zeroSpread)
    end

  #######the monthly file has more signficant differences
  # Some of the fields are of monthly resolution, so the dif between close and average is meaningless
  # These fields include :muni_high_rating_yield,  :corp_t_note_spread,
  # and :corp_industrial_yield

  taxDFm[:zeroSpread] = 0.0

  taxDFm[:fred_muni_yield_close] =
    fuseDataSeries(taxDFm, :muni_high_rating_yield, :muni_go_20_yield_close, viaFunction = viaFunction)

  taxDFm[:fred_muni_yield_avg] =
    fuseDataSeries(taxDFm, :muni_high_rating_yield, :muni_go_20_yield_avg, viaFunction = viaFunction)

    #combine the fred corporate yields and corporate spreads with the gfd treasury yields
  taxDFm[:fred_spread_close] = fuseDataSeries(Vector{Union{Float64, Missing}}(taxDFm[:corp_industrial_yield] .-
    taxDFm[:gfd_10y_gov_yield_close]), taxDFm[:corp_t_note_spread], viaFunction=viaFunction)

  taxDFm[:fred_spread_avg] = fuseDataSeries(Vector{Union{Float64, Missing}}(taxDFm[:corp_industrial_yield] .-
    taxDFm[:gfd_10y_gov_yield_avg]), taxDFm[:corp_t_note_spread], viaFunction=viaFunction)

  #we need to make the spreads from the corporate rates
  taxDFm[:gfd_credit_spread_close] =
    taxDFm[:gfd_lt_corp_yield_close] .- taxDFm[:gfd_20y_gov_yield_close]

  taxDFm[:gfd_credit_spread_avg] =
    taxDFm[:gfd_lt_corp_yield_avg] .- taxDFm[:gfd_20y_gov_yield_avg]

  #now calculate the muni implied tax rates
  for i::Int ∈ 1: length(muniTaxSymbols)
    taxDFm[muniTaxSymbols[i]] = getMuniTaxRate(taxDFm, muniSymbols[i], rfrSymbols[i], spreadSymbols[i])
    taxDFm[Symbol(muniTaxSymbols[i],"_0Credit")] = getMuniTaxRate(taxDFm, muniSymbols[i], rfrSymbols[i], :zeroSpread)
  end

  return taxDFm, taxDFq, taxDFy

end


#this is a formulaic method for getting the personal tax rate from the capital gains rate
getCumexPersonalTax(egStats::Vector{Union{Float64, Missing}},
  capitalGains::Vector{Union{Float64, Missing}})::Vector{Union{Float64, Missing}} =
  1.0 .- egStats .* (1.0 .- capitalGains)

#helper method which uses symbols instead of datavectors
getCumexPersonalTax(df::DataFrame, egStatsSym::Symbol,
  capitalGainsSym::Symbol)::Vector{Union{Float64, Missing}} =
  getCumexPersonalTax(df[egStatsSym], df[capitalGainsSym])

#this is a formulaic method for getting the personal tax rate from the capital gains rate
getCumexTotalTax(egStats::Vector{Union{Float64, Missing}},
  capitalGains::Vector{Union{Float64, Missing}})::Vector{Union{Float64, Missing}} =
  1.0 .- egStats

#helper method which uses symbols instead of datavectors
getCumexTotalTax(df::DataFrame, egStatsSym::Symbol,
  capitalGainsSym::Symbol)::Vector{Union{Float64, Missing}} =
  getCumexTotalTax(df[egStatsSym], df[capitalGainsSym])

function processCumExData!(taxDFm::DataFrame, taxDFq::DataFrame, taxDFy::DataFrame)::NTuple{3,DataFrame}

  #this was used to flatten the rate data, but the info was updated a bit- hence data pulled
  # using the below commented code would be slightly off
  #=taxDFm[:top_lt_gains_rate] =
    flattenYears(taxDFy[:year], taxDFy[:top_lt_gains_rate], taxDFm[:date])
  taxDFq[:top_lt_gains_rate] =
    flattenYears(taxDFy[:year], taxDFy[:top_lt_gains_rate], taxDFq[:date])=#


  capGainsSymbol::Symbol = :top_lt_gains_rate

  #get the ordinary tax rate from each cumex dividend time series
  for df::DataFrame ∈  [taxDFm, taxDFq, taxDFy]
    df[:cumex_ew_ord_tax] = getCumexPersonalTax(df, :eg_cumex_equal_weights, capGainsSymbol)
    df[:cumex_div_ord_tax] = getCumexPersonalTax(df, :eg_cumex_div, capGainsSymbol)
    df[:cumex_marketcap_ord_tax] = getCumexPersonalTax(df, :eg_cumex_marketcap, capGainsSymbol)

    df[:cumex_ew_tot_tax] = getCumexTotalTax(df, :eg_cumex_equal_weights, capGainsSymbol)
    df[:cumex_div_tot_tax] = getCumexTotalTax(df, :eg_cumex_div, capGainsSymbol)
    df[:cumex_marketcap_tot_tax] = getCumexTotalTax(df, :eg_cumex_marketcap, capGainsSymbol)

  end

  return taxDFm, taxDFq, taxDFy
end

#This funciton gets ex-post inflation over an arbitrary period
function getExPostFromIndex(curIndex::Vector{Union{Float64, Missing}},
    periodsExPost::Int, periodsPerYear::Int = 1)::Vector{Union{Float64, Missing}}

  N::Int = length(curIndex)
  outVec::Vector{Union{Float64, Missing}} = similar(curIndex)
  outVec .= missing

  for i::Int ∈ 1:(N-periodsExPost)
    outVec[i] =
      (curIndex[i+periodsExPost] / curIndex[i])^(Float64(periodsPerYear) / Float64(periodsExPost)) - 1.0
  end

  return outVec
end

#This funciton gets the geometric average of the index return from T-2 to T+1
function get3YSmoothedFromIndex(curIndex::Vector{Union{Float64, Missing}},
    beginOffset::Int = - 2, endOffset::Int = 1,
    periodsPerYear::Int = 1)::Vector{Union{Float64, Missing}}


  curIndexNoNA = collect(skipmissing(curIndex))
  N::Int = length(curIndexNoNA)

  outVecNoNA::Vector{Union{Float64, Missing}} = deepcopy(curIndexNoNA)
  outVecNoNA .= missing

  for i::Int ∈ 1:N
    beginIndex = max(1, i + beginOffset)
    endIndex = min(N, i + endOffset)
    outVecNoNA[i] =
      (curIndexNoNA[endIndex] / curIndexNoNA[beginIndex])^(Float64(periodsPerYear) /
      Float64(endIndex-beginIndex)) - 1.0
  end

  outVec::Vector{Union{Float64, Missing}} = deepcopy(curIndex)
  outVec .= missing
  ctr = 1
  for i::Int ∈ 1:length(curIndex)
    if !ismissing(curIndex[i])
      outVec[i] = outVecNoNA[ctr]
      ctr += 1
    end
  end


  #println(curIndex)
  #println(outVec)

  return outVec
end

#This funciton gets realized inflation over an arbitrary period
function getRealizedFromIndex(curIndex::Vector{Union{Float64, Missing}},
    periodsRealized::Int, periodsPerYear::Int = 1)::Vector{Union{Float64, Missing}}

  N::Int = length(curIndex)
  outVec::Vector{Union{Float64, Missing}} = similar(curIndex)
  outVec .= missing

  for i::Int ∈ ((1+periodsRealized):N)
    outVec[i] =
      (curIndex[i] / curIndex[i-periodsRealized])^(Float64(periodsPerYear) / Float64(periodsRealized)) - 1.0
  end

  return outVec
end

sigRound(x::Float64, maxPlaces::Int)::Float64 = round(signif(x,figures),maxPlaces)

formNum(x::Float64, decimals::Int; scaleHurdle::Float64 = 10.)::String =
  (x≥scaleHurdle) ? "$(Int(round(x)))" : format("{:.$(decimals)f}",x)

formPercentile(x::Float64, CDFFunc, decimals::Int)::String = formNum(CDFFunc(x), decimals)

function writeSelectedDataTable(df::DataFrame;
      fields::Vector{Symbol}=[:munitax, :munitax0credit, :interesttax, :topord, :tnote1, :tnote1r, :tnote1at, :tnote20, :tnote20r, :tnote20at],
  fieldNames::Vector{String}=["munitax", "0credit", "interesttax", "topord", "Nominal", "Real", "Txd-Real", "Nominal", "Real", "Txd-Real"],
  decimals::Int = 1)

  #table width
  K::Int = length(fields)
  N::Int = size(df, 1)

  colNames::Vector{Vector{String}} =
    [["Tax Rates", "1-Year T-Note", "20-Year T-Bond"], fieldNames]

  widthColNames::Vector{Vector{Int}} = [[4,3,3],ones(Int,K)]

  #get the content
  rows::Vector{Vector{String}} = [Vector{String}(undef, K) for r::Int ∈ 1:N]
  for r::Int ∈ 1:N
    for c::Int ∈ 1:K
      val::T where T<:Union{Float64,Missing} = df[r, fields[c]]
      rows[r][c] = ismissing(val) ? "" : "$(formNum(val*100, decimals, scaleHurdle=9999.))"
    end
  end

  rowNames::Vector{String} = (string).(df[:year])

  tableText::String = texTable("Tax Rate Data",
      "See tex file", #caption
      colNames, #colNames
      Vector{String}(),#contentRowNames
      Vector{Matrix{String}}(), #content
      rowNames, #descRowNames
      rows, #descContent
      Vector{String}(), #notes
      widthColNames = widthColNames)

  writeTables2File([tableText], HEADER_NAME, FOOTER_NAME,
    path=OUTPUT_PATH, outName = "selectedData.tex")
end

#this creates a summary table
function writeSummaryTable(df::DataFrame, fields::Vector{Symbol}, fieldNames::Vector{String};
  percentiles::Vector{Float64} = [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0],
  decimals::Int=3, σsAbove::Float64 = 2.0, σsAboveYear::Int = 2016,
  outName::String = "TableSummary.tex")::String

  #label the columns
  numRows::Int = length(fieldNames)
  yearIndex::Int = findfirst(isequal(σsAboveYear),df[:year])


  #label the rows
  summaryColNames::Vector{String} = ["N"; "mean"; "\$\\sigma\$"]
  nonPercentileNumCols = length(summaryColNames)
  summaryColNames = [summaryColNames; ["$(Int(p*100)) \\%" for p::Float64 ∈ percentiles];
    "\$\\%_{$σsAboveYear}\$"; "\$>x_{$σsAboveYear}+$(Int(σsAbove))\\sigma\$"]
  numCols::Int = length(summaryColNames)

  #make the content matrix
  summaryContent::Vector{Vector{String}} = [Vector{String}(undef, numCols) for i::Int ∈ 1:numRows]


  for r::Int ∈ 1:numRows
    summaryContent[r][1] = string(size(df,1)-sum(ismissing.(df[fields[r]])))
    summaryContent[r][2] = formNum(mean(skipmissing(df[fields[r]])),decimals)
    summaryContent[r][3] = formNum(std(skipmissing(df[fields[r]])),decimals)

    #get the percentiles
    for c ∈ 1:length(percentiles)
      summaryContent[r][c+nonPercentileNumCols] =
        formNum(quantile(Vector{Float64}(collect(skipmissing(df[fields[r]]))), percentiles[c]),decimals)
    end


    summaryContent[r][end-1] = ismissing(df[yearIndex,fields[r]]) ? "" :
        formPercentile(df[yearIndex,fields[r]],
          StatsBase.ecdf(Vector{Float64}(collect(skipmissing(df[fields[r]])))), decimals)

    numGreater::T where T <: Union{Int, Missing} = (sum(collect(skipmissing(df[fields[r]])) .≥
      df[df[:year].==σsAboveYear,fields[r]] .+ std(skipmissing(df[fields[r]])) * σsAbove))

    summaryContent[r][end] = ismissing(numGreater) ? "" : (
      numGreater==0 ? "\\text{none}" : string(numGreater))

  end

  TableSummarySD::String = texTable( "Summary of ST and LT Interest Yields",
    """See Tex File""", #caption
    [["Summary Statistics", "Percentiles"], summaryColNames], #colNames
    Vector{String}(),#contentRowNames
    Vector{Matrix{String}}(), #content
    fieldNames, #descRowNames
    summaryContent, #descContent
    Vector{String}(),
    #columnSepPt = -100, #NOTE its possible we will want these back. Look at the old IO version to use it.
    widthColNames = [[nonPercentileNumCols, length(percentiles), length(summaryColNames)-length(percentiles)-nonPercentileNumCols],
      ones(Int, length(summaryColNames))],
    alignmentColNames = [["c", "c"], ["r" for i::Int ∈ 1:length(summaryColNames)]],
    lineSpacer = "\\\\"
    )

    #println(df[:gfd_20_muni_tax_avg_at_gfd_20y_gov_yield_avg])
    if :tipsinflation ∈ fields
      ltecdf::StatsBase.ECDF = ecdf(Vector{Float64}(collect(skipmissing(df[1:(yearIndex-1),:tnote20at]))))
      ltrate::Float64 = (1+df[yearIndex,:tnote20] * (1-df[yearIndex,:munitax]))/
        (1+df[yearIndex,:tipsinflation]) - 1
      stecdf::StatsBase.ECDF = ecdf(Vector{Float64}(collect(skipmissing(df[1:(yearIndex-1),:tnote1at]))))
      strate::Float64 = (1+df[yearIndex,:tnote1] * (1-df[yearIndex,:munitax]))/
        (1+df[yearIndex,:tipsinflation]) - 1
      println("Quantile of $σsAboveYear breakeven 20yr: ", ltecdf(ltrate))
      println("Quantile of $σsAboveYear breakeven 1yr: ", stecdf(strate))
    end

    writeTables2File([TableSummarySD],
      HEADER_NAME, FOOTER_NAME, path=OUTPUT_PATH,
      outName = outName)

  return TableSummarySD
end

function writeDataAndSummarize(df::DataFrame)
    #interlude: Calculate the municipal bond less the after-tax credit spread
    df[:muni20MinusSpread] = df[:gfd_20y_muni_yield_avg] .-
      (df[:gfd_lt_corp_yield_avg] .- df[:gfd_20y_gov_yield_avg]) .*
      (1.0 .- df[:gfd_20_muni_tax_avg])

     println("...df filtered")
    #This vector of tuples is of the format oldFields, newField
    fieldsOut::Vector{Tuple{Symbol, Symbol}} = [(:year, :year),
      (:cpi_avg, :cpi), (:inflation, :inflation), (:inflation_3y_smoothed_avg, :inflation3y),
      (:gfd_20y_gov_yield_avg, :tnote20),
      (:gfd_1y_gov_yield_avg, :tnote1), (:gfd_10y_gov_yield_avg, :tnote10), (:gfd_20y_muni_yield_avg, :muni20),
      (:gfd_lt_corp_yield_avg, :corplt), (:us_max_rate, :topord), (:marg_rate, :avgord),
      (:mean_us_avg_rate, :avgordtotal),
      (:gfd_20_muni_tax_avg_0Credit, :munitax0credit),  (:gfd_20_muni_tax_avg, :munitax), (:us_lt_gain_marg_rate, :gainstax),
      (:gfd_1y_gov_yield_avg_r, :tnote1r), (:gfd_20_muni_tax_avg_at_gfd_1y_gov_yield_avg, :tnote1at),
      (:gfd_20y_gov_yield_avg_r, :tnote20r), (:gfd_20_muni_tax_avg_at_gfd_20y_gov_yield_avg, :tnote20at),
      (:us_interest_marg_rate, :interesttax), (:us_st_interest_marg_rate, :interesttaxstate),
      (:t_debt_depository_b, :tbankamt), (:t_debt_insurance_b, :tinsuranceamt),
      (:t_debt_pension_local_gov_b, :tpensionlocalamt), (:t_debt_pension_private, :tpensionprivateamt),
      (:t_debt_local_gov_b, :tlocalgovamt), (:t_debt_privately_owned_b, :ttotalprivateamt),
      #(:vwesx_10y_realized, :vanguardig), (:vwesx_10y_realized_real, :realvanguardig),
      #(:vfiix_10y_realized, :vanguardagency), (:vfiix_10y_realized_real, :realvanguardagency),
      (:tips_breakeven_inflation_avg,:tipsinflation), (:highTaxHoldings, :hightaxholdings),
      (:taxExemptHoldings, :taxexemptholdings), (:residualHoldings, :residualholdings),
      (:muni20MinusSpread, :muni20minusspread)]


    oldNames::Vector{Symbol} = [f[1] for f::Tuple{Symbol,Symbol} ∈ fieldsOut]
    newNames::Vector{Symbol} = [f[2] for f::Tuple{Symbol,Symbol} ∈ fieldsOut]

    #only get the fields we need
    df = df[ [oldNames;]]
    names!(df, newNames)
    println("...renaming complete")
    #oStream::IOStream = open("plot.csv", "w")
    println("...opened")
    #uCSV.write("plot.csv", df)
    CSV.write("plot.csv", df)
    println("...plot.csv written")

    truncatedDF::DataFrame = df[df[:year] .≥ 1950, :]

    writeSummaryTable(truncatedDF, newNames[2:end], (string).(newNames[2:end]),
      σsAbove = 2.0)
    println("...summary table written")

    writeSelectedDataTable(truncatedDF)
    println("...selected summary table written")

    fieldsForResultsTable::Vector{Symbol} = [:inflation3y, :interesttax, :munitax,
      :tnote1, :tnote1r, :tnote1at, :tnote20, :tnote20r, :tnote20at]
    writeSummaryTable(truncatedDF[:,[:year; fieldsForResultsTable]],
      fieldsForResultsTable, (string).(fieldsForResultsTable), outName="resultsTable.tex")
end


function makeReturnGraphs!(taxDFm::DataFrame, taxDFq::DataFrame, taxDFy::DataFrame;
  viaFunction::Function = pickMid, exPostYears::Int = 10, beginOffset=-2,
  endOffset=1)::NTuple{3,DataFrame}

  #########First make derived fields
  ###Can write the code for fields which are common accross the data frames in a compact manner
  for df::DataFrame ∈ [taxDFm, taxDFq, taxDFy]
    df[:gfd_ls_t_spread_close] = df[:gfd_10y_gov_yield_close] .- df[:gfd_st_gov_yield_close]
    df[:gfd_ls_t_spread_avg] = df[:gfd_10y_gov_yield_avg] .- df[:gfd_st_gov_yield_avg]
    df[:gfd_10_1_t_spread_close] = df[:gfd_10y_gov_yield_close] .- df[:gfd_1y_gov_yield_close]
    df[:gfd_10_1_t_spread_avg] = df[:gfd_10y_gov_yield_avg] .- df[:gfd_1y_gov_yield_avg]
    df[:gfd_10_2_t_spread_close] = df[:gfd_10y_gov_yield_close] .- df[:gfd_2y_gov_yield_close]
    df[:gfd_10_2_t_spread_avg] = df[:gfd_10y_gov_yield_avg] .- df[:gfd_2y_gov_yield_avg]
    df[:gfd_10_20_t_spread_close] = df[:gfd_10y_gov_yield_close] .- df[:gfd_20y_gov_yield_close]
    df[:gfd_10_20_t_spread_avg] = df[:gfd_10y_gov_yield_avg] .- df[:gfd_20y_gov_yield_avg]
    df[:gfd_20_10_2_t_spread_close] = df[:gfd_20y_gov_yield_close] .-
      df[:gfd_10y_gov_yield_close].- df[:gfd_2y_gov_yield_close]
    df[:gfd_20_10_2_t_spread_avg] = df[:gfd_20y_gov_yield_avg] .-
      df[:gfd_10y_gov_yield_avg] .- df[:gfd_2y_gov_yield_avg]
    df[:gfd_20_10_t_spread_close] = df[:gfd_20y_gov_yield_close] .- df[:gfd_10y_gov_yield_close]
    df[:gfd_20_10_t_spread_avg] = df[:gfd_20y_gov_yield_avg] .- df[:gfd_10y_gov_yield_avg]
    df[:gfd_20_2_t_spread_close] = df[:gfd_20y_gov_yield_close] .- df[:gfd_2y_gov_yield_close]
    df[:gfd_20_2_t_spread_avg] = df[:gfd_20y_gov_yield_avg] .- df[:gfd_2y_gov_yield_avg]
    df[:gfd_muni_spread_close] = df[:gfd_10y_gov_yield_close] .- df[:gfd_20y_muni_yield_close]
    df[:gfd_muni_spread_avg] = df[:gfd_10y_gov_yield_avg] .- df[:gfd_20y_muni_yield_avg]
    df[:gfd_20_muni_spread_close] = df[:gfd_20y_gov_yield_close] .- df[:gfd_20y_muni_yield_close]
    df[:gfd_20_muni_spread_avg] = df[:gfd_20y_gov_yield_avg] .- df[:gfd_20y_muni_yield_avg]
    df[:gfd_20_muni_cr_close] = df[:gfd_20y_muni_yield_close] .+
      (df[:gfd_lt_corp_yield_close] - df[:gfd_20y_gov_yield_close]) .* (1.0 .- df[:gfd_20_muni_tax_close])
    df[:gfd_20_muni_cr_avg] = df[:gfd_20y_muni_yield_avg] .+
      (df[:gfd_lt_corp_yield_avg] .- df[:gfd_20y_gov_yield_avg]) .* (1.0 .- df[:gfd_20_muni_tax_avg])
    df[:gfd_20y_muni_spread_cr_close] = df[:gfd_20y_gov_yield_close]  .- df[:gfd_20_muni_cr_close]
    df[:gfd_20y_muni_spread_cr_avg] = df[:gfd_20y_gov_yield_avg] .- df[:gfd_20_muni_cr_avg]
    df[:gfd_muni_spread_cr_close] = df[:gfd_10y_gov_yield_close] .-
      (df[:gfd_20y_muni_yield_close] .+ df[:gfd_credit_spread_close] .* (1.0 .- df[:gfd_muni_tax_close]))
    df[:gfd_muni_spread_cr_avg] = df[:gfd_10y_gov_yield_avg] .-
      (df[:gfd_20y_muni_yield_avg] .+ df[:gfd_credit_spread_avg] .* (1.0 .- df[:gfd_muni_tax_avg]))

  end

  ##create derived monthly DF series
  periodsPerYear::Int = 12
  taxDFm[:fred_t_st] = fuseDataSeries(taxDFm, :t_st_yield, :t_3mo_yield, viaFunction=viaFunction)
  taxDFm[:gfd_10_t_cpiyoy_spread_close] = taxDFm[:gfd_10y_gov_yield_close] .- taxDFm[:cpi_inflation]
  taxDFm[:gfd_10_t_cpiyoy_spread_avg] = taxDFm[:gfd_10y_gov_yield_avg] .- taxDFm[:cpi_inflation]
  taxDFm[ :inflation] = taxDFm[:cpi_inflation]
  taxDFm[ :inflation_10y_expost] = getExPostFromIndex(taxDFm[:cpi], exPostYears * periodsPerYear, periodsPerYear)
  taxDFm[ :inflation_10y_realized] = getRealizedFromIndex(taxDFm[:cpi], exPostYears * periodsPerYear, periodsPerYear)
  taxDFm[ :inflation_3y_smoothed] = get3YSmoothedFromIndex(taxDFm[:cpi],
    beginOffset * periodsPerYear, endOffset * periodsPerYear, periodsPerYear)
  taxDFm[:gfd_10_corp_spread_close] = taxDFm[:gfd_10y_corp_yield] .- taxDFm[:gfd_10y_gov_yield_close]
  taxDFm[:gfd_10_corp_spread_avg] = taxDFm[:gfd_10y_corp_yield] .- taxDFm[:gfd_10y_gov_yield_avg]

  ##create derived quarterly DF series
  periodsPerYear = 4
  taxDFq[:fred_t_st_close] = fuseDataSeries(taxDFq, :t_st_yield_close, :t_3mo_yield_close, viaFunction=viaFunction)
  taxDFq[:fred_t_st_avg] = fuseDataSeries(taxDFq, :t_st_yield_avg, :t_3mo_yield_avg, viaFunction=viaFunction)
  taxDFq[:gfd_10_t_cpiyoy_spread_close] = taxDFq[:gfd_10y_gov_yield_close].-taxDFq[:cpi_inflation_close]
  taxDFq[:gfd_10_t_cpiyoy_spread_avg] = taxDFq[:gfd_10y_gov_yield_avg].-taxDFq[:cpi_inflation_avg]
  taxDFq[ :inflation_close] = taxDFq[:cpi_inflation_close]
  taxDFq[ :inflation_avg] = taxDFq[:cpi_inflation_avg]
  taxDFq[ :inflation] = taxDFq[:cpi_inflation_avg] #!!!!
  taxDFq[ :inflation_10y_expost_close] = getExPostFromIndex(taxDFq[:cpi_close], exPostYears .* periodsPerYear, periodsPerYear)
  taxDFq[ :inflation_10y_expost_avg] = getExPostFromIndex(taxDFq[:cpi_avg], exPostYears .* periodsPerYear, periodsPerYear)
  taxDFq[ :inflation_10y_realized_close] = getRealizedFromIndex(taxDFq[:cpi_close], exPostYears .* periodsPerYear, periodsPerYear)
  taxDFq[ :inflation_10y_realized_avg] = getRealizedFromIndex(taxDFq[:cpi_avg], exPostYears .* periodsPerYear, periodsPerYear)
  taxDFq[ :inflation_3y_smoothed_close] = get3YSmoothedFromIndex(taxDFq[:cpi_close],
    beginOffset * periodsPerYear, endOffset * periodsPerYear, periodsPerYear)
  taxDFq[ :inflation_3y_smoothed_avg] = get3YSmoothedFromIndex(taxDFq[:cpi_avg],
    beginOffset * periodsPerYear, endOffset * periodsPerYear, periodsPerYear)
  taxDFq[:gfd_10_corp_spread_close] = taxDFq[:gfd_10y_corp_yield_close] .- taxDFq[:gfd_10y_gov_yield_close]
  taxDFq[:gfd_10_corp_spread_avg] = taxDFq[:gfd_10y_corp_yield_avg] .- taxDFq[:gfd_10y_gov_yield_avg]

  ##create derived yearly DF series
  periodsPerYear = 1
  taxDFy[:fred_t_st_close] = fuseDataSeries(taxDFy, :t_st_yield_close, :t_3mo_yield_close, viaFunction=viaFunction)
  taxDFy[:fred_t_st_avg] = fuseDataSeries(taxDFy, :t_st_yield_avg, :t_3mo_yield_avg, viaFunction=viaFunction)
  taxDFy[:gfd_10_t_cpiyoy_spread_close] = taxDFy[:gfd_10y_gov_yield_close].-taxDFy[:cpi_inflation_close]
  taxDFy[:gfd_10_t_cpiyoy_spread_avg] = taxDFy[:gfd_10y_gov_yield_avg].-taxDFy[:cpi_inflation_avg]
  taxDFy[ :inflation_close] = taxDFy[:cpi_inflation_close]
  taxDFy[ :inflation_avg] = taxDFy[:cpi_inflation_avg]
  taxDFy[ :inflation] = taxDFy[:cpi_inflation_avg]
  taxDFy[ :inflation_10y_expost_close] = getExPostFromIndex(taxDFy[:cpi_close], exPostYears)
  taxDFy[ :inflation_10y_expost_avg] = getExPostFromIndex(taxDFy[:cpi_avg], exPostYears)
  taxDFy[ :inflation_10y_realized_close] = getRealizedFromIndex(taxDFy[:cpi_close], exPostYears)
  taxDFy[ :inflation_10y_realized_avg] = getRealizedFromIndex(taxDFy[:cpi_avg], exPostYears)
  taxDFy[ :inflation_3y_smoothed_close] = get3YSmoothedFromIndex(taxDFy[:cpi_close],
    beginOffset * periodsPerYear, endOffset * periodsPerYear, periodsPerYear)
  taxDFy[ :inflation_3y_smoothed_avg] = get3YSmoothedFromIndex(taxDFy[:cpi_avg],
    beginOffset * periodsPerYear, endOffset * periodsPerYear, periodsPerYear)
  taxDFy[:gfd_10_corp_spread_close] = taxDFy[:gfd_10y_corp_yield_close] .- taxDFy[:gfd_10y_gov_yield_close]
  taxDFy[:gfd_10_corp_spread_avg] = taxDFy[:gfd_10y_corp_yield_avg] .- taxDFy[:gfd_10y_gov_yield_avg]

  #make the tax ratio fields
  taxDFy[:tax_ratio_ew_cumex] = 1.0 .- (1.0 .- taxDFy[:cumex_ew_tot_tax]) ./
    (1.0 .- taxDFy[:us_lt_gain_marg_rate])
  taxDFy[:tax_ratio_muni_avg] = 1.0 .- (1.0 .- taxDFy[:fred_muni_tax_avg]) ./
    (1.0 .- taxDFy[:us_lt_gain_marg_rate])

  #make the vanguard fields
  taxDFy[ :vwesx_10y_expost] = getExPostFromIndex(taxDFy[ :vwesx], exPostYears)
  taxDFy[ :vwesx_10y_realized] = getRealizedFromIndex(taxDFy[ :vwesx], exPostYears)
  taxDFy[ :vwesx_10y_expost_real] = (1.0 .+ taxDFy[ :vwesx_10y_expost]) ./
      (1.0 .+ getRealizedFromIndex(taxDFy[ :cpi_avg], exPostYears)) .- 1.0
  taxDFy[ :vwesx_10y_realized_real] = (1.0 .+ taxDFy[ :vwesx_10y_realized]) ./
      (1.0 .+ getRealizedFromIndex(taxDFy[ :cpi_avg], exPostYears)) .- 1.0

  taxDFy[ :vfiix_10y_expost] = getExPostFromIndex(taxDFy[ :vfiix], exPostYears)
  taxDFy[ :vfiix_10y_realized] = getRealizedFromIndex(taxDFy[ :vfiix], exPostYears)
  taxDFy[ :vfiix_10y_expost_real] = (1.0 .+ taxDFy[ :vfiix_10y_expost]) ./
      (1.0 .+ getRealizedFromIndex(taxDFy[ :cpi_avg], exPostYears)) .- 1.0
  taxDFy[ :vfiix_10y_realized_real] = (1.0 .+ taxDFy[ :vfiix_10y_realized]) ./
      (1.0 .+ getRealizedFromIndex(taxDFy[ :cpi_avg], exPostYears)) .- 1.0


  graphFrames::Vector{DataFrame} = Vector{Vector{Symbol}}()
  plotNames::Vector{String} =
  #comment out inflation rates other than 3 years
  [  "Real After Tax Long Term Treasury Return, 3-year Smoothed Inflation-1",
    "Real After Tax Long Term Treasury Return, 3-year Smoothed Inflation-2",

    "Real After Tax Short Term Treasury Return, 3-year Smoothed Inflation-1",
    "Real After Tax Short Term Treasury Return, 3-year Smoothed Inflation-2",

    "Real After Tax Corporate Bond Return, 3-year Smoothed Inflation-1",
    "Real After Tax Corporate Bond Return, 3-year Smoothed Inflation-2",

    "Real After Tax Corporate Bond Fund Return (VWESX), 3-year Smoothed Inflation-1",
    "Real After Tax Corporate Bond Fund Return (VWESX), 3-year Smoothed Inflation-2",

    "Real After Tax Agency Fund Return (VFIIX), 3-year Smoothed Inflation-1",
    "Real After Tax Agency Fund Return (VFIIX), 3-year Smoothed Inflation-2",

    "Nominal Interest Rate Series-1", "Nominal Interest Rate Series-2", "Nominal Interest Rate Series-3",
    "Real (ex-post) Interest Rate Series-1", "Real (ex-post) Interest Rate Series-2",
    "Implied and Statutory Tax Rates-1", "Implied and Statutory Tax Rates-2", "Implied and Statutory Tax Rates-3",

    "Muni and Cumex Implied Rates-1", "Muni and Cumex Implied Rates-2"
    ]

  #this is so we don't make a ton of extraneous fields
  plotDF = deepcopy(taxDFy)


  #set the marginal rate
  plotDF[ :marg_rate_1000] = plotDF[ :us_wage_marg_rate_1000]
  plotDF[ :marg_rate] = plotDF[:mean_us_marg_rate]
  plotDF[ :avg_rate] = plotDF[:mean_us_avg_rate]
  plotDF[ :marg_interest_rate] = plotDF[ :us_interest_marg_rate]
  outDF::DataFrame = deepcopy(plotDF)

  #set the fields we are interested in
  taxFields::Vector{Symbol} = [:fred_muni_tax_avg, :gfd_20_muni_tax_avg,
    :tax_ratio_muni_avg, :avg_rate, :marg_rate, :marg_rate_1000, :us_max_rate]
  returnFields::Vector{Symbol} = [:gfd_20y_gov_yield_avg, :gfd_1y_gov_yield_avg,
    :gfd_lt_corp_yield_avg, :vwesx_10y_expost, :vfiix_10y_expost]
  inflationFields::Vector{Symbol} = [#=:inflation_10y_expost_avg, :inflation_10y_realized_avg, =#:inflation_3y_smoothed_avg]
  #this symbol corresponds with the resulting after tax statistic
  statNames::Vector{Symbol} = ((s::Symbol)->Symbol(s,"_at")).(taxFields)

  Gadfly.push_theme(:default)
  plots::Vector{PlotContainer} = Vector{PlotContainer}()
  palette::Vector{Union{String, Gadfly.ColorTypes.RGB}} =
    ["SandyBrown", "Black", "DodgerBlue", colorant"#000001", colorant"#1F90FF"] #work around to have two of same color

   println("...Graph data processing complete")
  for iInd ∈ 1:length(inflationFields)
    for rInd ∈ 1:length(returnFields)
        #calculate (1+r*(1-t)/(1+π))
      for tInd ∈ 1:length(taxFields)
        plotDF[statNames[tInd]]  = ((1.0 .+ plotDF[returnFields[rInd]] .*
          (1.0 .- plotDF[taxFields[tInd]])) ./ (1.0 .+ plotDF[inflationFields[iInd]]) .- 1.0)
        newField::Symbol = Symbol(statNames[tInd], "_$(returnFields[rInd])")

        #record the data
        outDF[newField] = missings(Float64, size(outDF,1))
        outDF[newField] .= missing
        outDF[newField]  .= plotDF[statNames[tInd]]

      end

      #longDF::DataFrame = stackdf(plotDF[[:year; statNames]], statNames, [:year])
      #completecases!(longDF)

      minYear::Int = (floor(minimum(plotDF[completecases(plotDF[
        [inflationFields[iInd], returnFields[rInd]]]), :year])/10.0))*10
      maxYear::Int = (ceil(maximum(plotDF[completecases(plotDF[
        [inflationFields[iInd], returnFields[rInd]]]), :year])/10.0))*10

      push!(plots, plot(
          layer(plotDF[completecases(plotDF[[:year,statNames[1]]]),:], x="year", y=statNames[1], Geom.line,Theme(line_width=1.5pt, default_color=palette[1])),
          layer(plotDF[completecases(plotDF[[:year,statNames[2]]]),:], x="year", y=statNames[2], Geom.line,Theme(line_width=1.5pt, default_color=palette[2])),
          layer(plotDF[completecases(plotDF[[:year,statNames[3]]]),:], x="year", y=statNames[3], Geom.line,Theme(line_width=1.5pt, default_color=palette[3])),
          layer(yintercept=[0.0], Geom.hline(color=colorant"black", size=0.75pt)),
          Guide.ylabel("Return"), Guide.xlabel(nothing),
          Guide.xticks(ticks=collect(minYear:10:maxYear)),
          Guide.yticks(ticks=collect(-0.1:0.05:0.1)),
          Coord.cartesian(xmin=minYear, xmax=maxYear, ymin=-0.1, ymax=0.1),
          style(key_position = :right, key_max_columns = 6, key_label_font_size=8pt, key_label_font="courier"),
          Guide.title("Implied Real Rates of Return, Tax Rates from Financial Spreads"),
          Guide.manual_color_key("Legend",
            ["FRED Muni Tax Rt", "GFD Tax Rt", "Muni Tax/Cap Gain"],
            palette[1:3])))
      push!(plots, plot(
        layer(plotDF[completecases(plotDF[[:year,statNames[4]]]),:], x="year",
                y=statNames[4], Geom.line,Theme(line_width=1.5pt, default_color=palette[1])) ,
        layer(plotDF[completecases(plotDF[[:year,statNames[5]]]),:], x="year",
                y=statNames[5], Geom.line,Theme(line_width=1.5pt, default_color=palette[2])) ,
        layer(plotDF[completecases(plotDF[[:year,statNames[6]]]),:], x="year",
                y=statNames[6], Geom.line,Theme(line_width=2.5pt, default_color=palette[4]#=,
          line_style=:dash=#)), #the dashed line didn't work last I tried it
        layer(plotDF[completecases(plotDF[[:year,statNames[7]]]),:], x="year", y=statNames[7], Geom.line,Theme(line_width=1.5pt, default_color=palette[3])),
        layer(yintercept=[0.0], Geom.hline(color=colorant"black", size=0.75pt)),
        Guide.ylabel("Return"), Guide.xlabel(nothing),
        Guide.xticks(ticks=collect(minYear:10:maxYear)),
        Guide.yticks(ticks=collect(-0.1:0.05:0.1)),
        Coord.cartesian(xmin=minYear, xmax=maxYear, ymin=-0.1, ymax=0.1),
        style(key_position = :right, key_max_columns = 6, key_label_font_size=8pt, key_label_font="courier"),
        Guide.title("Implied Real Rates of Return, Based on NBER Statutory Tax Rates"),
        Guide.manual_color_key("Legend",
          ["Avg Tax","Marg Tax", "Marg >1MM (dashed)", "Max Tax"],
          palette[[1,2,4,3]])))
      push!(graphFrames, deepcopy(plotDF[[:year; statNames[1:3]]]))
      push!(graphFrames, deepcopy(plotDF[[:year; statNames[4:end]]]))
    end
  end

  println("...ror graphs complete")
  ###########################begin work on summary graphs

  ##Absolute Rates
  summaryRates::Vector{Symbol} = [:gfd_1y_gov_yield_avg, :gfd_10y_gov_yield_avg,  :gfd_20y_gov_yield_avg,
    :gfd_20_muni_cr_avg, :gfd_lt_corp_yield_avg, :gfd_10y_corp_yield_avg, :gfd_20y_muni_yield_avg,
    :gfd_credit_spread_avg, :gfd_10_corp_spread_avg, :gfd_20_muni_spread_avg, :gfd_20y_muni_spread_cr_avg]
  push!(plots, plot(
      layer(plotDF[completecases(plotDF[[:year,summaryRates[1]]]),:], x="year", y=summaryRates[1], Geom.line,Theme(line_width = 1.5pt, default_color=palette[1])),
      layer(plotDF[completecases(plotDF[[:year,summaryRates[2]]]),:], x="year", y=summaryRates[2], Geom.line,Theme(line_width = 1.5pt, default_color=palette[2])),
      layer(plotDF[completecases(plotDF[[:year,summaryRates[3]]]),:], x="year", y=summaryRates[3], Geom.line,Theme(line_width = 1.5pt, default_color=palette[3])),
      layer(plotDF[completecases(plotDF[[:year,summaryRates[4]]]),:], x="year", y=summaryRates[4], Geom.line,Theme(line_width = 2.5pt,
        default_color=palette[5]#=, line_style=:dot=#)),
      layer(yintercept=[0.0], Geom.hline(color=colorant"black", size=0.75pt)),
      Guide.ylabel("Rate"), Guide.xlabel(nothing),
      Guide.xticks(ticks=collect(1910:20:2020)),
      Guide.yticks(ticks=collect(0.0:0.05:0.15)),
      Guide.title("Treasury Interest Rates"),
      style(key_position = :right, key_max_columns = 1, key_label_font_size=8pt, key_label_font="courier"),
        Guide.manual_color_key("Legend", ["1Y T-Note", "10Y T-Note", "20Y T-Bond", "Muni+(1-t)Cr"],
        palette[[1,2,3,5]])))
    push!(plots, plot(
      layer(plotDF[completecases(plotDF[[:year,summaryRates[5]]]),:], x="year", y=summaryRates[5], Geom.line,Theme(line_width = 1.5pt, default_color=palette[1])),
      layer(plotDF[completecases(plotDF[[:year,summaryRates[6]]]),:], x="year", y=summaryRates[6], Geom.line,Theme(line_width = 1.5pt, default_color=palette[2])),
      layer(plotDF[completecases(plotDF[[:year,summaryRates[7]]]),:], x="year", y=summaryRates[7], Geom.line,Theme(line_width = 1.5pt, default_color=palette[3])),
      layer(yintercept=[0.0], Geom.hline(color=colorant"black", size=0.75pt)),
      Guide.ylabel("Rate"), Guide.xlabel(nothing),
      Guide.xticks(ticks=collect(1910:20:2020)),
      Guide.yticks(ticks=collect(0.0:0.05:0.15)),
      Guide.title("Other Nominal Interest Rates"),
      style(key_position = :right, key_max_columns = 1, key_label_font_size=8pt, key_label_font="courier"),
        Guide.manual_color_key("Legend", ["LT AAA Corp", "10Y A-AAA Corp", "20Y Muni           "],
        palette[1:3])))
    push!(plots,plot(
      layer(plotDF[completecases(plotDF[[:year,summaryRates[8]]]),:], x="year", y=summaryRates[8], Geom.line,Theme(line_width = 1.5pt, default_color=palette[1])),
      layer(plotDF[completecases(plotDF[[:year,summaryRates[9]]]),:], x="year", y=summaryRates[9], Geom.line,Theme(line_width = 1.5pt, default_color=palette[2])),
      layer(plotDF[completecases(plotDF[[:year,summaryRates[10]]]),:], x="year", y=summaryRates[10], Geom.line,Theme(line_width = 1.5pt, default_color=palette[3])),
      layer(plotDF[completecases(plotDF[[:year,summaryRates[11]]]),:], x="year", y=summaryRates[11], Geom.line,Theme(line_width = 2.5pt,
        default_color=palette[5]#=,   line_style=:dot=#)),
      layer(yintercept=[0.0], Geom.hline(color=colorant"black", size=0.75pt)),
      Guide.ylabel("Rate"), Guide.xlabel(nothing),
      Guide.xticks(ticks=collect(1910:20:2020)),
      Guide.yticks(ticks=collect(-0.02:0.02:0.02)),
      Guide.title("Interest Rate Spreads"),
      style(key_position = :right, key_max_columns = 1, key_label_font_size=8pt, key_label_font="courier"),
        Guide.manual_color_key("Legend", ["LT Corp-10Y T-Note", "10Y Corp-10Y T-Note",
          "20Y T-Bond-20Y Muni", "T-Note-(Muni+(1-t)Cr)"],
        palette[[1,2,3,5]])))

    push!(graphFrames, deepcopy(plotDF[[:year; summaryRates[1:4]]]))
    push!(graphFrames, deepcopy(plotDF[[:year; summaryRates[5:7]]]))
    push!(graphFrames, deepcopy(plotDF[[:year; summaryRates[8:end]]]))

    println("...Spreads and rates graphs complete (nominal)")
    ##Real Rates
    summaryRates = [:gfd_1y_gov_yield_avg, :gfd_10y_gov_yield_avg,  :gfd_20y_gov_yield_avg,
      :gfd_20_muni_cr_avg, :gfd_lt_corp_yield_avg, :gfd_10y_corp_yield_avg, :gfd_20y_muni_yield_avg]

    #calculate the real rates
    for i::Int ∈ 1:length(summaryRates)
      newField::Symbol = Symbol(summaryRates[i],"_r")
      plotDF[newField] =
        (1.0 .+ plotDF[summaryRates[i]]) ./ (1.0 .+ plotDF[:inflation_3y_smoothed_avg]) .- 1.0
      summaryRates[i] = newField

      outDF[newField] = missings(Float64, size(outDF,1))
      outDF[newField] .= missing
      outDF[newField] .= plotDF[newField]
    end

    push!(plots, plot(
        layer(plotDF[completecases(plotDF[[:year,summaryRates[1]]]),:], x="year", y=summaryRates[1], Geom.line,Theme(line_width = 1.5pt, default_color=palette[1])),
        layer(plotDF[completecases(plotDF[[:year,summaryRates[2]]]),:], x="year", y=summaryRates[2], Geom.line,Theme(line_width = 1.5pt, default_color=palette[2])),
        layer(plotDF[completecases(plotDF[[:year,summaryRates[3]]]),:], x="year", y=summaryRates[3], Geom.line,Theme(line_width = 1.5pt, default_color=palette[3])),
        layer(plotDF[completecases(plotDF[[:year,summaryRates[4]]]),:], x="year", y=summaryRates[4], Geom.line,Theme(line_width = 2.5pt,
          default_color=palette[5]#=, line_style=:dot=#)),
        layer(yintercept=[0.0], Geom.hline(color=colorant"black", size=0.75pt)),
        Guide.ylabel("Rate"), Guide.xlabel(nothing),
        Guide.xticks(ticks=collect(1910:20:2020)),
        Guide.yticks(ticks=collect(-0.05:0.05:0.1)),
        Guide.title("Treasury Interest Rates"),
        style(key_position = :right, key_max_columns = 1, key_label_font_size=8pt, key_label_font="courier"),
          Guide.manual_color_key("Legend", ["1Y T-Note", "10Y T-Note", "20Y T-Bond","Muni+(1-t)Cr"],
          palette[[1,2,3,5]])))
      push!(plots, plot(
        layer(plotDF[completecases(plotDF[[:year,summaryRates[5]]]),:], x="year", y=summaryRates[5], Geom.line,Theme(line_width = 1.5pt, default_color=palette[1])),
        layer(plotDF[completecases(plotDF[[:year,summaryRates[6]]]),:], x="year", y=summaryRates[6], Geom.line,Theme(line_width = 1.5pt, default_color=palette[2])),
        layer(plotDF[completecases(plotDF[[:year,summaryRates[7]]]),:], x="year", y=summaryRates[7], Geom.line,Theme(line_width = 1.5pt, default_color=palette[3])),
        layer(yintercept=[0.0], Geom.hline(color=colorant"black", size=0.75pt)),
        Guide.ylabel("Rate"), Guide.xlabel(nothing),
        Guide.xticks(ticks=collect(1910:20:2020)),
        Guide.yticks(ticks=collect(-0.05:0.05:0.1)),
        Guide.title("Other Real Interest Rates"),
        style(key_position = :right, key_max_columns = 1, key_label_font_size=8pt, key_label_font="courier"),
          Guide.manual_color_key("Legend", ["LT AAA Corp", "10Y A-AAA Corp", "20Y Muni"],
          palette[1:3])))

      push!(graphFrames, deepcopy(plotDF[[:year; summaryRates[1:4]]]))
      push!(graphFrames, deepcopy(plotDF[[:year; summaryRates[5:end]]]))
  println("...Spreads and rates graphs complete (real)")
  ## Tax Rates


  summaryTaxRates::Vector{Symbol} = [:fred_muni_tax_avg, :gfd_muni_tax_avg, :gfd_20_muni_tax_avg,
    :fred_muni_tax_avg_0Credit, :gfd_muni_tax_avg_0Credit, :gfd_20_muni_tax_avg_0Credit,
    :us_max_rate, :marg_rate, :marg_interest_rate, :us_lt_gain_marg_rate]
  push!(plots, plot(
    layer(plotDF[completecases(plotDF[[:year,summaryTaxRates[1]]]),:], x="year", y=summaryTaxRates[1], Geom.line,Theme(line_width = 1.5pt, default_color=palette[1])),
    layer(plotDF[completecases(plotDF[[:year,summaryTaxRates[2]]]),:], x="year", y=summaryTaxRates[2], Geom.line,Theme(line_width = 1.5pt, default_color=palette[2])),
    layer(plotDF[completecases(plotDF[[:year,summaryTaxRates[3]]]),:], x="year", y=summaryTaxRates[3], Geom.line,Theme(line_width = 1.5pt, default_color=palette[3])),
    layer(yintercept=[0.0], Geom.hline(color=colorant"black", size=0.75pt)),
    Guide.ylabel("Tax Rate"), Guide.xlabel(nothing),
    Guide.xticks(ticks=collect(1910:20:2020)),
    Guide.yticks(ticks=collect(0.0:0.2:1.0)),
    Coord.cartesian(xmin=1900, xmax=2020, ymin=0.0, ymax=1.0),
    Guide.title("Financial Securities Implied Tax Rates"),
    style(key_position = :right, key_max_columns = 7, key_label_font_size=8pt, key_label_font="courier"),
    Guide.manual_color_key("", ["Fred Muni", "GFD Muni", "20Y Muni"],
      palette[1:3])))
  push!(plots, plot(
    layer(plotDF[completecases(plotDF[[:year,summaryTaxRates[4]]]),:], x="year", y=summaryTaxRates[4], Geom.line,Theme(line_width = 1.5pt, default_color=palette[1])),
    layer(plotDF[completecases(plotDF[[:year,summaryTaxRates[5]]]),:], x="year", y=summaryTaxRates[5], Geom.line,Theme(line_width = 1.5pt, default_color=palette[2])),
    layer(plotDF[completecases(plotDF[[:year,summaryTaxRates[6]]]),:], x="year", y=summaryTaxRates[6], Geom.line,Theme(line_width = 1.5pt, default_color=palette[3])),
    layer(yintercept=[0.0], Geom.hline(color=colorant"black", size=0.75pt)),
    Guide.ylabel("Tax Rate"), Guide.xlabel(nothing),
    Guide.xticks(ticks=collect(1910:20:2020)),
    Guide.yticks(ticks=collect(0.0:0.2:1.0)),
    Coord.cartesian(xmin=1900, xmax=2020, ymin=0.0, ymax=1.0),
    Guide.title("Financial Securities Implied Tax Rates (No Credit Spread)"),
    style(key_position = :right, key_max_columns = 7, key_label_font_size=8pt, key_label_font="courier"),
    Guide.manual_color_key("", ["Fred Muni", "GFD Muni", "20Y Muni"],
      palette[1:3])))
  push!(plots, plot(
      layer(plotDF[completecases(plotDF[[:year,summaryTaxRates[7]]]),:], x="year", y=summaryTaxRates[7], Geom.line,Theme(line_width = 1.5pt, default_color=palette[1])),
      layer(plotDF[completecases(plotDF[[:year,summaryTaxRates[8]]]),:], x="year", y=summaryTaxRates[8], Geom.line,Theme(line_width = 1.5pt, default_color=palette[2])),
      layer(plotDF[completecases(plotDF[[:year,summaryTaxRates[9]]]),:], x="year", y=summaryTaxRates[9], Geom.line,Theme(line_width = 1.5pt,
        default_color=palette[2]#=, line_style=:dash=#)),
      layer(plotDF[completecases(plotDF[[:year,summaryTaxRates[10]]]),:], x="year", y=summaryTaxRates[10], Geom.line,Theme(line_width = 1.5pt, default_color=palette[3])),
      layer(yintercept=[0.0], Geom.hline(color=colorant"black", size=0.75pt)),
      Guide.ylabel("Tax Rate"), Guide.xlabel(nothing),
      Guide.title("NBER and Other Statutory Tax Rates"),
      Guide.xticks(ticks=collect(1910:20:2020)),
      Guide.yticks(ticks=collect(0.0:0.2:1.0)),
      Coord.cartesian(xmin=1900, xmax=2020, ymin=0.0, ymax=1.0),
      style(key_position = :right, key_max_columns = 7, key_label_font_size=8pt, key_label_font="courier"),
      Guide.manual_color_key("", ["Top Ordinary ", "Avg Marg (Ord.)",  "Avg Marg (Interest)", "Top LT Gains "],
        palette[[1,2,4,3]])))

    push!(graphFrames, deepcopy(plotDF[[:year; summaryRates[1:3]]]))
    push!(graphFrames, deepcopy(plotDF[[:year; summaryRates[4:6]]]))
    push!(graphFrames, deepcopy(plotDF[[:year; summaryRates[7:end]]]))
  println("...tax graphs complete")
  ###Statutory taxes for dividend analysis
  ## Tax Rates

  #make the ribbons for the ribbon plot
  plotDF[ :cumex_ew_tot_tax_h] =
    (plotDF[:div_levels] .- minimum(skipmissing(plotDF[:div_levels]))) ./
    (maximum(skipmissing(plotDF[:div_levels])) .- minimum(skipmissing(plotDF[:div_levels]))) .* .1
  plotDF[:cumex_ew_tot_tax_l] = plotDF[ :cumex_ew_tot_tax] .- plotDF[ :cumex_ew_tot_tax_h]
  plotDF[ :cumex_ew_tot_tax_h] += plotDF[ :cumex_ew_tot_tax]

  #note this is not the fully synced 20y muni
  summaryTaxRates = [:gfd_muni_tax_avg, :cumex_ew_tot_tax,
    :us_max_rate, :top_div_rate, :marg_rate, :us_lt_gain_marg_rate]
  push!(plots, plot(
    layer(plotDF[completecases(plotDF[[:year,summaryTaxRates[1]]]),:], x="year", y=summaryTaxRates[1], Geom.line,Theme(line_width = 1.5pt, default_color=palette[1])),
    layer(plotDF[completecases(plotDF[[:year,summaryTaxRates[2]]]),:], x="year", y=summaryTaxRates[2],
      ymin=:cumex_ew_tot_tax_l, ymax=:cumex_ew_tot_tax_h, Geom.line,Geom.ribbon,
      Theme(line_width = 1pt, default_color=palette[2], lowlight_color= (x)->color(colorant"Green"), lowlight_opacity = 0.5)),
    layer(yintercept=[0.0], Geom.hline(color=colorant"black", size=0.75pt)),
    Guide.ylabel("Tax Rate"), Guide.xlabel(nothing),
    Guide.xticks(ticks=collect(1910:20:2020)),
    Guide.yticks(ticks=collect(0.0:0.2:1.0)),
    Coord.cartesian(xmin=1900, xmax=2020, ymin=0.0, ymax=1.0),
    Guide.title("Financial Securities Implied Tax Rates"),
    style(key_position = :right, key_max_columns = 7, key_label_font_size=8pt, key_label_font="courier"),
    Guide.manual_color_key("", ["GFD Muni Rate", "Cumex Rate"],
      palette[1:2])))
  push!(plots, plot(
      layer(plotDF[completecases(plotDF[[:year,summaryTaxRates[3]]]),:], x="year", y=summaryTaxRates[3], Geom.line,Theme(line_width = 1.5pt,
        default_color=palette[2])),
      layer(plotDF[completecases(plotDF[[:year,summaryTaxRates[4]]]),:], x="year", y=summaryTaxRates[4], Geom.line,Theme(line_width = 2.5pt,
        default_color=palette[4]#=, line_style=:dash=#)),
      layer(plotDF[completecases(plotDF[[:year,summaryTaxRates[5]]]),:], x="year", y=summaryTaxRates[5], Geom.line,Theme(line_width = 1.5pt,
        default_color=palette[1])),
      layer(plotDF[completecases(plotDF[[:year,summaryTaxRates[6]]]),:], x="year", y=summaryTaxRates[6], Geom.line,Theme(line_width = 1.5pt,
        default_color=palette[3])),
      layer(yintercept=[0.0], Geom.hline(color=colorant"black", size=0.75pt)),
      Guide.ylabel("Tax Rate"), Guide.xlabel(nothing),
      Guide.title("NBER and Other Statutory Tax Rates"),
      Guide.xticks(ticks=collect(1910:20:2020)),
      Guide.yticks(ticks=collect(0.0:0.2:1.0)),
      Coord.cartesian(xmin=1900, xmax=2020, ymin=0.0, ymax=1.0),
      style(key_position = :right, key_max_columns = 7, key_label_font_size=8pt, key_label_font="courier"),
      Guide.manual_color_key("", ["Top Ordinary ", "Top Dividend (bold)", "Avg Marg (Ord.)",  "LT Gains"],
        palette[[2,4,1,3]])))
      push!(graphFrames, deepcopy(plotDF[[:year; summaryRates[1:2]]]))
      push!(graphFrames, deepcopy(plotDF[[:year; summaryRates[3:6]]]))
  println("...dividend graphs complete")

  for i ∈ 1:length(plots) #NOTE: Gadfly is too buggy to use as of 9/15/18, try again later
    #draw(PNG("$OUTPUT_PATH\\$(plotNames[i]).png", 8inch, 2inch),plots[i])
    CSV.write("$OUTPUT_PATH\\$(plotNames[i]).csv", graphFrames[i])
    println("CSV $(plotNames[i]) written.")
  end


  writeDataAndSummarize(outDF)

  return taxDFm, taxDFq, taxDFy
end

#this structure keeps track of the different fields used in the regressions
mutable struct TaxSpec
  indName::String
  depName::String

  taxField::CTSym
  muniField::CTSym
  muniTField::CTSym
  corpField::CTSym
  corpTField::CTSym

  muniSpread::CTSym
  corpSpread::CTSym
  TAT::CTSym
  CSAT::CTSym
end

function getFieldNames(t::TaxSpec)::Vector{Symbol}
  return [t.taxField, t.muniField, t.muniTField, t.corpField, t.corpTField,
    t.muniSpread, t.corpSpread, t.TAT, t.CSAT]
end

#keyword wrapper for tax specs
function TaxSpec(; indName::String = "", depName::String = "",
  taxField::CTSym = nothing,
  muniField::CTSym = DEFAULT_MUNI_FIELD, muniTField::CTSym = nothing,
  corpField::CTSym = nothing, corpTField::CTSym = nothing)

  #taxBonus::Symbol = Symbol(taxField, "_tb")
  muniSpread::Symbol = Symbol(taxField, "_ms")
  corpSpread::Symbol = Symbol(taxField, "_cs")
  TAT::Symbol = Symbol(taxField, "_tat")
  CSAT::Symbol = Symbol(taxField, "_csat")

  return TaxSpec(indName, depName, taxField, muniField, muniTField, corpField,
    corpTField, muniSpread, corpSpread, TAT, CSAT)
end


#this prepares a dataframe for regression. The dataframe is copied prior to modification
function prepRegression(df::DataFrame, taxSpecs::Vector{TaxSpec}, accountForStates::Bool = false)::DataFrame

  df = deepcopy(df) #we arn't interested in modifying the data here



  #rename some variables for convenience
  df[ :marg_rate_1000] = df[ :us_wage_marg_rate_1000]
  df[ :marg_rate] = df[:mean_us_marg_rate]
  df[ :avg_rate] = df[:mean_us_avg_rate]
  #df[ :marg_interest_rate] = df[ :us_interest_marg_rate]

  #we want to calculate the impact of the tax in units of yield
  for t::TaxSpec ∈ taxSpecs
    df[t.muniSpread] = df[ t.muniTField] .- df[ t.muniField]
    df[t.corpSpread] =
      (df[ t.corpTField] .- df[ t.corpField]) #.* (1.0 .+ df[t.taxField])

    df[t.TAT] = df[ t.muniTField] .* df[t.taxField]
    df[t.CSAT] = df[ t.corpField] .* df[t.taxField]


    #df[t.taxBonus] = (df[t.muniTField] .+ df[ t.corpSpread]) .* df[t.taxField]
    #println(skipmissing(df[t.corpSpread]))
  end
  #writetable("$OUTPUT_PATH\\df_test.csv",df, nastring="")

  for c ∈ 1:size(df,2)

    if eltype(df[c])==Any || eltype(df[c])==Union{Any,Missing}
      df[c] = Vector{Union{Float64,Missing}}(df[c])
    end
  end

  return df
end

#this function constructs a standard regression table
function standardReg(df::DataFrame, taxSpecs::Vector{TaxSpec},
  titleCaption::String, caption::String;
  otherSeriesStr::String = DEFAULT_OTHER_SERIES,
  otherSeriesSyms::Vector{Symbol} = DEFAULT_OTHER_SERIES_SYMS)
  otherSeriesNames::Vector{String} = DEFAULT_OTHER_SERIES_NAMES

  YSpecs::Vector{Symbol} = Vector{Symbol}()
  XSpecs::Vector{CTExpr} = Vector{CTExpr}()
  XNames::Vector{Vector{Symbol}} = Vector{Vector{Symbol}}()
  models::Vector{CTLM} = Vector{CTLM}()
  colNames::Vector{String} = Vector{String}()


  rowNames::Vector{String} =
    [ "20Y T-Bond Minus AAA Corp"; ["ibid*"*t.indName for t::TaxSpec ∈ taxSpecs];
    ["20Y T-Bond*"*t.indName for t::TaxSpec ∈ taxSpecs]; otherSeriesNames; "Intercept"]

  for t::TaxSpec ∈ taxSpecs
    push!(XSpecs, Meta.parse("$(t.TAT) + $(t.corpSpread) + $(t.CSAT)"))
    push!(XNames, [:intercept; t.taxField; :corpSpread; Symbol(t.taxField, "_CSAT")])
    push!(YSpecs, t.muniSpread)
    push!(colNames, t.depName)

    push!(XSpecs, Meta.parse("$(t.TAT) + $(t.corpSpread) + $(t.CSAT) + $(otherSeriesStr)"))
    push!(XNames, [:intercept; t.taxField; :corpSpread; Symbol(t.taxField, "_CSAT"); otherSeriesSyms])
    push!(YSpecs, t.muniSpread)
    push!(colNames, t.depName)
  end

  #run the regression
  for i::Int ∈ 1:length(XSpecs)
    println("Running spec $(XSpecs[i])")
    push!(models, CTLM(df, XSpecs[i], YSpecs[i], XNames=XNames[i], YName = :t_muni_spread))
  end

  descRowNames::Vector{String} = ["Num Years", "\$R^{2}\$"]

  descContent::Vector{Vector{String}} =
    [Vector{String}(undef, length(models)) for i∈ 1:length(descRowNames)]

  for i ∈ 1:length(XSpecs)
    descContent[1][i] = "$(models[i].N)"
    descContent[2][i] = "$(round(getR(models[i])^2.0,digits=2))"
  end

  #println(getFieldNames(taxSpecs[1]))
  taxSyms::Vector{Symbol} = [t.taxField for t::TaxSpec ∈ taxSpecs]
  CSATSyms::Vector{Symbol} = [Symbol(t.taxField, "_CSAT") for t::TaxSpec ∈ taxSpecs]

  TableText::String = texTable(models, getModWhiteΣ!, [:corpSpread; CSATSyms; taxSyms; otherSeriesSyms; :intercept],
    titleCaption = titleCaption,
    #colNames = [["($i)" for i∈1:length(models)]],
    colNames = [["Dependent Variable: Yield Spread (T-Bond--Muni)"],
        ["($i)" for i::Int ∈ 1:length(XSpecs)], colNames],
    contentRowNames = rowNames,
    descRowNames = descRowNames,
    descContent = descContent,
    stars=false,
    decimalDigits = 1,
    #columnSepPt = 0, #NOTE its possible we will want these back. Look at the old IO version to use it.
    starStrings = OVERRIDE_STAR_STRINGS,
    scaling = [100.0 for i::Int ∈ 1:length(rowNames)],
    notes=Vector{String}(),
    #clearMem = USE_AGGRESSIVE_GC,
    caption = caption,
    widthColNames = [[length(XNames)], ones(Int, length(XNames)), ones(Int, length(XNames))],
    alignmentColNames = [["c"]::Vector{String},
      ["r" for i::Int ∈ 1:length(XSpecs)]::Vector{String},
      ["r" for i::Int ∈ 1:length(XSpecs)]::Vector{String}])

    return TableText
end

#this function first-differences a data vector
function difference!(vec::Vector{Union{Float64, Missing}}, putForward::Bool = true)::Nothing

  if putForward
    vec[2:end] = vec[2:end] .- vec[1:(end-1)]
    vec[1] = missing
  else
    vec[1:(end-1)] = vec[2:end] .- vec[1:(end-1)]
    vec[end] = missing
  end

  return nothing
end

#this script function first first-differences the series before running the apporpriate regressions
#the data frame is copied before first differencing to avoid modifciation
function firstDifReg(dfΔ::DataFrame, taxSpecs::Vector{TaxSpec}, titleCaption::String, caption::String;
  otherSeriesStr::String = DEFAULT_OTHER_SERIES,
  otherSeriesSyms::Vector{Symbol} = DEFAULT_OTHER_SERIES_SYMS)

  dfΔ = deepcopy(dfΔ)

  #get the list of field names
  fieldNames::Vector{Symbol} = Vector{Symbol}()
  for taxSpec ∈ taxSpecs
    fieldNames = [fieldNames; getFieldNames(taxSpec); DEFAULT_OTHER_SERIES_SYMS]
  end



  #dedup the fields
  fieldNames = unique(fieldNames)

  #difference the appropriate fields

  for s::Symbol ∈ fieldNames

    if eltype(dfΔ[s]) ≠ Union{Float64, Missing}
      dfΔ[s]=Vector{Union{Float64, Missing}}(dfΔ[s])
    end

    difference!(dfΔ[s])
  end

  #run the regression as before
  return standardReg(dfΔ, taxSpecs, titleCaption, caption,
    otherSeriesStr = otherSeriesStr, otherSeriesSyms = otherSeriesSyms)
end


function getSDTable(df::DataFrame, fields::Vector{Symbol}, fieldNames::Vector{String};
  frequency::Float64=1.0, writeSDTable::Bool=false, label::Symbol=:full)::String
  N::Int = size(df,1)
  K::Int = length(fields)

  dfDA::DataFrame = deepcopy(df[fields])

  #de-annualize
  for s::Symbol ∈ fields
    dfDA[s] = (1.0 .+ dfDA[s]) .^ (1.0/frequency) .- 1.0
  end

  #calculate mean and standard deviation
  σsDA::Vector{Float64} = ((s::Symbol)->std(skipmissing(dfDA[s]))).(fields)
  μsDA::Vector{Float64} = ((s::Symbol)->mean(skipmissing(dfDA[s]))).(fields)
  cdfs::Vector{StatsBase.ECDF} = ((s::Symbol)->ecdf(collect(skipmissing(dfDA[s])))).(fields)

  σs::Vector{Float64} = σsDA .* frequency^0.5
  μs::Vector{Float64} = (1.0 .+ μsDA) .^ frequency .- 1.0

  σsDA .*= 100.0
  μsDA .*= 100.0
  σs .*= 100.0
  μs .*= 100.0

  colNames::Vector{Vector{String}} =
    [["$(fieldNames[i])" for i::Int ∈ 1:K],
    ["\$(\\mu=$(round(μs[i],digits=2)), \\sigma=$(round(σs[i],digits=2)))\$" for i::Int ∈ 1:K],
    vcat([["rate", "Z", "percentile"] for f::String ∈ fieldNames]...)]

  widthColNames::Vector{Vector{Int}} = [[3 for i::Int ∈ 1:K],[3 for i::Int ∈ 1:K],ones(Int,3*K)]

  #get the content
  rows::Vector{Vector{String}} = [Vector{String}(undef, 3*K) for r::Int ∈ 1:N]
  for r::Int ∈ 1:N
    for c::Int ∈ 1:K
      val::T where T<:Union{Missing,Float64} = df[r,fields[c]]
      valDA::T where T<:Union{Missing,Float64} = dfDA[r,fields[c]]
      rows[r][c*3-2] = ismissing(val) ? "" : "$(round(val*100.0,digits=2))"
      rows[r][c*3-1] = ismissing(valDA) ? "" : "$(round((valDA*100.0 - μsDA[c])/σsDA[c],digits=2))"
      rows[r][c*3] = ismissing(valDA) ? "" : "$(round(cdfs[c](valDA)*100.0,digits=2))"
    end
  end


  if frequency == 1
    rowNames::Vector{String} = string.(df[:year])
  else
    rowNames = ((d::Date)->"$(Dates.year(d))-$(Dates.monthname(d)[1:3])").(df[:date])
  end



  tableText::String = texTable("Table of Zs",
      "Standard deviations of returns", #caption
      colNames, #colNames
      Vector{String}(),#contentRowNames
      Vector{Matrix{String}}(), #content
      rowNames, #descRowNames
      rows, #descContent
      Vector{String}(), #notes
      widthColNames = widthColNames,
      #columnSepPt=25 #NOTE its possible we will want these back. Look at the old IO version to use it.
      )

  if writeSDTable
    colSyms::Vector{Symbol} = vcat([[Symbol("$(f)_rate"), Symbol("$(f)_Z"), Symbol("$(f)_perc")] for f::Symbol ∈ fields]...)
    for c::Int ∈ 1:length(colSyms)
      df[colSyms[c]] = ((r::Int)->rows[r][c]).(1:N)
    end
    CSV.write("$OUTPUT_PATH\\resultsSD_$label.csv",df[[:year; colSyms]])
  end




  return tableText

end


function getSDbyYear!(taxDFy::DataFrame)::Vector{String}
  taxFields::Vector{Symbol} = [:gfd_20_muni_tax_avg]
  returnFields::Vector{Symbol} = [:gfd_20y_gov_yield_avg, :gfd_1y_gov_yield_avg]
  inflationFields::Vector{Symbol} = [:inflation_3y_smoothed_avg]
  statNames::Vector{Symbol} = ((s::Symbol)->Symbol(s,"_at")).(taxFields)

  df = deepcopy(taxDFy[taxDFy[:year].≥1950,:])

  for iInd ∈ 1:length(inflationFields)
    for rInd ∈ 1:length(returnFields)
        #calculate (1+r*(1-t)/(1+π))
      for tInd ∈ 1:length(taxFields)
        df[statNames[tInd]]  = ((1.0 .+ df[returnFields[rInd]] .*
          (1.0 .- df[taxFields[tInd]])) ./ (1.0 .+ df[inflationFields[iInd]]) .- 1.0)
        newField::Symbol = Symbol(statNames[tInd], "_$(returnFields[rInd])")

        #record the data
        df[newField] = missings(Float64, size(df,1))
        df[newField] .= missing
        df[newField]  .= df[statNames[tInd]]
      end
    end
  end

  for i::Int ∈ 1:length(returnFields)
    newField::Symbol = Symbol(returnFields[i],"_r")
    df[newField] =
      (1.0 .+ df[returnFields[i]]) ./ (1.0 .+ df[:inflation_3y_smoothed_avg]) .- 1.0
    returnFields[i] = newField

  end
  #this symbol corresponds with the resulting after tax statistic

  fieldsST::Vector{Symbol} =
    [:gfd_1y_gov_yield_avg, :gfd_1y_gov_yield_avg_r, :gfd_20_muni_tax_avg_at_gfd_1y_gov_yield_avg]
  fieldNamesST::Vector{String} = ["1YR", "1YR-Real", "1YR-Real AT"]
  STTable::String = getSDTable(df, fieldsST, fieldNamesST, writeSDTable=true, label=:ST)

  fieldsLT::Vector{Symbol} =
    [:gfd_20y_gov_yield_avg,:gfd_20y_gov_yield_avg_r, :gfd_20_muni_tax_avg_at_gfd_20y_gov_yield_avg]
  fieldNamesLT::Vector{String} = ["20YR", "20YR-Real", "20YR-Real AT"]
  LTTable::String = getSDTable(df, fieldsLT, fieldNamesLT, writeSDTable=true, label=:LT)

  #=TableSummarySD = getSDSummaryTable(df, [fieldsST; fieldsLT;
    :vwesx_10y_realized; :vwesx_10y_realized_real; :gfd_lt_corp_yield_avg; :gfd_20y_muni_yield_avg],
    ["Nominal", "Real", "+Taxed","Nominal", "Real", "+Taxed",
    "10Y Act","Real", "20Y Corp", "20Y Muni"])=#

 # writetable("$DATA_PATH\\testSD.csv",df[[:year; fieldsST; fieldsLT; ]], nastring="")

  return [STTable, LTTable]
end


function getRegressionTables!(taxDFm::DataFrame, taxDFq::DataFrame, taxDFy::DataFrame)::Nothing

  #build the vector of tax specifications
  taxSpecs::Vector{TaxSpec} = Vector{TaxSpec}()


  push!(taxSpecs, TaxSpec(indName = "NBER Interest Tax", depName = "20y Minus Muni", taxField = :us_interest_marg_rate,
    muniTField = :gfd_20y_gov_yield_avg, corpField = :gfd_lt_corp_yield_avg, corpTField = :gfd_20y_gov_yield_avg))
  #push!(taxSpecs, TaxSpec(indName = "Max Stat. Tax", depName = "20y--Muni", taxField = :us_max_rate,
  #    muniTField = :gfd_20y_gov_yield_avg, corpField = :gfd_lt_corp_yield_avg, corpTField = :gfd_20y_gov_yield_avg))
  push!(taxSpecs, TaxSpec(indName = "NBER Gains Tax",depName = "20y Minus Muni", taxField = :us_lt_gain_marg_rate,
    muniTField = :gfd_20y_gov_yield_avg, corpField = :gfd_lt_corp_yield_avg, corpTField = :gfd_20y_gov_yield_avg))

  #flatten the appropriate variables
  for s ∈ [:top_div_rate, :us_wage_marg_rate_1000, :us_interest_marg_rate,
    :mean_us_avg_rate, :us_max_rate, :us_interest_marg_rate_1000]

    taxDFm[s] = flattenYears(taxDFy[:year],
      Vector{Union{T where T<:Real, Missing}}(taxDFy[s]), taxDFm[:date])
    taxDFq[s] = flattenYears(taxDFy[:year],
      Vector{Union{T where T<:Real, Missing}}(taxDFy[s]), taxDFq[:date])

  end

  df::DataFrame = prepRegression(taxDFy[taxDFy[:year] .≥ 1950,:], taxSpecs)

  titleCaption::String = "T-Bills -- Muni Bonds vs Assorted Tax Rates"
  caption = "See Compiled File."


  standardStr::String = standardReg(df, taxSpecs, titleCaption, caption)


  titleCaption = "T-Bills -- Muni Bonds vs Assorted Tax Rates (Differenced)"
  caption = "See Compiled File."
  difStr::String = firstDifReg(df, taxSpecs, titleCaption, caption)

  sdString::Vector{String} = getSDbyYear!(taxDFy)

  writeTables2File([standardStr; difStr; sdString],
    HEADER_NAME, FOOTER_NAME, path=OUTPUT_PATH,
    outName = "RegTables.tex")

  return nothing
end

#this function aggregates the treasury holding data
function meanHoldings!(taxDFq::DataFrame, taxDFy::DataFrame;
    fieldList::Vector{Symbol} =
    [:t_debt_tot_b, :t_debt_internal_b, :t_debt_privately_owned_b, :t_debt_depository_b,
    :t_debt_savings_bonds_b, :t_debt_pension_private, :t_debt_pension_local_gov_b,
    :t_debt_insurance_b, :t_debt_mutual_b, :t_debt_local_gov_b,
    :t_debt_intl_b, :t_debt_other_b])::DataFrame

  #so we can make temporary fields without clutter
  taxDFq = taxDFq[completecases(taxDFq[[:date; fieldList]]),[:date; fieldList]]
  taxDFq[:year] = (Dates.year).(taxDFq[:date])

  #index so we can access the dataframe by year
  yearInd::Dict = Dict(taxDFy[i,:year]=>i for i::Int ∈ 1:size(taxDFy,1))

  #make the new fields
  for f::Symbol ∈ fieldList
    taxDFy[f] = Vector{Union{Float64, Missing}}(undef, size(taxDFy,1))
    taxDFy[f] = missing
  end

  #split the higher resolution DF by year, and sum
  by(taxDFq, :year) do taxDFqSub::SubDataFrame
    yr::Int = taxDFqSub[1,:year]
    for f::Symbol ∈ fieldList
      if sum((!).((ismissing).(taxDFqSub[f]))) > 0
        taxDFy[yearInd[yr], f] = mean(skipmissing(taxDFqSub[f]))
      end
    end
  end

  #writetable("holdingsTest.csv", taxDFy[[:year; fieldList]])
  return taxDFy
end

#this categorizes the holdings
function processHoldings!(taxDFq::DataFrame, taxDFy::DataFrame)
  meanHoldings!(taxDFq, taxDFy)

  #fill the taxable holdings categories
  taxDFy[:highTaxHoldings] = taxDFy[:t_debt_depository_b] .+ taxDFy[:t_debt_insurance_b]
  taxDFy[:taxExemptHoldings] = taxDFy[:t_debt_pension_local_gov_b] .+ taxDFy[:t_debt_pension_private] .+
   taxDFy[:t_debt_local_gov_b]

  taxDFy[:residualHoldings] = taxDFy[ :t_debt_privately_owned_b] .-
    taxDFy[:highTaxHoldings] .- taxDFy[:taxExemptHoldings]

  #rewrite as percentage of holdings
  taxDFy[:highTaxHoldings] ./= taxDFy[ :t_debt_privately_owned_b]
  taxDFy[:taxExemptHoldings] ./= taxDFy[ :t_debt_privately_owned_b]
  taxDFy[:residualHoldings] ./= taxDFy[ :t_debt_privately_owned_b]

  return taxDFy
end


#this scripting function reads in the data and performs pre-processing activities
function processTaxData(;refreshWorkingData = true, replaceYahooData=true)::Nothing

  dtFormat::DateFormat = DateFormat(DATE_FORMAT_STR) ## a date format object

  if !isfile("$WORKING_PATH\\$DATA_FILE_NAME.jls") || refreshWorkingData
    writeTaxDataToDF()
  end

  iStream::IOStream = open("$WORKING_PATH\\$DATA_FILE_NAME.jls")
  taxDFm::DataFrame, taxDFq::DataFrame, taxDFy::DataFrame = deserialize(iStream)

  #adjust the types
  taxDFm[:date] = ((x::Int)->Date("$x",dtFormat)).(taxDFm[:date])
  taxDFq[:date] = ((x::Int)->Dates.firstdayofquarter(Date("$x",dtFormat))).(taxDFq[:date])

  sort!(taxDFm, [:date])
  sort!(taxDFq, [:date])
  sort!(taxDFy, [:year])

  println("Data loaded")
  if replaceYahooData
    processYahooData!(taxDFm, taxDFq, taxDFy)
  end
  #put these in a collection for convenience

  println("Yahoo processed.")
  scaleTaxData!(taxDFm, taxDFq, taxDFy)
  println("Tax data scaled")

  processMuniData!(taxDFm, taxDFq, taxDFy, viaFunction=pickMid)
  println("muniData processed")
  processCumExData!(taxDFm, taxDFq, taxDFy)
  println("Cumex processed")
  processHoldings!(taxDFq, taxDFy)
  println("Holdings processed")
  makeReturnGraphs!(taxDFm, taxDFq, taxDFy, viaFunction=pickMid)
  println("Return Graphs Made")
  getRegressionTables!(taxDFm, taxDFq, taxDFy)
  println("Regression tables made")

  CSV.write("$DATA_PATH\\$(DATA_FILE_NAME)-y.csv",taxDFy)
  CSV.write("$DATA_PATH\\$(DATA_FILE_NAME)-q.csv",taxDFq)
  CSV.write("$DATA_PATH\\$(DATA_FILE_NAME)-m.csv",taxDFm)

  return nothing
end

#top level global scripting
@time begin
  processTaxData(refreshWorkingData = true, replaceYahooData =  true)
end

end
