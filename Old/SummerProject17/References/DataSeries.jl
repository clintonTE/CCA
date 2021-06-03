module DataSeries

using  DataFrames, Distributions, StatsBase, Base.Dates

export coarseToFine, preProcessYahoo, flattenYears, flattenQuarters, getYahooData,
  dailyToEOM, dailyToEOQ, dailyToEOY

const YAHOO_DATE_FORMAT = "m/d/y"

const DEBUG_DATASERIES = true



function preProcessYahoo(path::String, name::String; makeReturns::Bool = true)::Void
  dateFormat::String = YAHOO_DATE_FORMAT
  DF::DataFrame = readtable("$path\\$name.csv")

  #fix the dates
  rename!(DF, :Date, :DATE)

  if makeReturns
    DF[:,:DATE] = ((s::String)->Date(s, dateFormat)::Date).(DF[:,:DATE])::Vector{Union{Date, Missing}}
    sort!(DF, cols = [:DATE])
    DF[:,:returns] = similar(DF[:,:Adj_Close])
    DF[2:end,:returns] = log.(DF[2:end,:Adj_Close] ./  DF[1:(end-1),:Adj_Close])
  end
  #showcols(DF)

  stream::IOStream = open("$path\\$name.jls", "w")
  serialize(stream, DF)
  close(stream)
  return nothing
end

#this function processes yahoo data, with an optional shortcut of the file read
function getYahooData(path::String, name::String; refreshData::Bool = true,
  makeReturns::Bool = true)::DataFrame

  if !isfile("$path\\$name.jls") || refreshData
    preProcessYahoo(path, name; makeReturns = makeReturns)
  end

  stream::IOStream = open("$path\\$name.jls")
    DF::DataFrame = deserialize(stream)
  close(stream)

  return DF
end

#this function gets the data point for the last available day in a month
function convertToMaxMonth(coarseDate::Date, fineDates::Vector{Union{Date, Missing}},
  fineData::Vector{Union{T, Missing}}, convertedFineDates::Vector{Union{Date, Missing}})::T where T<:Real

  maxDate::Date = maximum(fineDates[convertedFineDates.==coarseDate])
  return fineData[findfirst(fineDates, maxDate)]
end

dailyToEOM(coarseDates::Vector{Union{Date, Missing}},
  fineDates::Vector{Union{Date, Missing}}, fineData::Vector{Union{T, Missing}}) where T<: Real= fineToCoarse(coarseDates,
    fineDates, fineData, firstdayofmonth.(fineDates), convertToMaxMonth)

dailyToEOQ(coarseDates::Vector{Union{Date, Missing}},
  fineDates::Vector{Union{Date, Missing}}, fineData::Vector{Union{T, Missing}}) where T<: Real= fineToCoarse(coarseDates,
  fineDates, fineData, firstdayofquarter.(fineDates), convertToMaxMonth)

dailyToEOY(years::Vector{Union{Int, Missing}},
  fineDates::Vector{Union{Date, Missing}}, fineData::Vector{Union{T, Missing}}) where T<: Real= fineToCoarse(
  Vector{Union{Date, Missing}}(firstdayofyear.(Date.(years))),
  fineDates, fineData, firstdayofyear.(fineDates), convertToMaxMonth)

#fineToCoarse(coarseDF::DataFrame, fineDF::DataFrame, coarseField = :DATE, fine)

#this function takes in fine data and aggregates it upward
function fineToCoarse(coarseDates::Vector{Union{Date, Missing}},
  fineDates::Vector{Union{Date, Missing}}, fineData::Vector{Union{T, Missing}},
  convertedFineDates::Vector{Union{Date, Missing}},
  conversionFunc::Function)::Vector{Union{T, Missing}} where T<:Real

  #pre-allcoate the output vector
  out::Vector{Union{T, Missing}} = Vector{Union{T, Missing}}(Vector{T}(length(coarseDates)))

  for i::Int ∈ 1:length(coarseDates)
    if coarseDates[i] ∈ convertedFineDates
#      println("ctf test: ", coarseDates[i])
      out[i] = conversionFunc(coarseDates[i], fineDates, fineData, convertedFineDates)
    else
      out[i] = NA
    end
  end

  return out
end

#This function takes lower resolution time series and paste values over a higher resolution period
#IN: Data vectors for coarse dates, coarse data, higher resolution dates. Also requires a many to one
# mapping function of fine dates to coarse dates
#OUT: A data vector of the coarse data mapped to the higher reoslution fine dates
function coarseToFine(coarseDates::Vector{Union{Date, Missing}}, coarseData::Vector{Union{T, Missing}},
  fineDates::Vector{V}, conversionFunc::Function
  )::Vector where {T<:Real, V<:Union{Date,Missing}}



  #pre-allcoate the output vector
  out::Vector{Union{T, Missing}} = Vector{Union{T, Missing}}(Vector{Union{T, Missing}}(length(fineDates)))

  #create a lookup table for extacting coarse data given a coarse date
  dateToDataTable::Dict = Dict(coarseDates[i]=>coarseData[i] for i::Int ∈ 1:length(coarseDates))

  #get the data from the look-up table using the mapping function and the fine date input
  for i::Int ∈ 1:length(fineDates)
    out[i] = dateToDataTable[conversionFunc(fineDates[i])]
  end

  return out
end

flattenQuarters(coarseDates::Vector{Union{Date, Missing}}, coarseData::Vector{Union{T, Missing}},
  fineDates::Vector{Union{Date, Missing}}) where T<:Real = coarseToFine(coarseDates,
  coarseData, fineDates, Dates.firstdayofquarter)::Vector{Union{T, Missing}}

#wrapper function for coarseToFine which flattens yearly data to finer resolution
flattenYears(coarseDates::Vector{Union{Date, Missing}}, coarseData::Vector{Union{T, Missing}},
  fineDates::Vector{Date}) where T<:Real = coarseToFine(coarseDates,
  coarseData, fineDates, Dates.firstdayofyear)::Vector{Union{T, Missing}}

#helper method which takes years as an integer
flattenYears(years::Vector{Union{Int, Missing}}, coarseData::Vector{Union{T, Missing}},
  fineDates::Vector{Date}) where T<:Real = flattenYears(Vector{Union{Date, Missing}}(Date.(years)),
  coarseData, fineDates)::Vector{Union{T, Missing}}




end
