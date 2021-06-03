function minTest1()
  Rect3x2::Matrix{Float64} = [0.1 0.2; 0.3 0.4;0.5 0.6]
  Rect2x3::Matrix{Float64} = zeros(2,3)
  Rect2x3' .= Rect3x2
  println("Still Zeros: ", Rect2x3)
  Rect2x3 .= Rect3x2'
  println("Seems to Work: ", Rect2x3)
end

function foo(A, x, b=Vector{Float64}(length(x)))
  b .= A * x
  return b
end

println(foo([1 2; 3 4], [1,2]))
----------------
#this tests for a notice period negative bias to performance
#this version 
function NNP_1D(HFLong::DataFrame; dependentVar::Symbol = :performance1D,
    titleText::String =  "Performance Near Withdawels: Impact of Notice Period")::String
    thresholdLS::RedemptionThresholds = RedemptionThresholds(SIG_REDEMPTION_THRESHOLDS)

    modelsNNP_1D::Vector{CTLM} = Vector{CTLM}()

    XNNP_1DSpecs::Vector{CTExpr} = Vector{CTExpr}()
    XNNP_1DNames::Vector{Vector{Symbol}}  = Vector{Vector{Symbol}}()

    push!(XNNP_1DSpecs, parse("noticePeriod"))
    push!(XNNP_1DNames, [:intercept; :noticePeriod])

    push!(XNNP_1DSpecs, parse("noticePeriod + main_strategy"))
    push!(XNNP_1DNames, [:intercept; :noticePeriod])

    push!(XNNP_1DSpecs, parse("noticePeriod + fRedemptionDate + main_strategy"))
    push!(XNNP_1DNames, [:intercept; :noticePeriod])

    push!(XNNP_1DSpecs, parse("noticePeriod + fDate + main_strategy"))
    push!(XNNP_1DNames, [:intercept; :noticePeriod])

    push!(XNNP_1DSpecs, parse("noticePeriod&negLFlows"))
    push!(XNNP_1DNames, [:intercept; :noticePeriod_negLFlows])

    push!(XNNP_1DSpecs, parse("noticePeriod&negLFlows + main_strategy"))
    push!(XNNP_1DNames, [:intercept; :noticePeriod_negLFlows])

    push!(XNNP_1DSpecs, parse("noticePeriod&negLFlows + fRedemptionDate + main_strategy"))
    push!(XNNP_1DNames, [:intercept; :noticePeriod_negLFlows])

    push!(XNNP_1DSpecs, parse("noticePeriod&negLFlows + fDate + main_strategy"))
    push!(XNNP_1DNames, [:intercept; :noticePeriod_negLFlows])



    #push the specifications (use a common name to make a clean table)
    #=for s::Symbol ∈ thresholdLS.negLFlows

        push!(XNNP_1DSpecs, parse("noticePeriod&$s"))
        push!(XNNP_1DNames, [:intercept; :noticePeriod])

        push!(XNNP_1DSpecs, parse("noticePeriod&$s + fDate + main_strategy"))
        push!(XNNP_1DNames, [:intercept; :noticePeriod])

    end=#

    #run the model
    for iX ∈ 1:length(XNNP_1DSpecs)
        push!(modelsNNP_1D, CTLM(HFLong, XNNP_1DSpecs[iX],
            dependentVar, XNames=XNNP_1DNames[iX], YName = :performance))
    end

    ######IO Code for performance bias analysis

    descRowNamesNNP_1D::Vector{String} = ["Date F.E.",  "Redemp. Date F.E.","Strategy F.E.",#="Outflow \\% >", =# "N ('000s)"]

    #need to print the descriptive rows
    descContentNNP_1D::Vector{Vector{String}} =
        [Vector{String}(length(modelsNNP_1D)) for i∈ 1:length(descRowNamesNNP_1D)]

    #Determine the threshold levels (searches for the threshold level in each spec)
    thresholdVec::Vector{Int} = ((i::Int)->maximum(
        ((j::Int)->contains(XNNP_1DSpecs[i],"$(thresholdLS.negLFlows[j])")?
            Int(100.0*collect(SIG_REDEMPTION_THRESHOLDS)[j]):0).(1:length(thresholdLS.negLFlows))
        )).(1:length(XNNP_1DSpecs))

    for i ∈ 1:length(XNNP_1DSpecs)
        descContentNNP_1D[1][i] = "$(contains(XNNP_1DSpecs[i], "fDate")?"X":"")"
        descContentNNP_1D[2][i] = "$(contains(XNNP_1DSpecs[i], "fRedemptionDate")?"X":"")"
        descContentNNP_1D[3][i] = "$(contains(XNNP_1DSpecs[i], "main_strategy")?"X":"")"
        #descContentNNP_1D[3][i] = "$(contains(XNNP_1DSpecs[i], "negLFlows")?"X":"")"
        #descContentNNP_1D[4][i] = "$(thresholdVec[i])\\%"
        descContentNNP_1D[end][i] = "$(round(modelsNNP_1D[i].N/1000.0,1))"
    end

    #rows and variable names for the regression
    rowsNNP_1D::Vector{Symbol} = [:intercept; :noticePeriod; :noticePeriod_negLFlows]

    TableTextNNP_1D::String = texTable(modelsNNP_1D, getHomoskedΣ!, rowsNNP_1D,
        titleCaption = "Performance Flow-Bias Specifications",
        colNames = [["($i)" for i∈1:length(modelsNNP_1D)]],
        contentRowNames = ["intercept", "Notice Period", "Notice X Outflow"],
        descRowNames = descRowNamesNNP_1D,
        descContent = descContentNNP_1D,
        decimalDigits = 2,
        columnSepPt = -20,
        scaling = fill(100.0,length(rowsNNP_1D))::Vector{Float64},
        clearMem = USE_AGGRESSIVE_GC,
        caption = """NNP_1D specification compares performance prior to notice,
            with performance after notice. In some specifications, the notice
            period dummy is interacted with the amount of the outflow. Threshold
            refers to the minimum outflow included, thus some versions include
            only significant outflows. Data is from alive and dead hedge funds
            2000-2017. Interpret a value of 1.0 for the non-interacted
            specification as the discontinuity increasing performance by 1\\%.
            For the interacted specification, a 1.0 implies that outflows of
            10\\% corresponding to a 0.1\\% drop in monthly performance.""")

    return TableTextNNP_1D
end
