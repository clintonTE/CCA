const category = Dict(:A => :Arts,
    :B => :Education,
    :C => :Environmental,
    :D => :AnimalWelfare,
    :E => :Health,
    :F => :MentalHealth,
    :G => :MedicalDiseases,
    :H => :MedicalResearch,
    :I => :CrimeAndLegal,
    :J => :Employment,
    :K => :FoodAndAgriculture,
    :L => :Housing,
    :M => :PublicSafety,
    :N => :Recreation,
    :O => :YouthDevelopment,
    :P => :HumanServices,
    :Q => :ForeignAffairs,
    :R => :SocialAction,
    :S => :CommunityImprovement,
    :T => :GrantmakingFoundations,
    :U => :ScienceandTechResearch,
    :V => :SocialScienceResearch,
    :W => :SocietyBenefitAndMultipurpose,
    :X => :ReligionRelated,
    :Y => :MutualBenefitOrg,
    :Z => :Unknown,
    :ZZ => :AllOthers
)

const charitytype = Dict(:pc => :privatecharity,
    :pf => :privatefoundation,
    :co => :othercharity
)

const affiliationcodes = Dict(
    Symbol(0) => missing,
    Symbol(1) => :CentralNoGroup,
    Symbol(2) => :Intermediate,
    Symbol(3) => :Independent,
    Symbol(6) => :CentralOfGroupNoChurch,
    Symbol(7) => :Intermediate,
    Symbol(8) => :CentralOfGroupChurch,
    Symbol(9) => :Subordinate,
    missing => missing,
    Symbol(0.0) => missing,
    :None => missing,
    Symbol(1.0) => :CentralNoGroup,
    Symbol(2.0) => :Intermediate,
    Symbol(3.0) => :Independent,
    Symbol(6.0) => :CentralOfGroupNoChurch,
    Symbol(7.0) => :Intermediate,
    Symbol(8.0) => :CentralOfGroupChurch,
    Symbol(9.0) => :Subordinate
)

const conformedaffiliationcodes = Dict(
    missing => :noncentral, #Unknown
    :Unknown => :noncentral, #Unknown
    :CentralNoGroup => :central, #CentralNoGroup
    :Intermediate => :noncentral, #Intermediate
    :Independent => :central, #Independent
    :CentralOfGroupNoChurch => :central, #CentralOfGroupNoChurch
    :Intermediate => :noncentral, #Intermediate
    :CentralOfGroupChurch => :central, #CentralOfGroupChurch
    :Subordinate => :noncentral #Subordinate
)



const naicsCodes = Dict(
    Symbol("11") => :Agriculture,
    Symbol("21") => :Mining,
    Symbol("22") => :Utilities,
    Symbol("23") => :Construction,
    Symbol("31") => :Manufacturing,
    Symbol("32") => :Manufacturing,
    Symbol("33") => :Manufacturing,
    Symbol("42") => :WholesaleTrade,
    Symbol("44") => :RetailTrade,
    Symbol("45") => :RetailTrade,
    Symbol("48") => :Transportation,
    Symbol("49") => :Transportation,
    Symbol("51") => :Information,
    Symbol("52") => :Finance,
    Symbol("53") => :RealEstate,
    Symbol("54") => :ProfServices,
    Symbol("55") => :Management,
    Symbol("56") => :AdminAndWaste,
    Symbol("61") => :Education,
    Symbol("62") => :HealthCare,
    Symbol("71") => :Entertainment,
    Symbol("72") => :Hospitality,
    Symbol("81") => :OtherServices,
    Symbol("92") => :PublicAdmin,
    Symbol("99") => :Unknown,
    missing => :Unknown,
    :No => :Unknown)

    #=NAICS unshortened (2017, from US Census)
    11	Agriculture, Forestry, Fishing and Hunting
    21	Mining, Quarrying, and Oil and Gas Extraction
    22	Utilities
    23	Construction
    31-33	Manufacturing
    42	Wholesale Trade
    44-45	Retail Trade
    48-49	Transportation and Warehousing
    51	Information
    52	Finance and Insurance
    53	Real Estate and Rental and Leasing
    54	Professional, Scientific, and Technical Services
    55	Management of Companies and Enterprises
    56	Administrative and Support and Waste Management and Remediation Services
    61	Educational Services
    62	Health Care and Social Assistance
    71	Arts, Entertainment, and Recreation
    72	Accommodation and Food Services
    81	Other Services (except Public Administration)
    92	Public Administration
    missing => :other,
    :No => :other
    =#

const ntee2naics2names = Dict(
    :AnimalWelfare => :OtherServices,
    :Arts => :Entertainment,
    :CommunityImprovement => :OtherServices,
    :CrimeAndLegal => :ProfServices,
    :Education => :Education,
    :Employment => :OtherServices,
    :Environmental => :OtherServices,
    :FoodAndAgriculture => :Agriculture,
    :ForeignAffairs => :OtherServices,
    :GrantmakingFoundations => :OtherServices,
    :Health => :HealthCare,
    :Housing => :OtherServices,
    :HumanServices => :HealthCare,
    :MedicalDiseases => :OtherServices,
    :MedicalResearch => :ProfServices,
    :MentalHealth => :HealthCare,
    :MutualBenefitOrg => :OtherServices,
    :PublicSafety => :HealthCare,
    :Recreation => :Entertainment,
    :ReligionRelated => :OtherServices,
    :ScienceandTechResearch => :ProfServices,
    :SocialAction => :OtherServices,
    :SocialScienceResearch => :ProfServices,
    :SocietyBenefitAndMultipurpose => :Utilities,
    :YouthDevelopment => :OtherServices,
    :Unknown=>:Unknown,
    missing=>:Unknown,
    :No => :Unknown,
    :AllOthers=>:Unknown
)

const LBM_XFIELDS = Dict(
  :bmsp500=>[:lsp500t30],
  :bmff3=>[:lmkt3t30, :lSMB3, :lHML3],
  :bmff5=>[:lmkt5t30, :lSMB5, :lHML5, :lRMW5, :lCMA5])

const FF_SUFFIX_FIELDS = Dict(
"3"=>[:Mkt_RF, :SMB, :HML, :RF],
"5"=>[:Mkt_RF, :SMB, :HML, :RMW, :CMA, :RF])

#misc dictionary for shortening field names
const SHORT_NAME_INDEX = Dict(
  :qlreturn3yr=>:ret1y,
  :qlreturn3yr=>:ret2y,
  :qlreturn3yr=>:ret3y,
  :qlreturn4yr=>:ret4y,
  :qlreturn5yr=>:ret5y,
  :qpprogramexpense=>:Exp,
  :qpcontributions=>:Contrib,
  :qDpcontributions=>:DifContrib,
  :pprogramexpense=>:Exp,
  :pprogramexpense_avg=>:Exp,
  :b1lRetCon => :b1lRetCon,
  :qb1lRetCon => :qb1lRetCon,
  :b1lRetLCon=>:b1lRetLCon,
  :qb1lRetLCon=>:qb1lRetLCon,
  :b1lRet3yrCon => :b1lRet3yrCon,
  :qb1lRetCon => :qb1lRetCon,
  :qb1lRet3yrCon => :qb1lRet3yrCon,
  :b1lRet3yrNetSP500Con=>:b1lRet3yrNetSP500Con,
  :qb1lRet3yrNetSP500Con=>:qb1lRet3yrNetSP500Con,
  :qpprogramexpense_avg=>:Exp,

)
