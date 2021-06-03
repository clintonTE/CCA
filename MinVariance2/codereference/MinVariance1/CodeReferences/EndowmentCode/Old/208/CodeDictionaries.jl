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
