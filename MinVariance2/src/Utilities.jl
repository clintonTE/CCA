#use this to clear a dictionary
function clearcollection(c::AbstractDict)
  for i ∈ 1:length(c)
    pop!(c)
  end
end
