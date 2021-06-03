
function capacityio(;
  refreshio = PARAM[:refreshio])
  #provides an easy way to turn off writing the output
  if !PARAM[:refreshanyio]
    map(keys(refreshio) |> collect) do k
      refreshio[k] = false
    end
  end

  refreshio[:comomentumtable] && comomentumtables()
  refreshio[:fundtable] && fundtables()
  refreshio[:lltable] && lltables()
  #refreshio[:comomentumgraphs] && comomentumgraphs()
end
