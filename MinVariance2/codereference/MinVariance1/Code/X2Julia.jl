const W1SOL1 = ("(S2.*(S3P.*ws1+(-1).*S1P.*ws3).^2+S2P.^2.*(S3.*ws1.^2+ws3.*((-2).* ...
  S13.*ws1+S1.*ws3))+(-2).*S2P.*(S1P.*S3.*ws1.*ws2+(-1).*S13.*S3P.* ...
  ws1.*ws2+(-1).*S12.*S3P.*ws1.*ws3+(-1).*S13.*S1P.*ws2.*ws3+S1.* ...
  S3P.*ws2.*ws3+S12.*S1P.*ws3.^2+S23.*ws1.*(S3P.*ws1+(-1).*S1P.*ws3) ...
  )+ws2.*(S3P.^2.*((-2).*S12.*ws1+S1.*ws2)+2.*S1P.*S3P.*(S23.*ws1+( ...
  -1).*S13.*ws2+S12.*ws3)+S1P.^2.*(S3.*ws2+(-2).*S23.*ws3))).^(-1).* ...
  (S2P.^2.*S3.*ws1+(-2).*S23.*S2P.*S3P.*ws1+S2.*S3P.^2.*ws1+(-1).* ...
  S1P.*S2P.*S3.*ws2+S1P.*S23.*S3P.*ws2+(-1).*S2P.*S3.*SGP.*ws1.*ws2+ ...
  S23.*S3P.*SGP.*ws1.*ws2+S1P.*S3.*SGP.*ws2.^2+S1P.*S23.*S2P.*ws3+( ...
  -1).*S1P.*S2.*S3P.*ws3+S23.*S2P.*SGP.*ws1.*ws3+(-1).*S2.*S3P.* ...
  SGP.*ws1.*ws3+(-2).*S1P.*S23.*SGP.*ws2.*ws3+S1P.*S2.*SGP.*ws3.^2+ ...
  S13.*S2P.*(S3P.*ws2+(-1).*S2P.*ws3)+(-1).*S12.*S3P.*(S3P.*ws2+(-1) ...
  .*S2P.*ws3)+(-1).*S13.*SGP.*ws2.*(S3P.*ws2+(-1).*S2P.*ws3)+S12.* ...
  SGP.*ws3.*(S3P.*ws2+(-1).*S2P.*ws3)+(-1).*(S3P.*ws2+(-1).*S2P.* ...
  ws3).^2.*((S3P.*ws2+(-1).*S2P.*ws3).^(-2).*((-1).*S1.*S2P.^2.*S3+ ...
  2.*S1.*S23.*S2P.*S3P+S12.^2.*S3P.^2+(-1).*S1.*S2.*S3P.^2+(-2).* ...
  S12.*S2P.*S3.*SGP.*ws1+2.*S12.*S23.*S3P.*SGP.*ws1+S2P.^2.*S3.*SG.* ...
  ws1.^2+(-2).*S23.*S2P.*S3P.*SG.*ws1.^2+S2.*S3P.^2.*SG.*ws1.^2+ ...
  S23.^2.*SGP.^2.*ws1.^2+(-1).*S2.*S3.*SGP.^2.*ws1.^2+2.*S1.*S2P.* ...
  S3.*SGP.*ws2+(-2).*S1.*S23.*S3P.*SGP.*ws2+(-2).*S12.*S3P.^2.*SG.* ...
  ws1.*ws2+2.*S12.*S3.*SGP.^2.*ws1.*ws2+S1.*S3P.^2.*SG.*ws2.^2+(-1) ...
  .*S1.*S3.*SGP.^2.*ws2.^2+S13.^2.*(S2P+(-1).*SGP.*ws2).^2+(-2).* ...
  S1.*S23.*S2P.*SGP.*ws3+(-2).*S12.^2.*S3P.*SGP.*ws3+2.*S1.*S2.* ...
  S3P.*SGP.*ws3+2.*S12.*S2P.*S3P.*SG.*ws1.*ws3+(-2).*S12.*S23.* ...
  SGP.^2.*ws1.*ws3+(-2).*S1.*S2P.*S3P.*SG.*ws2.*ws3+2.*S1.*S23.* ...
  SGP.^2.*ws2.*ws3+S1.*S2P.^2.*SG.*ws3.^2+S12.^2.*SGP.^2.*ws3.^2+( ...
  -1).*S1.*S2.*SGP.^2.*ws3.^2+S1P.^2.*(S23.^2+(-1).*S2.*S3+S3.*SG.* ...
  ws2.^2+(-2).*S23.*SG.*ws2.*ws3+S2.*SG.*ws3.^2)+(-2).*S1P.*(S13.*(( ...
  -1).*S2.*S3P+S3P.*SG.*ws2.^2+S23.*(S2P+(-1).*SGP.*ws2)+S2.*SGP.* ...
  ws3+(-1).*S2P.*SG.*ws2.*ws3)+S12.*((-1).*S2P.*S3+S23.*S3P+S3.* ...
  SGP.*ws2+(-1).*S23.*SGP.*ws3+(-1).*S3P.*SG.*ws2.*ws3+S2P.*SG.* ...
  ws3.^2)+ws1.*(S23.^2.*SGP+(-1).*S2.*S3.*SGP+S2P.*S3.*SG.*ws2+S2.* ...
  S3P.*SG.*ws3+(-1).*S23.*SG.*(S3P.*ws2+S2P.*ws3)))+(-2).*S13.*( ...
  S12.*(S2P+(-1).*SGP.*ws2).*(S3P+(-1).*SGP.*ws3)+ws1.*(S23.*SGP.*(( ...
  -1).*S2P+SGP.*ws2)+S2P.*SG.*((-1).*S3P.*ws2+S2P.*ws3)+S2.*SGP.*( ...
  S3P+(-1).*SGP.*ws3))))).^(1/2))")


const W1SOL2 = "(S2.*(S3P.*ws1+(-1).*S1P.*ws3).^2+S2P.^2.*(S3.*ws1.^2+ws3.*((-2).*
 S13.*ws1+S1.*ws3))+(-2).*S2P.*(S1P.*S3.*ws1.*ws2+(-1).*S13.*S3P.*
 ws1.*ws2+(-1).*S12.*S3P.*ws1.*ws3+(-1).*S13.*S1P.*ws2.*ws3+S1.*
 S3P.*ws2.*ws3+S12.*S1P.*ws3.^2+S23.*ws1.*(S3P.*ws1+(-1).*S1P.*ws3)
 )+ws2.*(S3P.^2.*((-2).*S12.*ws1+S1.*ws2)+2.*S1P.*S3P.*(S23.*ws1+(
 -1).*S13.*ws2+S12.*ws3)+S1P.^2.*(S3.*ws2+(-2).*S23.*ws3))).^(-1).*
 (S2P.^2.*S3.*ws1+(-2).*S23.*S2P.*S3P.*ws1+S2.*S3P.^2.*ws1+(-1).*
 S1P.*S2P.*S3.*ws2+S1P.*S23.*S3P.*ws2+(-1).*S2P.*S3.*SGP.*ws1.*ws2+
 S23.*S3P.*SGP.*ws1.*ws2+S1P.*S3.*SGP.*ws2.^2+S1P.*S23.*S2P.*ws3+(
 -1).*S1P.*S2.*S3P.*ws3+S23.*S2P.*SGP.*ws1.*ws3+(-1).*S2.*S3P.* SGP.*ws1.*ws3+(-2).*S1P.*S23.*SGP.*ws2.*ws3+S1P.*S2.*SGP.*ws3.^2+
 S13.*S2P.*(S3P.*ws2+(-1).*S2P.*ws3)+(-1).*S12.*S3P.*(S3P.*ws2+(-1)
 .*S2P.*ws3)+(-1).*S13.*SGP.*ws2.*(S3P.*ws2+(-1).*S2P.*ws3)+S12.*
 SGP.*ws3.*(S3P.*ws2+(-1).*S2P.*ws3)+(S3P.*ws2+(-1).*S2P.*ws3).^2.*
 ((S3P.*ws2+(-1).*S2P.*ws3).^(-2).*((-1).*S1.*S2P.^2.*S3+2.*S1.*
 S23.*S2P.*S3P+S12.^2.*S3P.^2+(-1).*S1.*S2.*S3P.^2+(-2).*S12.*S2P.*
 S3.*SGP.*ws1+2.*S12.*S23.*S3P.*SGP.*ws1+S2P.^2.*S3.*SG.*ws1.^2+(
 -2).*S23.*S2P.*S3P.*SG.*ws1.^2+S2.*S3P.^2.*SG.*ws1.^2+S23.^2.* SGP.^2.*ws1.^2+(-1).*S2.*S3.*SGP.^2.*ws1.^2+2.*S1.*S2P.*S3.*SGP.*
 ws2+(-2).*S1.*S23.*S3P.*SGP.*ws2+(-2).*S12.*S3P.^2.*SG.*ws1.*ws2+
 2.*S12.*S3.*SGP.^2.*ws1.*ws2+S1.*S3P.^2.*SG.*ws2.^2+(-1).*S1.*S3.*
 SGP.^2.*ws2.^2+S13.^2.*(S2P+(-1).*SGP.*ws2).^2+(-2).*S1.*S23.* S2P.*SGP.*ws3+(-2).*S12.^2.*S3P.*SGP.*ws3+2.*S1.*S2.*S3P.*SGP.*
 ws3+2.*S12.*S2P.*S3P.*SG.*ws1.*ws3+(-2).*S12.*S23.*SGP.^2.*ws1.*
 ws3+(-2).*S1.*S2P.*S3P.*SG.*ws2.*ws3+2.*S1.*S23.*SGP.^2.*ws2.*ws3+
 S1.*S2P.^2.*SG.*ws3.^2+S12.^2.*SGP.^2.*ws3.^2+(-1).*S1.*S2.* SGP.^2.*ws3.^2+S1P.^2.*(S23.^2+(-1).*S2.*S3+S3.*SG.*ws2.^2+(-2).*
 S23.*SG.*ws2.*ws3+S2.*SG.*ws3.^2)+(-2).*S1P.*(S13.*((-1).*S2.*S3P+
 S3P.*SG.*ws2.^2+S23.*(S2P+(-1).*SGP.*ws2)+S2.*SGP.*ws3+(-1).*S2P.*
 SG.*ws2.*ws3)+S12.*((-1).*S2P.*S3+S23.*S3P+S3.*SGP.*ws2+(-1).* S23.*SGP.*ws3+(-1).*S3P.*SG.*ws2.*ws3+S2P.*SG.*ws3.^2)+ws1.*( S23.^2.*SGP+(-1).*S2.*S3.*SGP+S2P.*S3.*SG.*ws2+S2.*S3P.*SG.*ws3+(
 -1).*S23.*SG.*(S3P.*ws2+S2P.*ws3)))+(-2).*S13.*(S12.*(S2P+(-1).*
 SGP.*ws2).*(S3P+(-1).*SGP.*ws3)+ws1.*(S23.*SGP.*((-1).*S2P+SGP.*
 ws2)+S2P.*SG.*((-1).*S3P.*ws2+S2P.*ws3)+S2.*SGP.*(S3P+(-1).*SGP.*
 ws3))))).^(1/2))"

# To move from Mathematica, call ToMatlab on the equation after installing the appropriate package
#also need to find-replace out the backslashes
#the only purpose of this function is to handle the above monstrosity
function matlab2julia(input::String = W1SOL2, minLineBreak = 80)
  notationOLD::Vector{Pair} =   [
    "..." => "",
    ".*" => " * ",
    ".^" => "^",
    "+" => " + ",
    "S1P" => "SGP[1]",
    "S2P" => "SGP[2]",
    "S3P" => "SGP[3]",
    "S12" => "TG.S[1,2]",
    "S23" => "TG.S[2,3]",
    "S13" => "TG.S[1,3]",
    "S1" => "TG.S[1,1]",
    "S2" => "TG.S[2,2]",
    "S3" => "TG.S[3,3]",
    "ws1" => "TG.wₛ[1]",
    "ws2" => "TG.wₛ[2]",
    "ws3" => "TG.wₛ[3]",
    "SG" => "Θ.SG",
    "SGP" => "Θ.SGP",
    "Θ.Θ" => "Θ",
    "Θ.SGP[" => "SGP[",
    ";" => "",
    "\r\n" => "",
    "\n" => "",
    "  " => " ",
    "  " => " ",
    "  " => " "
    ]

    notation::Vector{Pair} =   [
        "..." => "",
        ".*" => " * ",
        ".^" => "^",
        "+" => " + ",
        "S1P" => "CSGP[1]",
        "S2P" => "CSGP[2]",
        "S3P" => "CSGP[3]",
        "S12" => "CS[1,2]",
        "S23" => "CS[2,3]",
        "S13" => "CS[1,3]",
        "S1" => "CS[1,1]",
        "S2" => "CS[2,2]",
        "S3" => "CS[3,3]",
        "ws1" => "Cwₛ[1]",
        "ws2" => "Cwₛ[2]",
        "ws3" => "Cwₛ[3]",
        "SG" => "Θ.SG",
        "SGP" => "Θ.SGP",
        "Θ.Θ" => "Θ",
        "Θ.SGP[" => "SGP[",
        ";" => "",
        "\r\n" => "",
        "\n" => "",
        "  " => " ",
        "  " => " ",
        "  " => " "
        ]

    out::String = input
    for i ∈ 1:length(notation)
      out = replace(out, notation[i])
    end

    outarray::Vector{Char} = Vector{Char}(out)

    #insert line breaks for readability (although it won't help much)
    widthctr::Int = 0
    for i::Int ∈ 1:length(outarray)
      widthctr += 1
      if widthctr > minLineBreak && outarray[i] == ' '
        outarray[i] = '\n'
        widthctr = 0
      end
    end

    out = String(outarray)



  return "\n---------------------------\n(\n$out\n)"
end
#################################################
#################################################

const W1SOL1RAW = "(S2P^2*S3*ws1 - 2*S23*S2P*S3P*ws1 + S2*S3P^2*ws1 - S1P*S2P*S3*ws2 +
  S1P*S23*S3P*ws2 - S2P*S3*SGP*ws1*ws2 + S23*S3P*SGP*ws1*ws2 + S1P*S3*SGP*ws2^2 +
  S1P*S23*S2P*ws3 - S1P*S2*S3P*ws3 + S23*S2P*SGP*ws1*ws3 - S2*S3P*SGP*ws1*ws3 -
  2*S1P*S23*SGP*ws2*ws3 + S1P*S2*SGP*ws3^2 + S13*S2P*(S3P*ws2 - S2P*ws3) -
  S12*S3P*(S3P*ws2 - S2P*ws3) - S13*SGP*ws2*(S3P*ws2 - S2P*ws3) + S12*SGP*ws3*(S3P*ws2 - S2P*ws3) -
  (S3P*ws2 - S2P*ws3)^2* Sqrt[(1/(S3P*ws2 - S2P*ws3)^2)*((-S1)*S2P^2*S3 + 2*S1*S23*S2P*S3P +
  S12^2*S3P^2 - S1*S2*S3P^2 - 2*S12*S2P*S3*SGP*ws1 + 2*S12*S23*S3P*SGP*ws1 + S2P^2*S3*SG*ws1^2 -
  2*S23*S2P*S3P*SG*ws1^2 + S2*S3P^2*SG*ws1^2 + S23^2*SGP^2*ws1^2 - S2*S3*SGP^2*ws1^2 +
  2*S1*S2P*S3*SGP*ws2 - 2*S1*S23*S3P*SGP*ws2 - 2*S12*S3P^2*SG*ws1*ws2 + 2*S12*S3*SGP^2*ws1*ws2 +
  S1*S3P^2*SG*ws2^2 - S1*S3*SGP^2*ws2^2 + S13^2*(S2P - SGP*ws2)^2 - 2*S1*S23*S2P*SGP*ws3 -
  2*S12^2*S3P*SGP*ws3 + 2*S1*S2*S3P*SGP*ws3 + 2*S12*S2P*S3P*SG*ws1*ws3 - 2*S12*S23*SGP^2*ws1*ws3 -
  2*S1*S2P*S3P*SG*ws2*ws3 + 2*S1*S23*SGP^2*ws2*ws3 + S1*S2P^2*SG*ws3^2 + S12^2*SGP^2*ws3^2 -
  S1*S2*SGP^2*ws3^2 + S1P^2*(S23^2 - S2*S3 + S3*SG*ws2^2 - 2*S23*SG*ws2*ws3 + S2*SG*ws3^2) -
  2*S1P*(S13*((-S2)*S3P + S3P*SG*ws2^2 + S23*(S2P - SGP*ws2) + S2*SGP*ws3 - S2P*SG*ws2*ws3) +
  S12*((-S2P)*S3 + S23*S3P + S3*SGP*ws2 - S23*SGP*ws3 - S3P*SG*ws2*ws3 + S2P*SG*ws3^2) +
  ws1*(S23^2*SGP - S2*S3*SGP + S2P*S3*SG*ws2 + S2*S3P*SG*ws3 - S23*SG*(S3P*ws2 + S2P*ws3))) -
  2*S13*(S12*(S2P - SGP*ws2)*(S3P - SGP*ws3) + ws1*(S23*SGP*(-S2P + SGP*ws2) +
  S2P*SG*((-S3P)*ws2 + S2P*ws3) + S2*SGP*(S3P - SGP*ws3))))])/ (S2*(S3P*ws1 - S1P*ws3)^2 +
  S2P^2*(S3*ws1^2 + ws3*(-2*S13*ws1 + S1*ws3)) - 2*S2P*(S1P*S3*ws1*ws2 - S13*S3P*ws1*ws2 -
  S12*S3P*ws1*ws3 - S13*S1P*ws2*ws3 + S1*S3P*ws2*ws3 + S12*S1P*ws3^2 + S23*ws1*(S3P*ws1 - S1P*ws3)) +
  ws2*(S3P^2*(-2*S12*ws1 + S1*ws2) + 2*S1P*S3P*(S23*ws1 - S13*ws2 + S12*ws3) +
  S1P^2*(S3*ws2 - 2*S23*ws3)))"

const W1SOL2RAW = "(S2P^2*S3*ws1 - 2*S23*S2P*S3P*ws1 + S2*S3P^2*ws1 - S1P*S2P*S3*ws2 +
  S1P*S23*S3P*ws2 - S2P*S3*SGP*ws1*ws2 + S23*S3P*SGP*ws1*ws2 +  S1P*S3*SGP*ws2^2 +
  S1P*S23*S2P*ws3 - S1P*S2*S3P*ws3 +  S23*S2P*SGP*ws1*ws3 - S2*S3P*SGP*ws1*ws3 -
  2*S1P*S23*SGP*ws2*ws3 + S1P*S2*SGP*ws3^2 +  S13*S2P*(S3P*ws2 - S2P*ws3) - S12*S3P*(S3P*ws2 -
  S2P*ws3) -  S13*SGP*ws2*(S3P*ws2 - S2P*ws3) +  S12*SGP*ws3*(S3P*ws2 - S2P*ws3) +
  (S3P*ws2 - S2P*ws3)^2*  Sqrt[(1/(S3P*ws2 - S2P*ws3)^2)*((-S1)*S2P^2*S3 +  2*S1*S23*S2P*S3P +
  S12^2*S3P^2 - S1*S2*S3P^2 -  2*S12*S2P*S3*SGP*ws1 +  2*S12*S23*S3P*SGP*ws1 + S2P^2*S3*SG*ws1^2 -
  2*S23*S2P*S3P*SG*ws1^2 + S2*S3P^2*SG*ws1^2 +  S23^2*SGP^2*ws1^2 - S2*S3*SGP^2*ws1^2 +
  2*S1*S2P*S3*SGP*ws2 -  2*S1*S23*S3P*SGP*ws2 -  2*S12*S3P^2*SG*ws1*ws2 + 2*S12*S3*SGP^2*ws1*ws2 +
  S1*S3P^2*SG*ws2^2 - S1*S3*SGP^2*ws2^2 +  S13^2*(S2P - SGP*ws2)^2 - 2*S1*S23*S2P*SGP*ws3 -
  2*S12^2*S3P*SGP*ws3 + 2*S1*S2*S3P*SGP*ws3 +  2*S12*S2P*S3P*SG*ws1*ws3 - 2*S12*S23*SGP^2*ws1*ws3 -
  2*S1*S2P*S3P*SG*ws2*ws3 + 2*S1*S23*SGP^2*ws2*ws3 +  S1*S2P^2*SG*ws3^2 + S12^2*SGP^2*ws3^2 -
  S1*S2*SGP^2*ws3^2 +   S1P^2*(S23^2 - S2*S3 + S3*SG*ws2^2 - 2*S23*SG*ws2*ws3 +  S2*SG*ws3^2) -
  2*S1P*(S13*((-S2)*S3P + S3P*SG*ws2^2 + S23*(S2P - SGP*ws2) +  S2*SGP*ws3 - S2P*SG*ws2*ws3) +
  S12*((-S2P)*S3 + S23*S3P + S3*SGP*ws2 - S23*SGP*ws3 -  S3P*SG*ws2*ws3 + S2P*SG*ws3^2) +
  ws1*(S23^2*SGP - S2*S3*SGP + S2P*S3*SG*ws2 +  S2*S3P*SG*ws3 - S23*SG*(S3P*ws2 + S2P*ws3))) -
  2*S13*(S12*(S2P - SGP*ws2)*(S3P - SGP*ws3) +  ws1*(S23*SGP*(-S2P + SGP*ws2) +
  S2P*SG*((-S3P)*ws2 + S2P*ws3) +  S2*SGP*(S3P - SGP*ws3))))])/  (S2*(S3P*ws1 - S1P*ws3)^2 +
  S2P^2*(S3*ws1^2 + ws3*(-2*S13*ws1 + S1*ws3)) -  2*S2P*(S1P*S3*ws1*ws2 - S13*S3P*ws1*ws2 -
  S12*S3P*ws1*ws3 -  S13*S1P*ws2*ws3 + S1*S3P*ws2*ws3 + S12*S1P*ws3^2 +
  S23*ws1*(S3P*ws1 - S1P*ws3)) +  ws2*(S3P^2*(-2*S12*ws1 + S1*ws2) +  2*S1P*S3P*(S23*ws1 -
  S13*ws2 + S12*ws3) +  S1P^2*(S3*ws2 - 2*S23*ws3)))"

function mathematica2juliaW(input::String = W1SOL2RAW, minLineBreak = 80)
  notation::Vector{Pair} =   [
    "S1P" => "CSGP[1]",
    "S2P" => "CSGP[2]",
    "S3P" => "CSGP[3]",
    "S12" => "CS[1,2]",
    "S23" => "CS[2,3]",
    "S13" => "CS[1,3]",
    "S1" => "CS[1,1]",
    "S2" => "CS[2,2]",
    "S3" => "CS[3,3]",
    "ws1" => "Cwₛ[1]",
    "ws2" => "Cwₛ[2]",
    "ws3" => "Cwₛ[3]",
    "Sqrt[" => "sqrt(",
    ")])" => ")))",
    "SG" => "Θ.SG",
    "SGP" => "Θ.SGP",
    "Θ.Θ" => "Θ",
    "Θ.SGP[" => "SGP[",
    "\r\n" => "",
    "\n" => "",
    "  " => " ",
    "  " => " ",
    "  " => " "
      ]

    out::String = input
    for i ∈ 1:length(notation)
      out = replace(out, notation[i])
    end

    outarray::Vector{Char} = Vector{Char}(out)

    #insert line breaks for readability (although it won't help much)
    widthctr::Int = 0
    for i::Int ∈ 1:length(outarray)
      widthctr += 1
      if widthctr > minLineBreak && outarray[i] == ' '
        outarray[i] = '\n'
        widthctr = 0
      end
    end

    out = String(outarray)



  return "\n---------------------------\n(\n$out\n)"
end

const eqn1 = "μG=((-SGP)*(-z11 + z12) - SG*(z12 - z22))/(z11 - 2*z12 + z22)\n
γG2 = (-z12^2 + z11*z22)/(z11 - 2*z12 + z22)"

function mathematica2juliaEq(input::String = eqn1, minLineBreak = 80)
  notation::Vector{Pair} =   [
    "z11" => "Θ.ζ²G",
    "z22" => "Θ.ζ²P",
    "z12" => "Θ.ζGP",
    "sg" => "Θ.σ²G",
    "SG" => "Θ.SG",
    "SGP" => "Θ.SGP",
    "Θ.Θ" => "Θ",
    "Θ.SGP[" => "SGP[",
    "  " => " ",
    "  " => " ",
    "  " => " "
      ]

    out::String = input
    for i ∈ 1:length(notation)
      out = replace(out, notation[i])
    end

    outarray::Vector{Char} = Vector{Char}(out)

    #insert line breaks for readability (although it won't help much)
    widthctr::Int = 0
    for i::Int ∈ 1:length(outarray)
      widthctr += 1
      if widthctr > minLineBreak && outarray[i] == ' '
        outarray[i] = '\n'
        widthctr = 0
      end
    end

    out = String(outarray)



  return "\n---------------------------\n(\n$out\n)"
end
