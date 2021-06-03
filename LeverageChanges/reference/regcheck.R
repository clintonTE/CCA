require(sandwich)

d = read.csv("plot.csv")
d = subset(d, year > 1946)

d["munispread"] = d["tnote20"] - d["muni20"]
d["creditspread"] = d["tnote20"] - d["corplt"]
d["taxtnote20"] = d["tnote20"] * d["interesttax"]
d["tnote10M20"] = d["tnote10"] - d["tnote20"]

s1 = as.formula("munispread ~ creditspread + taxtnote20")
s2 = as.formula("munispread ~ creditspread + taxtnote20 + tnote10M20 + inflation")

runreg = function(s, d) {
    reg = lm(s, d)
    print(reg)
    diag(NeweyWest(reg, lag = 3, prewhite = FALSE, adjust = TRUE)) ^ 0.5
}

runreg(s1,d)
runreg(s2,d)