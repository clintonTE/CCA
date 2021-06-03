#require iaw

setwd("C:\\Users\\Clinton Tepper\\Dropbox\\Projects\\histtaxes-R2\\histtaxes-R2")
d <- read.csv("plot.csv")
d <- subset(d, year > 1946)


scale <- 0.75

#print("Plot1")
pdf(file = "DONOTUSEInterestTaxesvWageTaxes.pdf", width = 12*scale, height = 6*scale)
with(d, {
    plot(year, topord, xaxs = "i", yaxs = "i", xlab = "Year", ylab = "Tax Rate", type = "n",
        cex.axis = 0.8, ylim = c(0, 0.5), xlim = c(1955, 2011))

    lines(year, avgord, col = "black", lwd = 3)
    points(year, avgord, col = "black", cex = 0.5)

    lines(year, interesttax, col = "blue", lwd = 3)
    points(year, interesttax, col = "black", cex = 0.5)

    text(1965, 0.30, "Ordinary interest tax rates\n from filed tax returns (NBER)", col = "blue", cex = 0.8, srt = 0)
    text(1965, 0.19, "Average tax rates from\n filed tax returns (NBER)", col = "black", cex = 0.8, srt = 0)
})
dev.off()


#print("Plot2")
#updated credit rates graph
pdf(file = "creditrates2.pdf", width = 12*scale, height = 6*scale)
with(d, {
    plot(year, xaxs = "i", yaxs = "i", xlab = "Year", ylab = "Nominal Interest Rate", type = "n", cex.axis = 0.8, ylim = c(0, 0.15), xlim = c(1950, 2018))
    lines(year, corplt, col = "blue", lwd = 3)
    lines(year, tnote20)
    lines(year, muni20, col = "blue", lwd = 3)

    points(year, corplt, col = "blue", cex = 0.5)
    points(year, tnote20, cex = 0.5, pty = 5)
    points(year, muni20, col = "blue", cex = 0.5)

    text(1968, 0.085, "LT AAA Corporates", col = "blue", cex = 0.9, srt = 30)
    text(1974.7, 0.0675, "LT Treasuries", col = "black", cex = 0.6, srt = 30)
    text(1978, 0.045, "LT AAA Municipals", col = "blue", cex = 0.9, srt = 30)

})
dev.off()

#print("Plot3")
pdf(file = "InterestTaxes.pdf", width = 12*scale, height = 6*scale)
with(d, {
    plot(year, interesttax, xaxs = "i", yaxs = "i", xlab = "Year", ylab = "Tax Rate", type = "n",
        cex.axis = 0.8, ylim = c(0, 0.5), xlim = c(1955, 2016))

    lines(year, avgord, col = "black", lwd = 3)
    points(year, avgord, col = "black", cex = 0.5)

    lines(year, interesttax, col = "blue", lwd = 3)
    points(year, interesttax, col = "blue", cex = 0.5)

    text(1975, 0.35, "Average marginal *ordinary income* tax rates from\n filed tax returns (NBER)", col = "black", cex =1 , srt = 0)
    text(1975, 0.17, "Average marginal *interest* tax rates\n from filed tax returns (NBER)", col = "blue", cex = 1, srt = 0)
})
dev.off()

#print("Plot4")
pdf(file = "FedWithStateTaxes.pdf", width = 12*scale, height = 6*scale)
with(d, {
    plot(year, topord, xaxs = "i", yaxs = "i", xlab = "Year", ylab = "Tax Rate", type = "n",
        cex.axis = 0.8, ylim = c(0, 0.5), xlim = c(1955, 2015))

    lines(year, interesttax, col = "black", lwd = 3)
    points(year, interesttax, col = "black", cex = 0.5)

    lines(year, interesttaxstate, col = "blue", lwd = 3)
    points(year, interesttaxstate, col = "blue", cex = 0.5)

    text(1990, 0.35, "Interest tax rates w/ state data\n from filed tax returns (NBER)", col = "blue", cex = 1, srt = 0)
    text(1965, 0.18, "Interest tax rates w/o state data\n from filed tax returns (NBER)", col = "black", cex = 1, srt = 0)
})
dev.off()

#print("Plot5")
pdf(file = "MarginalAndAverage.pdf", width = 12*scale, height = 6*scale)
with(d, {
    plot(year, avgordtotal, xaxs = "i", yaxs = "i", xlab = "Year", ylab = "Tax Rate", type = "n",
        cex.axis = 0.8, ylim = c(0, 0.5), xlim = c(1955, 2016))

    lines(year, avgordtotal, lwd = 3)
    points(year, avgordtotal, col = "black", cex = 0.5)

    lines(year, avgord, col = "blue", lwd = 3)
    points(year, avgord, col = "black", cex = 0.5)

    text(1980, 0.35, "Average *marginal* tax rates\n from filed tax returns (NBER)", col = "blue", cex = 1, srt = 0)
    text(1980, 0.08, "Average *average* tax rates\n  from filed tax returns (NBER)", col = "black", cex = 1, srt = 0)
})
dev.off()


#print("Plot7")
pdf(file = "tipsInflation.pdf", width = 12*scale, height = 6*scale)
with(d, {
    plot(year, inflation, xaxs = "i", yaxs = "i", xlab = "Year", ylab = "Tax Rate", type = "n",
        cex.axis = 0.8, ylim = c(0, 0.15), xlim = c(2003, 2018))

    lines(year, inflation, lwd = 3)
    points(year, inflation, col = "black", cex = 0.5)

    lines(year, tipsinflation, col = "blue", lwd = 3)
    points(year, tipsinflation, col = "black", cex = 0.5)

    text(2006.5, 0.05, "CPI 3Y Smoothed Inflation", col = "black", cex = 1, srt = 0)
    text(2005.5, 0.01, "TIPS Inflation", col = "blue", cex = 1, srt = 0)
})
dev.off()

pdf(file = "holdings.pdf", width = 12 * scale, height = 6 * scale)
with(d, {
    plot(year, hightaxholdings, xaxs = "i", yaxs = "i", xlab = "Year", ylab = "Percent of Total Private Treasury Holdings", type = "n",
        cex.axis = 0.8, ylim = c(0, 0.5), xlim = c(1989, 2018))

    lines(year, hightaxholdings, lwd = 3)
    points(year, hightaxholdings, col = "black", cex = 0.5)

    lines(year, taxexemptholdings, col = "blue", lwd = 3)
    points(year, taxexemptholdings, col = "black", cex = 0.5)

    lines(year, residualholdings, col = "black", lty = 3, lwd = 1)
    points(year, residualholdings, col = "black", cex = 0.5)

    lines(year, munitax, col = "red", lty = 3, lwd = 1)
    points(year, munitax, col = "black", cex = 0.5)

    #lines(year, interesttax, col = "red", lwd = 3)
    #points(year, interesttax, col = "black", cex = 0.5)

    text(1994, 0.12, "Taxable Holdings", col = "black", cex = 1, srt = 0)
    text(1993, 0.345, "Tax Exempt Holdings", col = "blue", cex = 1, srt = 0)
    text(1994, 0.78, "Unclassified Holdings", col = "black", cex = 1, srt = 0)
    text(2001, 0.32, "Muni Implied Tax Rate", col = "red", cex = 1, srt = 0)

})
dev.off()

pdf(file = "muniCrossover.pdf", width = 12 * scale, height = 6 * scale)
with(d, {
    plot(year, avgordtotal, xaxs = "i", yaxs = "i", xlab = "Year", ylab = "Tax Rate", type = "n",
        cex.axis = 0.8, ylim = c(0, 0.1), xlim = c(2000, 2018))

    lines(year, muni20, lwd = 3)
    points(year, muni20, col = "black", cex = 0.5)

    lines(year, tnote20, col = "blue", lwd = 3)
    points(year, tnote20, col = "black", cex = 0.5)

    lines(year, muni20minusspread, col = "red", lty = 3, lwd = 2)
    points(year, muni20minusspread, col = "black", cex = 0.5)

    #lines(year, muni20 - (corplt - tnote20) * (1 - avgord), col = "red", lty = 3, lwd = 2)
    #points(year, muni20 - (corplt - tnote20) * (1 - avgord), col = "black", cex = 0.5)

    text(2012.25, 0.0475, "20 Yr Municipal Bonds", col = "black", cex = 1, srt = 0)
    text(2013.5, 0.01675, "Credit Adj Municipal", col = "red", cex = 1, srt = 0)
    text(2001.75, 0.065, "20 Yr Treasuries", col = "blue", cex = 1, srt = 0)
})
dev.off()

pdf(file = "muniCrossover CHECK-DO NOT USE.pdf", width = 12 * scale, height = 6 * scale)
with(d, {
    plot(year, avgordtotal, xaxs = "i", yaxs = "i", xlab = "Year", ylab = "Tax Rate", type = "n",
        cex.axis = 0.8, ylim = c(0, 0.1), xlim = c(2000, 2018))

    lines(year, muni20, lwd = 3)
    points(year, muni20, col = "black", cex = 0.5)

    lines(year, tnote20, col = "blue", lwd = 3)
    points(year, tnote20, col = "black", cex = 0.5)

    lines(year, muni20 - (corplt - tnote20) * (1 - interesttax), col = "red", lty = 3, lwd = 2)
    points(year, muni20 - (corplt - tnote20) * (1 - interesttax), col = "black", cex = 0.5)

    text(2012.25, 0.0475, "20 Yr Municipal Bonds", col = "black", cex = 1, srt = 0)
    text(2013.5, 0.01675, "Credit Adj Municipal", col = "red", cex = 1, srt = 0)
    text(2001.75, 0.065, "20 Yr Treasuries", col = "blue", cex = 1, srt = 0)
})
dev.off()
