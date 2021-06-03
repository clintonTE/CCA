
setwd("C:\\Users\\Clinton Tepper\\Dropbox\\Projects\\histtaxes-R2\\histtaxes-R2")
d <- read.csv("PlotRep.csv", check.names = FALSE)
#print(colnames(d))
d <- subset(d, year > 1950)

scale <- 0.75

pdf(file = "JC IAW Rep Graph 20Y.pdf", width = 10 * scale, height = 8 * scale)
with(d, {
    plot(year, avgordtotal, xaxs = "i", yaxs = "i", xlab = "Year", ylab = "Rate", type = "n",
        cex.axis = 0.8, ylim = c(-0.05, 0.09), xlim = c(1950, 2017))


    lines(year, dp, col="blue", lwd = 3)
    points(year, dp, col = "blue", cex = 0.5)

    lines(year, tnote20r, col = "black", lwd = 3)
    points(year, tnote20r, col = "black", cex = 0.5)

    lines(year, tnote20at, col = "red", lwd = 3)
    points(year, tnote20at, col = "red", cex = 0.5)

    text(2013, 0.03, "Div. Yield", col = "blue", cex = 1, srt = 0)
    text(1996, 0.06, "20Y Real (3Y Inflation)", col = "black", cex = 1, srt = 0)
    text(1990, 0.008, "20Y Real After Tax", col = "red", cex = 1, srt = 0)
})
dev.off()

pdf(file = "JC IAW Rep Graph 10Y DP.pdf", width = 10 * scale, height = 8 * scale)
with(d, {
    plot(year, avgordtotal, xaxs = "i", yaxs = "i", xlab = "Year", ylab = "Rate", type = "n",
        cex.axis = 0.8, ylim = c(-0.01, 0.09), xlim = c(1982, 2017))


    lines(year, dp, col = "blue", lwd = 3)
    points(year, dp, col = "blue", cex = 0.5)

    lines(year, upperdp, col = "blue", lty = 3, lwd = 2)
    points(year, upperdp, col = "blue", cex = 0.5)

    lines(year, lowerdp, col = "blue", lty = 3, lwd = 2)
    points(year, lowerdp, col = "blue", cex = 0.5)

    lines(year, clefed, col = "black", lwd = 3)
    points(year, clefed, col = "black", cex = 0.5)

    lines(year, clefedat, col = "red", lwd = 3)
    points(year, clefedat, col = "red", cex = 0.5)

    text(2013, 0.03, "Div. Yield", col = "blue", cex = 1, srt = 0)
    text(1990.5, 0.05, "10Y Real", col = "black", cex = 1, srt = 0)
    text(1990.5, 0.01, "10Y Real After Tax", col = "red", cex = 1, srt = 0)
})
dev.off()

pdf(file = "JC IAW Rep Graph 10Y DP DL.pdf", width = 10 * scale, height = 8 * scale)
with(d, {
    plot(year, avgordtotal, xaxs = "i", yaxs = "i", xlab = "Year", ylab = "Rate", type = "n",
        cex.axis = 0.8, ylim = c(-0.01, 0.09), xlim = c(1982, 2017))


    lines(year, dp, col = "blue", lwd = 3)
    points(year, dp, col = "blue", cex = 0.5)

    lines(year, upperdpdl, col = "blue", lty = 3, lwd = 2)
    points(year, upperdpdl, col = "blue", cex = 0.5)

    lines(year, lowerdpdl, col = "blue", lty = 3, lwd = 2)
    points(year, lowerdpdl, col = "blue", cex = 0.5)

    lines(year, clefed, col = "black", lwd = 3)
    points(year, clefed, col = "black", cex = 0.5)

    lines(year, clefedat, col = "red", lwd = 3)
    points(year, clefedat, col = "red", cex = 0.5)

    text(2013, 0.03, "Div. Yield", col = "blue", cex = 1, srt = 0)
    text(1990.5, 0.05, "10Y Real", col = "black", cex = 1, srt = 0)
    text(1990.5, 0.01, "10Y Real After Tax", col = "red", cex = 1, srt = 0)
})
dev.off()

pdf(file = "JC IAW Rep Graph 10Y EP.pdf", width = 10 * scale, height = 8 * scale)
with(d, {
    plot(year, avgordtotal, xaxs = "i", yaxs = "i", xlab = "Year", ylab = "Rate", type = "n",
        cex.axis = 0.8, ylim = c(-0.01, 0.09), xlim = c(1982, 2017))


    lines(year, ep, col = "blue", lwd = 3)
    points(year, ep, col = "blue", cex = 0.5)

    lines(year, upperep, col = "blue", lty = 3, lwd = 2)
    points(year, upperep, col = "blue", cex = 0.5)

    lines(year, lowerep, col = "blue", lty = 3, lwd = 2)
    points(year, lowerep, col = "blue", cex = 0.5)

    lines(year, clefed, col = "black", lwd = 3)
    points(year, clefed, col = "black", cex = 0.5)

    lines(year, clefedat, col = "red", lwd = 3)
    points(year, clefedat, col = "red", cex = 0.5)

    text(2006, 0.07, "Earnings Yield", col = "blue", cex = 1, srt = 0)
    text(1987.5, 0.05, "10Y Real", col = "black", cex = 1, srt = 0)
    text(1986.5, 0.01, "10Y Real After Tax", col = "red", cex = 1, srt = 0)
})
dev.off()

pdf(file = "JC IAW Rep Graph 10Y EP DL.pdf", width = 10 * scale, height = 8 * scale)
with(d, {
    plot(year, avgordtotal, xaxs = "i", yaxs = "i", xlab = "Year", ylab = "Rate", type = "n",
        cex.axis = 0.8, ylim = c(-0.01, 0.09), xlim = c(1982, 2017))


    lines(year, ep, col = "blue", lwd = 3)
    points(year, ep, col = "blue", cex = 0.5)

    lines(year, upperepdl, col = "blue", lty = 3, lwd = 2)
    points(year, upperepdl, col = "blue", cex = 0.5)

    lines(year, lowerepdl, col = "blue", lty = 3, lwd = 2)
    points(year, lowerepdl, col = "blue", cex = 0.5)

    lines(year, clefed, col = "black", lwd = 3)
    points(year, clefed, col = "black", cex = 0.5)

    lines(year, clefedat, col = "red", lwd = 3)
    points(year, clefedat, col = "red", cex = 0.5)

    text(2006, 0.07, "Earnings Yield", col = "blue", cex = 1, srt = 0)
    text(1987.5, 0.05, "10Y Real", col = "black", cex = 1, srt = 0)
    text(1986.5, 0.01, "10Y Real After Tax", col = "red", cex = 1, srt = 0)
})
dev.off()