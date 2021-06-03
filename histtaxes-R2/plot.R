
d <- read.csv("plot.csv")
d <- subset(d, year > 1946)

# Treasuries	tnote1, tnote20
# taxrates	munitax, munitax0credit, avgord, topord
# shortterm	tnote1 tnote1r tnote1at
# longterm	tnote20 tnote20r tnote20at

# appendix creditrates
# appendix InterestTaxes

# pset <- subset(d, TRUE, select= c(year, munitax,munitax0credit, avgord, topord, tnote1, tnote1r, tnote1at, tnote20, tnote20r, tnote20at))
# for (i in 2:ncol(pset)) {  pset[,i] <- sprintf("%.1f", pset[,i]*100) }

# options(max.print=999999999)
# print(pset, sep=" & ")

iaw$pdf.start("Treasuries.pdf", wd=12, ht=6)

with( d, {
    plot( year, tnote1, xaxs="i", yaxs="i", xlab="Year", ylab="Nominal Interest Rate", type="n", cex.axis=0.8, ylim=c(0,0.15) )
    lines( year, tnote1 );
    lines( year, tnote20, col="blue", lwd=2 );
    points( year, tnote1, cex=0.5 );
    points( year, tnote20, col="blue", cex=0.5 );

    iaw$hline( tail( tnote20, 1 ), lty=2, col="blue" )
    iaw$hline( tail( tnote1, 1 ), lty=2, col="black" )

    text( 1995, 0.08, "LT Treasuries", col="blue", cex=0.8, srt=-25 )
    text( 1988.2, 0.046, "ST Treasuries", col="black", cex=0.8, srt=-38 )

    iaw$vline( 1981, lty=3 )
    text( 1981, 0.05, "1981: Highest\nInterest Rates", cex=0.5, col="gray", srt=90 )
})

iaw$pdf.end()


iaw$pdf.start("creditrates.pdf", wd=12, ht=6)
with( d, {
    plot( year, corplt, xaxs="i", yaxs="i", xlab="Year", ylab="Nominal Interest Rate", type="n", cex.axis=0.8, ylim=c(0,0.15) )
    lines( year, corplt );
    lines( year, tnote20, col="blue", lwd=3 );
    lines( year, muni20, col="blue" );

    points( year, corplt, cex=0.5 );
    points( year, tnote20, col="blue", cex=0.5 );
    points( year, muni20, col="blue", cex=0.8, pty=5 );
    iaw$hline( tail( tnote20, 1 ), lty=2 )

    text( 1968, 0.085, "LT AAA Corporates", col="black", cex=0.8, srt=30 )
    text( 1978, 0.05, "LT Municipals", col="blue", cex=0.8, srt=30 )

})
iaw$pdf.end()


iaw$pdf.start("taxrates.pdf", wd=12, ht=6)

with( d, {
    plot( year, topord, xaxs="i", yaxs="i", xlab="Year", ylab="Tax Rate", type="n", cex.axis=0.8, ylim=c(0,1) )

    lines( year, munitax, col="blue", lwd=3 );
    points( year, munitax, col="black", cex=0.5 )

    lines( year, munitax0credit, col="blue", lty=3, lwd=1 );

    lines(year, interesttax, col = "black", lwd = 3)
    points(year, interesttax, col = "black", cex = 0.5)

#    lines( year, avgint, col="black", lwd=3, lty=2 );

    lines( year, topord, col="red", lty=3, lwd=1 );
#    lines( year, topgains, col="red", lwd=1, lty=2 );

    text( 1971, 0.82, "Top statutory rate\non ordinary income", col="red", cex=0.6, srt=0 )
    text( 1955.5, 0.58, "Implied effective tax rate\nfrom credit-adjusted\nMunis - Treasuries", col="blue", cex=0.8, srt=30 )
    text( 1961, 0.36, "Average tax rates\nfrom filed tax\nreturns (NBER)", col="black", cex=0.8, srt=30 )

    text( 1997, 0.05, "Unadjusted  Muni\n-Treasury tax rate", col="blue", cex=0.6, srt=0 )

    iaw$vline(1987, lty=2); text( 1988, 0.7, "Reagan Tax Reform", cex=0.7, srt=90);

})

iaw$pdf.end()



iaw$pdf.start("longterm.pdf", wd=12, ht=6)

with( d, {
    plot( year, tnote20, xaxs="i", yaxs="i", xlab="Year", ylab="Long-Term Interest Rate", type="n", cex.axis=0.8, ylim=c(-0.05,0.15) )
    lines( year, tnote20, col="gray", lty=2 );
    lines( year, tnote20r, col="gray", lty=3 );
    lines( year, tnote20at, col="black", lwd=3 );

    points( year, tnote20, col="gray", cex=0.5 );
    points( year, tnote20r, col="gray", cex=0.5 );
    points( year, tnote20at, col="black", cex=0.5 );

#    iaw$p.arrows( 1979, 0.09, 1979, 0 );
#    iaw$p.arrows( 1979, 0.0, 1979, -0.04 ); text(1979, -0.015, "Avg\nTax", srt=90, cex=0.6)

#    iaw$p.arrows( 1981, 0.13, 1981, 0.04, col="gray" ); text(1981, 0.08, "9% CPI\nInflation", srt=90, cex=0.6, col="gray")
#    iaw$p.arrows( 1981, 0.04, 1981, 0.0, col="gray" ); text(1980, 0.023, "4% Tax\nPayment", srt=90, cex=0.5, col="gray")

#    iaw$hline( tail( tnote20at, 1 ), lty=2, col="blue" )
    iaw$hline( 0, lty=2, col="gray" )

    text( 1995, 0.08, "Nominal Rate", col="gray", cex=0.8, srt=-20 )
    text( 1995, 0.053, "Untaxed Real Rate", col="gray", cex=0.8, srt=-10 )
    text( 1993, 0.01, "Taxed Real Rate", col="black", cex=0.8, srt=0 )
})

iaw$pdf.end()


iaw$pdf.start("shortterm.pdf", wd=12, ht=6)

with( d, {
    plot( year, tnote1, xaxs="i", yaxs="i", xlab="Year", ylab="Short-Term Interest Rate", type="n", cex.axis=0.8, ylim=c(-0.05,0.15) )
    lines( year, tnote1, col="gray", lty=2 );
    lines( year, tnote1r, col="gray", lty=3 );
    lines( year, tnote1at, col="black", lwd=3 );

    points( year, tnote1, col="gray", cex=0.5 );
    points( year, tnote1r, col="gray", cex=0.5 );
    points( year, tnote1at, col="black", cex=0.5 );

#    iaw$p.arrows( 1979, 0.09, 1979, 0 );
#    iaw$p.arrows( 1979, 0.0, 1979, -0.04 ); text(1979, -0.015, "Avg\nTax", srt=90, cex=0.6)

#    iaw$p.arrows( 1981, 0.13, 1981, 0.04, col="gray" ); text(1981, 0.08, "9% CPI\nInflation", srt=90, cex=0.6, col="gray")
#    iaw$p.arrows( 1981, 0.04, 1981, 0.0, col="gray" ); text(1980, 0.023, "4% Tax\nPayment", srt=90, cex=0.5, col="gray")

#    iaw$hline( tail( tnote1at, 1 ), lty=2, col="blue" )
    iaw$hline( 0, lty=2, col="gray" )

    text( 1997, 0.068, "Nominal Rate", col="gray", cex=0.8, srt=-1 )
    text( 1975, 0.01, "Taxed Rate", col="gray", cex=0.8, srt=0 )
    text( 1993, -0.02, "Post-Tax Nominal Rate", col="black", cex=0.8, srt=0 )

})

iaw$pdf.end()




