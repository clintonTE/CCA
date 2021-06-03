#allhigh - level scripts
#functions and lower level code are in the below file

rm(list = ls());

require(aod)
require(AER);
require(DataAnalytics);
require(devtools);
require(foreign);
require(ggplot2);
require(Matrix);
require(plyr);
require(reshape2);
require(sandwich);
require(stargazer);
require(xtable);
require(mfx);
#require(compiler);



source("Rep Functions.R");

######switch between debugging/testing to final output
DEBUG = TRUE;
ORIGINAL_WC_METHOD = TRUE;
CUSTOM_ROBUST_ERRORS = FALSE;
UNSOLD_ESTIMATE_DISCOUNT = 0.4;


if (DEBUG == TRUE) {
    STARTYPE = "text";
    TEST = TRUE;
} else {
    STARTYPE = "latex";
    TEST = FALSE;
#    enableJIT(3);
}

#' Title
#'
#' @return
#' @export
#'
#' @examples
#' Title
#'
#' @return
#' @export
#'
#' @examples
main = function() {

    dfC = read.dta("Data\\contempp.dta");

    monthCutoff = 42;

    #########prep for contemporary art initial regression and summaries

    #drop some missing values
    dfC = dfC[!(is.na(dfC$low_est) | is.na(dfC$high_est) | is.na(dfC$price) | (dfC$price == 0) |
        is.na(dfC$date_ptg) | is.na(dfC$len) | is.na(dfC$wid) | is.na(dfC$med) | (dfC$med == "") |
        is.na(dfC$artist) | is.na(dfC$date)),];

        #initial sort
    dfC = dfC[order(dfC$date, dfC$lot),];
    iCurLength = nrow(dfC);

        #scale data appropriately
    dfC[, c("low_est", "high_est", "price")] = dfC[, c("low_est", "high_est", "price")] * 1000;

        #do a bunch of other transformations that they do in their paper
    dfC[, c("lest", "hest", "lprice", "lwid", "llen", "luscpi", "lukcpi", "lav", "age")] =
            cbind(log(dfC[, c("low_est", "high_est", "price", "wid", "len", "uscpi", "ukcpi")]),
            log((dfC$low_est + dfC$high_est) / 2),
            dfC[, "date_ptg"]);
    dfC[, c("llow", "lhigh")] = dfC[, c("lest", "hest")];

        #number the lots
    dfC = dfC[order(dfC$date, dfC$lot, dfC$id),];

        #this is their fudge factor for items which are not sold
    dfC[, "price"] = ifelse(dfC[, "sold"] == 0, dfC$low_est * (1 - UNSOLD_ESTIMATE_DISCOUNT), dfC$price);
    dfC[, "lprice"] = ifelse(dfC[, "sold"] == 0, dfC$llow + log(1 - UNSOLD_ESTIMATE_DISCOUNT), dfC$lprice);
    dfC[, "count"] = sapply(dfC$artist, function(x) sum(dfC$artist == x, na.rm = TRUE));
    dfC[, "lage"] = log(dfC$age - min(dfC$age) + 1);

    dfC = dfC[!(dfC$count < 2),];
    iCurLength = nrow(dfC);

        #date processing section- begin by rordering
    dfC = dfC[order(dfC$id, dfC$date),]

    #core code to identify the sample set
    dfC[, "repSale"] = (dfC[, "id"] == c(NA, dfC[1:(iCurLength - 1), "id"])) & c(NA, dfC[1:(iCurLength - 1), "sold"] > 0);
    dfC[, "repSale"] = ifelse(is.na(dfC[, "repSale"]), FALSE, dfC[, "repSale"]);
    dfC[, "nextIsRepSale"] = c(dfC[2:iCurLength, "repSale"], FALSE);
    dfC[, "nextIsRepSaleSold"] = c(dfC[2:iCurLength, "repSale"], FALSE) & c(dfC[2:iCurLength, "sold"] == 1, FALSE);

    dfC[dfC$repSale, "difdate"] = (dfC[dfC$repSale, "date"] - dfC[dfC$nextIsRepSale, "date"]) / 30;
    dfC[dfC$repSale, "difus"] = -dfC[dfC$repSale, "luscpi"] + dfC[dfC$nextIsRepSale, "luscpi"];
    dfC[dfC$repSale, "difuk"] = -dfC[dfC$repSale, "lukcpi"] + dfC[dfC$nextIsRepSale, "lukcpi"];

    ####################################################output############################

    stargazer(dfC[dfC$sold > 0, c("low_est", "high_est", "price")],
            dfC[, c("low_est", "high_est", "price")],
            dfC[dfC$repSale & dfC$sold > 0, c("low_est", "high_est", "price", "difdate")],
            dfC[dfC$repSale, c("low_est", "high_est", "price", "difdate")],
            type = STARTYPE, single.row = TRUE,
            digits = 0,
            style = "aer", align = TRUE, no.space = TRUE,
            title = c("Summary: Paintings Sold", "Summary: All Paintings", "Summary: Repeated Sales", "Summary: Repeated Listings"),
            omit.summary.stat = c("max", "min"),
            flip = FALSE)

    ##run the first stage regressions
    strHal2Thru25 = paste("hal", 2:25, sep = "");
    sCCovarsStg1 = append(c("lage", "llen", "lwid", "artist", "med"), strHal2Thru25)
    frmCStg1 = formula(paste("lprice ~ ", paste(sCCovarsStg1, collapse = " + ")));

    regCSoldStg1 = lm(frmCStg1, data = dfC[dfC$sold == 1,]);
    regCAllStg1 = lm(frmCStg1, data = dfC);

    #efficiently get the predicted values
    dfC[, "residAll"] = regCAllStg1$residuals;
    dfC[dfC$sold == 1, "residSold"] = regCSoldStg1$residuals;
    dfC[, "ppriceAll"] = dfC[, "lprice"] - regCAllStg1$residuals;
    dfC[dfC$sold == 1, "ppriceSold"] = dfC[dfC$sold == 1, "lprice"] - regCSoldStg1$residuals;

    ##get some lagged stats
    dfC[dfC$repSale, "lresidAll"] = dfC[dfC$nextIsRepSale, "residAll"];
    dfC[dfC$repSale, "llpriceAll"] = dfC[dfC$nextIsRepSale, "lprice"];
    dfC[dfC$repSale, "lppriceAll"] = dfC[dfC$nextIsRepSale, "ppriceAll"];

    dfC[dfC$sold == 1 & dfC$repSale, "lresidSold"] = dfC[dfC$nextIsRepSaleSold, "residSold"];
    dfC[dfC$sold == 1 & dfC$repSale, "llpriceSold"] = dfC[dfC$nextIsRepSaleSold, "lprice"];
    dfC[dfC$sold == 1 & dfC$repSale, "lppriceSold"] = dfC[dfC$nextIsRepSaleSold, "ppriceSold"];

    dfC[dfC$repSale, "lossGainAll"] =
        dfC[dfC$repSale, "llpriceAll"] - dfC[dfC$repSale, "ppriceAll"];
    dfC[dfC$sold == 1 & dfC$repSale, "lossGainSold"] =
        dfC[dfC$sold == 1 & dfC$repSale, "llpriceSold"] - dfC[dfC$sold == 1 & dfC$repSale, "ppriceSold"];

    #now run the second stage stats
    sCCovarStg2Sold = c("lossGainSold", "ppriceSold", "lresidSold", "difdate");
    sCCovarStg2All = c("lossGainAll", "ppriceAll", "lresidAll", "difdate");

    frmCSoldStg2Price = formula(paste("lprice ~ ", paste(sCCovarStg2Sold, collapse = " + ")));
    frmCSoldStg2Low = formula(paste("llow ~ ", paste(sCCovarStg2Sold, collapse = " + ")));
    frmCAllStg2Price = formula(paste("lprice ~ ", paste(sCCovarStg2All, collapse = " + ")));
    frmCAllStg2Low = formula(paste("llow ~ ", paste(sCCovarStg2All, collapse = " + ")));

    regCSoldStg2Price = lm(frmCSoldStg2Price, data = dfC[(dfC$sold == 1) & (dfC$repSale),]);
    regCSoldStg2Low = lm(frmCSoldStg2Low, data = dfC[dfC$sold == 1 & dfC$repSale,]);
    regCSoldStg2PriceCut = lm(frmCSoldStg2Price, data = dfC[dfC$sold == 1 & dfC$repSale & dfC$difdate < monthCutoff,]);
    regCSoldStg2LowCut = lm(frmCSoldStg2Low, data = dfC[dfC$sold == 1 & dfC$repSale & dfC$difdate < monthCutoff,]);

    regCAllStg2Price = lm(frmCAllStg2Price, data = dfC[dfC$repSale,]);
    regCAllStg2Low = lm(frmCAllStg2Low, data = dfC[dfC$repSale,]);
    regCAllStg2PriceCut = lm(frmCAllStg2Price, data = dfC[dfC$repSale & dfC$difdate < monthCutoff,]);
    regCAllStg2LowCut = lm(frmCAllStg2Low, data = dfC[dfC$repSale & dfC$difdate < monthCutoff,]);

    if (CUSTOM_ROBUST_ERRORS) {
    errCSoldStg2Price = GetMWhiteErrors(dfC[(dfC$sold == 1) & (dfC$repSale), sCCovarStg2Sold], regCSoldStg2Price$residuals)
    errCSoldStg2Low = GetMWhiteErrors(dfC[(dfC$sold == 1) & (dfC$repSale), sCCovarStg2Sold], regCSoldStg2Low$residuals)
    errCSoldStg2PriceCut = GetMWhiteErrors(dfC[dfC$sold == 1 & dfC$repSale & dfC$difdate < monthCutoff, sCCovarStg2Sold], regCSoldStg2PriceCut$residuals)
    errCSoldStg2LowCut = GetMWhiteErrors(dfC[dfC$sold == 1 & dfC$repSale & dfC$difdate < monthCutoff, sCCovarStg2Sold], regCSoldStg2LowCut$residuals)

    errCAllStg2Price = GetMWhiteErrors(dfC[dfC$repSale, sCCovarStg2All], regCAllStg2Price$residuals)
    errCAllStg2Low = GetMWhiteErrors(dfC[dfC$repSale, sCCovarStg2All], regCAllStg2Low$residuals)
    errCAllStg2PriceCut = GetMWhiteErrors(dfC[dfC$repSale & dfC$difdate < monthCutoff, sCCovarStg2All], regCAllStg2PriceCut$residuals)
    errCAllStg2LowCut = GetMWhiteErrors(dfC[dfC$repSale & dfC$difdate < monthCutoff, sCCovarStg2All], regCAllStg2LowCut$residuals)
    } else { #test using sandwich package
        errCSoldStg2Price = diag(vcovHC(regCSoldStg2Price, type = "HC1"))^.5;
        errCSoldStg2Low = diag(vcovHC(regCSoldStg2Low, type = "HC1")) ^ .5;
        errCSoldStg2PriceCut = diag(vcovHC(regCSoldStg2PriceCut, type = "HC1")) ^ .5;
        errCSoldStg2LowCut = diag(vcovHC(regCSoldStg2LowCut, type = "HC1")) ^ .5;

        errCAllStg2Price = diag(vcovHC(regCAllStg2Price, type = "HC1")) ^ .5;
        errCAllStg2Low = diag(vcovHC(regCAllStg2Low, type = "HC1")) ^ .5;
        errCAllStg2PriceCut = diag(vcovHC(regCAllStg2PriceCut, type = "HC1")) ^ .5;
        errCAllStg2LowCut = diag(vcovHC(regCAllStg2LowCut, type = "HC1")) ^ .5;
    }

    #Probit

    frmCAllStg2PB = formula(paste("sold ~ ", paste(sCCovarStg2All, collapse = " + ")));
    regCAllStg2PB = glm(frmCAllStg2PB, data = dfC[dfC$repSale,], family = binomial(link = "probit"));
    regCAllStg2PBCut = glm(frmCAllStg2PB, data = dfC[dfC$repSale & dfC$difdate < monthCutoff,], family = binomial(link = "probit"));

    print(wald.test(b = coef(regCAllStg2PB), Sigma = vcovHC(regCAllStg2PB, "HC1"), Terms = 1:5, verbose = TRUE));
    #print(wald.test(b = coef(regCAllStg2PBCut), Sigma = errCAllStg2PBCut, Terms = 1:5));


    ##I need to count the dummy variables, as stargazer can't do it for me
    lDummyStats = list(c("Half-Year Fixed Effects",
        sum(apply(dfC[dfC$sold == 1, c("hal1", strHal2Thru25)], 2, max)),
        sum(apply(dfC[, c("hal1", strHal2Thru25)], 2, max))),
    c("Artist Fixed Effects",
            NumTypes(dfC[dfC$sold == 1, "artist"]),
            NumTypes(dfC[, "artist"])),
    c("Medium Fixed Effects",
            NumTypes(dfC[dfC$sold == 1, "med"]),
            NumTypes(dfC[, "med"]))
        );


    stargazer(regCSoldStg1, regCAllStg1, type = STARTYPE,
        summary.logical = FALSE,
        omit = c("hal", "artist", "med"),
        covariate.labels = c("Log(Age)", "Log(Length)", "Log(Width)"),
        column.labels = c("Sold", "All"),
        dep.var.labels = c("Log(Price)"),
        omit.stat = c("f", "rsq"),
        #omit.labels = c("Year Fixed Effects", "Artist Fixed Effects", "Medium Fixed Effects"),
        add.lines = lDummyStats,
        title = "Anchoring Effects: Contemporary Art 1st Stage",
        style = "aer", align = TRUE, no.space = TRUE,
        star.char = "",
        df = FALSE);

    stargazer(regCSoldStg2Price, regCSoldStg2PriceCut, regCSoldStg2Low, regCSoldStg2LowCut,
        type = STARTYPE,
        se = list(errCSoldStg2Price, errCSoldStg2PriceCut, errCSoldStg2Low, errCSoldStg2LowCut),
        covariate.labels = c("Anchoring Eff", "Predicted Price", "Sale Residual", "Months Since Sale"),
        column.labels = c("Sold", "Sold <42M", "Sold", "Sold <42M"),
        dep.var.labels = c("Log(Price)", "Log(Low Est.)"),
        omit.stat = c("f", "rsq"),
        star.char = "",
        notes = "Robust errors",
        title = "Anchoring Effects: Contemporary Art- Sold",
        style = "aer", align = TRUE, no.space = TRUE,
        df = FALSE);

    stargazer(regCAllStg2Price, regCAllStg2PriceCut, regCAllStg2Low, regCAllStg2LowCut,
        type = STARTYPE,
        se = list(errCAllStg2Price, errCAllStg2PriceCut, errCAllStg2Low, errCAllStg2LowCut),
        covariate.labels = c("Anchoring Eff", "Predicted Price", "Sale Residual", "Months Since Sale"),
        column.labels = c("Listed", "Listed <42M", "Listed", "Listed <42M"),
        dep.var.labels = c("Log(Price)", "Log(Low Est.)"),
        omit.stat = c("f", "rsq"),
        star.char = "",
        notes = "Prices of unsold calculated at 80% of low estimate, robust errors",
        title = "Anchoring Effects: Contemporary Art- All Listed",
        style = "aer", align = TRUE, no.space = TRUE,
        df = FALSE);


    ##############################Switch to impressionist art##################
    dfI = read.dta("Data\\impresss.dta");

    #initial data transformations
    dfI[, c("llow", "lhigh", "lprice", "sold", "house", "loc")] = cbind(
        log(dfI[, c("LOW_EST", "HIGH_EST", "S_PRICE")]),
        ifelse(!is.na(dfI$S_PRICE), 1, 0),
        ifelse((!is.na(dfI$HOUSE)) & dfI$HOUSE == "SO", 1, 0),
        ifelse((!is.na(dfI$AUCT_LOC)) & dfI$AUCT_LOC == "NY", 1, 0));

    dfI[, "newsold"] = sapply(dfI$AUCTION, function(x) mean(dfI[dfI$AUCTION == x, "sold"]));

    dfI = dfI[order(dfI$id, dfI$date, dfI$HOUSE, dfI$AUCT_LOC),];
    iCurLength = nrow(dfI);

    #need some early logic for lags and leads
    dfI[, "idEqIDm1"] = dfI[, "id"] == c(NA, dfI[1:(iCurLength - 1), "id"]);
    dfI[, "idEqIDm1"] = ifelse(!is.na(dfI[, "idEqIDm1"]), dfI[, "idEqIDm1"], FALSE);

    dfI[, "idEqIDp1"] = dfI[, "id"] == c(dfI[2:iCurLength, "id"], NA);
    dfI[, "idEqIDp1"] = ifelse(!is.na(dfI[, "idEqIDp1"]), dfI[, "idEqIDp1"], FALSE);


    #some DQ code in the stata file
    dfI[dfI$idEqIDm1, "ARTIST"] = ifelse(!is.na(dfI[dfI$idEqIDm1, "ARTIST"]), dfI[dfI$idEqIDm1, "ARTIST"], dfI[dfI$idEqIDp1, "ARTIST"]);
    dfI[dfI$idEqIDm1, "ART_MED"] = ifelse(!is.na(dfI[dfI$idEqIDm1, "ART_MED"]), dfI[dfI$idEqIDm1, "ART_MED"], dfI[dfI$idEqIDp1, "ART_MED"]);
    dfI[dfI$idEqIDm1, "DIM_A"] = ifelse(!is.na(dfI[dfI$idEqIDm1, "DIM_A"]), dfI[dfI$idEqIDm1, "DIM_A"], dfI[dfI$idEqIDp1, "DIM_A"]);
    dfI[dfI$idEqIDm1, "DIM_B"] = ifelse(!is.na(dfI[dfI$idEqIDm1, "DIM_B"]), dfI[dfI$idEqIDm1, "DIM_B"], dfI[dfI$idEqIDp1, "DIM_B"]);
    dfI[dfI$idEqIDm1, "SIGNED"] = ifelse(!is.na(dfI[dfI$idEqIDm1, "SIGNED"]), dfI[dfI$idEqIDm1, "SIGNED"], dfI[dfI$idEqIDp1, "SIGNED"]);
    dfI[dfI$idEqIDm1, "DATE_PTG"] = ifelse(!is.na(dfI[dfI$idEqIDm1, "DATE_PTG"]), dfI[dfI$idEqIDm1, "DATE_PTG"], dfI[dfI$idEqIDp1, "DATE_PTG"]);

    dfI[dfI$idEqIDp1, "ARTIST"] = ifelse(!is.na(dfI[dfI$idEqIDp1, "ARTIST"]), dfI[dfI$idEqIDp1, "ARTIST"], dfI[dfI$idEqIDm1, "ARTIST"]);
    dfI[dfI$idEqIDp1, "ART_MED"] = ifelse(!is.na(dfI[dfI$idEqIDp1, "ART_MED"]), dfI[dfI$idEqIDp1, "ART_MED"], dfI[dfI$idEqIDm1, "ART_MED"]);
    dfI[dfI$idEqIDp1, "DIM_A"] = ifelse(!is.na(dfI[dfI$idEqIDp1, "DIM_A"]), dfI[dfI$idEqIDp1, "DIM_A"], dfI[dfI$idEqIDm1, "DIM_A"]);
    dfI[dfI$idEqIDp1, "DIM_B"] = ifelse(!is.na(dfI[dfI$idEqIDp1, "DIM_B"]), dfI[dfI$idEqIDp1, "DIM_B"], dfI[dfI$idEqIDm1, "DIM_B"]);
    dfI[dfI$idEqIDp1, "SIGNED"] = ifelse(!is.na(dfI[dfI$idEqIDp1, "SIGNED"]), dfI[dfI$idEqIDp1, "SIGNED"], dfI[dfI$idEqIDm1, "SIGNED"]);
    dfI[dfI$idEqIDp1, "DATE_PTG"] = ifelse(!is.na(dfI[dfI$idEqIDp1, "DATE_PTG"]), dfI[dfI$idEqIDp1, "DATE_PTG"], dfI[dfI$idEqIDm1, "DATE_PTG"]);

    dfI[, "AUCTION"] = ifelse(!is.na(dfI[, "AUCTION"]), dfI[, "AUCTION"], dfI[, "date"]);

    #data conformity/more DQ
    dfI[, "med"] = as.factor(dfI$ART_MED);
    levels(dfI$med) = c(levels(dfI$med), "1001", "1002", "1003")
    dfI[is.na(dfI$med), "med"] = "1001";

    ####NOTE: THE FOLLOWING IS A CODING ERROR ON THEIR PART- MEDIUM 39 IS LISTED AS  WATERCOLOR
    #### BUT THEY ALSO CONFORM IT TO OIL- THIS IS WRONG- HENCE WE WILL DO THIS TWICE
    if (ORIGINAL_WC_METHOD) {
        dfI[, "wcDum"] = dfI$med == "39";
        dfI[dfI$med %in% c("18", "21", "39"), "med"] = "1002";
    } else {

        dfI[dfI$med %in% c("18", "21"), "med"] = "1002";

    }

    dfI[dfI$med %in% c("1", "7", "10", "13", "14",
                        "16", "17", "19", "20", "22",
                        "26", "28", "29", "31", "35",
                        "37", "38", "61", "700"), "med"] = "1001";

    dfI[, "artist"] = dfI$ANAME;
    dfI[is.na(dfI$artist), "artist"] = "UNK";
    dfI[dfI$artist %in% c("ERTEROMADET", "ERTEROMANDET"), "artist"] = "Erte";

    dfI[, "signed"] = as.factor(dfI[, "SIGNED"]);
    levels(dfI$signed) = c(levels(dfI$signed), "0");
    dfI[is.na(dfI$signed) | !(dfI$signed %in% c("1", "2", "3")), "signed"] = "0";

    dfI$DATE_PTG = dfI$DATE_PTG - min(dfI$DATE_PTG, na.rm = TRUE) + 1;
    dfI[, c("ldate", "llen", "lwid")] = c(log(dfI$DATE_PTG), log(dfI$DIM_A), log(dfI$DIM_B));

    dfI[, "LO"] = ifelse((!is.na(dfI$AUCT_LOC)) & dfI$AUCT_LOC == "LO", 1, 0);
    dfI[, "CH"] = ifelse((!is.na(dfI$AUCT_LOC)) & dfI$AUCT_LOC == "CH", 1, 0);

    dfI[, "lav"] = log((dfI$LOW_EST + dfI$HIGH_EST) / 2);

    #did they split out the sold form unsold correctly?
    dfI[dfI$sold == 0, "lprice"] = dfI[dfI$sold == 0, "llow"] + log(1 - UNSOLD_ESTIMATE_DISCOUNT);
    dfI[dfI$sold == 0, "S_PRICE"] = dfI[dfI$sold == 0, "LOW_EST"] * (1 - UNSOLD_ESTIMATE_DISCOUNT);
    #*************
    dfI = dfI[!(is.na(dfI$ldate) | is.na(dfI$llen) | is.na(dfI$lwid) | (dfI$hal23 == 1) |
        is.na(dfI$lprice) | dfI$DIM_A == 0 | is.na(dfI$DIM_A) | dfI$DIM_B == 0 | is.na(dfI$DIM_B)),];
    iCurLength = nrow(dfI);

    ###start work on the first regression
    strHal2Thru23 = paste("hal", 2:23, sep = "");

    if (ORIGINAL_WC_METHOD) {
        sICovarsStg1 = append(c("ldate", "llen", "lwid", "signed", "med", "artist", "wcDum"), strHal2Thru23)
    } else {
        sICovarsStg1 = append(c("ldate", "llen", "lwid", "signed", "med", "artist"), strHal2Thru23)
    }
    frmIStg1 = formula(paste("lprice ~ ", paste(sICovarsStg1, collapse = " + ")));


    regISoldStg1LO = lm(frmIStg1, data = dfI[dfI$sold == 1 & dfI$LO == 1,]);
    regISoldStg1NY = lm(frmIStg1, data = dfI[dfI$sold == 1 & dfI$LO == 0,]);
    regIAllStg1LO = lm(frmIStg1, data = dfI[dfI$LO == 1,]);
    regIAllStg1NY = lm(frmIStg1, data = dfI[dfI$LO == 0,]);

    #currency conversion
    dfI[dfI$LO == 1, c("LOW_EST", "HIGH_EST", "S_PRICE")] =
    dfI[dfI$LO == 1, c("LOW_EST", "HIGH_EST", "S_PRICE")] * dfI[dfI$LO == 1, "CNV_RATE"];

    #get the predicted values
    dfI[dfI$sold == 1 & dfI$LO == 1, "residSoldLO"] = regISoldStg1LO$residuals;
    dfI[dfI$sold == 1 & dfI$LO == 0, "residSoldNY"] = regISoldStg1NY$residuals;
    dfI[dfI$LO == 1, "residAllLO"] = regIAllStg1LO$residuals
    dfI[dfI$LO == 0, "residAllNY"] = regIAllStg1NY$residuals;

    dfI[dfI$sold == 1, "residSold"] = ifelse(dfI[dfI$sold == 1, "LO"] == 1, dfI[dfI$sold == 1, "residSoldLO"], dfI[dfI$sold == 1, "residSoldNY"]);
    dfI[, "residAll"] = ifelse(dfI[, "LO"] == 1, dfI[, "residAllLO"], dfI[, "residAllNY"]);

    dfI[dfI$sold == 1, "ppriceSold"] = dfI[dfI$sold == 1, "lprice"] - dfI[dfI$sold == 1, "residSold"];
    dfI[, "ppriceAll"] = dfI[, "lprice"] - dfI[, "residAll"];

    #JIC sort
    dfI = dfI[order(dfI$id, dfI$date, dfI$HOUSE, dfI$AUCT_LOC),];

    #Core code to identify the sampel set (note locations must be equal
    dfI[, "repSale"] = (dfI[, "id"] == c(NA, dfI[1:(iCurLength - 1), "id"])) & c(NA, dfI[1:(iCurLength - 1), "sold"] > 0) &
        dfI[, "LO"] == c(NA, dfI[1:(iCurLength - 1), "LO"]);
    dfI[, "repSale"] = ifelse(is.na(dfI[, "repSale"]), FALSE, dfI[, "repSale"]);
    dfI[, "nextIsRepSale"] = c(dfI[2:iCurLength, "repSale"], FALSE);
    dfI[, "nextIsRepSaleSold"] = c(dfI[2:iCurLength, "repSale"], FALSE) & c(dfI[2:iCurLength, "sold"] == 1, FALSE);

    ##generate our contingent variables
    dfI[dfI$repSale, "difdate"] = (dfI[dfI$repSale, "date"] - dfI[dfI$nextIsRepSale, "date"]) / 30;

    dfI[dfI$repSale, "lresidSold"] = dfI[dfI$nextIsRepSale, "residSold"];
    dfI[dfI$repSale, "lresidAll"] = dfI[dfI$nextIsRepSale, "residAll"];

    dfI[dfI$repSale, "lossGainSold"] = dfI[dfI$nextIsRepSale, "lprice"] - dfI[dfI$repSale, "ppriceSold"];
    dfI[dfI$repSale, "lossGainAll"] = dfI[dfI$nextIsRepSale, "lprice"] - dfI[dfI$repSale, "ppriceAll"];


    sICovarStg2Sold = c("lossGainSold", "ppriceSold", "lresidSold", "difdate");
    sICovarStg2All = c("lossGainAll", "ppriceAll", "lresidAll", "difdate");


    frmISoldStg2Price = formula(paste("lprice ~ ", paste(sICovarStg2Sold, collapse = " + ")));
    frmISoldStg2Low = formula(paste("llow ~ ", paste(sICovarStg2Sold, collapse = " + ")));
    frmIAllStg2Price = formula(paste("lprice ~ ", paste(sICovarStg2All, collapse = " + ")));
    frmIAllStg2Low = formula(paste("llow ~ ", paste(sICovarStg2All, collapse = " + ")));

    regISoldStg2Price = lm(frmISoldStg2Price, data = dfI[(dfI$sold == 1) & (dfI$repSale),]);
    regISoldStg2Low = lm(frmISoldStg2Low, data = dfI[dfI$sold == 1 & dfI$repSale,]);
    regISoldStg2PriceCut = lm(frmISoldStg2Price, data = dfI[dfI$sold == 1 & dfI$repSale & dfI$difdate < monthCutoff,]);
    regISoldStg2LowCut = lm(frmISoldStg2Low, data = dfI[dfI$sold == 1 & dfI$repSale & dfI$difdate < monthCutoff,]);

    regIAllStg2Price = lm(frmIAllStg2Price, data = dfI[dfI$repSale,]);
    regIAllStg2Low = lm(frmIAllStg2Low, data = dfI[dfI$repSale,]);
    regIAllStg2PriceCut = lm(frmIAllStg2Price, data = dfI[dfI$repSale & dfI$difdate < monthCutoff,]);
    regIAllStg2LowCut = lm(frmIAllStg2Low, data = dfI[dfI$repSale & dfI$difdate < monthCutoff,]);

    if (CUSTOM_ROBUST_ERRORS) { 

        errISoldStg2Price = GetMWhiteErrors(dfI[(dfI$sold == 1) & (dfI$repSale), sICovarStg2Sold], regISoldStg2Price$residuals)
        errISoldStg2Low = GetMWhiteErrors(dfI[(dfI$sold == 1) & (dfI$repSale), sICovarStg2Sold], regISoldStg2Low$residuals)
        errISoldStg2PriceCut = GetMWhiteErrors(dfI[dfI$sold == 1 & dfI$repSale & dfI$difdate < monthCutoff, sICovarStg2Sold], regISoldStg2PriceCut$residuals)
        errISoldStg2LowCut = GetMWhiteErrors(dfI[dfI$sold == 1 & dfI$repSale & dfI$difdate < monthCutoff, sICovarStg2Sold], regISoldStg2LowCut$residuals)

        errIAllStg2Price = GetMWhiteErrors(dfI[dfI$repSale, sICovarStg2All], regIAllStg2Price$residuals)
        errIAllStg2Low = GetMWhiteErrors(dfI[dfI$repSale, sICovarStg2All], regIAllStg2Low$residuals)
        errIAllStg2PriceCut = GetMWhiteErrors(dfI[dfI$repSale & dfI$difdate < monthCutoff, sICovarStg2All], regIAllStg2PriceCut$residuals)
        errIAllStg2LowCut = GetMWhiteErrors(dfI[dfI$repSale & dfI$difdate < monthCutoff, sICovarStg2All], regIAllStg2LowCut$residuals)
    } else {
        errISoldStg2Price = diag(vcovHC(regISoldStg2Price, type = "HC1")) ^ .5;
        errISoldStg2Low = diag(vcovHC(regISoldStg2Low, type = "HC1")) ^ .5;
        errISoldStg2PriceCut = diag(vcovHC(regISoldStg2PriceCut, type = "HC1")) ^ .5;
        errISoldStg2LowCut = diag(vcovHC(regISoldStg2LowCut, type = "HC1")) ^ .5;

        errIAllStg2Price = diag(vcovHC(regIAllStg2Price, type = "HC1")) ^ .5;
        errIAllStg2Low = diag(vcovHC(regIAllStg2Low, type = "HC1")) ^ .5;
        errIAllStg2PriceCut = diag(vcovHC(regIAllStg2PriceCut, type = "HC1")) ^ .5;
        errIAllStg2LowCut = diag(vcovHC(regIAllStg2LowCut, type = "HC1")) ^ .5;
    }



    stargazer(dfI[dfI$sold == 1 & dfI$LO == 1, c("LOW_EST", "HIGH_EST", "S_PRICE")],
        dfI[dfI$sold == 1 & dfI$LO == 0, c("LOW_EST", "HIGH_EST", "S_PRICE")],
        dfI[dfI$LO == 1, c("LOW_EST", "HIGH_EST", "S_PRICE")],
        dfI[dfI$LO == 0, c("LOW_EST", "HIGH_EST", "S_PRICE")],
        dfI[dfI$repSale & dfI$sold == 1 & dfI$LO == 1, c("LOW_EST", "HIGH_EST", "S_PRICE", "difdate")],
        dfI[dfI$repSale & dfI$sold == 1 & dfI$LO == 0, c("LOW_EST", "HIGH_EST", "S_PRICE", "difdate")],
        dfI[dfI$repSale & dfI$LO == 1, c("LOW_EST", "HIGH_EST", "S_PRICE", "difdate")],
        dfI[dfI$repSale & dfI$LO == 0, c("LOW_EST", "HIGH_EST", "S_PRICE", "difdate")],
        type = STARTYPE, single.row = TRUE, digits = 0, style = "aer",
        align = TRUE, no.space = TRUE,
        title = c("Summary: London, Sold", "Summary: NY, Sold", "Summary: London, All", "Summary: NY, All",
        "Summary: London, Rep. Sales", "Summary: NY, Rep. Sales", "Summary: London, Rep. Listings", "Summary: NY, Rep. Listings"),
        omit.summary.stat = c("max", "min"), flip = FALSE);

lDummyStatsI = list(c("Half-Year Fixed Effects",
        sum(apply(dfI[dfI$sold == 1 & dfI$LO == 1, c("hal1", strHal2Thru23)], 2, max)),
        sum(apply(dfI[dfI$sold == 1 & dfI$LO == 0, c("hal1", strHal2Thru23)], 2, max)),
        sum(apply(dfI[dfI$LO == 1, c("hal1", strHal2Thru23)], 2, max)),
        sum(apply(dfI[dfI$LO == 0, c("hal1", strHal2Thru23)], 2, max))),
    c("Artist Fixed Effects",
            NumTypes(dfI[dfI$sold == 1 & dfI$LO == 1, "artist"]),
            NumTypes(dfI[dfI$sold == 1 & dfI$LO == 0, "artist"]),
            NumTypes(dfI[dfI$LO == 1, "artist"]),
            NumTypes(dfI[dfI$LO == 0, "artist"])),
    c("Medium Fixed Effects",
            NumTypes(dfI[dfI$sold == 1 & dfI$LO == 1, "med"]),
            NumTypes(dfI[dfI$sold == 1 & dfI$LO == 0, "med"]),
            NumTypes(dfI[dfI$LO == 1, "med"]),
            NumTypes(dfI[dfI$LO == 0, "med"]))
        );



stargazer(regISoldStg1LO, regISoldStg1NY, regIAllStg1LO, regIAllStg1NY,
    type = STARTYPE,
    summary.logical = FALSE,
    omit = c("hal", "artist", "med", "wcD"),
    covariate.labels = c("Log(Age)", "Log(Length)", "Log(Width)", "Signed", "Monagrammed", "Stamped"),
    column.labels = c("LO: Sold", "NY: Sold", "LO: All", "NY: All"),
    dep.var.labels = c("Log(Price)"),
    omit.stat = c("f"),
    #omit.labels = c("Year Fixed Effects", "Artist Fixed Effects", "Medium Fixed Effects"),
    add.lines = lDummyStatsI,
    title = "Anchoring Effects: Impressionist Art 1st Stage",
    style = "aer", align = TRUE, no.space = TRUE,
    star.char = "",
    df = FALSE);
    
    stargazer(regISoldStg2Price, regISoldStg2PriceCut, regISoldStg2Low, regISoldStg2LowCut,
            type = STARTYPE,
            se = list(errISoldStg2Price, errISoldStg2PriceCut, errISoldStg2Low, errISoldStg2LowCut),
            covariate.labels = c("Anchoring Eff", "Predicted Price", "Sale Residual", "Months Since Sale"),
            column.labels = c("Sold", "Sold <42M", "Sold", "Sold <42M"),
            dep.var.labels = c("Log(Price)", "Log(Low Est.)"),
            omit.stat = c("f", "rsq"),
            star.char = "",
            notes = "Robust errors",
            title = "Anchoring Effects: Impressionist Art- Sold",
            style = "aer", align = TRUE, no.space = TRUE,
            df = FALSE);

    stargazer(regIAllStg2Price, regIAllStg2PriceCut, regIAllStg2Low, regIAllStg2LowCut,
            type = STARTYPE,
            se = list(errIAllStg2Price, errIAllStg2PriceCut, errIAllStg2Low, errIAllStg2LowCut),
            covariate.labels = c("Anchoring Eff", "Predicted Price", "Sale Residual", "Months Since Sale"),
            column.labels = c("Listed", "Listed <42M", "Listed", "Listed <42M"),
            dep.var.labels = c("Log(Price)", "Log(Low Est.)"),
            omit.stat = c("f", "rsq"),
            star.char = "",
            notes = "Prices of unsold calculated at 80% of low estimate, robust errors",
            title = "Anchoring Effects: Impressionist Art- All Listed",
            style = "aer", align = TRUE, no.space = TRUE,
            df = FALSE);
    #Probit

    frmIAllStg2PB = formula(paste("sold ~ ", paste(sICovarStg2All, collapse = " + ")));
    regIAllStg2PB = glm(frmIAllStg2PB, data = dfI[dfI$repSale,], family = binomial(link = "probit"));
    regIAllStg2PBCut = glm(frmIAllStg2PB, data = dfI[dfI$repSale & dfI$difdate < monthCutoff,],
        family = binomial(link = "probit"));

   
    print(wald.test(b = coef(regIAllStg2PB), Sigma = vcovHC(regIAllStg2PB, "HC3"), Terms = 1:5, verbose = TRUE));
}


main();








