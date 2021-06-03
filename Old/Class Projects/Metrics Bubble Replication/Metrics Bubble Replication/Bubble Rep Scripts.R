#all high-level scripts
TEST = TRUE;

source("Bubble Rep Functions.R")

dat = ReadBubbleData()

if (TEST = TRUE) {
    print(head(dat$dfDaily))
    print(head(dat$dfInstVolume))
    print(head(dat$dfIntraday))
    print(head(dat$dfMaturity))
}

