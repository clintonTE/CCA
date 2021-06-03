import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as sm
import os

def practicepd():
    locdir = r'C:\Users\Clinton\Dropbox\Projects\pythonprep'
    pth=os.path.join(locdir, r'data\crspm.csv')
    crsp = pd.read_csv(pth, sep=',',dtype={
        "PERMNO" : int, 
        "date" : "string", 
        "DLRET" : "string", 
        "RET": "string",
        "SICCD": "category",
        "SHRCD": "category",
        "EXCHCD": "category",
        "DISTCD" : "category",        
        "PRC": "string",
        "DLPRC": "string",
        "SHROUT": "string",
        "ALTPRC": "string",
        "ALTPRCDT": "string",
        })

    crsp.columns = crsp.columns.str.strip().str.lower()

    for c in ["dlprc", "ret", "prc", "dlprc", "shrout", "altprc"]:
        crsp[c] = pd.to_numeric(crsp[c], errors="coerce")  

    
    crsp['date'] = pd.to_datetime(crsp['date'])
    crsp['altprcdt'] = pd.to_datetime(crsp['altprcdt'])

    crsp = crsp.dropna(subset=["ret", "prc", "shrcd", "exchcd", "date"])

    print(crsp.dtypes)
    crsp = crsp[crsp['exchcd'].isin(['1','2','3']) & crsp['shrcd'].isin(['10','11'])]

    crsp['prc1'] = 1/crsp['prc']

    l = sm.ols(formula="ret ~ prc1 + C(exchcd)", data=crsp).fit()
    print(l.params)

    print(crsp.dtypes)


practicepd()
    


