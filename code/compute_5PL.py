import pandas as pd
import numpy as np
from scipy.optimize import least_squares
import re

def get_standard(data_df, col_df, run="20211021", gene='M-CSF', dilution=4):
    '''
    returns the standard for that run and the respective gene
    '''
    
    # retrieve only the controls from the raw data
    control_df = data_df.loc[data_df['Type'].fillna("").str.match(r"^[SC][1-9]?"),:]
    
    run = int(run)
    # retrieve standards and control from control_df
    s = control_df.query('Run == @run and Gene == @gene')
    ss = s.loc[s['Type'].str.match(r"^S[1-8]$"),:]
    sc = s.loc[s['Type'].str.match("^C[12]$"),:]
    
    # get the starting concentration for that dilution
    conc = col_df.query('Run == @run and Gene == @gene')['S1'].iloc[0]
    # fill the dilution series with the last being 0
    ss.loc[:, 'conc'] = conc / np.power(dilution, ss.loc[:, 'Type'].str.extract("S([1-8])").astype(int)-1)
    ss.loc[ss['Type'] == "S8", 'conc'] = 0.001
    return ss.loc[:, ['Run', 'Gene', 'conc', 'FI']]


def PL5(conc,params):
    '''
    equation for regression
    f(conc) = FI
    '''
    
    A,B,C,D,E = params
    return A + (B-A)/np.power(1+np.power(conc/C, D), E)

def retro_5PL(fi, params):
    '''
    equation for deriving conc from FI
    '''
    
    A,B,C,D,E = params
    conc = C * np.power(np.power((B-A)/(fi-A), 1/E)-1, 1/D).fillna(0)

    return np.round(conc,2)
    
    
def residuals(params, x, y):
    return y - PL5(x,params)


def r_squared(params, s):
    '''
    regression coefficient
    '''
    
    y_mean = np.mean(s['FI'])
    
    res_ss = np.sum(residuals(params, s['conc'], s['FI']) ** 2)
    tot_ss = np.sum((s['FI'] - y_mean)**2)
    return 1 - (res_ss / tot_ss)
    

def fit_standard(s):
    '''
    takes a standard with columns conc and FI 
    and returns the params for the 5LP regression
    '''
    
    # init the params to sensible values (domain-specific)
    p0 = [10, 1000, 10000, -1, 1]
    
    # fit using leastsq
    plsq = least_squares(residuals, p0, jac="2-point", method="trf", args=(s['conc'], s['FI']),bounds=([-np.inf,-np.inf,-np.inf,-1,0.1],[np.inf, np.inf, np.inf,0,5]))
    
    params = list(plsq['x'])
    
    return params, r_squared(params, s)


def fit_curve(s):
    '''
    returns the fit curve for plotting
    '''
    params, R = fit_standard(s)
    conc = np.power(10, np.linspace(-3, 5, 10000))
    fit = PL5(conc, params)

    return conc, fit, R, params


def get_params_from_string(string):
    '''
    little helper
    extracts the params from the luminex curve fit string
    "Std. Curve: FI = 30,2261 + (14539,8 - 30,2261) / ((1 + (Conc / 678,757)^-1,04439))^1,87946"
    --> 
    '''
    
    A,B,_,C,D,E = [float(f) for f in re.findall(r"-?[0-9]+\.[0-9]+", string.replace(",", "."))]
    return [A,B,C,D,E]