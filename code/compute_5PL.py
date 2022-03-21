import pandas as pd
import numpy as np
from scipy.optimize import least_squares


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


def compute_samples(df, params):
    '''
    calculate the expected controls/samples from 5PL fit and compare to bounds from
    luminex params
    Coff is the relative distance of computed conc from range limited by Cmean and Cbound
    '''

    # compute the concentration from the fit for all non-standard samples
    df.loc[df['conc'] != df['conc'], "conc"] = retro_5PL(df['FI'], params)
    df.loc[:, "bound"] = np.log(df['Cmax']/df['Cmin']) / 2
    df['Coff'] = (np.log((df['conc'] + .1) / df['Cmin']) - df['bound']) / df['bound']
    return df


def get_confidence(params, fraction=0.9):
    '''
    returns Cmin and Cmax (as onfidence range around that mean for given 5PL-params)
    '''
    Frange = (params[1] - params[0])
    FoffSet = Frange * (1-fraction) / 2
    Fmin = params[0] + FoffSet
    Fmax = params[1] - FoffSet
    return list(retro_5PL(pd.DataFrame([Fmin, Fmax]), params)[0])
