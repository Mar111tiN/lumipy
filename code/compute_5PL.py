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
    # for values above Fmax set the highest value
    conc.loc[conc<0] = conc.max()
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
    

def fit_standard(s, B_bound=np.inf, **kwargs):
    '''
    takes a standard with columns conc and FI 
    and returns the params for the 5LP regression
    '''
    
    # init the params to sensible values (domain-specific)
    p0 = [10, 1000, 10000, -1, 1]
    # fit using leastsq
    plsq = least_squares(residuals, p0, jac="2-point", method="trf", args=(s['conc'], s['FI']), bounds=([-np.inf,-np.inf,-np.inf,-1,0.1],[np.inf, B_bound, np.inf,0,5]))
    
    params = list(plsq['x'])
    
    return params, r_squared(params, s)


def compute_conc(df, standard_row, conc_col_suff=""):
    '''
    calculate the expected controls/samples from 5PL fit and compare to bounds from
    luminex params
    computed values are stored in conc_col
    Fpos is the relative distance between Fmin and Fmax as a measure of confidence
    '''


    MINVALUE=0.01

    conc_col = "conc" + conc_col_suff
    Fpos_col = "Fpos" + conc_col_suff
    # extract the params from the standard_row params string
    params = [float(p) for p in standard_row['params'].split(" | ")]
    df.loc[:, conc_col] = retro_5PL(df['FI'], params)
    # extract Fmin and Fmax
    Fmin, Fmax = standard_row.loc[['Fmin', 'Fmax']]
    df.loc[:, Fpos_col] = np.round((df['FI'] - Fmin) / (Fmax - Fmin),3)
    # upgrade 0 values to 0.001
    df.loc[df[conc_col] == 0, conc_col] = MINVALUE
    # distances in the C space should be log-linear
    # Cbound = np.log(Cmin/Cmax) / 2
    # df['Coff'] = (np.log((df[conc_col] + .1) / Cmin) - Cbound) / Cbound
    return df


def get_confidence(params, fraction=0.9):
    '''
    returns a pd.Series of Fmin, Fmax, Cmin and Cmax (as onfidence range around that mean for given 5PL-params)
    this is calculated from the definition of the params
    ?? should the FI confidence be determined from log2 ?
    '''
    Frange = (params[1] - params[0])
    FoffSet = Frange * (1-fraction) / 2
    Fmin = params[0] + FoffSet
    Fmax = params[1] - FoffSet
    conf = retro_5PL(pd.DataFrame([Fmin, Fmax]), params)[0]
    conf.index = ['ConcMin', 'ConcMax']
    # add Fmin and Fmax to the pd.Series
    return pd.concat([pd.Series(dict(Fmin=Fmin, Fmax=Fmax)), conf])


def analyse_standard(standard_row, s, dilution=4, confidence=0.9, **kwargs):
    '''
    take a standard_row and and samples reduced to one protein
    return (all squeeced into the returned standard_row):
        - ConcMin, ConcMax and Fmin and Fmax and C2pos as a confidence interval
        - Fpos as a measure of how well the samples fit into the sigmoidal curve
            (should be between 0 and 1 (optimally around 0.5))
        - StMax as maximum Fpos in the standard dilution series as a measure of control suitability
            this is a measure of the reach of the maximal standard concentrations
            StMax < 0.6 mean the sigmoidal curve is largely extrapolated 
        - the standard_df (ss)
    '''

    # retrieve standard rows from s
    ss = s.loc[s['Type'].str.match(r"^[SB][1-8]?$"),:]
    ss.loc[ss['Type'] == "B", "Type"] = "S8"

    # get the starting concentration for that dilution
    starting_conc = standard_row['S1']

    # fill the dilution series with the last being 0
    ss.loc[:, 'conc'] = starting_conc / np.power(dilution, ss['Type'].str.extract(r"S([1-8])", expand=False).astype(int) -1)
    ss.loc[ss['Type'] == "S8", 'conc'] = 0
    
    # fit the params
    params, R = fit_standard(ss, **kwargs)
    # add the params and fit to the standard_row
    fit_series = pd.Series([" | ".join([str(round(p,3)) for p in params]), round(R, 6)], index=['params','R^2'])
    # add the ConcMin and ConcMax to standard_row
    conf_series = get_confidence(params, fraction=confidence)

    standard_row = pd.concat([standard_row, fit_series, conf_series])

    # do not compute the fit concentrations!!! 
    # ss = compute_conc(ss, standard_row)
    # compute StMax as maximum Fpos in the standard dilution series
    # this is a measure of the reach of the maximal standard concentrations
    # StMax < 0.6 mean the sigmoidal curve is largely extrapolated 
    Fmin, Fmax = standard_row.loc[['Fmin', 'Fmax']]
    ss.loc[:, 'Fpos'] = np.round((ss['FI'] - Fmin) / (Fmax - Fmin),3)
    standard_row["StMax"] = np.round(ss['Fpos'].max(),2)

    # fix if StMax is very small Fmax needs to

    standard_row['ss'] = ss
    return standard_row

def analyse_control(standard_row, s):
    '''
    take a standard_row and and samples reduced to one protein
    return (all squeeced into the returned standard_row):
        - C1pos and C2pos as a measure of control suitability (should be between 0 and 1 (optimally around 0.5))
        - C1fit and C1fit as a measure of fit of known concentration to estimated conc (0 < 0.5 < 1)
        - the control_df (sc)
    '''
    
    # extract the controls and load the controls into sc
    sc = s.loc[s['Type'].str.match("^C[12]$"),:]
    
    # skip if controls are not included
    if not len(sc.index):
        standard_row['sc'] = None
        return standard_row
    
    C_extract_pattern = r"(?P<Cmin>[0-9]+(?:\.[0-9])?) ?[-–] ?(?P<Cmax>[0-9]+(?:\.[0-9])?)"
    sc = sc.merge(standard_row.loc[['C1', 'C2']].str.extract(C_extract_pattern).astype(float), left_on='Type', right_index=True)
    
    # calculate the conc estimated from params
    sc = compute_conc(sc, standard_row)
    
    control_pos = sc.groupby('Type')['Fpos'].mean().rename({'C1':'C1Pos', 'C2':'C2pos'})
    # compute the fit of the control samples
    sc.loc[:, 'Cfit'] = np.log(sc['conc'] / sc['Cmin']) / (np.log(sc['Cmax'] / sc['Cmin']))
    control_fit = sc.groupby('Type')['Cfit'].mean().rename({'C1':'C1fit', 'C2':'C2fit'})
    
    # add these metrices to the standard_row
    standard_row = pd.concat([standard_row,control_pos, control_fit])
    # load control_df (sc) into standard_row
    standard_row['sc'] = sc
    return standard_row


def compute_external_fit(standard_row, df=pd.DataFrame()):
    '''
    for all samples, goes through the standard rows and applies where applicable
    '''
    
    # extract run and protein for df extraction
    run, protein = standard_row.loc[['Run', 'Protein']]
    cols = [f"conc{run}", f"Fpos{run}"]
    
    df.loc[df['Protein'] == protein, cols] = compute_conc(df.loc[df['Protein'] == protein, :], standard_row, conc_col_suff = run).loc[:, cols]
    return df


def mean_row(row, standard_df=pd.DataFrame(), external_mean_method="arithmetic", minFpos=0, **kwargs):
    '''
    get the mean of all these runs and standards
    '''
    conc = []
    Fpos = [] 
    # cycle through conc and Fpos for all standard_runs
    for run in standard_df['Run'].unique():
        ccol = f"conc{run}"
        fcol = f"Fpos{run}"
        # only append Fpos if not NaN
        if row[fcol] == row[fcol]:
            Fpos.append(row[fcol])
        # only append conc if Fpos > minFpos
        if row[fcol] > minFpos:
            conc.append(row[ccol])
    if (l := len(conc)):
        if external_mean_method=="arithmetic":
            return pd.Series(
                dict(
                    concMean=sum(conc) / len(conc),
                    concStd=np.std(conc),
                    FposMean=sum(Fpos) / len(Fpos)             
                    ))
        elif external_mean_method=="geometric":
            geo_mean = np.power(np.prod(np.abs(conc)), 1/l)
            return pd.Series(
                dict(
                    concMean=round(geo_mean, 3),
                    concStd=np.std(conc),
                    FposMean=sum(Fpos) / len(Fpos)
                )
            )