from math import *
import numpy as np
import piecewise_regression
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import piecewise_regression
import pandas as pd

# Constants
Na = 6.02214e23 # Avogadro constant
NuH2O=30e-30    # 30 A^3
kB=1.3806e-23 # J/K
R=8.3145 # J/K/mol --> kB * Na

def loadCsvToArrWrapper(args):
   return loadCsvToArr(*args)

def loadCsvToArr(file, skip):
    
    try:
        if skip != 0:
            dfTmp=pd.read_csv(file, sep=' ', header=1, escapechar="#", skiprows= lambda a : a % skip != 0 and a != 1)
        else:
            dfTmp=pd.read_csv(file, sep=' ', header=1, escapechar="#")
    
    except Exception as e:
        print(f"problem with {file} error: {e}")
        return None

    print("file ", file, " imported successfully")
    return dfTmp

def plot_cube(ax, roots, tops, color, axisMin=None, axisMax=None):
    # Define the eight vertices of the cube
    vertices = [
        (roots[0], roots[1], roots[2]),
        (tops[0], roots[1], roots[2]),
        (tops[0], tops[1], roots[2]),
        (roots[0], tops[1], roots[2]),
        (roots[0], roots[1], tops[2]),
        (tops[0], roots[1], tops[2]),
        (tops[0], tops[1], tops[2]),
        (roots[0], tops[1], tops[2]),
    ]

    # Define the six faces of the cube using the vertices
    faces = [
        [vertices[0], vertices[1], vertices[2], vertices[3]],
        [vertices[4], vertices[5], vertices[6], vertices[7]],
        [vertices[0], vertices[1], vertices[5], vertices[4]],
        [vertices[2], vertices[3], vertices[7], vertices[6]],
        [vertices[1], vertices[2], vertices[6], vertices[5]],
        [vertices[0], vertices[3], vertices[7], vertices[4]],
    ]

    # Plot the cube
    ax.add_collection3d(Poly3DCollection(faces, facecolors=color, linewidths=1, edgecolors='r', alpha=.15))
    # Set the axis limits
    if axisMin == None:
        axisMin=np.min(roots)
        axisMax=np.max(tops)
    ax.set_xlim(axisMin, axisMax)
    ax.set_ylim(axisMin, axisMax)
    ax.set_zlim(axisMin, axisMax)

    # Set axis labels
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')


def concentrationsH1TnNa(Ntensio, meanWaterInVc, Vc, Nm, Rho=3):
    # Control volume
    NNA=Ntensio # Same number of Na+ ions
    NbeadsVc=Rho*Vc
    NWB=meanWaterInVc
    NH2O=NWB*3
    NH2O+=NNA*3
    MH2O=18.01528e-03 # kg/mol
    Cm=Ntensio/(NH2O*MH2O) # Avogadro constant simplified
    C=Ntensio/(Na*1e3*Vc*(Rho*Nm*NuH2O))
    phi_vol = Ntensio / meanWaterInVc

    return {"C (mol/L)": C, "C (mol/kg)": Cm, "phi_vol": phi_vol}

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a, ddof=0)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return [m, m-h, m+h, h]

def rollAvg(data, window_size):
    rolls=data.rolling(window_size, center=True,step=window_size,min_periods=window_size,closed="right")

    muBlocs = rolls.mean()
    muBlocs.dropna()
    nbBlocs=len(muBlocs)

    muMean = muBlocs.mean()
    muBlocsStd=muBlocs.std()
    muStd=muBlocsStd/sqrt(nbBlocs-1)
    #print(f"The osmotic pressure is {muMean} \u00B1 {muStd}")
    return [nbBlocs,muMean, muStd]

def statsPosmoH1TnNa(df, Vc=None, rho=None,rc=6.4633E-10,
 NuH2O=30e-30, T=298, Nm=3, verbose=True, Nions = 0, Cions=0):

    VcDeclared = False if Vc == None else True
    df["Phi"] = 0.0
    dfFinal=pd.DataFrame(columns=["NTensio","runNb", "C (mol/L)","Cm (mol/kg)","Pos mean","Pos error", "Pos error #blocs", "Pos error blocLen", "Pos mean (bar)","Pos error (bar)", "Phi mean"]).set_index(["NTensio", "runNb"])
    
    # Check that all simulations have the same length
    first = True
    lenTmp=0
    diffLen = False
    for NT in df.index.unique(0):
        for rn in df.loc[NT].index.unique(0):
            work = df.loc[NT, rn]
            currLen = len(work)
            #print(f'{NT} / {rn} / {currLen}')
            if lenTmp != currLen and first == False:
                print(f"The length of simulation {NT} / {rn} / {currLen} is different from {lenTmp}")
                diffLen = True
            lenTmp = currLen
            first = False
    if diffLen == False:
        print(f"All the simulations have the same length with {lenTmp} frames")

    for NTensio in df.index.unique(0):
        for runNb in df.loc[NTensio].index.unique(0):

            # Osmotic pressure rolling averages and Standard Error of the Mean

            # Wether to calculate the Vc or use the defined one
            if VcDeclared == False:
                if 'v_VcVol' in df.columns:
                    Vc=df.loc[(NTensio, runNb)]["v_VcVol"].iloc[0]
                if 'Vc' in df.columns:
                    Vc=df.loc[(NTensio, runNb)]["Vc"].iloc[0]
                
            if Vc == None:
                    raise Exception("Sorry, no control volume found")
                    
            
            prodRun=df.loc[(NTensio, runNb)].reset_index()["v_osmoticPressure"]
            meanPos=prodRun.mean()
            meanPosBar=meanPos * (kB*T/rc**3) * 1E-05
            resultsList=[]
            for windowSize in range(1,500,1) :
                adding=[windowSize]
                adding.extend(rollAvg(prodRun, windowSize))
                resultsList+=[adding]
            results=pd.DataFrame(resultsList, columns=["windowSize","#blocs","mean","std"])
            
            # Find the convergence of standard deviation when the second derivative change of sign
            plateau = results[20:].loc[results["std"].diff().diff().apply(np.sign).diff().ne(0)].iloc[0]
            #print(f'For mu= {mu} , the mean is: {plateau["mean"]} \u00B1 {plateau["std"]} using {plateau["#blocs"]} blocs')
            
            # Wether to use the meanRho or the defined rho 
            if rho != None :
                meanRho = rho
            else :
                meanRho = df.loc[(NTensio, runNb),"v_sysdensity"].mean()

            meanWaterInVc = df.loc[(NTensio, runNb),"v_waterInVc"].mean()

            CCm = concentrationsH1TnNa(NTensio, meanWaterInVc, meanRho, Vc, Nm)
            
            dfFinal.loc[(NTensio, runNb),"C (mol/L)"] = CCm["C (mol/L)"]
            dfFinal.loc[(NTensio, runNb),"log(C (mol/L))"] = np.log(CCm["C (mol/L)"])
            dfFinal.loc[(NTensio, runNb),"Cm (mol/kg)"] = CCm["C (mol/kg)"]
            dfFinal.loc[(NTensio, runNb),"phi_vol"] = CCm["phi_vol"]
            dfFinal.loc[(NTensio, runNb),"Pos mean"] = meanPos
            dfFinal.loc[(NTensio, runNb),"Pos error"] = plateau["std"]*2
            dfFinal.loc[(NTensio, runNb),"Pos error #blocs"] = plateau["#blocs"]
            dfFinal.loc[(NTensio, runNb),"Pos error blocLen"] = plateau["windowSize"]

            dfFinal.loc[(NTensio, runNb),"Pos mean (bar)"] = meanPosBar
            dfFinal.loc[(NTensio, runNb),"Pos error (bar)"] = plateau["std"]*2 * (kB*T/rc**3) * 1E-05

            # -----------------------------------------------------------------------------------
            #######  Osmotic coefficient rolling averages and Standard Error of the Mean  #######
            # -----------------------------------------------------------------------------------

            idealPos = (NTensio+Nions) * 2 / Vc if rho != None else (NTensio+Nions) * 2 / meanWaterInVc * meanRho

            idealPosBar = (CCm["C (mol/L)"]+Cions) * 2 * R * T * 1E03 * 1E-05
            dfFinal.loc[(NTensio, runNb),"Pos Ideal (bar)"] = idealPosBar
            df_filtered = df.loc[(NTensio, runNb),"v_osmoticPressure"].values
            df_filtered /= idealPos
            df.loc[(NTensio, runNb),"Phi"] = df_filtered

            phiError=rollAvg(df.loc[(NTensio, runNb)].reset_index()["Phi"], plateau["windowSize"].astype(int))
            
            
            # Osmotic coefficient mean
            dfFinal.loc[(NTensio, runNb),"Phi mean"] = meanPos / idealPos
            dfFinal.loc[(NTensio, runNb),"Phi error"] = phiError[2]*2    


    dfFinal=dfFinal.reset_index().set_index(["C (mol/L)", "Cm (mol/kg)", "NTensio"])
    dfFinal.sort_index(level=0, inplace=True)
    if verbose == True : display(dfFinal)

    tmp=[]
    for NTensio in dfFinal.index.get_level_values("NTensio").unique():
        workingDf=dfFinal.xs(NTensio, level="NTensio").reset_index()
        
        meanAllRuns=workingDf["Pos mean (bar)"].mean()
        PhimeanAllRuns=workingDf["Phi mean"].mean()
        nbRuns=len(workingDf["Pos mean (bar)"])

        if nbRuns==1:
            SEMAllRuns=workingDf["Pos error (bar)"].iloc[0]
            PhiSEMAllRuns=workingDf["Phi error"].iloc[0]
        else:
            # mean of variances technic :
            meanVarPos = workingDf["Pos error (bar)"].apply(lambda s: s**2).mean()
            SEMAllRuns = np.sqrt(meanVarPos)
            meanVarPhi = workingDf["Phi error"].apply(lambda s: s**2).mean()
            PhiSEMAllRuns = np.sqrt(meanVarPhi)

            # standard error of the mean technic :
            #sdAllRuns=workingDf["Pos mean (bar)"].std()
            #SEMAllRuns=sdAllRuns/sqrt(nbRuns) * 2
            #PhisdAllRuns=workingDf["Phi mean"].std()
            #PhiSEMAllRuns=PhisdAllRuns/sqrt(nbRuns) * 2

        minCI = meanAllRuns - SEMAllRuns
        maxCI = meanAllRuns + SEMAllRuns


        PhiminCI = PhimeanAllRuns - PhiSEMAllRuns
        PhimaxCI = PhimeanAllRuns + PhiSEMAllRuns


        tmp+=[[workingDf["C (mol/L)"].mean(), workingDf["log(C (mol/L))"].mean(), workingDf["Cm (mol/kg)"].mean(),
         NTensio, meanAllRuns, SEMAllRuns, maxCI, minCI, workingDf["Phi mean"].mean(), PhiSEMAllRuns, PhiminCI, PhimaxCI]]

    dfFinalMerged=pd.DataFrame(tmp, columns=["C (mol/L)", "log(C (mol/L))", "Cm (mol/kg)", "NTensio", "Pos mean (bar)", "Pos error",
     "Pos +", "Pos -", "Phi mean", "Phi error", "Phi -", "Phi +"])
    if verbose == True : display(dfFinalMerged)

    return dfFinal, dfFinalMerged

# Return the power of ten of a string in scientific notation
# "3.42e3" --> 3
def get_power_of_ten(scientific_notation):
    # Split the scientific notation into coefficient and exponent parts
    coefficient, exponent = scientific_notation.split('e')
    return int(exponent)

# Returns the value with 00 below the error position
def round_on_error(val, err):
    val_power_ten = get_power_of_ten(f'{val:.10e}')
    err_power_ten = get_power_of_ten(f'{err:.10e}')
    significant_numbers = abs(err_power_ten - val_power_ten)
    significant_numbers_err = min(3,significant_numbers)

    formatted_result = "{f_val:.{SN}f}e{f_exp:+03} ± {f_err:.{SNerr}f}e{f_exp_err:+03}".format(
    f_val=val*10**(-val_power_ten), SN=significant_numbers, SNerr= significant_numbers_err,
    f_exp = val_power_ten, f_exp_err = err_power_ten, f_err=err*10**(-err_power_ten))
    
    formatted_value = "{f_val:.{SN}f}e{f_exp:+03}".format(
    f_val=val*10**(-val_power_ten), SN=significant_numbers, f_exp = val_power_ten)
    
    formatted_error = "{f_err:.{SNerr}f}e{f_exp_err:+03}".format(
    SNerr= significant_numbers_err, f_exp_err = err_power_ten, f_err=err*10**(-err_power_ten))

    return formatted_result, formatted_value, formatted_error


def plotCMC(dfFinal, dfFinalMerged, xC="C (mol/L)", yC="Pos mean (bar)", yCerr="Pos error",
 xLabel=None , yLabel=None, VH=True, bp=6, bpN=1, down=0.0, up=100.0, T=298, CNACL=0.0, ax=None, exp=None):
    dfFinal.sort_index(inplace=True)
    wanted=dfFinal.loc[down:up].reset_index()
    x=wanted[xC].astype(float).values
    y=wanted[yC].astype(float).values
    
    dfFinalMerged = dfFinalMerged.set_index(xC)
    dfFinalMerged.sort_index(inplace=True)
    wanted2=dfFinalMerged.loc[down:up].reset_index()
    xM=wanted2[xC].values
    yM=wanted2[yC].values
    yMerr=wanted2[yCerr].values if yCerr in wanted2.columns else None


    # Fit the data
    ms = piecewise_regression.ModelSelection(x, y, max_breakpoints=bp)
    nbb=pd.DataFrame(data=ms.model_summaries)
    bestbb=int(nbb.iloc[nbb["bic"].idxmin()]["n_breakpoints"])
    print(f"The best number accuracy is obtained using {bestbb} breakpoint(s)")
    pw_fit = piecewise_regression.Fit(x, y, n_breakpoints=bestbb if bestbb > 0 else 1, verbose=False)
    pw_fit.summary()

    # Plot the data, fit, breakpoints and confidence intervals
    plt.style.use('default')
    if ax == None:
        fig, ax = plt.subplots(nrows=1, ncols=1)

    ax.errorbar(xM,yM,yMerr, fmt='.', color="blue", label="DPD")

    # Ideal Van't Hoff osmotic pressure
    # Don't go higher than the maximum measured osmotic pressure
    if VH == True :
        posIdealX=np.linspace(x.min(),abs(x.max()-x.min())+x.min(),100)
        posIdealY= 2 * R*T * (posIdealX + CNACL) * 1E03 * 1E-05
        posIdeal = np.column_stack((posIdealX, posIdealY))
        max_y=np.max(y)
        sortedPosIdeal = posIdeal[posIdeal[:, 1] <= max_y*1.1]
        ax.plot(sortedPosIdeal[:,0], sortedPosIdeal[:,1], "--", color="purple", label="Van't Hoff")

    pw_fit.plot_fit(color="red", linewidth=0.5, label="Fit")
    ax.axvline(pw_fit.best_muggeo.best_fit.next_breakpoints[bpN-1],color="green", label="Breakpoint")
    cibb=pw_fit.best_muggeo.best_fit.estimates[f"breakpoint{bpN}"]["confidence_interval"]
    ax.axvspan(cibb[0], cibb[1], color="green", alpha=0.1)

    CMC = pw_fit.get_results()["estimates"][f"breakpoint{bpN}"]["estimate"]
    CMCErr = pw_fit.get_results()["estimates"][f"breakpoint{bpN}"]["se"]*2
    formatted_CMC = round_on_error(CMC, CMCErr)
    display(f"The CMC found is {formatted_CMC[0]} mol/L")
    arrowY=abs(ax.get_ylim()[1]-ax.get_ylim()[0])*0.5+ax.get_ylim()[0]
    arrowX=abs(ax.get_xlim()[1]-ax.get_xlim()[0])*0.6+ax.get_xlim()[0]

    plt.annotate(f'CMC={formatted_CMC[1]} \n\u00B1 {formatted_CMC[2]} mol/L', xy=(CMC,arrowY), xytext=(arrowX,arrowY),
                arrowprops=dict(facecolor='black', arrowstyle="]->", alpha=1), ha='left')

    # Experimental data:
    if exp != None and isinstance(exp, list) and len(exp) == 2:
        ax.axvspan(*exp, color="red", alpha=0.3, label="Exp")
        
    xLabel = xC if xLabel == None else xLabel
    yLabel = yC if yLabel == None else yLabel
    
    ax.set_xlabel(xLabel)
    ax.set_ylabel(yLabel)
    ax.legend()

    return CMC, CMCErr
