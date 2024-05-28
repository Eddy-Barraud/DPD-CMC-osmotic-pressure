# This script gather all the simulations run and collect data inside a common DataFrame
# It saves the DataFrame to disk with the root folder name passed as first argument

# import libraries
import pandas as pd
import numpy as np
from math import *
import os
import re
import sys
import multiprocessing as mp
import workers
pool = mp.Pool(8) # Parallel processing of simulation data

regex1="^N([0-9.]*)_([0-9]*)$" # Get back the number of surfactants and the run number with the folder name
regex2="^4.ave.csv$" # The production runs save thermo data inside these files

filesToImport=[]
root=str(sys.argv[1])+"/"
for folder in os.listdir(root):
    matches = re.match(regex1, folder)
    if matches:
        prod = False
        NTensio = int(matches[1])
        runNb = int(matches[2])
        path=root+folder
        #print(f"Found run {runNb} with {NTensio} surfactants in directory {path}")
        for file in os.listdir(path):
            matches2 = re.match(regex2, file)
            path2 = path + "/" + file
            if matches2:
                prod=True
                #print(f"found {file} at path {path2}")
                filesToImport+=[[NTensio, runNb, path2]]
        if prod == False:
            print(f"problem with {NTensio} run {runNb}")

dfsList=pool.map(workers.loadCsvToArrWrapper, [(row[2], 0) for row in filesToImport])
df=pd.concat(dfsList, keys=[(row[0], row[1]) for row in filesToImport], names=["NTensio", "runNb", "index"])
df=df.sort_index(level=0)
df.to_csv("df."+str(sys.argv[1])+".saved.csv.gz", compression="gzip")
filesToImport=np.array(filesToImport)
np.save('filesToImport.'+str(sys.argv[1])+'.npy', filesToImport)
pool.close()
exit()