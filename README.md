# Steps to run the osmotic pressure simulations for K-PFOS (H1T4 molecules)

1. Generate simulation cells coordinates for the concentration range wanted using the Jupyter Notebook `0.KPFOS.genCmds.ipynb`
2. Carry out the simulations using the LAMMPS software by running the command `bash run.bunch.sh` in the folder `H1T4-PFOS-K-1`. It will create folders for each concentration and run number with a name consisting of N-{#surfactants}-{run#}
3. Merge the results using the command `python 1.1.merge.results.py H1T4-PFOS-K-1`. It will save a csv file: `df.H1T4-PFOS-K-1.saved.csv.gz`
4. Analyze the results and determine the CMC with the Jupyter Notebook `1.2.Post.Analysis.H1T4.PFOS-K.ipynb`