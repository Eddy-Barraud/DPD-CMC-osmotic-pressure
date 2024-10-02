#!/bin/bash

NTensio=$1
runNb=$2
lenTensio=$3
xLen=$4
yLen=$5
zLen=$6
wall1Pos=$7
wall2Pos=$8

fName="N"$NTensio"_"$runNb

mkdir $fName

sbatch <<EOT
#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node 36
#SBATCH -J Lammps
#SBATCH --exclusive

#SBATCH -e $fName/job.`date +%y.%m.%d-%Hh%M`.%A.err
#SBATCH -o $fName/job.`date +%y.%m.%d-%Hh%M`.%A.out
#SBATCH --time=2-0:0:0
#SBATCH -p standard


set -e
# Error handling function
handle_error() {
    echo "An error occurred with NTensio $NTensio run $runNb . Exiting..."
    echo "An error occurred with NTensio $NTensio run $runNb . Exiting..." >> errors.txt
    exit 1
}
# Register the error handling function
trap handle_error ERR

if [ ! -f $fName"/3.eq.done.restart" ]; then
    mpirun -ppn 36 lmp -var seed `bash -c 'echo $((1 + $RANDOM))'` -var runNb $runNb -var NTensio $NTensio -var lenTensio $lenTensio  -var xLen $xLen -v yLen $yLen -v zLen $zLen -v wall1Pos $wall1Pos -v wall2Pos $wall2Pos -log $fName/eq.logs -in 0.eq.lmps -sf hybrid gpu omp -pk gpu 0
    mpirun -ppn 36 lmp -var seed `bash -c 'echo $((1 + $RANDOM))'` -var runNb $runNb -var NTensio $NTensio -var lenTensio $lenTensio  -var xLen $xLen -v yLen $yLen -v zLen $zLen -v wall1Pos $wall1Pos -v wall2Pos $wall2Pos -log $fName/prod.logs -in 0.prod.lmps -sf hybrid gpu omp -pk gpu 0
else
    echo "equilibration already done going directly to production"
    mpirun -ppn 36 lmp -var seed `bash -c 'echo $((1 + $RANDOM))'` -var runNb $runNb -var NTensio $NTensio -var lenTensio $lenTensio  -var xLen $xLen -v yLen $yLen -v zLen $zLen -v wall1Pos $wall1Pos -v wall2Pos $wall2Pos -log $fName/prod.logs -in 0.prod.lmps -sf hybrid gpu omp -pk gpu 0
fi

exit 0
EOT
