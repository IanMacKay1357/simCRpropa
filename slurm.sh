#!/bin/bash -l
#SBATCH --job-name=crpropa # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<>   # Where to send mail	
#SBATCH --time=1-00:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --output=/home/iamackay/simCRpropa/output/logs/slurm_%A.log  # Standard output and error log
#SBATCH --mem=1G                    
#SBATCH --array=0-499
#SBATCH --partition=windfall

echo "Config file used:"
echo $1

module load crpropa
export PYTHONPATH=$PYTHONPATH:/home/iamackay/simCRpropa/simCRpropa:/home/iamackay/fermiAnalysis
export LSCRATCH="/home/iamackay/simCRpropa/simCRpropa/output/_scratch"
python /home/iamackay/simCRpropa/simCRpropa/scripts/run_crpropa_em_cascade.py --conf $1 

