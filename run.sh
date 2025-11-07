#$ -q all.q
#$ -j y                      # Merge stdout and stderr
#$ -cwd                      # Use the current working directory
#$ -S /bin/bash              # Use Bash shell
#$ -terse                    # Output only job ID
#$ -pe smp 8                # Request 24 cores
#$ -N slingshot       # Name the job
#$ -l h_vmem=150G              # Request 8GB RAM
#$ -l h_rt=24:00:00          # Request 24 hour runtime
#$ -m beas                   # Email at the beginning and end of the job
#$ -M june.yoon@nationwidechildrens.org
#$ -o /igm/home/hxy008/Slingshot_Tracy/slingshot.log           # Save combined stdout/stderr to this log file


# Load the R module (adjust for your environment)
module purge
module load R-4.3.3 ## change this to curret R version
Rscript /igm/home/hxy008/Slingshot_Tracy/run_sh.R ## change the R script path of yours
