#!/bin/bash
#PBS -N fwang1_cnvplot
#PBS -l nodes=1:ppn=1
#PBS -l mem=10gb
#PBS -l walltime=09:59:59
#PBS -o /group/lyang-lab/fan/pcawg_data/log/fwang1_cnvplot.stdout.log
#PBS -e /group/lyang-lab/fan/pcawg_data/log/fwang1_cnvplot.stderr.log

module load gcc/6.2.0
module load R/3.5.0

/apps/software/gcc-6.2.0/R/3.5.0/bin/Rscript /group/lyang-lab/fan/pcawg_data/PCAWG_TSS_breakpoint/cnvplot.R 1 2 3
