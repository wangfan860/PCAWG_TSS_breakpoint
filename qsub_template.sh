#!/bin/bash
#PBS -N fwang1_XX_CHUNK_XX
#PBS -l nodes=1:ppn=1
#PBS -l mem=4gb
#PBS -l walltime=05:59:59
#PBS -o /group/lyang-lab/fan/pcawg_data/fwang1_XX_CHUNK_XX.stdout.log
#PBS -e /group/lyang-lab/fan/pcawg_data/fwang1_XX_CHUNK_XX.stderr.log

/apps/software/gcc-6.2.0/python/3.6.0/bin/python3 /group/lyang-lab/fan/pcawg_data/PCAWG_TSS_breakpoint/distance_table.py -g XX_CHUNK_XX
