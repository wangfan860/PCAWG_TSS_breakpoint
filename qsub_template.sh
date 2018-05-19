#PBS -N fwang1_XX_CHUNK_XX
#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb
#PBS -S /bin/bash
#PBS -o /group/lyang-lab/fan/pcawg_data/fwang1_XX_CHUNK_XX.stdout.log
#PBS -e /group/lyang-lab/fan/pcawg_data/fwang1_XX_CHUNK_XX.stderr.log

module load gcc/6.2.0 python/3.5.3

python /group/lyang-lab/fan/pcawg_data/PCAWG_TSS_breakpoint/distance_table.py -g XX_CHUNK_XX
