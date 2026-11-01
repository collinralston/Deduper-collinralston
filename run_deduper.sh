#!/bin/bash
#!/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00

#/usr/bin/time -v pixi run samtools sort -O sam -o temp.sam /projects/bgmp/shared/deduper/C1_SE_uniqAlign.sam
/usr/bin/time -v ./Ralston_deduper.py -f temp.sam -o 2026-10-28-deduper_test6.sam -u STL96.txt
