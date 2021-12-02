#!/bin/bash
#SBATCH -t 10:00:00
#SBATCH --mem=100G
#SBATCH -o log.out
#SBATCH -e errlog.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=b.lin@umcutrecht.nl

dir=/hpc/hers_en/blin/UkBiobank/PE_MR/data/TwosampleGSMR
gcta=/hpc/hers_en/blin/program/gcta_1.93.0beta/gcta64

$gcta \
--mbfile $dir/ref.file.list \
--gsmr-file  $dir/f2040 $dir/outcome_PE.list \
--gsmr-snp-min 1 \
--clump-r2 0.01 \
--effect-plot \
--gsmr-direction 0 \
--out $dir/36_PE/2040_PE \
--gwas-thresh 1e-6 \
--thread-num 10

$gcta \
--mbfile $dir/ref.file.list \
--gsmr-file  $dir/f2040 $dir/outcome_SCZ.list \
--gsmr-snp-min 1 \
--clump-r2 0.01 \
--effect-plot \
--gsmr-direction 0 \
--out $dir/36_PE/2040_SCZ \
--gwas-thresh 1e-6 \
--thread-num 10


