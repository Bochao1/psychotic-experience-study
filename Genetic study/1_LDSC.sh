
#!/bin/bash
#SBATCH -t 4:00:00
#SBATCH --mem=10G
#SBATCH --mail-type=END
#SBATCH --mail-user=b.lin@umcutrecht.nl
#SBATCH -o log.out
#SBATCH -e errlog.out


cd /hpc/hers_en/blin/UkBiobank/GWAS/LDSC/ldsc

dir=/hpc/hers_en/blin/UkBiobank/PE_MR/data/LDSC/

module load python

#### Munge Data

while IFS= read -r line; do
./munge_sumstats.py --sumstats ${dir}$line".txt.gz" --out   ${dir}$line --merge-alleles w_hm3.noMHC.snplist 
done < /hpc/hers_en/blin/UkBiobank/PE_MR/data/riskfactor/list_small


./munge_sumstats.py \
--sumstats ${dir}PE.txt.gz \
--out   ${dir}PE \
--merge-alleles w_hm3.noMHC.snplist

module load python 
cd /hpc/hers_en/blin/UkBiobank/GWAS/LDSC/ldsc
./munge_sumstats.py \
--sumstats /hpc/hers_en/blin/UkBiobank/GWAS/PE/old/PE_Ad.txt.gz \
--out   ${dir}PE_Ad \
--merge-alleles w_hm3.noMHC.snplist

### LD Score Regression to estimate SNP heritability 

cd /hpc/hers_en/blin/UkBiobank/GWAS/LDSC/ldsc
while IFS= read -r line; do
./ldsc.py --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out ${dir}${line} \
 --h2 ${dir}${line}.sumstats.gz \
done < /hpc/hers_en/blin/UkBiobank/PE_MR/data/riskfactor/list

./ldsc.py --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out ${dir}PE \
 --h2 ${dir}PE.sumstats.gz 

./ldsc.py --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out ${dir}PE_Ad \
 --h2 ${dir}PE_Ad.sumstats.gz 


## LD Score Regression to estimate  Genetic Correlation 


while IFS= read -r line; do
./ldsc.py --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/  --out ${dir}${line}_PE_Ad \
--rg ${dir}${line}.sumstats.gz,${dir}PE_Ad.sumstats.gz \
--no-intercept \
done < /hpc/hers_en/blin/UkBiobank/PE_MR/data/riskfactor/list

cd /hpc/hers_en/blin/UkBiobank/GWAS/LDSC/ldsc
while IFS= read -r line; do
./ldsc.py --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/  --out ${dir}OUT/${line}_PE \
--rg ${dir}PE.sumstats.gz,${dir}${line}.sumstats.gz \
done < /hpc/hers_en/blin/UkBiobank/PE_MR/data/riskfactor/list



