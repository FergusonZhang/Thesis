conda activate thesis
cd Thesis/
git pull

cd Data_raw/
gzip -d Cgrand_scaffold_1_shapeit4.vcf.gz &
gzip -d Cgrand_scaffold_2_shapeit4.vcf.gz &
gzip -d Cgrand_scaffold_3_shapeit4.vcf.gz &
gzip -d Cgrand_scaffold_4_shapeit4.vcf.gz &
gzip -d Cgrand_scaffold_5_shapeit4.vcf.gz &
gzip -d Cgrand_scaffold_6_shapeit4.vcf.gz &
gzip -d Cgrand_scaffold_7_shapeit4.vcf.gz &
gzip -d Cgrand_scaffold_8_shapeit4.vcf.gz &
cd ..

python code_parse_the_chromosome.py Data_raw/Cgrand_scaffold_1_shapeit4.vcf 100 > Cgrand_scaffold_1_shapeit4.vcf_result.txt &
python code_parse_the_chromosome.py Data_raw/Cgrand_scaffold_2_shapeit4.vcf 100 > Cgrand_scaffold_2_shapeit4.vcf_result.txt &
python code_parse_the_chromosome.py Data_raw/Cgrand_scaffold_3_shapeit4.vcf 100 > Cgrand_scaffold_3_shapeit4.vcf_result.txt &
python code_parse_the_chromosome.py Data_raw/Cgrand_scaffold_4_shapeit4.vcf 100 > Cgrand_scaffold_4_shapeit4.vcf_result.txt &
python code_parse_the_chromosome.py Data_raw/Cgrand_scaffold_5_shapeit4.vcf 100 > Cgrand_scaffold_5_shapeit4.vcf_result.txt &
python code_parse_the_chromosome.py Data_raw/Cgrand_scaffold_6_shapeit4.vcf 100 > Cgrand_scaffold_6_shapeit4.vcf_result.txt &
python code_parse_the_chromosome.py Data_raw/Cgrand_scaffold_7_shapeit4.vcf 100 > Cgrand_scaffold_7_shapeit4.vcf_result.txt &
python code_parse_the_chromosome.py Data_raw/Cgrand_scaffold_8_shapeit4.vcf 100 > Cgrand_scaffold_8_shapeit4.vcf_result.txt &
python code_analyze_the_scores_Tajima.py 50 &
python code_pick_the_genes.py Tajima > Candidates.txt &

vcftools --vcf Data_raw/Cgrand_scaffold_1_shapeit4.vcf --plink --out scaffold_1_plink
vcftools --vcf Data_raw/Cgrand_scaffold_2_shapeit4.vcf --plink --out scaffold_2_plink
vcftools --vcf Data_raw/Cgrand_scaffold_3_shapeit4.vcf --plink --out scaffold_3_plink
vcftools --vcf Data_raw/Cgrand_scaffold_4_shapeit4.vcf --plink --out scaffold_4_plink
vcftools --vcf Data_raw/Cgrand_scaffold_5_shapeit4.vcf --plink --out scaffold_5_plink
vcftools --vcf Data_raw/Cgrand_scaffold_6_shapeit4.vcf --plink --out scaffold_6_plink
vcftools --vcf Data_raw/Cgrand_scaffold_7_shapeit4.vcf --plink --out scaffold_7_plink
vcftools --vcf Data_raw/Cgrand_scaffold_8_shapeit4.vcf --plink --out scaffold_8_plink
python code_change_the_columns.py &

selscan --ihh12 --vcf Data_raw/Cgrand_scaffold_1_shapeit4.vcf --map Data_plink/Cgrand_scaffold_1_shapeit4.vcf_re_plink.map --out Cgrand_scaffold_1_shapeit4.vcf
selscan --ihh12 --vcf Data_raw/Cgrand_scaffold_2_shapeit4.vcf --map Data_plink/Cgrand_scaffold_2_shapeit4.vcf_re_plink.map --out Cgrand_scaffold_2_shapeit4.vcf
selscan --ihh12 --vcf Data_raw/Cgrand_scaffold_3_shapeit4.vcf --map Data_plink/Cgrand_scaffold_3_shapeit4.vcf_re_plink.map --out Cgrand_scaffold_3_shapeit4.vcf
selscan --ihh12 --vcf Data_raw/Cgrand_scaffold_4_shapeit4.vcf --map Data_plink/Cgrand_scaffold_4_shapeit4.vcf_re_plink.map --out Cgrand_scaffold_4_shapeit4.vcf
selscan --ihh12 --vcf Data_raw/Cgrand_scaffold_5_shapeit4.vcf --map Data_plink/Cgrand_scaffold_5_shapeit4.vcf_re_plink.map --out Cgrand_scaffold_5_shapeit4.vcf
selscan --ihh12 --vcf Data_raw/Cgrand_scaffold_6_shapeit4.vcf --map Data_plink/Cgrand_scaffold_6_shapeit4.vcf_re_plink.map --out Cgrand_scaffold_6_shapeit4.vcf
selscan --ihh12 --vcf Data_raw/Cgrand_scaffold_7_shapeit4.vcf --map Data_plink/Cgrand_scaffold_7_shapeit4.vcf_re_plink.map --out Cgrand_scaffold_7_shapeit4.vcf
selscan --ihh12 --vcf Data_raw/Cgrand_scaffold_8_shapeit4.vcf --map Data_plink/Cgrand_scaffold_8_shapeit4.vcf_re_plink.map --out Cgrand_scaffold_8_shapeit4.vcf
python code_analyze_the_scores_iHH12.py 50 &

selscan --ehh scaffold_1:5617863 --vcf Data_raw/Cgrand_scaffold_1_shapeit4.vcf --maf 0.001 --map Data_plink/Cgrand_scaffold_1_shapeit4.vcf_re_plink.map --out Locus_1
selscan --ehh scaffold_1:6734760 --vcf Data_raw/Cgrand_scaffold_1_shapeit4.vcf --maf 0.001 --map Data_plink/Cgrand_scaffold_1_shapeit4.vcf_re_plink.map --out Locus_2
selscan --ehh scaffold_1:18910362 --vcf Data_raw/Cgrand_scaffold_1_shapeit4.vcf --maf 0.001 --map Data_plink/Cgrand_scaffold_1_shapeit4.vcf_re_plink.map --out Locus_3
selscan --ehh scaffold_2:11831023 --vcf Data_raw/Cgrand_scaffold_2_shapeit4.vcf --maf 0.001 --map Data_plink/Cgrand_scaffold_2_shapeit4.vcf_re_plink.map --out Locus_4
selscan --ehh scaffold_3:13306930 --vcf Data_raw/Cgrand_scaffold_3_shapeit4.vcf --maf 0.001 --map Data_plink/Cgrand_scaffold_3_shapeit4.vcf_re_plink.map --out Locus_5
selscan --ehh scaffold_4:1753949 --vcf Data_raw/Cgrand_scaffold_4_shapeit4.vcf --maf 0.001 --map Data_plink/Cgrand_scaffold_4_shapeit4.vcf_re_plink.map --out Locus_6
selscan --ehh scaffold_4:7292520 --vcf Data_raw/Cgrand_scaffold_4_shapeit4.vcf --maf 0.001 --map Data_plink/Cgrand_scaffold_4_shapeit4.vcf_re_plink.map --out Locus_7
selscan --ehh scaffold_5:1877966 --vcf Data_raw/Cgrand_scaffold_5_shapeit4.vcf --maf 0.001 --map Data_plink/Cgrand_scaffold_5_shapeit4.vcf_re_plink.map --out Locus_8
selscan --ehh scaffold_5:2187764 --vcf Data_raw/Cgrand_scaffold_5_shapeit4.vcf --maf 0.001 --map Data_plink/Cgrand_scaffold_5_shapeit4.vcf_re_plink.map --out Locus_9
selscan --ehh scaffold_5:7836339 --vcf Data_raw/Cgrand_scaffold_5_shapeit4.vcf --maf 0.001 --map Data_plink/Cgrand_scaffold_5_shapeit4.vcf_re_plink.map --out Locus_10
selscan --ehh scaffold_5:8801691 --vcf Data_raw/Cgrand_scaffold_5_shapeit4.vcf --maf 0.001 --map Data_plink/Cgrand_scaffold_5_shapeit4.vcf_re_plink.map --out Locus_11
selscan --ehh scaffold_6:1398068 --vcf Data_raw/Cgrand_scaffold_6_shapeit4.vcf --maf 0.001 --map Data_plink/Cgrand_scaffold_6_shapeit4.vcf_re_plink.map --out Locus_12
selscan --ehh scaffold_6:3898014 --vcf Data_raw/Cgrand_scaffold_6_shapeit4.vcf --maf 0.001 --map Data_plink/Cgrand_scaffold_6_shapeit4.vcf_re_plink.map --out Locus_13
selscan --ehh scaffold_6:4742683 --vcf Data_raw/Cgrand_scaffold_6_shapeit4.vcf --maf 0.001 --map Data_plink/Cgrand_scaffold_6_shapeit4.vcf_re_plink.map --out Locus_14
selscan --ehh scaffold_6:8025196 --vcf Data_raw/Cgrand_scaffold_6_shapeit4.vcf --maf 0.001 --map Data_plink/Cgrand_scaffold_6_shapeit4.vcf_re_plink.map --out Locus_15
selscan --ehh scaffold_6:13784834 --vcf Data_raw/Cgrand_scaffold_6_shapeit4.vcf --maf 0.001 --map Data_plink/Cgrand_scaffold_6_shapeit4.vcf_re_plink.map --out Locus_16
selscan --ehh scaffold_6:14094248 --vcf Data_raw/Cgrand_scaffold_6_shapeit4.vcf --maf 0.001 --map Data_plink/Cgrand_scaffold_6_shapeit4.vcf_re_plink.map --out Locus_17
selscan --ehh scaffold_6:14859274 --vcf Data_raw/Cgrand_scaffold_6_shapeit4.vcf --maf 0.001 --map Data_plink/Cgrand_scaffold_6_shapeit4.vcf_re_plink.map --out Locus_18
selscan --ehh scaffold_7:8200301 --vcf Data_raw/Cgrand_scaffold_7_shapeit4.vcf --maf 0.001 --map Data_plink/Cgrand_scaffold_7_shapeit4.vcf_re_plink.map --out Locus_19
selscan --ehh scaffold_7:14582699 --vcf Data_raw/Cgrand_scaffold_7_shapeit4.vcf --maf 0.001 --map Data_plink/Cgrand_scaffold_7_shapeit4.vcf_re_plink.map --out Locus_20
selscan --ehh scaffold_7:16044114 --vcf Data_raw/Cgrand_scaffold_7_shapeit4.vcf --maf 0.001 --map Data_plink/Cgrand_scaffold_7_shapeit4.vcf_re_plink.map --out Locus_21
selscan --ehh scaffold_8:1190253 --vcf Data_raw/Cgrand_scaffold_8_shapeit4.vcf --maf 0.001 --map Data_plink/Cgrand_scaffold_8_shapeit4.vcf_re_plink.map --out Locus_22
selscan --ehh scaffold_8:6178744 --vcf Data_raw/Cgrand_scaffold_8_shapeit4.vcf --maf 0.001 --map Data_plink/Cgrand_scaffold_8_shapeit4.vcf_re_plink.map --out Locus_23
selscan --ehh scaffold_8:7032339 --vcf Data_raw/Cgrand_scaffold_8_shapeit4.vcf --maf 0.001 --map Data_plink/Cgrand_scaffold_8_shapeit4.vcf_re_plink.map --out Locus_24
selscan --ehh scaffold_8:7579458 --vcf Data_raw/Cgrand_scaffold_8_shapeit4.vcf --maf 0.001 --map Data_plink/Cgrand_scaffold_8_shapeit4.vcf_re_plink.map --out Locus_25
selscan --ehh scaffold_8:7706991 --vcf Data_raw/Cgrand_scaffold_8_shapeit4.vcf --maf 0.001 --map Data_plink/Cgrand_scaffold_8_shapeit4.vcf_re_plink.map --out Locus_26
selscan --ehh scaffold_8:8766997 --vcf Data_raw/Cgrand_scaffold_8_shapeit4.vcf --maf 0.001 --map Data_plink/Cgrand_scaffold_8_shapeit4.vcf_re_plink.map --out Locus_27
selscan --ehh scaffold_8:9745220 --vcf Data_raw/Cgrand_scaffold_8_shapeit4.vcf --maf 0.001 --map Data_plink/Cgrand_scaffold_8_shapeit4.vcf_re_plink.map --out Locus_28

cd Data_raw/
gzip Cgrand_scaffold_1_shapeit4.vcf.gz &
gzip Cgrand_scaffold_2_shapeit4.vcf.gz &
gzip Cgrand_scaffold_3_shapeit4.vcf.gz &
gzip Cgrand_scaffold_4_shapeit4.vcf.gz &
gzip Cgrand_scaffold_5_shapeit4.vcf.gz &
gzip Cgrand_scaffold_6_shapeit4.vcf.gz &
gzip Cgrand_scaffold_7_shapeit4.vcf.gz &
gzip Cgrand_scaffold_8_shapeit4.vcf.gz &
cd ..