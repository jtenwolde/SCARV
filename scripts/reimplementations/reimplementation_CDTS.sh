# reimplementation CDTS

cd /rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/
mkdir reimplementing_CDTS_gnomad_hg38
cd reimplementing_CDTS_gnomad_hg38

# create bed file called: "hg38_chrXPAR.bed" with following lines,
#chrX    10001   2781479
#chrX    155701383   156030895

# get heptamers
mkdir heptamers

for chr in {1..22}
do
awk -v chrom="chr"${chr} '$1==chrom' /rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/quality_filtering/reliable_sites.bed |\
    bedops --chop 1 - | bedops --range 3 --everything - | bedtools getfasta -tab -fi /rds/project/who1000-1/rds-who1000-cbrc/ref/UCSC/hg38/hg38.fa -bed stdin -fo stdout |\
    awk -F[":-"] '{print $1,$2,toupper($3)}' OFS='\t' - | bedops --range -3 --everything - > heptamers/reliable_heptamers_chr${chr}.bed &
done 

chr=XPAR
awk '$1=="chrX"' /rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/quality_filtering/reliable_sites.bed |\
    bedops --intersect - hg38_chrXPAR.bed | bedops --chop 1 - | bedops --range 3 --everything - | bedtools getfasta -tab -fi /rds/project/who1000-1/rds-who1000-cbrc/ref/UCSC/hg38/hg38.fa -bed stdin -fo stdout |\
    awk -F[":-"] '{print $1,$2,toupper($3)}' OFS='\t' - | bedops --range -3 --everything - > heptamers/reliable_heptamers_chr${chr}.bed &

chr=XnonPAR
awk '$1=="chrX"' /rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/quality_filtering/reliable_sites.bed |\
    bedops --difference - hg38_chrXPAR.bed | bedops --chop 1 - | bedops --range 3 --everything - | bedtools getfasta -tab -fi /rds/project/who1000-1/rds-who1000-cbrc/ref/UCSC/hg38/hg38.fa -bed stdin -fo stdout |\
    awk -F[":-"] '{print $1,$2,toupper($3)}' OFS='\t' - | bedops --range -3 --everything - > heptamers/reliable_heptamers_chr${chr}.bed &

wait


# get snvs by site
mkdir annot_with_snvs

for chr in {1..22}
do
    awk -v chrom="chr"${chr} '$1==chrom && $6/$7 > 0.0001 {print $1,$2,$3}' OFS='\t' /rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/variants/pass_snvs.bed |\
        uniq | bedmap --echo --indicator heptamers/reliable_heptamers_chr${chr}.bed - | sed 's/|/\t/g' - > annot_with_snvs/reliable_heptamers_chr${chr}_annot_with_SNVs_present.bed &
done 

chr=XPAR
awk '$6/$7 > 0.0001 {print $1,$2,$3}' OFS='\t' /rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/variants/pass_snvs.bed |\
    uniq | bedmap --echo --indicator heptamers/reliable_heptamers_chr${chr}.bed - | sed 's/|/\t/g' - > annot_with_snvs/reliable_heptamers_chr${chr}_annot_with_SNVs_present.bed &

chr=XnonPAR
awk '$6/$7 > 0.000075 {print $1,$2,$3}' OFS='\t' /rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/variants/pass_snvs.bed |\
    uniq | bedmap --echo --indicator heptamers/reliable_heptamers_chr${chr}.bed - | sed 's/|/\t/g' - > annot_with_snvs/reliable_heptamers_chr${chr}_annot_with_SNVs_present.bed &

wait


# get proportion snv per k-mer
mkdir counts

for chr in {1..22}
do
    awk '{arr[$4]+=1; if ($5==1) {arr2[$4]+=1;} else {arr2[$4]+=0;}} END {for (key in arr) {print key,arr[key],arr2[key];}}' annot_with_snvs/reliable_heptamers_chr${chr}_annot_with_SNVs_present.bed > counts/heptamers_counts_and_hits_chr${chr}.txt &
done 

chr=XPAR
awk '{arr[$4]+=1; if ($5==1) {arr2[$4]+=1;} else {arr2[$4]+=0;}} END {for (key in arr) {print key,arr[key],arr2[key];}}' annot_with_snvs/reliable_heptamers_chr${chr}_annot_with_SNVs_present.bed > counts/heptamers_counts_and_hits_chr${chr}.txt &

chr=XnonPAR
awk '{arr[$4]+=1; if ($5==1) {arr2[$4]+=1;} else {arr2[$4]+=0;}} END {for (key in arr) {print key,arr[key],arr2[key];}}' annot_with_snvs/reliable_heptamers_chr${chr}_annot_with_SNVs_present.bed > counts/heptamers_counts_and_hits_chr${chr}.txt &

wait


# combine counts
cat counts/heptamers_counts_and_hits_chr* | awk '{print $1}' | sort | uniq | grep -v N - > counts/unique_hepts.txt
join -1 1 -2 1 counts/unique_hepts.txt <(sort counts/heptamers_counts_and_hits_chr1.txt) -a1 |\
    awk '{if (NF==1) {print $0,0,0;} else {print $0;}}' OFS='\t' - > counts/counts_mgd.txt

for chr in {2..22}
do
    join -1 1 -2 1 <(sort counts/counts_mgd.txt) <(sort counts/heptamers_counts_and_hits_chr${chr}.txt) -a1 |\
        awk '{if (NF==3) {print $0,0,0;} else {print $0;}}' OFS='\t' > counts/counts_mgd_tmp.txt
    mv counts/counts_mgd_tmp.txt counts/counts_mgd.txt 
done 

chr=XPAR
join -1 1 -2 1 <(sort counts/counts_mgd.txt) <(sort counts/heptamers_counts_and_hits_chr${chr}.txt) -a1 |\
        awk '{if (NF==3) {print $0,0,0;} else {print $0;}}' OFS='\t' > counts/counts_mgd_tmp.txt
mv counts/counts_mgd_tmp.txt counts/counts_mgd.txt 

awk '{N=0; k=0; for (i=2; i<=NF; i++) {if (i%2==0) {N+=$i;} else {k+=$i;};} print $1,k/N}' OFS='\t' counts/counts_mgd.txt > rate_by_hepta_autosomes_and_chrXPAR.txt


chr=XnonPAR
sort -k1 counts/heptamers_counts_and_hits_chr${chr}.txt | awk '{print $1,$3/$2}' OFS='\t' - > rate_by_hepta_chrXnonPAR.txt


# map k-mers to proportions
mkdir annot_with_rate

for chr in {1..22}
do
    awk 'NR==FNR{maps[$1]=$2} NR!=FNR{print $0, maps[$4]}' OFS='\t' rate_by_hepta_autosomes_and_chrXPAR.txt annot_with_snvs/reliable_heptamers_chr${chr}_annot_with_SNVs_present.bed |\
    grep -v N - > annot_with_rate/reliable_heptamers_chr${chr}_annot_incl_rate.bed &
done  

chr=XPAR
awk 'NR==FNR{maps[$1]=$2} NR!=FNR{print $0, maps[$4]}' OFS='\t' rate_by_hepta_autosomes_and_chrXPAR.txt annot_with_snvs/reliable_heptamers_chr${chr}_annot_with_SNVs_present.bed |\
    grep -v N -  > annot_with_rate/reliable_heptamers_chr${chr}_annot_incl_rate.bed &

chr=XnonPAR
awk 'NR==FNR{maps[$1]=$2} NR!=FNR{print $0, maps[$4]}' OFS='\t' rate_by_hepta_chrXnonPAR.txt annot_with_snvs/reliable_heptamers_chr${chr}_annot_with_SNVs_present.bed |\
    grep -v N -  > annot_with_rate/reliable_heptamers_chr${chr}_annot_incl_rate.bed &

wait


# sliding window
mkdir CDTS

window=550
sliding=10

for chr in {1..22}
do
    grep -w "chr"${chr} ../hg38.genome | \
        awk -v sld=$sliding -v wd=$window 'BEGIN{OFS=FS="\t"}{for(i=0; i<=($2); i=i+sld) {stt=i;end=i+wd; print $1, stt,end}}' | \
        bedmap --echo --bases --sum - <(cut -f1-3,5-6 annot_with_rate/reliable_heptamers_chr${chr}_annot_incl_rate.bed) | \
        sed 's/|/\t/g' | awk -v wd=$window 'BEGIN{FS=OFS="\t"}{if ($4 >= wd*0.9) print $0}' | \
        bedmap --echo --count - <(awk '$5==1' annot_with_rate/reliable_heptamers_chr${chr}_annot_incl_rate.bed) | sed 's/|/\t/g' | \
        awk 'BEGIN{OFS=FS="\t"}{print $1, $2+270, $3-270, $6, $5, $6-$5}' > CDTS/CDTS_diff_chr${chr}.bed & # diff stands for difference
done

chr=XPAR
grep -w "chrX" ../hg38.genome | \
    awk -v sld=$sliding -v wd=$window 'BEGIN{OFS=FS="\t"}{for(i=0; i<=($2); i=i+sld) {stt=i;end=i+wd; print $1, stt,end}}' | \
    bedmap --echo --bases --sum - <(cut -f1-3,5-6 annot_with_rate/reliable_heptamers_chr${chr}_annot_incl_rate.bed) | \
    sed 's/|/\t/g' | awk -v wd=$window 'BEGIN{FS=OFS="\t"}{if ($4 >= wd*0.9) print $0}' | \
    bedmap --echo --count - <(awk '$5==1' annot_with_rate/reliable_heptamers_chr${chr}_annot_incl_rate.bed) | sed 's/|/\t/g' | \
    awk 'BEGIN{OFS=FS="\t"}{print $1, $2+270, $3-270, $6, $5, $6-$5}' > CDTS/CDTS_diff_chr${chr}.bed & 

chr=XnonPAR
grep -w "chrX" ../hg38.genome | \
    awk -v sld=$sliding -v wd=$window 'BEGIN{OFS=FS="\t"}{for(i=0; i<=($2); i=i+sld) {stt=i;end=i+wd; print $1, stt,end}}' | \
    bedmap --echo --bases --sum - <(cut -f1-3,5-6 annot_with_rate/reliable_heptamers_chr${chr}_annot_incl_rate.bed) | \
    sed 's/|/\t/g' | awk -v wd=$window 'BEGIN{FS=OFS="\t"}{if ($4 >= wd*0.9) print $0}' | \
    bedmap --echo --count - <(awk '$5==1' annot_with_rate/reliable_heptamers_chr${chr}_annot_incl_rate.bed) | sed 's/|/\t/g' | \
    awk 'BEGIN{OFS=FS="\t"}{print $1, $2+270, $3-270, $6, $5, $6-$5}' > CDTS/CDTS_diff_chr${chr}.bed & 

wait



wc -l CDTS/* | tail -1 | awk '{print $1}' -
nb=2628352 #== 262835190/100 rounded

cat CDTS/CDTS_diff_chr*.bed | sort -k6,6g -k4,4n --temporary-directory=./ | \
split --lines=${nb} --numeric-suffixes=1 --suffix-length=3 - CDTS/CDTS_diff_perc.bed.


touch CDTS_diff_perc.bed
i=1
FILES=`ls CDTS/CDTS_diff_perc.bed.*`
for f in $FILES
do
    echo "doing " $f
    cat $f | awk -v ivar=${i} -v nbLine=$nb 'BEGIN{OFS=FS="\t"}{printf "%s\t%i\t%i\t%s\t%s\t%s\t%.6f\n", \
        $1, $2, $3, $4, $5, $6, ivar-1+(NR/nbLine)}' >> CDTS/CDTS_diff_perc.bed
    i=$((${i}+1))
done

sort-bed CDTS_diff_perc.bed | gzip - > CDTS_diff_perc_coordsorted_gnomAD_N32299_hg38.bed.gz
rm CDTS_diff_perc.bed






