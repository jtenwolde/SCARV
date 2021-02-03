# reimplementation Orion
source activate scars_env # required for ipython

cd /rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/
mkdir reimplementing_Orion_gnomad_hg38
cd reimplementing_Orion_gnomad_hg38

# get trinic
mkdir trinucs

for chr in {1..22}
do
awk -v chrom="chr"${chr} '$1==chrom' /rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/quality_filtering/reliable_sites.bed |\
    bedops --chop 1 - | bedops --range 1 --everything - | bedtools getfasta -tab -fi /rds/project/who1000-1/rds-who1000-cbrc/ref/UCSC/hg38/hg38.fa -bed stdin -fo stdout |\
    awk -F[":-"] '{print $1,$2,toupper($3)}' OFS='\t' - > trinucs/reliable_trinucs_chr${chr}.bed &
done 

wait


# sum mut rates by trinuc
awk 'NR!=1 {mut_rate[$1]+=$3;} END {for (trinuc in mut_rate) {print trinuc,mut_rate[trinuc];}}' OFS='\t' mut_rate_table.txt | sort -k1 > mut_rate_table_collapsed.txt

# map k-mers to mut_rates
mkdir annot_with_rate

for chr in {1..22}
do
    awk 'NR==FNR{maps[$1]=$2} NR!=FNR{print $0, maps[$4]}' OFS='\t' mut_rate_table_collapsed.txt trinucs/reliable_trinucs_chr${chr}.bed |\
        grep -v N - > annot_with_rate/reliable_trinucs_chr${chr}_annot_incl_rate.bed &
done  

wait


# sliding window
mkdir Orion

window=501
sliding=1
pop_size=32299

for chr in {1..22}
do
    grep -w "chr"${chr} ../hg38.genome | \
        awk -v sld=$sliding -v wd=$window 'BEGIN{OFS=FS="\t"}{for(i=0; i<=($2); i=i+sld) {stt=i;end=i+wd; print $1, stt,end}}' | \
        bedmap --echo --bases --sum --sci - <(cut -f1-3,4-5 annot_with_rate/reliable_trinucs_chr${chr}_annot_incl_rate.bed | bedops --range -1 --everything  -) | \
        sed 's/|/\t/g' | awk -v wd=$window 'BEGIN{FS=OFS="\t"}{if ($4 >= wd*0.9) print $0}' | \
        ipython calculate_Orion.py $pop_size "chr"$chr /rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/variants/pass_snvs_ac_by_chr/ac_by_pos_chr${chr}.txt > Orion/Orion_chr${chr}.bed & 
done

# map to percentiles with python script
ipython map_score_to_percentiles.py 



