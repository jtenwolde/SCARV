# script to set up training and test set for pathogenicity classifier

cd /rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/

mkdir classifier
cd classifier

mkdir variants
cd variants

# download pathogenic variants as provided by ncER (hg19)
mkdir ncER_SNVs
wget https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6868241/bin/41467_2019_13212_MOESM3_ESM.txt

mv 41467_2019_13212_MOESM3_ESM.txt ncER_SNVs

awk 'NR!=1 && $5=="training" {print $1,$2-1,$2}' OFS='\t' ncER_SNVs/41467_2019_13212_MOESM3_ESM.txt | \
    /home/jwt44/software/liftOver stdin /home/jwt44/software/hg19ToHg38.over.chain ncER_SNVs/ncER_training_SNVs_hg38.bed ncER_SNVs/ncER_training_SNVs_liftOver_unmapped.txt
awk 'NR!=1 && $5=="testing" {print $1,$2-1,$2}' OFS='\t' 41467_2019_13212_MOESM3_ESM.txt | \
    /home/jwt44/software/liftOver stdin /home/jwt44/software/hg19ToHg38.over.chain ncER_SNVs/ncER_testing_SNVs_hg38.bed ncER_SNVs/ncER_testing_SNVs_liftOver_unmapped.txt


# gnomAD variants with AF > 1% across populations
awk '$6/$7>0.01' ../../scarv_pipeline_gnomad_hg38/nfe/variants/pass_snvs.bed | \
    bedmap --echo --echo-map --skip-unmapped - <(awk '$6/$7>0.01' ../../scarv_pipeline_gnomad_hg38/afr/variants/pass_snvs.bed) | \
    awk '$4==$10 && $5==$11 {print $1,$2,$3,$4,$5}' OFS='\t' - | \
    bedmap --echo --skip-unmapped - ../../scarv_pipeline_gnomad_hg38/nfe/quality_filtering/covered_loci.bed > \
    gnomAD_hg38_covered_SNVs_AFgt0p01_AFRandNFE.bed


# regulatory elements
cd ..
mkdir regulatory_build
cd regulatory_build

wget ftp://ftp.ensembl.org/pub/release-99/regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.gff.gz

zcat homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.gff.gz | awk 'length($1)<3 {print "chr"$1,$4-1,$5,$3}' OFS='\t' - | \
    sort-bed - > homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329_simplified.bed

awk '$4=="CTCF_binding_site"' homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329_simplified.bed > homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329_CTCF_simplified.bed
awk '$4=="TF_binding_site"' homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329_simplified.bed > homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329_TFBS_simplified.bed
awk '$4=="promoter"' homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329_simplified.bed > homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329_Promoter_simplified.bed
awk '$4=="enhancer"' homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329_simplified.bed > homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329_Enhancer_simplified.bed
awk '$4=="promoter_flanking_region"' homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329_simplified.bed > homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329_Promoter_Flanking_simplified.bed
awk '$4=="open_chromatin_region"' homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329_simplified.bed > homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329_Open_Chromatin_Region_simplified.bed


# gencode annotations
cd ..
mkdir gencode_v27
cd gencode_v27

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz


zcat gencode.v27.annotation.gtf.gz  | grep -v '#' | \
    awk '$3=="gene" {print $1,$4-1,$5}' OFS='\t' - | sort-bed - | bedops -m - > gencode_genes_hg38_v27.bed

zcat gencode.v27.annotation.gtf.gz  | grep -v '#' | \
    awk '$3=="exon" {print $1,$4-1,$5}' OFS='\t' - | sort-bed - > gencode_exons_hg38_v27.bed

awk '{printf "%s\t%d\t%d\n%s\t%d\t%d\n", $1,$2-1,$2,$1,$3,$3+1}' gencode_exons_hg38_v27.bed | \
    bedops --intersect - gencode_genes_hg38_v27.bed > gencode_splice_sites_hg38_v27.bed

zcat gencode.v27.annotation.gtf.gz | grep 'lincRNA\|ncRNA\|snoRNA' - | awk '{print $1,$4-1,$5}' OFS='\t' - | \
    sort-bed - | bedops -m - | awk '{print $0,"Exon-ncRNA"}' OFS='\t' - > gencode_v27_ncRNAs_simplified.bed

zcat gencode.v27.annotation.gtf.gz | awk '$3=="UTR" {print $1,$4-1,$5}' OFS='\t' - | sort-bed - | bedops -m - | \
    awk '{print $0, "UTR"}' OFS='\t' - > gencode_v27_UTR_simplified.bed

zcat gencode.v27.annotation.gtf.gz | awk '$3=="CDS" {print $1,$4-1,$5}' OFS='\t' - |  sort-bed - | bedops -m - | \
    awk '{print $0, "CDS"}' OFS='\t' - > gencode_v27_CDS_simplified.bed

bedops --difference gencode_v27_introns_simplified.bed <(bedops --range 10 --everything gencode_splice_sites_hg38_v27.bed) | \
    bedops -m - | awk '{print $0, "Intron-Distal"}' OFS='\t' - > gencode_v27_introns_distal_simplified.bed

bedops --intersect gencode_v27_introns_simplified.bed <(bedops --range 10 --everything gencode_splice_sites_hg38_v27.bed) | \
    bedops -m - | awk '{print $0, "Intron-Cis"}' OFS='\t' - > gencode_v27_introns_cis_simplified.bed





