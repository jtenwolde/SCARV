paste <(zcat /rds/project/who1000-1/rds-who1000-cbrc/user/dg333/joep-explained-cases-vars3.txt.gz) <(zcat /rds/project/who1000-1/rds-who1000-cbrc/user/dg333/joep-ird2.txt.gz) | \
    awk '{print "chr"$1,$2-1,$2,$3","$4, $5","$6","$7","$8","$9","$10","$11","$12","$17","$18","$19","$20}' OFS='\t' - | \
    sort-bed - | /home/jwt44/software/liftOver stdin /home/jwt44/software/hg19ToHg38.over.chain stdout /home/jwt44/liftOver_unmapped | sort-bed - | \
    bedmap --echo --skip-unmapped - <(zcat /rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/opr/minOPR_gt0.99_hg38_liftedOver_sorted.bed.gz) | \
    gzip > /rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/humpy_analysis/WGS10K_rare_alleles/by_readlength/joep-explained-cases-vars3_hg38_minOPR_fltrd.bed.gz


# remove exons 
zcat /rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/humpy_analysis/WGS10K_rare_alleles/by_readlength/joep-explained-cases-vars3_hg38_minOPR_fltrd.bed.gz | \
    bedops --difference - /rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/gencode_v27/gencode_exons_hg38_v27.bed | \
    bedmap --echo --skip-unmapped <(zcat /rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/humpy_analysis/WGS10K_rare_alleles/by_readlength/joep-explained-cases-vars3_hg38_minOPR_fltrd.bed.gz) - | \
    gzip > /rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/humpy_analysis/WGS10K_rare_alleles/by_readlength/joep-explained-cases-vars3_hg38_minOPR_fltrd_non_exonic.bed.gz


# map to humpies
zcat /rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/humpy_analysis/WGS10K_rare_alleles/by_readlength/joep-explained-cases-vars3_hg38_minOPR_fltrd_non_exonic.bed.gz | \
    bedmap --echo --echo-map-id --skip-unmapped - /home/jwt44/MK_humpies_excl_merged_and_chopped_and_annotated_sorted.bed | \
    sed 's/|/\t/g' - | gzip > /rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/humpy_analysis/WGS10K_rare_alleles/by_readlength/MK_humpies_excl_joep-explained-cases-vars3_hg38_minOPR_fltrd_non_exonic_annot_with_SCARV_clf.bed.gz
