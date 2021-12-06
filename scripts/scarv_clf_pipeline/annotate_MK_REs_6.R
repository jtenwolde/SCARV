# annotate expanded humpies with max-atac, then liftOver to hg38 

require(rtracklayer)
options(scipen = 999)

chain = import.chain("/home/jwt44/software/hg19ToHg38.over.chain")

lifted_over = list()
for (cell_type in c("MK", "EB", "B", "aCD4", "rCD4", "MON")) {
	filename = paste0("/rds/project/who1000/rds-who1000-wgs10k/analysis/regulome/humpy/", cell_type, "-humpy400_-1_47-expansion.bed")
	gr = import(filename)
	gr_hg38 = unlist(liftOver(gr, chain))
	lifted_over[cell_type] = gr_hg38
} 

# concatenate and merge all the liftOver humpies
MK_humpies = lifted_over$MK
non_MK_humpies_merged = reduce(Reduce('c', lifted_over[2:length(lifted_over)]))
MK_only_humpies = setdiff(MK_humpies, non_MK_humpies_merged)

df = data.table(chrom=as.vector(seqnames(MK_only_humpies)),
				start=start(MK_only_humpies)-1,
				end=end(MK_only_humpies))
fwrite(df, sep='\t', col.names = FALSE, file='/home/jwt44/MK_humpies_excl_merged_unsorted.bed')
system("sort-bed /home/jwt44/MK_humpies_excl_merged_unsorted.bed > /home/jwt44/MK_humpies_excl_merged_sorted.bed")

cmd = "tabix /rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/tracks/scarv_clf.bed.gz -R /home/jwt44/MK_humpies_excl_merged_sorted.bed"
SCARV_scores_raw = system(cmd, intern = TRUE)


SCARV_scores_split = strsplit(SCARV_scores_raw, split='\t')
SCARV_scores_gr = GRanges(unlist(lapply(SCARV_scores_split, function(x) {x[1]})),
					IRanges(as.numeric(unlist(lapply(SCARV_scores_split, function(x) {x[2]}))),
							as.numeric(unlist(lapply(SCARV_scores_split, function(x) {x[3]})))))
SCARV_scores_gr$SCARV_clf = unlist(lapply(SCARV_scores_split, function(x) {x[4]}))
							

# split and annotate each individual site with SCARV-clf 
df = data.table(chrom=as.vector(seqnames(SCARV_scores_gr)),
				start=start(SCARV_scores_gr),
				end=end(SCARV_scores_gr),
				score=SCARV_scores_gr$SCARV_clf)
fwrite(df, sep='\t', col.names=FALSE, file='/home/jwt44/MK_humpies_excl_merged_and_chopped_and_annotated_unsorted.bed')
system("sort-bed /home/jwt44/MK_humpies_excl_merged_and_chopped_and_annotated_unsorted.bed > /home/jwt44/MK_humpies_excl_merged_and_chopped_and_annotated_sorted.bed")





