require(data.table)
require(gsEasy)
require(ggplot2)
require(ggpubr)
require(survcomp)


# read data frame
f_in = "zcat /rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/humpy_analysis/WGS10K_rare_alleles/by_readlength/MK_humpies_excl_joep-explained-cases-vars3_hg38_minOPR_fltrd_non_exonic_annot_with_SCARV_clf.bed.gz"
df = fread(cmd = f_in, sep='\t', col.names=c("chrom", "start", "end", "ref_alt", "counts", "SCARV_clf"))


# extract reference and alternate allele 
ref_alt_split = t(df[, strsplit(ref_alt, ",")])
df[, ref := unname(ref_alt_split)[,1]]
df[, alt := unname(ref_alt_split)[,2]]
df[, ref_alt := NULL]


# split comma separated column of allele counts into multiple columns
counts_split = t(df[, strsplit(counts, ",")])

df[, AC_explained_125 := as.numeric(unname(counts_split)[,1])]
df[, AC_controls_and_unexp_125 := as.numeric(unname(counts_split)[,2])]

df[, AC_platelet_cases_125 := as.numeric(unname(counts_split)[,3])]
df[, AC_platelet_controls_125 := as.numeric(unname(counts_split)[,4])]

df[, AC_explained_150 := as.numeric(unname(counts_split)[,5])]
df[, AC_controls_and_unexp_150 := as.numeric(unname(counts_split)[,6])]

df[, AC_platelet_cases_150 := as.numeric(unname(counts_split)[,7])]
df[, AC_platelet_controls_150 := as.numeric(unname(counts_split)[,8])]

df[, AC_ird_cases_125 := as.numeric(unname(counts_split)[,9])]
df[, AC_ird_controls_125 := as.numeric(unname(counts_split)[,10])]

df[, AC_ird_cases_150 := as.numeric(unname(counts_split)[,11])]
df[, AC_ird_controls_150 := as.numeric(unname(counts_split)[,12])]

df[, counts := NULL]


# expand SCARV-clf percentile scores according to allele count in each group
explained_125_scores = df[, rep(SCARV_clf, AC_explained_125)]
controls_and_unexp_125_scores = df[, rep(SCARV_clf, AC_controls_and_unexp_125)]

explained_150_scores = df[, rep(SCARV_clf, AC_explained_150)]
controls_and_unexp_150_scores = df[, rep(SCARV_clf, AC_controls_and_unexp_150)]

platelet_cases_125_scores = df[, rep(SCARV_clf, AC_platelet_cases_125)]
platelet_controls_125_scores = df[, rep(SCARV_clf, AC_platelet_controls_125)]

platelet_cases_150_scores = df[, rep(SCARV_clf, AC_platelet_cases_150)]
platelet_controls_150_scores = df[, rep(SCARV_clf, AC_platelet_controls_150)]

ird_cases_125_scores = df[, rep(SCARV_clf, AC_ird_cases_125)]
ird_controls_125_scores = df[, rep(SCARV_clf, AC_ird_controls_125)]

ird_cases_150_scores = df[, rep(SCARV_clf, AC_ird_cases_150)]
ird_controls_150_scores = df[, rep(SCARV_clf, AC_ird_controls_150)]


# sum the number of alleles in each category, used later for meta analysis
weight_platelet_125 = length(platelet_cases_125_scores) + length(platelet_controls_125_scores)
weight_platelet_150 = length(platelet_cases_150_scores) + length(platelet_controls_150_scores)
weights_platelet = c(weight_platelet_125, weight_platelet_150)

weight_explained_125 = length(explained_125_scores) + length(controls_and_unexp_125_scores)
weight_explained_150 = length(explained_150_scores) + length(controls_and_unexp_150_scores)
weights_explained = c(weight_explained_125, weight_explained_150)

weight_ird_125 = length(ird_cases_125_scores) + length(ird_controls_125_scores)
weight_ird_150 = length(ird_cases_150_scores) + length(ird_controls_150_scores)
weights_ird = c(weight_ird_125, weight_ird_150)


dot_plot_data_ranksum = data.table(GROUP=c("Expl_125", "Expl_150", "Plat_125", "Plat_150", "IRD_125", "IRD_150", "Expl_meta", "Plat_meta", "IRD_meta"),
																	 PVAL=rep(-1, 9))


# conduct rank sum test and use Stouffer's Z-score to combine independent tests
dot_plot_data_ranksum[GROUP=="Expl_125", PVAL:=wilcox.test(explained_125_scores, controls_and_unexp_125_scores, alternative="greater")$p.value]
dot_plot_data_ranksum[GROUP=="Expl_150", PVAL:=wilcox.test(explained_150_scores, controls_and_unexp_150_scores, alternative="greater")$p.value]

dot_plot_data_ranksum[GROUP=="Plat_125", PVAL:=wilcox.test(platelet_cases_125_scores, platelet_controls_125_scores, alternative="greater")$p.value]
dot_plot_data_ranksum[GROUP=="Plat_150", PVAL:=wilcox.test(platelet_cases_150_scores, platelet_controls_150_scores, alternative="greater")$p.value]

dot_plot_data_ranksum[GROUP=="IRD_125", PVAL:=wilcox.test(ird_cases_125_scores, ird_controls_125_scores, alternative="greater")$p.value]
dot_plot_data_ranksum[GROUP=="IRD_150", PVAL:=wilcox.test(ird_cases_150_scores, ird_controls_150_scores, alternative="greater")$p.value]

dot_plot_data_ranksum[GROUP=="Expl_meta", PVAL:=combine.test(p=c(dot_plot_data_ranksum[GROUP=="Expl_125", PVAL], dot_plot_data_ranksum[GROUP=="Expl_150", PVAL]),
														weight=weights_explained, method="z.transform")]
dot_plot_data_ranksum[GROUP=="Plat_meta", PVAL:=combine.test(p=c(dot_plot_data_ranksum[GROUP=="Plat_125", PVAL], dot_plot_data_ranksum[GROUP=="Plat_150", PVAL]),
														weight=weights_platelet, method="z.transform")]
dot_plot_data_ranksum[GROUP=="IRD_meta", PVAL:=combine.test(p=c(dot_plot_data_ranksum[GROUP=="IRD_125", PVAL], dot_plot_data_ranksum[GROUP=="IRD_150", PVAL]),
														weight=weights_ird, method="z.transform")]


dot_plot <- ggplot(dot_plot_data_ranksum) + geom_dotplot(aes(x=GROUP, y=PVAL), dotsize=0.8, binaxis='y', stackdir='center') + 
	theme_bw() + ylim(0, 1) + scale_y_continuous(trans='log10')

pdf("/home/jwt44/dot_plot.pdf")
dot_plot
dev.off()


############ logistic test (thresholding) ##############

# platelet

batch = c(rep("125", length(platelet_cases_125_scores) + length(platelet_controls_125_scores)),
				  rep("150", length(platelet_cases_150_scores) + length(platelet_controls_150_scores)))

case_control_status = c(rep(1, length(platelet_cases_125_scores)),
												rep(0, length(platelet_controls_125_scores)),
												rep(1, length(platelet_cases_150_scores)),
												rep(0, length(platelet_controls_150_scores)))

scores = c(platelet_cases_125_scores, platelet_controls_125_scores,
					 platelet_cases_150_scores, platelet_controls_150_scores)

pvals = numeric(99)
thresholds = seq(1, 99, 1)

for (threshold in thresholds) {
	fit = glm(case_control_status ~ factor(scores>threshold) + factor(batch), family=binomial)
	pvals[threshold] = summary(fit)$coefficients[2, 4]
}

min_threshold = which.min(pvals)
min_pval = pvals[min_threshold]

fit_min = glm(case_control_status ~ factor(scores>min_threshold) + factor(batch), family=binomial)
M = confint(fit_min)
CI_lower = exp(M[2,1])
CI_upper = exp(M[2,2])
OR_pred_min = exp(summary(fit_min)$coefficients[2,1])

print(paste(min_pval, CI_lower, OR_pred_min, CI_upper))


# explained

batch = c(rep("125", length(explained_125_scores) + length(controls_and_unexp_125_scores)),
				  rep("150", length(explained_150_scores) + length(controls_and_unexp_150_scores)))

case_control_status = c(rep(1, length(explained_125_scores)),
												rep(0, length(controls_and_unexp_125_scores)),
												rep(1, length(explained_150_scores)),
												rep(0, length(controls_and_unexp_150_scores)))

scores = c(explained_125_scores, controls_and_unexp_125_scores,
					 explained_150_scores, controls_and_unexp_150_scores)

fit_expl = glm(case_control_status ~ factor(scores>min_threshold) + factor(batch), family=binomial)
pval_expl = summary(fit_expl)$coefficients[2, 4]
M_expl = confint(fit_expl)
CI_lower_expl = exp(M_expl[2,1])
CI_upper_expl = exp(M_expl[2,2])
OR_pred_expl = exp(summary(fit_expl)$coefficients[2,1])

print(paste(pval_expl, CI_lower_expl, OR_pred_expl, CI_upper_expl))



# IRD

batch = c(rep("125", length(ird_cases_125_scores) + length(ird_controls_125_scores)),
				  rep("150", length(ird_cases_150_scores) + length(ird_controls_150_scores)))

case_control_status = c(rep(1, length(ird_cases_125_scores)),
												rep(0, length(ird_controls_125_scores)),
												rep(1, length(ird_cases_150_scores)),
												rep(0, length(ird_controls_150_scores)))

scores = c(ird_cases_125_scores, ird_controls_125_scores,
					 ird_cases_150_scores, ird_controls_150_scores)

pvals = numeric(99)
thresholds = seq(1, 99, 1)

fit_ird = glm(case_control_status ~ factor(scores>min_threshold) + factor(batch), family=binomial)
pval_ird = summary(fit_ird)$coefficients[2, 4]
M_ird = confint(fit_ird)
CI_lower_ird = exp(M_ird[2,1])
CI_upper_ird = exp(M_ird[2,2])
OR_pred_ird = exp(summary(fit_ird)$coefficients[2,1])

print(paste(pval_ird, CI_lower_ird, OR_pred_ird, CI_upper_ird))



ci_plot_data = data.table(case_group = c("platelet", "explained", "ird"), 
					 								pred_or = c(OR_pred_min, OR_pred_expl, OR_pred_ird),
					 								ci_lower = c(CI_lower, CI_lower_expl, CI_lower_ird),
					 								ci_upper = c(CI_upper, CI_upper_expl, CI_upper_ird))

CI_plot <- ggplot(ci_plot_data, aes(x = case_group, y = pred_or, colour = case_group)) + geom_point(size = 4) +
  geom_errorbar(aes(ymax = ci_upper, ymin = ci_lower, colour = case_group)) + guides(colour=FALSE) + 
  xlab("Case group") + ylab("Exponentiated coefficient") + theme_bw() + ylim(0.89, 1.11)

pdf("/home/jwt44/ci_plot_logistic.pdf")
CI_plot
dev.off()

