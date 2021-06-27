
require(ggplot2)
require(ggpubr)
require(data.table)


dt = fread("/home/jwt44/case_control_regression_data.tsv")
dt = dt[!is.na(median_SCARS_clf)]

coefs = c()
CI_lower = c()
CI_upper = c()
pvals = c()
thresholds = seq(95, 99.9, 0.1)

for (threshold in thresholds) {
	dt[, SCARS_id := factor(median_SCARS_clf>threshold)]
	
	fit = with(dt, glm(case_control ~ SCARS_id, family=binomial))
	coefs = c(coefs, exp(summary(fit)$coefficients[2, 1]))
	pvals = c(pvals, summary(fit)$coefficients[2, 4])

	M = confint(fit)
	CI_lower = c(CI_lower, exp(M[2,1]))
	CI_upper = c(CI_upper, exp(M[2,2]))
}


df <- data.frame(x = thresholds,
                 F = coefs,
                 L = CI_lower,
                 U = CI_upper,
                 P = pvals,
				 S = pvals < 0.05)


CI_plot <- ggplot(df, aes(x = x, y = F, colour = S)) + geom_point(size = 4) +
  geom_errorbar(aes(ymax = U, ymin = L, colour = S)) + guides(colour=FALSE) + 
  xlab("Median SCARS-clf threshold") + ylab("Exponentiated coefficient") + theme_bw() + coord_equal()

PVAL_plot <- ggplot(df, aes(x = x, y = P)) + geom_line() + 
	xlab("Median SCARS-clf threshold") + ylab("P-value") + theme_bw() +  
	scale_y_continuous(trans='log10') + coord_equal()

pdf("/home/jwt44/regression_CI_analysis.pdf", height=5, width=10)
ggarrange(CI_plot, PVAL_plot, ncol=2, nrow=1)
dev.off()
