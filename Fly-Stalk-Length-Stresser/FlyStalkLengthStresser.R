load('fly.RData')
mu_stressor1 <- mean(fly_stressor1)
mu_stressor2 <- mean(fly_stressor2)
effect_size_stressors <- mu_stressor2 - mu_stressor1
t.test(fly_stressor1,fly_stressor2, paired = FALSE, var.equal = TRUE)
  #p-value is 0.00000000666 which is very small indicating that the results are statistically significant and there is strong evidence against the null hypothesis.

wilcox.test(fly_stressor1,fly_stressor2, paired = FALSE, var.equal = TRUE)
#the non-parametric test also gives a very small of 8.857e-07 indicates a statistically significant results and strong evidence against the null
#Hence even if the assumptions of the parametric test donot hold, the wilcoxon rank gives a satisfying p-value

par(mfrow=c(2,2))
#Assumptions for parametric test:
  "Variances of the two samples are equal. The samples are normally distributed"
#In order to check these assumptions:
hist(fly_stressor1, col = "purple")
hist(fly_stressor2, col = "green")
#The histograms are appproximately normally distributed with the fly_stressor2 having a peak at around 1.5 
#the fly_stressor2 has a peak at around 2.0
#the histogram shows visually that an approximate normal distribution exists with the a peak and then symmetry on both sides
stripchart(fly_stressor1, col = "purple")
stripchart(fly_stressor2, col = "green")
#The stripcharts show that the spread of the values is approximately similar because the scatter looks similar.
#this indicates that the two samples have an approximate equal variance
#therefore the assumptions are met for this parametric test
