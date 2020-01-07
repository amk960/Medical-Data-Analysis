a <- read.table(file = "q3_calcium_regression.tsv", header = TRUE, col.names = c("density", "metabolite"))
plot(a)
#it is a scatter plot that indicates a negative correlation between density and metabolite
#since this is visually linear, the best fitting line equation will be y = mx + c 
#in this equation the m is the slope and the c is the intercept
abline(lm(a[["metabolite"]]~a[["density"]]))
summary(lm(a[["metabolite"]]~a[["density"]]))
#Since the adjusted R-squared value is 0.7554, it is a value smaller than the perfect fit 1. It is however close to 1.
# Because the R-squared is greater than 0.5 hence there is a strong correlation between the variables. 
#the R-squared value indicates that the 75% of the variations in the data values collected is explained
#the p-value is 0.0003152 which indicates a strong evidence against the null and the results are statistically significant
#there is strong evidence of a strong association

#The assumptions for a linear regresssion are:
" -> linear relationship between X and Y
  -> dependent variable metabolite is contiouous
  -> residual is normally distributed
  -> variance of residual is same across X
"
#Runnig diagnostics
par(mfrow=c(2,2))
plot(lm(a), which=1:4)

# the residual vs. fitted shows a shows a spread around y=0 without any distinst pattern. 
#Hence there is not any non-linear relationship between the two variables not adressed in the regression
#The normal q-q plot checks of the residuals are normally distributed
#because there is a straight line observed in the normal q-q plot, the residuals are normally distributed
#The scale-location line shows if the variance of residual is same across X
#Because the line is not nearly horizontal, the variance of residual is not quite the same across X
#The Cook's distance plot shows which observation in the data has the most influence on the regression line
#it is shown in the cook's distance plot to be at 9.

#Hence there is no non-linear relation between the variables and the residuals are normally distributed. 
#The observation at X=9 is an outlier but otherwise the data has strong evidence of a positive linear correlation
summary(lm(a[["metabolite"]]~a[["density"]]), data=a, subset=-9)
#the summary indicates that even after removing 9, the p_value is still 0.0003152 which is less than 0.005.
#there is still strong evidence for a linear relationship between the variables
