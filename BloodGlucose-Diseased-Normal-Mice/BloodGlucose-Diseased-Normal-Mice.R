library(ggplot2)

load( 'midterm_mice.RData')
attach(mtcars)
par(mfrow=c(2,2))

x <- list("normal" = normal , "diseased" = diseased)

#strip-chart
stripchart(x,
           main="Strip-chart",
           xlab="Blood Glucose levels (mmol/L)",
           ylab="Population",
           method="jitter",
           col=c("darkturquoise","firebrick1"),
           pch=16
)
#box-plot
boxplot(x, main = "Box-plot", pch = 16, col = c("darkturquoise","firebrick1"),  ylab="Blood Glucose levels (mmol/L)",
        xlab="Population")
#The elements of a box-plot include the median that is the thick black line and the line blow it is the first quartile and third quartile is above it.
#The median tells the middle value from the data set.
#The first quartile tells the value of the 25th pertentile and the third quartile is the value of the 75th percentile
#The ends of a box-plot show the maximum and the minimum value. 
#The advantages of using box-plot is that it points out where the median lies and how many values lie between the quartiles.
#It is used if the statistics of the data-set are to be exactly demostrated using actual values

#A strip plot can be used to visualize the spread of data.
#It can be used to tell where most of the scatter of one data set lies as compared to another.

#A stacked histogram can tell the difference of the values in the actual data-sets by observed height of stacked bars
#Instead of giving two separate datasets and visualizing them to tell the difference, stacked histogram is just one figure
#It can be used to tell which data-set contains values that are more extreme

#stacked histogram
dd <- data.frame(Blood_Glucose_Levels = c(normal, diseased),dist = c(rep("normal", length(normal)), rep("diseased", length(diseased))))
ggplot(dd, aes(x = Blood_Glucose_Levels, fill = dist)) +
  geom_histogram(position = "stack", binwidth = 0.5, col = 'black')

#the effect-size is the difference of the means
effect_size <- mean(normal) - mean(diseased)

#Comparison of data
#As illustrated in the strip chart the scatter for diseased population lies in the higher blood glucose region around 8-10mmol/L
#While the scatter for normal population lies in the 4-6mmol/L region
#As illustrated in the box-plot, the the median for normal is around 4 with a minimum and maximum around -2 and 11 respectively
#Whereas the median for the diseased is around 7 with a minimuma and a maximum around 0 and 14
#therefore these plots can be helpful in visually interpreting the statistics of the data


#vector containing the data from both the diseased and nomal population
combine_vec <- c()

for (i in 1:length(normal)){
  combine_vec[i] = normal[i]
}
ind <- length(normal) +1
for (i in diseased){
  combine_vec[ind] = i
  ind <- ind+1
}

#PART 1
# parametric test
# two-sample (pooled variance) t-test
df1 = length(normal)-1
df2 = length(diseased)-1
pooled_variance <- (df1*var(normal)+df2*var(diseased))/(df1+df2)
SEM <- sqrt(pooled_variance) * sqrt(((1/length(normal))+(1/length(diseased))))
test.stat <- (mean(normal)-mean(diseased))/SEM
test.stat

# p.value
df = df1+df2
p_value = pt(test.stat,df=df)*2 
#(assuming two-sided t-test)
# p = 4.2206 e-06
#statistically significant results because p<0.005. There is strong evidence against the null. 

#the assumptions made in this parametric test is that the sample is collected randomly from the population. (unbiased)
#Another assumption is that the population has a normal distribution (Gaussian) for the parameter
#Another assumption is that the variations of both the data sets are equal
#this means that the limits of the method include that it cannot be used for population data that is not normally distributed.
#Another limitation is that it can not be used for small sample sizes.
#The advantages are such that it is a very powerful test that can give accurate results if the limitations are met.

#checking the assumptions:
attach(mtcars)
par(mfrow=c(3,2))
hist(combine_vec, col = "orange")
stripchart(combine_vec, col ="orange")
boxplot(combine_vec, col = "orange")
#the assumptions for the parametric test are met as the histogram of the combine_vec does have a symmetry and a normal distribution.
#other assumption is that the variances are equal
#the strip chart shows condensed scatter in the middle and then fewer spots around. Hence the data is approximately gaussian distributed.
hist(diseased, col = "red")
hist(normal, col = "blue")
stripchart(diseased, col = "red")
stripchart(normal, col = "blue")
#these pair of histograms and the stripcharts can show that the variances of both the data sets are almost equal
#as the spread of values aroung the symmetry in the histogram and the strip chart is visually equal


#PART 2
# two-sample Wilcoxon 2-sample rank-sum
wilcox.test(normal, diseased, paired = FALSE)
# p-value = 0.000004508
#statistically significant results, p<0.005, strong evidence against the null

#the assumptions made in this non-parametric test is that the sample is collected randomly from the population (unbiased)
#the limits of the method include that they are not very powerful and donot always give trustworthy results
#therefore it is better to perform a parametric test if their assumptions are met
#the advantages are such that it does not assume a gaussian (normal) distribution so can be used for data with any distrubution
#Moreover it is robust to outliers because the ranks are compared not the mean


#PART 3
numb_extreme <- 0
#a vector containing the mean differences of between the diseased and normal population of 1000 permutations
mean_diff_vector <- c()
for (i in 1:1000){
  new_vec <- sample(combine_vec, replace = FALSE)
  normal_half <- new_vec[1:(length(normal))]
  diseased_half <- new_vec[(length(normal)+1):length(new_vec)]
  mu_normal <- mean(normal_half)
  mu_diseased <- mean(diseased_half)
  mean_diff_vector[i] <- (mu_normal-mu_diseased)
}

for (i in mean_diff_vector){
  if (i >= ( abs(mean(normal) - mean(diseased) ) ) ){
    numb_extreme <- numb_extreme + 1
  }
}
#fraction of extreme values from the 1000 iterations
p_value <- numb_extreme/1000
#p_value = 0. This is a statistically significant result, and there is strong evidence against the null.

#this histogram can demonstrate the p_value in the graph of all the permutations of the sample
hist(mean_diff_vector, breaks = 40, xlim = c(-3,3), col = "blue")
abline(v = effect_size, col = "red")
#the area (probability) of getting result observed or more extreme is observed to be zero

#the assumptions made in this non-parametric test is that the sample is collected randomly from the population (unbiased)
#the limits of the method include that they are not very powerful and donot always give trustworthy results
#therefore it is better to perform a parametric test if their assumptions are met
#the advantages are such that it does not assume a gaussian (normal) distribution so can be used for data with any distrubution
#Moreover it can be used for any test-statistic



