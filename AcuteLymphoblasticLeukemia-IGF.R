
library("ISwR")
juul
str(juul)
#This shows the structure that there are 1339 observations and 6 variables: age as numeric, menarche as int, sex as numeric, igf1 as numeric, tanner as int, testvol as int.
head(juul)
#This shows the first six observations, where most of the values are not recorded


par(mfrow = c(2,2))
#igf1 is already numeric, while tanner is int so it is converted using the as.factor()
plot(lm(formula = igf1 ~as.factor(tanner) , data = juul), which= 1:4)
#The residuals vs. fitted shows that the spread of the variables across the line y=0 is even,
#Therefore there is no evidence of any non-linear relation between the two variables.

#The normal q-q plot shows that the residuals are almost normally distributed.
#There is an almost straight line observed, therefore the residuals are normally distributed

#The Scale-location plot shows that the variance of the residuals is the same across the explanatory variable
#Because there is an almost horizontal line, the residuals have the same variance across the explanatory variable

#The Cook's distance plot shows that the values at observation numbers which have the most effect on the relation
#918, 355 and 335 observation numbers have the most effect on the linear relationship of the response and the explanatory variable.

anova(lm(igf1~as.factor(tanner),data=juul))
#without using the multiple testing correction, the omnibust p-value is very small 2.2e-16.
#hence there is strong evidence of a correlation between the response and the explanatory variable.
#Then usin the pairwise t.test and using FWER (both Holm and Bonferroni)
pairwise.t.test(juul$igf1,as.factor(juul$tanner),data=juul, p.adj = "bonferroni")
#The Bonferroni and the Holm both give similar results. The p-value is very small for comparing the factor1 with all other factors.
#Both in Bonferroni and Holm it is 2e-16, also factor2 with other factors is very small in the e-08 exponents. 
#Therefore the factor2 and factor1 are very different from other factors and have a difference of means that is large
#factor 3 however has a very high p-value of 1 in Bonferroni and 0.4 in Holm. Therefore this is strong evidence that the factor 3 is not different from the other factors
#The factor 4 and 5 also have a relatively different mean.
pairwise.t.test(juul$igf1,as.factor(juul$tanner),data=juul, p.adj = "holm")

#Then the Tukey HSD is used. This is to get the 95% confidence intervals for the difference of means of the factors.
val <- aov(lm(igf1~as.factor(tanner),data=juul))
TukeyHSD(val)

plot(TukeyHSD(val))
#This is illustrated in the graph as an interval.
#The intervals that cross the dotted line for zero indicate that the pairs have similar means. 
#The factors 3-1, 4-1, 5-1 are very far apart from the dotted line
#The 5-4, 4-3, 5-3 have very high p-values and cross the dotted line for 0, therefore their means are similar.

#======================================================================================================
#It is an observational study because just the gene expression of the people suffering from acute lymphoblastic leukemia and normal people was studied.
#There was no other extermal variable introduced
library("multtest","qvalue")
load("ALLsfilt.RData")
str(ALLsfilt)
gene_names <- dimnames(ALLsfilt)[1]
#The structure that it has 8812 rows of different genes and 79 columns of 79 different individuals
mol.biol
#initializing the cancer vectors and the control vectors to null
cancer_vector = c()
control_vector = c()
#running the loop to check which rows correspond to the control and which to the cancer
for (i in 1:length(mol.biol)){
  if (mol.biol[i] == "NEG"){
    control_vector <- append(x=control_vector,values=i)
  }
  else{
    cancer_vector <- append(x=cancer_vector,values=i)
  }
}
#print out the vectors to see which has what indeces of the row numbers
cancer_vector
control_vector
#initializing the result vector to null
result <- NULL
#running the loop where the t-test is conducted on the values of that gene at the normal and the defective level
for (i in 1:nrow(ALLsfilt)){
  result[i]<- t.test(ALLsfilt[i,cancer_vector],ALLsfilt[i,control_vector], var.equal = TRUE, paired = F)$p.value
}

#creating a histogram of the p-values collected in the result vector
hist(result,col = "blue", breaks = 30 )
#the features of the distribution include that it has a background of null distribution but a peak close to zero. 
#This shows that there are about 1100 genes that are deferentially expressed in a normal and cancer person
#All the other genes are divided into 400-500 sets that have a uniform p-value distribution and come form the null distribution

library(multtest)
library(qvalue)
#Thhis mt is to check all the corrections: Bonferroni, Holm and BH.
mt <- mt.rawp2adjp(result, proc=c("Bonferroni","Holm","BH"))
#this is to get the q-values, adjusted p-values, from the q-values func in the library
q.storey <- qvalue(sort(result))$qvalues
#This is to create a data frame of the adjuted p-values obtained from all the corrections
d_null_fdr <- data.frame(mt$adjp, q.storey)
#This is just to check the first three rows of the matrix and check for the adjusted p-values.
d_null_fdr[1:3,]
#Then the number of genes that are differentially expressed and have less tha 20% false positives are calculated using the apply method below
apply(d_null_fdr, 2, function(x) sum(x < 0.2))
# this shows that the FWER Bonferroni and Holm correction is very conservative and gives only 42 differentially expressed genes with p-value less than 0.2
#This is because the error is bound to the whole family of tests 
#The raw-p-values are very small and give a huge 2713 differentially expressed genes with multiple testing correction
#Therefore the FDR are less conservaive and gives meaningful error control with large numbers of genes tested
#The adjusted result below is therefore set to the FDR Benjamini-Hockberg column of p-values.
adjusted_result <- d_null_fdr$BH
adjusted_result


#set the list of differentially expressed genes to be empty
list_differential_pval <- c()
#Run a loop and if the p-value is less than 0.2, add it to the list of differentially expressed
for (i in 1:length(adjusted_result)){
  if (adjusted_result[i]<0.2){
    list_differential_pval<-append(list_differential_pval,adjusted_result[i])
  }
}
#The list of p-values that are less than 0.2
list_differential_pval
length(list_differential_pval)
#The length of this vector is 593 which is the same number of genes showed by the apply() method as differentially expressed


raw_pvalue <- result
gene_names
set_name <- data.frame(raw_pvalue,gene_names)
set_name = order(set_name)
#The first 593 numbers are the row numbers in the ordered set, that give the differentially expressed genes ordered in ascending order
list_gene_index = NULL
for (i in 1:593){
  list_gene_index = append(list_gene_index,set_name[i])
}
list_gene_index
#This list contains the list of indexes of the genes in the ALLsfilt data that are differentially expressed.




