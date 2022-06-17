#Autophagy random forest model classification

setwd("/media/mohaksharda/MyPassport/Somatic_Cancer_analysis_wo_merging_technical_replicate_bams_17_may_2021/scripts/autophagy_analysis/")

library(dplyr)
library(ggplot2)
library(randomForest)
library(ROCR)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)

#Create the dataframe
autophagy_table <- as.data.frame(read.table("merged_transposed_expression_with_smoking_status_1_2_only_autophagy_analysis.tsv",sep="\t",header=T))
autophagy_table<-select(autophagy_table, -Hugo_Symbol)
summary(autophagy_table)
autophagy_table <- transform(autophagy_table,
                             smoking_status=as.factor(smoking_status))
autophagy_table[,-370] <- log(1+autophagy_table[-370],2) #take log of all columns except the last one with smoking_status (370th column)
autophagy_table[,-370] <- scale(autophagy_table[-370])
head(autophagy_table)
write.csv(x = autophagy_table,file = "autophagy_genes_expression_converted_to_zscores.csv")
#split data
# Set random seed to make results reproducible:
set.seed(17)
# Calculate the size of each of the data sets:
data_set_size <- floor(3*nrow(autophagy_table)/4)
# Generate a random sample of "data_set_size" indexes
indexes <- sample(1:nrow(autophagy_table), size = data_set_size)
# Assign the data to the correct sets
autophagy_table_train <- autophagy_table[indexes,]
autophagy_table_test <- autophagy_table[-indexes,]
autophagy_table_train

#Random Forest Modelling
rf_classifier1 = randomForest(smoking_status ~ ., data=autophagy_table_train, ntree=1500, importance=TRUE)

rf_classifier1

library(e1071)
x<-as.data.frame(rf_classifier1$importance)
str(x)
x[order(-x$MeanDecreaseGini),]
x[order(-x$MeanDecreaseAccuracy),]
varImpPlot(rf_classifier1) #feature importance plot

prediction_for_table <- predict(rf_classifier1,autophagy_table_test[,-370])
table(observed=autophagy_table_test[,370],predicted=prediction_for_table)
cm_prediction <- confusionMatrix(prediction_for_table,autophagy_table_test[,370],positive="2")
cm_prediction #prints all metrics

# Validation set assessment #2: ROC curves and AUC

# Calculate the probability of new observations belonging to each class
# prediction_for_roc_curve will be a matrix with dimensions data_set_size x number_of_classes
prediction_for_roc_curve <- predict(rf_classifier1,autophagy_table_test[,-370],type="prob")
prediction_for_roc_curve

# Use pretty colours:
pretty_colours <- c("#00BA38","#F8766D") #green:non-smokers, red:smokers

# Specify the different classes 
classes <- levels(autophagy_table_test$smoking_status)

# For each class
for (i in 1:2)
{
  # Define which observations belong to class[i]
  true_values <- ifelse(autophagy_table_test[,370]==classes[i],1,0)
  # Assess the performance of classifier for class[i]
  pred <- prediction(prediction_for_roc_curve[,i],true_values)
  perf <- performance(pred, "tpr", "fpr")
  if (i==1)
  {
    plot(perf,main="ROC Curve",col=pretty_colours[i]) 
  }
  else
  {
    plot(perf,main="ROC Curve",col=pretty_colours[i],add=TRUE) 
  }
  # Calculate the AUC and print it to screen
  auc.perf <- performance(pred, measure = "auc")
  print(auc.perf@y.values)
}

#selected top autophagy genes boxplot and wilcoxon test to compare expression between smokers and non-smokers

#autophagy_table <- autophagy_table[, c('BIRC5', 'BNIP3', 'CLEC16A',
#'DNM1L','EIF4EBP1','FEZ1','MUL1',
#'PRKD1','RIPK2','SESN1','SESN2','TNFSF10',
#'TRIM22','smoking_status')]

plot_gene_plot <- autophagy_table %>% ggplot(aes(x=smoking_status,y=BNIP3,fill=factor(smoking_status))) + 
  geom_boxplot(width=0.5,lwd=1) + 
  geom_jitter(width=0.15)+ 
  labs(subtitle="BIRC5 expression")

plot_gene_plot + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black"))

autophagy_table %>% 
  select(TNFSF10,BIRC5, BNIP3, CLEC16A,DNM1L,EIF4EBP1,FEZ1,MUL1,PRKD1,RIPK2,SESN1,SESN2,TRIM22,smoking_status) %>% 
  pivot_longer(.,cols=c(TNFSF10,BIRC5, BNIP3, CLEC16A, DNM1L,EIF4EBP1,FEZ1,MUL1,PRKD1,RIPK2,SESN1,SESN2,TRIM22), names_to = "Top_Autophagy_genes", values_to = "Expression_values") %>% 
  ggplot(aes(x= Top_Autophagy_genes, y= Expression_values, fill = smoking_status)) + geom_boxplot() + stat_compare_means()

autophagy_table %>% 
  select(TNFSF10,BIRC5,smoking_status) %>% 
  pivot_longer(.,cols=c(TNFSF10,BIRC5), names_to = "Top_Autophagy_genes", values_to = "Expression_values") %>% 
  ggplot(aes(x= Top_Autophagy_genes, y= Expression_values, fill = smoking_status)) + geom_boxplot() + stat_compare_means()

wilcox.test(TRIM22 ~ smoking_status, data = autophagy_table,
            exact = FALSE)