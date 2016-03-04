###################################################################################
### Evaluation: GEA and PRR for drug classes
###################################################################################
#
# This file will create ROC Plots for each OMOP outcome:
# Additional File 1 in BMC publication 
#
# This file will evaluate performance of GEA and PRR
# at different levels and thresholds:
# Additional File 2 in BMC publication
#
library(GEA)
library(ROCR)


##################################################################################
### Helper function for evaluating TP, TN, FP, FN, precison, recall, F1
##################################################################################
getPerformance <- function(results, level, thr, reverse=FALSE) {
  
  omop_result_test <- results[[level]]
  
  if (reverse==TRUE) {
    TP = sum(omop_result_test[omop_result_test$set=="POS","pvalue"]<thr )
    TN = sum(omop_result_test[omop_result_test$set=="NEG","pvalue"]>=thr )
    FP = sum(omop_result_test[omop_result_test$set=="NEG","pvalue"]<thr )
    FN = sum(omop_result_test[omop_result_test$set=="POS","pvalue"]>=thr )
  } else {
    TP = sum(omop_result_test[omop_result_test$set=="POS","pvalue"]>thr )
    TN = sum(omop_result_test[omop_result_test$set=="NEG","pvalue"]<=thr )
    FP = sum(omop_result_test[omop_result_test$set=="NEG","pvalue"]>thr )
    FN = sum(omop_result_test[omop_result_test$set=="POS","pvalue"]<=thr )   
  }
  
  PR = TP/(TP+FP)
  RC = TP/(TP+FN)
  F1 = 2*(PR*RC)/(PR+RC)
  
  c(TP,TN,FP,FN,PR,RC,F1)
}

##################################################################################
### Helper function for ROC Plots for each OMOP outcome
##################################################################################
# (x-axis: fpr, y-axis: tpr)
plotROC <- function(plottype = "AMI", title=plottype) {
  
  getPerformanceObject <- function(results, level, plottype, reverse=FALSE) {
    
    omop_result_test <- results[[level]]
    
    if (reverse==TRUE) {
      pred710 <- prediction( -1*omop_result_test[omop_result_test$Outcome==plottype,"pvalue"], omop_result_test[omop_result_test$Outcome==plottype,"set_bin"])
    } else {
      pred710 <- prediction( omop_result_test[omop_result_test$Outcome==plottype,"pvalue"], omop_result_test[omop_result_test$Outcome==plottype,"set_bin"])  
    }
    
    tprfpr <- performance(pred710,"tpr","fpr")
    auc    <- round(performance(pred710,"auc")@y.values[[1]],3)
    
    return(list("tprfpr"=tprfpr, "auc"=auc))
  }
  
  pred710 <- getPerformanceObject(results, "7-10",  plottype, reverse=TRUE)
  pred457 <- getPerformanceObject(results, "4.5-7", plottype, reverse=TRUE)
  pred155 <- getPerformanceObject(results, "1-5.5", plottype, reverse=TRUE)
  pred2nd <- getPerformanceObject(results, "2nd",   plottype, reverse=TRUE)
  
  predprr710 <- getPerformanceObject(results, "PRR7-10",  plottype)
  predprr457 <- getPerformanceObject(results, "PRR4.5-7", plottype)
  predprr155 <- getPerformanceObject(results, "PRR1-5.5", plottype)
  predprr2nd <- getPerformanceObject(results, "PRR2nd",   plottype)
  
  plot( pred710$tprfpr,                       col = "red",  lwd = 2, main=title)
  plot( pred457$tprfpr,    add = TRUE, lty=2, col = "red",  lwd = 2)
  #plot( pred155$tprfpr,    add = TRUE, lty=3, col = "red",  lwd = 2)
  #plot( pred2nd$tprfpr,    add = TRUE, lty=4, col = "red",  lwd = 2)
  
  plot( predprr710$tprfpr, add = TRUE,        col = "blue", lwd = 2)
  #plot( predprr457$tprfpr, add = TRUE, lty=2, col = "blue", lwd = 2)
  #plot( predprr155$tprfpr, add = TRUE, lty=3, col = "blue", lwd = 2)
  plot( predprr2nd$tprfpr, add = TRUE, lty=4, col = "blue", lwd = 2)
  
  abline(a=0, b= 1)
  legend(0.5,0.2,  lwd=c(2.5,2.5), 
         c(paste0("gEA 7-10  - AUC ", pred710$auc   , collapse = NULL),
           paste0("gEA 4.5-7 - AUC ", pred457$auc   , collapse = NULL),
           #paste0("gEA 1-5.5 - AUC ", pred155$auc   , collapse = NULL),
           #paste0("gEA 2nd   - AUC ", pred2nd$auc   , collapse = NULL),
           paste0("PRR 7-10  - AUC ", predprr710$auc, collapse = NULL),
           #paste0("PRR 4.5-7 - AUC ", predprr457$auc, collapse = NULL),
           #paste0("PRR 1-5.5 - AUC ", predprr155$auc, collapse = NULL),
           paste0("PRR 2nd   - AUC ", predprr2nd$auc, collapse = NULL)
         ),
         col=c("red","red","blue","blue" ), lty=1:4)
}


##################################################################################
### Start Evaluation 
##################################################################################

##################################################################################
### Load drug classes and adverse event terms
##################################################################################

drugs <- Drugs("drugs_ATC_MEDLINE_360k.txt", classLevel=TRUE)
events <- AdverseEvents("manifestations_MEDLINE_360k.txt")


##################################################################################
### Create abstraction levels
##################################################################################

abs710    <- getAbstractionLayer(drugs,events,c(7,10))
abs457    <- getAbstractionLayer(drugs,events,c(4.5,7))
abs155    <- getAbstractionLayer(drugs,events,c(1,5.5))
abs2nd    <- getAbstractionLayer(drugs,events,fix=2)

absprr710 <- getAbstractionLayer(drugs,events,c(7,10))
absprr457 <- getAbstractionLayer(drugs,events,c(4.5,7))
absprr155 <- getAbstractionLayer(drugs,events,c(1,5.5))
absprr2nd <- getAbstractionLayer(drugs,events,fix=2)

##################################################################################
### Create Gold Standard (will read OMOP data from flat file)
##################################################################################

gold <- GS("OMOP_GS_drug_classes.csv")


################################################################################## 
# Perform Enrichment analysis and return pvalue for each drug - outcome pair 
# in Gold Standard
# EA will be performed at 4 different abstraction levels
################################################################################## 

omop_result_all710     <- data.frame(gold,  pvalue=apply( gold , 1, function(x) { getPValue(x[6],x[2],abs710)}))
omop_result_all457     <- data.frame(gold,  pvalue=apply( gold , 1, function(x) { getPValue(x[6],x[3],abs457)}))
omop_result_all155     <- data.frame(gold,  pvalue=apply( gold , 1, function(x) { getPValue(x[6],x[4],abs155)}))
omop_result_all2nd     <- data.frame(gold,  pvalue=apply( gold , 1, function(x) { getPValue(x[6],x["X2nd"],abs2nd)}))

################################################################################## 
# Calculate PRR for each drug - outcome pair in Gold Standard
# PRR will be calculated at 4 different abstraction levels
##################################################################################

omop_result_prr_all710 <- data.frame(gold,  pvalue=apply( gold , 1, function(x) { getPRR(x[6],x[2],absprr710)}))
omop_result_prr_all457 <- data.frame(gold,  pvalue=apply( gold , 1, function(x) { getPRR(x[6],x[3],absprr457)}))
omop_result_prr_all155 <- data.frame(gold,  pvalue=apply( gold , 1, function(x) { getPRR(x[6],x[4],absprr155)}))
omop_result_prr_all2nd <- data.frame(gold,  pvalue=apply( gold , 1, function(x) { getPRR(x[6],x["X2nd"],absprr2nd)}))


##################################################################################
# Combine individual results and plot
##################################################################################

results<- list(   "7-10" = omop_result_all710,       "4.5-7" = omop_result_all457,        "1-5.5" = omop_result_all155,        "2nd"= omop_result_all2nd, 
                  "PRR7-10" = omop_result_prr_all710,"PRR4.5-7" = omop_result_prr_all457, "PRR1-5.5" = omop_result_prr_all155, "PRR2nd"= omop_result_prr_all2nd)

##################################################################################
# Plot curves for each outcome (Additional File 1)
# Of not, for GEA levels 1-5.5 and 2nd are currently uncommented in plot function
# For PRR levels 4.5-7 and 1-1.5 are uncommented
##################################################################################

pdf("plot.pdf" , height=12, width=12,pointsize=13)
par(mfrow=c(2,2),mar=c(5.1,5.1,4.1,2.1))
plottype = "AKI"
plotROC(plottype, "Acute Kidney Injury")

plottype = "ALI"
plotROC(plottype, "Acute Liver Injury")

plottype = "AMI"
plotROC(plottype, "Acute Myocardial Infarction")

plottype = "GI"
plotROC(plottype, "GI Bleed")
dev.off()


##################################################################################
# Create performance table for each method, layer, and threshold
# (Additional Table 2)
##################################################################################

level = "7-10"
tableAdditional2<- cbind(        getPerformance(results, level, 0.1, reverse=TRUE), getPerformance(results, level, 1.00e-7, reverse=TRUE), getPerformance(results, level, 1.00e-50, reverse=TRUE))
level = "4.5-7"
tableAdditional2 <- cbind(tableAdditional2, getPerformance(results, level, 0.1, reverse=TRUE), getPerformance(results, level, 1.00e-7, reverse=TRUE), getPerformance(results, level, 1.00e-16, reverse=TRUE))
level = "1-5.5"
tableAdditional2 <- cbind(tableAdditional2, getPerformance(results, level, 0.1, reverse=TRUE), getPerformance(results, level, 1.00e-7, reverse=TRUE), getPerformance(results, level, 1.00e-16, reverse=TRUE))
level = "2nd"
tableAdditional2 <- cbind(tableAdditional2, getPerformance(results, level, 0.1, reverse=TRUE), getPerformance(results, level, 1.00e-7, reverse=TRUE), getPerformance(results, level, 1.00e-16, reverse=TRUE))
colnames(tableAdditional2)<-c("7-10@0.1","7-10@1.00e-7","7-10@1.00e-50","4.5-7@0.1","4.5-7@1.00e-7","4.5-7@1.00e-16","1-5.5@0.1","1-5.5@1.00e-7","1-5.5@1.00e-16","2nd@0.1","2nd@1.00e-7","2nd@1.00e-16")
rownames(tableAdditional2)<-c("TP","TN","FP","FN","Precision","Recall","F1-Measure")
tableAdditional2

level = "PRR7-10"
tableAdditional2b <- cbind(         getPerformance(results, level, 1), getPerformance(results, level, 1.5), getPerformance(results, level, 5))
level = "PRR4.5-7"
tableAdditional2b <- cbind(tableAdditional2b, getPerformance(results, level, 1), getPerformance(results, level, 1.5), getPerformance(results, level, 5))
level = "PRR1-5.5"
tableAdditional2b <- cbind(tableAdditional2b, getPerformance(results, level, 1), getPerformance(results, level, 1.5), getPerformance(results, level, 5))
level = "PRR2nd"
tableAdditional2b <- cbind(tableAdditional2b, getPerformance(results, level, 1), getPerformance(results, level, 1.5), getPerformance(results, level, 5))
colnames(tableAdditional2b)<-c("7-10@1","7-10@1.5","7-10@5","4.5-7@1","4.5-7@1.5","4.5-7@5","1-5.5@1","1-5.5@1.5","1-5.5@5","2nd@1","2nd@1.5","2nd@5")
rownames(tableAdditional2b)<-c("TP","TN","FP","FN","Precision","Recall","F1-Measure")
tableAdditional2b