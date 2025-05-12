MultiFunR<-function(MultiOmicsSummary,PredictionData,PredictionModel){
  # MultiOmicsSummary<-MultiOmics_summary
  # PredictionData<-Annotation
  # PredictionModel<-"randomforest"
  
  #Preload the necessary packages and functions#####
  library(data.table)
  library(dplyr)
  library(randomForest)
  library(pROC)
  library(mixKernel)
  library(RMKL)
  source('/public/home/test_blank/tyd/Applicataion/Code/MultiFun_Function_20241206.R')

  #Obtain causal probability of multi-omics based on constructing causal SNP prediction models with summary statistics and prediction data (annotation data) of different omics.
  PriorProbability<-NULL
  for(i in 1:length(MultiOmicsSummary)){
    PriorProbability[[i]]<-PriorCausalProbability(MultiOmicsSummary[[i]],PredictionData,PredictionModel)
  }
  Omics_PriorProbability<-data.frame(PredictionData[,c(1:4)])
  for(i in 1:length(MultiOmicsSummary)){
    temp<-names(MultiOmicsSummary)[i]
    Omics_PriorProbability[,temp]<-PriorProbability[[i]][["Probability"]]
  }
  # AUC<-NULL
  # for(i in 1:length(MultiOmicsSummary)){
  #   AUC[i]<-paste0(names(MultiOmicsSummary)[i],"'s AUC:",PriorProbability[[i]][["ROC"]]$auc)
  # }
  # print(AUC)
  
  #Obtain weight of multi-omics
  MultiOmics_Weight<-MultiOmics_Weight_MKL(MultiOmicsSummary,PredictionData)
  
  #MultiFun-MKL:
  Omics_PriorProbability$PriorCausalProbability_MKL<-0
  for(i in 1:length(MultiOmicsSummary)){
    temp<-names(MultiOmicsSummary)[i]
    Omics_PriorProbability$PriorCausalProbability_MKL<-Omics_PriorProbability$PriorCausalProbability_MKL+
      MultiOmics_Weight[i]*Omics_PriorProbability[,temp]
  }
  
  #MultiFun-MAX:
  Omics_PriorProbability$PriorCausalProbability_MAX<-0
  Omics_PriorProbability$PriorCausalProbability_MAX<-apply(Omics_PriorProbability[,names(MultiOmicsSummary)],1,max)
  
  Omics_Result<-list(MultiOmics_Weight,Omics_PriorProbability)
  names(Omics_Result)<-c("MultiOmics_Weight","MultiOmics_PriorProbability")
  return(Omics_Result)
}

