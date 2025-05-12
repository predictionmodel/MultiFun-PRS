#MultiFun Function
library(dplyr)
library(randomForest)
library(caret)
library(RMKL)
library(mixKernel) 
library(data.table)

#######################################################################
#1. Obtain causal probability of different omics####

# Helper function to balance positive and negative samples
balance_samples <- function(df, pos_limit = 1000, neg_ratio = 4) {
  df_pos <- df %>% filter(Y == 1)
  df_neg <- df %>% filter(Y == 0)
  
  if (nrow(df_pos) > pos_limit) {
    set.seed(20230204)
    df_pos <- slice_sample(df_pos, n = pos_limit)
  }
  #Negative samples were randomly downsampled to achieve 4:1 positive versus negative ratio.
  if (nrow(df_neg) > nrow(df_pos) * neg_ratio) {
    set.seed(20230204)
    df_neg <- slice_sample(df_neg, n = nrow(df_pos) * neg_ratio)
  }
  bind_rows(df_pos, df_neg)
}


# Obtain causal probability of multi-omics based on constructing causal SNP prediction models with summary statistics and prediction data (annotation data) of different omics.
#Y_Data：Obtain from summary statistics of specific-omics
#X_Data：Funtional annotations
#Prediction_Model：Machine learning algorithms for constructing predictive models
PriorCausalProbability<-function(Y_Data,X_Data,Prediction_Model){
  # Y_Data<-MultiOmicsSummary[[i]]
  # X_Data<-PredictionData
  # Prediction_Model<-PredictionModel
  
  Y_Data$Y<-Y_Data$Index
  #table(Y_Data[,"Y"])
  
  # Join Y and X data
  data <- merge(X_Data, as.data.frame(Y_Data[, c("POS", "Y")]), by = "POS", all.x = TRUE)
  
  # Select columns for training and testing
  testdata<-data[,-which(colnames(data) %in% c("POS","CHR","SNP","CM"))]
  traindata <- data %>% filter(!is.na(Y))
  traindata<-traindata[,-which(colnames(traindata) %in% c("POS","CHR","SNP","CM"))]
  #table(traindata$Y)
  
  # Balance the dataset
  traindata <- balance_samples(traindata)
  #table(traindata$Y)
  
  # Convert Y to factor
  traindata$Y <- as.factor(traindata$Y)
  testdata$Y <- as.factor(testdata$Y)
  
  # Train prediction model and obtain predictive probability
  PredictModel <- switch(tolower(Prediction_Model),
                         logisticregression = Logistic,
                         randomforest = RandomForest,
                         svm = SVM,
                         gradientboosing = GBM,
                         adaboost = AdaBoost)
  
  predict_result <- PredictModel(Y,traindata, testdata)
  names(predict_result)<-c("Probability","ROC")
  return(predict_result)
}

#RandomForest
RandomForest<-function(Y,traindata,testdata){
  set.seed(20230204)
  RandomForest<-randomForest(Y~.,data=traindata,importance=TRUE,na.action = na.pass)
  
  testdata2<-testdata[,-which(colnames(testdata)=="Y")]
  RandomForest_Probability<-predict(RandomForest,testdata2,type="prob")[,2]
  RandomForest_roc<-roc(testdata$Y,RandomForest_Probability)
  list<-list(RandomForest_Probability,RandomForest_roc)
  return(list)
}

#######################################################################
#2.Obtain weight of multi-omics######

#obtained weights of each omics based on the MKL algorithm
#Multi_Y_Data：The list from of summary statistics of each omics
#X_Data：Funtional annotations
MultiOmics_Weight_MKL<-function(Multi_Y_Data,X_Data){
  # Multi_Y_Data<-MultiOmicsSummary
  # X_Data<-PredictionData
  
  Traindata_MultiOmics <- list()
  ImportanceVar_MultiOmics <- list()
  # Process each Y data set
  for (i in seq_along(Multi_Y_Data)) {
    summarystatistic <- Multi_Y_Data[[i]]
    ##Y：Define whether a causal SNP is a Y,which is by summary statistics
    summarystatistic$Y <- summarystatistic$Index
    
    data <- left_join(X_Data, summarystatistic[, c("POS", "Y")], by = "POS")
    
    traindata <- data %>%
      filter(!is.na(Y)) %>%
      balance_samples()
    
    Traindata_MultiOmics[[i]] <- traindata
    traindata<-traindata[,-which(colnames(traindata) %in% c("POS","CHR","SNP","CM"))]
    ImportanceVar_MultiOmics[[i]] <- RandomForest_ImportanceVar(Y, traindata)
  }
  
  # Further sample balancing for the combined dataset
  traindata_all<-NULL
  for (i in seq_along(Multi_Y_Data)) {
    temp <- as.data.frame(Traindata_MultiOmics[[i]]) %>%
      filter(Y == 1 | Y == 0) %>%
      balance_samples(pos_limit = round(1000 / (2 * length(Multi_Y_Data))), neg_ratio = 1)
    
    traindata_all <- bind_rows(traindata_all, temp)
  }
  
  traindata_all <- distinct(traindata_all, POS, .keep_all = TRUE)
  
  # Convert Y labels from 0 to -1
  traindata_all <- mutate(traindata_all, Y = ifelse(Y == 0, -1, Y))
  
  # Compute kernel matrices and calculate weights
  SNPMatrix <- list()
  for (i in seq_along(Multi_Y_Data)) {
    temp <- traindata_all %>%
      dplyr::select(all_of(ImportanceVar_MultiOmics[[i]])) %>%
      dplyr::select_if(~!is.character(.))
    temp <- temp[, colSums(is.na(temp) | temp == 0 | temp == "")<nrow(temp)]#Delete the columns that are all na/0/""
    kernel <- compute.kernel(temp, kernel.func = "linear")
    SNPMatrix[[i]] <- kernel$kernel
    #names(SNPMatrix)[[i]] <- paste0(strsplit(names(Multi_Y_Data)[[i]], "_")[[1]][1], "_SNPMatrix")
  }
  
  C <- 1
  Weight <- SEMKL.classification(SNPMatrix, traindata_all$Y, C)$gamma
  return(Weight)
}

##Screened important variables of each omics from random forests algorithm
RandomForest_ImportanceVar<-function(Y,traindata){
  set.seed(20230204)
  traindata$Y<-as.factor(traindata$Y)
  RandomForest<-randomForest(Y~.,data=traindata,importance=TRUE,na.action = na.pass)
  Select_Number<-100 #Select top 100
  importance_train <- data.frame(importance(RandomForest), check.names = FALSE)
  importance_train <- importance_train[order(importance_train$MeanDecreaseAccuracy, decreasing = TRUE), ]
  importance_train_select <- rownames(importance_train[1:Select_Number, ])
  return(importance_train_select)
}

#######################################################################
#3. Construct MultiFun#####
MultiFunR<-function(MultiOmicsSummary,PredictionData,PredictionModel){
  # MultiOmicsSummary<-MultiOmics_summary
  # PredictionData<-Annotation
  # PredictionModel<-"randomforest"
  
  #1. Obtain causal probability of different omics
  PriorProbability<-NULL
  for(i in 1:length(MultiOmicsSummary)){
    PriorProbability[[i]]<-PriorCausalProbability(MultiOmicsSummary[[i]],PredictionData,PredictionModel)
  }
  Omics_PriorProbability<-data.frame(PredictionData[,c(1:4)])
  for(i in 1:length(MultiOmicsSummary)){
    temp<-names(MultiOmicsSummary)[i]
    Omics_PriorProbability[,temp]<-PriorProbability[[i]][["Probability"]]
  }
  AUC<-NULL
  for(i in 1:length(MultiOmicsSummary)){
    AUC[i]<-paste0(names(MultiOmicsSummary)[i],"'s AUC:",PriorProbability[[i]][["ROC"]]$auc)
  }
  print(AUC)
  
  #2.Obtain weight of multi-omics
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

