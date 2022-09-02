#------FEATURE PRE-FILTERING AND SELECTION FOR PREDICTION OF PATHOLOGICAL STAGE IN PROSTATE CANCER------#

#----Feature Selection and Pre-Filtering techniques on Clinical Variables--#

#----Libraries----#\
library(caret)
library(pROC)
library(randomForest)
library(glmnet)
library(ROCR)

#----1. LASSO Feature Selection with Kruskal-Wallis and Chi-Squared----#
new_dat <- pros_new_dat
new_dat$PathStage = relevel(new_dat$PathStage, ref='Advanced_Stage')


R = 5
accuracy.outer = numeric(R)
accuracy.avg.inner = numeric(R)
auc.outer = numeric(R)
roc.outer = list(R)
precision = numeric(R)
recall = numeric(R)
f1.score = numeric(R)
standard_error = numeric(R)
features_selected = c()
feature_frequency_matrix = matrix(data=0,nrow = 1, ncol = ncol(pros_new_dat)-1)
colnames(feature_frequency_matrix) = colnames(pros_new_dat)[c(1:12)]
feature_stability_list = c()

sensitivities = c()
specificities = c()
thresholds = c()


set.seed(6090)
for(r in 1:R){
  n = nrow(new_dat)
  i.cv = sample(1:n, replace=FALSE)
  pc_copy = new_dat[i.cv,]
  pc.i <- createDataPartition(pc_copy$PathStage,p=0.70, list=F)
  pc_train = pc_copy[pc.i,]
  pc_test = pc_copy[-pc.i,]
  
  K = 10
  length = length(pc.i)
  folds = cut(1:length, K, labels=FALSE)
  acc.rf = numeric(K)
  
  feature.matrix = matrix(data=0,nrow = K, ncol = ncol(pc_train))
  feature.matrix[,1:10]
  colnames(feature.matrix) = colnames(pc_train)
  feature.matrix = feature.matrix[,-c(13)]
  
  set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    pros_train = pc_train[i.train,]
    pros_valid = pc_train[-i.train,]
    
    pros_train_numeric <- pros_train[,c(3:7,13)]
    pros_train_categorical <- pros_train[,c(1,2,8:12,13)]
    
    numeric_features = c()
    categorical_features = c()
    features = c()
    
    for(col in 1:(ncol(pros_train_numeric)-1)){
      kruskal_test <- kruskal.test(pros_train_numeric[,col]~pros_train_numeric$PathStage)
      if(kruskal_test$p.value<0.05){
        numeric_features = c(numeric_features,colnames(pros_train_numeric)[col])
      }
    }
    for(col in 1:(ncol(pros_train_categorical)-1)){
      pros_train_categorical[,col] = as.numeric(pros_train_categorical[,col])
      if(length(levels(as.factor(pros_train_categorical[,col])))>1){
        pros_train_categorical[,col] = as.factor(pros_train_categorical[,col])
        chisq_test <- chisq.test(pros_train_categorical[,col], pros_train_categorical$PathStage, correct = FALSE)
        if(round(chisq_test$p.value,2)<=0.05){
          categorical_features = c(categorical_features,colnames(pros_train_categorical)[col])
        }
      }
    }
    numerical_Pval = p.adjust(numerical_pval,method = 'fdr')
    for(i in 1:length(numerical_pval_adj)){
      if(numerical_pval[i]<=0.05){
        numeric_features = c(numeric_features,numerical_pval_features[i])
      }
    }
    
    categorical_Pval = p.adjust(categorical_pval,method = 'fdr')
    for(i in 1:length(categorical_pval_adj)){
      if(categorical_pval[i]<=0.05){
        categorical_features = c(categorical_features,catgorical_pval_features[i])
      }
    }
    
    features = c(numeric_features, categorical_features)
    features = c(categorical_features, numeric_features,'PathStage')
    pros_train = pros_train[,features]
    pros_valid = pros_valid[,features]
    
    x.train = pros_train[,-c(ncol(pros_train))]
    y.train = pros_train$PathStage
    x.valid = pros_valid[,-c(ncol(pros_valid))]
    y.valid = pros_valid$PathStage
    xm.train = model.matrix(PathStage~., data=pros_train)[,-1]
    xm.valid = model.matrix(PathStage~., data=pros_valid)[,-1]
    
    lam = cv.glmnet(xm.train, y.train, data=pros_train, alpha=1, family='binomial')
    lasso = glmnet(xm.train, y.train, data=pros_train, lambda=lam$lambda.min, alpha=1, family='binomial')
    lasso_features = which(coef(lasso)[-1] != 0)
    
    feature_stability_list = c(feature_stability_list, length(which(coef(lasso)[-1] != 0)))
    if(length(which(coef(lasso)[-1] != 0)) != 0){
      lasso_feature_names = colnames(xm.train)[lasso_features]
      feature.matrix[k,lasso_feature_names] = 1
      pros_train = pros_train[,c(lasso_feature_names,'PathStage')]
      pros_valid = pros_valid[,c(lasso_feature_names,'PathStage')]
      
      rf.o = randomForest(PathStage~., data=pros_train)
      rf.p = predict(rf.o, newdata=pros_valid)
      tb.rf = table(rf.p, pros_valid$PathStage)
      acc.rf[k] = sum(diag(tb.rf)) / sum(tb.rf)
      
      feature_frequency_matrix[,lasso_feature_names] = feature_frequency_matrix[,lasso_feature_names] + 1
    }
  }
  accuracy.avg.inner[r] = round(mean(acc.rf),2)
  
  fs = c()
  for(i in 1:ncol(feature.matrix)){
    if((sum(feature.matrix[,i])/nrow(feature.matrix)) >= 0.1){
      fs = c(fs, colnames(feature.matrix)[i])
    }
  }
  fs = c(fs, "PathStage")
  features_selected[r] = list(fs)
  pc_train = pc_train[,fs] 
  pc_test = pc_test[,fs]
  rf.fit = randomForest(PathStage~., data=pc_train)
  rf.pred = predict(rf.fit, newdata=pc_test)
  table.rf = confusionMatrix(rf.pred, pc_test$PathStage, mode = 'everything')
  accuracy.outer[r] = round(table.rf$overall[1],2)
  precision[r] = round(table.rf$byClass[5],2)
  recall[r] = round(table.rf$byClass[6],2)
  f1.score[r] = round(table.rf$byClass[7],2)
  rf.p = predict(rf.fit, pc_test, type="prob")[,2]
  roc.rf = roc(response=pc_test$PathStage, predictor=rf.p)
  roc.outer[[r]] = roc.rf
  auc.outer[r] = round(roc.rf$auc,2)
  p.hat = length(which(pc_test$PathStage=='Early_Stage'))/(nrow(pc_test))
  standard_error[r] = sqrt((p.hat*(1-p.hat))/nrow(pc_train))
  
  sensitivities[r] = list(roc.rf$sensitivities)
  specificities[r] = list(roc.rf$specificities)
  thresholds[r] = list(roc.rf$thresholds)
}


features_with_frequency = c()
for(i in 1:ncol(feature_frequency_matrix)){
  if(feature_frequency_matrix[,i] >= 0){
    features_with_frequency = c(features_with_frequency,colnames(feature_frequency_matrix)[i])
  }
}
KWc_L_RF_Freq_feature_frequency_matrix = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(KWc_L_RF_Freq_feature_frequency_matrix) = features_with_frequency


max_len_sensitivity = max(lengths(sensitivities))

for(i in 1:length(sensitivities)){
  if(length(sensitivities[[i]])<max_len_sensitivity){
    repeats = max_len_sensitivity-length(sensitivities[[i]])
    sensitivities[[i]] = append(sensitivities[[i]], rep(sensitivities[[i]][length(sensitivities[[i]])],times=repeats))
  }
}
max_len_specificity = max(lengths(specificities))

for(i in 1:length(specificities)){
  if(length(specificities[[i]])<max_len_specificity){
    repeats = max_len_specificity-length(specificities[[i]])
    specificities[[i]] = append(specificities[[i]], rep(specificities[[i]][length(specificities[[i]])],times=repeats))
  }
}
new_sensitivities = c()

for(i in 1:max_len_sensitivity){
  mean_sensitivity_i = (sensitivities[[1]][i]+sensitivities[[2]][i]+sensitivities[[3]][i]+sensitivities[[4]][i]+sensitivities[[5]][i])/5
  new_sensitivities = c(new_sensitivities,mean_sensitivity_i)
}
new_specificities = c()
for(i in 1:max_len_specificity){
  mean_specificity_i = (specificities[[1]][i]+specificities[[2]][i]+specificities[[3]][i]+specificities[[4]][i]+specificities[[5]][i])/5
  new_specificities = c(new_specificities,mean_specificity_i)
}

KWc_L_RF_Freq_new_sensitivities_1 = new_sensitivities
KWc_L_RF_Freq_new_specificities_1 = new_specificities

kw_chi_lasso_fs = feature_stability_list


#----2. LASSO Feature Selection with Kruskal-Wallis and One Rule----#
new_dat <- pros_new_dat
new_dat$PathStage = relevel(new_dat$PathStage, ref='Advanced_Stage')

R = 5
accuracy.outer = numeric(R)
accuracy.avg.inner = numeric(R)
auc.outer = numeric(R)
roc.outer = list(R)
precision = numeric(R)
recall = numeric(R)
f1.score = numeric(R)
standard_error = numeric(R)
features_selected = c()
feature_frequency_matrix = matrix(data=0,nrow = 1, ncol = ncol(pros_new_dat)-1)
colnames(feature_frequency_matrix) = colnames(pros_new_dat)[c(1:12)]
feature_stability_list = c()

sensitivities = c()
specificities = c()
thresholds = c()

set.seed(6090)
for(r in 1:R){
  n = nrow(new_dat)
  i.cv = sample(1:n, replace=FALSE)
  pc_copy = new_dat[i.cv,]
  pc.i <- createDataPartition(pc_copy$PathStage,p=0.70, list=F)
  pc_train = pc_copy[pc.i,]
  pc_test = pc_copy[-pc.i,]
  
  K = 10
  length = length(pc.i)
  folds = cut(1:length, K, labels=FALSE)
  acc.rf = numeric(K)
  
  feature.matrix = matrix(data=0,nrow = K, ncol = ncol(pc_train))
  feature.matrix[,1:10]
  colnames(feature.matrix) = colnames(pc_train)
  feature.matrix = feature.matrix[,-c(13)]
  
  set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    pros_train = pc_train[i.train,]
    pros_valid = pc_train[-i.train,]
    
    pros_train_numeric <- pros_train[,c(3:7,13)]
    pros_train_categorical <- pros_train[,c(1,2,8:12,13)]
    numeric_features = c()
    categorical_features = c()
    features = c()
    
    for(col in 1:(ncol(pros_train_numeric)-1)){
      kruskal_test <- kruskal.test(pros_train_numeric[,col]~pros_train_numeric$PathStage)
      if(kruskal_test$p.value<0.05){
        numeric_features = c(numeric_features,colnames(pros_train_numeric)[col])
      }
    }
    for(col in 1:(ncol(pros_train_categorical)-1)){
      oner_test <- oner_Type <- OneR::OneR(pros_train_categorical$PathStage~pros_train_categorical[,col])
      if((oner_test$correct_instances/oner_test$total_instances)>0.6){
        categorical_features = c(categorical_features,colnames(pros_train_categorical)[col])
      }
    }
    
    features = c(numeric_features, categorical_features)
    features = c(categorical_features, numeric_features,'PathStage')
    pros_train = pros_train[,features]
    pros_valid = pros_valid[,features]
    
    x.train = pros_train[,-c(ncol(pros_train))]
    y.train = pros_train$PathStage
    x.valid = pros_valid[,-c(ncol(pros_valid))]
    y.valid = pros_valid$PathStage
    xm.train = model.matrix(PathStage~., data=pros_train)[,-1]
    xm.valid = model.matrix(PathStage~., data=pros_valid)[,-1]
    
    lam = cv.glmnet(xm.train, y.train, data=pros_train, alpha=1, family='binomial')
    lasso = glmnet(xm.train, y.train, data=pros_train, lambda=lam$lambda.min, alpha=1, family='binomial')
    lasso_features = which(coef(lasso)[-1] != 0)
    
    feature_stability_list = c(feature_stability_list, length(which(coef(lasso)[-1] != 0)))
    if(length(which(coef(lasso)[-1] != 0)) != 0){
      lasso_feature_names = colnames(xm.train)[lasso_features]
      feature.matrix[k,lasso_feature_names] = 1
      pros_train = pros_train[,c(lasso_feature_names,'PathStage')]
      pros_valid = pros_valid[,c(lasso_feature_names,'PathStage')]
      
      rf.o = randomForest(PathStage~., data=pros_train)
      rf.p = predict(rf.o, newdata=pros_valid)
      tb.rf = table(rf.p, pros_valid$PathStage)
      acc.rf[k] = sum(diag(tb.rf)) / sum(tb.rf)
      
      feature_frequency_matrix[,lasso_feature_names] = feature_frequency_matrix[,lasso_feature_names] + 1
    }
  }
  accuracy.avg.inner[r] = round(mean(acc.rf),2)
  
  fs = c()
  for(i in 1:ncol(feature.matrix)){
    if((sum(feature.matrix[,i])/nrow(feature.matrix)) >= 0.1){
      fs = c(fs, colnames(feature.matrix)[i])
    }
  }
  fs = c(fs, "PathStage")
  features_selected[r] = list(fs)
  pc_train = pc_train[,fs] 
  pc_test = pc_test[,fs]
  rf.fit = randomForest(PathStage~., data=pc_train)
  rf.pred = predict(rf.fit, newdata=pc_test)
  table.rf = confusionMatrix(rf.pred, pc_test$PathStage, mode = 'everything')
  accuracy.outer[r] = round(table.rf$overall[1],2)
  precision[r] = round(table.rf$byClass[5],2)
  recall[r] = round(table.rf$byClass[6],2)
  f1.score[r] = round(table.rf$byClass[7],2)
  rf.p = predict(rf.fit, pc_test, type="prob")[,2]
  roc.rf = roc(response=pc_test$PathStage, predictor=rf.p)
  roc.outer[[r]] = roc.rf
  auc.outer[r] = round(roc.rf$auc,2)
  p.hat = length(which(pc_test$PathStage=='Early_Stage'))/(nrow(pc_test))
  standard_error[r] = sqrt((p.hat*(1-p.hat))/nrow(pc_train))
  
  sensitivities[r] = list(roc.rf$sensitivities)
  specificities[r] = list(roc.rf$specificities)
  thresholds[r] = list(roc.rf$thresholds)
}

features_with_frequency = c()
for(i in 1:ncol(feature_frequency_matrix)){
  if(feature_frequency_matrix[,i] > 0){
    features_with_frequency = c(features_with_frequency,colnames(feature_frequency_matrix)[i])
  }
}

KWo_L_RF_Freq_feature_frequency_matrix = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(KWo_L_RF_Freq_feature_frequency_matrix) = features_with_frequency


max_len_sensitivity = max(lengths(sensitivities))

for(i in 1:length(sensitivities)){
  if(length(sensitivities[[i]])<max_len_sensitivity){
    repeats = max_len_sensitivity-length(sensitivities[[i]])
    sensitivities[[i]] = append(sensitivities[[i]], rep(sensitivities[[i]][length(sensitivities[[i]])],times=repeats))
  }
}

max_len_specificity = max(lengths(specificities))

for(i in 1:length(specificities)){
  if(length(specificities[[i]])<max_len_specificity){
    repeats = max_len_specificity-length(specificities[[i]])
    specificities[[i]] = append(specificities[[i]], rep(specificities[[i]][length(specificities[[i]])],times=repeats))
  }
}
new_sensitivities = c()
for(i in 1:max_len_sensitivity){
  mean_sensitivity_i = (sensitivities[[1]][i]+sensitivities[[2]][i]+sensitivities[[3]][i]+sensitivities[[4]][i]+sensitivities[[5]][i])/5
  new_sensitivities = c(new_sensitivities,mean_sensitivity_i)
}
new_specificities = c()

for(i in 1:max_len_specificity){
  mean_specificity_i = (specificities[[1]][i]+specificities[[2]][i]+specificities[[3]][i]+specificities[[4]][i]+specificities[[5]][i])/5
  new_specificities = c(new_specificities,mean_specificity_i)
}

KWo_L_RF_Freq_new_sensitivities_2 = new_sensitivities
KWo_L_RF_Freq_new_specificities_2 = new_specificities

kw_oneR_lasso_fs = feature_stability_list


#----3. LASSO Feature Selection with Kruskal-Wallis and Information Gain----#

new_dat <- pros_new_dat
new_dat$PathStage = relevel(new_dat$PathStage, ref='Advanced_Stage')

R = 5
accuracy.outer = numeric(R)
accuracy.avg.inner = numeric(R)
auc.outer = numeric(R)
roc.outer = list(R)
precision = numeric(R)
recall = numeric(R)
f1.score = numeric(R)
standard_error = numeric(R)
features_selected = c()
feature_frequency_matrix = matrix(data=0,nrow = 1, ncol = ncol(pros_new_dat)-1)
colnames(feature_frequency_matrix) = colnames(pros_new_dat)[c(1:12)]
feature_stability_list = c()

sensitivities = c()
specificities = c()
thresholds = c()


set.seed(6090)
for(r in 1:R){
  n = nrow(new_dat)
  i.cv = sample(1:n, replace=FALSE)
  pc_copy = new_dat[i.cv,]
  pc.i <- createDataPartition(pc_copy$PathStage,p=0.70, list=F)
  pc_train = pc_copy[pc.i,]
  pc_test = pc_copy[-pc.i,]
  
  K = 10
  length = length(pc.i)
  folds = cut(1:length, K, labels=FALSE)
  acc.rf = numeric(K)
  
  feature.matrix = matrix(data=0,nrow = K, ncol = ncol(pc_train))
  feature.matrix[,1:10]
  colnames(feature.matrix) = colnames(pc_train)
  feature.matrix = feature.matrix[,-c(13)]
  
  set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    pros_train = pc_train[i.train,]
    pros_valid = pc_train[-i.train,]
    
    pros_train_numeric <- pros_train[,c(3:7,13)]
    pros_train_categorical <- pros_train[,c(1,2,8:12,13)]
    numeric_features = c()
    categorical_features = c()
    features = c()
    
    for(col in 1:(ncol(pros_train_numeric)-1)){
      kruskal_test <- kruskal.test(pros_train_numeric[,col]~pros_train_numeric$PathStage)
      if(kruskal_test$p.value<0.05){
        numeric_features = c(numeric_features,colnames(pros_train_numeric)[col])
      }
    }
    for(col in 1:(ncol(pros_train_categorical)-1)){
      pros_train_categorical[,col] = as.factor(pros_train_categorical[,col])
      pros_train_categorical_temp = pros_train_categorical[,c(col,ncol(pros_train_categorical))]
      infogain_test <- information_gain(PathStage~., data=pros_train_categorical_temp, type="infogain")
      if(round(infogain_test$importance,2)>0){
        categorical_features = c(categorical_features,colnames(pros_train_categorical)[col])
      }
      pros_train_categorical[,col] = as.numeric(pros_train_categorical[,col])
    }
    
    features = c(numeric_features, categorical_features)
    features = c(categorical_features, numeric_features,'PathStage')
    pros_train = pros_train[,features]
    pros_valid = pros_valid[,features]
    
    x.train = pros_train[,-c(ncol(pros_train))]
    y.train = pros_train$PathStage
    x.valid = pros_valid[,-c(ncol(pros_valid))]
    y.valid = pros_valid$PathStage
    xm.train = model.matrix(PathStage~., data=pros_train)[,-1]
    xm.valid = model.matrix(PathStage~., data=pros_valid)[,-1]
    
    lam = cv.glmnet(xm.train, y.train, data=pros_train, alpha=1, family='binomial')
    lasso = glmnet(xm.train, y.train, data=pros_train, lambda=lam$lambda.min, alpha=1, family='binomial')
    lasso_features = which(coef(lasso)[-1] != 0)
    
    feature_stability_list = c(feature_stability_list, length(which(coef(lasso)[-1] != 0)))
    if(length(which(coef(lasso)[-1] != 0)) != 0){
      lasso_feature_names = colnames(xm.train)[lasso_features]
      feature.matrix[k,lasso_feature_names] = 1
      pros_train = pros_train[,c(lasso_feature_names,'PathStage')]
      pros_valid = pros_valid[,c(lasso_feature_names,'PathStage')]
      
      rf.o = randomForest(PathStage~., data=pros_train)
      rf.p = predict(rf.o, newdata=pros_valid)
      tb.rf = table(rf.p, pros_valid$PathStage)
      acc.rf[k] = sum(diag(tb.rf)) / sum(tb.rf)
      
      feature_frequency_matrix[,lasso_feature_names] = feature_frequency_matrix[,lasso_feature_names] + 1
    }
  }
  accuracy.avg.inner[r] = round(mean(acc.rf),2)
  
  fs = c()
  for(i in 1:ncol(feature.matrix)){
    if((sum(feature.matrix[,i])/nrow(feature.matrix)) >= 0.1){
      fs = c(fs, colnames(feature.matrix)[i])
    }
  }
  fs = c(fs, "PathStage")
  features_selected[r] = list(fs)
  pc_train = pc_train[,fs] 
  pc_test = pc_test[,fs]
  rf.fit = randomForest(PathStage~., data=pc_train)
  rf.pred = predict(rf.fit, newdata=pc_test)
  table.rf = confusionMatrix(rf.pred, pc_test$PathStage, mode = 'everything')
  accuracy.outer[r] = round(table.rf$overall[1],2)
  precision[r] = round(table.rf$byClass[5],2)
  recall[r] = round(table.rf$byClass[6],2)
  f1.score[r] = round(table.rf$byClass[7],2)
  rf.p = predict(rf.fit, pc_test, type="prob")[,2]
  roc.rf = roc(response=pc_test$PathStage, predictor=rf.p)
  roc.outer[[r]] = roc.rf
  auc.outer[r] = round(roc.rf$auc,2)
  p.hat = length(which(pc_test$PathStage=='Early_Stage'))/(nrow(pc_test))
  standard_error[r] = sqrt((p.hat*(1-p.hat))/nrow(pc_train))
  
  sensitivities[r] = list(roc.rf$sensitivities)
  specificities[r] = list(roc.rf$specificities)
  thresholds[r] = list(roc.rf$thresholds)
}

features_with_frequency = c()
for(i in 1:ncol(feature_frequency_matrix)){
  if(feature_frequency_matrix[,i] > 0){
    features_with_frequency = c(features_with_frequency,colnames(feature_frequency_matrix)[i])
  }
}

kwi_L_RF_Freq_feature_frequency_matrix = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(kwi_L_RF_Freq_feature_frequency_matrix) = features_with_frequency


max_len_sensitivity = max(lengths(sensitivities))

for(i in 1:length(sensitivities)){
  if(length(sensitivities[[i]])<max_len_sensitivity){
    repeats = max_len_sensitivity-length(sensitivities[[i]])
    sensitivities[[i]] = append(sensitivities[[i]], rep(sensitivities[[i]][length(sensitivities[[i]])],times=repeats))
  }
}

max_len_specificity = max(lengths(specificities))

for(i in 1:length(specificities)){
  if(length(specificities[[i]])<max_len_specificity){
    repeats = max_len_specificity-length(specificities[[i]])
    specificities[[i]] = append(specificities[[i]], rep(specificities[[i]][length(specificities[[i]])],times=repeats))
  }
}


new_sensitivities = c()

for(i in 1:max_len_sensitivity){
  mean_sensitivity_i = (sensitivities[[1]][i]+sensitivities[[2]][i]+sensitivities[[3]][i]+sensitivities[[4]][i]+sensitivities[[5]][i])/5
  new_sensitivities = c(new_sensitivities,mean_sensitivity_i)
}
new_sensitivities

new_specificities = c()

for(i in 1:max_len_specificity){
  mean_specificity_i = (specificities[[1]][i]+specificities[[2]][i]+specificities[[3]][i]+specificities[[4]][i]+specificities[[5]][i])/5
  new_specificities = c(new_specificities,mean_specificity_i)
}
new_specificities

kwi_L_RF_Freq_new_sensitivities_3 = new_sensitivities
kwi_L_RF_Freq_new_specificities_3 = new_specificities


kw_infogain_lasso_fs = feature_stability_list


#----4. LASSO Feature Selection with Kruskal-Wallis and Symmetrical Uncertainty----#

new_dat <- pros_new_dat
new_dat$PathStage = relevel(new_dat$PathStage, ref='Advanced_Stage')
levels(new_dat$PathStage)
nrow(new_dat)
ncol(new_dat)
#--correlation features--#
library(caret)
library(pROC)
library(randomForest)
library(glmnet)
library(ROCR)

R = 5
accuracy.outer = numeric(R)
accuracy.avg.inner = numeric(R)
auc.outer = numeric(R)
roc.outer = list(R)
precision = numeric(R)
recall = numeric(R)
f1.score = numeric(R)
standard_error = numeric(R)
features_selected = c()
feature_frequency_matrix = matrix(data=0,nrow = 1, ncol = ncol(pros_new_dat)-1)
colnames(feature_frequency_matrix) = colnames(pros_new_dat)[c(1:12)]
feature_stability_list = c()

sensitivities = c()
specificities = c()
thresholds = c()


set.seed(6090)
for(r in 1:R){
  n = nrow(new_dat)
  i.cv = sample(1:n, replace=FALSE)
  pc_copy = new_dat[i.cv,]
  pc.i <- createDataPartition(pc_copy$PathStage,p=0.70, list=F)
  pc_train = pc_copy[pc.i,]
  pc_test = pc_copy[-pc.i,]
  
  K = 10
  length = length(pc.i)
  folds = cut(1:length, K, labels=FALSE)
  acc.rf = numeric(K)
  
  feature.matrix = matrix(data=0,nrow = K, ncol = ncol(pc_train))
  feature.matrix[,1:10]
  colnames(feature.matrix) = colnames(pc_train)
  feature.matrix = feature.matrix[,-c(13)]
  
  set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    pros_train = pc_train[i.train,]
    pros_valid = pc_train[-i.train,]
    
    pros_train_numeric <- pros_train[,c(3:7,13)]
    pros_train_categorical <- pros_train[,c(1,2,8:12,13)]
    numeric_features = c()
    categorical_features = c()
    features = c()
    
    for(col in 1:(ncol(pros_train_numeric)-1)){
      kruskal_test <- kruskal.test(pros_train_numeric[,col]~pros_train_numeric$PathStage)
      if(kruskal_test$p.value<0.05){
        numeric_features = c(numeric_features,colnames(pros_train_numeric)[col])
      }
    }
    for(col in 1:(ncol(pros_train_categorical)-1)){
      pros_train_categorical[,col] = as.factor(pros_train_categorical[,col])
      pros_train_categorical_temp = pros_train_categorical[,c(col,ncol(pros_train_categorical))]
      infogain_test <- information_gain(PathStage~., data=pros_train_categorical_temp, type="symuncert")
      if(round(infogain_test$importance,2)>0){
        categorical_features = c(categorical_features,colnames(pros_train_categorical)[col])
      }
      pros_train_categorical[,col] = as.numeric(pros_train_categorical[,col])
    }
    
    features = c(numeric_features, categorical_features)
    features = c(categorical_features, numeric_features,'PathStage')
    pros_train = pros_train[,features]
    pros_valid = pros_valid[,features]
    
    x.train = pros_train[,-c(ncol(pros_train))]
    y.train = pros_train$PathStage
    x.valid = pros_valid[,-c(ncol(pros_valid))]
    y.valid = pros_valid$PathStage
    xm.train = model.matrix(PathStage~., data=pros_train)[,-1]
    xm.valid = model.matrix(PathStage~., data=pros_valid)[,-1]
    
    lam = cv.glmnet(xm.train, y.train, data=pros_train, alpha=1, family='binomial')
    lasso = glmnet(xm.train, y.train, data=pros_train, lambda=lam$lambda.min, alpha=1, family='binomial')
    lasso_features = which(coef(lasso)[-1] != 0)
    
    feature_stability_list = c(feature_stability_list, length(which(coef(lasso)[-1] != 0)))
    if(length(which(coef(lasso)[-1] != 0)) != 0){
      lasso_feature_names = colnames(xm.train)[lasso_features]
      feature.matrix[k,lasso_feature_names] = 1
      pros_train = pros_train[,c(lasso_feature_names,'PathStage')]
      pros_valid = pros_valid[,c(lasso_feature_names,'PathStage')]
      
      rf.o = randomForest(PathStage~., data=pros_train)
      rf.p = predict(rf.o, newdata=pros_valid)
      tb.rf = table(rf.p, pros_valid$PathStage)
      acc.rf[k] = sum(diag(tb.rf)) / sum(tb.rf)
      
      feature_frequency_matrix[,lasso_feature_names] = feature_frequency_matrix[,lasso_feature_names] + 1
    }
  }
  accuracy.avg.inner[r] = round(mean(acc.rf),2)
  
  fs = c()
  for(i in 1:ncol(feature.matrix)){
    if((sum(feature.matrix[,i])/nrow(feature.matrix)) >= 0.1){
      fs = c(fs, colnames(feature.matrix)[i])
    }
  }
  fs = c(fs, "PathStage")
  features_selected[r] = list(fs)
  pc_train = pc_train[,fs] 
  pc_test = pc_test[,fs]
  rf.fit = randomForest(PathStage~., data=pc_train)
  rf.pred = predict(rf.fit, newdata=pc_test)
  table.rf = confusionMatrix(rf.pred, pc_test$PathStage, mode = 'everything')
  accuracy.outer[r] = round(table.rf$overall[1],2)
  precision[r] = round(table.rf$byClass[5],2)
  recall[r] = round(table.rf$byClass[6],2)
  f1.score[r] = round(table.rf$byClass[7],2)
  rf.p = predict(rf.fit, pc_test, type="prob")[,2]
  roc.rf = roc(response=pc_test$PathStage, predictor=rf.p)
  roc.outer[[r]] = roc.rf
  auc.outer[r] = round(roc.rf$auc,2)
  p.hat = length(which(pc_test$PathStage=='Early_Stage'))/(nrow(pc_test))
  standard_error[r] = sqrt((p.hat*(1-p.hat))/nrow(pc_train))
  
  sensitivities[r] = list(roc.rf$sensitivities)
  specificities[r] = list(roc.rf$specificities)
  thresholds[r] = list(roc.rf$thresholds)
}


features_with_frequency = c()
for(i in 1:ncol(feature_frequency_matrix)){
  if(feature_frequency_matrix[,i] > 0){
    features_with_frequency = c(features_with_frequency,colnames(feature_frequency_matrix)[i])
  }
}

kws_L_RF_Freq_feature_frequency_matrix= matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(kws_L_RF_Freq_feature_frequency_matrix) = features_with_frequency

accuracy.outer
round(mean(accuracy.outer),2)
round(sd(accuracy.outer),2)

accuracy.avg.inner
round(mean(accuracy.avg.inner),2)

precision
round(mean(precision),2)
round(sd(precision),2)

recall
round(mean(recall),2)
round(sd(recall),2)

f1.score
round(mean(f1.score),2)
round(sd(f1.score),2)

auc.outer
round(mean(auc.outer),2)
round(sd(auc.outer),2)


boxplot(feature_stability_list)

par(mfrow=c(1,2))
plot(roc.outer[[1]], col='black')
plot(roc.outer[[2]], add=TRUE, col='red')
plot(roc.outer[[3]], add=TRUE, col='skyblue')
plot(roc.outer[[4]], add=TRUE, col='green')
plot(roc.outer[[5]], add=TRUE, col='orange')
legend(x = 'bottomright', legend=c("Iter 1", "Iter 2", "Iter 3","Iter 4","Iter 5"), 
       col = c("black","red",'skyblue', 'green','orange'), lty=1, cex=1.0, title = 'Kruskal-Wallis - Chi-Sq'
)

max_len_sensitivity = max(lengths(sensitivities))

for(i in 1:length(sensitivities)){
  if(length(sensitivities[[i]])<max_len_sensitivity){
    repeats = max_len_sensitivity-length(sensitivities[[i]])
    sensitivities[[i]] = append(sensitivities[[i]], rep(sensitivities[[i]][length(sensitivities[[i]])],times=repeats))
  }
}

max_len_specificity = max(lengths(specificities))

for(i in 1:length(specificities)){
  if(length(specificities[[i]])<max_len_specificity){
    repeats = max_len_specificity-length(specificities[[i]])
    specificities[[i]] = append(specificities[[i]], rep(specificities[[i]][length(specificities[[i]])],times=repeats))
  }
}


new_sensitivities = c()

for(i in 1:max_len_sensitivity){
  mean_sensitivity_i = (sensitivities[[1]][i]+sensitivities[[2]][i]+sensitivities[[3]][i]+sensitivities[[4]][i]+sensitivities[[5]][i])/5
  new_sensitivities = c(new_sensitivities,mean_sensitivity_i)
}
new_sensitivities

new_specificities = c()

for(i in 1:max_len_specificity){
  mean_specificity_i = (specificities[[1]][i]+specificities[[2]][i]+specificities[[3]][i]+specificities[[4]][i]+specificities[[5]][i])/5
  new_specificities = c(new_specificities,mean_specificity_i)
}
new_specificities

kws_L_RF_Freq_new_sensitivities_4 = new_sensitivities
kws_L_RF_Freq_new_specificities_4 = new_specificities

kw_symun_lasso_fs = feature_stability_list

#----5. Forward Stepwise Feature Selection with Kruskal-Wallis and Chi-Squared----#

new_dat <- pros_new_dat
new_dat$PathStage = relevel(new_dat$PathStage, ref='Advanced_Stage')

R = 5
accuracy.outer = numeric(R)
accuracy.avg.inner = numeric(R)
auc.outer = numeric(R)
roc.outer = list(R)
precision = numeric(R)
recall = numeric(R)
f1.score = numeric(R)
standard_error = numeric(R)
features_selected = c()
feature_frequency_matrix = matrix(data=0,nrow = 1, ncol = ncol(pros_new_dat)-1)
colnames(feature_frequency_matrix) = colnames(pros_new_dat)[c(1:12)]
feature_stability_list = c()

sensitivities = c()
specificities = c()
thresholds = c()


set.seed(6090)
for(r in 1:R){
  n = nrow(new_dat)
  i.cv = sample(1:n, replace=FALSE)
  pc_copy = new_dat[i.cv,]
  pc.i <- createDataPartition(pc_copy$PathStage,p=0.70, list=F)
  pc_train = pc_copy[pc.i,]
  pc_test = pc_copy[-pc.i,]
  
  K = 10
  length = length(pc.i)
  folds = cut(1:length, K, labels=FALSE)
  acc.rf = numeric(K)
  
  feature.matrix = matrix(data=0,nrow = K, ncol = ncol(pc_train))
  feature.matrix[,1:10]
  colnames(feature.matrix) = colnames(pc_train)
  feature.matrix = feature.matrix[,-c(13)]
  
  set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    pros_train = pc_train[i.train,]
    pros_valid = pc_train[-i.train,]
    
    pros_train_numeric <- pros_train[,c(3:7,13)]
    pros_train_categorical <- pros_train[,c(1,2,8:12,13)]
    numeric_features = c()
    categorical_features = c()
    features = c()
    
    for(col in 1:(ncol(pros_train_numeric)-1)){
      kruskal_test <- kruskal.test(pros_train_numeric[,col]~pros_train_numeric$PathStage)
      if(kruskal_test$p.value<0.05){
        numeric_features = c(numeric_features,colnames(pros_train_numeric)[col])
      }
    }
    for(col in 1:(ncol(pros_train_categorical)-1)){
      pros_train_categorical[,col] = as.numeric(pros_train_categorical[,col])
      if(length(levels(as.factor(pros_train_categorical[,col])))>1){
        pros_train_categorical[,col] = as.factor(pros_train_categorical[,col])
        chisq_test <- chisq.test(pros_train_categorical[,col], pros_train_categorical$PathStage, correct = FALSE)
        if(round(chisq_test$p.value,2)<=0.05){
          categorical_features = c(categorical_features,colnames(pros_train_categorical)[col])
        }
      }
    }
    
    features = c(numeric_features, categorical_features)
    features = c(categorical_features, numeric_features,'PathStage')
    pros_train = pros_train[,features]
    pros_valid = pros_valid[,features]
    
    reg.fwd = regsubsets(PathStage~., data=pros_train, method="forward")
    R2adj = summary(reg.fwd)$adjr2
    R2adj.index = which.max(R2adj)
    mod = summary(reg.fwd)$which[R2adj.index,]
    vars = rownames(data.frame(which(mod)[-1]))
    feature.matrix[k,vars] = 1
    pros_train = pros_train[,c(vars,'PathStage')]
    pros_valid = pros_valid[,c(vars,'PathStage')]
    
    rf.o = randomForest(PathStage~., data=pros_train)
    rf.p = predict(rf.o, newdata=pros_valid)
    tb.rf = table(rf.p, pros_valid$PathStage)
    acc.rf[k] = sum(diag(tb.rf)) / sum(tb.rf)
    
    feature_frequency_matrix[,vars] = feature_frequency_matrix[,vars] + 1
    feature_stability_list = c(feature_stability_list, length(vars))
    
  }
  accuracy.avg.inner[r] = round(mean(acc.rf),2)
  
  fs = c()
  for(i in 1:ncol(feature.matrix)){
    if((sum(feature.matrix[,i])/nrow(feature.matrix)) >= 0.1){
      fs = c(fs, colnames(feature.matrix)[i])
    }
  }
  fs = c(fs, "PathStage")
  features_selected[r] = list(fs)
  pc_train = pc_train[,fs] 
  pc_test = pc_test[,fs]
  rf.fit = randomForest(PathStage~., data=pc_train)
  rf.pred = predict(rf.fit, newdata=pc_test)
  table.rf = confusionMatrix(rf.pred, pc_test$PathStage, mode = 'everything')
  accuracy.outer[r] = round(table.rf$overall[1],2)
  precision[r] = round(table.rf$byClass[5],2)
  recall[r] = round(table.rf$byClass[6],2)
  f1.score[r] = round(table.rf$byClass[7],2)
  rf.p = predict(rf.fit, pc_test, type="prob")[,2]
  roc.rf = roc(response=pc_test$PathStage, predictor=rf.p)
  roc.outer[[r]] = roc.rf
  auc.outer[r] = round(roc.rf$auc,2)
  p.hat = length(which(pc_test$PathStage=='Early_Stage'))/(nrow(pc_test))
  standard_error[r] = sqrt((p.hat*(1-p.hat))/nrow(pc_train))
  
  sensitivities[r] = list(roc.rf$sensitivities)
  specificities[r] = list(roc.rf$specificities)
  thresholds[r] = list(roc.rf$thresholds)
}

features_with_frequency = c()
for(i in 1:ncol(feature_frequency_matrix)){
  if(feature_frequency_matrix[,i] > 0){
    features_with_frequency = c(features_with_frequency,colnames(feature_frequency_matrix)[i])
  }
}

kC_F_RF_Freq_feature_frequency_matrix = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(kC_F_RF_Freq_feature_frequency_matrix) = features_with_frequency



max_len_sensitivity = max(lengths(sensitivities))

for(i in 1:length(sensitivities)){
  if(length(sensitivities[[i]])<max_len_sensitivity){
    repeats = max_len_sensitivity-length(sensitivities[[i]])
    sensitivities[[i]] = append(sensitivities[[i]], rep(sensitivities[[i]][length(sensitivities[[i]])],times=repeats))
  }
}

max_len_specificity = max(lengths(specificities))

for(i in 1:length(specificities)){
  if(length(specificities[[i]])<max_len_specificity){
    repeats = max_len_specificity-length(specificities[[i]])
    specificities[[i]] = append(specificities[[i]], rep(specificities[[i]][length(specificities[[i]])],times=repeats))
  }
}


new_sensitivities = c()

for(i in 1:max_len_sensitivity){
  mean_sensitivity_i = (sensitivities[[1]][i]+sensitivities[[2]][i]+sensitivities[[3]][i]+sensitivities[[4]][i]+sensitivities[[5]][i])/5
  new_sensitivities = c(new_sensitivities,mean_sensitivity_i)
}
new_sensitivities

new_specificities = c()

for(i in 1:max_len_specificity){
  mean_specificity_i = (specificities[[1]][i]+specificities[[2]][i]+specificities[[3]][i]+specificities[[4]][i]+specificities[[5]][i])/5
  new_specificities = c(new_specificities,mean_specificity_i)
}
new_specificities

kC_F_RF_Freq_new_sensitivities_1 = new_sensitivities
kC_F_RF_Freq_new_specificities_1 = new_specificities

kw_chi_forward_fs = feature_stability_list

#----6. Forward Stepwise Feature Selection with Kruskal-Wallis and One Rule----#
new_dat <- pros_new_dat
new_dat$PathStage = relevel(new_dat$PathStage, ref='Advanced_Stage')

R = 5
accuracy.outer = numeric(R)
accuracy.avg.inner = numeric(R)
auc.outer = numeric(R)
roc.outer = list(R)
precision = numeric(R)
recall = numeric(R)
f1.score = numeric(R)
standard_error = numeric(R)
features_selected = c()
feature_frequency_matrix = matrix(data=0,nrow = 1, ncol = ncol(pros_new_dat)-1)
colnames(feature_frequency_matrix) = colnames(pros_new_dat)[c(1:12)]
feature_stability_list = c()

sensitivities = c()
specificities = c()
thresholds = c()


set.seed(6090)
for(r in 1:R){
  n = nrow(new_dat)
  i.cv = sample(1:n, replace=FALSE)
  pc_copy = new_dat[i.cv,]
  pc.i <- createDataPartition(pc_copy$PathStage,p=0.70, list=F)
  pc_train = pc_copy[pc.i,]
  pc_test = pc_copy[-pc.i,]
  
  K = 10
  length = length(pc.i)
  folds = cut(1:length, K, labels=FALSE)
  acc.rf = numeric(K)
  
  feature.matrix = matrix(data=0,nrow = K, ncol = ncol(pc_train))
  feature.matrix[,1:10]
  colnames(feature.matrix) = colnames(pc_train)
  feature.matrix = feature.matrix[,-c(13)]
  
  set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    pros_train = pc_train[i.train,]
    pros_valid = pc_train[-i.train,]
    
    pros_train_numeric <- pros_train[,c(3:7,13)]
    pros_train_categorical <- pros_train[,c(1,2,8:12,13)]
    numeric_features = c()
    categorical_features = c()
    features = c()
    
    for(col in 1:(ncol(pros_train_numeric)-1)){
      kruskal_test <- kruskal.test(pros_train_numeric[,col]~pros_train_numeric$PathStage)
      if(kruskal_test$p.value<0.05){
        numeric_features = c(numeric_features,colnames(pros_train_numeric)[col])
      }
    }
    for(col in 1:(ncol(pros_train_categorical)-1)){
      oner_test <- oner_Type <- OneR::OneR(pros_train_categorical$PathStage~pros_train_categorical[,col])
      if((oner_test$correct_instances/oner_test$total_instances)>=0.6){
        categorical_features = c(categorical_features,colnames(pros_train_categorical)[col])
      }
    }
    
    features = c(numeric_features, categorical_features)
    features = c(categorical_features, numeric_features,'PathStage')
    pros_train = pros_train[,features]
    pros_valid = pros_valid[,features]
    
    reg.fwd = regsubsets(PathStage~., data=pros_train, method="forward")
    R2adj = summary(reg.fwd)$adjr2
    R2adj.index = which.max(R2adj)
    mod = summary(reg.fwd)$which[R2adj.index,]
    vars = rownames(data.frame(which(mod)[-1]))
    feature.matrix[k,vars] = 1
    pros_train = pros_train[,c(vars,'PathStage')]
    pros_valid = pros_valid[,c(vars,'PathStage')]
    
    rf.o = randomForest(PathStage~., data=pros_train)
    rf.p = predict(rf.o, newdata=pros_valid)
    tb.rf = table(rf.p, pros_valid$PathStage)
    acc.rf[k] = sum(diag(tb.rf)) / sum(tb.rf)
    
    feature_frequency_matrix[,vars] = feature_frequency_matrix[,vars] + 1
    feature_stability_list = c(feature_stability_list, length(vars))
    
  }
  accuracy.avg.inner[r] = round(mean(acc.rf),2)
  
  fs = c()
  for(i in 1:ncol(feature.matrix)){
    if((sum(feature.matrix[,i])/nrow(feature.matrix)) >= 0.1){
      fs = c(fs, colnames(feature.matrix)[i])
    }
  }
  fs = c(fs, "PathStage")
  features_selected[r] = list(fs)
  pc_train = pc_train[,fs] 
  pc_test = pc_test[,fs]
  rf.fit = randomForest(PathStage~., data=pc_train)
  rf.pred = predict(rf.fit, newdata=pc_test)
  table.rf = confusionMatrix(rf.pred, pc_test$PathStage, mode = 'everything')
  accuracy.outer[r] = round(table.rf$overall[1],2)
  precision[r] = round(table.rf$byClass[5],2)
  recall[r] = round(table.rf$byClass[6],2)
  f1.score[r] = round(table.rf$byClass[7],2)
  rf.p = predict(rf.fit, pc_test, type="prob")[,2]
  roc.rf = roc(response=pc_test$PathStage, predictor=rf.p)
  roc.outer[[r]] = roc.rf
  auc.outer[r] = round(roc.rf$auc,2)
  p.hat = length(which(pc_test$PathStage=='Early_Stage'))/(nrow(pc_test))
  standard_error[r] = sqrt((p.hat*(1-p.hat))/nrow(pc_train))
  
  sensitivities[r] = list(roc.rf$sensitivities)
  specificities[r] = list(roc.rf$specificities)
  thresholds[r] = list(roc.rf$thresholds)
}

features_with_frequency = c()
for(i in 1:ncol(feature_frequency_matrix)){
  if(feature_frequency_matrix[,i] > 0){
    features_with_frequency = c(features_with_frequency,colnames(feature_frequency_matrix)[i])
  }
}

ko_F_RF_Freq_feature_frequency_matrix = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(ko_F_RF_Freq_feature_frequency_matrix) = features_with_frequency





max_len_sensitivity = max(lengths(sensitivities))

for(i in 1:length(sensitivities)){
  if(length(sensitivities[[i]])<max_len_sensitivity){
    repeats = max_len_sensitivity-length(sensitivities[[i]])
    sensitivities[[i]] = append(sensitivities[[i]], rep(sensitivities[[i]][length(sensitivities[[i]])],times=repeats))
  }
}

max_len_specificity = max(lengths(specificities))

for(i in 1:length(specificities)){
  if(length(specificities[[i]])<max_len_specificity){
    repeats = max_len_specificity-length(specificities[[i]])
    specificities[[i]] = append(specificities[[i]], rep(specificities[[i]][length(specificities[[i]])],times=repeats))
  }
}


new_sensitivities = c()

for(i in 1:max_len_sensitivity){
  mean_sensitivity_i = (sensitivities[[1]][i]+sensitivities[[2]][i]+sensitivities[[3]][i]+sensitivities[[4]][i]+sensitivities[[5]][i])/5
  new_sensitivities = c(new_sensitivities,mean_sensitivity_i)
}
new_sensitivities

new_specificities = c()

for(i in 1:max_len_specificity){
  mean_specificity_i = (specificities[[1]][i]+specificities[[2]][i]+specificities[[3]][i]+specificities[[4]][i]+specificities[[5]][i])/5
  new_specificities = c(new_specificities,mean_specificity_i)
}
new_specificities

ko_F_RF_Freq_new_sensitivities_2 = new_sensitivities
ko_F_RF_Freq_new_specificities_2 = new_specificities



kw_oneR_forward_fs = feature_stability_list

#----7. Forward Stepwise Feature Selection with Kruskal-Wallis and Information Gain----#
new_dat <- pros_new_dat
new_dat$PathStage = relevel(new_dat$PathStage, ref='Advanced_Stage')

R = 5
accuracy.outer = numeric(R)
accuracy.avg.inner = numeric(R)
auc.outer = numeric(R)
roc.outer = list(R)
precision = numeric(R)
recall = numeric(R)
f1.score = numeric(R)
standard_error = numeric(R)
features_selected = c()
feature_frequency_matrix = matrix(data=0,nrow = 1, ncol = ncol(pros_new_dat)-1)
colnames(feature_frequency_matrix) = colnames(pros_new_dat)[c(1:12)]
feature_stability_list = c()

sensitivities = c()
specificities = c()
thresholds = c()


set.seed(6090)
for(r in 1:R){
  n = nrow(new_dat)
  i.cv = sample(1:n, replace=FALSE)
  pc_copy = new_dat[i.cv,]
  pc.i <- createDataPartition(pc_copy$PathStage,p=0.70, list=F)
  pc_train = pc_copy[pc.i,]
  pc_test = pc_copy[-pc.i,]
  
  K = 10
  length = length(pc.i)
  folds = cut(1:length, K, labels=FALSE)
  acc.rf = numeric(K)
  
  feature.matrix = matrix(data=0,nrow = K, ncol = ncol(pc_train))
  feature.matrix[,1:10]
  colnames(feature.matrix) = colnames(pc_train)
  feature.matrix = feature.matrix[,-c(13)]
  
  set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    pros_train = pc_train[i.train,]
    pros_valid = pc_train[-i.train,]
    
    pros_train_numeric <- pros_train[,c(3:7,13)]
    pros_train_categorical <- pros_train[,c(1,2,8:12,13)]
    numeric_features = c()
    categorical_features = c()
    features = c()
    
    for(col in 1:(ncol(pros_train_numeric)-1)){
      kruskal_test <- kruskal.test(pros_train_numeric[,col]~pros_train_numeric$PathStage)
      if(kruskal_test$p.value<0.05){
        numeric_features = c(numeric_features,colnames(pros_train_numeric)[col])
      }
    }
    for(col in 1:(ncol(pros_train_categorical)-1)){
      pros_train_categorical[,col] = as.factor(pros_train_categorical[,col])
      pros_train_categorical_temp = pros_train_categorical[,c(col,ncol(pros_train_categorical))]
      infogain_test <- information_gain(PathStage~., data=pros_train_categorical_temp, type="infogain")
      if(round(infogain_test$importance,2)>0){
        categorical_features = c(categorical_features,colnames(pros_train_categorical)[col])
      }
      pros_train_categorical[,col] = as.numeric(pros_train_categorical[,col])
    }
    
    features = c(numeric_features, categorical_features)
    features = c(categorical_features, numeric_features,'PathStage')
    pros_train = pros_train[,features]
    pros_valid = pros_valid[,features]
    
    reg.fwd = regsubsets(PathStage~., data=pros_train, method="forward")
    R2adj = summary(reg.fwd)$adjr2
    R2adj.index = which.max(R2adj)
    mod = summary(reg.fwd)$which[R2adj.index,]
    vars = rownames(data.frame(which(mod)[-1]))
    feature.matrix[k,vars] = 1
    pros_train = pros_train[,c(vars,'PathStage')]
    pros_valid = pros_valid[,c(vars,'PathStage')]
    
    rf.o = randomForest(PathStage~., data=pros_train)
    rf.p = predict(rf.o, newdata=pros_valid)
    tb.rf = table(rf.p, pros_valid$PathStage)
    acc.rf[k] = sum(diag(tb.rf)) / sum(tb.rf)
    
    feature_frequency_matrix[,vars] = feature_frequency_matrix[,vars] + 1
    feature_stability_list = c(feature_stability_list, length(vars))
    
  }
  accuracy.avg.inner[r] = round(mean(acc.rf),2)
  
  fs = c()
  for(i in 1:ncol(feature.matrix)){
    if((sum(feature.matrix[,i])/nrow(feature.matrix)) >= 0.1){
      fs = c(fs, colnames(feature.matrix)[i])
    }
  }
  fs = c(fs, "PathStage")
  print("Features after Forward Stepwise")
  print(fs)
  features_selected[r] = list(fs)
  pc_train = pc_train[,fs] 
  pc_test = pc_test[,fs]
  rf.fit = randomForest(PathStage~., data=pc_train)
  rf.pred = predict(rf.fit, newdata=pc_test)
  table.rf = confusionMatrix(rf.pred, pc_test$PathStage, mode = 'everything')
  accuracy.outer[r] = round(table.rf$overall[1],2)
  precision[r] = round(table.rf$byClass[5],2)
  recall[r] = round(table.rf$byClass[6],2)
  f1.score[r] = round(table.rf$byClass[7],2)
  rf.p = predict(rf.fit, pc_test, type="prob")[,2]
  roc.rf = roc(response=pc_test$PathStage, predictor=rf.p)
  roc.outer[[r]] = roc.rf
  auc.outer[r] = round(roc.rf$auc,2)
  p.hat = length(which(pc_test$PathStage=='Early_Stage'))/(nrow(pc_test))
  standard_error[r] = sqrt((p.hat*(1-p.hat))/nrow(pc_train))
  
  sensitivities[r] = list(roc.rf$sensitivities)
  specificities[r] = list(roc.rf$specificities)
  thresholds[r] = list(roc.rf$thresholds)
}

features_with_frequency = c()
for(i in 1:ncol(feature_frequency_matrix)){
  if(feature_frequency_matrix[,i] > 0){
    features_with_frequency = c(features_with_frequency,colnames(feature_frequency_matrix)[i])
  }
}


ki_F_RF_Freq_feature_frequency_matrix = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(feature_frequency_matrix_3) = features_with_frequency


max_len_sensitivity = max(lengths(sensitivities))

for(i in 1:length(sensitivities)){
  if(length(sensitivities[[i]])<max_len_sensitivity){
    repeats = max_len_sensitivity-length(sensitivities[[i]])
    sensitivities[[i]] = append(sensitivities[[i]], rep(sensitivities[[i]][length(sensitivities[[i]])],times=repeats))
  }
}

max_len_specificity = max(lengths(specificities))

for(i in 1:length(specificities)){
  if(length(specificities[[i]])<max_len_specificity){
    repeats = max_len_specificity-length(specificities[[i]])
    specificities[[i]] = append(specificities[[i]], rep(specificities[[i]][length(specificities[[i]])],times=repeats))
  }
}


new_sensitivities = c()

for(i in 1:max_len_sensitivity){
  mean_sensitivity_i = (sensitivities[[1]][i]+sensitivities[[2]][i]+sensitivities[[3]][i]+sensitivities[[4]][i]+sensitivities[[5]][i])/5
  new_sensitivities = c(new_sensitivities,mean_sensitivity_i)
}
new_sensitivities

new_specificities = c()

for(i in 1:max_len_specificity){
  mean_specificity_i = (specificities[[1]][i]+specificities[[2]][i]+specificities[[3]][i]+specificities[[4]][i]+specificities[[5]][i])/5
  new_specificities = c(new_specificities,mean_specificity_i)
}
new_specificities

ki_F_RF_Freq_new_sensitivities_3 = new_sensitivities
ki_F_RF_Freq_new_specificities_3 = new_specificities

kw_infogain_forward_fs = feature_stability_list

#----8. Forward Stepwise Feature Selection with Kruskal-Wallis and Symmetrical Uncertainty----#
new_dat <- pros_new_dat
new_dat$PathStage = relevel(new_dat$PathStage, ref='Advanced_Stage')

R = 5
accuracy.outer = numeric(R)
accuracy.avg.inner = numeric(R)
auc.outer = numeric(R)
roc.outer = list(R)
precision = numeric(R)
recall = numeric(R)
f1.score = numeric(R)
standard_error = numeric(R)
features_selected = c()
feature_frequency_matrix = matrix(data=0,nrow = 1, ncol = ncol(pros_new_dat)-1)
colnames(feature_frequency_matrix) = colnames(pros_new_dat)[c(1:12)]
feature_stability_list = c()

sensitivities = c()
specificities = c()
thresholds = c()


set.seed(6090)
for(r in 1:R){
  n = nrow(new_dat)
  i.cv = sample(1:n, replace=FALSE)
  pc_copy = new_dat[i.cv,]
  pc.i <- createDataPartition(pc_copy$PathStage,p=0.70, list=F)
  pc_train = pc_copy[pc.i,]
  pc_test = pc_copy[-pc.i,]
  
  K = 10
  length = length(pc.i)
  folds = cut(1:length, K, labels=FALSE)
  acc.rf = numeric(K)
  
  feature.matrix = matrix(data=0,nrow = K, ncol = ncol(pc_train))
  feature.matrix[,1:10]
  colnames(feature.matrix) = colnames(pc_train)
  feature.matrix = feature.matrix[,-c(13)]
  
  set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    pros_train = pc_train[i.train,]
    pros_valid = pc_train[-i.train,]
    
    pros_train_numeric <- pros_train[,c(3:7,13)]
    pros_train_categorical <- pros_train[,c(1,2,8:12,13)]
    numeric_features = c()
    categorical_features = c()
    features = c()
    
    for(col in 1:(ncol(pros_train_numeric)-1)){
      kruskal_test <- kruskal.test(pros_train_numeric[,col]~pros_train_numeric$PathStage)
      if(kruskal_test$p.value<0.05){
        numeric_features = c(numeric_features,colnames(pros_train_numeric)[col])
      }
    }
    for(col in 1:(ncol(pros_train_categorical)-1)){
      pros_train_categorical[,col] = as.factor(pros_train_categorical[,col])
      pros_train_categorical_temp = pros_train_categorical[,c(col,ncol(pros_train_categorical))]
      infogain_test <- information_gain(PathStage~., data=pros_train_categorical_temp, type="symuncert")
      if(round(infogain_test$importance,2)>0){
        categorical_features = c(categorical_features,colnames(pros_train_categorical)[col])
      }
      pros_train_categorical[,col] = as.numeric(pros_train_categorical[,col])
    }
    
    features = c(numeric_features, categorical_features)
    features = c(categorical_features, numeric_features,'PathStage')
    pros_train = pros_train[,features]
    pros_valid = pros_valid[,features]
    
    reg.fwd = regsubsets(PathStage~., data=pros_train, method="forward")
    R2adj = summary(reg.fwd)$adjr2
    R2adj.index = which.max(R2adj)
    mod = summary(reg.fwd)$which[R2adj.index,]
    vars = rownames(data.frame(which(mod)[-1]))
    feature.matrix[k,vars] = 1
    pros_train = pros_train[,c(vars,'PathStage')]
    pros_valid = pros_valid[,c(vars,'PathStage')]
    
    rf.o = randomForest(PathStage~., data=pros_train)
    rf.p = predict(rf.o, newdata=pros_valid)
    tb.rf = table(rf.p, pros_valid$PathStage)
    acc.rf[k] = sum(diag(tb.rf)) / sum(tb.rf)
    
    feature_frequency_matrix[,vars] = feature_frequency_matrix[,vars] + 1
    feature_stability_list = c(feature_stability_list, length(vars))
    
  }
  accuracy.avg.inner[r] = round(mean(acc.rf),2)
  
  fs = c()
  for(i in 1:ncol(feature.matrix)){
    if((sum(feature.matrix[,i])/nrow(feature.matrix)) >= 0.1){
      fs = c(fs, colnames(feature.matrix)[i])
    }
  }
  fs = c(fs, "PathStage")
  print("Features after Forward Stepwise")
  print(fs)
  features_selected[r] = list(fs)
  pc_train = pc_train[,fs] 
  pc_test = pc_test[,fs]
  rf.fit = randomForest(PathStage~., data=pc_train)
  rf.pred = predict(rf.fit, newdata=pc_test)
  table.rf = confusionMatrix(rf.pred, pc_test$PathStage, mode = 'everything')
  accuracy.outer[r] = round(table.rf$overall[1],2)
  precision[r] = round(table.rf$byClass[5],2)
  recall[r] = round(table.rf$byClass[6],2)
  f1.score[r] = round(table.rf$byClass[7],2)
  rf.p = predict(rf.fit, pc_test, type="prob")[,2]
  roc.rf = roc(response=pc_test$PathStage, predictor=rf.p)
  roc.outer[[r]] = roc.rf
  auc.outer[r] = round(roc.rf$auc,2)
  p.hat = length(which(pc_test$PathStage=='Early_Stage'))/(nrow(pc_test))
  standard_error[r] = sqrt((p.hat*(1-p.hat))/nrow(pc_train))
  
  sensitivities[r] = list(roc.rf$sensitivities)
  specificities[r] = list(roc.rf$specificities)
  thresholds[r] = list(roc.rf$thresholds)
}

features_with_frequency = c()
for(i in 1:ncol(feature_frequency_matrix)){
  if(feature_frequency_matrix[,i] > 0){
    features_with_frequency = c(features_with_frequency,colnames(feature_frequency_matrix)[i])
  }
}


ks_F_RF_Freq = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(feature_frequency_matrix_4) = features_with_frequency



max_len_sensitivity = max(lengths(sensitivities))

for(i in 1:length(sensitivities)){
  if(length(sensitivities[[i]])<max_len_sensitivity){
    repeats = max_len_sensitivity-length(sensitivities[[i]])
    sensitivities[[i]] = append(sensitivities[[i]], rep(sensitivities[[i]][length(sensitivities[[i]])],times=repeats))
  }
}

max_len_specificity = max(lengths(specificities))

for(i in 1:length(specificities)){
  if(length(specificities[[i]])<max_len_specificity){
    repeats = max_len_specificity-length(specificities[[i]])
    specificities[[i]] = append(specificities[[i]], rep(specificities[[i]][length(specificities[[i]])],times=repeats))
  }
}


new_sensitivities = c()

for(i in 1:max_len_sensitivity){
  mean_sensitivity_i = (sensitivities[[1]][i]+sensitivities[[2]][i]+sensitivities[[3]][i]+sensitivities[[4]][i]+sensitivities[[5]][i])/5
  new_sensitivities = c(new_sensitivities,mean_sensitivity_i)
}
new_sensitivities

new_specificities = c()

for(i in 1:max_len_specificity){
  mean_specificity_i = (specificities[[1]][i]+specificities[[2]][i]+specificities[[3]][i]+specificities[[4]][i]+specificities[[5]][i])/5
  new_specificities = c(new_specificities,mean_specificity_i)
}
new_specificities

ks_F_RF_Freq_new_sensitivities_4 = new_sensitivities
ks_F_RF_Freq_new_specificities_4 = new_specificities



kw_symun_forward_fs = feature_stability_list


#----9. Backward Stepwise Feature Selection with Kruskal-Wallis and Chi-Squared----#
new_dat <- pros_new_dat
new_dat$PathStage = relevel(new_dat$PathStage, ref='Advanced_Stage')

R = 5
accuracy.outer = numeric(R)
accuracy.avg.inner = numeric(R)
auc.outer = numeric(R)
roc.outer = list(R)
precision = numeric(R)
recall = numeric(R)
f1.score = numeric(R)
standard_error = numeric(R)
features_selected = c()
feature_frequency_matrix = matrix(data=0,nrow = 1, ncol = ncol(pros_new_dat)-1)
colnames(feature_frequency_matrix) = colnames(pros_new_dat)[c(1:12)]
feature_stability_list = c()

sensitivities = c()
specificities = c()
thresholds = c()


set.seed(6090)
for(r in 1:R){
  n = nrow(new_dat)
  i.cv = sample(1:n, replace=FALSE)
  pc_copy = new_dat[i.cv,]
  pc.i <- createDataPartition(pc_copy$PathStage,p=0.70, list=F)
  pc_train = pc_copy[pc.i,]
  pc_test = pc_copy[-pc.i,]
  
  K = 10
  length = length(pc.i)
  folds = cut(1:length, K, labels=FALSE)
  acc.rf = numeric(K)
  
  feature.matrix = matrix(data=0,nrow = K, ncol = ncol(pc_train))
  feature.matrix[,1:10]
  colnames(feature.matrix) = colnames(pc_train)
  feature.matrix = feature.matrix[,-c(13)]
  
  set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    pros_train = pc_train[i.train,]
    pros_valid = pc_train[-i.train,]
    
    pros_train_numeric <- pros_train[,c(3:7,13)]
    pros_train_categorical <- pros_train[,c(1,2,8:12,13)]
    numeric_features = c()
    categorical_features = c()
    features = c()
    
    for(col in 1:(ncol(pros_train_numeric)-1)){
      kruskal_test <- kruskal.test(pros_train_numeric[,col]~pros_train_numeric$PathStage)
      if(kruskal_test$p.value<0.05){
        numeric_features = c(numeric_features,colnames(pros_train_numeric)[col])
      }
    }
    for(col in 1:(ncol(pros_train_categorical)-1)){
      pros_train_categorical[,col] = as.numeric(pros_train_categorical[,col])
      if(length(levels(as.factor(pros_train_categorical[,col])))>1){
        pros_train_categorical[,col] = as.factor(pros_train_categorical[,col])
        chisq_test <- chisq.test(pros_train_categorical[,col], pros_train_categorical$PathStage, correct = FALSE)
        if(round(chisq_test$p.value,2)<=0.05){
          categorical_features = c(categorical_features,colnames(pros_train_categorical)[col])
        }
      }
    }
    
    features = c(numeric_features, categorical_features)
    features = c(categorical_features, numeric_features,'PathStage')
    pros_train = pros_train[,features]
    pros_valid = pros_valid[,features]
    
    reg.fwd = regsubsets(PathStage~., data=pros_train, method="backward")
    R2adj = summary(reg.fwd)$adjr2
    R2adj.index = which.max(R2adj)
    mod = summary(reg.fwd)$which[R2adj.index,]
    vars = rownames(data.frame(which(mod)[-1]))
    feature.matrix[k,vars] = 1
    pros_train = pros_train[,c(vars,'PathStage')]
    pros_valid = pros_valid[,c(vars,'PathStage')]
    
    rf.o = randomForest(PathStage~., data=pros_train)
    rf.p = predict(rf.o, newdata=pros_valid)
    tb.rf = table(rf.p, pros_valid$PathStage)
    acc.rf[k] = sum(diag(tb.rf)) / sum(tb.rf)
    
    feature_frequency_matrix[,vars] = feature_frequency_matrix[,vars] + 1
    feature_stability_list = c(feature_stability_list, length(vars))
    
  }
  accuracy.avg.inner[r] = round(mean(acc.rf),2)
  
  fs = c()
  for(i in 1:ncol(feature.matrix)){
    if((sum(feature.matrix[,i])/nrow(feature.matrix)) >= 0.1){
      fs = c(fs, colnames(feature.matrix)[i])
    }
  }
  fs = c(fs, "PathStage")
  features_selected[r] = list(fs)
  pc_train = pc_train[,fs] 
  pc_test = pc_test[,fs]
  rf.fit = randomForest(PathStage~., data=pc_train)
  rf.pred = predict(rf.fit, newdata=pc_test)
  table.rf = confusionMatrix(rf.pred, pc_test$PathStage, mode = 'everything')
  accuracy.outer[r] = round(table.rf$overall[1],2)
  precision[r] = round(table.rf$byClass[5],2)
  recall[r] = round(table.rf$byClass[6],2)
  f1.score[r] = round(table.rf$byClass[7],2)
  rf.p = predict(rf.fit, pc_test, type="prob")[,2]
  roc.rf = roc(response=pc_test$PathStage, predictor=rf.p)
  roc.outer[[r]] = roc.rf
  auc.outer[r] = round(roc.rf$auc,2)
  p.hat = length(which(pc_test$PathStage=='Early_Stage'))/(nrow(pc_test))
  standard_error[r] = sqrt((p.hat*(1-p.hat))/nrow(pc_train))
  
  sensitivities[r] = list(roc.rf$sensitivities)
  specificities[r] = list(roc.rf$specificities)
  thresholds[r] = list(roc.rf$thresholds)
}

features_with_frequency = c()
for(i in 1:ncol(feature_frequency_matrix)){
  if(feature_frequency_matrix[,i] > 0){
    features_with_frequency = c(features_with_frequency,colnames(feature_frequency_matrix)[i])
  }
}

kwc_B_RF_Freq_feature_frequency_matrix = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(kwc_B_RF_Freq_feature_frequency_matrix) = features_with_frequency




max_len_sensitivity = max(lengths(sensitivities))

for(i in 1:length(sensitivities)){
  if(length(sensitivities[[i]])<max_len_sensitivity){
    repeats = max_len_sensitivity-length(sensitivities[[i]])
    sensitivities[[i]] = append(sensitivities[[i]], rep(sensitivities[[i]][length(sensitivities[[i]])],times=repeats))
  }
}

max_len_specificity = max(lengths(specificities))

for(i in 1:length(specificities)){
  if(length(specificities[[i]])<max_len_specificity){
    repeats = max_len_specificity-length(specificities[[i]])
    specificities[[i]] = append(specificities[[i]], rep(specificities[[i]][length(specificities[[i]])],times=repeats))
  }
}


new_sensitivities = c()

for(i in 1:max_len_sensitivity){
  mean_sensitivity_i = (sensitivities[[1]][i]+sensitivities[[2]][i]+sensitivities[[3]][i]+sensitivities[[4]][i]+sensitivities[[5]][i])/5
  new_sensitivities = c(new_sensitivities,mean_sensitivity_i)
}
new_sensitivities

new_specificities = c()

for(i in 1:max_len_specificity){
  mean_specificity_i = (specificities[[1]][i]+specificities[[2]][i]+specificities[[3]][i]+specificities[[4]][i]+specificities[[5]][i])/5
  new_specificities = c(new_specificities,mean_specificity_i)
}
new_specificities

kwc_B_RF_Freq_new_sensitivities_1 = new_sensitivities
kwc_B_RF_Freq_new_specificities_1 = new_specificities




kw_chi_backward_fs = feature_stability_list

#----10. Backward Stepwise Feature Selection with Kruskal-Wallis and One Rule----#
new_dat <- pros_new_dat
new_dat$PathStage = relevel(new_dat$PathStage, ref='Advanced_Stage')

R = 5
accuracy.outer = numeric(R)
accuracy.avg.inner = numeric(R)
auc.outer = numeric(R)
roc.outer = list(R)
precision = numeric(R)
recall = numeric(R)
f1.score = numeric(R)
standard_error = numeric(R)
features_selected = c()
feature_frequency_matrix = matrix(data=0,nrow = 1, ncol = ncol(pros_new_dat)-1)
colnames(feature_frequency_matrix) = colnames(pros_new_dat)[c(1:12)]
feature_stability_list = c()

sensitivities = c()
specificities = c()
thresholds = c()


set.seed(6090)
for(r in 1:R){
  n = nrow(new_dat)
  i.cv = sample(1:n, replace=FALSE)
  pc_copy = new_dat[i.cv,]
  pc.i <- createDataPartition(pc_copy$PathStage,p=0.70, list=F)
  pc_train = pc_copy[pc.i,]
  pc_test = pc_copy[-pc.i,]
  
  K = 10
  length = length(pc.i)
  folds = cut(1:length, K, labels=FALSE)
  acc.rf = numeric(K)
  
  feature.matrix = matrix(data=0,nrow = K, ncol = ncol(pc_train))
  feature.matrix[,1:10]
  colnames(feature.matrix) = colnames(pc_train)
  feature.matrix = feature.matrix[,-c(13)]
  
  set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    pros_train = pc_train[i.train,]
    pros_valid = pc_train[-i.train,]
    
    pros_train_numeric <- pros_train[,c(3:7,13)]
    pros_train_categorical <- pros_train[,c(1,2,8:12,13)]
    numeric_features = c()
    categorical_features = c()
    features = c()
    
    for(col in 1:(ncol(pros_train_numeric)-1)){
      kruskal_test <- kruskal.test(pros_train_numeric[,col]~pros_train_numeric$PathStage)
      if(kruskal_test$p.value<0.05){
        numeric_features = c(numeric_features,colnames(pros_train_numeric)[col])
      }
    }
    for(col in 1:(ncol(pros_train_categorical)-1)){
      oner_test <- oner_Type <- OneR::OneR(pros_train_categorical$PathStage~pros_train_categorical[,col])
      if((oner_test$correct_instances/oner_test$total_instances)>=0.6){
        categorical_features = c(categorical_features,colnames(pros_train_categorical)[col])
      }
    }
    
    features = c(numeric_features, categorical_features)
    features = c(categorical_features, numeric_features,'PathStage')
    pros_train = pros_train[,features]
    pros_valid = pros_valid[,features]
    
    reg.fwd = regsubsets(PathStage~., data=pros_train, method="backward")
    R2adj = summary(reg.fwd)$adjr2
    R2adj.index = which.max(R2adj)
    mod = summary(reg.fwd)$which[R2adj.index,]
    vars = rownames(data.frame(which(mod)[-1]))
    feature.matrix[k,vars] = 1
    pros_train = pros_train[,c(vars,'PathStage')]
    pros_valid = pros_valid[,c(vars,'PathStage')]
    
    rf.o = randomForest(PathStage~., data=pros_train)
    rf.p = predict(rf.o, newdata=pros_valid)
    tb.rf = table(rf.p, pros_valid$PathStage)
    acc.rf[k] = sum(diag(tb.rf)) / sum(tb.rf)
    
    feature_frequency_matrix[,vars] = feature_frequency_matrix[,vars] + 1
    feature_stability_list = c(feature_stability_list, length(vars))
    
  }
  accuracy.avg.inner[r] = round(mean(acc.rf),2)
  
  fs = c()
  for(i in 1:ncol(feature.matrix)){
    if((sum(feature.matrix[,i])/nrow(feature.matrix)) >= 0.1){
      fs = c(fs, colnames(feature.matrix)[i])
    }
  }
  fs = c(fs, "PathStage")
  features_selected[r] = list(fs)
  pc_train = pc_train[,fs] 
  pc_test = pc_test[,fs]
  rf.fit = randomForest(PathStage~., data=pc_train)
  rf.pred = predict(rf.fit, newdata=pc_test)
  table.rf = confusionMatrix(rf.pred, pc_test$PathStage, mode = 'everything')
  accuracy.outer[r] = round(table.rf$overall[1],2)
  precision[r] = round(table.rf$byClass[5],2)
  recall[r] = round(table.rf$byClass[6],2)
  f1.score[r] = round(table.rf$byClass[7],2)
  rf.p = predict(rf.fit, pc_test, type="prob")[,2]
  roc.rf = roc(response=pc_test$PathStage, predictor=rf.p)
  roc.outer[[r]] = roc.rf
  auc.outer[r] = round(roc.rf$auc,2)
  p.hat = length(which(pc_test$PathStage=='Early_Stage'))/(nrow(pc_test))
  standard_error[r] = sqrt((p.hat*(1-p.hat))/nrow(pc_train))
  
  sensitivities[r] = list(roc.rf$sensitivities)
  specificities[r] = list(roc.rf$specificities)
  thresholds[r] = list(roc.rf$thresholds)
}

features_with_frequency = c()
for(i in 1:ncol(feature_frequency_matrix)){
  if(feature_frequency_matrix[,i] > 0){
    features_with_frequency = c(features_with_frequency,colnames(feature_frequency_matrix)[i])
  }
}


kwo_B_RF_Freq_feature_frequency_matrix = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(kwo_B_RF_Freq_feature_frequency_matrix) = features_with_frequency


max_len_sensitivity = max(lengths(sensitivities))

for(i in 1:length(sensitivities)){
  if(length(sensitivities[[i]])<max_len_sensitivity){
    repeats = max_len_sensitivity-length(sensitivities[[i]])
    sensitivities[[i]] = append(sensitivities[[i]], rep(sensitivities[[i]][length(sensitivities[[i]])],times=repeats))
  }
}

max_len_specificity = max(lengths(specificities))

for(i in 1:length(specificities)){
  if(length(specificities[[i]])<max_len_specificity){
    repeats = max_len_specificity-length(specificities[[i]])
    specificities[[i]] = append(specificities[[i]], rep(specificities[[i]][length(specificities[[i]])],times=repeats))
  }
}


new_sensitivities = c()

for(i in 1:max_len_sensitivity){
  mean_sensitivity_i = (sensitivities[[1]][i]+sensitivities[[2]][i]+sensitivities[[3]][i]+sensitivities[[4]][i]+sensitivities[[5]][i])/5
  new_sensitivities = c(new_sensitivities,mean_sensitivity_i)
}
new_sensitivities

new_specificities = c()

for(i in 1:max_len_specificity){
  mean_specificity_i = (specificities[[1]][i]+specificities[[2]][i]+specificities[[3]][i]+specificities[[4]][i]+specificities[[5]][i])/5
  new_specificities = c(new_specificities,mean_specificity_i)
}
new_specificities

kwo_B_RF_Freq_new_sensitivities_2 = new_sensitivities
kwo_B_RF_Freq_new_specificities_2 = new_specificities

kw_oneR_backward_fs = feature_stability_list

#----11. Backward Stepwise Feature Selection with Kruskal-Wallis and Information Gain----#
new_dat <- pros_new_dat
new_dat$PathStage = relevel(new_dat$PathStage, ref='Advanced_Stage')

R = 5
accuracy.outer = numeric(R)
accuracy.avg.inner = numeric(R)
auc.outer = numeric(R)
roc.outer = list(R)
precision = numeric(R)
recall = numeric(R)
f1.score = numeric(R)
standard_error = numeric(R)
features_selected = c()
feature_frequency_matrix = matrix(data=0,nrow = 1, ncol = ncol(pros_new_dat)-1)
colnames(feature_frequency_matrix) = colnames(pros_new_dat)[c(1:12)]
feature_stability_list = c()

sensitivities = c()
specificities = c()
thresholds = c()


set.seed(6090)
for(r in 1:R){
  n = nrow(new_dat)
  i.cv = sample(1:n, replace=FALSE)
  pc_copy = new_dat[i.cv,]
  pc.i <- createDataPartition(pc_copy$PathStage,p=0.70, list=F)
  pc_train = pc_copy[pc.i,]
  pc_test = pc_copy[-pc.i,]
  
  K = 10
  length = length(pc.i)
  folds = cut(1:length, K, labels=FALSE)
  acc.rf = numeric(K)
  
  feature.matrix = matrix(data=0,nrow = K, ncol = ncol(pc_train))
  feature.matrix[,1:10]
  colnames(feature.matrix) = colnames(pc_train)
  feature.matrix = feature.matrix[,-c(13)]
  
  set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    pros_train = pc_train[i.train,]
    pros_valid = pc_train[-i.train,]
    
    pros_train_numeric <- pros_train[,c(3:7,13)]
    pros_train_categorical <- pros_train[,c(1,2,8:12,13)]
    numeric_features = c()
    categorical_features = c()
    features = c()
    
    for(col in 1:(ncol(pros_train_numeric)-1)){
      kruskal_test <- kruskal.test(pros_train_numeric[,col]~pros_train_numeric$PathStage)
      if(kruskal_test$p.value<0.05){
        numeric_features = c(numeric_features,colnames(pros_train_numeric)[col])
      }
    }
    for(col in 1:(ncol(pros_train_categorical)-1)){
      pros_train_categorical[,col] = as.factor(pros_train_categorical[,col])
      pros_train_categorical_temp = pros_train_categorical[,c(col,ncol(pros_train_categorical))]
      infogain_test <- information_gain(PathStage~., data=pros_train_categorical_temp, type="infogain")
      if(round(infogain_test$importance,2)>0){
        categorical_features = c(categorical_features,colnames(pros_train_categorical)[col])
      }
      pros_train_categorical[,col] = as.numeric(pros_train_categorical[,col])
    }
    
    features = c(numeric_features, categorical_features)
    features = c(categorical_features, numeric_features,'PathStage')
    pros_train = pros_train[,features]
    pros_valid = pros_valid[,features]
    
    reg.fwd = regsubsets(PathStage~., data=pros_train, method="backward")
    R2adj = summary(reg.fwd)$adjr2
    R2adj.index = which.max(R2adj)
    mod = summary(reg.fwd)$which[R2adj.index,]
    vars = rownames(data.frame(which(mod)[-1]))
    feature.matrix[k,vars] = 1
    pros_train = pros_train[,c(vars,'PathStage')]
    pros_valid = pros_valid[,c(vars,'PathStage')]
    
    rf.o = randomForest(PathStage~., data=pros_train)
    rf.p = predict(rf.o, newdata=pros_valid)
    tb.rf = table(rf.p, pros_valid$PathStage)
    acc.rf[k] = sum(diag(tb.rf)) / sum(tb.rf)
    
    feature_frequency_matrix[,vars] = feature_frequency_matrix[,vars] + 1
    feature_stability_list = c(feature_stability_list, length(vars))
    
  }
  accuracy.avg.inner[r] = round(mean(acc.rf),2)
  
  fs = c()
  for(i in 1:ncol(feature.matrix)){
    if((sum(feature.matrix[,i])/nrow(feature.matrix)) >= 0.1){
      fs = c(fs, colnames(feature.matrix)[i])
    }
  }
  fs = c(fs, "PathStage")
  features_selected[r] = list(fs)
  pc_train = pc_train[,fs] 
  pc_test = pc_test[,fs]
  rf.fit = randomForest(PathStage~., data=pc_train)
  rf.pred = predict(rf.fit, newdata=pc_test)
  table.rf = confusionMatrix(rf.pred, pc_test$PathStage, mode = 'everything')
  accuracy.outer[r] = round(table.rf$overall[1],2)
  precision[r] = round(table.rf$byClass[5],2)
  recall[r] = round(table.rf$byClass[6],2)
  f1.score[r] = round(table.rf$byClass[7],2)
  rf.p = predict(rf.fit, pc_test, type="prob")[,2]
  roc.rf = roc(response=pc_test$PathStage, predictor=rf.p)
  roc.outer[[r]] = roc.rf
  auc.outer[r] = round(roc.rf$auc,2)
  p.hat = length(which(pc_test$PathStage=='Early_Stage'))/(nrow(pc_test))
  standard_error[r] = sqrt((p.hat*(1-p.hat))/nrow(pc_train))
  
  sensitivities[r] = list(roc.rf$sensitivities)
  specificities[r] = list(roc.rf$specificities)
  thresholds[r] = list(roc.rf$thresholds)
}

features_with_frequency = c()
for(i in 1:ncol(feature_frequency_matrix)){
  if(feature_frequency_matrix[,i] > 0){
    features_with_frequency = c(features_with_frequency,colnames(feature_frequency_matrix)[i])
  }
}

kwi_B_RF_Freq_feature_frequency_matrix = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(feature_frequency_matrix_3) = features_with_frequency


max_len_sensitivity = max(lengths(sensitivities))

for(i in 1:length(sensitivities)){
  if(length(sensitivities[[i]])<max_len_sensitivity){
    repeats = max_len_sensitivity-length(sensitivities[[i]])
    sensitivities[[i]] = append(sensitivities[[i]], rep(sensitivities[[i]][length(sensitivities[[i]])],times=repeats))
  }
}

max_len_specificity = max(lengths(specificities))

for(i in 1:length(specificities)){
  if(length(specificities[[i]])<max_len_specificity){
    repeats = max_len_specificity-length(specificities[[i]])
    specificities[[i]] = append(specificities[[i]], rep(specificities[[i]][length(specificities[[i]])],times=repeats))
  }
}


new_sensitivities = c()

for(i in 1:max_len_sensitivity){
  mean_sensitivity_i = (sensitivities[[1]][i]+sensitivities[[2]][i]+sensitivities[[3]][i]+sensitivities[[4]][i]+sensitivities[[5]][i])/5
  new_sensitivities = c(new_sensitivities,mean_sensitivity_i)
}
new_sensitivities

new_specificities = c()

for(i in 1:max_len_specificity){
  mean_specificity_i = (specificities[[1]][i]+specificities[[2]][i]+specificities[[3]][i]+specificities[[4]][i]+specificities[[5]][i])/5
  new_specificities = c(new_specificities,mean_specificity_i)
}
new_specificities

kwi_B_RF_Freq_new_sensitivities_3 = new_sensitivities
kwi_B_RF_Freq_new_specificities_3 = new_specificities


kw_infogain_backward_fs = feature_stability_list

#----12. Backward Stepwise Feature Selection with Kruskal-Wallis and Symmetrical Uncertainty----#
new_dat <- pros_new_dat
new_dat$PathStage = relevel(new_dat$PathStage, ref='Advanced_Stage')

R = 5
accuracy.outer = numeric(R)
accuracy.avg.inner = numeric(R)
auc.outer = numeric(R)
roc.outer = list(R)
precision = numeric(R)
recall = numeric(R)
f1.score = numeric(R)
standard_error = numeric(R)
features_selected = c()
feature_frequency_matrix = matrix(data=0,nrow = 1, ncol = ncol(pros_new_dat)-1)
colnames(feature_frequency_matrix) = colnames(pros_new_dat)[c(1:12)]
feature_stability_list = c()

sensitivities = c()
specificities = c()
thresholds = c()


set.seed(6090)
for(r in 1:R){
  n = nrow(new_dat)
  i.cv = sample(1:n, replace=FALSE)
  pc_copy = new_dat[i.cv,]
  pc.i <- createDataPartition(pc_copy$PathStage,p=0.70, list=F)
  pc_train = pc_copy[pc.i,]
  pc_test = pc_copy[-pc.i,]
  
  K = 10
  length = length(pc.i)
  folds = cut(1:length, K, labels=FALSE)
  acc.rf = numeric(K)
  
  feature.matrix = matrix(data=0,nrow = K, ncol = ncol(pc_train))
  feature.matrix[,1:10]
  colnames(feature.matrix) = colnames(pc_train)
  feature.matrix = feature.matrix[,-c(13)]
  
  set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    pros_train = pc_train[i.train,]
    pros_valid = pc_train[-i.train,]
    
    pros_train_numeric <- pros_train[,c(3:7,13)]
    pros_train_categorical <- pros_train[,c(1,2,8:12,13)]
    numeric_features = c()
    categorical_features = c()
    features = c()
    
    for(col in 1:(ncol(pros_train_numeric)-1)){
      kruskal_test <- kruskal.test(pros_train_numeric[,col]~pros_train_numeric$PathStage)
      if(kruskal_test$p.value<0.05){
        numeric_features = c(numeric_features,colnames(pros_train_numeric)[col])
      }
    }
    for(col in 1:(ncol(pros_train_categorical)-1)){
      pros_train_categorical[,col] = as.factor(pros_train_categorical[,col])
      pros_train_categorical_temp = pros_train_categorical[,c(col,ncol(pros_train_categorical))]
      infogain_test <- information_gain(PathStage~., data=pros_train_categorical_temp, type="symuncert")
      if(round(infogain_test$importance,2)>0){
        categorical_features = c(categorical_features,colnames(pros_train_categorical)[col])
      }
      pros_train_categorical[,col] = as.numeric(pros_train_categorical[,col])
    }
    
    features = c(numeric_features, categorical_features)
    features = c(categorical_features, numeric_features,'PathStage')
    pros_train = pros_train[,features]
    pros_valid = pros_valid[,features]
    
    reg.fwd = regsubsets(PathStage~., data=pros_train, method="backward")
    R2adj = summary(reg.fwd)$adjr2
    R2adj.index = which.max(R2adj)
    mod = summary(reg.fwd)$which[R2adj.index,]
    vars = rownames(data.frame(which(mod)[-1]))
    feature.matrix[k,vars] = 1
    pros_train = pros_train[,c(vars,'PathStage')]
    pros_valid = pros_valid[,c(vars,'PathStage')]
    
    rf.o = randomForest(PathStage~., data=pros_train)
    rf.p = predict(rf.o, newdata=pros_valid)
    tb.rf = table(rf.p, pros_valid$PathStage)
    acc.rf[k] = sum(diag(tb.rf)) / sum(tb.rf)
    
    feature_frequency_matrix[,vars] = feature_frequency_matrix[,vars] + 1
    feature_stability_list = c(feature_stability_list, length(vars))
    
  }
  accuracy.avg.inner[r] = round(mean(acc.rf),2)
  
  fs = c()
  for(i in 1:ncol(feature.matrix)){
    if((sum(feature.matrix[,i])/nrow(feature.matrix)) >= 0.1){
      fs = c(fs, colnames(feature.matrix)[i])
    }
  }
  fs = c(fs, "PathStage")
  features_selected[r] = list(fs)
  pc_train = pc_train[,fs] 
  pc_test = pc_test[,fs]
  rf.fit = randomForest(PathStage~., data=pc_train)
  rf.pred = predict(rf.fit, newdata=pc_test)
  table.rf = confusionMatrix(rf.pred, pc_test$PathStage, mode = 'everything')
  accuracy.outer[r] = round(table.rf$overall[1],2)
  precision[r] = round(table.rf$byClass[5],2)
  recall[r] = round(table.rf$byClass[6],2)
  f1.score[r] = round(table.rf$byClass[7],2)
  rf.p = predict(rf.fit, pc_test, type="prob")[,2]
  roc.rf = roc(response=pc_test$PathStage, predictor=rf.p)
  roc.outer[[r]] = roc.rf
  auc.outer[r] = round(roc.rf$auc,2)
  p.hat = length(which(pc_test$PathStage=='Early_Stage'))/(nrow(pc_test))
  standard_error[r] = sqrt((p.hat*(1-p.hat))/nrow(pc_train))
  
  sensitivities[r] = list(roc.rf$sensitivities)
  specificities[r] = list(roc.rf$specificities)
  thresholds[r] = list(roc.rf$thresholds)
}

features_with_frequency = c()
for(i in 1:ncol(feature_frequency_matrix)){
  if(feature_frequency_matrix[,i] > 0){
    features_with_frequency = c(features_with_frequency,colnames(feature_frequency_matrix)[i])
  }
}

kws_B_RF_Freq_feature_frequency_matrix = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(kws_B_RF_Freq_feature_frequency_matrix) = features_with_frequency




max_len_sensitivity = max(lengths(sensitivities))

for(i in 1:length(sensitivities)){
  if(length(sensitivities[[i]])<max_len_sensitivity){
    repeats = max_len_sensitivity-length(sensitivities[[i]])
    sensitivities[[i]] = append(sensitivities[[i]], rep(sensitivities[[i]][length(sensitivities[[i]])],times=repeats))
  }
}

max_len_specificity = max(lengths(specificities))

for(i in 1:length(specificities)){
  if(length(specificities[[i]])<max_len_specificity){
    repeats = max_len_specificity-length(specificities[[i]])
    specificities[[i]] = append(specificities[[i]], rep(specificities[[i]][length(specificities[[i]])],times=repeats))
  }
}


new_sensitivities = c()

for(i in 1:max_len_sensitivity){
  mean_sensitivity_i = (sensitivities[[1]][i]+sensitivities[[2]][i]+sensitivities[[3]][i]+sensitivities[[4]][i]+sensitivities[[5]][i])/5
  new_sensitivities = c(new_sensitivities,mean_sensitivity_i)
}
new_sensitivities

new_specificities = c()

for(i in 1:max_len_specificity){
  mean_specificity_i = (specificities[[1]][i]+specificities[[2]][i]+specificities[[3]][i]+specificities[[4]][i]+specificities[[5]][i])/5
  new_specificities = c(new_specificities,mean_specificity_i)
}
new_specificities

kws_B_RF_Freq_new_sensitivities_4 = new_sensitivities
kws_B_RF_Freq_new_specificities_4 = new_specificities



kw_symun_backward_fs = feature_stability_list


#----13. RFE Feature Selection with Kruskal-Wallis and Chi-Squared----#

new_dat <- pros_new_dat
new_dat$PathStage = relevel(new_dat$PathStage, ref='Advanced_Stage')

R = 5
accuracy.outer = numeric(R)
accuracy.avg.inner = numeric(R)
auc.outer = numeric(R)
roc.outer = list(R)
precision = numeric(R)
recall = numeric(R)
f1.score = numeric(R)
standard_error = numeric(R)
features_selected = c()
feature_frequency_matrix = matrix(data=0,nrow = 1, ncol = ncol(pros_new_dat)-1)
colnames(feature_frequency_matrix) = colnames(pros_new_dat)[c(1:12)]
feature_stability_list = c()

sensitivities = c()
specificities = c()
thresholds = c()


set.seed(6090)
for(r in 1:R){
  n = nrow(new_dat)
  i.cv = sample(1:n, replace=FALSE)
  pc_copy = new_dat[i.cv,]
  pc.i <- createDataPartition(pc_copy$PathStage,p=0.70, list=F)
  pc_train = pc_copy[pc.i,]
  pc_test = pc_copy[-pc.i,]
  
  K = 10
  length = length(pc.i)
  folds = cut(1:length, K, labels=FALSE)
  acc.rf = numeric(K)
  
  feature.matrix = matrix(data=0,nrow = K, ncol = ncol(pc_train))
  feature.matrix[,1:10]
  colnames(feature.matrix) = colnames(pc_train)
  feature.matrix = feature.matrix[,-c(13)]
  
  set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    pros_train = pc_train[i.train,]
    pros_valid = pc_train[-i.train,]
    
    pros_train_numeric <- pros_train[,c(3:7,13)]
    pros_train_categorical <- pros_train[,c(1,2,8:12,13)]
    numeric_features = c()
    categorical_features = c()
    features = c()
    
    for(col in 1:(ncol(pros_train_numeric)-1)){
      kruskal_test <- kruskal.test(pros_train_numeric[,col]~pros_train_numeric$PathStage)
      if(kruskal_test$p.value<0.05){
        numeric_features = c(numeric_features,colnames(pros_train_numeric)[col])
      }
    }
    for(col in 1:(ncol(pros_train_categorical)-1)){
      pros_train_categorical[,col] = as.numeric(pros_train_categorical[,col])
      if(length(levels(as.factor(pros_train_categorical[,col])))>1){
        pros_train_categorical[,col] = as.factor(pros_train_categorical[,col])
        chisq_test <- chisq.test(pros_train_categorical[,col], pros_train_categorical$PathStage, correct = FALSE)
        if(round(chisq_test$p.value,2)<=0.05){
          categorical_features = c(categorical_features,colnames(pros_train_categorical)[col])
        }
      }
    }
    
    features = c(numeric_features, categorical_features)
    features = c(categorical_features, numeric_features,'PathStage')
    pros_train = pros_train[,features]
    pros_valid = pros_valid[,features]
    
    x.train = pros_train[,-c(ncol(pros_train))]
    y.train = pros_train$PathStage
    x.valid = pros_valid[,-c(ncol(pros_valid))]
    y.valid = pros_valid$PathStage
    
    subsets <- c(1:ncol(x.train))
    ctrl <- rfeControl(functions = rfFuncs,
                       method = "cv",
                       number = 10,
                       verbose = FALSE)
    rf.rfe <- rfe(x.train, y.train,
                  sizes = subsets,
                  rfeControl = ctrl)
    
    rfe_features = match(rf.rfe$optVariables,colnames(pros_train))
    feature.matrix[k,rf.rfe$optVariables] = 1
    pros_train = pros_train[,c(rf.rfe$optVariables,'PathStage')]
    pros_valid = pros_valid[,c(rf.rfe$optVariables,'PathStage')]
    
    rf.o = randomForest(PathStage~., data=pros_train)
    rf.p = predict(rf.o, newdata=pros_valid)
    tb.rf = table(rf.p, pros_valid$PathStage)
    acc.rf[k] = sum(diag(tb.rf)) / sum(tb.rf)
    
    feature_frequency_matrix[,rf.rfe$optVariables] = feature_frequency_matrix[,rf.rfe$optVariables] + 1
    feature_stability_list = c(feature_stability_list, length(rf.rfe$optVariables))
  }
  accuracy.avg.inner[r] = round(mean(acc.rf),2)
  
  fs = c()
  for(i in 1:ncol(feature.matrix)){
    if((sum(feature.matrix[,i])/nrow(feature.matrix)) >= 0.1){
      fs = c(fs, colnames(feature.matrix)[i])
    }
  }
  fs = c(fs, "PathStage")
  features_selected[r] = list(fs)
  pc_train = pc_train[,fs] 
  pc_test = pc_test[,fs]
  rf.fit = randomForest(PathStage~., data=pc_train)
  rf.pred = predict(rf.fit, newdata=pc_test)
  table.rf = confusionMatrix(rf.pred, pc_test$PathStage, mode = 'everything')
  accuracy.outer[r] = round(table.rf$overall[1],2)
  precision[r] = round(table.rf$byClass[5],2)
  recall[r] = round(table.rf$byClass[6],2)
  f1.score[r] = round(table.rf$byClass[7],2)
  rf.p = predict(rf.fit, pc_test, type="prob")[,2]
  roc.rf = roc(response=pc_test$PathStage, predictor=rf.p)
  roc.outer[[r]] = roc.rf
  auc.outer[r] = round(roc.rf$auc,2)
  p.hat = length(which(pc_test$PathStage=='Early_Stage'))/(nrow(pc_test))
  standard_error[r] = sqrt((p.hat*(1-p.hat))/nrow(pc_train))
  
  sensitivities[r] = list(roc.rf$sensitivities)
  specificities[r] = list(roc.rf$specificities)
  thresholds[r] = list(roc.rf$thresholds)
}

features_with_frequency = c()
for(i in 1:ncol(feature_frequency_matrix)){
  if(feature_frequency_matrix[,i] > 0){
    features_with_frequency = c(features_with_frequency,colnames(feature_frequency_matrix)[i])
  }
}



kwc_RFE_RF_Freq_feature_frequency_matrix = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(kwc_RFE_RF_Freq_feature_frequency_matrix) = features_with_frequency


max_len_sensitivity = max(lengths(sensitivities))

for(i in 1:length(sensitivities)){
  if(length(sensitivities[[i]])<max_len_sensitivity){
    repeats = max_len_sensitivity-length(sensitivities[[i]])
    sensitivities[[i]] = append(sensitivities[[i]], rep(sensitivities[[i]][length(sensitivities[[i]])],times=repeats))
  }
}

max_len_specificity = max(lengths(specificities))

for(i in 1:length(specificities)){
  if(length(specificities[[i]])<max_len_specificity){
    repeats = max_len_specificity-length(specificities[[i]])
    specificities[[i]] = append(specificities[[i]], rep(specificities[[i]][length(specificities[[i]])],times=repeats))
  }
}


new_sensitivities = c()

for(i in 1:max_len_sensitivity){
  mean_sensitivity_i = (sensitivities[[1]][i]+sensitivities[[2]][i]+sensitivities[[3]][i]+sensitivities[[4]][i]+sensitivities[[5]][i])/5
  new_sensitivities = c(new_sensitivities,mean_sensitivity_i)
}
new_sensitivities

new_specificities = c()

for(i in 1:max_len_specificity){
  mean_specificity_i = (specificities[[1]][i]+specificities[[2]][i]+specificities[[3]][i]+specificities[[4]][i]+specificities[[5]][i])/5
  new_specificities = c(new_specificities,mean_specificity_i)
}
new_specificities

kwc_RFE_RF_Freq_new_sensitivities_1 = new_sensitivities
kwc_RFE_RF_Freq_new_specificities_1 = new_specificities

kw_chi_rfe_fs = feature_stability_list

#----14. RFE Feature Selection with Kruskal-Wallis and One Rule----#
new_dat <- pros_new_dat
new_dat$PathStage = relevel(new_dat$PathStage, ref='Advanced_Stage')

R = 5
accuracy.outer = numeric(R)
accuracy.avg.inner = numeric(R)
auc.outer = numeric(R)
roc.outer = list(R)
precision = numeric(R)
recall = numeric(R)
f1.score = numeric(R)
standard_error = numeric(R)
features_selected = c()
feature_frequency_matrix = matrix(data=0,nrow = 1, ncol = ncol(pros_new_dat)-1)
colnames(feature_frequency_matrix) = colnames(pros_new_dat)[c(1:12)]
feature_stability_list = c()

sensitivities = c()
specificities = c()
thresholds = c()


set.seed(6090)
for(r in 1:R){
  n = nrow(new_dat)
  i.cv = sample(1:n, replace=FALSE)
  pc_copy = new_dat[i.cv,]
  pc.i <- createDataPartition(pc_copy$PathStage,p=0.70, list=F)
  pc_train = pc_copy[pc.i,]
  pc_test = pc_copy[-pc.i,]
  
  K = 10
  length = length(pc.i)
  folds = cut(1:length, K, labels=FALSE)
  acc.rf = numeric(K)
  
  feature.matrix = matrix(data=0,nrow = K, ncol = ncol(pc_train))
  feature.matrix[,1:10]
  colnames(feature.matrix) = colnames(pc_train)
  feature.matrix = feature.matrix[,-c(13)]
  
  set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    pros_train = pc_train[i.train,]
    pros_valid = pc_train[-i.train,]
    
    pros_train_numeric <- pros_train[,c(3:7,13)]
    pros_train_categorical <- pros_train[,c(1,2,8:12,13)]
    numeric_features = c()
    categorical_features = c()
    features = c()
    
    for(col in 1:(ncol(pros_train_numeric)-1)){
      kruskal_test <- kruskal.test(pros_train_numeric[,col]~pros_train_numeric$PathStage)
      if(kruskal_test$p.value<0.05){
        numeric_features = c(numeric_features,colnames(pros_train_numeric)[col])
      }
    }
    for(col in 1:(ncol(pros_train_categorical)-1)){
      oner_test <- oner_Type <- OneR::OneR(pros_train_categorical$PathStage~pros_train_categorical[,col])
      if((oner_test$correct_instances/oner_test$total_instances)>0.6){
        categorical_features = c(categorical_features,colnames(pros_train_categorical)[col])
      }
    }
    
    features = c(numeric_features, categorical_features)
    features = c(categorical_features, numeric_features,'PathStage')
    pros_train = pros_train[,features]
    pros_valid = pros_valid[,features]
    
    x.train = pros_train[,-c(ncol(pros_train))]
    y.train = pros_train$PathStage
    x.valid = pros_valid[,-c(ncol(pros_valid))]
    y.valid = pros_valid$PathStage
    
    subsets <- c(1:ncol(x.train))
    ctrl <- rfeControl(functions = rfFuncs,
                       method = "cv",
                       number = 10,
                       verbose = FALSE)
    rf.rfe <- rfe(x.train, y.train,
                  sizes = subsets,
                  rfeControl = ctrl)
    
    rfe_features = match(rf.rfe$optVariables,colnames(pros_train))
    feature.matrix[k,rf.rfe$optVariables] = 1
    pros_train = pros_train[,c(rf.rfe$optVariables,'PathStage')]
    pros_valid = pros_valid[,c(rf.rfe$optVariables,'PathStage')]
    
    rf.o = randomForest(PathStage~., data=pros_train)
    rf.p = predict(rf.o, newdata=pros_valid)
    tb.rf = table(rf.p, pros_valid$PathStage)
    acc.rf[k] = sum(diag(tb.rf)) / sum(tb.rf)
    
    feature_frequency_matrix[,rf.rfe$optVariables] = feature_frequency_matrix[,rf.rfe$optVariables] + 1
    feature_stability_list = c(feature_stability_list, length(rf.rfe$optVariables))
  }
  accuracy.avg.inner[r] = round(mean(acc.rf),2)
  
  fs = c()
  for(i in 1:ncol(feature.matrix)){
    if((sum(feature.matrix[,i])/nrow(feature.matrix)) >= 0.1){
      fs = c(fs, colnames(feature.matrix)[i])
    }
  }
  fs = c(fs, "PathStage")
  features_selected[r] = list(fs)
  pc_train = pc_train[,fs] 
  pc_test = pc_test[,fs]
  rf.fit = randomForest(PathStage~., data=pc_train)
  rf.pred = predict(rf.fit, newdata=pc_test)
  table.rf = confusionMatrix(rf.pred, pc_test$PathStage, mode = 'everything')
  accuracy.outer[r] = round(table.rf$overall[1],2)
  precision[r] = round(table.rf$byClass[5],2)
  recall[r] = round(table.rf$byClass[6],2)
  f1.score[r] = round(table.rf$byClass[7],2)
  rf.p = predict(rf.fit, pc_test, type="prob")[,2]
  roc.rf = roc(response=pc_test$PathStage, predictor=rf.p)
  roc.outer[[r]] = roc.rf
  auc.outer[r] = round(roc.rf$auc,2)
  p.hat = length(which(pc_test$PathStage=='Early_Stage'))/(nrow(pc_test))
  standard_error[r] = sqrt((p.hat*(1-p.hat))/nrow(pc_train))
  
  sensitivities[r] = list(roc.rf$sensitivities)
  specificities[r] = list(roc.rf$specificities)
  thresholds[r] = list(roc.rf$thresholds)
}

features_with_frequency = c()
for(i in 1:ncol(feature_frequency_matrix)){
  if(feature_frequency_matrix[,i] > 0){
    features_with_frequency = c(features_with_frequency,colnames(feature_frequency_matrix)[i])
  }
}


kwo_RFE_RF_Freq_feature_frequency_matrix = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(kwo_RFE_RF_Freq_feature_frequency_matrix) = features_with_frequency




max_len_sensitivity = max(lengths(sensitivities))

for(i in 1:length(sensitivities)){
  if(length(sensitivities[[i]])<max_len_sensitivity){
    repeats = max_len_sensitivity-length(sensitivities[[i]])
    sensitivities[[i]] = append(sensitivities[[i]], rep(sensitivities[[i]][length(sensitivities[[i]])],times=repeats))
  }
}

max_len_specificity = max(lengths(specificities))

for(i in 1:length(specificities)){
  if(length(specificities[[i]])<max_len_specificity){
    repeats = max_len_specificity-length(specificities[[i]])
    specificities[[i]] = append(specificities[[i]], rep(specificities[[i]][length(specificities[[i]])],times=repeats))
  }
}


new_sensitivities = c()

for(i in 1:max_len_sensitivity){
  mean_sensitivity_i = (sensitivities[[1]][i]+sensitivities[[2]][i]+sensitivities[[3]][i]+sensitivities[[4]][i]+sensitivities[[5]][i])/5
  new_sensitivities = c(new_sensitivities,mean_sensitivity_i)
}
new_sensitivities

new_specificities = c()

for(i in 1:max_len_specificity){
  mean_specificity_i = (specificities[[1]][i]+specificities[[2]][i]+specificities[[3]][i]+specificities[[4]][i]+specificities[[5]][i])/5
  new_specificities = c(new_specificities,mean_specificity_i)
}
new_specificities

kwo_RFE_RF_Freq_new_sensitivities_2 = new_sensitivities
kwo_RFE_RF_Freq_new_specificities_2 = new_specificities

kw_oner_rfe_fs = feature_stability_list

#----15. RFE Feature Selection with Kruskal-Wallis and Information Gain----#
new_dat <- pros_new_dat
new_dat$PathStage = relevel(new_dat$PathStage, ref='Advanced_Stage')

R = 5
accuracy.outer = numeric(R)
accuracy.avg.inner = numeric(R)
auc.outer = numeric(R)
roc.outer = list(R)
precision = numeric(R)
recall = numeric(R)
f1.score = numeric(R)
standard_error = numeric(R)
features_selected = c()
feature_frequency_matrix = matrix(data=0,nrow = 1, ncol = ncol(pros_new_dat)-1)
colnames(feature_frequency_matrix) = colnames(pros_new_dat)[c(1:12)]
feature_stability_list = c()

sensitivities = c()
specificities = c()
thresholds = c()


set.seed(6090)
for(r in 1:R){
  n = nrow(new_dat)
  i.cv = sample(1:n, replace=FALSE)
  pc_copy = new_dat[i.cv,]
  pc.i <- createDataPartition(pc_copy$PathStage,p=0.70, list=F)
  pc_train = pc_copy[pc.i,]
  pc_test = pc_copy[-pc.i,]
  
  K = 10
  length = length(pc.i)
  folds = cut(1:length, K, labels=FALSE)
  acc.rf = numeric(K)
  
  feature.matrix = matrix(data=0,nrow = K, ncol = ncol(pc_train))
  feature.matrix[,1:10]
  colnames(feature.matrix) = colnames(pc_train)
  feature.matrix = feature.matrix[,-c(13)]
  
  set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    pros_train = pc_train[i.train,]
    pros_valid = pc_train[-i.train,]
    
    pros_train_numeric <- pros_train[,c(3:7,13)]
    pros_train_categorical <- pros_train[,c(1,2,8:12,13)]
    numeric_features = c()
    categorical_features = c()
    features = c()
    
    for(col in 1:(ncol(pros_train_numeric)-1)){
      kruskal_test <- kruskal.test(pros_train_numeric[,col]~pros_train_numeric$PathStage)
      if(kruskal_test$p.value<0.05){
        numeric_features = c(numeric_features,colnames(pros_train_numeric)[col])
      }
    }
    for(col in 1:(ncol(pros_train_categorical)-1)){
      pros_train_categorical[,col] = as.factor(pros_train_categorical[,col])
      pros_train_categorical_temp = pros_train_categorical[,c(col,ncol(pros_train_categorical))]
      infogain_test <- information_gain(PathStage~., data=pros_train_categorical_temp, type="infogain")
      if(round(infogain_test$importance,2)>0){
        categorical_features = c(categorical_features,colnames(pros_train_categorical)[col])
      }
      pros_train_categorical[,col] = as.numeric(pros_train_categorical[,col])
    }
    
    features = c(numeric_features, categorical_features)
    features = c(categorical_features, numeric_features,'PathStage')
    pros_train = pros_train[,features]
    pros_valid = pros_valid[,features]
    
    x.train = pros_train[,-c(ncol(pros_train))]
    y.train = pros_train$PathStage
    x.valid = pros_valid[,-c(ncol(pros_valid))]
    y.valid = pros_valid$PathStage
    
    subsets <- c(1:ncol(x.train))
    ctrl <- rfeControl(functions = rfFuncs,
                       method = "cv",
                       number = 10,
                       verbose = FALSE)
    rf.rfe <- rfe(x.train, y.train,
                  sizes = subsets,
                  rfeControl = ctrl)
    
    rfe_features = match(rf.rfe$optVariables,colnames(pros_train))
    feature.matrix[k,rf.rfe$optVariables] = 1
    pros_train = pros_train[,c(rf.rfe$optVariables,'PathStage')]
    pros_valid = pros_valid[,c(rf.rfe$optVariables,'PathStage')]
    
    rf.o = randomForest(PathStage~., data=pros_train)
    rf.p = predict(rf.o, newdata=pros_valid)
    tb.rf = table(rf.p, pros_valid$PathStage)
    acc.rf[k] = sum(diag(tb.rf)) / sum(tb.rf)
    
    feature_frequency_matrix[,rf.rfe$optVariables] = feature_frequency_matrix[,rf.rfe$optVariables] + 1
    feature_stability_list = c(feature_stability_list, length(rf.rfe$optVariables))
  }
  accuracy.avg.inner[r] = round(mean(acc.rf),2)
  
  fs = c()
  for(i in 1:ncol(feature.matrix)){
    if((sum(feature.matrix[,i])/nrow(feature.matrix)) >= 0.1){
      fs = c(fs, colnames(feature.matrix)[i])
    }
  }
  fs = c(fs, "PathStage")
  features_selected[r] = list(fs)
  pc_train = pc_train[,fs] 
  pc_test = pc_test[,fs]
  rf.fit = randomForest(PathStage~., data=pc_train)
  rf.pred = predict(rf.fit, newdata=pc_test)
  table.rf = confusionMatrix(rf.pred, pc_test$PathStage, mode = 'everything')
  accuracy.outer[r] = round(table.rf$overall[1],2)
  precision[r] = round(table.rf$byClass[5],2)
  recall[r] = round(table.rf$byClass[6],2)
  f1.score[r] = round(table.rf$byClass[7],2)
  rf.p = predict(rf.fit, pc_test, type="prob")[,2]
  roc.rf = roc(response=pc_test$PathStage, predictor=rf.p)
  roc.outer[[r]] = roc.rf
  auc.outer[r] = round(roc.rf$auc,2)
  p.hat = length(which(pc_test$PathStage=='Early_Stage'))/(nrow(pc_test))
  standard_error[r] = sqrt((p.hat*(1-p.hat))/nrow(pc_train))
  
  sensitivities[r] = list(roc.rf$sensitivities)
  specificities[r] = list(roc.rf$specificities)
  thresholds[r] = list(roc.rf$thresholds)
}

features_with_frequency = c()
for(i in 1:ncol(feature_frequency_matrix)){
  if(feature_frequency_matrix[,i] > 0){
    features_with_frequency = c(features_with_frequency,colnames(feature_frequency_matrix)[i])
  }
}

kwi_RFE_RF_Freq_feature_frequency_matrix = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(kwi_RFE_RF_Freq_feature_frequency_matrix) = features_with_frequency


max_len_sensitivity = max(lengths(sensitivities))

for(i in 1:length(sensitivities)){
  if(length(sensitivities[[i]])<max_len_sensitivity){
    repeats = max_len_sensitivity-length(sensitivities[[i]])
    sensitivities[[i]] = append(sensitivities[[i]], rep(sensitivities[[i]][length(sensitivities[[i]])],times=repeats))
  }
}

max_len_specificity = max(lengths(specificities))

for(i in 1:length(specificities)){
  if(length(specificities[[i]])<max_len_specificity){
    repeats = max_len_specificity-length(specificities[[i]])
    specificities[[i]] = append(specificities[[i]], rep(specificities[[i]][length(specificities[[i]])],times=repeats))
  }
}


new_sensitivities = c()

for(i in 1:max_len_sensitivity){
  mean_sensitivity_i = (sensitivities[[1]][i]+sensitivities[[2]][i]+sensitivities[[3]][i]+sensitivities[[4]][i]+sensitivities[[5]][i])/5
  new_sensitivities = c(new_sensitivities,mean_sensitivity_i)
}
new_sensitivities

new_specificities = c()

for(i in 1:max_len_specificity){
  mean_specificity_i = (specificities[[1]][i]+specificities[[2]][i]+specificities[[3]][i]+specificities[[4]][i]+specificities[[5]][i])/5
  new_specificities = c(new_specificities,mean_specificity_i)
}
new_specificities

kwi_RFE_RF_Freq_new_sensitivities_3 = new_sensitivities
kwi_RFE_RF_Freq_new_specificities_3 = new_specificities

kw_infogain_rfe_fs = feature_stability_list

#----16. RFE Feature Selection with Kruskal-Wallis and Symmetrical Uncertainty----#

new_dat <- pros_new_dat
new_dat$PathStage = relevel(new_dat$PathStage, ref='Advanced_Stage')

R = 5
accuracy.outer = numeric(R)
accuracy.avg.inner = numeric(R)
auc.outer = numeric(R)
roc.outer = list(R)
precision = numeric(R)
recall = numeric(R)
f1.score = numeric(R)
standard_error = numeric(R)
features_selected = c()
feature_frequency_matrix = matrix(data=0,nrow = 1, ncol = ncol(pros_new_dat)-1)
colnames(feature_frequency_matrix) = colnames(pros_new_dat)[c(1:12)]
feature_stability_list = c()

sensitivities = c()
specificities = c()
thresholds = c()


set.seed(6090)
for(r in 1:R){
  n = nrow(new_dat)
  i.cv = sample(1:n, replace=FALSE)
  pc_copy = new_dat[i.cv,]
  pc.i <- createDataPartition(pc_copy$PathStage,p=0.70, list=F)
  pc_train = pc_copy[pc.i,]
  pc_test = pc_copy[-pc.i,]
  
  K = 10
  length = length(pc.i)
  folds = cut(1:length, K, labels=FALSE)
  acc.rf = numeric(K)
  
  feature.matrix = matrix(data=0,nrow = K, ncol = ncol(pc_train))
  feature.matrix[,1:10]
  colnames(feature.matrix) = colnames(pc_train)
  feature.matrix = feature.matrix[,-c(13)]
  
  set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    pros_train = pc_train[i.train,]
    pros_valid = pc_train[-i.train,]
    
    pros_train_numeric <- pros_train[,c(3:7,13)]
    pros_train_categorical <- pros_train[,c(1,2,8:12,13)]
    numeric_features = c()
    categorical_features = c()
    features = c()
    
    for(col in 1:(ncol(pros_train_numeric)-1)){
      kruskal_test <- kruskal.test(pros_train_numeric[,col]~pros_train_numeric$PathStage)
      if(kruskal_test$p.value<0.05){
        numeric_features = c(numeric_features,colnames(pros_train_numeric)[col])
      }
    }
    for(col in 1:(ncol(pros_train_categorical)-1)){
      pros_train_categorical[,col] = as.factor(pros_train_categorical[,col])
      pros_train_categorical_temp = pros_train_categorical[,c(col,ncol(pros_train_categorical))]
      infogain_test <- information_gain(PathStage~., data=pros_train_categorical_temp, type="symuncert")
      if(round(infogain_test$importance,2)>0){
        categorical_features = c(categorical_features,colnames(pros_train_categorical)[col])
      }
      pros_train_categorical[,col] = as.numeric(pros_train_categorical[,col])
    }
    
    features = c(numeric_features, categorical_features)
    features = c(categorical_features, numeric_features,'PathStage')
    pros_train = pros_train[,features]
    pros_valid = pros_valid[,features]
    
    x.train = pros_train[,-c(ncol(pros_train))]
    y.train = pros_train$PathStage
    x.valid = pros_valid[,-c(ncol(pros_valid))]
    y.valid = pros_valid$PathStage
    
    subsets <- c(1:ncol(x.train))
    ctrl <- rfeControl(functions = rfFuncs,
                       method = "cv",
                       number = 10,
                       verbose = FALSE)
    rf.rfe <- rfe(x.train, y.train,
                  sizes = subsets,
                  rfeControl = ctrl)
    
    rfe_features = match(rf.rfe$optVariables,colnames(pros_train))
    feature.matrix[k,rf.rfe$optVariables] = 1
    pros_train = pros_train[,c(rf.rfe$optVariables,'PathStage')]
    pros_valid = pros_valid[,c(rf.rfe$optVariables,'PathStage')]
    
    rf.o = randomForest(PathStage~., data=pros_train)
    rf.p = predict(rf.o, newdata=pros_valid)
    tb.rf = table(rf.p, pros_valid$PathStage)
    acc.rf[k] = sum(diag(tb.rf)) / sum(tb.rf)
    
    feature_frequency_matrix[,rf.rfe$optVariables] = feature_frequency_matrix[,rf.rfe$optVariables] + 1
    feature_stability_list = c(feature_stability_list, length(rf.rfe$optVariables))
  }
  accuracy.avg.inner[r] = round(mean(acc.rf),2)
  
  fs = c()
  for(i in 1:ncol(feature.matrix)){
    if((sum(feature.matrix[,i])/nrow(feature.matrix)) >= 0.1){
      fs = c(fs, colnames(feature.matrix)[i])
    }
  }
  fs = c(fs, "PathStage")
  features_selected[r] = list(fs)
  pc_train = pc_train[,fs] 
  pc_test = pc_test[,fs]
  rf.fit = randomForest(PathStage~., data=pc_train)
  rf.pred = predict(rf.fit, newdata=pc_test)
  table.rf = confusionMatrix(rf.pred, pc_test$PathStage, mode = 'everything')
  accuracy.outer[r] = round(table.rf$overall[1],2)
  precision[r] = round(table.rf$byClass[5],2)
  recall[r] = round(table.rf$byClass[6],2)
  f1.score[r] = round(table.rf$byClass[7],2)
  rf.p = predict(rf.fit, pc_test, type="prob")[,2]
  roc.rf = roc(response=pc_test$PathStage, predictor=rf.p)
  roc.outer[[r]] = roc.rf
  auc.outer[r] = round(roc.rf$auc,2)
  p.hat = length(which(pc_test$PathStage=='Early_Stage'))/(nrow(pc_test))
  standard_error[r] = sqrt((p.hat*(1-p.hat))/nrow(pc_train))
  
  sensitivities[r] = list(roc.rf$sensitivities)
  specificities[r] = list(roc.rf$specificities)
  thresholds[r] = list(roc.rf$thresholds)
}

features_with_frequency = c()
for(i in 1:ncol(feature_frequency_matrix)){
  if(feature_frequency_matrix[,i] > 0){
    features_with_frequency = c(features_with_frequency,colnames(feature_frequency_matrix)[i])
  }
}


kws_RFE_RF_Freq_feature_frequency_matrix = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(kws_RFE_RF_Freq_feature_frequency_matrix) = features_with_frequency


max_len_sensitivity = max(lengths(sensitivities))

for(i in 1:length(sensitivities)){
  if(length(sensitivities[[i]])<max_len_sensitivity){
    repeats = max_len_sensitivity-length(sensitivities[[i]])
    sensitivities[[i]] = append(sensitivities[[i]], rep(sensitivities[[i]][length(sensitivities[[i]])],times=repeats))
  }
}

max_len_specificity = max(lengths(specificities))

for(i in 1:length(specificities)){
  if(length(specificities[[i]])<max_len_specificity){
    repeats = max_len_specificity-length(specificities[[i]])
    specificities[[i]] = append(specificities[[i]], rep(specificities[[i]][length(specificities[[i]])],times=repeats))
  }
}


new_sensitivities = c()

for(i in 1:max_len_sensitivity){
  mean_sensitivity_i = (sensitivities[[1]][i]+sensitivities[[2]][i]+sensitivities[[3]][i]+sensitivities[[4]][i]+sensitivities[[5]][i])/5
  new_sensitivities = c(new_sensitivities,mean_sensitivity_i)
}
new_sensitivities

new_specificities = c()

for(i in 1:max_len_specificity){
  mean_specificity_i = (specificities[[1]][i]+specificities[[2]][i]+specificities[[3]][i]+specificities[[4]][i]+specificities[[5]][i])/5
  new_specificities = c(new_specificities,mean_specificity_i)
}
new_specificities

kws_RFE_RF_Freq_new_sensitivities_4 = new_sensitivities
kws_RFE_RF_Freq_new_specificities_4 = new_specificities


kw_symun_rfe_fs = feature_stability_list

#----Feature Frequency Plots----#
KWc_L_RF_Freq <- barplot(KWc_L_RF_Freq_feature_frequency_matrix[1,], xaxt="n", beside=TRUE, cex.names=0.55, las=3, col=c('red'), main='Kruskal-Wallis - Chi-Sq - LASSO')
axis(2,at=seq(0,50,5), lwd=2.0)
xlabs <- paste(colnames(KWc_L_RF_Freq_feature_frequency_matrix))
text(KWc_L_RF_Freq, par("usr")[3] - 0.025, srt = 45, adj = 1, 
     labels = xlabs, 
     xpd = TRUE, font = 2)

KWo_L_RF_Freq <- barplot(KWo_L_RF_Freq_feature_frequency_matrix[1,], xaxt="n", beside=TRUE, cex.names=0.55, las=3, col=c('blue'), main='Kruskal-Wallis - oneR - LASSO')
axis(2,at=seq(0,50,5), lwd=2.0)
xlabs <- paste(colnames(KWo_L_RF_Freq_feature_frequency_matrix))
text(KWo_L_RF_Freq, par("usr")[3] - 0.025, srt = 45, adj = 1, 
     labels = xlabs, 
     xpd = TRUE, font = 2)

kwi_L_RF_Freq <- barplot(kwi_L_RF_Freq_feature_frequency_matrix[1,], xaxt="n", beside=TRUE, cex.names=0.55, las=3, col=c('sea green'), main='Kruskal-Wallis - InfoGain - LASSO')
axis(2,at=seq(0,50,5), lwd=2.0)
xlabs <- paste(colnames(kwi_L_RF_Freq_feature_frequency_matrix))
text(kwi_L_RF_Freq, par("usr")[3] - 0.025, srt = 45, adj = 1, 
     labels = xlabs, 
     xpd = TRUE, font = 2)

kws_L_RF_Freq <- barplot(kws_L_RF_Freq_feature_frequency_matrix[1,], xaxt="n", beside=TRUE, cex.names=0.55, las=3, col=c('black'), main='Kruskal-Wallis - SymUncert - LASSO')
axis(2,at=seq(0,50,5), lwd=2.0)
xlabs <- paste(colnames(kws_L_RF_Freq_feature_frequency_matrix))
text(kws_L_RF_Freq, par("usr")[3] - 0.025, srt = 45, adj = 1, 
     labels = xlabs, 
     xpd = TRUE, font = 2)

kC_F_RF_Freq <- barplot(kC_F_RF_Freq_feature_frequency_matrix[1,], xaxt="n", beside=TRUE, cex.names=0.55, las=3, col=c('red'), main='Kruskal-Wallis - Chi-Sq - Forward')
axis(2,at=seq(0,50,5), lwd=2.0)
xlabs <- paste(colnames(kC_F_RF_Freq_feature_frequency_matrix))
text(kC_F_RF_Freq, par("usr")[3] - 0.025, srt = 45, adj = 1, 
     labels = xlabs, 
     xpd = TRUE, font = 2)

ko_F_RF_Freq <- barplot(ko_F_RF_Freq_feature_frequency_matrix[1,], xaxt="n", beside=TRUE, cex.names=0.55, las=3, col=c('blue'), main='Kruskal-Wallis - oneR - Forward')
axis(2,at=seq(0,50,5), lwd=2.0)
xlabs <- paste(colnames(ko_F_RF_Freq_feature_frequency_matrix))
text(ko_F_RF_Freq, par("usr")[3] - 0.025, srt = 45, adj = 1, 
     labels = xlabs, 
     xpd = TRUE, font = 2)


ki_F_RF_Freq <- barplot(ki_F_RF_Freq_feature_frequency_matrix[1,], xaxt="n", beside=TRUE, cex.names=0.55, las=3, col=c('sea green'), main='Kruskal-Wallis - InfoGain - Forward')
axis(2,at=seq(0,50,5), lwd=2.0)
xlabs <- paste(colnames(ki_F_RF_Freq_feature_frequency_matrix))
text(ki_F_RF_Freq, par("usr")[3] - 0.025, srt = 45, adj = 1, 
     labels = xlabs, 
     xpd = TRUE, font = 2)


ks_F_RF_Freq <- barplot(feature_frequency_matrix_4[1,], xaxt="n", beside=TRUE, cex.names=0.55, las=3, col=c('black'), main='Kruskal-Wallis - SymUncert - Forward')
axis(2,at=seq(0,50,5), lwd=2.0)
xlabs <- paste(colnames(feature_frequency_matrix_4))
text(x, par("usr")[3] - 0.025, srt = 45, adj = 1, 
     labels = xlabs, 
     xpd = TRUE, font = 2)

kwc_B_RF_Freq <- barplot(kwc_B_RF_Freq_feature_frequency_matrix[1,], xaxt="n", beside=TRUE, cex.names=0.55, las=3, col=c('red'), main='Kruskal-Wallis - Chi-Sq - Backward')
axis(2,at=seq(0,50,5), lwd=2.0)
xlabs <- paste(colnames(kwc_B_RF_Freq_feature_frequency_matrix))
text(kwc_B_RF_Freq, par("usr")[3] - 0.025, srt = 45, adj = 1, 
     labels = xlabs, 
     xpd = TRUE, font = 2)


kwo_B_RF_Freq <- barplot(feature_frequency_matrix_2[1,], xaxt="n", beside=TRUE, cex.names=0.55, las=3, col=c('blue'), main='Kruskal-Wallis - oneR - Backward')
axis(2,at=seq(0,50,5), lwd=2.0)
xlabs <- paste(colnames(feature_frequency_matrix_2))
text(kwo_B_RF_Freq, par("usr")[3] - 0.025, srt = 45, adj = 1, 
     labels = xlabs, 
     xpd = TRUE, font = 2)


kwi_B_RF_Freq <- barplot(feature_frequency_matrix_3[1,], xaxt="n", beside=TRUE, cex.names=0.55, las=3, col=c('sea green'), main='Kruskal-Wallis - InfoGain - Backward')
axis(2,at=seq(0,50,5), lwd=2.0)
xlabs <- paste(colnames(feature_frequency_matrix_3))
text(kwi_B_RF_Freq, par("usr")[3] - 0.025, srt = 45, adj = 1, 
     labels = xlabs, 
     xpd = TRUE, font = 2)


kws_B_RF_Freq <- barplot(kws_B_RF_Freq_feature_frequency_matrix[1,], xaxt="n", beside=TRUE, cex.names=0.55, las=3, col=c('black'), main='Kruskal-Wallis - SymUncert - Backward')
axis(2,at=seq(0,50,5), lwd=2.0)
xlabs <- paste(colnames(kws_B_RF_Freq_feature_frequency_matrix))
text(kws_B_RF_Freq, par("usr")[3] - 0.025, srt = 45, adj = 1, 
     labels = xlabs, 
     xpd = TRUE, font = 2)

kwc_RFE_RF_Freq <- barplot(kwc_RFE_RF_Freq_feature_frequency_matrix[1,], xaxt="n", beside=TRUE, cex.names=0.55, las=3, col=c('red'), main='Kruskal-Wallis - Chi-Sq - RFE')
axis(2,at=seq(0,50,5), lwd=2.0)
xlabs <- paste(colnames(kwc_RFE_RF_Freq_feature_frequency_matrix))
text(kwc_RFE_RF_Freq, par("usr")[3] - 0.025, srt = 45, adj = 1, 
     labels = xlabs, 
     xpd = TRUE, font = 2)

kwo_RFE_RF_Freq <- barplot(kwo_RFE_RF_Freq_feature_frequency_matrix[1,], xaxt="n", beside=TRUE, cex.names=0.55, las=3, col=c('blue'), main='Kruskal-Wallis - oneR - RFE')
axis(2,at=seq(0,50,5), lwd=2.0)
xlabs <- paste(colnames(kwo_RFE_RF_Freq_feature_frequency_matrix))
text(kwo_RFE_RF_Freq, par("usr")[3] - 0.025, srt = 45, adj = 1, 
     labels = xlabs, 
     xpd = TRUE, font = 2)


kwi_RFE_RF_Freq <- barplot(kwi_RFE_RF_Freq_feature_frequency_matrix[1,], xaxt="n", beside=TRUE, cex.names=0.55, las=3, col=c('sea green'), main='Kruskal-Wallis - InfoGain - RFE')
axis(2,at=seq(0,50,5), lwd=2.0)
xlabs <- paste(colnames(kwi_RFE_RF_Freq_feature_frequency_matrix))
text(kwi_RFE_RF_Freq, par("usr")[3] - 0.025, srt = 45, adj = 1, 
     labels = xlabs, 
     xpd = TRUE, font = 2)


kws_RFE_RF_Freq <- barplot(kws_RFE_RF_Freq_feature_frequency_matrix[1,], xaxt="n", beside=TRUE, cex.names=0.55, las=3, col=c('black'), main='Kruskal-Wallis - SymUncert - RFE')
axis(2,at=seq(0,50,5), lwd=2.0)
xlabs <- paste(colnames(kws_RFE_RF_Freq_feature_frequency_matrix))
text(kws_RFE_RF_Freq, par("usr")[3] - 0.025, srt = 45, adj = 1, 
     labels = xlabs, 
     xpd = TRUE, font = 2)



#----ROC-AUC Plots----#
par(pty='s')
plot(1-KWc_L_RF_Freq_new_specificities_1, KWc_L_RF_Freq_new_sensitivities_1, type='l', xlab='1-Specificity', ylab='Sensitivity', col='red', lwd=2.0, 
     main = 'Kruskal-Wallis - Chi-Sq - LASSO')
abline(coef=c(0,1),lwd=2.0)


par(pty='s')
plot(1-KWo_L_RF_Freq_new_specificities_2, KWo_L_RF_Freq_new_sensitivities_2, type='l', xlab='1-Specificity', ylab='Sensitivity', col='blue', lwd=2.0, 
     main = 'Kruskal-Wallis - oneR - LASSO')
abline(coef=c(0,1),lwd=2.0)



par(pty='s')
plot(1-kwi_L_RF_Freq_new_specificities_3, kwi_L_RF_Freq_new_sensitivities_3, type='l', xlab='1-Specificity', ylab='Sensitivity', col='sea green', lwd=2.0, 
     main = 'Kruskal-Wallis - InfoGain - LASSO')
abline(coef=c(0,1),lwd=2.0)


par(pty='s')
plot(1-kws_L_RF_Freq_new_specificities_4, kws_L_RF_Freq_new_sensitivities_4, type='l', xlab='1-Specificity', ylab='Sensitivity', col='black', lwd=2.0, 
     main = 'Kruskal-Wallis - SymUncert - LASSO')
abline(coef=c(0,1),lwd=2.0)


par(pty='s')
plot(1-kC_F_RF_Freq_new_specificities_1, kC_F_RF_Freq_new_sensitivities_1, type='l', xlab='1-Specificity', ylab='Sensitivity', col='red', lwd=2.0, 
     main = 'Kruskal-Wallis - Chi-Sq - Forward')
abline(coef=c(0,1),lwd=2.0)



par(pty='s')
plot(1-ko_F_RF_Freq_new_specificities_2, ko_F_RF_Freq_new_sensitivities_2, type='l', xlab='1-Specificity', ylab='Sensitivity', col='blue', lwd=2.0, 
     main = 'Kruskal-Wallis - oneR - Forward')
abline(coef=c(0,1),lwd=2.0)


par(pty='s')
plot(1-ki_F_RF_Freq_new_specificities_3, ki_F_RF_Freq_new_sensitivities_3, type='l', xlab='1-Specificity', ylab='Sensitivity', col='sea green', lwd=2.0, 
     main = 'Kruskal-Wallis - InfoGain - Forward')
abline(coef=c(0,1),lwd=2.0)

par(pty='s')
plot(1-ks_F_RF_Freq_new_specificities_4, ks_F_RF_Freq_new_sensitivities_4, type='l', xlab='1-Specificity', ylab='Sensitivity', col='black', lwd=2.0, 
     main = 'Kruskal-Wallis - SymUncert - Forward')
abline(coef=c(0,1),lwd=2.0)


par(pty='s')
plot(1-kwc_B_RF_Freq_new_specificities_1, kwc_B_RF_Freq_new_sensitivities_1, type='l', xlab='1-Specificity', ylab='Sensitivity', col='red', lwd=2.0, 
     main = 'Kruskal-Wallis - Chi-Sq - Backward')
abline(coef=c(0,1),lwd=2.0)

par(pty='s')
plot(1-kwi_B_RF_Freq_new_specificities_3, kwi_B_RF_Freq_new_sensitivities_3, type='l', xlab='1-Specificity', ylab='Sensitivity', col='sea green', lwd=2.0, 
     main = 'Kruskal-Wallis - InfoGain - Backward')
abline(coef=c(0,1),lwd=2.0)

par(pty='s')
plot(1-kwc_RFE_RF_Freq_new_specificities_1, kwc_RFE_RF_Freq_new_sensitivities_1, type='l', xlab='1-Specificity', ylab='Sensitivity', col='red', lwd=2.0, 
     main = 'Kruskal-Wallis - Chi-Sq - RFE')
abline(coef=c(0,1),lwd=2.0)

par(pty='s')
plot(1-kwo_B_RF_Freq_new_specificities_2, kwo_B_RF_Freq_new_sensitivities_2, type='l', xlab='1-Specificity', ylab='Sensitivity', col='blue', lwd=2.0, 
     main = 'Kruskal-Wallis - oneR - Backward')
abline(coef=c(0,1),lwd=2.0)

par(pty='s')
plot(1-kws_B_RF_Freq_new_specificities_4, kws_B_RF_Freq_new_sensitivities_4, type='l', xlab='1-Specificity', ylab='Sensitivity', col='black', lwd=2.0, 
     main = 'Kruskal-Wallis - SymUncert - Backward')
abline(coef=c(0,1),lwd=2.0)


par(pty='s')
plot(1-kwo_RFE_RF_Freq_new_specificities_2, kwo_RFE_RF_Freq_new_sensitivities_2, type='l', xlab='1-Specificity', ylab='Sensitivity', col='blue', lwd=2.0, 
     main = 'Kruskal-Wallis - oneR - RFE')
abline(coef=c(0,1),lwd=2.0)

par(pty='s')
plot(1-kwi_RFE_RF_Freq_new_specificities_3, kwi_RFE_RF_Freq_new_sensitivities_3, type='l', xlab='1-Specificity', ylab='Sensitivity', col='sea green', lwd=2.0, 
     main = 'Kruskal-Wallis - InfoGain - RFE')
abline(coef=c(0,1),lwd=2.0)


par(pty='s')
plot(1-kws_RFE_RF_Freq_new_specificities_4, kws_RFE_RF_Freq_new_sensitivities_4, type='l', xlab='1-Specificity', ylab='Sensitivity', col='black', lwd=2.0, 
     main = 'Kruskal-Wallis - SymUncert - RFE')
abline(coef=c(0,1),lwd=2.0)

#---Feature Stability plot for all models--#

feature_stability_clinical_vars <- data.frame(matrix(0,nrow=50,ncol=16))
feature_stability_clinical_vars[,1] = kw_chi_lasso_fs
feature_stability_clinical_vars[,2] = kw_oneR_lasso_fs
feature_stability_clinical_vars[,3] = kw_infogain_lasso_fs
feature_stability_clinical_vars[,4] = kw_symun_lasso_fs
feature_stability_clinical_vars[,5] = kw_chi_forward_fs
feature_stability_clinical_vars[,6] = kw_oneR_forward_fs
feature_stability_clinical_vars[,7] = kw_infogain_forward_fs
feature_stability_clinical_vars[,8] = kw_symun_forward_fs
feature_stability_clinical_vars[,9] = kw_chi_backward_fs
feature_stability_clinical_vars[,10] = kw_oneR_backward_fs
feature_stability_clinical_vars[,11] = kw_infogain_backward_fs
feature_stability_clinical_vars[,12] = kw_symun_backward_fs
feature_stability_clinical_vars[,13] = kw_chi_rfe_fs
feature_stability_clinical_vars[,14] = kw_oner_rfe_fs
feature_stability_clinical_vars[,15] = kw_infogain_rfe_fs
feature_stability_clinical_vars[,16] = kw_symun_rfe_fs

colnames(feature_stability_clinical_vars) = c('KW - Chi-Sq - LASSO','KW - oneR - LASSO','KW - InfoGain - LASSO','KW - SymUncert - LASSO','KW - Chi-Sq - Forward','KW - oneR - Forward','KW - InfoGain - Forward','KW - SymUncert - Forward','KW - Chi-Sq - Backward','KW - oneR - Backward','KW - InfoGain - Backward','KW - SymUncert - Backward','KW - Chi-Sq - RFE','KW - oneR - RFE','KW - InfoGain - RFE','KW - SymUncert - RFE')

ggplot(stack(feature_stability_clinical_vars), aes(x = ind, y = values)) +
  geom_boxplot( fill=rep(c('#F8766D','#619CFF','#00BA38','grey'),4)) +
  theme(axis.text.x = element_text(face = 'bold',angle = 90, vjust = 0.5, hjust=0.5, size=12), axis.text.y = element_text(face='bold', size=12), plot.title = element_text(hjust = 0.5), axis.title = element_text(face='bold')) +
  labs(title='Feature Stability Plot for Feature Selection Models on Clinical Variables' , x = 'Models') +
  scale_y_continuous(breaks = seq(0, 10, by = 1))