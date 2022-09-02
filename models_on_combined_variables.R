#------FEATURE PRE-FILTERING AND SELECTION FOR PREDICTION OF PATHOLOGICAL STAGE IN PROSTATE CANCER------#

#----Feature Selection and Pre-Filtering techniques on Combined Variables--#

#----Libraries----#\
library(caret)
library(pROC)
library(randomForest)
library(glmnet)
library(ROCR)

#----1. LASSO Feature Selection with Mutual Information and Chi-Squared----#
new_dat <- pros_clinical_mrna_dat_copy
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
feature_frequency_matrix = matrix(data=0,nrow = 1, ncol = ncol(new_dat)-1)
colnames(feature_frequency_matrix) = colnames(new_dat)[c(1:18576)]
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
  feature.matrix = feature.matrix[,-c(18577)]
  
  set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    pros_train = pc_train[i.train,]
    pros_valid = pc_train[-i.train,]
    
    pros_train_numeric <- pros_train[,c(1:18564,18567,18571,18577)]
    pros_train_categorical <- pros_train[,c(18565,18566,18568:18570,18572:18577)]
    numeric_features = c()
    categorical_features = c()
    features = c()
    
    for(col in 1:(ncol(pros_train_numeric)-1)){
      cns <- as.matrix(as.numeric(pros_train_numeric[,col]))
      disc <- data.frame(pros_train_numeric$PathStage)
      mi_test <- mmi(cns,disc)
      if(round(mi_test$mi,2) > 0.13){
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
    xm.train = model.matrix(PathStage~., data=pros_train)[,-1]
    xm.valid = model.matrix(PathStage~., data=pros_valid)[,-1]
    
    lam = cv.glmnet(xm.train, y.train, data=pros_train, alpha=1, family='binomial')
    lasso = glmnet(xm.train, y.train, data=pros_train, lambda=lam$lambda.min, alpha=1, family='binomial')
    lasso_features = which(coef(lasso)[-1] != 0)
    
    feature_stability_list = c(feature_stability_list, length(which(coef(lasso)[-1] != 0)))
    if(length(which(coef(lasso)[-1] != 0)) != 0){
      lasso_feature_names = colnames(xm.train)[lasso_features]
      print('LASSO features:')
      print(lasso_feature_names)
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

features_selected

features_with_frequency = c()
for(i in 1:ncol(feature_frequency_matrix)){
  if(feature_frequency_matrix[,i] > 0){
    features_with_frequency = c(features_with_frequency,colnames(feature_frequency_matrix)[i])
  }
}
feature_frequency_matrix_1 = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(feature_frequency_matrix_1) = features_with_frequency
feature_frequency_matrix_1

mic_L_RF_Freq_feature_frequency_matrix = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(mic_L_RF_Freq_feature_frequency_matrix) = features_with_frequency
mic_L_RF_Freq_feature_frequency_matrix


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

mic_L_RF_Freq_feature_frequency_matrix
mic_L_RF_Freq <- barplot(mic_L_RF_Freq_feature_frequency_matrix[1,], xaxt="n", beside=TRUE, cex.names=0.55, las=3, col=c('red'), main='Mutual Info - Chi-Sq - LASSO')
axis(2,at=seq(0,50,5), lwd=2.0)
xlabs <- paste(colnames(mic_L_RF_Freq_feature_frequency_matrix))
text(mic_L_RF_Freq, par("usr")[3] - 0.025, srt = 90, adj = 1, 
     labels = xlabs, 
     xpd = TRUE, font = 2,cex=0.7)


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


mic_L_RF_new_sensitivities_1 = new_sensitivities
mic_L_RF_new_specificities_1 = new_specificities




par(pty='s')
plot(1-mic_L_RF_new_specificities_1, mic_L_RF_new_sensitivities_1, type='l', xlab='1-Specificity', ylab='Sensitivity', col='red', lwd=2.0, 
     main = 'Mutual Info - Chi-Sq - LASSO')
abline(coef=c(0,1),lwd=2.0)

mi_chi_lasso_fs = feature_stability_list

#----2. LASSO Feature Selection with Mutual Information and One Rule----#

new_dat <- pros_clinical_mrna_dat_copy
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
feature_frequency_matrix = matrix(data=0,nrow = 1, ncol = ncol(new_dat)-1)
colnames(feature_frequency_matrix) = colnames(new_dat)[c(1:18576)]
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
  feature.matrix = feature.matrix[,-c(18577)]
  
  set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    pros_train = pc_train[i.train,]
    pros_valid = pc_train[-i.train,]
    
    pros_train_numeric <- pros_train[,c(1:18564,18567,18571,18577)]
    pros_train_categorical <- pros_train[,c(18565,18566,18568:18570,18572:18577)]
    numeric_features = c()
    categorical_features = c()
    features = c()
    
    for(col in 1:(ncol(pros_train_numeric)-1)){
      cns <- as.matrix(as.numeric(pros_train_numeric[,col]))
      disc <- data.frame(pros_train_numeric$PathStage)
      mi_test <- mmi(cns,disc)
      if(round(mi_test$mi,2) > 0.13){
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
      print('LASSO features:')
      print(lasso_feature_names)
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

features_selected

features_with_frequency = c()
for(i in 1:ncol(feature_frequency_matrix)){
  if(feature_frequency_matrix[,i] > 0){
    features_with_frequency = c(features_with_frequency,colnames(feature_frequency_matrix)[i])
  }
}
feature_frequency_matrix_2 = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(feature_frequency_matrix_2) = features_with_frequency
feature_frequency_matrix_2

mio_L_RF_Freq_feature_frequency_matrix = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(mio_L_RF_Freq_feature_frequency_matrix) = features_with_frequency
mio_L_RF_Freq_feature_frequency_matrix


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

mio_L_RF_Freq_feature_frequency_matrix
mio_L_RF_Freq <- barplot(mio_L_RF_Freq_feature_frequency_matrix[1,], xaxt="n", beside=TRUE, cex.names=0.55, las=3, col=c('blue'), main='Mutual Info - oneR - LASSO')
axis(2,at=seq(0,50,5), lwd=2.0)
xlabs <- paste(colnames(mio_L_RF_Freq_feature_frequency_matrix))
text(mio_L_RF_Freq, par("usr")[3] - 0.025, srt = 90, adj = 1, 
     labels = xlabs, 
     xpd = TRUE, font = 2,cex=0.7)



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


mio_L_RF_new_sensitivities_2 = new_sensitivities
mio_L_RF_new_specificities_2 = new_specificities




par(pty='s')
plot(1-mio_L_RF_new_specificities_2, mio_L_RF_new_sensitivities_2, type='l', xlab='1-Specificity', ylab='Sensitivity', col='blue', lwd=2.0, 
     main = 'Mutual Info - oneR - LASSO')
abline(coef=c(0,1),lwd=2.0)

mi_oneR_lasso_fs = feature_stability_list

#----3. LASSO Feature Selection with Mutual Information and Information Gain----#
new_dat <- pros_clinical_mrna_dat_copy
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
feature_frequency_matrix = matrix(data=0,nrow = 1, ncol = ncol(new_dat)-1)
colnames(feature_frequency_matrix) = colnames(new_dat)[c(1:18576)]
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
  feature.matrix = feature.matrix[,-c(18577)]
  
  set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    pros_train = pc_train[i.train,]
    pros_valid = pc_train[-i.train,]
    
    pros_train_numeric <- pros_train[,c(1:18564,18567,18571,18577)]
    pros_train_categorical <- pros_train[,c(18565,18566,18568:18570,18572:18577)]
    numeric_features = c()
    categorical_features = c()
    features = c()
    
    for(col in 1:(ncol(pros_train_numeric)-1)){
      cns <- as.matrix(as.numeric(pros_train_numeric[,col]))
      disc <- data.frame(pros_train_numeric$PathStage)
      mi_test <- mmi(cns,disc)
      if(round(mi_test$mi,2) > 0.13){
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
      print('LASSO features:')
      print(lasso_feature_names)
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

features_selected

features_with_frequency = c()
for(i in 1:ncol(feature_frequency_matrix)){
  if(feature_frequency_matrix[,i] > 0){
    features_with_frequency = c(features_with_frequency,colnames(feature_frequency_matrix)[i])
  }
}
feature_frequency_matrix_3 = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(feature_frequency_matrix_3) = features_with_frequency
feature_frequency_matrix_3

mii_L_RF_Freq_feature_frequency_matrix = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(mii_L_RF_Freq_feature_frequency_matrix) = features_with_frequency
mii_L_RF_Freq_feature_frequency_matrix


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

mii_L_RF_Freq_feature_frequency_matrix
mii_L_RF_Freq <- barplot(mii_L_RF_Freq_feature_frequency_matrix[1,], xaxt="n", beside=TRUE, cex.names=0.55, las=3, col=c('sea green'), main='Mutual Info - Info Gain - LASSO')
axis(2,at=seq(0,50,5), lwd=2.0)
xlabs <- paste(colnames(mii_L_RF_Freq_feature_frequency_matrix))
text(mii_L_RF_Freq, par("usr")[3] - 0.025, srt = 90, adj = 1, 
     labels = xlabs, 
     xpd = TRUE, font = 2,cex=0.7)



boxplot(feature_stability_list)


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


mii_L_RF_new_sensitivities_3 = new_sensitivities
mii_L_RF_new_specificities_3 = new_specificities




par(pty='s')
plot(1-mii_L_RF_new_specificities_3, mii_L_RF_new_sensitivities_3, type='l', xlab='1-Specificity', ylab='Sensitivity', col='sea green', lwd=2.0, 
     main = 'Mutual Info - Info Gain - LASSO')
abline(coef=c(0,1),lwd=2.0)

mi_infogain_lasso_fs = feature_stability_list

#----4. LASSO Feature Selection with Mutual Information and Symmetrical Uncertainty----#
new_dat <- pros_clinical_mrna_dat_copy
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
feature_frequency_matrix = matrix(data=0,nrow = 1, ncol = ncol(new_dat)-1)
colnames(feature_frequency_matrix) = colnames(new_dat)[c(1:18576)]
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
  feature.matrix = feature.matrix[,-c(18577)]
  
  set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    pros_train = pc_train[i.train,]
    pros_valid = pc_train[-i.train,]
    
    pros_train_numeric <- pros_train[,c(1:18564,18567,18571,18577)]
    pros_train_categorical <- pros_train[,c(18565,18566,18568:18570,18572:18577)]
    numeric_features = c()
    categorical_features = c()
    features = c()
    
    for(col in 1:(ncol(pros_train_numeric)-1)){
      cns <- as.matrix(as.numeric(pros_train_numeric[,col]))
      disc <- data.frame(pros_train_numeric$PathStage)
      mi_test <- mmi(cns,disc)
      if(round(mi_test$mi,2) > 0.13){
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
      print('LASSO features:')
      print(lasso_feature_names)
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

features_selected

features_with_frequency = c()
for(i in 1:ncol(feature_frequency_matrix)){
  if(feature_frequency_matrix[,i] > 0){
    features_with_frequency = c(features_with_frequency,colnames(feature_frequency_matrix)[i])
  }
}
feature_frequency_matrix_4 = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(feature_frequency_matrix_4) = features_with_frequency
feature_frequency_matrix_4

mis_L_RF_Freq_feature_frequency_matrix = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(mis_L_RF_Freq_feature_frequency_matrix) = features_with_frequency
mis_L_RF_Freq_feature_frequency_matrix


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

mis_L_RF_Freq_feature_frequency_matrix
mis_L_RF_Freq <- barplot(feature_frequency_matrix_4[1,], xaxt="n", beside=TRUE, cex.names=0.55, las=3, col=c('black'), main='Mutual Info - Sym Uncert - LASSO')
axis(2,at=seq(0,50,5), lwd=2.0)
xlabs <- paste(colnames(feature_frequency_matrix_4))
text(mis_L_RF_Freq, par("usr")[3] - 0.025, srt = 90, adj = 1, 
     labels = xlabs, 
     xpd = TRUE, font = 2,cex=0.7)


boxplot(feature_stability_list)


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


mis_L_RF_new_sensitivities_4 = new_sensitivities
mis_L_RF_new_specificities_4 = new_specificities




par(pty='s')
plot(1-mis_L_RF_new_specificities_4, mis_L_RF_new_sensitivities_4, type='l', xlab='1-Specificity', ylab='Sensitivity', col='black', lwd=2.0, 
     main = 'Mutual Info - Sym Uncert - LASSO')
abline(coef=c(0,1),lwd=2.0)

mi_symun_lasso_fs = feature_stability_list

#----5. Forward Stepwise Feature Selection with Mutual Information and Chi-Squared----#
new_dat <- pros_clinical_mrna_dat_copy
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
feature_frequency_matrix = matrix(data=0,nrow = 1, ncol = ncol(new_dat)-1)
colnames(feature_frequency_matrix) = colnames(new_dat)[c(1:18576)]
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
  feature.matrix = feature.matrix[,-c(18577)]
  
  set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    pros_train = pc_train[i.train,]
    pros_valid = pc_train[-i.train,]
    
    pros_train_numeric <- pros_train[,c(1:18564,18567,18571,18577)]
    pros_train_categorical <- pros_train[,c(18565,18566,18568:18570,18572:18577)]
    numeric_features = c()
    categorical_features = c()
    features = c()
    
    for(col in 1:(ncol(pros_train_numeric)-1)){
      cns <- as.matrix(as.numeric(pros_train_numeric[,col]))
      disc <- data.frame(pros_train_numeric$PathStage)
      mi_test <- mmi(cns,disc)
      if(round(mi_test$mi,2) > 0.13){
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
    
    #print(pros_train)
    
    reg.fwd = regsubsets(PathStage~., data=pros_train, method="forward")
    R2adj = summary(reg.fwd)$adjr2
    R2adj.index = which.max(R2adj)
    mod = summary(reg.fwd)$which[R2adj.index,]
    vars = rownames(data.frame(which(mod)[-1]))
    print("Forward Features:")
    print(vars)
    
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
feature_frequency_matrix_1 = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(feature_frequency_matrix_1) = features_with_frequency

mic_F_RF_Freq_feature_frequency_matrix = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(mic_F_RF_Freq_feature_frequency_matrix) = features_with_frequency

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



mic_F_RF_Freq_feature_frequency_matrix
mic_F_RF_Freq <- barplot(mic_F_RF_Freq_feature_frequency_matrix[1,], xaxt="n", beside=TRUE, cex.names=0.55, las=3, col=c('red'), main='Mutual Info - Chi-Sq - Forward')
axis(2,at=seq(0,50,5))
xlabs <- paste(colnames(mic_F_RF_Freq_feature_frequency_matrix))
text(mic_F_RF_Freq, par("usr")[3] - 0.025, srt = 90, adj = 1, 
     labels = xlabs, 
     xpd = TRUE, font = 2, cex=0.5)



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


mic_F_RF_Freq_new_sensitivities_1 = new_sensitivities
mic_F_RF_Freq_new_specificities_1 = new_specificities



par(pty='s')
plot(1-mic_F_RF_Freq_new_specificities_1, mic_F_RF_Freq_new_sensitivities_1, type='l', xlab='1-Specificity', ylab='Sensitivity', col='red', lwd=2.0, 
     main = 'Mutual Info - Chi-Sq - Forward')
abline(coef=c(0,1),lwd=0.5)

mi_chi_forward_fs = feature_stability_list

#----6. Forward Stepwise Feature Selection with Mutual Information and One Rule----#
new_dat <- pros_clinical_mrna_dat_copy
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
feature_frequency_matrix = matrix(data=0,nrow = 1, ncol = ncol(new_dat)-1)
colnames(feature_frequency_matrix) = colnames(new_dat)[c(1:18576)]
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
  feature.matrix = feature.matrix[,-c(18577)]
  
  set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    pros_train = pc_train[i.train,]
    pros_valid = pc_train[-i.train,]
    
    pros_train_numeric <- pros_train[,c(1:18564,18567,18571,18577)]
    pros_train_categorical <- pros_train[,c(18565,18566,18568:18570,18572:18577)]
    numeric_features = c()
    categorical_features = c()
    features = c()
    
    for(col in 1:(ncol(pros_train_numeric)-1)){
      cns <- as.matrix(as.numeric(pros_train_numeric[,col]))
      disc <- data.frame(pros_train_numeric$PathStage)
      mi_test <- mmi(cns,disc)
      if(round(mi_test$mi,2) > 0.13){
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
    
    #print(pros_train)
    
    reg.fwd = regsubsets(PathStage~., data=pros_train, method="forward")
    R2adj = summary(reg.fwd)$adjr2
    R2adj.index = which.max(R2adj)
    mod = summary(reg.fwd)$which[R2adj.index,]
    vars = rownames(data.frame(which(mod)[-1]))
    print("Forward Features:")
    print(vars)
    
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
feature_frequency_matrix_2 = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(feature_frequency_matrix_2) = features_with_frequency
feature_frequency_matrix_2

mio_F_RF_Freq_feature_frequency_matrix = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(mio_F_RF_Freq_feature_frequency_matrix) = features_with_frequency
mio_F_RF_Freq_feature_frequency_matrix

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



mio_F_RF_Freq_feature_frequency_matrix
mio_F_RF_Freq <- barplot(mio_F_RF_Freq_feature_frequency_matrix[1,], xaxt="n", beside=TRUE, cex.names=0.55, las=3, col=c('blue'), main='Mutual Info - OneR - Forward')
axis(2,at=seq(0,50,5))
xlabs <- paste(colnames(mio_F_RF_Freq_feature_frequency_matrix))
text(mio_F_RF_Freq, par("usr")[3] - 0.025, srt = 90, adj = 1, 
     labels = xlabs, 
     xpd = TRUE, font = 2, cex=0.5)


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


mio_F_RF_Freq_new_sensitivities_2 = new_sensitivities
mio_F_RF_Freq_new_specificities_2 = new_specificities



par(pty='s')
plot(1-mio_F_RF_Freq_new_specificities_2, mio_F_RF_Freq_new_sensitivities_2, type='l', xlab='Specificity', ylab='Sensitivity', col='blue', lwd=2.0, 
     main = 'Mutual Info - oneR - Forward')
abline(coef=c(0,1),lwd=0.5)

mi_oneR_forward_fs = feature_stability_list

#----7. Forward Stepwise Feature Selection with Mutual Information and Information Gain----#
new_dat <- pros_clinical_mrna_dat_copy
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
feature_frequency_matrix = matrix(data=0,nrow = 1, ncol = ncol(new_dat)-1)
colnames(feature_frequency_matrix) = colnames(new_dat)[c(1:18576)]
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
  feature.matrix = feature.matrix[,-c(18577)]
  
  set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    pros_train = pc_train[i.train,]
    pros_valid = pc_train[-i.train,]
    
    pros_train_numeric <- pros_train[,c(1:18564,18567,18571,18577)]
    pros_train_categorical <- pros_train[,c(18565,18566,18568:18570,18572:18577)]
    numeric_features = c()
    categorical_features = c()
    features = c()
    
    for(col in 1:(ncol(pros_train_numeric)-1)){
      cns <- as.matrix(as.numeric(pros_train_numeric[,col]))
      disc <- data.frame(pros_train_numeric$PathStage)
      mi_test <- mmi(cns,disc)
      if(round(mi_test$mi,2) > 0.13){
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
    
    #print(pros_train)
    
    reg.fwd = regsubsets(PathStage~., data=pros_train, method="forward")
    R2adj = summary(reg.fwd)$adjr2
    R2adj.index = which.max(R2adj)
    mod = summary(reg.fwd)$which[R2adj.index,]
    vars = rownames(data.frame(which(mod)[-1]))
    print("Forward Features:")
    print(vars)
    
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
feature_frequency_matrix_3 = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(feature_frequency_matrix_3) = features_with_frequency
feature_frequency_matrix_3

mii_F_RF_Freq_feature_frequency_matrix = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(mii_F_RF_Freq_feature_frequency_matrix) = features_with_frequency
mii_F_RF_Freq_feature_frequency_matrix

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



mii_F_RF_Freq_feature_frequency_matrix
mii_F_RF_Freq <- barplot(mii_F_RF_Freq_feature_frequency_matrix[1,], xaxt="n", beside=TRUE, cex.names=0.55, las=3, col=c('sea green'), main='Mutual Info - InfoGain - Forward')
axis(2,at=seq(0,50,5))
xlabs <- paste(colnames(mii_F_RF_Freq_feature_frequency_matrix))
text(mii_F_RF_Freq, par("usr")[3] - 0.025, srt = 90, adj = 1, 
     labels = xlabs, 
     xpd = TRUE, font = 2, cex=0.5)


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

mii_F_RF_Freq_new_sensitivities_3 = new_sensitivities
mii_F_RF_Freq_new_specificities_3 = new_specificities



par(pty='s')
plot(1-mii_F_RF_Freq_new_specificities_3, mii_F_RF_Freq_new_sensitivities_3, type='l', xlab='1-Specificities', ylab='Sensitivities', col='sea green', lwd=2.0, 
     main = 'Mutual Info - InfoGain - Forward')
abline(coef=c(0,1),lwd=0.5)

mi_infogain_forward_fs = feature_stability_list

#----8. Forward Stepwise Feature Selection with Mutual Information and Symmetrical Uncertainty----#
new_dat <- pros_clinical_mrna_dat_copy
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
feature_frequency_matrix = matrix(data=0,nrow = 1, ncol = ncol(new_dat)-1)
colnames(feature_frequency_matrix) = colnames(new_dat)[c(1:18576)]
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
  feature.matrix = feature.matrix[,-c(18577)]
  
  set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    pros_train = pc_train[i.train,]
    pros_valid = pc_train[-i.train,]
    
    pros_train_numeric <- pros_train[,c(1:18564,18567,18571,18577)]
    pros_train_categorical <- pros_train[,c(18565,18566,18568:18570,18572:18577)]
    numeric_features = c()
    categorical_features = c()
    features = c()
    
    for(col in 1:(ncol(pros_train_numeric)-1)){
      cns <- as.matrix(as.numeric(pros_train_numeric[,col]))
      disc <- data.frame(pros_train_numeric$PathStage)
      mi_test <- mmi(cns,disc)
      if(round(mi_test$mi,2) > 0.13){
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
    
    #print(pros_train)
    
    reg.fwd = regsubsets(PathStage~., data=pros_train, method="forward")
    R2adj = summary(reg.fwd)$adjr2
    R2adj.index = which.max(R2adj)
    mod = summary(reg.fwd)$which[R2adj.index,]
    vars = rownames(data.frame(which(mod)[-1]))
    print("Forward Features:")
    print(vars)
    
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
feature_frequency_matrix_4 = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(feature_frequency_matrix_4) = features_with_frequency
feature_frequency_matrix_4

mis_F_RF_Freq_feature_frequency_matrix = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(mis_F_RF_Freq_feature_frequency_matrix) = features_with_frequency
mis_F_RF_Freq_feature_frequency_matrix

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



mis_F_RF_Freq_feature_frequency_matrix
mis_F_RF_Freq <- barplot(mis_F_RF_Freq_feature_frequency_matrix[1,], xaxt="n", beside=TRUE, cex.names=0.55, las=3, col=c('black'), main='Mutual Info - SymUncert - Forward')
axis(2,at=seq(0,50,5))
xlabs <- paste(colnames(mis_F_RF_Freq_feature_frequency_matrix))
text(mis_F_RF_Freq, par("usr")[3] - 0.025, srt = 90, adj = 1, 
     labels = xlabs, 
     xpd = TRUE, font = 2, cex=0.5)



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

mis_F_RF_Freq_new_sensitivities_4 = new_sensitivities
mis_F_RF_Freq_new_specificities_4 = new_specificities


par(pty='s')
plot(1-mis_F_RF_Freq_new_specificities_4, mis_F_RF_Freq_new_sensitivities_4, type='l', xlab='1-Specificities', ylab='Sensitivities', col='black', lwd=2.0, 
     main = 'Mutual Info - SymUncert - Forward')
abline(coef=c(0,1),lwd=0.5)

mi_symun_forward_fs = feature_stability_list

#----9. Backward Stepwise Feature Selection with Mutual Information and Chi-Squared----#
new_dat <- pros_clinical_mrna_dat_copy
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
feature_frequency_matrix = matrix(data=0,nrow = 1, ncol = ncol(new_dat)-1)
colnames(feature_frequency_matrix) = colnames(new_dat)[c(1:18576)]
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
  feature.matrix = feature.matrix[,-c(18577)]
  
  set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    pros_train = pc_train[i.train,]
    pros_valid = pc_train[-i.train,]
    
    pros_train_numeric <- pros_train[,c(1:18564,18567,18571,18577)]
    pros_train_categorical <- pros_train[,c(18565,18566,18568:18570,18572:18577)]
    numeric_features = c()
    categorical_features = c()
    features = c()
    
    for(col in 1:(ncol(pros_train_numeric)-1)){
      cns <- as.matrix(as.numeric(pros_train_numeric[,col]))
      disc <- data.frame(pros_train_numeric$PathStage)
      mi_test <- mmi(cns,disc)
      if(round(mi_test$mi,2) > 0.13){
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
    print("Backward Features:")
    print(vars)
    
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
feature_frequency_matrix_1 = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(feature_frequency_matrix_1) = features_with_frequency
feature_frequency_matrix_1

mic_B_RF_Freq_feature_frequency_matrix = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(mic_B_RF_Freq_feature_frequency_matrix) = features_with_frequency
mic_B_RF_Freq_feature_frequency_matrix




mic_B_RF_Freq_feature_frequency_matrix
mic_B_RF_Freq <- barplot(mic_B_RF_Freq_feature_frequency_matrix[1,], xaxt="n", beside=TRUE, cex.names=0.55, las=3, col=c('red'), main='Mutual Info - Chi-Sq - Backward')
axis(2,at=seq(0,50,5))
xlabs <- paste(colnames(mic_B_RF_Freq_feature_frequency_matrix))
text(mic_B_RF_Freq, par("usr")[3] - 0.025, srt = 90, adj = 1, 
     labels = xlabs, 
     xpd = TRUE, font = 2, cex=0.5)


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


mic_B_RF_Freq_new_sensitivities_1 = new_sensitivities
mic_B_RF_Freq_new_specificities_1 = new_specificities



par(pty='s')
plot(1-mic_B_RF_Freq_new_specificities_1, mic_B_RF_Freq_new_sensitivities_1, type='l', xlab='1-Specificity', ylab='Sensitivity', col='red', lwd=2.0, 
     main = 'Mutual Info - Chi-Sq - Backward')
abline(coef=c(0,1),lwd=0.5)

mi_chi_backward_fs = feature_stability_list

#----10. Backward Stepwise Feature Selection with Mutual Information and One Rule----#
new_dat <- pros_clinical_mrna_dat_copy
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
feature_frequency_matrix = matrix(data=0,nrow = 1, ncol = ncol(new_dat)-1)
colnames(feature_frequency_matrix) = colnames(new_dat)[c(1:18576)]
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
  feature.matrix = feature.matrix[,-c(18577)]
  
  set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    pros_train = pc_train[i.train,]
    pros_valid = pc_train[-i.train,]
    
    pros_train_numeric <- pros_train[,c(1:18564,18567,18571,18577)]
    pros_train_categorical <- pros_train[,c(18565,18566,18568:18570,18572:18577)]
    numeric_features = c()
    categorical_features = c()
    features = c()
    
    for(col in 1:(ncol(pros_train_numeric)-1)){
      cns <- as.matrix(as.numeric(pros_train_numeric[,col]))
      disc <- data.frame(pros_train_numeric$PathStage)
      mi_test <- mmi(cns,disc)
      if(round(mi_test$mi,2) > 0.13){
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
    
    #print(pros_train)
    
    reg.fwd = regsubsets(PathStage~., data=pros_train, method="backward")
    R2adj = summary(reg.fwd)$adjr2
    R2adj.index = which.max(R2adj)
    mod = summary(reg.fwd)$which[R2adj.index,]
    vars = rownames(data.frame(which(mod)[-1]))
    print("Backward Features:")
    print(vars)
    
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
feature_frequency_matrix_2 = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(feature_frequency_matrix_2) = features_with_frequency
feature_frequency_matrix_2

mio_B_RF_Freq_feature_frequency_matrix = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(mio_B_RF_Freq_feature_frequency_matrix) = features_with_frequency
mio_B_RF_Freq_feature_frequency_matrix




mio_B_RF_Freq_feature_frequency_matrix
mio_B_RF_Freq <- barplot(mio_B_RF_Freq_feature_frequency_matrix[1,], xaxt="n", beside=TRUE, cex.names=0.55, las=3, col=c('blue'), main='Mutual Info - OneR - Backward')
axis(2,at=seq(0,50,5))
xlabs <- paste(colnames(mio_B_RF_Freq_feature_frequency_matrix))
text(mio_B_RF_Freq, par("usr")[3] - 0.025, srt = 90, adj = 1, 
     labels = xlabs, 
     xpd = TRUE, font = 2, cex=0.5)



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


mio_B_RF_Freq_new_sensitivities_2 = new_sensitivities
mio_B_RF_Freq_new_specificities_2 = new_specificities



par(pty='s')
plot(1-mio_B_RF_Freq_new_specificities_2, mio_B_RF_Freq_new_sensitivities_2, type='l', xlab='1-Specificity', ylab='Sensitivity', col='blue', lwd=2.0, 
     main = 'Mutual Info - oneR - Backward')
abline(coef=c(0,1),lwd=0.5)

mi_oneR_backward_fs = feature_stability_list


#----11. Backward Stepwise Feature Selection with Mutual Information and Information Gain----#
new_dat <- pros_clinical_mrna_dat_copy
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
feature_frequency_matrix = matrix(data=0,nrow = 1, ncol = ncol(new_dat)-1)
colnames(feature_frequency_matrix) = colnames(new_dat)[c(1:18576)]
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
  feature.matrix = feature.matrix[,-c(18577)]
  
  set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    pros_train = pc_train[i.train,]
    pros_valid = pc_train[-i.train,]
    
    pros_train_numeric <- pros_train[,c(1:18564,18567,18571,18577)]
    pros_train_categorical <- pros_train[,c(18565,18566,18568:18570,18572:18577)]
    numeric_features = c()
    categorical_features = c()
    features = c()
    
    for(col in 1:(ncol(pros_train_numeric)-1)){
      cns <- as.matrix(as.numeric(pros_train_numeric[,col]))
      disc <- data.frame(pros_train_numeric$PathStage)
      mi_test <- mmi(cns,disc)
      if(round(mi_test$mi,2) > 0.13){
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
    print("Backward Features:")
    print(vars)
    
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
feature_frequency_matrix_3 = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(feature_frequency_matrix_3) = features_with_frequency
feature_frequency_matrix_3

mii_B_RF_Freq_feature_frequency_matrix = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(mii_B_RF_Freq_feature_frequency_matrix) = features_with_frequency
mii_B_RF_Freq_feature_frequency_matrix



mii_B_RF_Freq_feature_frequency_matrix
mii_B_RF_Freq <- barplot(mii_B_RF_Freq_feature_frequency_matrix[1,], xaxt="n", beside=TRUE, cex.names=0.55, las=3, col=c('sea green'), main='Mutual Info - InfoGain - Backward')
axis(2,at=seq(0,50,5))
xlabs <- paste(colnames(mii_B_RF_Freq_feature_frequency_matrix))
text(mii_B_RF_Freq, par("usr")[3] - 0.025, srt = 90, adj = 1, 
     labels = xlabs, 
     xpd = TRUE, font = 2, cex=0.5)



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

mii_B_RF_Freq_new_sensitivities_3 = new_sensitivities
mii_B_RF_Freq_new_specificities_3 = new_specificities



par(pty='s')
plot(1-mii_B_RF_Freq_new_specificities_3, mii_B_RF_Freq_new_sensitivities_3, type='l', xlab='1-Specificity', ylab='Sensitivity', col='sea green', lwd=2.0, 
     main = 'Mutual Info - InfoGain - Backward')
abline(coef=c(0,1),lwd=0.5)

mi_infogain_backward_fs = feature_stability_list

#----12. Backward Stepwise Feature Selection with Mutual Information and Symmetrical Uncertainty----#
new_dat <- pros_clinical_mrna_dat_copy
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
feature_frequency_matrix = matrix(data=0,nrow = 1, ncol = ncol(new_dat)-1)
colnames(feature_frequency_matrix) = colnames(new_dat)[c(1:18576)]
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
  feature.matrix = feature.matrix[,-c(18577)]
  
  set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    pros_train = pc_train[i.train,]
    pros_valid = pc_train[-i.train,]
    
    pros_train_numeric <- pros_train[,c(1:18564,18567,18571,18577)]
    pros_train_categorical <- pros_train[,c(18565,18566,18568:18570,18572:18577)]
    numeric_features = c()
    categorical_features = c()
    features = c()
    
    for(col in 1:(ncol(pros_train_numeric)-1)){
      cns <- as.matrix(as.numeric(pros_train_numeric[,col]))
      disc <- data.frame(pros_train_numeric$PathStage)
      mi_test <- mmi(cns,disc)
      if(round(mi_test$mi,2) > 0.13){
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
    
    print("Numerical Features:")
    print(numeric_features)
    print("Categorical Feature:")
    print(categorical_features)
    
    features = c(numeric_features, categorical_features)
    features = c(categorical_features, numeric_features,'PathStage')
    pros_train = pros_train[,features]
    pros_valid = pros_valid[,features]
    
    #print(pros_train)
    
    reg.fwd = regsubsets(PathStage~., data=pros_train, method="backward")
    R2adj = summary(reg.fwd)$adjr2
    R2adj.index = which.max(R2adj)
    mod = summary(reg.fwd)$which[R2adj.index,]
    vars = rownames(data.frame(which(mod)[-1]))
    print("Backward Features:")
    print(vars)
    
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
feature_frequency_matrix_4 = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(feature_frequency_matrix_4) = features_with_frequency
feature_frequency_matrix_4

mis_B_RF_Freq_feature_frequency_matrix = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(mis_B_RF_Freq_feature_frequency_matrix) = features_with_frequency
mis_B_RF_Freq_feature_frequency_matrix



mis_B_RF_Freq_feature_frequency_matrix
mis_B_RF_Freq <- barplot(mis_B_RF_Freq_feature_frequency_matrix[1,], xaxt="n", beside=TRUE, cex.names=0.55, las=3, col=c('black'), main='Mutual Info - SymUncert - Backward')
axis(2,at=seq(0,50,5))
xlabs <- paste(colnames(mis_B_RF_Freq_feature_frequency_matrix))
text(mis_B_RF_Freq, par("usr")[3] - 0.025, srt = 90, adj = 1, 
     labels = xlabs, 
     xpd = TRUE, font = 2, cex=0.5)



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

mis_B_RF_Freq_new_sensitivities_4 = new_sensitivities
mis_B_RF_Freq_new_specificities_4 = new_specificities



par(pty='s')
plot(1-mis_B_RF_Freq_new_specificities_4, mis_B_RF_Freq_new_sensitivities_4, type='l', xlab='1-Specificity', ylab='Sensitivity', col='black', lwd=2.0, 
     main = 'Mutual Info - SymUncert - Backward')
abline(coef=c(0,1),lwd=0.5)

mi_symun_backward_fs = feature_stability_list

#----13. RFE Feature Selection with Mutual Information and Chi-Squared----#
new_dat <- pros_clinical_mrna_dat_copy
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
feature_frequency_matrix = matrix(data=0,nrow = 1, ncol = ncol(new_dat)-1)
colnames(feature_frequency_matrix) = colnames(new_dat)[c(1:18576)]
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
  feature.matrix = feature.matrix[,-c(18577)]
  
  set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    pros_train = pc_train[i.train,]
    pros_valid = pc_train[-i.train,]
    
    pros_train_numeric <- pros_train[,c(1:18564,18567,18571,18577)]
    pros_train_categorical <- pros_train[,c(8565,18566,18568:18570,18572:18577)]
    numeric_features = c()
    categorical_features = c()
    features = c()
    
    for(col in 1:(ncol(pros_train_numeric)-1)){
      cns <- as.matrix(as.numeric(pros_train_numeric[,col]))
      disc <- data.frame(pros_train_numeric$PathStage)
      mi_test <- mmi(cns,disc)
      if(round(mi_test$mi,2) > 0.13){
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
    
    print('RFE Features:')
    print(rf.rfe$optVariables)
    
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
feature_frequency_matrix_1 = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(feature_frequency_matrix_1) = features_with_frequency
feature_frequency_matrix_1

mic_RFE_RF_Freq_feature_frequency_matrix = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(mic_RFE_RF_Freq_feature_frequency_matrix) = features_with_frequency
mic_RFE_RF_Freq_feature_frequency_matrix



mic_RFE_RF_Freq_feature_frequency_matrix
mic_RFE_RF_Freq <- barplot(mic_RFE_RF_Freq_feature_frequency_matrix[1,], xaxt="n", beside=TRUE, cex.names=0.55, las=3, col=c('red'), main='Mutual Info - Chi-Sq - RFE')
axis(2,at=seq(0,50,5))
xlabs <- paste(colnames(mic_RFE_RF_Freq_feature_frequency_matrix))
text(mic_RFE_RF_Freq, par("usr")[3] - 0.025, srt = 90, adj = 1, 
     labels = xlabs, 
     xpd = TRUE, font = 2, cex=0.6)


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


mic_RFE_RF_Freq_new_sensitivities_1 = new_sensitivities
mic_RFE_RF_Freq_new_specificities_1 = new_specificities



par(pty='s')
plot(1-mic_RFE_RF_Freq_new_specificities_1, mic_RFE_RF_Freq_new_sensitivities_1, type='l', xlab='1-Specificity', ylab='Sensitivity', col='red', lwd=2.0, 
     main = 'Mutual Info - Chi-Sq - RFE')
abline(coef=c(0,1),lwd=0.5)

mi_chi_rfe_fs = feature_stability_list

#----14. RFE Feature Selection with Mutual Information and One Rule----#
new_dat <- pros_clinical_mrna_dat_copy
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
feature_frequency_matrix = matrix(data=0,nrow = 1, ncol = ncol(new_dat)-1)
colnames(feature_frequency_matrix) = colnames(new_dat)[c(1:18576)]
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
  feature.matrix = feature.matrix[,-c(18577)]
  
  set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    pros_train = pc_train[i.train,]
    pros_valid = pc_train[-i.train,]
    
    pros_train_numeric <- pros_train[,c(1:18564,18567,18571,18577)]
    pros_train_categorical <- pros_train[,c(8565,18566,18568:18570,18572:18577)]
    numeric_features = c()
    categorical_features = c()
    features = c()
    
    for(col in 1:(ncol(pros_train_numeric)-1)){
      cns <- as.matrix(as.numeric(pros_train_numeric[,col]))
      disc <- data.frame(pros_train_numeric$PathStage)
      mi_test <- mmi(cns,disc)
      if(round(mi_test$mi,2) > 0.13){
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
    
    print('RFE Features:')
    print(rf.rfe$optVariables)
    
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
feature_frequency_matrix_2 = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(feature_frequency_matrix_2) = features_with_frequency
feature_frequency_matrix_2


mio_RFE_RF_Freq_feature_frequency_matrix = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(mio_RFE_RF_Freq_feature_frequency_matrix) = features_with_frequency
mio_RFE_RF_Freq_feature_frequency_matrix


mio_RFE_RF_Freq_feature_frequency_matrix
mio_RFE_RF_Freq <- barplot(mio_RFE_RF_Freq_feature_frequency_matrix[1,], xaxt="n", beside=TRUE, cex.names=0.55, las=3, col=c('blue'), main='Mutual Info - oneR - RFE')
axis(2,at=seq(0,50,5))
xlabs <- paste(colnames(mio_RFE_RF_Freq_feature_frequency_matrix))
text(mio_RFE_RF_Freq, par("usr")[3] - 0.025, srt = 90, adj = 1, 
     labels = xlabs, 
     xpd = TRUE, font = 2, cex=0.5)


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


mio_RFE_RF_Freq_new_sensitivities_2 = new_sensitivities
mio_RFE_RF_Freq_new_specificities_2 = new_specificities



par(pty='s')
plot(1-mio_RFE_RF_Freq_new_specificities_2, mio_RFE_RF_Freq_new_sensitivities_2, type='l', xlab='1-Specificity', ylab='Sensitivity', col='blue', lwd=2.0, 
     main = 'Mutual Info - oneR - RFE')
abline(coef=c(0,1),lwd=0.5)

mi_oneR_rfe_fs = feature_stability_list

#----15. RFE Feature Selection with Mutual Information and Information Gain----#
new_dat <- pros_clinical_mrna_dat_copy
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
feature_frequency_matrix = matrix(data=0,nrow = 1, ncol = ncol(new_dat)-1)
colnames(feature_frequency_matrix) = colnames(new_dat)[c(1:18576)]
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
  feature.matrix = feature.matrix[,-c(18577)]
  
  set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    pros_train = pc_train[i.train,]
    pros_valid = pc_train[-i.train,]
    
    pros_train_numeric <- pros_train[,c(1:18564,18567,18571,18577)]
    pros_train_categorical <- pros_train[,c(8565,18566,18568:18570,18572:18577)]
    numeric_features = c()
    categorical_features = c()
    features = c()
    
    for(col in 1:(ncol(pros_train_numeric)-1)){
      cns <- as.matrix(as.numeric(pros_train_numeric[,col]))
      disc <- data.frame(pros_train_numeric$PathStage)
      mi_test <- mmi(cns,disc)
      if(round(mi_test$mi,2) > 0.13){
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
    
    print('RFE Features:')
    print(rf.rfe$optVariables)
    
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
feature_frequency_matrix_3 = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(feature_frequency_matrix_3) = features_with_frequency
feature_frequency_matrix_3

mii_RFE_RF_Freq_feature_frequency_matrix = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(mii_RFE_RF_Freq_feature_frequency_matrix) = features_with_frequency
mii_RFE_RF_Freq_feature_frequency_matrix


mii_RFE_RF_Freq_feature_frequency_matrix
mii_RFE_RF_Freq <- barplot(mii_RFE_RF_Freq_feature_frequency_matrix[1,], xaxt="n", beside=TRUE, cex.names=0.55, las=3, col=c('sea green'), main='Mutual Info - Info Gain - RFE')
axis(2,at=seq(0,50,5))
xlabs <- paste(colnames(mii_RFE_RF_Freq_feature_frequency_matrix))
text(mii_RFE_RF_Freq, par("usr")[3] - 0.025, srt = 90, adj = 1, 
     labels = xlabs, 
     xpd = TRUE, font = 2, cex=0.6)



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


mii_RFE_RF_Freq_new_sensitivities_3 = new_sensitivities
mii_RFE_RF_Freq_new_specificities_3 = new_specificities



par(pty='s')
plot(1-mii_RFE_RF_Freq_new_specificities_3, mii_RFE_RF_Freq_new_sensitivities_3, type='l', xlab='1-Specificity', ylab='Sensitivity', col='sea green', lwd=2.0, 
     main = 'Mutual Info - Info Gain - RFE')
abline(coef=c(0,1),lwd=0.5)

mi_infogain_rfe_fs = feature_stability_list

#----16. RFE Feature Selection with Mutual Information and Symmetrical Uncertainty----#
new_dat <- pros_clinical_mrna_dat_copy
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
feature_frequency_matrix = matrix(data=0,nrow = 1, ncol = ncol(new_dat)-1)
colnames(feature_frequency_matrix) = colnames(new_dat)[c(1:18576)]
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
  feature.matrix = feature.matrix[,-c(18577)]
  
  set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    pros_train = pc_train[i.train,]
    pros_valid = pc_train[-i.train,]
    
    pros_train_numeric <- pros_train[,c(1:18564,18567,18571,18577)]
    pros_train_categorical <- pros_train[,c(8565,18566,18568:18570,18572:18577)]
    numeric_features = c()
    categorical_features = c()
    features = c()
    
    for(col in 1:(ncol(pros_train_numeric)-1)){
      cns <- as.matrix(as.numeric(pros_train_numeric[,col]))
      disc <- data.frame(pros_train_numeric$PathStage)
      mi_test <- mmi(cns,disc)
      if(round(mi_test$mi,2) > 0.13){
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
    
    print('RFE Features:')
    print(rf.rfe$optVariables)
    
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
feature_frequency_matrix_4 = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(feature_frequency_matrix_4) = features_with_frequency
feature_frequency_matrix_4

mis_RFE_RF_Freq_feature_frequency_matrix = matrix(unlist(feature_frequency_matrix[,features_with_frequency]),nrow=1, ncol=length(features_with_frequency))
colnames(mis_RFE_RF_Freq_feature_frequency_matrix) = features_with_frequency
mis_RFE_RF_Freq_feature_frequency_matrix



mis_RFE_RF_Freq_feature_frequency_matrix
mis_RFE_RF_Freq <- barplot(feature_frequency_matrix_4[1,], xaxt="n", beside=TRUE, cex.names=0.55, las=3, col=c('black'), main='Mutual Info - Sym Uncert - RFE')
axis(2,at=seq(0,50,5))
xlabs <- paste(colnames(feature_frequency_matrix_4))
text(mis_RFE_RF_Freq, par("usr")[3] - 0.025, srt = 90, adj = 1, 
     labels = xlabs, 
     xpd = TRUE, font = 2, cex=0.6)


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


mis_RFE_RF_Freq_new_sensitivities_4 = new_sensitivities
mis_RFE_RF_Freq_new_specificities_4 = new_specificities



par(pty='s')
plot(1-mis_RFE_RF_Freq_new_specificities_4, mis_RFE_RF_Freq_new_sensitivities_4, type='l', xlab='1-Specificity', ylab='Sensitivity', col='black', lwd=2.0, 
     main = 'Mutual Info - Sym Uncert - RFE')
abline(coef=c(0,1),lwd=0.5)

mi_symun_rfe_fs = feature_stability_list

#---Feature Stability plot for all models--#

feature_stability_combined_vars <- data.frame(matrix(0,nrow=50,ncol=16))
feature_stability_combined_vars[,1] = mi_chi_lasso_fs
feature_stability_combined_vars[,2] = mi_oneR_lasso_fs
feature_stability_combined_vars[,3] = mi_infogain_lasso_fs
feature_stability_combined_vars[,4] = mi_symun_lasso_fs
feature_stability_combined_vars[,5] = mi_chi_forward_fs
feature_stability_combined_vars[,6] = mi_oneR_forward_fs
feature_stability_combined_vars[,7] = mi_infogain_forward_fs
feature_stability_combined_vars[,8] = mi_symun_forward_fs
feature_stability_combined_vars[,9] = mi_chi_backward_fs
feature_stability_combined_vars[,10] = mi_oneR_backward_fs
feature_stability_combined_vars[,11] = mi_infogain_backward_fs
feature_stability_combined_vars[,12] = mi_symun_backward_fs
feature_stability_combined_vars[,13] = mi_chi_rfe_fs
feature_stability_combined_vars[,14] = mi_oneR_rfe_fs
feature_stability_combined_vars[,15] = mi_infogain_rfe_fs
feature_stability_combined_vars[,16] = mi_symun_rfe_fs

colnames(feature_stability_combined_vars) = c('MI - Chi-Sq - LASSO','MI - oneR - LASSO','MI - InfoGain - LASSO','MI - SymUncert - LASSO','MI - Chi-Sq - Forward','MI - oneR - Forward','MI - InfoGain - Forward','MI - SymUncert - Forward','MI - Chi-Sq - Backward','MI - oneR - Backward','MI - InfoGain - Backward','MI - SymUncert - Backward','MI - Chi-Sq - RFE','MI - oneR - RFE','MI - InfoGain - RFE','MI - SymUncert - RFE')

ggplot(stack(feature_stability_combined_vars), aes(x = ind, y = values)) +
  geom_boxplot( fill=rep(c('#F8766D','#619CFF','#00BA38','grey'),4)) +
  theme(axis.text.x = element_text(face = 'bold',angle = 90, vjust = 0.5, hjust=0.5, size=12), axis.text.y = element_text(face='bold', size=12), plot.title = element_text(hjust = 0.5), axis.title = element_text(face='bold')) +
  labs(title='Feature Stability Plot for Feature Selection Models on Combined Variables' , x = 'Models') +
  scale_y_continuous(breaks = seq(0, 60, by = 2))