#Loading of important packages

library(dplyr)
library(caTools)
library(randomForest)
library(pROC)
library(MASS)
library(caret)

#Data loading

data <- read.csv("C:\\Users\\mehak\\Downloads\\csv_result-bone-marrow.csv") 

#Data pre-processing

data$Disease <- recode(data$Disease, 'ALL' = '4', 'lymphoma' = '0', 'AML' = '1', 'chronic' = '2', 'nonmalignant' = '3')
data_numeric <- data.frame(lapply(data, function(x) as.numeric(x)))  #converting the values to numeric
data_clean <- na.omit(data_numeric)  #removing the NA values
correlation <- cor(data_clean)
data_clean <- subset(data_clean, select = -id)

#Removal of highly correlated features: survival_time and relapse

data_clean <- subset(data_clean, select = -survival_time)
data_clean <- subset(data_clean, select = -Relapse)
View(data_clean)
#Feature selection using rfe: recursive feature elimination

data_clean$survival_status <- as.factor(data_clean$survival_status)

set.seed(2)

control <- rfeControl(functions = rfFuncs,
                          method = "cv",
                          number = 10)

results <- rfe(x = data_clean[-35], 
                   y = data_clean$survival_status, 
                   sizes = c(1:34), 
                   rfeControl = control)
print(results)
print("10 variables got selected for best accuracy out of which top 5 features are: Riskgroup, Txpostrelapse, extcGvHD, RecipientRh and Disease")
plot(results)
selected_features <- predictors(results)
cat("The 10 selected features are:", selected_features)

#Supervised models on 2 types of data: whole data set and reduced data set with top 5 features selected after feature selection

#Supervised model 1: Random Forest using 5 fold cross validation

#Whole data:

set.seed(2)
num_splits <- 5
par(mfrow=c(3, 2))
train_data_list <- list()
test_data_list <- list()
accuracy_list_train <- list()
accuracy_list_test <- list()
auc_list <- list()
tpr_list <- list()
fpr_list <- list()
color_list <- c("red", "blue", "green", "pink", "purple")

for(i in 1:num_splits) {
  split <- sample.split(data_clean$survival_status, SplitRatio = 0.7)
  train_data <- subset(data_clean, split == TRUE)
  test_data <- subset(data_clean, split == FALSE)
  
  train_data_list[[i]] <- train_data
  test_data_list[[i]] <- test_data
  
  # Random Forest model
  model_RF <- randomForest(x = train_data[-35], 
                           y = train_data$survival_status, 
                           ntree = 500)
  
  # Predictions on train data
  RF_pred_train <- predict(model_RF, newdata = train_data)
  
  # Predictions on test data
  test_data <- na.omit(test_data)
  RF_pred_test <- predict(model_RF, newdata = test_data)
  
  #Calculate accuracy
  accuracy_train <- sum(diag(table(train_data[, 35], RF_pred_train)))/length(train_data$survival_status)
  accuracy_test <- sum(diag(table(test_data[, 35], RF_pred_test)))/length(test_data$survival_status)
  
  accuracy_list_train[[i]] <- accuracy_train
  accuracy_list_test[[i]] <- accuracy_test
  
  #ROC curve
  RF_probs <- predict(model_RF, newdata = test_data, type = "prob")
  positive_class_probs <- RF_probs[, 1]
  roc_obj <- roc(test_data$survival_status, positive_class_probs)
  auc_list[[i]] <- auc(roc_obj)
  tpr_list[[i]] <- roc_obj$sensitivities
  fpr_list[[i]] <- 1 - roc_obj$specificities
  plot(fpr_list[[i]], tpr_list[[i]], type = "l", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "False Positive Rate", ylab = "True Positive Rate",
       main = paste("ROC Curve for test split", i), col = color_list[i], lwd = 2)
  abline(a = 0, b = 1, lty = 2, col = "gray")
}
print(accuracy_list_train)
print(accuracy_list_test)
print(auc_list)
mean(unlist(accuracy_list_train))
mean(unlist(accuracy_list_test))
mean(unlist(auc_list))
sd(unlist(accuracy_list_test))
sd(unlist(auc_list))
par(mfrow=c(3, 2))

#Reduced data:

data_new <- data.frame(data_clean$Riskgroup, data_clean$Txpostrelapse, data_clean$extcGvHD, data_clean$RecipientRh, data_clean$Disease, data_clean$survival_status)

set.seed(2)
par(mfrow=c(3, 2))
train_data_list_new <- list()
test_data_list_new <- list()
accuracy_list_new_train <- list()
accuracy_list_new_test <- list()
auc_list_new <- list()
tpr_list_new <- list()
fpr_list_new <- list()

for(i in 1:num_splits) {
  split_new <- sample.split(data_new$data_clean.survival_status, SplitRatio = 0.7)
  train_data_new <- subset(data_new, split_new == TRUE)
  test_data_new <- subset(data_new, split_new == FALSE)
  
  train_data_list_new[[i]] <- train_data_new
  test_data_list_new[[i]] <- test_data_new
  
  # Random Forest model
  model_RF_new <- randomForest(x = train_data_new[-6], 
                           y = train_data_new$data_clean.survival_status, 
                           ntree = 500)
  
  # Predictions on train data
  RF_pred_new_train <- predict(model_RF_new, newdata = train_data_new)
  
  # Predictions on test data
  RF_pred_new_test <- predict(model_RF_new, newdata = test_data_new)
  
  #Calculate accuracy
  accuracy_new_train <- sum(diag(table(train_data_new[, 6], RF_pred_new_train)))/length(train_data_new$data_clean.survival_status)
  accuracy_new_test <- sum(diag(table(test_data_new[, 6], RF_pred_new_test)))/length(test_data_new$data_clean.survival_status)
  accuracy_list_new_train[[i]] <- accuracy_new_train
  accuracy_list_new_test[[i]] <- accuracy_new_test
  
  #ROC curve
  RF_probs_new <- predict(model_RF_new, newdata = test_data_new, type = "prob")
  positive_class_probs_new <- RF_probs_new[, 1]
  roc_obj_new <- roc(test_data_new$data_clean.survival_status, positive_class_probs_new)
  auc_list_new[[i]] <- auc(roc_obj_new)
  tpr_list_new[[i]] <- roc_obj_new$sensitivities
  fpr_list_new[[i]] <- 1 - roc_obj_new$specificities
  plot(fpr_list_new[[i]], tpr_list_new[[i]], type = "l", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "False Positive Rate", ylab = "True Positive Rate",
       main = paste("ROC Curve for test split", i), col = color_list[i], lwd = 2)
  abline(a = 0, b = 1, lty = 2, col = "gray")
}
print(accuracy_list_new_train)
print(accuracy_list_new_test)
print(auc_list_new)

mean(unlist(accuracy_list_new_train))
mean(unlist(accuracy_list_new_test))
mean(unlist(auc_list_new))
sd(unlist(accuracy_list_new_test))
sd(unlist(auc_list_new))
par(mfrow=c(3, 2))

#Supervised model 2: LDA using 5 fold cross validation

#Whole data:

set.seed(2)
par(mfrow=c(3, 2))
accuracy_list_lda_train <- list()
accuracy_list_lda_test <- list()
auc_list_lda <- list()
tpr_list_lda <- list()
fpr_list_lda <- list()

for(i in 1:num_splits) {
  split <- sample.split(data_clean$survival_status, SplitRatio = 0.7)
  train_data <- subset(data_clean, split == TRUE)
  test_data <- subset(data_clean, split == FALSE)
  #LDA model
  lda_model <- lda(survival_status ~ ., data = train_data)
  
  #Predictions on train data
  LDA_pred_train <- predict(lda_model, newdata = train_data) 

  #Predictions on test data
  LDA_pred_test <- predict(lda_model, newdata = test_data)
  
  #Calculate accuracy
  accuracy_lda_train <- sum(diag(table(train_data$survival_status, LDA_pred_train$class))) / nrow(train_data)
  accuracy_lda_test <- sum(diag(table(test_data$survival_status, LDA_pred_test$class))) / nrow(test_data)
  accuracy_list_lda_train[[i]] <- accuracy_lda_train
  accuracy_list_lda_test[[i]] <- accuracy_lda_test
  
  #ROC curve
  roc_obj_lda <- roc(test_data$survival_status, LDA_pred_test$posterior[,1])
  auc_list_lda[[i]] <- auc(roc_obj_lda)
  tpr_list_lda[[i]] <- roc_obj_lda$sensitivities
  fpr_list_lda[[i]] <- 1 - roc_obj_lda$specificities
  plot(fpr_list_lda[[i]], tpr_list_lda[[i]], type = "l", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "False Positive Rate", ylab = "True Positive Rate",
       main = paste("ROC Curve for test split", i), col = color_list[i], lwd = 2)
  abline(a = 0, b = 1, lty = 2, col = "gray")
}
accuracy_list_lda_train
accuracy_list_lda_test
auc_list_lda

mean(unlist(accuracy_list_lda_train))
mean(unlist(accuracy_list_lda_test))
mean(unlist(auc_list_lda))
sd(unlist(accuracy_list_lda_test))
sd(unlist(auc_list_lda))
par(mfrow=c(3, 2))

#Reduced data set

set.seed(2)
par(mfrow=c(3, 2))
accuracy_list_lda_train_new <- list()
accuracy_list_lda_test_new <- list()
auc_list_lda_new <- list()
tpr_list_lda_new <- list()
fpr_list_lda_new <- list()

for(i in 1:num_splits) {
  split_new <- sample.split(data_new$data_clean.survival_status, SplitRatio = 0.7)
  train_data_new <- subset(data_new, split_new == TRUE)
  test_data_new <- subset(data_new, split_new == FALSE)
  
  #LDA model
  lda_model_new <- lda(data_clean.survival_status ~ ., data = train_data_new)
  
  #Predictions on train data
  LDA_pred_train_new <- predict(lda_model_new, newdata = train_data_new) 
  
  #Predictions on test data
  LDA_pred_test_new <- predict(lda_model_new, newdata = test_data_new)
  
  #Calculate accuracy
  accuracy_lda_train_new <- sum(diag(table(train_data_new$data_clean.survival_status, LDA_pred_train_new$class))) / nrow(train_data_new)
  accuracy_lda_test_new <- sum(diag(table(test_data_new$data_clean.survival_status, LDA_pred_test_new$class))) / nrow(test_data_new)
  accuracy_list_lda_train_new[[i]] <- accuracy_lda_train_new
  accuracy_list_lda_test_new[[i]] <- accuracy_lda_test_new
  
  #ROC curve
  roc_obj_lda_new <- roc(test_data_new$data_clean.survival_status, LDA_pred_test_new$posterior[,1])
  auc_list_lda_new[[i]] <- auc(roc_obj_lda_new)
  tpr_list_lda_new[[i]] <- roc_obj_lda_new$sensitivities
  fpr_list_lda_new[[i]] <- 1 - roc_obj_lda_new$specificities
  plot(fpr_list_lda_new[[i]], tpr_list_lda_new[[i]], type = "l", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "False Positive Rate", ylab = "True Positive Rate",
       main = paste("ROC Curve for test split", i), col = color_list[i], lwd = 2)
  abline(a = 0, b = 1, lty = 2, col = "gray")
}
accuracy_list_lda_train_new
accuracy_list_lda_test_new
auc_list_lda_new

mean(unlist(accuracy_list_lda_train_new))
mean(unlist(accuracy_list_lda_test_new))
mean(unlist(auc_list_lda_new))
sd(unlist(accuracy_list_lda_test_new))
sd(unlist(auc_list_lda_new))
par(mfrow=c(3, 2))

#Supervised model 3: Logistic regression using 5 fold cross validation

#Whole data:

set.seed(2)
par(mfrow=c(3, 2))
accuracy_list_lr_train <- list()
accuracy_list_lr_test <- list()
auc_list_lr <- list()
tpr_list_lr <- list()
fpr_list_lr <- list()

for(i in 1:num_splits) {
  split <- sample.split(data_clean$survival_status, SplitRatio = 0.7)
  train_data <- subset(data_clean, split == TRUE)
  test_data <- subset(data_clean, split == FALSE)
  
  #lr model
  lr_model <- glm(survival_status ~. , family="binomial", data = train_data)
  
  #Predictions on train data
  LR_pred_train <- pred <- predict(lr_model, newdata = train_data, type = "response")
  y_pred_train <- ifelse(LR_pred_train > 0.5, "pos", "neg")
  
  #Predictions on test data
  LR_pred_test <- predict(lr_model, newdata = test_data, type = "response")
  y_pred_test <- ifelse(LR_pred_test > 0.5, "pos", "neg")
  
  #Calculate accuracy
  accuracy_lr_train <- sum(diag(table(y_pred_train, train_data$survival_status)))/length(train_data$survival_status)
  accuracy_lr_test <- sum(diag(table(y_pred_test, test_data$survival_status)))/length(test_data$survival_status)
  accuracy_list_lr_train[[i]] <- accuracy_lr_train
  accuracy_list_lr_test[[i]] <- accuracy_lr_test
  
  #ROC curve
  roc_obj_lr <- roc(test_data$survival_status, LR_pred_test)
  auc_list_lr[[i]] <- auc(roc_obj_lr)
  tpr_list_lr[[i]] <- roc_obj_lr$sensitivities
  fpr_list_lr[[i]] <- 1 - roc_obj_lr$specificities
  plot(fpr_list_lr[[i]], tpr_list_lr[[i]], type = "l", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "False Positive Rate", ylab = "True Positive Rate",
       main = paste("ROC Curve for test split", i), col = color_list[i], lwd = 2)
  abline(a = 0, b = 1, lty = 2, col = "gray")
}
accuracy_list_lr_train
accuracy_list_lr_test
auc_list_lr

mean(unlist(accuracy_list_lr_train))
mean(unlist(accuracy_list_lr_test))
mean(unlist(auc_list_lr))
sd(unlist(accuracy_list_lr_test))
sd(unlist(auc_list_lr))
par(mfrow=c(3, 2))

#Reduced data

set.seed(2)
par(mfrow=c(3, 2))
accuracy_list_lr_train_new <- list()
accuracy_list_lr_test_new <- list()
auc_list_lr_new <- list()
tpr_list_lr_new <- list()
fpr_list_lr_new <- list()

for(i in 1:num_splits) {
  split_new <- sample.split(data_new$data_clean.survival_status, SplitRatio = 0.7)
  train_data_new <- subset(data_new, split_new == TRUE)
  test_data_new <- subset(data_new, split_new == FALSE)
  
  #lr model
  lr_model_new <- glm(data_clean.survival_status ~. , family="binomial", data = train_data_new)
  
  #Predictions on train data
  LR_pred_train_new <- predict(lr_model_new, newdata = train_data_new, type = "response")
  y_pred_train_new <- ifelse(LR_pred_train_new > 0.5, "pos", "neg")
  
  #Predictions on test data
  LR_pred_test_new <- predict(lr_model_new, newdata = test_data_new, type = "response")
  y_pred_test_new <- ifelse(LR_pred_test_new > 0.5, "pos", "neg")
  
  #Calculate accuracy
  accuracy_lr_train_new <- sum(diag(table(y_pred_train_new, train_data_new$data_clean.survival_status)))/length(train_data_new$data_clean.survival_status)
  accuracy_lr_test_new <- sum(diag(table(y_pred_test_new, test_data_new$data_clean.survival_status)))/length(test_data_new$data_clean.survival_status)
  accuracy_list_lr_train_new[[i]] <- accuracy_lr_train_new
  accuracy_list_lr_test_new[[i]] <- accuracy_lr_test_new
  
  #ROC curve
  roc_obj_lr_new <- roc(test_data_new$data_clean.survival_status, LR_pred_test_new)
  auc_list_lr_new[[i]] <- auc(roc_obj_lr_new)
  tpr_list_lr_new[[i]] <- roc_obj_lr_new$sensitivities
  fpr_list_lr_new[[i]] <- 1 - roc_obj_lr_new$specificities
  plot(fpr_list_lr_new[[i]], tpr_list_lr_new[[i]], type = "l", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "False Positive Rate", ylab = "True Positive Rate",
       main = paste("ROC Curve for test split", i), col = color_list[i], lwd = 2)
  abline(a = 0, b = 1, lty = 2, col = "gray")
}
accuracy_list_lr_train_new
accuracy_list_lr_test_new
auc_list_lr_new

mean(unlist(accuracy_list_lr_train_new))
mean(unlist(accuracy_list_lr_test_new))
mean(unlist(auc_list_lr_new))
sd(unlist(accuracy_list_lr_test_new))
sd(unlist(auc_list_lr_new))
par(mfrow=c(3, 2))