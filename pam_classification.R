#load library for PAM
library(pamr);

#load library for error bar charts
library(Hmisc);

#set random seed
set.seed(174457795);

#number of bootstrap, aka number of times to repeat the experiment
n_repeat = 500;
#max_threshold of cross validation
max_threshold = 10;
#number of thresholds to look at in cross validation
n_threshold = 30;

#NOTE FOR pamr: need to convert to data.frame then to matrix

#read in design matrix
X = read.table("gene146.csv");
#convert to data frame, convert to matrix
X = as.matrix(as.data.frame(X));
n = ncol(X);

#read in response vector (class labels);
Y = read.table("y146gene")
#convert numericals to factors
Y[Y==1] = "Label 1"; #factor 1
Y[Y==2] = "Label 2"; #factor 2
Y[Y==3] = "Label 3"; #factor 3
#convert to data frame
Y = as.data.frame(Y);
#convert to matrix
Y = as.matrix(Y);

#make a copy of X and Y
X_bootstrap = X;
Y_bootstrap = Y;

#the thresholds to look at in cross validation
thresholdCV_array = seq(from=0,to=max_threshold,length.out=n_threshold);

test_error_array = rep(0,n_repeat); #array of test error for every repeat of the experiment
threshold_array = rep(0,n_repeat); #array of optimal threshold for every repeat of the experiment
n_genes_array = rep(0,n_repeat); #array of number of genes which survive the optimal threshold for every repeat of the experiment
missclassification_array = matrix(0,nrow=n_repeat,ncol=n_threshold); #matrix of missclassification for every experiment (rows) for every threshold (column)

#for n_repeat times
for (i in 1:n_repeat){

  #if this is not the 1st run, bootstrap the data
  if (i >= 2){
    index = sample(1:n,n,TRUE); 
    X_bootstrap = X[,index];
    Y_bootstrap = Y[index];
  }
  
  #TRAIN THE CLASSIFIER
  train = pamr.train(list(x = X_bootstrap, y = Y_bootstrap),threshold = thresholdCV_array);
  #GET CROSS VALIDATION ERROR
  results = pamr.cv(train, list(x = X_bootstrap, y = Y_bootstrap), nfold = 10);
  
  #SAVE THE MINIMUM VALIDATION ERROR
  test_error_array[i] = min(results$error);
  #WORK OUT THE INDEX OF THE OPTIMAL THRESHOLD: the largest threshold with the minimum error 
  optimal_index = max(which(results$error==min(results$error)));
  #SAVE THE OPTIMAL THRESHOLD
  threshold_array[i] = results$threshold[optimal_index];
  #SAVE THE NUMBER OF GENES SURVIVE
  n_genes_array[i] = results$size[optimal_index]
  #SAVE THE MISSCLASSIFICATION ERROR
  missclassification_array[i,] = results$error;
  
  #if this is the first run, get the CV plot
  if (i == 1){
    pamr.plotcv(results);
  }
  #save the gene names
  #pamr.listgenes(train, list(x = X_bootstrap, y = Y_bootstrap), threshold=0.1, genenames=TRUE);
}

#print results of the test error
print("Mean and standard deviation of the test error (%)");
print( paste( mean(test_error_array)*100, "+/-", sd(test_error_array)*100)  );
#state the estimated standard deviation of the standard deviation in order to judge the number of significant figures to be used
print("Estimated standard deviation of standard deviation (%)");
print(sd(test_error_array)*100/sqrt(2*(n_repeat-1)));

print("");

#print results if the optimal threshold
print("Mean and standard deviation of the optimal threshold (AU)");
print( paste( mean(threshold_array), "+/-", sd(threshold_array) ) );
#state the estimated standard deviation of the standard deviation in order to judge the number of significant figures to be used
print("Estimated standard deviation of standard deviation (AU)");
print(sd(threshold_array)/sqrt(2*(n_repeat-1)));

print("");

#print results of the number of surviving genes
print("Mean and standard deviation of the number of surviving genes (count)");
print( paste( mean(n_genes_array), "+/-", sd(n_genes_array)) );
#state the estimated standard deviation of the standard deviation in order to judge the number of significant figures to be used
print("Estimated standard deviation of standard deviation (count)");
print(sd(n_genes_array)/sqrt(2*(n_repeat-1)));

#plot the error for each threshold
missclassifcation_mean = colMeans(missclassification_array)*100; #calculate the mean
missclassification_sd = apply(missclassification_array,2,sd)*100; #calculate the std (error bar)
errbar(thresholdCV_array, missclassifcation_mean, missclassifcation_mean+missclassification_sd, missclassifcation_mean-missclassification_sd,xlab="Threshold",ylab="Validation Error (%)");

#save the variables
save(list = ls(all.names = TRUE), file = "pam_classification_results_500.RData", envir = .GlobalEnv);