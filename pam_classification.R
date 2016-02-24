#load library for PAM
library(pamr);
#load library for cum. sd.
library(TTR);

#NOTE FOR pamr: need to convert to data.frame then to matrix

#read in design matrix
X = read.table("gene146.csv");
#convert to data frame, convert to matrix
X = as.matrix(as.data.frame(X));
n = ncol(X);

#read in response vector (class labels);
Y = read.table("y146gene")
#convert numericals to factors
Y[Y==1] = "CCS1"; #factor 1
Y[Y==2] = "CCS2"; #factor 2
Y[Y==3] = "CCS3"; #factor 3
#convert to data frame
Y = as.data.frame(Y);
#convert to matrix
Y = as.matrix(Y);

#make a copy of X and Y
X_bootstrap = X;
Y_bootstrap = Y;

#for n_repeat times
n_repeat = 500;
test_error_array = rep(0,n_repeat); #array of test error for every repeat of the experiment
for (i in 1:n_repeat){

  #if this is not the 1st run, bootstrap the data
  if (i >= 2){
    index = sample(1:n,n,TRUE); 
    X_bootstrap = X[,index];
    Y_bootstrap = Y[index];
  }
  
  #TRAIN THE CLASSIFIER
  train = pamr.train(list(x = X_bootstrap, y = Y_bootstrap));
  #GET CROSS VALIDATION ERROR
  results = pamr.cv(train, list(x = X_bootstrap, y = Y_bootstrap), nfold = 10);
  #SAVE THE MINIMUM VALIDATION ERROR
  test_error_array[i] = min(results$error);
  
  #save the gene names
  #pamr.listgenes(train, list(x = X_bootstrap, y = Y_bootstrap), threshold=0.1, genenames=TRUE);
}

#print results
print("Mean and standard deviation of the test error (%)");
print(mean(test_error_array)*100);
print(sd(test_error_array)*100);

#plot running standard deviation to judge how many significant figures to use
#runSD(test_error_array, n=1, cumulative=TRUE)[(n_repeat/2):n_repeat]
plot(runSD(test_error_array, n=1, cumulative=TRUE)[(n_repeat/2):n_repeat],type="l")
hist(runSD(test_error_array, n=1, cumulative=TRUE)[(n_repeat/2):n_repeat])
print("Standard deviation of 2nd half of running standard deviations (%)");
print(sd(runSD(test_error_array, n=1, cumulative=TRUE)[(n_repeat/2):n_repeat])*100);
print("Estimate standard deviation of standard deviation (%)");
print(sd(test_error_array)*100/sqrt(2*(n_repeat-1)));

#save the variables
save(list = ls(all.names = TRUE), file = "pam_classification_results_500.RData", envir = .GlobalEnv);