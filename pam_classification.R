#load library for PAM
library(pamr);

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
n_repeat = 100;
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
}

#print results
print("Mean and standard deviation of the test error (%)");
print(mean(test_error_array)*100);
print(sd(test_error_array)*100);