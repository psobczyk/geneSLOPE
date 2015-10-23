# auxilliary functions

#fast p-value computation for simple marginal lm fit test
pValComp <- function(x,y,n,suma){
  a <- lm.fit(cbind(x,1),y)
  b <- sum(a$residuals^2)
  1-pf((suma-b)/b*n,1, n)
}


# Function to replace missing values with mean for that col
replace_na_with_mean <- function(x) {
  x_bar <- mean(x, na.rm = TRUE)
  ifelse(is.na(x), x_bar, x)
}

# recoding snps as was described in
# http://pngu.mgh.harvard.edu/~purcell/plink/dataman.shtml#recode
recodeAD <- function(x){
  sapply(seq(1, length(x), 2), function(i){
    if(any(is.na(x[c(i,i+1)])))
      return(c(NA,NA))
    if(any(x[c(i,i+1)]==0))
      return(c(NA,NA))
    if(x[i]==1 & x[i+1]==1){
      c(0,0)
    } else if(x[i]==2 & x[i+1]==2){
      c(2,0)
    } else{
      c(1,1)
    }
  })
}
