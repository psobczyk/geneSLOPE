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
