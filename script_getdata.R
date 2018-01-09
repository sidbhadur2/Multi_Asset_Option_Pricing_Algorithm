setwd("/Users/siddharthbhaduri/Desktop/Work/Spring-2017/IE-525-Numerical_Methods/Multi_Asset_Option_Project/")

library(fImport)
ETFs= read.csv("/Users/siddharthbhaduri/Desktop/Work/Spring-2017/IE-525-Numerical_Methods/Multi_Asset_Option_Project/Test_ETFs.csv")

date_from <- "11/11/2016"
date_to <- "01/01/2017"

for(i in 1:length(ETFs[,1])){
  Ticker = as.character(ETFs[i,1])
  ETF.data <- yahooSeries(Ticker, from = date_from, to= date_to)
  if(i == 1) {
    ETFs.data <- ETF.data[,4]
  } else {
    ETFs.data <- cbind(ETFs.data, ETF.data[,4])
  }
}

#Initial_volatility.data <- as.numeric(ETFs.data[1,])

Initial_price.data <- as.numeric(ETFs.data[1,])
input_table <- data.frame(ETFs$Symbol[1:length(ETFs[,1])], Initial_price.data, 0.02, 0.0213, 1)
colnames(input_table) <-c("Ticker", "S_price", "Volatility", "r" , "T")

sprintf("Complete")
View(ETFs)
View(ETFs.data)
View(input_table)