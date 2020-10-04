library(dplyr)
library(lubridate)
library(xts)
library(imputeTS)
library(purrr)
fit<-read.csv("/Users/ceciliaxia/Downloads/Fitbit code/fitbit.csv")
fit_split<-split.data.frame(fit,fit$ID)

evaluate_peaks_coef<-function(subject,lag,threshold,influence){
ThresholdingAlgo <- function(y,lag,threshold,influence) {
  signals <- rep(0,length(y))
  filteredY <- y[0:lag]
  avgFilter <- NULL
  stdFilter <- NULL
  avgFilter[lag] <- mean(y[0:lag], na.rm=TRUE)
  stdFilter[lag] <- sd(y[0:lag], na.rm=TRUE)
  for (i in (lag+1):length(y)){
    if (abs(y[i]-avgFilter[i-1]) > threshold*stdFilter[i-1]) {
      if (y[i] > avgFilter[i-1]) {
        signals[i] <- 1;
      } else {
        signals[i] <- -1;
      }
      filteredY[i] <- influence*y[i]+(1-influence)*filteredY[i-1]
    } else {
      signals[i] <- 0
      filteredY[i] <- y[i]
    }
    avgFilter[i] <- mean(filteredY[(i-lag):i], na.rm=TRUE)
    stdFilter[i] <- sd(filteredY[(i-lag):i], na.rm=TRUE)
  }
  return(list("signals"=signals,"avgFilter"=avgFilter,"stdFilter"=stdFilter))
}

  x<-fit_split[[subject]]
  
  month<-c("JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")
  monnum<-c("01","02","03","04","05","06","07","08","09","10","11","12")
  date<-x$date
  
  for(i in 1:12){
    date<-gsub(month[i],monnum[i],date)
  }
  
  x$date<-as.Date(date,format="%d%m%Y")
  days<-split.data.frame(x,x$date)
  date_day<-unique(x$date)
  dataset2<-cbind.data.frame(date_day,rep(0,length(date_day)),rep(0,length(date_day)))
  colnames(dataset2)<-c("day","is_peak","is_not_peak")
  
  
  
  
  for(i in 1:length(date_day)){
    
  z <- seq.POSIXt(as.POSIXct(date_day[i])+7*60*60, as.POSIXct(date_day[i]+1)+7*60*60, by = "1 min")
  tzone(z)<-"America/Los_angeles"
  minute<-hour(z)*60+minute(z)
  newdata<-data.frame(minute=minute,DATE=z)
  days[[i]]<-inner_join(newdata,days[[i]],by="minute")
  days[[i]]<-days[[i]][,c(2,4,5)]
  day_ts_min<-zoo(days[[i]][,-1],order.by=days[[i]]$DATE)
  df<- merge(day_ts_min,zoo(,seq(start(day_ts_min),end(day_ts_min),by='1 min')), all=TRUE)
  df1<-data.frame(DATE=seq(start(day_ts_min),end(day_ts_min),by='1 min'),steps=as.numeric(df$steps))
  df1$steps<-na_replace(df1$steps,fill=0)
  y<-df1$steps
  
  tryCatch({
  result <- ThresholdingAlgo(y,lag,threshold,influence)
  dataset1<-cbind.data.frame(df1$DATE,y,result$signals,1:length(y))
  colnames(dataset1)<-c("Date","steps","signal","Absolute_Date")
  regression<-lm(steps~Absolute_Date+signal,data=dataset1)
  
  dataset2$is_peak[i]<-regression$coefficients[1]+regression$coefficients[3]
  dataset2$is_not_peak[i]<-regression$coefficients[1]}, error=function(e){})
}
  par(mfrow=c(2,1))
  hist(dataset2$is_peak,main="Histogram of peaks coefficient")
  hist(dataset2$is_not_peak,main="Histogram of non_peaks coefficient")

  }
evaluate_peaks_coef(6,10,3.5,0.1)
  
  


