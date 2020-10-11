library(dplyr)
library(lubridate)
library(xts)
library(imputeTS)
library(purrr)
library(dtwclust)
fit<-read.csv("/Users/ceciliaxia/Downloads/Fitbit code/fitbit.csv")
fit_split<-split.data.frame(fit,fit$ID)

cluster_plot<-function(subject,day,lag,threshold,influence){
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
  
  i<-day
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
  y<-df1$steps
  result <- ThresholdingAlgo(y,lag,threshold,influence)
  original_signals<-result$signals
  back_up_signals<-result$signals
  new_signals<-rep(0,length(result$signals))
  for(i in (lag+1):(length(result$signals)-2)){
    
    if(original_signals[i]==0){
    front<-original_signals[(i+1):(i+2)]
    back<-original_signals[(i-2):(i-1)]
    sum<-sum(back==rep(1,2))+sum(front==rep(1,2))
    if(sum>1|sum==1){back_up_signals[i]<-1}
    }
  }
  
  
  rle<-rle(back_up_signals)
  order<-length(rle$lengths)
  store_point<-rep(0,length(df1$DATE))
  store_list<-list()
  ts_list<-list()
  sum=0
  for(i in 1:order){
    sum=sum+rle$lengths[i]
    if(i%%2==0){
      a=sum+1-rle$lengths[i]
      b=sum
      new_data<-df1[a:b,]
      store_list[(i/2)]<-list(new_data)
      ts_list[(i/2)]<-list(y[a:b])
    }
  }
  cluster_result<-tsclust(series=ts_list,k=3)
  cluster<-cluster_result@cluster
  cluster_result<-tsclust(series=ts_list,k=3)
  cluster<-cluster_result@cluster
  plot(cluster_result)
  plot(df1$DATE,y,type="l",ylab="",xlab="")
  for(i in c(which(cluster==1))){
    t_data<-store_list[[i]]
    lines(t_data$DATE,t_data$steps,col=2,lwd=2)
  }
  for(i in c(which(cluster==2))){
    t_data<-store_list[[i]]
    lines(t_data$DATE,t_data$steps,col=3,lwd=2)
  }
  for(i in c(which(cluster==3))){
    t_data<-store_list[[i]]
    lines(t_data$DATE,t_data$steps,col=4,lwd=2)
  }
  legend("topright", legend=c("type 1", "type 2","type 3"),
         col=c(2,3,4), lty=1, cex=0.8,lwd=2,
         title="cluster types", text.font=4, bg='lightblue')
}
 
 
 
