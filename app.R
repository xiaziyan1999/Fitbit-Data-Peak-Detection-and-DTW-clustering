#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
options(shiny.maxRequestSize = 3000*1024^2)
library(shiny)
library(shinythemes)
library(dplyr)
library(lubridate)
library(xts)
library(imputeTS)
library(purrr)
library(dtwclust)

# Define UI for application that draws a histogram
ui <- fluidPage(
    theme = shinytheme("superhero"),
    
    # Application title
    titlePanel("Characterize Intra-subject Variability"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            fileInput("file", "Choose FITBIT File",
                      multiple = TRUE,
                      accept = c("text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")),
            selectInput("y",
                        "Subject:",
                        choices = c(1:104),
                        selected = 5 ),
            sliderInput("d",
                        label = "Day:",
                        min = 1,
                        max = 1000,
                        step=1,
                        value =1 ),
            sliderInput("lag",
                        label = "Lag:",
                        min = 1,
                        max = 100,
                        step=1,
                        value =10 ),
            sliderInput("thresh",
                        label = "Threshold:",
                        min = 1,
                        max = 5,
                        step=0.5,
                        value =3.5),
            sliderInput("inf",
                        label = "Influence:",
                        min = 0,
                        max = 1,
                        step=0.01,
                        value =0.1 ),
            radioButtons("co", "Choose a plot",
                         choices = c("Peak Detection results"= "arima",
                                     "Cluster results" = "struct",
                                     "Cluster on Time Series"="comp")
            ),
        ),
        # Show a plot of the generated distribution
        mainPanel(
            plotOutput("plot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output,session) {
    
    xdat<-reactive({
        req(input$file)
        fit<-read.csv(input$file$datapath)
        fit_split<-split.data.frame(fit,fit$ID)
        subject<-paste(input$y)
        dat_ori<-fit_split[subject]
        x<-dat_ori[[1]]
        month<-c("JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")
        monnum<-c("01","02","03","04","05","06","07","08","09","10","11","12")
        date<-x$date
        for(i in 1:12){
            date<-gsub(month[i],monnum[i],date)
        }
        x$date<-as.Date(date,format="%d%m%Y")
        x
    })
    
    output$plot <- renderPlot({
        subject<-paste(input$y)
        days<-split.data.frame(xdat(),xdat()$date)
        date_day<-unique(xdat()$date)
        i<-input$d
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
        Lag<-input$lag
        threshold<-input$thresh
        influence<-input$inf
        result <- ThresholdingAlgo(y,Lag,threshold,influence)
        original_signals<-result$signals
        back_up_signals<-result$signals
        new_signals<-rep(0,length(result$signals))
        for(i in (Lag+1):(length(result$signals)-2)){
            
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
        if(input$co=="arima"){
        # Plot result
        par(mfrow = c(2,1),oma = c(2,2,0,0) + 0.1,mar = c(0,0,2,1) + 0.2)
        plot(df1$DATE,y,type="l",ylab="",xlab="",main = paste("Segmentation Results for",paste("subject",subject),"on",date_day[i])) 
        lines(df1$DATE,result$avgFilter,type="l",col="cyan",lwd=2)
        lines(df1$DATE,result$avgFilter+threshold*result$stdFilter,type="l",col="green",lwd=2)
        lines(df1$DATE,result$avgFilter-threshold*result$stdFilter,type="l",col="green",lwd=2)
        plot(df1$DATE,result$signals,type="S",col="red",ylab="",xlab="",ylim=c(-1.5,1.5),lwd=2)
        }else if(input$co=="struct"){
            plot(df1$DATE,y,type="l",ylab="",xlab="",,main = paste("Clustering Results for",paste("subject",subject),"on",date_day[i]))
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
        }else{plot(cluster_result)
            
        }
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
