
options(shiny.maxRequestSize = 3000*1024^2)
library(shiny)
library(shinythemes)
library(dplyr)
library(lubridate)
library(forecast)
library(xts)
library(imputeTS)
library(writexl)
library(purrr)

# Define UI for data upload app ----
ui <- fluidPage(
    theme = shinytheme("superhero"),
    # Application title
    titlePanel("Fitbit Data Forecasting"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            fileInput("file1", "Choose FITBIT File",
                      multiple = TRUE,
                      accept = c("text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")),
            fileInput("file2", "Choose Temperature File",
                      multiple = TRUE,
                      accept = c("text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")),
            # Horizontal line ----
            tags$hr(),
            
            # Input: Checkbox if file has header ----
            checkboxInput("scale", "Standardize", TRUE),
            # Input: Select time scale ----
            radioButtons("sep", "Time Scale",
                         choices = c("Weekly"="week",
                                     "Daily"="day",
                                     "15 minutes"="15 minutes")
            ),
            # Input: Select subject ----
            selectInput("y",
                        "Subject:",
                        choices = c(1:104),
                        selected = 5 ),
            sliderInput("rate",
                        label = "Trainging Set spiliting Rate:",
                        min = 0.5,
                        max = 1,
                        step=0.01,
                        value =0.9 ),
            radioButtons("summ", "Summarizing Figures and Evaluation Figures",
                         choices = c("ARIMA results"= "arima",
                                     "StrutTS results" = "struct",
                                     "Comparison of Two Models"="comp",
                                     "ACF and PACF plots of ARIMA model"="arac",
                                     "ACF and PACF plots of StructTS model"="strac",
                                     "Barplot of MSE of Nested cross validation"="mse")
            ),
            
            
            
            
            
        ),
        
        # Show a plot of the generated data
        mainPanel(
            plotOutput("plot")
        )
    )
)
# Define server logic to read selected file ----
server <- function(input, output,session) {
    
    xdat<-reactive({
        req(input$file1)
        req(input$file2)
        fit<-read.csv(input$file1$datapath)
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
    
    Temp<-reactive({
        Temp<-read.csv(input$file2$datapath)
        Temp$X<-as.Date(Temp$X)
        colnames(Temp)[1]<-"DATE"
        Temp })
    
    
    output$plot <- renderPlot({
        
        
        
        Sys.setlocale("LC_TIME", "English")
        
        
        
        
        
        if(input$sep=="week"){
            day_count<-xdat() %>%group_by(date) %>%summarize(steps=sum(steps))
            day_count<-day_count[order(day_count$date),]
            day_ts<-zoo(day_count[,-1],order.by=as.Date(day_count$date, "%Y/%m/%d"))
            df<- merge(day_ts,zoo(,seq(start(day_ts),end(day_ts),by="day")), all=TRUE)
            index<-seq(start(day_ts),end(day_ts),by="day")
            df1<-data.frame(DATE=seq(start(day_ts),end(day_ts),by="day"),steps=df$steps)
            
            df2<-inner_join(df1,Temp(),by="DATE")
            
            df2$steps<-na_kalman(df2$steps)
            
            nth<-wday(start(day_ts))
            df3<-df2 %>% group_by(Week = floor_date(df2$DATE-nth+1, unit="week")) %>% 
                summarize(week_steps=sum(steps), min_date = min(DATE), max_date = max(DATE),temp=mean(TAVG))
            df3$nth_week<-1:length(df3$Week)
            
            if(input$scale==FALSE){
                start_week<-ceiling(input$rate*length(df3$Week))
                end_week<-length(df3$Week)
                n<-end_week-start_week+1
                trainset<-df3[1:(start_week-1),]
                testset<-df3[start_week:end_week,]
                
                model_lm<-lm(data=trainset,week_steps~nth_week+temp)
                summary(model_lm)
                forecast_lsstep<-predict(model_lm,testset)
                
                autofit_arima<-auto.arima(model_lm$residuals)
                summary(autofit_arima)
                forecast_arimastep<-forecast(autofit_arima,h=n)
                forecast_arima<-forecast_arimastep$mean+forecast_lsstep
                
                autofit_structTS<-StructTS(model_lm$residuals)
                summary(autofit_structTS)
                forecast_StructTSstep<-forecast(autofit_structTS,h=n)
                forecast_StructTS<-forecast_StructTSstep$mean+forecast_lsstep
                if(input$summ=="arima"){
                    plot(df3$Week,df3$week_steps,type="l",col="grey54",xlab="Week",ylab="Weekly Steps",main=paste("Weekly Step Count for subject",input$y,sep=" "))
                    lines(df3$Week[start_week:end_week],forecast_arima,col="blue",lwd=1)
                    lines(df3$Week[start_week:end_week],forecast_lsstep+forecast_arimastep$upper[,1],col="green",lwd=1,lty=1)
                    lines(df3$Week[start_week:end_week],forecast_lsstep+forecast_arimastep$lower[,1],col="green",lwd=1,lty=1)
                    lines(df3$Week[start_week:end_week],forecast_lsstep+forecast_arimastep$upper[,2],col="pink",lwd=1,lty=1)
                    lines(df3$Week[start_week:end_week],forecast_lsstep+forecast_arimastep$lower[,2],col="pink",lwd=1,lty=1)
                    legend("topright", legend = c("true values", "forecasting values mean","forecasting value 80% confidence interval line","forecasting value 90% confidence interval line"), 
                           col = c("grey54", "blue","green","pink"), pch = c(46,46,46,46), lty = c(1,1,1,1),cex=0.6)
                }else if(input$summ=="struct"){
                    plot(df3$Week,df3$week_steps,type="l",col="grey54",xlab="Week",ylab="Weekly Steps",main=paste("Weekly Step Count for subject",input$y,sep=" "))
                    lines(df3$Week[start_week:end_week],forecast_StructTS,col="blue",lwd=1)
                    lines(df3$Week[start_week:end_week],forecast_lsstep+forecast_StructTSstep$upper[,1],col="green",lwd=1,lty=1)
                    lines(df3$Week[start_week:end_week],forecast_lsstep+forecast_StructTSstep$lower[,1],col="green",lwd=1,lty=1)
                    lines(df3$Week[start_week:end_week],forecast_lsstep+forecast_StructTSstep$upper[,2],col="pink",lwd=1,lty=1)
                    lines(df3$Week[start_week:end_week],forecast_lsstep+forecast_StructTSstep$lower[,2],col="pink",lwd=1,lty=1)
                    legend("topright", legend = c("true values", "forecasting values mean","forecasting value 80% confidence interval line","forecasting value 90% confidence interval line"), 
                           col = c("grey54", "blue","green","pink"), pch = c(46,46,46,46), lty = c(1,1,1,1),cex=0.6)
                    
                }else if(input$summ=="comp"){
                    
                    plot(df3$Week,df3$week_steps,type="l",col="grey54",xlab="Week",ylab="Weekly Steps",main=paste("Weekly Step Count for subject",input$y,sep=" "))
                    lines(df3$Week[start_week:end_week],forecast_arima,col="green",lwd=1)
                    lines(df3$Week[start_week:end_week],forecast_StructTS,col="red",lwd=1)
                    legend("topright", legend = c("true values", "forecasting values of ARIMA ","forecasting values of StructTS model"), 
                           col = c("grey54", "green","red"), pch = c(46,46,46), lty = c(1,1,1),cex=0.6)
                }else if(input$summ=="arac"){
                    arima_residuals<-model_lm$residuals-autofit_arima$fitted
                    par(mfrow=c(1,2))
                    acf(arima_residuals)
                    pacf(arima_residuals)
                    
                }else if(input$summ=="strac"){
                    structTS_residuals<-autofit_structTS$residuals
                    par(mfrow=c(1,2))
                    acf(structTS_residuals)
                    pacf(structTS_residuals)
                }else{
                    df3$week_steps<-scale(df3$week_steps)
                    sum_squared_error<-c(0,0)
                    for(week in start_week:end_week){
                        #split train set and test set for every run
                        train_order<-window(1:end_week,end=week-0.001)
                        test_order<-window(1:end_week,start=week,end=week+1-0.01)
                        train_set<-df3[train_order,]
                        test_set<-df3[test_order,]
                        #fit least squares linear model
                        model_lm<-lm(data=train_set,week_steps~nth_week+temp)
                        #fit arima model on residuals
                        autofit_arima<-auto.arima(model_lm$residuals)
                        #fit StructTS model on residuals
                        autofit_structTS<-StructTS(model_lm$residuals)
                        forecast_arimastep<-forecast(autofit_arima,h=1)
                        forecast_StructTSstep<-forecast(autofit_structTS,h=1)
                        lsstep<-predict(model_lm,test_set)
                        forecast_arima<-forecast_arimastep$mean+lsstep
                        forecast_StructTS<-forecast_StructTSstep$mean+lsstep
                        #calcuate the sum of squared errors
                        sum_squared_error[1]<-sum_squared_error[1]+sum((as.matrix(forecast_arima)-test_set$week_steps)^2)
                        sum_squared_error[2]<-sum_squared_error[2]+sum((as.matrix(forecast_StructTS)-test_set$week_steps)^2)}
                    
                    
                    mse<-sum_squared_error/n
                    barplot(c(ARIMA=mse[1],StructTS=mse[2]))
                }
                
            }else{
                df3$week_steps<-scale(df3$week_steps)
                start_week<-ceiling(input$rate*length(df3$Week))
                end_week<-length(df3$Week)
                n<-end_week-start_week+1
                trainset<-df3[1:(start_week-1),]
                testset<-df3[start_week:end_week,]
                
                model_lm<-lm(data=trainset,week_steps~nth_week+temp)
                summary(model_lm)
                forecast_lsstep<-predict(model_lm,testset)
                
                autofit_arima<-auto.arima(model_lm$residuals)
                summary(autofit_arima)
                forecast_arimastep<-forecast(autofit_arima,h=n)
                forecast_arima<-forecast_arimastep$mean+forecast_lsstep
                
                autofit_structTS<-StructTS(model_lm$residuals)
                summary(autofit_structTS)
                forecast_StructTSstep<-forecast(autofit_structTS,h=n)
                forecast_StructTS<-forecast_StructTSstep$mean+forecast_lsstep
                if(input$summ=="arima"){
                    plot(df3$Week,df3$week_steps,type="l",col="grey54",xlab="Week",ylab="Weekly Steps",main=paste("Weekly Step Count for subject",input$y,sep=" "))
                    lines(df3$Week[start_week:end_week],forecast_arima,col="blue",lwd=1)
                    lines(df3$Week[start_week:end_week],forecast_lsstep+forecast_arimastep$upper[,1],col="green",lwd=1,lty=1)
                    lines(df3$Week[start_week:end_week],forecast_lsstep+forecast_arimastep$lower[,1],col="green",lwd=1,lty=1)
                    lines(df3$Week[start_week:end_week],forecast_lsstep+forecast_arimastep$upper[,2],col="pink",lwd=1,lty=1)
                    lines(df3$Week[start_week:end_week],forecast_lsstep+forecast_arimastep$lower[,2],col="pink",lwd=1,lty=1)
                    legend("topright", legend = c("true values", "forecasting values mean","forecasting value 80% confidence interval line","forecasting value 90% confidence interval line"), 
                           col = c("grey54", "blue","green","pink"), pch = c(46,46,46,46), lty = c(1,1,1,1),cex=0.6)
                }else if(input$summ=="struct"){
                    plot(df3$Week,df3$week_steps,type="l",col="grey54",xlab="Week",ylab="Weekly Steps",main=paste("Weekly Step Count for subject",input$y,sep=" "))
                    lines(df3$Week[start_week:end_week],forecast_StructTS,col="blue",lwd=1)
                    lines(df3$Week[start_week:end_week],forecast_lsstep+forecast_StructTSstep$upper[,1],col="green",lwd=1,lty=1)
                    lines(df3$Week[start_week:end_week],forecast_lsstep+forecast_StructTSstep$lower[,1],col="green",lwd=1,lty=1)
                    lines(df3$Week[start_week:end_week],forecast_lsstep+forecast_StructTSstep$upper[,2],col="pink",lwd=1,lty=1)
                    lines(df3$Week[start_week:end_week],forecast_lsstep+forecast_StructTSstep$lower[,2],col="pink",lwd=1,lty=1)
                    legend("topright", legend = c("true values", "forecasting values mean","forecasting value 80% confidence interval line","forecasting value 90% confidence interval line"), 
                           col = c("grey54", "blue","green","pink"), pch = c(46,46,46,46), lty = c(1,1,1,1),cex=0.6)
                    
                }else if(input$summ=="comp"){
                    
                    plot(df3$Week,df3$week_steps,type="l",col="grey54",xlab="Week",ylab="Weekly Steps",main=paste("Weekly Step Count for subject",input$y,sep=" "))
                    lines(df3$Week[start_week:end_week],forecast_arima,col="green",lwd=1)
                    lines(df3$Week[start_week:end_week],forecast_StructTS,col="red",lwd=1)
                    legend("topright", legend = c("true values", "forecasting values of ARIMA ","forecasting values of StructTS model"), 
                           col = c("grey54", "green","red"), pch = c(46,46,46), lty = c(1,1,1),cex=0.6)
                }else if(input$summ=="arac"){
                    arima_residuals<-model_lm$residuals-autofit_arima$fitted
                    par(mfrow=c(1,2))
                    acf(arima_residuals)
                    pacf(arima_residuals)
                    
                }else if(input$summ=="strac"){
                    structTS_residuals<-autofit_structTS$residuals
                    par(mfrow=c(1,2))
                    acf(structTS_residuals)
                    pacf(structTS_residuals)
                }else{
                    
                    sum_squared_error<-c(0,0)
                    for(week in start_week:end_week){
                        #split train set and test set for every run
                        train_order<-window(1:end_week,end=week-0.001)
                        test_order<-window(1:end_week,start=week,end=week+1-0.01)
                        train_set<-df3[train_order,]
                        test_set<-df3[test_order,]
                        #fit least squares linear model
                        model_lm<-lm(data=train_set,week_steps~nth_week+temp)
                        #fit arima model on residuals
                        autofit_arima<-auto.arima(model_lm$residuals)
                        #fit StructTS model on residuals
                        autofit_structTS<-StructTS(model_lm$residuals)
                        forecast_arimastep<-forecast(autofit_arima,h=1)
                        forecast_StructTSstep<-forecast(autofit_structTS,h=1)
                        lsstep<-predict(model_lm,test_set)
                        forecast_arima<-forecast_arimastep$mean+lsstep
                        forecast_StructTS<-forecast_StructTSstep$mean+lsstep
                        #calcuate the sum of squared errors
                        sum_squared_error[1]<-sum_squared_error[1]+sum((as.matrix(forecast_arima)-test_set$week_steps)^2)
                        sum_squared_error[2]<-sum_squared_error[2]+sum((as.matrix(forecast_StructTS)-test_set$week_steps)^2)}
                    
                    
                    mse<-sum_squared_error/n
                    barplot(c(ARIMA=mse[1],StructTS=mse[2]))
                }
                
                
            }
            
        }else if(input$sep=="day"){
            
            
            
            day_count<-xdat() %>%group_by(date) %>%summarize(steps=sum(steps))
            day_count<-day_count[order(day_count$date),]
            
            
            day_ts<-zoo(day_count[,-1],order.by=as.Date(day_count$date, "%Y/%m/%d"))
            df<- merge(day_ts,zoo(,seq(start(day_ts),end(day_ts),by="day")), all=TRUE)
            index<-seq(start(day_ts),end(day_ts),by="day")
            df1<-data.frame(DATE=seq(start(day_ts),end(day_ts),by="day"),steps=df$steps)
            
            
            df2<-inner_join(df1,Temp(),by="DATE")
            
            
            df2$wday<-wday(df2$DATE)
            df2$t<-1:length(df2$wday)
            df2$is_not_weekend<-rep(0,length(df2$wday))
            df2$is_not_weekend[which(df2$wday==1|df2$wday==7)]<-1
            
            
            
            
            if(input$scale==FALSE){
                df2_with_all_NA<-df2
                df2_with_all_NA$steps<-scale(df2_with_all_NA$steps)
                start_day<-ceiling(input$rate*length(df2$wday))
                end_day<-length(df2$wday)
                n<-end_day-start_day+1
                id.na <- which(is.na(df2$steps))
                id.na1<-id.na[which(id.na<start_day)]
                df2$steps[id.na1]<-na_kalman(df2$steps)[id.na1]
                trainset<-df2[1:(start_day-1),]
                testset<-df2[start_day:end_day,]
                df2_with_all_NA$steps<-na_kalman(df2_with_all_NA$steps)
                model_lm<-lm(data=trainset,steps~t+wday+is_not_weekend+TAVG+TMAX+TMIN)
                
                
                forecast_lsstep<-predict(model_lm,df2[start_day:end_day,])
                
                autofit_arima<-auto.arima(model_lm$residuals)
                summary(autofit_arima)
                forecast_arimastep<-forecast(autofit_arima,h=n)
                forecast_arima<-forecast_arimastep$mean+forecast_lsstep
                
                autofit_structTS<-StructTS(model_lm$residuals)
                summary(autofit_structTS)
                forecast_StructTSstep<-forecast(autofit_structTS,h=n)
                forecast_StructTS<-forecast_StructTSstep$mean+forecast_lsstep
                if(input$summ=="arima"){
                    plot(seq(start(day_ts),end(day_ts),by="day"),df2$steps,type="l",col="grey54",xlab="day",ylab="steps",main=paste("Daily Step Count for subject",input$y,sep=" "))
                    points(index[id.na1],df2[id.na1,2],col="red",pch=1,lwd=0.1)
                    lines(index[start_day:end_day],forecast_arima,col="blue",lwd=1)
                    lines(index[start_day:end_day],forecast_lsstep+forecast_arimastep$upper[,1],col="green",lwd=1,lty=1)
                    lines(index[start_day:end_day],forecast_lsstep+forecast_arimastep$lower[,1],col="green",lwd=1,lty=1)
                    lines(index[start_day:end_day],forecast_lsstep+forecast_arimastep$upper[,2],col="pink",lwd=1,lty=1)
                    lines(index[start_day:end_day],forecast_lsstep+forecast_arimastep$lower[,2],col="pink",lwd=1,lty=1)
                    legend("topright", legend = c("true values", "imputed missing values","forecasting values mean","forecasting value 80% confidence interval line","forecasting value 90% confidence interval line"), 
                           col = c("grey54", "red","blue","green","pink"), pch = c(46, 1,46,46,46), lty = c(1, 0,1,1,1),cex=0.6)
                }else if(input$summ=="struct"){
                    plot(seq(start(day_ts),end(day_ts),by="day"),df2$steps,type="l",col="grey54",xlab="day",ylab="steps",main=paste("Daily Step Count for subject",input$y,sep=" "))
                    points(index[id.na1],df2[id.na1,2],col="red",pch=1,lwd=0.1)
                    lines(index[start_day:end_day],forecast_StructTS,col="blue",lwd=1)
                    lines(index[start_day:end_day],forecast_lsstep+forecast_StructTSstep$upper[,1],col="green",lwd=1,lty=1)
                    lines(index[start_day:end_day],forecast_lsstep+forecast_StructTSstep$lower[,1],col="green",lwd=1,lty=1)
                    lines(index[start_day:end_day],forecast_lsstep+forecast_StructTSstep$upper[,2],col="pink",lwd=1,lty=1)
                    lines(index[start_day:end_day],forecast_lsstep+forecast_StructTSstep$lower[,2],col="pink",lwd=1,lty=1)
                    legend("topright", legend = c("true values", "imputed missing values","forecasting values mean","forecasting value 80% confidence interval line","forecasting value 90% confidence interval line"), 
                           col = c("grey54", "red","blue","green","pink"), pch = c(46, 1,46,46,46), lty = c(1, 0,1,1,1),cex=0.6)
                    
                }else if(input$summ=="comp"){
                    
                    plot(seq(start(day_ts),end(day_ts),by="day"),df2$steps,type="l",col="grey54",xlab="day",ylab="steps",main=paste("Daily Step Count for subject",input$y,sep=" "))
                    points(index[id.na1],df2[id.na1,2],col="blue",pch=1,lwd=0.1)
                    lines(index[start_day:end_day],forecast_arima,col="green",lwd=1,lty=1)
                    lines(index[start_day:end_day],forecast_StructTS,col="red",lwd=1,lty=1)
                    legend("topright", legend = c("true values", "imputed missing values","forecasting values of ARIMA ","forecasting values of StructTS model"), 
                           col = c("grey54", "blue","green","red"), pch = c(46, 1,46,46), lty = c(1, 0,1,1),cex=0.6)
                }else if(input$summ=="arac"){
                    arima_residuals<-model_lm$residuals-autofit_arima$fitted
                    par(mfrow=c(1,2))
                    acf(arima_residuals)
                    pacf(arima_residuals)
                    
                }else if(input$summ=="strac"){
                    structTS_residuals<-autofit_structTS$residuals
                    par(mfrow=c(1,2))
                    acf(structTS_residuals)
                    pacf(structTS_residuals)
                }else{
                    
                    sum_squared_error<-c(0,0)
                    for(day in start_day:end_day){
                        #split train set and test set for every run
                        train_order<-window(1:end_day,end=day-0.001)
                        test_order<-window(1:end_day,start=day,end=day+1-0.01)
                        train_set<-df2_with_all_NA[train_order,]
                        test_set<-df2_with_all_NA[test_order,]
                        test_set
                        #fit least squares linear model
                        model_lm<-lm(data=train_set,steps~t+wday+is_not_weekend+TAVG+TMAX+TMIN)
                        #fit arima model on residuals
                        autofit_arima<-auto.arima(model_lm$residuals)
                        #fit StructTS model on residuals
                        autofit_structTS<-StructTS(model_lm$residuals)
                        forecast_arimastep<-forecast(autofit_arima,h=1)
                        forecast_StructTSstep<-forecast(autofit_structTS,h=1)
                        lsstep<-predict(model_lm,test_set)
                        forecast_arima<-forecast_arimastep$mean+lsstep
                        forecast_StructTS<-forecast_StructTSstep$mean+lsstep
                        #calcuate the sum of squared errors
                        sum_squared_error[1]<-sum_squared_error[1]+sum((as.matrix(forecast_arima)-test_set$steps)^2)
                        sum_squared_error[2]<-sum_squared_error[2]+sum((as.matrix(forecast_StructTS)-test_set$steps)^2)}
                    sum_squared_error/n
                    
                    mse<-sum_squared_error/n
                    
                    
                    mse<-sum_squared_error/n
                    barplot(c(ARIMA=mse[1],StructTS=mse[2]))
                }
                
                
                
            }else{
                df2$steps<-scale(df2$steps)
                df2_with_all_NA<-df2
                start_day<-ceiling(input$rate*length(df2$wday))
                end_day<-length(df2$wday)
                n<-end_day-start_day+1
                id.na <- which(is.na(df2$steps))
                id.na1<-id.na[which(id.na<start_day)]
                df2$steps[id.na1]<-na_kalman(df2$steps)[id.na1]
                trainset<-df2[1:(start_day-1),]
                testset<-df2[start_day:end_day,]
                df2_with_all_NA$steps<-na_kalman(df2_with_all_NA$steps)
                model_lm<-lm(data=trainset,steps~t+wday+is_not_weekend+TAVG+TMAX+TMIN)
                
                
                forecast_lsstep<-predict(model_lm,df2[start_day:end_day,])
                
                autofit_arima<-auto.arima(model_lm$residuals)
                summary(autofit_arima)
                forecast_arimastep<-forecast(autofit_arima,h=n)
                forecast_arima<-forecast_arimastep$mean+forecast_lsstep
                
                autofit_structTS<-StructTS(model_lm$residuals)
                summary(autofit_structTS)
                forecast_StructTSstep<-forecast(autofit_structTS,h=n)
                forecast_StructTS<-forecast_StructTSstep$mean+forecast_lsstep
                
                if(input$summ=="arima"){
                    plot(seq(start(day_ts),end(day_ts),by="day"),df2$steps,type="l",col="grey54",xlab="day",ylab="steps",main=paste("Daily Step Count for subject",input$y,sep=" "))
                    points(index[id.na1],df2[id.na1,2],col="red",pch=1,lwd=0.1)
                    lines(index[start_day:end_day],forecast_arima,col="blue",lwd=1)
                    lines(index[start_day:end_day],forecast_lsstep+forecast_arimastep$upper[,1],col="green",lwd=1,lty=1)
                    lines(index[start_day:end_day],forecast_lsstep+forecast_arimastep$lower[,1],col="green",lwd=1,lty=1)
                    lines(index[start_day:end_day],forecast_lsstep+forecast_arimastep$upper[,2],col="pink",lwd=1,lty=1)
                    lines(index[start_day:end_day],forecast_lsstep+forecast_arimastep$lower[,2],col="pink",lwd=1,lty=1)
                    legend("topright", legend = c("true values", "imputed missing values","forecasting values mean","forecasting value 80% confidence interval line","forecasting value 90% confidence interval line"), 
                           col = c("grey54", "red","blue","green","pink"), pch = c(46, 1,46,46,46), lty = c(1, 0,1,1,1),cex=0.6)
                }else if(input$summ=="struct"){
                    plot(seq(start(day_ts),end(day_ts),by="day"),df2$steps,type="l",col="grey54",xlab="day",ylab="steps",main=paste("Daily Step Count for subject",input$y,sep=" "))
                    points(index[id.na1],df2[id.na1,2],col="red",pch=1,lwd=0.1)
                    lines(index[start_day:end_day],forecast_StructTS,col="blue",lwd=1)
                    lines(index[start_day:end_day],forecast_lsstep+forecast_StructTSstep$upper[,1],col="green",lwd=1,lty=1)
                    lines(index[start_day:end_day],forecast_lsstep+forecast_StructTSstep$lower[,1],col="green",lwd=1,lty=1)
                    lines(index[start_day:end_day],forecast_lsstep+forecast_StructTSstep$upper[,2],col="pink",lwd=1,lty=1)
                    lines(index[start_day:end_day],forecast_lsstep+forecast_StructTSstep$lower[,2],col="pink",lwd=1,lty=1)
                    legend("topright", legend = c("true values", "imputed missing values","forecasting values mean","forecasting value 80% confidence interval line","forecasting value 90% confidence interval line"), 
                           col = c("grey54", "red","blue","green","pink"), pch = c(46, 1,46,46,46), lty = c(1, 0,1,1,1),cex=0.6)
                    
                }else if(input$summ=="comp"){
                    
                    plot(seq(start(day_ts),end(day_ts),by="day"),df2$steps,type="l",col="grey54",xlab="day",ylab="steps",main=paste("Daily Step Count for subject",input$y,sep=" "))
                    points(index[id.na1],df2[id.na1,2],col="blue",pch=1,lwd=0.1)
                    lines(index[start_day:end_day],forecast_arima,col="green",lwd=1,lty=1)
                    lines(index[start_day:end_day],forecast_StructTS,col="red",lwd=1,lty=1)
                    legend("topright", legend = c("true values", "imputed missing values","forecasting values of ARIMA ","forecasting values of StructTS model"), 
                           col = c("grey54", "blue","green","red"), pch = c(46, 1,46,46), lty = c(1, 0,1,1),cex=0.6)
                }else if(input$summ=="arac"){
                    arima_residuals<-model_lm$residuals-autofit_arima$fitted
                    par(mfrow=c(1,2))
                    acf(arima_residuals)
                    pacf(arima_residuals)
                    
                }else if(input$summ=="strac"){
                    structTS_residuals<-autofit_structTS$residuals
                    par(mfrow=c(1,2))
                    acf(structTS_residuals)
                    pacf(structTS_residuals)
                }else{
                    
                    sum_squared_error<-c(0,0)
                    for(day in start_day:end_day){
                        #split train set and test set for every run
                        train_order<-window(1:end_day,end=day-0.001)
                        test_order<-window(1:end_day,start=day,end=day+1-0.01)
                        train_set<-df2_with_all_NA[train_order,]
                        test_set<-df2_with_all_NA[test_order,]
                        test_set
                        #fit least squares linear model
                        model_lm<-lm(data=train_set,steps~t+wday+is_not_weekend+TAVG+TMAX+TMIN)
                        #fit arima model on residuals
                        autofit_arima<-auto.arima(model_lm$residuals)
                        #fit StructTS model on residuals
                        autofit_structTS<-StructTS(model_lm$residuals)
                        forecast_arimastep<-forecast(autofit_arima,h=1)
                        forecast_StructTSstep<-forecast(autofit_structTS,h=1)
                        lsstep<-predict(model_lm,test_set)
                        forecast_arima<-forecast_arimastep$mean+lsstep
                        forecast_StructTS<-forecast_StructTSstep$mean+lsstep
                        #calcuate the sum of squared errors
                        sum_squared_error[1]<-sum_squared_error[1]+sum((as.matrix(forecast_arima)-test_set$steps)^2)
                        sum_squared_error[2]<-sum_squared_error[2]+sum((as.matrix(forecast_StructTS)-test_set$steps)^2)}
                    sum_squared_error/n
                    
                    mse<-sum_squared_error/n
                    
                    
                    mse<-sum_squared_error/n
                    barplot(c(ARIMA=mse[1],StructTS=mse[2]))
                }
                
                
                
                
            }
            
        }else{
            Sys.setlocale("LC_TIME", "English")
            
            days<-split.data.frame(xdat(),xdat()$date)
            date_day<-unique(xdat()$date)
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
                
                days[[i]]<-df1 %>% group_by(groups = floor_date(df1$DATE, unit="15 min")) %>% summarize(group_steps=sum(steps))
            }
            
            dat<-days%>%reduce(rbind)
            dat$DATE<-as.Date(dat$groups)
            
            day_count<-xdat()%>%group_by(date) %>%summarize(steps=sum(steps))
            day_count<-day_count[order(day_count$date),]
            day_ts<-zoo(day_count[,-1],order.by=as.Date(day_count$date, "%Y/%m/%d"))
            df<- merge(day_ts,zoo(,seq(start(day_ts),end(day_ts),by="day")), all=TRUE)
            index<-seq(start(day_ts),end(day_ts),by="day")
            df1<-data.frame(DATE=seq(start(day_ts),end(day_ts),by="day"),steps=df$steps)
            df2<-inner_join(df1,Temp(),by="DATE")
            df2$wday<-wday(df2$DATE)
            df2$t<-1:length(df2$wday)
            dat<-inner_join(dat,df2,by="DATE")
            dat_new<-dat
            
            
            if(input$scale==FALSE){
                start<-ceiling(input$rate*length(dat_new$groups))
                end<-length(dat_new$groups)
                n<-end-start+1
                dat_new$wday<-as.character(dat_new$wday)
                trainset<-dat_new[1:(start-1),]
                testset<-dat_new[start:end,]
                
                model_lm<-lm(data=trainset,group_steps~t+wday+TAVG+TMAX+TMIN)
                
                
                forecast_lsstep<-predict(model_lm,testset)
                
                autofit_arima<-auto.arima(model_lm$residuals)
                summary(autofit_arima)
                forecast_arimastep<-forecast(autofit_arima,h=n)
                forecast_arima<-forecast_arimastep$mean+forecast_lsstep
                
                autofit_structTS<-StructTS(model_lm$residuals)
                summary(autofit_structTS)
                forecast_StructTSstep<-forecast(autofit_structTS,h=n)
                forecast_StructTS<-forecast_StructTSstep$mean+forecast_lsstep
                if(input$summ=="arima"){
                    plot(testset$groups,testset$group_steps,type="l",   xlab="15 minutes",ylab="15 minutes Steps",main=paste("15 minutes Step Count for subject",input$y,sep=" "))
                    lines(testset$groups,forecast_arima,col="blue",lwd=1)
                    lines(testset$groups,forecast_lsstep+forecast_arimastep$upper[,1],col="green",lwd=1,lty=1)
                    lines(testset$groups,forecast_lsstep+forecast_arimastep$lower[,1],col="green",lwd=1,lty=1)
                    lines(testset$groups,forecast_lsstep+forecast_arimastep$upper[,2],col="pink",lwd=1,lty=1)
                    lines(testset$groups,forecast_lsstep+forecast_arimastep$lower[,2],col="pink",lwd=1,lty=1)
                    legend("topright", legend = c("true values", "forecasting values mean","forecasting value 80% confidence interval line","forecasting value 90% confidence interval line"), 
                           col = c("grey54", "blue","green","pink"), pch = c(46,46,46,46), lty = c(1,1,1,1),cex=0.6)
                }else if(input$summ=="struct"){
                    plot(testset$groups,testset$group_steps,type="l",   xlab="15 minutes",ylab="15 minutes Steps",main=paste("15 minutes Step Count for subject",input$y,sep=" "))
                    lines(testset$groups,forecast_StructTS,col="blue",lwd=1)
                    lines(testset$groups,forecast_lsstep+forecast_StructTSstep$upper[,1],col="green",lwd=1,lty=1)
                    lines(testset$groups,forecast_lsstep+forecast_StructTSstep$lower[,1],col="green",lwd=1,lty=1)
                    lines(testset$groups,forecast_lsstep+forecast_StructTSstep$upper[,2],col="pink",lwd=1,lty=1)
                    lines(testset$groups,forecast_lsstep+forecast_StructTSstep$lower[,2],col="pink",lwd=1,lty=1)
                    legend("topright", legend = c("true values", "forecasting values mean","forecasting value 80% confidence interval line","forecasting value 90% confidence interval line"), 
                           col = c("grey54", "blue","green","pink"), pch = c(46,46,46,46), lty = c(1,1,1,1),cex=0.6)
                    
                }else if(input$summ=="comp"){
                    
                    plot(testset$groups,testset$group_steps,type="l",col="grey54",xlab="Day",ylab="Weekly Steps",main=paste("Weekly Step Count for subject",input$y,sep=" "))
                    lines(testset$groups,forecast_arima,col="green",lwd=1)
                    lines(testset$groups,forecast_StructTS,col="red",lwd=1)
                    legend("topright", legend = c("true values", "forecasting values of ARIMA ","forecasting values of StructTS model"), 
                           col = c("grey54", "green","red"), pch = c(46,46,46), lty = c(1,1,1),cex=0.6)
                }else if(input$summ=="arac"){
                    arima_residuals<-model_lm$residuals-autofit_arima$fitted
                    par(mfrow=c(1,2))
                    acf(arima_residuals)
                    pacf(arima_residuals)
                }else if(input$summ=="strac"){
                    structTS_residuals<-autofit_structTS$residuals
                    par(mfrow=c(1,2))
                    acf(structTS_residuals)
                    pacf(structTS_residuals)
                }else{
                    sum_squared_error<-c(0,0)
                    sum_squared_error[1]<-sum((testset$group_steps-forecast_arima)^2)
                    sum_squared_error[2]<-sum((testset$group_steps-forecast_StructTS)^2)
                    mse<-sum_squared_error/n
                    barplot(c(ARIMA=mse[1],StructTS=mse[2]))
                    
                }
                
                
                
                
                
            }else{
                dat_new$group_steps<-scale(dat_new$group_steps)
                start<-ceiling(input$rate*length(dat_new$groups))
                end<-length(dat_new$groups)
                n<-end-start+1
                dat_new$wday<-as.character(dat_new$wday)
                trainset<-dat_new[1:(start-1),]
                testset<-dat_new[start:end,]
                
                model_lm<-lm(data=trainset,group_steps~t+wday+TAVG+TMAX+TMIN)
                
                
                forecast_lsstep<-predict(model_lm,testset)
                
                autofit_arima<-auto.arima(model_lm$residuals)
                summary(autofit_arima)
                forecast_arimastep<-forecast(autofit_arima,h=n)
                forecast_arima<-forecast_arimastep$mean+forecast_lsstep
                
                autofit_structTS<-StructTS(model_lm$residuals)
                summary(autofit_structTS)
                forecast_StructTSstep<-forecast(autofit_structTS,h=n)
                forecast_StructTS<-forecast_StructTSstep$mean+forecast_lsstep
                if(input$summ=="arima"){
                    plot(testset$groups,testset$group_steps,type="l",   xlab="15 minutes",ylab="15 minutes Steps",main=paste("15 minutes Step Count for subject",input$y,sep=" "))
                    lines(testset$groups,forecast_arima,col="blue",lwd=1)
                    lines(testset$groups,forecast_lsstep+forecast_arimastep$upper[,1],col="green",lwd=1,lty=1)
                    lines(testset$groups,forecast_lsstep+forecast_arimastep$lower[,1],col="green",lwd=1,lty=1)
                    lines(testset$groups,forecast_lsstep+forecast_arimastep$upper[,2],col="pink",lwd=1,lty=1)
                    lines(testset$groups,forecast_lsstep+forecast_arimastep$lower[,2],col="pink",lwd=1,lty=1)
                    legend("topright", legend = c("true values", "forecasting values mean","forecasting value 80% confidence interval line","forecasting value 90% confidence interval line"), 
                           col = c("grey54", "blue","green","pink"), pch = c(46,46,46,46), lty = c(1,1,1,1),cex=0.6)
                }else if(input$summ=="struct"){
                    plot(testset$groups,testset$group_steps,type="l",   xlab="15 minutes",ylab="15 minutes Steps",main=paste("15 minutes Step Count for subject",input$y,sep=" "))
                    lines(testset$groups,forecast_StructTS,col="blue",lwd=1)
                    lines(testset$groups,forecast_lsstep+forecast_StructTSstep$upper[,1],col="green",lwd=1,lty=1)
                    lines(testset$groups,forecast_lsstep+forecast_StructTSstep$lower[,1],col="green",lwd=1,lty=1)
                    lines(testset$groups,forecast_lsstep+forecast_StructTSstep$upper[,2],col="pink",lwd=1,lty=1)
                    lines(testset$groups,forecast_lsstep+forecast_StructTSstep$lower[,2],col="pink",lwd=1,lty=1)
                    legend("topright", legend = c("true values", "forecasting values mean","forecasting value 80% confidence interval line","forecasting value 90% confidence interval line"), 
                           col = c("grey54", "blue","green","pink"), pch = c(46,46,46,46), lty = c(1,1,1,1),cex=0.6)
                    
                }else if(input$summ=="comp"){
                    
                    plot(testset$groups,testset$group_steps,type="l",col="grey54",xlab="Day",ylab="Weekly Steps",main=paste("Weekly Step Count for subject",input$y,sep=" "))
                    lines(testset$groups,forecast_arima,col="green",lwd=1)
                    lines(testset$groups,forecast_StructTS,col="red",lwd=1)
                    legend("topright", legend = c("true values", "forecasting values of ARIMA ","forecasting values of StructTS model"), 
                           col = c("grey54", "green","red"), pch = c(46,46,46), lty = c(1,1,1),cex=0.6)
                }else if(input$summ=="arac"){
                    arima_residuals<-model_lm$residuals-autofit_arima$fitted
                    par(mfrow=c(1,2))
                    acf(arima_residuals)
                    pacf(arima_residuals)
                }else if(input$summ=="strac"){
                    structTS_residuals<-autofit_structTS$residuals
                    par(mfrow=c(1,2))
                    acf(structTS_residuals)
                    pacf(structTS_residuals)
                }else{
                    sum_squared_error<-c(0,0)
                    sum_squared_error[1]<-sum((testset$group_steps-forecast_arima)^2)
                    sum_squared_error[2]<-sum((testset$group_steps-forecast_StructTS)^2)
                    mse<-sum_squared_error/n
                    barplot(c(ARIMA=mse[1],StructTS=mse[2]))
                    
                }
                
                
                
                
            }
        }
        
    })
    
}
# Run the app ----
shinyApp(ui, server)