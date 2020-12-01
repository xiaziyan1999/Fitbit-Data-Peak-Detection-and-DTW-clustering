#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)

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
            selectInput("ID",
                        "Subject:",
                        choices = c(1:98),
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
                        value =20 ),
            sliderInput("thresh",
                        label = "Threshold:",
                        min = 1,
                        max = 5,
                        step=0.5,
                        value =3),
            sliderInput("inf",
                        label = "Influence:",
                        min = 0,
                        max = 1,
                        step=0.01,
                        value =0.1 ),
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
        x<-fit_split[[input$y]]
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
        lag<-input$l
        threshold<-input$thresh
        influence<-input$inf
        result <- ThresholdingAlgo(y,lag,threshold,influence)
        
        # Plot result
        par(mfrow = c(2,1),oma = c(2,2,0,0) + 0.1,mar = c(0,0,2,1) + 0.2)
        plot(df1$DATE,y,type="l",ylab="",xlab="") 
        lines(df1$DATE,result$avgFilter,type="l",col="cyan",lwd=2)
        lines(df1$DATE,result$avgFilter+threshold*result$stdFilter,type="l",col="green",lwd=2)
        lines(df1$DATE,result$avgFilter-threshold*result$stdFilter,type="l",col="green",lwd=2)
        plot(df1$DATE,result$signals,type="S",col="red",ylab="",xlab="",ylim=c(-1.5,1.5),lwd=2)
    })

}

# Run the application 
shinyApp(ui = ui, server = server)
