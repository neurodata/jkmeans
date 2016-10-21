setwd("~/git/jkmeans/Rmd/")

library(shiny)
library(datasets)
library(ggplot2)

# Define server logic required to plot various variables against mpg
shinyServer(function(input, output) {
  
  loadResult<- reactive({
    K<- input$K
    balanced<- input$weight
    sigma<- input$sigma
    initial_truth<- input$initial_truth
    
    filename = paste("sims/result_",paste(K,balanced ,sigma , initial_truth,sep =  "_"),".Rda",sep="")
    
    K<- as.numeric(K)
    
    load(file = filename)
    result
  }
  )
  
  
  
  output$ratePlotMCE <- renderPlot({
    
    result<- loadResult()
    
    K<- as.numeric(input$K)
    
    nSeries<- c(seq(20,200,by = 20),300,400,500)
    
    plotDF<- data.frame( "n" = rep(nSeries, K*2),
                         "MC14Error"= c(result$MCEmean),
                         "MC14ErrorL"= c(result$MCEq25),
                         "MC14ErrorU"= c(result$MCEq975),
                         "J"= as.factor(rep(rep(c(1:K),2),each= length(nSeries))),
                         "Model"=as.factor( rep(c("fixW","flexW"),each=length(nSeries)*K)))
    
    plotRMSE<- data.frame( "n" = rep(nSeries, K*2),
                           "RMSE"= c(result$RMSEmean),
                           "RMSEL"= c(result$RMSEq25),
                           "RMSEU"= c(result$RMSEq975),
                           "J"= as.factor(rep(rep(c(1:K),2),each= length(nSeries))),
                           "Model"=as.factor( rep(c("fixW","flexW"),each=length(nSeries)*K)))
    
    
    if(!input$fixW){
      plotDF<- plotDF[plotDF$Model!= "fixW",]
      plotRMSE<- plotRMSE[plotRMSE$Model!= "fixW",]
      
    }
    
    if(!input$GMM){
      plotDF<- plotDF[plotDF$Model!= "flexW",]
      plotRMSE<- plotRMSE[plotRMSE$Model!= "flexW",]
    }
    
    
    p<- ggplot(data=plotDF, aes())
    p<- p+  geom_line(aes(x= n, y=MC14Error, colour=J,linetype=Model ))
    
    # pE<- p+  geom_line(aes(x= n, y=MC14ErrorL, colour=J),linetype=2)
    # pE<- pE+  geom_line(aes(x= n, y=MC14ErrorU, colour=J),linetype=2)
    
    p+ theme_bw() + geom_hline(yintercept = min(result$BayesError))
    
    
  })
  
  output$ratePlotRMSE <- renderPlot({
    
    result<- loadResult()
    
    K<- as.numeric(input$K)
    
    nSeries<- c(seq(20,200,by = 20),300,400,500)
    
    plotDF<- data.frame( "n" = rep(nSeries, K*2),
                         "MC14Error"= c(result$MCEmean),
                         "MC14ErrorL"= c(result$MCEq25),
                         "MC14ErrorU"= c(result$MCEq975),
                         "J"= as.factor(rep(rep(c(1:K),2),each= length(nSeries))),
                         "Model"=as.factor( rep(c("fixW","flexW"),each=length(nSeries)*K)))
    
    plotRMSE<- data.frame( "n" = rep(nSeries, K*2),
                           "RMSE"= c(result$RMSEmean),
                           "RMSEL"= c(result$RMSEq25),
                           "RMSEU"= c(result$RMSEq975),
                           "J"= as.factor(rep(rep(c(1:K),2),each= length(nSeries))),
                           "Model"=as.factor( rep(c("fixW","flexW"),each=length(nSeries)*K)))
    
    
    if(!input$fixW){
      plotDF<- plotDF[plotDF$Model!= "fixW",]
      plotRMSE<- plotRMSE[plotRMSE$Model!= "fixW",]
      
    }
    
    if(!input$GMM){
      plotDF<- plotDF[plotDF$Model!= "flexW",]
      plotRMSE<- plotRMSE[plotRMSE$Model!= "flexW",]
    }
    
    
    p<- ggplot(data=plotRMSE, aes())
    p<- p+  geom_line(aes(x= n, y=RMSE, colour=J,linetype=Model ))
    
    pE<- p+  geom_line(aes(x= n, y=RMSEL, colour=Model),linetype=2)
    pE<- pE+  geom_line(aes(x= n, y=RMSEU, colour=Model),linetype=2)
    
    p+ theme_bw()
  })
  
  
})