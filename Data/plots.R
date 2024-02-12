

## Example of the weighted function
w <- function(x){ # Pondering function (it is a logistic regression with a fix)
  b_1 = log(1000-1)/20 # param[5] is the delay
  b_0 = -b_1*50 # param[6] is the A/C boundary
  w <- 1/ (1 + exp(-b_0-b_1*x))
  return(w)
}


x<- seq(0,100,.1)

ws <- w(x)


plot(x,ws,type='l',col=rgb(0,0,0,.8),ylab ='Proportion' ,xlab='Depth',main='Weighted Function')
segments(x0 = 50,x1 = 50,y0 = -.1,y1=w(50),col=rgb(1,0,0,.8))
segments(x0 = 30,x1 = 30,y0 = -.1,y1=w(30),col=rgb(1,0,0,.8))
segments(x0 = 70,x1 = 70,y0 = -.1,y1=w(70),col=rgb(1,0,0,.8))
pdf('~/GitHub/Bayesian_Carbon_Acc/Manuscript/Figures/weighted_fun.pdf')
dev.off()