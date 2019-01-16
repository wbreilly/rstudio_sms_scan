#  demos

x=abs(rnorm(100*100,50,25))
x=matrix(x,nrow=100)

require(reshape2)
require(ggplot2)
x1=melt(x)
names(x1)=c("x","y","color")

x1$color=factor(x1$color>50)
levels(x1$color)=c("lessthan50","more than 50")

qplot(x, y, fill=color, data=x1,geom='tile')
