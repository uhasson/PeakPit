`get.visibility` <-
function(y,x=1:length(y),size=Inf,sets=F)
{
#Created by Thomas Jagger Version 1.0
#For support send email to tjagger@fsu.edu
#We assume that y and x are pairs, with x being the independent vector, and y the dependent vector
#We form slopes to each item and diffrerence the slopes from both sides to get the horizon.
#All points inside the horizen are kept.
#returns span. You can create the actual list from span as
#size: length to search in forward, back directions.
ly<-length(y)
ly1<-ly-1
tx<-x
ty<-y
val<-vector(mode="list",length=ly)
mss<-ly:2
if(size != Inf)
mss<-(size-mss)*(mss > size) + mss
for(i in 1:ly1)
 {
 ms<-mss[i]
 xx<-tx[2:ms]
 yy<-ty[2:ms]
 slopes<-(yy-ty[1])/(xx-tx[1])
 val[[i]]<-which(cummax(slopes)==slopes)+i
 tx<-tx[-1]
 ty<-ty[-1]
 }
  val.sparse<-cbind(row=rep(1:length(val),sapply(val,function(x) length(x))),col=unlist(val))
  val.sparse<-rbind(val.sparse,val.sparse[,c(2,1)])
  val.sparse<-val.sparse[order(val.sparse[,1]),]
  by.node<-split(val.sparse[,2],val.sparse[,1])
  k<-sapply(by.node,length)
  degree.dist<-table(k)
  degree.dist<-data.frame(k=as.numeric(names(degree.dist)),degree=as.vector(degree.dist))
  degree.dist<-cbind(degree.dist,P=degree.dist[,2]/ly)
  return(list(sm=val.sparse,node=by.node,pk=degree.dist))
}
