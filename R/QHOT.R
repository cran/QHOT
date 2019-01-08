#' @title QTL Hotspot Detection
#' @description This function produces both the numerical and graphical summaries of the QTL hotspot detection in the genomes that are available on the worldwide web including the flanking markers of QTLs.
#' Man-Hsia Yang, Dong-Hong Wu, Chen-Hung Kao. 2019. A Statistical Procedure for Genome-wide Detection of QTL Hotspots Using Public Databases with Application to Rice. G3-Genes Genom Genet DOI: 10.1534/g3.118.200922.
#' @param DataQTL a data-frame of values for QTL information including the trait names, which chromosomes localized, the left and right flanking marker positions of QTLs for the first to fourth columns, respectively.
#' @param DataCrop	a data-frame of values for chromosome information consisting of the names, center positions and lengths of chromosomes for the first to third columns, respectively.
#' @param ScanStep a value for the length of every bin.
#' @param NH a value for the number of spurious hotspots in the proposed method.
#' @param NP a value for permutation times to calculate the threshold.
#' @references Man-Hsia Yang, Dong-Hong Wu, Chen-Hung Kao. 2019. A Statistical Procedure for Genome-wide Detection of QTL Hotspots Using Public Databases with Application to Rice. G3-Genes Genom Genet DOI: 10.1534/g3.118.200922.
#' @examples
#' Trait<-paste("t",sample(1:9,100,replace=TRUE,prob=c(2,rep(1,8))/10),sep="")
#' chr<-1
#' L<-sample(seq(0,90,by=10),100,replace=TRUE,prob=c(0.5,0.5,5.5,rep(0.5,7))/10)
#' R<-L+sample(c(0.5,1,5,10,50),100,replace=TRUE)
#' R[R>100]<-100
#' DataQTL.t<-data.frame(Trait,chr,L,R)
#' DataCrop.t<-data.frame(chr=1,center=75,length=100)
#' QHOT(DataQTL.t, DataCrop.t, ScanStep=0.5, NH=1, NP=1000)
#' @import grDevices graphics
#' @export

QHOT<-function(DataQTL,DataCrop,ScanStep,NH,NP){



Data.freq<-function(DataQTL,ScanStep){

  Len.DataQTL<-nrow(DataQTL)
  Data.h<-c()
  Na.h<-c()


X.max<-ceiling(max(DataQTL)/ScanStep)*ScanStep

  Break<-seq(0,X.max,by=ScanStep)

  for(j in 1:Len.DataQTL){

    M1<-DataQTL[j,1]
    M2<-DataQTL[j,2]
    Break1<-Break[max(which(Break<=M1)):min(which(Break>=M2))]


      if(length(Break1)<=2){
         Data.h<-c(Data.h,1)
         Na.h<-c(Na.h,Break1[1])
      }else{
         Temp.break<-Break[min(which(Break>M1)):max(which(Break<M2))]
         Len.QTL<-diff(c(M1,M2))
         Data.h<-c(Data.h,diff(sort(c(Temp.break,M1,M2)))/Len.QTL)
         Na.h<-c(Na.h,Break1[1:length(Break1)-1])
      }

  }
  Data.h<-c(Data.h,rep(0,times=length(Break)-1))
  Na.h<-c(Na.h,Break[1:length(Break)-1])
  Freq.h<-tapply(Data.h,Na.h,sum)
  Break.h1<-levels(cut(1,breaks=Break,right=F))
  data.frame(BreakPoint=Break.h1,Left=Break[1:length(Break)-1],Right=Break[2:length(Break)],ExpectedNumber=Freq.h)
}

PHot<-function(Data,NP,NH){
  L<-dim(Data)[1]
  Trait.no<-dim(Data)[2]
  t1<-Sys.time()
  TV.p<-c()
  ps<-c()
  Data2<-Data
  for(k in 1:NP){
    for(j in 1:Trait.no){
        Index<-sample(1:L,L)
        Data2[,j]<- Data[Index,j]
    }
    tmp<-sort(apply(Data2,1,sum),decreasing=TRUE)
    ps<-c(ps,tmp)
  }
 ps2<-matrix(ps,byrow=T,nrow=NP)
 TV.p<-apply(ps2,2,sort)[round(NP*0.95),1:NH]

 t2<-Sys.time()
 Output<-list(time=t2-t1,TV.p=TV.p)
 invisible(Output)

}

  NoQTL<-table(DataQTL[,2])
  Na.chr<-names(NoQTL)

  Na.chr.crop<-factor(DataCrop[,1])
  Index2<-Na.chr.crop%in%Na.chr
  DataCrop2<-DataCrop[Index2,]
  Na.chr.crop2<-Na.chr.crop[Index2]
  Index3<-order(Na.chr.crop2)
  DataCrop3<-DataCrop2[Index3,]
  ChrL<-DataCrop3[,3]

  Na.trait<-levels(factor(as.vector(DataQTL[,1])))
  N.trait<-length(Na.trait)

  DataF<-vector("list",length(NoQTL))
  DataF.P<-c()
  DataF.Q<-c()
  L.ALL<-ceiling(DataCrop[,3]/ScanStep)

  L.ALL2<-sum(L.ALL)


for(k in 1:length(NoQTL)){

  DataQTL2<-DataQTL[DataQTL[,2]==Na.chr[k],3:4]
  DataF[[k]]<-Data.freq(DataQTL2,ScanStep)
  L.break<-dim(DataF[[k]])[1]
  for(t in 1:N.trait){

      DataQTL.t<-DataQTL[DataQTL[,2]==Na.chr[k]&DataQTL[,1]==Na.trait[t],3:4]

    if(nrow(DataQTL.t)!=0){
      Temp<-Data.freq(DataQTL.t,ScanStep)$ExpectedNumber
      L.diff<-L.break-length(Temp)
      Temp2<-c(Temp,rep(0,L.diff))}else{
      Temp2<-0
      }
      DataF[[k]]<-data.frame(DataF[[k]],Temp2)
      names(DataF[[k]])[ncol(DataF[[k]])]<-Na.trait[t]

  }



  for(tt in 1:dim(DataQTL2)[1]){
      Temp<-Data.freq(DataQTL2[tt,],ScanStep)$ExpectedNumber
      DataF.Q<-c(DataF.Q,Temp[Temp>0])
  }


  DataF.P<-rbind(DataF.P,DataF[[k]][,4:(4+N.trait)])
}


DataF.tmp<-data.frame(matrix(0,nrow=(L.ALL2-dim(DataF.P)[1]),ncol=(N.trait+1)))
names(DataF.tmp)<-names(DataF.P)
DataF.P2<-rbind(DataF.P,DataF.tmp)

Output<-PHot(DataF.P2[,2:(N.trait+1)],NP,NH)
TV.P2<-Output$TV.p

TV.Q<-c()
for(i in 1:NP){
  Index.Q<-sample(1:L.ALL2, length(DataF.Q),replace=T)
  Index.Q2<-c(Index.Q,1:L.ALL2)
  DataF.Q2<-c(DataF.Q,rep(0,L.ALL2))
  TV.Q<-c(TV.Q,max(tapply(DataF.Q2,Index.Q2,sum)))
}
TV.Q2<-sort(TV.Q)[ceiling(NP*0.95)]


DataHotspot<-vector("list", length(NoQTL))
TableH<-c()
for(i in 1:length(NoQTL)){

   Data.Freq<-DataF[[i]][,4]
   Data.Left<-DataF[[i]][,2]
   DataHotspot[[i]]<-vector("list",(NH+1))
   for(nh in 1:NH){
     Index.P<-Data.Freq>=TV.P2[nh]
     DataHotspot[[i]][[nh]]<-Data.Left[Index.P]
     TableH<-c(TableH,sum(Index.P))
   }
   Index.Q<-Data.Freq>=TV.Q2
   DataHotspot[[i]][[(NH+1)]]<-Data.Left[Index.Q]
   TableH<-c(TableH,sum(Index.Q))
}

TableH2<-matrix(TableH,byrow=F,nrow=(NH+1))

Ymax<-max(DataF.P2[,1])

dev.new(7,5,record=TRUE)
par(mar=rep(0.1,4))

x1<-max(ChrL)*(-0.05)
x2<-max(ChrL)*1.1
y1<-max(Ymax)*(-0.4)
y2<-max(Ymax)*1.1


for(i in 1:length(NoQTL)){

plot(0,1,type="n",bty='n',xaxt="n",yaxt="n",
      xlim=c(x1,x2),ylim=c(y1,y2),main=" ",xlab=" ",ylab=" ")

yw<-max(Ymax)*0.2

xb<-0
yb<--yw
xt<-ChrL[i]
yt<-0
rect(xb,yb,xt,yt,col="gray40",border=1)


  Hotspot<-DataHotspot[[i]]
  Y1<-(seq(yb*0.4,yb*0.6,length=2))
  color<-c("Navy","red")
  for(k in 1:2){
    Hot.k<-NH:(NH+1)
    X1<-Hotspot[[Hot.k[k]]]
    if(length(X1)!=0){
       X2<-X1+0.5*ScanStep
      Y2<-Y1[k]
      points(X2, rep(Y2, length(X2)),pch=20,col=color[k])
     }
   }


    Na<-paste("Chr",Na.chr[i],sep="")
    text(0,(yt+yb)*0.5,Na,pos=2,col=1,cex=1)

    ChrC<-DataCrop3[,2]
    if(is.na(ChrC[i])!=TRUE){
    points(ChrC[i],yb*1.1,col=1,pch=17)}

    by.x<-ceiling(diff(seq(0,xt,length=20))[1]/5)*5
    axis(1,seq(0,xt,by=by.x),labels =seq(0,xt,by=by.x),cex.axis=0.75,las=0,tck=-0.005,pos=yb,col=1,col.axis="gray40",mgp=c(3,0.005,0))

    MarkerPos<-c(as.matrix(DataQTL[DataQTL[,2]==Na.chr[i],3:4]))

    axis(1,MarkerPos,labels =FALSE,cex.axis=0.6,las=0,tck=-0.02,
             pos=0,col="green",col.axis=1,mgp=c(3,0.05,0))


    ylim<-ceiling(max(Ymax)/5)*5
    by.y<-ceiling(diff(seq(0,ylim,length=20))[1]/5)*5

    axis(2,seq(0,ylim,by=by.y),labels =seq(0,ylim,by=by.y),cex.axis=0.7,las=1,tck=-0.005,
           pos=0,col="gray40",col.axis="gray40",mgp=c(3,0.5,0))

  by.y.qtl<-seq(y2*0.05,y2*0.95,length=NoQTL[i])

  count.t<-0
  for(j in 1:N.trait){
   DataQTL2<-DataQTL[DataQTL[,2]==Na.chr[i]&DataQTL[,1]==Na.trait[j],3:4]

   if(nrow(DataQTL2)>0){
     Index.L<-order(DataQTL2[,1])
     DataQTL3<-DataQTL2[Index.L,]
     No.trait<-nrow(DataQTL3)
     count.t<-count.t+No.trait
     by.y2<-by.y.qtl[(count.t-No.trait+1):count.t]
     L3<-DataQTL3[,1]
     R3<-DataQTL3[,2]
     segments(L3,by.y2,R3,by.y2,col=rainbow(N.trait,start=0.1,end=0.9,v=0.8)[j],lty=j)
     points(L3,by.y2,pch=j,col=rainbow(N.trait,start=0.1,end=0.9,v=0.8)[j],lty=j)
     points(R3,by.y2,pch=j,col=rainbow(N.trait,start=0.1,end=0.9,v=0.8)[j],lty=j)

  }

  }

  DataQTL.chr<-factor(as.vector(DataQTL[DataQTL[,2]==Na.chr[i],1] ),levels=Na.trait)

  DataQTL.chr2<-table(DataQTL.chr)
  DataQTL.chr3<-round(DataQTL.chr2/length(DataQTL.chr),3)*100

  by.y.chr<-seq(yt,y2,length=(N.trait+2))
  y.pt<-by.y.chr[2:(N.trait+1)]
  x.pt<-rep(xt+(0.02*x2),N.trait)

  points(x.pt,y.pt,pch=1:N.trait,
         col=rainbow(N.trait,start=0.1,end=0.9,v=0.8)[1:N.trait])

  Na.x<-paste(DataQTL.chr2," ( ",DataQTL.chr3,"% )",sep="")
  axis(4,y.pt,labels =Na.x,cex.axis=1,las=1,tck=-0.005,
     pos=xt+(0.05*x2),col=par()$bg,col.axis="navy",mgp=c(3,0.005,0))

  axis(4,yt,labels ="Total",cex.axis=1,las=1,tck=-0.005,
     pos=xt+(0.05*x2),col=par()$bg,col.axis="navy",mgp=c(3,0.005,0))

  Na.x.total<-paste(NoQTL[i]," ( ",100,"% )",sep="")
  axis(4,0.5*(yb+yt),labels =Na.x.total,cex.axis=1,las=1,tck=-0.005,
     pos=xt+(0.05*x2),col=par()$bg,col.axis="navy",mgp=c(3,0.005,0))

Left<-DataF[[i]][,2]
Freq<-DataF[[i]][,4]

if((max(Left)+ScanStep)<xt){
x.last<-(max(Left)+ScanStep)
}else{
  x.last<-xt
}
lines(c(Left,x.last),c(Freq,0),type="s",lwd=1,col="gray20")

Na.submain<-paste("Draw all flanking markers of QTLs in the chromosome ",Na.chr[i],sep="")
axis(1,0.5*ChrL[i],labels =Na.submain,cex.axis=1,las=0,tck=-0.005,
     pos=0.5*(yb+y1),col=par()$bg,col.axis=1,mgp=c(3,0.005,0))

Hot.max<-sort(TableH2[,i],decreasing=TRUE)[1]

if(Hot.max>0){
Index.H<-order(TableH2[,i],decreasing=TRUE)[1]
Pos.hot<-DataHotspot[[i]][[Index.H]]

for(nHot in 1:Hot.max){

plot(0,1,type="n",bty='n',xaxt="n",yaxt="n",
      xlim=c(x1,x2),ylim=c(y1,y2),main=" ",xlab=" ",ylab=" ")

yw<-max(Ymax)*0.2

xb<-0
yb<--yw
xt<-ChrL[i]
yt<-0
rect(xb,yb,xt,yt,col="gray40",border=1)

rect(Pos.hot[nHot],yb,Pos.hot[nHot]+ScanStep,yt,col="gray50",border="gray80")
rect(Pos.hot[nHot],yt,Pos.hot[nHot]+ScanStep,y2,col="gray99",border="gray95")

  Hotspot<-DataHotspot[[i]]
  Y1<-(seq(yb*0.4,yb*0.6,length=2))
  color<-c("Navy","red")
  for(k in 1:2){
    Hot.k<-NH:(NH+1)
    X1<-Hotspot[[Hot.k[k]]]
    if(length(X1)!=0){
       X2<-X1+0.5*ScanStep
      Y2<-Y1[k]
      points(X2, rep(Y2, length(X2)),pch=20,col=color[k])
     }
   }

    Na<-paste("Chr",Na.chr[i],sep="")
    text(0,(yt+yb)*0.5,Na,pos=2,col=1,cex=1)
    ChrC<-DataCrop3[,2]
    if(is.na(ChrC[i])!=TRUE){
    points(ChrC[i],yb*1.1,col=1,pch=17)}
    by.x<-ceiling(diff(seq(0,xt,length=20))[1]/5)*5
    axis(1,seq(0,xt,by=by.x),labels =seq(0,xt,by=by.x),cex.axis=0.75,las=0,tck=-0.005,pos=yb,col=1,col.axis="gray40",mgp=c(3,0.005,0))
    segments(0,0,xt,0,col=1)
    ylim<-ceiling(max(Ymax)/5)*5
    by.y<-ceiling(diff(seq(0,ylim,length=20))[1]/5)*5
    axis(2,seq(0,ylim,by=by.y),labels =seq(0,ylim,by=by.y),cex.axis=0.7,las=1,tck=-0.005,
           pos=0,col="gray40",col.axis="gray40",mgp=c(3,0.5,0))

   DataQTL.nHot<-DataQTL[DataQTL[,2]==Na.chr[i],3:4]

   Index.LH<-DataQTL.nHot[,1]<=(Pos.hot[nHot]+ScanStep)
   Index.RH<-DataQTL.nHot[,2]>=Pos.hot[nHot]
   Index.LRH<-Index.LH&Index.RH

  axis(4,yt,labels ="Total",cex.axis=1,las=1,tck=-0.005,
     pos=xt+(0.01*x2),col=par()$bg,col.axis="navy",mgp=c(3,0.005,0))

  Na.x.total<-paste(sum(Index.LRH)," ( ",100,"% )",sep="")
  axis(4,0.5*(yb+yt),labels =Na.x.total,cex.axis=1,las=1,tck=-0.005,
     pos=xt+(0.05*x2),col=par()$bg,col.axis="navy",mgp=c(3,0.005,0))
   DataQTL.nHot2<-factor(as.vector(DataQTL[DataQTL[,2]==Na.chr[i],1][Index.LRH]),levels=Na.trait)

    DataQTL.nHot3<-table(DataQTL.nHot2)
    DataQTL.nHot4<-round(DataQTL.nHot3/sum(DataQTL.nHot3),3)*100
  points(x.pt,y.pt,pch=1:N.trait,
         col=rainbow(N.trait,start=0.1,end=0.9,v=0.8)[1:N.trait])
  Na.x<-paste(  DataQTL.nHot3," ( ", DataQTL.nHot4,"% )",sep="")
  axis(4,y.pt,labels =Na.x,cex.axis=1,las=1,tck=-0.005,
     pos=xt+(0.05*x2),col=par()$bg,col.axis="navy",mgp=c(3,0.005,0))


  by.y.qtl<-seq(y2*0.05,y2*0.95,length=sum(Index.LRH))

  count.t<-0
  for(j in 1:N.trait){

    DataQTL2<-DataQTL[DataQTL[,2]==Na.chr[i]&DataQTL[,1]==Na.trait[j],3:4]
   Index.LH<-DataQTL2$L<=(Pos.hot[nHot]+ScanStep)
   Index.RH<-DataQTL2$R>=Pos.hot[nHot]
   Index.LRH<-Index.LH&Index.RH
   if(sum(Index.LRH)>0){
     DataQTL3<-DataQTL2[Index.LRH,]
     Index.L<-order(DataQTL3[,1])
     DataQTL4<-DataQTL3[Index.L,]
     No.trait<-nrow(DataQTL4)
     count.t<-count.t+No.trait
     by.y2<-by.y.qtl[(count.t-No.trait+1):count.t]
     L3<-DataQTL3[,1]
     R3<-DataQTL3[,2]
     segments(L3,by.y2,R3,by.y2,col=rainbow(N.trait,start=0.1,end=0.9,v=0.8)[j],lty=j)
     points(L3,by.y2,pch=j,col=rainbow(N.trait,start=0.1,end=0.9,v=0.8)[j],lty=j)
     points(R3,by.y2,pch=j,col=rainbow(N.trait,start=0.1,end=0.9,v=0.8)[j],lty=j)

  }


  }


Left<-DataF[[i]][,2]
Freq<-DataF[[i]][,4]

if((max(Left)+ScanStep)<xt){
x.last<-(max(Left)+ScanStep)
}else{
  x.last<-xt
}
lines(c(Left,x.last),c(Freq,0),type="s",lwd=1,col="gray20")



Index.hh<-DataF[[i]]$Left==Pos.hot[nHot]

Na.submain<-paste(DataF[[i]]$BreakPoint[Index.hh],"Expected QTL frequency =",
                  round(DataF[[i]]$ExpectedNumber[Index.hh],2),sep=" ")
axis(1,0.5*ChrL[i],labels =Na.submain,cex.axis=1,las=0,tck=-0.005,
     pos=0.5*(yb+y1),col=par()$bg,col.axis=1,mgp=c(3,0.005,0))


}

}

}


layout(matrix(c(1,2),1,2,byrow=TRUE),widths=c(1,1))

plot(0,1,type="n",bty='n',xaxt="n",yaxt="n",main=" ",xlab=" ",ylab=" ")

Na.legd<-c(paste("Proposed method n= ",NH,sep=""),"Q method" )
legend("left",Na.legd,pch=20,col=color,bty = "n")



plot(0,1,type="n",bty='n',xaxt="n",yaxt="n",main=" ",xlab=" ",ylab=" ")

Na.legd<-legend("left",Na.trait,pch=1:N.trait,lty=1:N.trait,
         col=rainbow(N.trait,start=0.1,end=0.9,v=0.8)[1:N.trait],bty = "n")


names(DataF)<-paste("Chr",Na.chr,sep=" ")
names(TV.P2)<-paste("n= ",1:NH,sep="")
names(TV.Q2)<-c("Q method")
rownames(TableH2)<-c(paste("Proposed method n= ",1:NH,sep=""),"Q method")
colnames(TableH2)<-paste("Chr.",Na.chr,sep="")


print(list(Expected.QTL.frequency.in.every.bin=DataF,Proposed.Method.threshold=TV.P2,
         Q.method.threshold=TV.Q2,The.number.of.Hotspots=TableH2))

}

