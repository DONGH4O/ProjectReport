# ProjectReport
R code for the Imperial summer project
#Summer Project


############
#Simulation#
############


#Set up
d<-1024
eps<-c(0.001,0.005,0.01,0.1,0.25)
sig<-c(1,1,1)
amp<-c(0.5,1,3,5,7)
repl<-50


##make data simulations
data.sim<-function(Eps,Amp,d){
ind<-sample(1:d,d*eps[Eps],replace=FALSE)
thet<-(c(1:d) %in% ind)*amp[Amp]
rnorm(d,mean=0,sd=sig[Sig])+thet
}
#trace plot
d.100<-data.sim(3,3,100)
d.10000<-data.sim(3,3,10000)
eps.1<-data.sim(1,3,1000)
eps.5<-data.sim(5,3,1000)
amp.1<-data.sim(3,1,1000)
amp.5<-data.sim(3,5,1000)

par(mfrow=c(1,2))
plot(d.100,type = "l",ylab = "data value",main = "n=100")
plot(d.10000,type = "l",ylab = "data value",main = "n=10000")
plot(eps.1,type = "l",ylab = "data value",main = "epsilon=0.001")
plot(eps.5,type = "l",ylab = "data value",main = "epsilion=0.01")
plot(amp.1,type = "l",ylab = "data value",main = "amplitude=0.5")
plot(amp.5,type = "l",ylab = "data value",main = "amplitude=3")


#SURE
SURE<-function(Sam,Sig,Lam){
  #d*sig[Sig]^2+2*sig[Sig]^2*Lam*sum((Sam>Lam)-(Sam<(-Lam))+Sam*(abs(Sam)<=Lam))-sum(abs(Sam)<=Lam)
  d*sig[Sig]^2-2*sig[Sig]^2*sum((abs(Sam)<=Lam))+sum(Sam^2*(Sam^2<=Lam^2)+Lam^2*(Sam^2>Lam^2))
}

#Soft Threshold
mu.est<-function(Mu,Lam){
  (Mu<(-Lam))*(Mu+Lam)+(Mu>Lam)*(Mu-Lam)
}
#Lam
Lam.ind<-function(Mu.sort,d,Sig){
  sure<-sapply(1:d,function(i){
    lambda<-rep(abs(Mu.sort[i]),d)
    SURE(Mu.sort,Sig,lambda)
  })
  which.min(sure)
}


#estimation
estimat<-function(Eps,Amp,Sig,d){
  ind<-sample(1:d,d*eps[Eps],replace=FALSE)
  thet<-(c(1:d) %in% ind)*amp[Amp]
  mu<-rnorm(d,mean=0,sd=sig[Sig])+thet
  #SURE
  lamind<-Lam.ind(mu,d,Sig)
  lam<-rep(abs(mu[lamind]),d)
  su.est<-mu.est(mu,lam)
  #D-J
  lam.F<-rep(sqrt(2*log(d)),d)
  dj.est<-mu.est(mu,lam.F)
  list("mu"=mu,"su.est"=su.est,"dj.est"=dj.est,"thet"=thet)
}

samp.est<-estimat(4,3,1,50)
plot(unlist(samp.est[1]),type = "l",ylab = "values")
lines(unlist(samp.est[2]),col="red")
lines(unlist(samp.est[3]),col="blue")
lines(unlist(samp.est[4]),col="purple")
legend("topleft",legend = c("Data","SUREshrink","sqrt(2log n)","theta"),col = c("black","red","blue","purple"),lty = c(1,1,1,1))








#MSE comparison
sam.mak<-function(Eps,Sig,Amp){
  ind<-sample(1:d,d*eps[Eps],replace=FALSE)
  thet<-(c(1:d) %in% ind)*amp[Amp]
  mu<-rnorm(d,mean=0,sd=sig[Sig])+thet
  #SURE method
  mu.sort<-sort(mu)
  lamind<-Lam.ind(mu.sort,d,Sig)
  lam<-rep(abs(mu.sort[lamind]),d)
  Mse.S<-mean((mu.est(mu,lam)-thet)^2)
  #New method
  lam.F<-rep(sqrt(2*log(d)),d)
  Mse.F<-mean((mu.est(mu,lam.F)-thet)^2)
  list(Mse.s=Mse.S,Mse.f=Mse.F)
}


#First plot when level of sparsity changes
mse.15.3<-sapply(1:length(eps),function(i) replicate(repl,sam.mak(i,1,3)))



sf.15.1<-matrix(0,nrow = repl,ncol = 10)
for(i in 1:repl){
  for(j in 1:length(eps)){
    sf.15.1[i,2*j]<-sqrt(unlist(mse.15.1[2*i,j]))
    sf.15.1[i,2*j-1]<-sqrt(unlist(mse.15.1[i*2-1,j]))
  }
}
colnames(sf.15.1)<-c("S eps0.001","F eps0.001","S eps0.005","F eps0.005","S eps0.01","F eps0.01","S eps0.1","F eps0.1","S eps0.25","F eps0.25")
sf.eps<-data.frame(sf.15.1)
boxplot(sf.eps,col = c(0,1,0,1,0,1,0,1),main = "Contrast of MSE of two thresholdings with level of sparsity as variable",ylab="Square root of MSE",xlab = "Sparcity")
legend("topleft",legend = c("SUREshrink","Global rule"),fill = c(0,1))



#second plot when level of amplitude changes
mse.3.15<-sapply(1:length(amp),function(i) replicate(repl,sam.mak(3,1,i)))



sf.3.15<-matrix(0,nrow = repl,ncol = 10)
for(i in 1:repl){
  for(j in 1:length(amp)){
    sf.3.15[i,2*j]<-sqrt(unlist(mse.3.15[2*i,j]))
    sf.3.15[i,2*j-1]<-sqrt(unlist(mse.3.15[i*2-1,j]))
  }
}
colnames(sf.3.15)<-c("S amp0.5","F amp0.5","S amp1","F amp1","S amp3","F amp3","S amp5","F amp5","S amp7","F amp7")
sf.amp<-data.frame(sf.3.15)
boxplot(sf.amp,col = c(1,0,1,0,1,0,1,0,1,0),main = "Contrast of MSE of two thresholdings with amplitude as variable",xlab = "Amplitude")
legend("topleft",legend = c("SUREshrink","Global rule"),fill = c(1,0))


#################
##  FDR check  ##
#################
fdr<-function(thet,zeroind){
  if(length(zeroind)!=0){
    fal.dis<-which(thet[zeroind]==0)
    length(fal.dis)/length(zeroind)
  }else{
    0
  }
}

fdr.make<-function(Amp,Eps,Sig,d){
  #generating data
  ind<-sample(1:d,d*eps[Eps],replace=FALSE)
  thet<-(c(1:d) %in% ind)*amp[Amp]
  mu<-rnorm(d,mean=0,sd=sig[Sig])+thet
  #SURE method
  lamind<-Lam.ind(mu,d,Sig)
  lam<-rep(abs(mu[lamind]),d)
  su.est<-mu.est(mu,lam)
  su.zin<-which(su.est!=0)
  su.fd<-fdr(thet,su.zin)
  #D-J method
  lam.F<-rep(sqrt(2*log(d)),d)
  dj.est<-mu.est(mu,lam.F)
  dj.zin<-which(dj.est!=0)
  dj.fd<-fdr(thet,dj.zin)
  list("su.fd"=su.fd,"dj.fd"=dj.fd)
}

fdr.comp<-function(Amp,Eps,Sig,d,repli){
  sure.fd<-0
  djme.fd<-0
  fdr.result<-replicate(repli,fdr.make(Amp,Eps,Sig,d))
  list("su.fdr"=sum(as.numeric(fdr.result[1,]))/repli,"dj.fdr"=sum(as.numeric(fdr.result[2,]))/repli)
}



fdr.15.1<-sapply(1:length(eps),function(i) replicate(repl,fdr.make(3,i,1,1024)))

fmt.15.1<-matrix(0,nrow = repl,ncol = 10)
for(i in 1:repl){
  for(j in 1:length(eps)){
    fmt.15.1[i,2*j]<-unlist(fdr.15.1[2*i,j])
    fmt.15.1[i,2*j-1]<-unlist(fdr.15.1[i*2-1,j])
  }
}
colnames(fmt.15.1)<-c("S eps0.001","F eps0.001","S eps0.005","F eps0.005","S eps0.01","F eps0.01","S eps0.1","F eps0.1","S eps0.25","F eps0.25")
fmt.eps<-data.frame(fmt.15.1)
boxplot(fmt.eps,col = c(0,1,0,1,0,1,0,1),main = "Contrast of FDR of two thresholdings with level of sparsity as variable",ylab="False discovery rate",xlab = "Sparcity")
legend("topright",legend = c("SUREshrink","Global rule"),fill = c(0,1))



fdr.3.15<-sapply(1:length(amp),function(i) replicate(repl,fdr.make(i,3,1,1024)))

fmt.3.15<-matrix(0,nrow = repl,ncol = 10)
for(i in 1:repl){
  for(j in 1:length(amp)){
    fmt.3.15[i,2*j]<-unlist(fdr.3.15[2*i,j])
    fmt.3.15[i,2*j-1]<-unlist(fdr.3.15[i*2-1,j])
  }
}
colnames(fmt.3.15)<-c("S amp0.5","F amp0.5","S amp1","F amp1","S amp3","F amp3","S amp5","F amp5","S amp7","F amp7")
fmt.amp<-data.frame(fmt.3.15)
boxplot(fmt.amp,col = c(0,1,0,1,0,1,0,1,0,1),main = "Contrast of FDR of two thresholdings with amplitude as variable",ylab="False dicovery rate",xlab = "Amplitude")
legend(x=0.5,y=0.5,legend = c("SUREshrink","Global Rule"),fill = c(0,1))

#############################
##                         ##
##  MSE extreme situation  ##
##                         ##
#############################

mse.1.1<-replicate(repl,sam.mak(1,1,1))
mse.5.5<-replicate(repl,sam.mak(5,1,5))

aaa<-cbind(unlist(mse.1.1[1,]),unlist(mse.1.1[2,]))
colnames(aaa)<-c("SUREshrink","Globalrule")
bbb<-cbind(unlist(mse.5.5[1,]),unlist(mse.5.5[2,]))
colnames(bbb)<-c("SUREshrink","Globalrule")

par(mfrow=c(2,1))
boxplot(aaa,main = "Contrast of MSE with low level of sparsity and low amplitude",ylab="Mean Squared Error")
boxplot(bbb,main = "Contrast of MSE with high level of sparsity and high amplitude",ylab="Mean Squared Error")




mse.comb<-t(rbind(mse.1.1,mse.5.5))
colnames(mse.comb)<-c("S.small","F.small","S.large","F.large")
sf.1.1<-as.data.frame(mse.comb)
boxplot(mse.comb,col = c(0,1,0,1),main = "Contrast of FDR of two thresholdings with amplitude as variable",ylab="False dicovery rate",xlab = "Amplitude")


write.table(mse.comb,file = "E:/Documents/Imperial College London/Project/data/sf11.txt",quote = TRUE)


##############################
##Check for adaptive methods##
##############################


###Hybrid method
sdsq<-function(Sam,d){
  sum(Sam^2-1)/d
}
gam<-function(d){
  log(d,base = 2)^1.5/sqrt(d)
}
#adaptive method
adrep<-function(Theta,Sig,d){
  mu<-rnorm(d,mean=0,sd=sig[Sig])+Theta
  mu.sort<-sort(mu)
  lamind<-Lam.ind(mu.sort,d,Sig)
  lam<-rep(abs(mu.sort[lamind]),d)
  Mse.S<-mean((mu.est(mu,lam)-Theta)^2)
  #New method
  lam.F<-rep(sqrt(2*log(d)),d)
  Mse.F<-mean((mu.est(mu,lam.F)-Theta)^2)
  max(c(0,1)*c(Mse.S<=Mse.F,Mse.F<Mse.S))
}


##Test four new method
examin<-function(Amp,Eps,Sig,d,ad.repl){
  #generating data
  ind<-sample(1:d,d*eps[Eps],replace=FALSE)
  thet<-(c(1:d) %in% ind)*amp[Amp]
  mu<-rnorm(d,mean=0,sd=sig[Sig])+thet
  #SURE method
  lamind<-Lam.ind(mu,d,Sig)
  lam<-rep(abs(mu[lamind]),d)
  Mse.S<-mean((mu.est(mu,lam)-thet)^2)
  #New method
  lam.F<-rep(sqrt(2*log(d)),d)
  Mse.F<-mean((mu.est(mu,lam.F)-thet)^2)
  su.m<-Mse.S
  dj.m<-Mse.F
  #hybrid method
  benc<-sdsq(mu,d)
  thres<-gam(d)
  if(benc<=thres) {
    hy.m<-Mse.F
  }else{
    hy.m<-Mse.S
  }
  #adaptive method
  lam.ad<-rep(sqrt(2*log(d)),d)
  mu.ad<-mu.est(mu,lam.ad)
  ad.resu<-replicate(ad.repl,adrep(mu.ad,Sig,d))
  ad.cho<-round(mean(ad.resu),0)
  if(ad.cho<0.5){
    ad.m<-Mse.S
  }else{
    ad.m<-Mse.F
  }
  list("hy.m"=hy.m,"ad.m"=ad.m,"su.m"=su.m,"dj.m"=dj.m)
}


comp<-function(Amp,Eps,Sig,d,ad.repl,repli){
comp.result<-replicate(repli,examin(Amp,Eps,Sig,d,ad.repl))
list("hy.m"=sum(as.numeric(comp.result[1,])),"ad.m"=sum(as.numeric(comp.result[2,])),"su.m"=sum(as.numeric(comp.result[3,])),"dj.m"=sum(as.numeric(comp.result[4,])))
}



compare.result.13<-comp(1,3,1,1024,500,1000)
compare.result.23<-comp(2,3,1,1024,500,1000)
compare.result.33<-comp(3,3,1,1024,500,1000)
compare.result.43<-comp(4,3,1,1024,500,1000)
compare.result.53<-comp(5,3,1,1024,500,1000)

#compare.result.31<-comp(3,1,1,1024,500,1000)
#compare.result.32<-comp(3,2,1,1024,500,1000)
#compare.result.33<-comp(3,3,1,1024,500,1000)#already have
compare.result.34<-comp(3,4,1,1024,500,1000)
compare.result.35<-comp(3,5,1,1024,500,1000)

compare.result.31<-read.table("E:/Documents/Imperial College London/Project/data/compare.result.31.txt")
compare.result.32<-read.table("E:/Documents/Imperial College London/Project/data/compare.result.32.txt")

compare.result.31<-as.list(compare.result.31)
compare.result.32<-as.list(compare.result.32)



##########################
##                      ##
##  parallel computing  ##
##                      ##
##########################
library(parallel)
cores<-detectCores()-1
cl<-makeCluster(cores)
clusterExport(cl,c('fdr.find','fdr','fdr.re','eps','amp','sig','Lam.ind','mu.est','sdsq','gam','adrep',"SURE",'sam.mak','d','fdr.find.SUREbase','fdr.re.SUREbase','comp','examin'))
compare.result.13.3<-parSapply(cl,1:3,function(i) comp(i,3,1,1024,500,1000))
stopCluster(cl)




spar<-cbind(unlist(compare.result.31),unlist(compare.result.32),unlist(compare.result.33),unlist(compare.result.34),unlist(compare.result.35))
ampl<-cbind(unlist(compare.result.13),unlist(compare.result.23),unlist(compare.result.33),unlist(compare.result.43),unlist(compare.result.53))

write.table(spar,file = "E:/Documents/Imperial College London/Project/data/spar_newone.txt",quote=TRUE)
write.table(ampl,file = "E:/Documents/Imperial College London/Project/data/ampl_newone.txt",quote=TRUE)




#Graph
par(mfrow=c(2,2))

plot(eps,spar[2,],ylim=c(0,2000) ,main = "Total mse amoung four methods for sparsity changing (Global rule base)",xlab = "Proportion of sparsity",ylab = "Mse",type = "l")
lines(eps,spar[4,],col="blue")
lines(eps,spar[6,],col="green")
lines(eps,spar[8,],col="red")
legend("topleft",legend = c("Hybrid method","Adaptive method","Sure method","Globalrule method"),col = c("black","blue","green","red"),lty = c(1,1,1,1))

plot(amp,ampl[2,],ylim = c(0,150),main = "Total mse amoung four methods for amplitude changing (Global rule base)",xlab = "Amplitude",ylab = "Mse",type = "l")
lines(amp,ampl[4,],col="blue")
lines(amp,ampl[6,],col="green")
lines(amp,ampl[8,],col="red")
legend("topleft",legend = c("Hybrid method","Adaptive method","Sure method","Global rule method"),col = c("black","blue","green","red"),lty = c(1,1,1,1))




#########################
##False Discovery Rate###
#########################

fdr<-function(thet,zeroind){
  if(length(zeroind)!=0){
  fal.dis<-which(thet[zeroind]==0)
  length(fal.dis)/length(zeroind)
  }else{
  0
}
}
fdr.find<-function(Amp,Eps,Sig,d,ad.repl){
  #generating data
  ind<-sample(1:d,d*eps[Eps],replace=FALSE)
  thet<-(c(1:d) %in% ind)*amp[Amp]
  mu<-rnorm(d,mean=0,sd=sig[Sig])+thet
  #SURE method
  lamind<-Lam.ind(mu,d,Sig)
  lam<-rep(abs(mu[lamind]),d)
  su.est<-mu.est(mu,lam)
  su.zin<-which(su.est!=0)
  su.fd<-fdr(thet,su.zin)
  #D-J method
  lam.F<-rep(sqrt(2*log(d)),d)
  dj.est<-mu.est(mu,lam.F)
  dj.zin<-which(dj.est!=0)
  dj.fd<-fdr(thet,dj.zin)
  fdr.cho<-c(su.fd,dj.fd)
  #hybrid method
  benc<-sdsq(mu,d)
  thres<-gam(d)
  hy.ind<-(benc<=thres)+1
  hy.fd<-fdr.cho[hy.ind]
  #adaptive method
  lam.ad<-rep(sqrt(2*log(d)),d)
  mu.ad<-mu.est(mu,lam.ad)
  ad.resu<-replicate(ad.repl,adrep(mu.ad,Sig,d))
  ad.cho<-round(mean(ad.resu),0)+1
  ad.fd<-fdr.cho[ad.cho]
  list("su.fd"=su.fd,"dj.fd"=dj.fd,"hy.fd"=hy.fd,"ad.fd"=ad.fd)
}



fdr.re<-function(Amp,Eps,Sig,d,ad.repl,repli){
  fdr.result<-replicate(repli,fdr.find(Amp,Eps,Sig,d,ad.repl))
  list("su.fdr"=sum(as.numeric(fdr.result[1,]))/repli,"dj.fdr"=sum(as.numeric(fdr.result[2,]))/repli,"hy.fdr"=sum(as.numeric(fdr.result[3,]))/repli,"ad.fdr"=sum(as.numeric(fdr.result[4,]))/repli)
}

fdr.result.31<-fdr.re(3,1,1,1024,500,1000)
fdr.result.32<-fdr.re(3,2,1,1024,500,1000)
fdr.result.33<-fdr.re(3,3,1,1024,500,1000)
fdr.result.34<-fdr.re(3,4,1,1024,500,1000)
fdr.result.35<-fdr.re(3,5,1,1024,500,1000)

fdr.result.13<-fdr.re(1,3,1,1024,500,1000)
fdr.result.23<-fdr.re(2,3,1,1024,500,1000)
fdr.result.43<-fdr.re(4,3,1,1024,500,1000)
fdr.result.53<-fdr.re(5,3,1,1024,500,1000)#left uncalculated

##########################
##                      ##
##  parallel computing  ##
##                      ##
##########################
library(parallel)
cores<-detectCores()-1
cl<-makeCluster(cores)
clusterExport(cl,c('fdr.find','fdr','fdr.re','eps','amp','sig','Lam.ind','mu.est','sdsq','gam','adrep',"SURE",'sam.mak','d'))
fdr.result.14.3<-parSapply(cl,c(1,2,4),function(i) fdr.re(i,3,1,1024,500,1000))#already calculated.
stopCluster(cl)

write.table(fdr.result.14.3,file = "E:/Documents/Imperial College London/Project/data/fdr_result_14_3.txt",quote = TRUE)



spar.fdr<-cbind(unlist(fdr.result.31),unlist(fdr.result.32),unlist(fdr.result.33),unlist(fdr.result.34),unlist(fdr.result.35))
ampl.fdr<-cbind(unlist(fdr.result.13),unlist(fdr.result.23),unlist(fdr.result.33),unlist(fdr.result.43),unlist(fdr.result.53))

write.table(spar.fdr,file = "E:/Documents/Imperial College London/Project/data/spar_fdr.txt",quote=TRUE)
write.table(ampl.fdr,file = "E:/Documents/Imperial College London/Project/data/ampl_fdr.txt",quote=TRUE)

spar.fdr11<-read.table("E:/Documents/Imperial College London/Project/data/spar_fdr.txt",header = TRUE)
ampl.fdr11<-read.table("E:/Documents/Imperial College London/Project/data/ampl_fdr.txt",header = TRUE)

spar.fdr1<-as.list(spar.fdr11)
ampl.fdr1<-as.list(ampl.fdr11)

plot(eps,spar.fdr[1,],ylim = c(0,1),main = "False discovery rate amoung four methods for sparsity changing (Global rule base)",xlab = "Proportion of sparsity",ylab = "FDR",type = "l")
lines(eps,spar.fdr[2,],col="blue")
lines(eps,spar.fdr[3,],col="green")
lines(eps,spar.fdr[4,],col="red")
legend("topright",legend = c("Sure method","Globalrule method","Hybrid method","Adaptive method"),col = c("black","blue","green","red"),lty = c(1,1,1,1))


plot(amp,ampl.fdr[1,],ylim = c(0,1),main = "False discovery rate amoung four methods for amplitude changing(Global rule base)",xlab = "amplitude",ylab = "FDR",type = "l")
lines(amp,ampl.fdr[2,],col="blue")
lines(amp,ampl.fdr[3,],col="green")
lines(amp,ampl.fdr[4,],col="red")
legend(x=0.24,y=0.6,legend = c("Sure method","Globalrule method","Hybrid method","Adaptive method"),col = c("black","blue","green","red"),lty = c(1,1,1,1))


###################################
##new method with SURE as initial##
###################################

###########
##       ##
##  MSE  ##
##       ##
###########


###Hybrid method
sdsq<-function(Sam,d){
  sum(Sam^2-1)/d
}
gam<-function(d){
  log(d,base = 2)^1.5/sqrt(d)
}
#adaptive method
adrep<-function(Theta,Sig,d){
  mu<-rnorm(d,mean=0,sd=sig[Sig])+Theta
  mu.sort<-sort(mu)
  lamind<-Lam.ind(mu.sort,d,Sig)
  lam<-rep(abs(mu.sort[lamind]),d)
  Mse.S<-mean((mu.est(mu,lam)-Theta)^2)
  #New method
  lam.F<-rep(sqrt(2*log(d)),d)
  Mse.F<-mean((mu.est(mu,lam.F)-Theta)^2)
  max(c(0,1)*c(Mse.S<=Mse.F,Mse.F<Mse.S))
}


##Test four new method
examin.su<-function(Amp,Eps,Sig,d,ad.repl){
  #generating data
  ind<-sample(1:d,d*eps[Eps],replace=FALSE)
  thet<-(c(1:d) %in% ind)*amp[Amp]
  mu<-rnorm(d,mean=0,sd=sig[Sig])+thet
  #SURE method
  lamind<-Lam.ind(mu,d,Sig)
  lam<-rep(abs(mu[lamind]),d)
  Mse.S<-mean((mu.est(mu,lam)-thet)^2)
  #New method
  lam.F<-rep(sqrt(2*log(d)),d)
  Mse.F<-mean((mu.est(mu,lam.F)-thet)^2)
  
  su.m<-Mse.S
  dj.m<-Mse.F
  #hybrid method
  benc<-sdsq(mu,d)
  thres<-gam(d)
  if(benc<=thres) {
    hy.m<-Mse.F
  }else{
    hy.m<-Mse.S
  }
  #adaptive method
  mu.su<-mu.est(mu,lam)
  ad.resu<-replicate(ad.repl,adrep(mu.su,Sig,d))
  ad.cho<-round(mean(ad.resu),0)
  if(ad.cho<0.5){
    ad.m<-Mse.S
  }else{
    ad.m<-Mse.F
  }
  list("hy.m"=hy.m,"ad.m"=ad.m,"su.m"=su.m,"dj.m"=dj.m)
}


comp.su<-function(Amp,Eps,Sig,d,ad.repl,repli){
  comp.result<-replicate(repli,examin.su(Amp,Eps,Sig,d,ad.repl))
  list("hy.m"=sum(as.numeric(comp.result[1,])),"ad.m"=sum(as.numeric(comp.result[2,])),"su.m"=sum(as.numeric(comp.result[3,])),"dj.m"=sum(as.numeric(comp.result[4,])))
}


compare.result.31.su<-comp(3,1,1,1024,500,1000)
compare.result.32.su<-comp(3,2,1,1024,500,1000)
compare.result.33.su<-comp(3,3,1,1024,500,1000)
compare.result.34.su<-comp(3,4,1,1024,500,1000)
compare.result.35.su<-comp(3,5,1,1024,500,1000)


compare.result.13.su<-comp.su(1,3,1,1024,500,1000)#calculated
compare.result.23.su<-comp.su(2,3,1,1024,500,1000)#
#compare.result.33.su<-comp.su(3,3,1,1024,500,1000)
compare.result.43.su<-comp.su(4,3,1,1024,500,1000)#
compare.result.53.su<-comp.su(5,3,1,1024,500,1000)#


##########################
##                      ##
##  parallel computing  ##
##                      ##
##########################
library(parallel)
cores<-detectCores(logical = TRUE)-2
cl<-makeCluster(cores)
clusterExport(cl,c('fdr.find','fdr','fdr.re','eps','amp','sig','Lam.ind','mu.est','sdsq','gam','adrep',"SURE",'sam.mak','d','comp.su','examin.su'))
compare.result.3.45.su<-parSapply(cl,c(4,5),function(i) comp.su(3,i,1,1024,500,1000))
stopCluster(cl)




compare.result.15.3.su<-parSapply(cl,c(1,2,4,5),function(i) comp.su(i,3,1,1024,500,1000))
write.table(compare.result.15.3.su,file = "E:/Documents/Imperial College London/Project/data/compare_result_15_3.txt",quote = TRUE)
compare.result.3.13.su<-parSapply(cl,c(1,2,3),function(i) comp.su(3,i,1,1024,500,1000))

###########
##       ##
##  fdr  ##
##       ##
###########


######new method with SURE as initial
fdr.find.SUREbase<-function(Amp,Eps,Sig,d,ad.repl){
  #generating data
  ind<-sample(1:d,d*eps[Eps],replace=FALSE)
  thet<-(c(1:d) %in% ind)*amp[Amp]
  mu<-rnorm(d,mean=0,sd=sig[Sig])+thet
  #SURE method
  lamind<-Lam.ind(mu,d,Sig)
  lam<-rep(abs(mu[lamind]),d)
  su.est<-mu.est(mu,lam)
  su.zin<-which(su.est!=0)
  su.fd<-fdr(thet,su.zin)
  #D-J method
  lam.F<-rep(sqrt(2*log(d)),d)
  dj.est<-mu.est(mu,lam.F)
  dj.zin<-which(dj.est!=0)
  dj.fd<-fdr(thet,dj.zin)
  fdr.cho<-c(su.fd,dj.fd)
  #hybrid method
  benc<-sdsq(mu,d)
  thres<-gam(d)
  hy.ind<-(benc<=thres)+1
  hy.fd<-fdr.cho[hy.ind]
  #adaptive method
  mu.su<-su.est
  ad.resu<-replicate(ad.repl,adrep(mu.su,Sig,d))
  ad.cho<-round(mean(ad.resu),0)+1
  ad.fd<-fdr.cho[ad.cho]
  list("su.fd"=su.fd,"dj.fd"=dj.fd,"hy.fd"=hy.fd,"ad.fd"=ad.fd)
}



fdr.re.SUREbase<-function(Amp,Eps,Sig,d,ad.repl,repli){
  fdr.result<-replicate(repli,fdr.find.SUREbase(Amp,Eps,Sig,d,ad.repl))
  list("su.fdr"=sum(as.numeric(fdr.result[1,]))/repli,"dj.fdr"=sum(as.numeric(fdr.result[2,]))/repli,"hy.fdr"=sum(as.numeric(fdr.result[3,]))/repli,"ad.fdr"=sum(as.numeric(fdr.result[4,]))/repli)
}

fdr.result.su.31<-fdr.re.SUREbase(3,1,1,1024,500,1000)
fdr.result.su.32<-fdr.re.SUREbase(3,2,1,1024,500,1000)
fdr.result.su.33<-fdr.re.SUREbase(3,3,1,1024,500,1000)
fdr.result.su.34<-fdr.re.SUREbase(3,4,1,1024,500,1000)
fdr.result.su.35<-fdr.re.SUREbase(3,5,1,1024,500,1000)


fdr.result.su.13<-fdr.re.SUREbase(1,3,1,1024,500,1000)
fdr.result.su.23<-fdr.re.SUREbase(2,3,1,1024,500,1000)
#fdr.result.su.33<-fdr.re.SUREbase(3,3,1,1024,500,1000)
fdr.result.su.43<-fdr.re.SUREbase(4,3,1,1024,500,1000)
fdr.result.su.53<-fdr.re.SUREbase(5,3,1,1024,500,1000)#only one left



##########################
##                      ##
##  parallel computing  ##
##                      ##
##########################
library(parallel)
cores<-detectCores()-2
cl<-makeCluster(cores)
clusterExport(cl,c('fdr.find','fdr','fdr.re','eps','amp','sig','Lam.ind','mu.est','sdsq','gam','adrep',"SURE",'sam.mak','d','fdr.find.SUREbase','fdr.re.SUREbase'))
fdr.result.su.3.45<-parSapply(cl,c(4,5),function(i) fdr.re.SUREbase(3,i,1,1024,500,1000))
stopCluster(cl)


fdr.result.su.3.13<-parSapply(cl,1:3,function(i) fdr.re.SUREbase(3,i,1,1024,500,1000))
fdr.result.su.14.3<-parSapply(cl,c(1,2,4),function(i) fdr.re.SUREbase(i,3,1,1024,500,1000))
write.table(fdr.result.su.3.13,file = "E:/Documents/Imperial College London/Project/data/fdr_result_su_3_13.txt",quote = TRUE)
write.table(fdr.result.su.3.45,file = "E:/Documents/Imperial College London/Project/data/fdr_result_su_3_45.txt",quote = TRUE)
write.table(fdr.result.su.14.3,file = "E:/Documents/Imperial College London/Project/data/fdr_result_su_14_3.txt",quote = TRUE)

fdr.result.su.3.45<-read.table("E:/Documents/Imperial College London/Project/data/fdr_result_su_3_45.txt",header = TRUE)
fdr.result.su.14.3<-read.table("E:/Documents/Imperial College London/Project/data/fdr_result_su_14_3.txt",header = TRUE)

spar.fdr.su<-cbind(unlist(fdr.result.su.31),unlist(fdr.result.su.32),unlist(fdr.result.su.33),unlist(fdr.result.su.34),unlist(fdr.result.su.35))
ampl.fdr.su<-cbind(unlist(fdr.result.su.13),unlist(fdr.result.su.23),unlist(fdr.result.su.33),unlist(fdr.result.su.43),unlist(fdr.result.su.53))

plot(eps,spar.fdr.su[1,],ylim = c(0,1),main = "False discovery rate amoung four methods for sparsity changing(SURE base)",xlab = "Proportion of sparsity",ylab = "FDR",type = "l")
lines(eps,spar.fdr.su[2,],col="blue")
lines(eps,spar.fdr.su[3,],col="green")
lines(eps,spar.fdr.su[4,],col="red")
legend("topright",legend = c("Sure method","Globalrule method","Hybrid method","Adaptive method"),col = c("black","blue","green","red"),lty = c(1,1,1,1))


plot(amp,ampl.fdr.su[1,],ylim = c(0,1),main = "False discovery rate amoung four methods for amplitude changing(SURE base)",xlab = "amplitude",ylab = "FDR",type = "l")
lines(amp,ampl.fdr.su[2,],col="blue")
lines(amp,ampl.fdr.su[3,],col="green")
lines(amp,ampl.fdr.su[4,],col="red")
legend(x=4.17,y=0.6,legend = c("Sure method","Globalrule method","Hybrid method","Adaptive method"),col = c("black","blue","green","red"),lty = c(1,1,1,1))



plot(eps[1:2],spar.fdr.su[1,1:2],type = "l",ylim = c(0,1))

plot(eps[2],spar.fdr.su[1,2])














################
##            ##
##  New Aera  ##
##            ##
################

fdr.result.53<-read.table("E:/Documents/Imperial College London/Project/data/fdr_result_53.txt",header=TRUE)
spar.fdr<-read.table("E:/Documents/Imperial College London/Project/data/spar_fdr_new.txt",header=TRUE)
fdr.result.su.53<-read.table("E:/Documents/Imperial College London/Project/data/fdr_result_su_53.txt",header=TRUE)

spar.fdr<-as.data.frame(spar.fdr)

ampl.fdr<-cbind(fdr.result.14.3[,1:2],spar.fdr[,3],unlist(fdr.result.14.3[,3]),unlist(fdr.result.53))

spar.fdr.su<-cbind(unlist(fdr.result.su.3.13[,1]),unlist(fdr.result.su.3.13[,2]),unlist(fdr.result.su.3.13[,3]),unlist(fdr.result.su.3.45[,1]),unlist(fdr.result.su.3.45[,2]))

ampl.fdr.su<-cbind(fdr.result.su.14.3[,1:2],unlist(spar.fdr.su[,3]),unlist(fdr.result.su.14.3[,3]),unlist(fdr.result.su.53))



ampl.mse.su<-read.table("E:/Documents/Imperial College London/Project/data/ampl_mse_su.txt",header=TRUE)
spar.mse.su<-read.table("E:/Documents/Imperial College London/Project/data/spar_mse_su.txt",header=TRUE)

plot(eps,spar.mse.su[1,],ylim=c(0,2000),main = "Total mse amoung four methods for sparsity changing (SUREshrink base)",xlab = "Proportion of sparsity",ylab = "Mse",type = "l")
lines(eps,spar.mse.su[2,],col="blue")
lines(eps,spar.mse.su[3,],col="green")
lines(eps,spar.mse.su[4,],col="red")
legend("topleft",legend = c("Hybrid method","Adaptive method","SUREshrink method","Globalrule method"),col = c("black","blue","green","red"),lty = c(1,1,1,1))

plot(amp,ampl.mse.su[1,],ylim = c(0,150),main = "Total mse amoung four methods for amplitude changing (SUREshrink base)",xlab = "Amplitude",ylab = "Mse",type = "l")
lines(amp,ampl.mse.su[2,],col="blue")
lines(amp,ampl.mse.su[3,],col="green")
lines(amp,ampl.mse.su[4,],col="red")
legend("topleft",legend = c("Hybrid method","Adaptive method","SUREshrink method","Globalrule method"),col = c("black","blue","green","red"),lty = c(1,1,1,1))






#############choosing rate
cr<-function(Amp,Eps,Sig,d,ad.repl){
  ind<-sample(1:d,d*eps[Eps],replace=FALSE)
  thet<-(c(1:d) %in% ind)*amp[Amp]
  mu<-rnorm(d,mean=0,sd=sig[Sig])+thet
  #hybrid method
  benc<-sdsq(mu,d)
  thres<-gam(d)
  hy<-as.numeric(benc<=thres)
  #adaptive method Global base
  lam.ad<-rep(sqrt(2*log(d)),d)
  mu.ad<-mu.est(mu,lam.ad)
  ad.resu<-replicate(ad.repl,adrep(mu.ad,Sig,d))
  ad<-round(mean(ad.resu),0)
  #adaptive method SUREbase
  lamind.su<-Lam.ind(mu,d,Sig)
  lam.su<-rep(abs(mu[lamind.su]),d)
  mu.su<-mu.est(mu,lam.su)
  ad.resu.su<-replicate(ad.repl,adrep(mu.su,Sig,d))
  ad.su<-round(mean(ad.resu.su),0)
  #### 1=DJ, 0=SURE
  list('hy'=hy,'ad'=ad,'ad.su'=ad.su)
}


cr.result<-function(Amp,Eps,Sig,d,ad.repl,repli){
  result<-replicate(repli,cr(Amp,Eps,Sig,d,ad.repl))
  list("hy.cr"=mean(unlist(result[1,])),"ad.cr"=mean(unlist(result[2,])),"ad.su.cr"=mean(unlist(result[3,])))
}


library(parallel)
cores<-detectCores()-2
cl<-makeCluster(cores)
clusterExport(cl,c('fdr.find','fdr','fdr.re','eps','amp','sig','Lam.ind','mu.est','sdsq','gam','adrep',"SURE",'sam.mak','d','fdr.find.SUREbase','fdr.re.SUREbase','cr','cr.result'))
cr.result.12.3<-parSapply(cl,c(1,2),function(i) cr.result(i,3,1,1024,500,1000))
stopCluster(cl)





