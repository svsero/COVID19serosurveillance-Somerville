####code for stochastic simulation to estimate uncertainty 
###in serology-based cumulative incidence by sampling strategy
###to accompany: 
###"GPS-estimated foot traffic data and venue selection for COVID-19 serosurveillance studies"

library(reshape)

#####SET UP#####
# (1) Proportion of population by age group, Somerville
# Age bins match those for MA case reporting: 0-19,20-29,30-39,40-49,50-65,>65
d1 <- c(0.17,0.26,0.23,0.12,0.13,0.09)

# (2) MA age distribution of confirmed COVID cases
ma_cases_p <- c(0.06312,0.14897,0.15306,0.14164,0.22020,0.27302)
sville_cases_p <- c(0,0.18,0.20,0.18,0.17+0.06,0.08+0.06+0.07)

# (3) somerville participant distribution, GPS-estimated visitor distribution
# cases by ward, & testing effort
sville_d <- readRDS('somervilleData.rds')

spop1 <- t(sapply(sville_d$population,function(x){return(x*d1)}))
spop_p1 <- spop1/sum(spop1)

# (4) values to test over for number of tests per group if evenly distributed 
nn <- 5:50
# (5) values to test over for "multiplier" (factor by which true > reported cases)
mm <- 4:40
# (6) all combinations of nn and mm
vmn <- expand.grid(mm,nn)

####FUNCTIONS####
# (1) randomly sample prior infecteds in sample
#binomial distribution with probability i and number of tries n
rb <- Vectorize(function(n,i)(return(rbinom(1,n,i))))

# (2) apply vectorized function over n and i matrices to output weighted prevalence
npos1 <- function(n,i,d_){
  ps <- matrix(ncol=6,nrow=7,rb(n,i))
  
  tp <- matrix(ncol=6,nrow=7,rb(ps,0.90))
  fp <- matrix(ncol=6,nrow=7,rb(n-ps,0.005))
  
  p <- tp+fp
  pr <- p/n
  
  return(c(sum(p),sum(pr*spop_p1)))
}

# (3) wrapper function to obtain uncertainty estimate (width of confidence interval)
# for different example sampling strategies

f2 <- function(m_,n_) {
  #get "ground truth" cumulative incidence by group, i_
  #assume that incidence matches MA case distribution
  #cumulative incidence by ward assuming multiplier m_, distributed same by age distribution in each location
  
  i_ <- t(sapply(sville_d$per*m_,function(x)x*ma_cases_p))
  
  #null: equal sampling across all groups
  n_e <- matrix(ncol=6,nrow=7,n_)
  
  #demo: weighted by population demographics
  n_demo <- round(spop_p1*n_*42)
  re <- sum(n_e)-sum(n_demo);sre <- re/abs(re)
  cfx <- sample(1:length(n_demo),abs(re));n_demo[cfx] <- n_demo[cfx]+sre
  n_demo[which(n_demo==0)] <- 1
  
  #matched to observed distribution of incident cases
  n_g <- round(sville_d$per/sum(sville_d$per)*n_*7)
  re <- sum(n_e[,1])-sum(n_g);cfx <- sample(1:7,1);n_g[cfx] <- n_g[cfx]+re #fix rounding error
  n_g <- cbind(n_g,n_g,n_g,n_g,n_g,n_g)
  n_g[which(n_g==0)] <- 1
  
  #matched to observed distribution of cases and demographics
  n_gd<-round(t(sapply(apply(n_g,1,sum),function(x)x*d1)))
  re <- sum(n_e)-sum(n_gd);sre <- re/abs(re)
  cfx <- sample(1:length(n_gd),abs(re));n_gd[cfx] <- n_gd[cfx]+sre
  n_gd[which(n_gd==0)] <- 1
  
  #matched to observed participants at study site
  n_p <- round(sville_d$par1/sum(sville_d$par1)*n_*7)
  re <- sum(n_e[,1])-sum(n_p);cfx <- sample(1:7,1);n_p[cfx] <- n_p[cfx]+re #fix rounding error
  n_p <- cbind(n_p,n_p,n_p,n_p,n_p,n_p)
  n_p[which(n_p==0)] <- 1
  
  #matched to observed participants at study site and demographics
  n_pd<-round(t(sapply(apply(n_p,1,sum),function(x)x*d1)))
  re <- sum(n_e)-sum(n_pd);sre <- re/abs(re)
  cfx <- sample(1:length(n_pd),abs(re));n_pd[cfx] <- n_pd[cfx]+sre
  n_pd[which(n_pd==0)] <- 1
  
  #matched to estimated visitors to study site
  n_v1 <- round(sville_d$vis1_psv*n_*7)
  re <- sum(n_e[,1])-sum(n_v1);cfx <- sample(1:7,1);n_v1[cfx] <- n_v1[cfx]+re #fix rounding error
  n_v1 <- cbind(n_v1,n_v1,n_v1,n_v1,n_v1,n_v1)
  n_v1[which(n_v1==0)] <- 1
  
  #matched to estimated visitors to study site & age distribution
  n_vd1<-round(t(sapply(apply(n_v1,1,sum),function(x)x*d1)))
  re <- sum(n_e)-sum(n_vd1);sre <- re/abs(re)
  cfx <- sample(1:length(n_vd1),abs(re));n_vd1[cfx] <- n_vd1[cfx]+sre
  n_vd1[which(n_vd1==0)] <- 1
  
  #matched to estimated visitors to alternate site
  n_v2 <- round(sville_d$vis2_psv*n_*7)
  re <- sum(n_e[,1])-sum(n_v2);cfx <- sample(1:7,1);n_v2[cfx] <- n_v2[cfx]+re #fix rounding error
  n_v2 <- cbind(n_v2,n_v2,n_v2,n_v2,n_v2,n_v2)
  n_v2[which(n_v2==0)] <- 1
  
  #matched to estimated visitors to alternate site & age distribution
  n_vd2<-round(t(sapply(apply(n_v2,1,sum),function(x)x*d1)))
  re <- sum(n_e)-sum(n_vd2);sre <- re/abs(re)
  cfx <- sample(1:length(n_vd2),abs(re));n_vd2[cfx] <- n_vd1[cfx]+sre
  n_vd2[which(n_vd2==0)] <- 1
  
  e<-sapply(1:1000,npos1,n=n_e,i=i_)
  g<-sapply(1:1000,npos1,n=n_g,i=i_)
  gd<-sapply(1:1000,npos1,n=n_gd,i=i_)
  p<-sapply(1:1000,npos1,n=n_p,i=i_)
  pd<-sapply(1:1000,npos1,n=n_pd,i=i_)
  demo <- sapply(1:1000,npos1,n=n_demo,i=i_)
  
  v1<-sapply(1:1000,npos1,n=n_v1,i=i_)
  vd1<-sapply(1:1000,npos1,n=n_vd1,i=i_)
  v2<-sapply(1:1000,npos1,n=n_v2,i=i_)
  vd2<-sapply(1:1000,npos1,n=n_vd2,i=i_)
  
  #return(c(mean(e),var(e),mean(g),var(g),mean(gd),var(gd)))
  return(c(quantile(e[2,],0.95,na.rm=T)-quantile(e[2,],0.05,na.rm=T), #e
           quantile(g[2,],0.95,na.rm=T)-quantile(g[2,],0.05,na.rm=T), #g
           quantile(gd[2,],0.95,na.rm=T)-quantile(gd[2,],0.05,na.rm=T), #gd
           quantile(p[2,],0.95,na.rm=T)-quantile(p[2,],0.05,na.rm=T), #p
           quantile(pd[2,],0.95,na.rm=T)-quantile(pd[2,],0.05,na.rm=T), #pd
           quantile(demo[2,],0.95,na.rm=T)-quantile(demo[2,],0.05,na.rm=T), #demo
           quantile(v1[2,],0.95,na.rm=T)-quantile(v1[2,],0.05,na.rm=T), #v1
           quantile(vd1[2,],0.95,na.rm=T)-quantile(vd1[2,],0.05,na.rm=T), #vd1
           quantile(v2[2,],0.95,na.rm=T)-quantile(v2[2,],0.05,na.rm=T), #v2
           quantile(vd2[2,],0.95,na.rm=T)-quantile(vd2[2,],0.05,na.rm=T),
           mean(e[2,],na.rm=T), #e
           mean(g[2,],na.rm=T), #g
           mean(gd[2,],na.rm=T), #gd
           mean(p[2,],na.rm=T), #p
           mean(pd[2,],na.rm=T), #pd
           mean(demo[2,],na.rm=T), #demo
           mean(v1[2,],na.rm=T), #v1
           mean(vd1[2,],na.rm=T), #vd1
           mean(v2[2,],na.rm=T), #v2
           mean(vd2[2,])#vd2
  ))
  
}

vf2 <- Vectorize(f2)

####RUN####
#apply f2 over all combinations of mm and nn
chf2 <- vf2(vmn[,1],vmn[,2])

while (length(which(apply(chf2,2,function(x)any(is.na(x))))) > 0){
for (k in which(apply(chf2,2,function(x)any(is.na(x))))){chf2[,k]<-vf2(vmn[k,1],vmn[k,2])}
}

####REFORMAT###
#melt into mm by nn (rows x columns) matrices for each sampling strategy
#confidence interval width for estimated population-weighted cumulative incidence
w_e <- melt(matrix(ncol=length(nn),nrow=length(mm),t(chf2)[,1]))
w_g <- melt(matrix(ncol=length(nn),nrow=length(mm),t(chf2)[,2]))
w_gd <- melt(matrix(ncol=length(nn),nrow=length(mm),t(chf2)[,3]))
w_p <- melt(matrix(ncol=length(nn),nrow=length(mm),t(chf2)[,4]))
w_pd <- melt(matrix(ncol=length(nn),nrow=length(mm),t(chf2)[,5]))
w_demo <- melt(matrix(ncol=length(nn),nrow=length(mm),t(chf2)[,6]))
w_v1 <- melt(matrix(ncol=length(nn),nrow=length(mm),t(chf2)[,7]))
w_vd1 <- melt(matrix(ncol=length(nn),nrow=length(mm),t(chf2)[,8]))
w_v2 <- melt(matrix(ncol=length(nn),nrow=length(mm),t(chf2)[,9]))
w_vd2 <- melt(matrix(ncol=length(nn),nrow=length(mm),t(chf2)[,10]))

#mean population-weighted cumulative incidence
m_e <- melt(matrix(ncol=length(nn),nrow=length(mm),t(chf2)[,11]))
m_g <- melt(matrix(ncol=length(nn),nrow=length(mm),t(chf2)[,12]))
m_gd <- melt(matrix(ncol=length(nn),nrow=length(mm),t(chf2)[,13]))
m_p <- melt(matrix(ncol=length(nn),nrow=length(mm),t(chf2)[,14]))
m_pd <- melt(matrix(ncol=length(nn),nrow=length(mm),t(chf2)[,15]))
m_demo <- melt(matrix(ncol=length(nn),nrow=length(mm),t(chf2)[,16]))
m_v1 <- melt(matrix(ncol=length(nn),nrow=length(mm),t(chf2)[,17]))
m_vd1 <- melt(matrix(ncol=length(nn),nrow=length(mm),t(chf2)[,18]))
m_v2 <- melt(matrix(ncol=length(nn),nrow=length(mm),t(chf2)[,19]))
m_vd2 <- melt(matrix(ncol=length(nn),nrow=length(mm),t(chf2)[,20]))

####SAVE OUTPUT####
o_n <- c('_e','_g','_gd','_p','_pd','_demo','_v1','_vd1','_v2','_vd2')
for (i in 1:10){
  saveRDS(get(paste0('w',o_n[i])),file=paste0('CIwidth',o_n[i],'.rds'))
  saveRDS(get(paste0('m',o_n[i])),file=paste0('mean',o_n[i],'.rds'))
}

