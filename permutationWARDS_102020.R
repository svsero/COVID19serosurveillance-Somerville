####permutation testing to test significance of correlations between
####participant and visitor counts by geographic unit (electoral ward)
####to accompany:
####"GPS-estimated foot traffic data and venue selection for COVID-19 serosurveillance studies"
library(ggplot2)

####load data####
sville_d <- readRDS('somervilleData.rds')

####functions####
#function for correlation with single random permutation
###of geographic location (ward) for each variable

prm <- function(a,b,d){
  ap <- a[sample(1:length(a))]
  bp <- a[sample(1:length(b))]
  return(cor(ap,bp))
}

incidence <- sville_d$per/sum(sville_d$per)

corners <- data.frame(x=c(0,0.51,0.51),y=c(0,0,0.51))

####plots####
# (1) directly-recruited study participants vs GPS-estimated visitors (study site)
r_vS_pD <- round(cor(sv_joined$vis1_psv,sv_joined$par1_psv),2)
p_vS_pD <- round(length(which(r_vS_pD < sapply(1:10000,prm,a=sv_joined$vis1_psv,b=sv_joined$par1_psv)))/10000,2)

visS_vs_parD <- ggplot(corners) + geom_polygon(data=corners,aes(x=x,y=y),fill='grey95') + geom_point(data=sv_joined,aes(y=vis1_psv,x=par1_psv),color='blue') +
  geom_text(data=sv_joined,aes(y=vis1_psv,x=par1_psv),label=1:7,
            color='blue',vjust = "top", hjust = "right", nudge_y=0.01, nudge_x=-0.01) +
  theme_classic() + xlab(expression(p['direct'])) + xlim(0,0.51) + ylim(0,0.51) + ylab(expression(v['site'])) +
  annotate("text",  x=-Inf, y = Inf, label = paste0('Pearson\'s r = ',r_vS_pD), vjust=1,hjust=-0.1,size=4) +
  annotate("text",  x=-Inf, y = Inf, label = paste0('p = ',p_vS_pD), vjust=3,hjust=-0.2,size=4)

ggsave('directparticipantsVSvisitors.pdf',plot=visS_vs_parD,device = 'pdf',width=4,height=4)

# (2) directly-recruited study participants vs incidence of laboratory-confirmed infections
r_i_pD <- round(cor(incidence,sv_joined$par1_psv),2)
p_i_pD <- round(length(which(r_i_pD < sapply(1:10000,prm,a=incidence,b=sv_joined$par1_psv)))/10000,2)

i_vs_parD <- ggplot(corners) + geom_polygon(data=corners,aes(x=x,y=y),fill='grey95') + geom_point(data=sv_joined,aes(y=per/sum(per),x=par1_psv),color='red') +
  geom_text(data=sv_joined,aes(y=per/sum(per),x=par1_psv),label=1:7,
            color='red',vjust = "top", hjust = "right", nudge_y=0.01, nudge_x=-0.01) +
  theme_classic() + xlab(expression(p['direct'])) + xlim(0,0.51) + ylim(0,0.51) + ylab(expression(theta['ward,confirmed'])) +
  annotate("text",  x=-Inf, y = Inf, label = paste0('Pearson\'s r = ',r_i_pD), vjust=1,hjust=-0.1,size=4) +
  annotate("text",  x=-Inf, y = Inf, label = paste0('p = ',p_i_pD), vjust=3,hjust=-0.2,size=4)

ggsave('directparticipantsVSincidence.pdf',plot=i_vs_parD,device = 'pdf',width=4,height=4)

# (3) all study participants (directly-recruited and other) vs incidence of laboratory-confirmed infections
r_i_pA <- round(cor(incidence,sv_joined$par2_psv),2)
p_i_pA <- round(length(which(r_i_pA < sapply(1:10000,prm,a=incidence,b=sv_joined$par2_psv)))/10000,2)

i_vs_parA <- ggplot(corners) + geom_polygon(data=corners,aes(x=x,y=y),fill='grey95') + geom_point(data=sv_joined,aes(y=per/sum(per),x=par2_psv),color='red') +
  geom_text(data=sv_joined,aes(y=per/sum(per),x=par2_psv),label=1:7,
            color='red',vjust = "top", hjust = "right", nudge_y=0.01, nudge_x=-0.01) +
  theme_classic() + xlab(expression(p['all'])) + xlim(0,0.51) + ylim(0,0.51) + ylab(expression(theta['ward,confirmed'])) +
  annotate("text",  x=-Inf, y = Inf, label = paste0('Pearson\'s r = ',r_i_pA), vjust=1,hjust=-0.1,size=4) +
  annotate("text",  x=-Inf, y = Inf, label = paste0('p = ',p_i_pA), vjust=3,hjust=-0.2,size=4)

ggsave('allparticipantsVSincidence.pdf',plot=i_vs_parA,device = 'pdf',width=4,height=4)

# (4) GPS-estimated visitors (study site) vs incidence of laboratory-confirmed infections
r_i_pv1 <- round(cor(incidence,sv_joined$vis1_psv),2)
p_i_pv1 <- round(length(which(r_i_pv1 < sapply(1:10000,prm,a=incidence,b=sv_joined$vis1_psv)))/10000,2)

i_vs_visS <- ggplot(corners) + geom_polygon(data=corners,aes(x=x,y=y),fill='grey95') + geom_point(data=sv_joined,aes(y=per/sum(per),x=vis1_psv),color='red') +
  geom_text(data=sv_joined,aes(y=per/sum(per),x=vis1_psv),label=1:7,
            color='red',vjust = "top", hjust = "right", nudge_y=0.01, nudge_x=-0.01) +
  theme_classic() + xlab(expression(v['site'])) + xlim(0,0.51) + ylim(0,0.51) + ylab(expression(theta['ward,confirmed'])) +
  annotate("text",  x=-Inf, y = Inf, label = paste0('Pearson\'s r = ',r_i_pv1), vjust=1,hjust=-0.1,size=4) +
  annotate("text",  x=-Inf, y = Inf, label = paste0('p = ',p_i_pv1), vjust=3,hjust=-0.2,size=4)

ggsave('visitorsVSincidence.pdf',plot=i_vs_visS,device = 'pdf',width=4,height=4)

# (5) GPS-estimated visitors (alternate site) vs incidence of laboratory-confirmed infections
r_i_pv2 <- round(cor(incidence,sv_joined$vis2_psv),2)
p_i_pv2 <- round(length(which(r_i_pv2 < sapply(1:10000,prm,a=incidence,b=sv_joined$vis2_psv)))/10000,2)

i_vs_visA <- ggplot(corners) + geom_polygon(data=corners,aes(x=x,y=y),fill='grey95') + geom_point(data=sv_joined,aes(y=per/sum(per),x=vis2_psv),color='red') +
  geom_text(data=sv_joined,aes(y=per/sum(per),x=vis2_psv),label=1:7,
            color='red',vjust = "top", hjust = "right", nudge_y=0.01, nudge_x=-0.01) +
  theme_classic() + xlab(expression(v['alternate'])) + xlim(0,0.51) + ylim(0,0.51) + ylab(expression(theta['ward,confirmed'])) +
  annotate("text",  x=-Inf, y = Inf, label = paste0('Pearson\'s r = ',r_i_pv2), vjust=1,hjust=-0.1,size=4) +
  annotate("text",  x=-Inf, y = Inf, label = paste0('p = ',p_i_pv2), vjust=3,hjust=-0.2,size=4)

ggsave('visitorsalternateVSincidence.pdf',plot=i_vs_visA,device = 'pdf',width=4,height=4)

