rm(list=ls())
#First we set our working directory - this needs to be modified to wherever your data are located
setwd("")

###Create my own plots!
library(chron) #chron helps deal with dates and times in R
library(data.table)
library(dplyr)
library(lubridate)
library(ggplot2)
library(ggthemes)
library(ggfortify)
library(gridExtra)
library(RColorBrewer)
library(pracma)
library(patchwork)
library(FSA)
library(car)
library(fmsb)
library(ggridges)
library(viridis)
library(hrbrthemes)
library(plyr)
library(rlang)
library(ggpubr)
library(grid)

options(stringsAsFactors = FALSE)

#MIMS Rep Data N2:Ar and O2:Ar model fig
data<-read.csv(file="mimsrepdatall.newextratimepoint.csv")
data$Season<-factor(data$Season, levels = c("Spring","Summer","Fall"))
data$stream<-factor(data$stream, levels = c("Stream","River"))

p<-ggplot(data, ylim=c(0,50)) + 
  geom_point(data=data, aes(x=new.timepoint,y=N2toArMassRatioRep1,fill=Stream),shape=21,size=2,color="black")+ #,fill="black"
  geom_point(data=data, aes(x=new.timepoint,y=N2toArMassRatioRep2,fill=Stream),shape=21,size=2,color="black")+ #,fill="black"
  geom_point(data=data, aes(x=new.timepoint,y=N2toArMassRatioRep3,fill=Stream),shape=21,size=2,color="black")+ #,fill="black"
  geom_line(data=data, aes(x=new.timepoint,y=dnmodel),linetype="dashed",size=.8)+
  geom_point(data=data, aes(x=new.timepoint,y=NSat/ArSat),shape=21,size=2,color="black",fill="grey")+ 
  geom_line(data=data, aes(x=new.timepoint,y=nodnmodel,color="grey"),linetype="dashed",color="grey",size=.8)+
  #scale_fill_manual(values=c("#e4565b","#e4565b"))+
  theme_classic()+
  scale_fill_manual(values=c("#8FAADC","#586887"))+
  scale_x_continuous(expand = c(0, 0),limits = c(-1, 39),breaks=c(1.5,7.5,13.5,19.5,25.5,31.5,37.5),labels=c("18:00","00:00","06:00","12:00","18:00","00:00","06:00"))+
  #scale_x_discrete(breaks=c("1.5","7.5","13.5","25.5","31.5","37.5"),
  #                 labels=c("18:00","00:00","06:00","12:00","18:00","00:00"))+
  xlab(expression(paste(Time)))+
  #scale_y_continuous(position = "left")+
  #scale_x_continuous(expand=c(0,0))+
  ylab(expression(paste({N}[2], ":Ar")))+
  theme(axis.title.x=element_text(size=14,color="black",vjust=-0.1))+
  theme(axis.title.y=element_text(size=14,color="black",vjust=1.5))+
  theme(axis.text.y=element_text(size=13,color="black"))+
  facet_grid(Season~stream)+
  theme(axis.text.x=element_text(size=13,color="black"))+
  theme(strip.text = element_text(size=15,face = "bold"))+
  theme(legend.position="none")


#NOW TRY TO COLOR THE FACETS (Added the P above)
g <- ggplot_gtable(ggplot_build(p))
stripr <- which(grepl('strip-', g$layout$name))
fills <- c("#8FAADC","#586887","#788E45","#F9BC5D","#8E3600") ##B45B5E
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)

p<-ggplot(data, ylim=c(0,20)) + 
  geom_point(data=data, aes(x=new.timepoint,y=O2toArMassRatioRep1,fill=Stream),shape=21,size=2,color="black")+ 
  geom_point(data=data, aes(x=new.timepoint,y=O2toArMassRatioRep2,fill=Stream),shape=21,size=2,color="black")+ 
  geom_point(data=data, aes(x=new.timepoint,y=O2toArMassRatioRep3,fill=Stream),shape=21,size=2,color="black")+ 
  geom_line(data=data, aes(x=new.timepoint,y=mimsoxy),linetype="dashed",size=0.8)+
  #geom_point(data=data, aes(x=Timepoint,y=OSat/ArSat),shape=21,size=2,color="black",fill="grey")+ 
  #geom_line(data=data, aes(x=Timepoint,y=nodnmodel,color="grey"),linetype="dashed",color="grey")+
  #scale_fill_manual(values=c("#e4565b","#e4565b"))+
  theme_classic()+
  xlab(expression(paste(Time)))+
  scale_x_continuous(expand = c(0, 0),limits = c(-1, 39),breaks=c(1.5,7.5,13.5,19.5,25.5,31.5,37.5),labels=c("18:00","00:00","06:00","12:00","18:00","00:00","06:00"))+
  #scale_y_continuous(position = "left")+
  #scale_x_continuous(expand=c(0,0))+
  ylab(expression(paste({O}[2], ":Ar")))+
  scale_fill_manual(values=c("#8FAADC","#586887"))+
  theme(axis.title.x=element_text(size=14,color="black",vjust=-0.1))+
  theme(axis.title.y=element_text(size=14,color="black",vjust=1.5))+
  theme(axis.text.y=element_text(size=13,color="black"))+
  facet_grid(Season~stream)+
  #theme(
  #  panel.background = element_rect(fill='transparent'), #transparent panel bg
  #  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  #  panel.grid.major = element_blank(), #remove major gridlines
  #  panel.grid.minor = element_blank(), #remove minor gridlines
  #  legend.background = element_rect(fill='transparent'), #transparent legend bg
  #  legend.box.background = element_rect(fill='transparent') #transparent legend panel
  #)+
  theme(axis.text.x=element_text(size=13,color="black"))+
  theme(strip.text=element_text(size=15,face="bold"))+
  theme(legend.position="none")
#theme(
#  strip.background = element_blank(),
# strip.text.x = element_blank(),
# strip.text.y = element_blank(),
# legend.position="none"
#)

g <- ggplot_gtable(ggplot_build(p))
stripr <- which(grepl('strip-', g$layout$name))
fills <- c("#8FAADC","#586887","#788E45","#F9BC5D","#8E3600") ##B45B5E
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)

#Then again for Dissolved O2
data<-read.csv(file="oxydataall.csv")
data$Season<-factor(data$Season, levels = c("Spring","Summer","Fall"))
data$stream<-factor(data$stream, levels = c("Stream","River"))

#By time
p<-ggplot(data, ylim=c(0,50)) + 
  geom_point(data=data, aes(x=new.time,y=oxy,fill=Stream),shape=21,size=2,color="black")+ #,fill="grey"
  geom_line(data=data, aes(x=new.time,y=sondemetabmodel),linetype="dashed",size=0.9)+
  #scale_fill_manual(values=c("#e4565b","#e4565b"))+
  theme_classic()+
  xlab(expression(paste(Time)))+
  scale_fill_manual(values=c("#8FAADC","#586887"))+
  #scale_y_continuous(position = "left")+
  #scale_x_continuous(expand=c(0,0))+
  scale_x_continuous(expand = c(0, 0),limits = c(-1, 39),breaks=c(.9,6.9,12.9,18.9,24.9,30.9,36.4),labels=c("18:00","00:00","06:00","12:00","18:00","00:00","06:00"))+
  ylab(expression(paste(Dissolved~O[2]~(g~m^{"-3"}))))+
  theme(axis.title.x=element_text(size=14,color="black",vjust=-0.1))+
  theme(axis.title.y=element_text(size=14,color="black",vjust=1.5))+
  theme(axis.text.y=element_text(size=13,color="black"))+
  facet_grid(Season~stream)+
  #theme(
  #  panel.background = element_rect(fill='transparent'), #transparent panel bg
  #  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  #  panel.grid.major = element_blank(), #remove major gridlines
  #  panel.grid.minor = element_blank(), #remove minor gridlines
  #  legend.background = element_rect(fill='transparent'), #transparent legend bg
  #  legend.box.background = element_rect(fill='transparent') #transparent legend panel
  #)+
  theme(axis.text.x=element_text(size=13,color="black"))+
  theme(strip.text=element_text(size=15,face="bold"))+
  theme(legend.position="none")
#theme(
#  strip.background = element_blank(),
# strip.text.x = element_blank(),
# strip.text.y = element_blank(),
# legend.position="none"
#)

#NOW TRY TO COLOR THE FACETS (Added the P above)
g <- ggplot_gtable(ggplot_build(p))
stripr <- which(grepl('strip-', g$layout$name))
fills <- c("#8FAADC","#586887","#788E45","#F9BC5D","#8E3600") ##B45B5E
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)

#Then the metabolism Plot
dt<-read.csv(file="dielmasterdatLONG.csv")
dt$Season<-factor(dt$Season, levels = c("Spring","Summer", "Fall"))
dt$stream<-factor(dt$stream, levels = c("Stream","River"))
dt.sha<-subset(dt, stream=="Stream")
dt.tip<-subset(dt, stream=="River")
dt$Method<-factor(dt$Method, levels = c("Sonde","MIMS"))
# Example bar plot
cor.test(dt.tip$abser, dt.tip$denit.n, method = "pearson")

# Example bar plot
p<-ggplot(dt, aes(x = Season, y = gpp, fill = Season,group=Method)) +
  geom_bar(
    stat = "identity", color="black",
    position = position_dodge(width = .9)
  ) +
  geom_errorbar(data=dt, aes(ymin=low.gpp, ymax=high.gpp,group=Method), width=.2,
                position=position_dodge(width = .9))+
  scale_fill_manual(values = c("#788E45","#F9BC5D","#8E3600")) +
  facet_wrap(~ stream) +  # Optional: facet by watershed
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Season") +
  ylab("GPP (g O₂ m⁻² d⁻¹)") +
  scale_y_continuous(expand=c(0,0), limits=c(0,12))+
  theme(
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.position = "none",
    legend.title = element_text(size = 13,color="black"),
    legend.text = element_text(size = 12,color="black")
  )
p
g <- ggplot_gtable(ggplot_build(p))
stripr <- which(grepl('strip-', g$layout$name))
fills <- c("#8FAADC","#586887") ##B45B5E #,"#788E45","#F9BC5D","#8E3600"
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)

er<-ggplot(dt, aes(x = Season, y = abs(er), fill = Season,group=Method)) +
  geom_bar(
    stat = "identity", color="black",
    position = position_dodge(width = .9)
  ) +
  geom_errorbar(data=dt, aes(ymin=low.er, ymax=high.er,group=Method), width=.2,
                position=position_dodge(width = .9))+
  scale_fill_manual(values = c("#788E45","#F9BC5D","#8E3600")) +
  facet_wrap(~ stream) +  # Optional: facet by watershed
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Season") +
  ylab(expression(paste("|"~ER~"|"~(g~O["2"]~m^{"-2"}~d^{"-1"}))))+
  scale_y_continuous(expand=c(0,0), limits=c(0,35))+
  theme(
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.position = "none",
    legend.title = element_text(size = 13,color="black"),
    legend.text = element_text(size = 12,color="black")
  )

j<-p/er

#Then for denitrification
dt<-read.csv(file="dielmasterdat.csv")
#dt<-subset(data,use=="y")
dt.sdw<-subset(dt,Stream=="SHA")
dt.tip<-subset(dt,Stream=="TIP")
sdw.f<-subset(dt.sdw,Season=="Fall")
sdw.s<-subset(dt.sdw,Season=="Summer")
sdw.sp<-subset(dt.sdw,Season=="Spring")
tip.f<-subset(dt.tip,Season=="Fall")
tip.s<-subset(dt.tip,Season=="Summer")
tip.sp<-subset(dt.tip,Season=="Spring")
dt$Stream<-factor(dt$Stream, levels = c("SHA","TIP"))
dt$type<-factor(dt$type, levels = c("Stream","River"))
dt$Season <- as.factor(dt$Season)
dt$Season<-factor(dt$Season, levels = c("Spring","Summer", "Fall"))
dt.sdw$Season<-factor(dt.sdw$Season, levels = c("Spring","Summer", "Fall"))
dt.tip$Season<-factor(dt.tip$Season, levels = c("Spring","Summer", "Fall"))

facet_labels <- c("SHA" = "Stream", "TIP" = "River")

p1<-ggplot(dt, aes( y=denit.n, x=Season, fill=Season))+ 
  geom_bar(position="dodge", stat="identity", colour="black", width=0.8)+
  geom_errorbar(data=sdw.sp, aes(ymin=127.7, ymax=168.8), width=.2,
                position=position_dodge(0.05))+
  geom_errorbar(data=sdw.s, aes(ymin=69.3, ymax=81.6), width=.2,
                position=position_dodge(0.05))+
  geom_errorbar(data=sdw.f, aes(ymin=52.2, ymax=79), width=.2,
                position=position_dodge(0.05))+
  geom_errorbar(data=tip.sp, aes(ymin=2.7, ymax=8.8), width=.2,
                position=position_dodge(0.05))+
  geom_errorbar(data=tip.s, aes(ymin=0.07, ymax=2.2), width=.2,
                position=position_dodge(0.05))+
  geom_errorbar(data=tip.f, aes(ymin=13.1, ymax=24.9), width=.2,
                position=position_dodge(0.05))+
  theme_classic()+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(name="Season", values=c("#788E45","#F9BC5D","#8E3600"))+ ##F9F97C
  xlab(expression(Season))+ 
  ylab(expression(paste(Denitrification~Rate~(mg~N~m^{"-2"}~hr^{"-1"}))))+
  facet_grid(.~Stream, labeller=as_labeller(facet_labels))+
  theme(panel.spacing.y  = unit(.4, "cm"))+
  theme_classic()+
  theme(axis.text.x=element_text(color="black"))+
  theme(axis.text.y=element_text(color="black"))+
  #theme(text = element_text(size=12))+
  theme(axis.title.x=element_text(size=14,color="black",vjust=-0.1))+
  theme(axis.title.y=element_text(size=14,color="black",vjust=1.5))+
  theme(axis.text.y=element_text(size=13,color="black"))+
  theme(axis.text.x=element_text(size=13,color="black"))+
  scale_y_continuous(expand = c(0, 0),limits=c(0,175))+
  #scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  theme(legend.position="none")+
  theme(strip.text=element_text(size=15,face="bold"))
p
g <- ggplot_gtable(ggplot_build(p))
stripr <- which(grepl('strip-', g$layout$name))
fills <- c("#8FAADC","#586887") ##B45B5E #,"#788E45","#F9BC5D","#8E3600"
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)

#Denit Scaled to per km
p2<-ggplot(dt, aes( y=denit.n.g.km, x=Season, fill=Season))+ 
  geom_bar(position="dodge", stat="identity", colour="black", width=0.8)+
  geom_errorbar(data=sdw.sp, aes(ymin=280.94, ymax=371.36), width=.2,
                position=position_dodge(0.05))+
  geom_errorbar(data=sdw.s, aes(ymin=152.46, ymax=179.52), width=.2,
                position=position_dodge(0.05))+
  geom_errorbar(data=sdw.f, aes(ymin=114.84, ymax=173.8), width=.2,
                position=position_dodge(0.05))+
  geom_errorbar(data=tip.sp, aes(ymin=162, ymax=528), width=.2,
                position=position_dodge(0.05))+
  geom_errorbar(data=tip.s, aes(ymin=4.2, ymax=132), width=.2,
                position=position_dodge(0.05))+
  geom_errorbar(data=tip.f, aes(ymin=786, ymax=1494), width=.2,
                position=position_dodge(0.05))+
  theme_classic()+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(name="Season", values=c("#788E45","#F9BC5D","#8E3600"))+ ##F9F97C
  xlab(expression(Season))+ 
  ylab(expression(paste(Linear~Denitrification~Rate~(g~N~km^{"-1"}~hr^{"-1"}))))+
  facet_grid(.~Stream, labeller=as_labeller(facet_labels))+
  theme(panel.spacing.y  = unit(.4, "cm"))+
  theme_classic()+
  theme(axis.text.x=element_text(color="black"))+
  theme(axis.text.y=element_text(color="black"))+
  #theme(text = element_text(size=12))+
  theme(axis.title.x=element_text(size=14,color="black",vjust=-0.1))+
  theme(axis.title.y=element_text(size=14,color="black",vjust=1.5))+
  theme(axis.text.y=element_text(size=13,color="black"))+
  theme(axis.text.x=element_text(size=13,color="black"))+
  scale_y_continuous(expand = c(0, 0),limits=c(0,1750))+
  #scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  theme(legend.position="none")+
  theme(strip.text=element_text(size=15,face="bold"))
g <- ggplot_gtable(ggplot_build(p))
stripr <- which(grepl('strip-', g$layout$name))
fills <- c("#8FAADC","#586887") ##B45B5E #,"#788E45","#F9BC5D","#8E3600"
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)

g1 <- ggplot_gtable(ggplot_build(p1))  # Assuming 'p1' is your first ggplot
stripr1 <- which(grepl('strip-', g1$layout$name))
fills1 <- c("#8FAADC","#586887") 
k <- 1
for (i in stripr1) {
  j <- which(grepl('rect', g1$grobs[[i]]$grobs[[1]]$childrenOrder))
  g1$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills1[k]
  k <- k + 1
}

# Create the second plot and modify facet colors
g2 <- ggplot_gtable(ggplot_build(p2))  # Assuming 'p2' is your second ggplot
stripr2 <- which(grepl('strip-', g2$layout$name))
fills2 <- c("#8FAADC","#586887")  
k <- 1
for (i in stripr2) {
  j <- which(grepl('rect', g2$grobs[[i]]$grobs[[1]]$childrenOrder))
  g2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills2[k]
  k <- k + 1
}
grid.arrange(g1, g2, ncol = 2)

#Now compare the ER and GPP from the models
GPPooo<-ggplot(dt, aes(x=sonde.gpp, y=mims.gpp))+ #, ylim=c(0,80)
  geom_point(size=5,color="black", aes(fill=Season, shape=Stream))+
  scale_fill_manual(values=c("#788E45","#F9BC5D","#8E3600"))+
  scale_shape_manual(values = c(21, 22)) +
  theme_classic()+
  geom_errorbar(data=sdw.sp, aes(ymin=6.9, ymax=7.9), width=.2,
                position=position_dodge(0.05))+
  geom_errorbar(data=sdw.s, aes(ymin=4.8, ymax=5.4), width=.2,
                position=position_dodge(0.05))+
  geom_errorbar(data=sdw.f, aes(ymin=1.1, ymax=1.4), width=.2,
                position=position_dodge(0.05))+
  geom_errorbar(data=tip.sp, aes(ymin=1.7, ymax=2.4), width=.2,
                position=position_dodge(0.05))+
  geom_errorbar(data=tip.s, aes(ymin=3.3, ymax=3.9), width=.2,
                position=position_dodge(0.05))+
  geom_errorbar(data=tip.f, aes(ymin=.2, ymax=1.4), width=.2,
                position=position_dodge(0.05))+
  geom_errorbarh(data=sdw.sp, aes(xmin=8.9, xmax=10.0), width=.2,
                 position=position_dodge(0.05))+
  geom_errorbarh(data=sdw.s, aes(xmin=5.5, xmax=6.2), width=.2,
                 position=position_dodge(0.05))+
  geom_errorbarh(data=sdw.f, aes(xmin=.95, xmax=1.1), width=.2,
                 position=position_dodge(0.05))+
  geom_errorbarh(data=tip.sp, aes(xmin=2.4, xmax=2.8), width=.2,
                 position=position_dodge(0.05))+
  geom_errorbarh(data=tip.s, aes(xmin=5.0, xmax=5.3), width=.2,
                 position=position_dodge(0.05))+
  geom_errorbarh(data=tip.f, aes(xmin=.7, xmax=1.3), width=.2,
                 position=position_dodge(0.05))+
  #stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),label.x = 0, label.y = 155, size=4, color="black") +
  xlab(expression(paste(Sonde~GPP~(g~O["2"]~m^{"-2"}~d^{"-1"}))))+
  ylab(expression(paste(MIMS~GPP~(g~O["2"]~m^{"-2"}~d^{"-1"}))))+
  geom_abline(sdw.f, yintercept=28,slope=1,color="red",size=1)+
  #geom_smooth(method="lm",se=FALSE, color="black",size=0.8)+
  theme(axis.title.x=element_text(size=13,color="black",vjust=-0.1))+
  theme(axis.title.y=element_text(size=13,color="black",vjust=1.5))+
  theme(axis.text.y=element_text(size=12,color="black"))+
  theme(axis.text.x=element_text(size=12,color="black"))+
  scale_y_continuous(expand=c(0,0),limits=c(0,12))+
  scale_x_continuous(expand=c(0,0),limits=c(0,12))+
  theme(legend.position="none")
ERooo<-ggplot(dt, aes(x=abs(sonde.er), y=abs(mims.er)))+ #, ylim=c(0,80)
  geom_point(size=5,color="black", aes(fill=Season,shape=Stream))+
  scale_fill_manual(values=c("#788E45","#F9BC5D","#8E3600"))+
  scale_shape_manual(values = c(21, 22)) +
  theme_classic()+
  geom_errorbar(data=sdw.sp, aes(ymin=14.8, ymax=12.8), width=.2,
                position=position_dodge(0.05))+
  geom_errorbar(data=sdw.s, aes(ymin=14.4, ymax=13.1), width=.2,
                position=position_dodge(0.05))+
  geom_errorbar(data=sdw.f, aes(ymin=7.7, ymax=6.1), width=.2,
                position=position_dodge(0.05))+
  geom_errorbar(data=tip.sp, aes(ymin=2.1, ymax=1.5), width=.2,
                position=position_dodge(0.05))+
  geom_errorbar(data=tip.s, aes(ymin=9.9, ymax=8), width=.2,
                position=position_dodge(0.05))+
  geom_errorbar(data=tip.f, aes(ymin=20, ymax=12.4), width=.2,
                position=position_dodge(0.05))+
  geom_errorbarh(data=sdw.sp, aes(xmin=17.4, xmax=15.6), width=.2,
                 position=position_dodge(0.05))+
  geom_errorbarh(data=sdw.s, aes(xmin=14.5, xmax=13.0), width=.2,
                 position=position_dodge(0.05))+
  geom_errorbarh(data=sdw.f, aes(xmin=5.7, xmax=4.9), width=.2,
                 position=position_dodge(0.05))+
  geom_errorbarh(data=tip.sp, aes(xmin=2.9, xmax=2.6), width=.2,
                 position=position_dodge(0.05))+
  geom_errorbarh(data=tip.s, aes(xmin=12.7, xmax=11.75), width=.2,
                 position=position_dodge(0.05))+
  geom_errorbarh(data=tip.f, aes(xmin=34.0, xmax=15.0), width=.2,
                 position=position_dodge(0.05))+
  #stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),label.x = 0, label.y = 155, size=4, color="black") +
  xlab(expression(paste(Sonde~"|"~ER~"|"~(g~O["2"]~m^{"-2"}~d^{"-1"}))))+
  ylab(expression(paste(MIMS~"|"~ER~"|"~(g~O["2"]~m^{"-2"}~d^{"-1"}))))+
  geom_abline(sdw.f, yintercept=28,slope=1,color="red",size=1)+
  #geom_smooth(method="lm",se=FALSE, color="black",size=0.8)+
  theme(axis.title.x=element_text(size=13,color="black",vjust=-0.1))+
  theme(axis.title.y=element_text(size=13,color="black",vjust=1.5))+
  theme(axis.text.y=element_text(size=12,color="black"))+
  theme(axis.text.x=element_text(size=12,color="black"))+
  scale_y_continuous(expand=c(0,0),limits=c(0,35))+
  scale_x_continuous(expand=c(0,0),limits=c(0,35))+
  theme(legend.position="none")

GPPooo/ERooo


#Metanalysis
data<-read.csv(file="streammetadenit.csv")
data$study<-factor(data$study, levels = c("Mulholland et al. 2009", "Findlay et al. 2011","Roley et al. 2012","Mahl et al. 2015","Reisinger et al. 2016","Speir et al. 2020","This Study"))

#Metanalysis Figure
data<-read.csv(file="denitmeta.5.23.23.csv")
data$Paper<-factor(data$Paper, levels = c("Laursen and Seitzinger 2002","Laursen and Seitzinger 2004","McCutchan and Lewis 2008","Smith et al. 2008","Mulholland et al. 2009", "Ritz et al. 2018","Nifong et al. 2020","Reisinger et al. 2016","This Study-SHA","This Study-TIP"))
dt.study<-subset(data,Paper=="This Study-SHA"|Paper=="This Study-TIP") #,"Reisinger et al. 2016"

ggplot(data, aes(x = NO3.Conc.mgL, y = denit.rate.mgm2hr)) +
  # Plot points for "Other Studies"
  geom_point(
    data = subset(data), #, Group == "Other Studies"
    aes(fill = num, shape=Stream),
    size = 4, color="darkgrey"
  ) +
  # Plot points for "This Study"
  geom_point(
    data = subset(dt.study),
    aes(fill = Season, shape = Stream),
    size = 5,
    color = "black"
  ) +
  # Add regression line for "This Study"
  geom_smooth(
    data = dt.study,
    method = "lm", color="black",
    se = FALSE,
    size = 0.8
  ) +
  theme_classic() +
  xlab(expression(paste(NO["3"]^{"-"}-N~(mg~N~L^{-1}))))+
  ylab(expression(paste(Denitrification~Rate~(mg~N~m^{"-2"}~hr^{"-1"})))) +
  scale_x_continuous( limits = c(0, 15)) +
  scale_y_continuous( limits = c(0, 250)) +
  scale_fill_manual(
    name = "Season",
    values = c("Spring" = "#788E45", "Summer" = "#F9BC5D", "Fall" = "#8E3600", "Other" = "#DBDBDB"))+
  scale_shape_manual(name = "Stream", values = c(22, 21)) +
  theme(
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 13, color = "black"),
    axis.title.y = element_text(size = 13, color = "black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13),
    legend.position = "none")
