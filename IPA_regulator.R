setwd("C:/Users/yudonglin/Desktop/硕士论文写作/硕士论文/figure-整理/IPA分析")
library(ggplot2)

ipa<-read.csv("Upstream Regulators_E0.csv")
#ipa_2<-ipa[,c(1,4,5,7)]
ipa_2<-ipa[,c(3,6,7,11)]
head(ipa_2)
ipa_2[is.na(ipa_2)]<-0
colnames(ipa_2)
ipa_2$log_p<--log10(ipa_2$p.value.of.overlap)
ipa_2$Upstream.Regulator<-as.character(ipa_2$Upstream.Regulator)
ipa_2$NAME<-ifelse(abs(ipa_2$Activation.z.score)>=2 & ipa_2$p.value.of.overlap< 0.05,
                   ipa_2$Upstream.Regulator,"" )
ipa_2$NAME<-ifelse(abs(ipa_2$Activation.z.score)>=6 & ipa_2$p.value.of.overlap< 0.05,
                   ipa_2$Upstream.Regulator,"" )



#ipa_2$Predicted.Activation.State
ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(size=3.5)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("black","red","blue"))+
  geom_point()+geom_text(nudge_y = 1,size=3)+ylab("-log10(p-value)")

ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(size=.1)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("red","blue","black"))+
  geom_point()+ylab("-log10(p-value)")+geom_text(nudge_y = 1,size=3)





ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(aes(color = Predicted.Activation.State, size = log_p), show.legend = F)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("red","blue","black"))+
  geom_point()+ylab("-log10(p-value)")+geom_text(nudge_y = 1,size=3)

library(tidyverse)
library(ggrepel)
library(ggridges)
library(lubridate)
library(modelr)
library(patchwork)

ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(size=.1)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("red","blue","black"))+
  geom_point()+ylab("-log10(p-value)")+geom_text_repel(nudge_y = 1,size=3)
ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(aes(color = Predicted.Activation.State, size = log_p), show.legend = F)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("red","blue","black"))+
  geom_point()+ylab("-log10(p-value)")+geom_text_repel(nudge_y = 1,size=3)


pdf("E0_upstream.pdf",width = 6,height = 6)
ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(aes(color = Predicted.Activation.State, size = log_p), show.legend = F)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("red","blue","black"))+
  geom_point()+ylab("-log10(p-value)")+geom_text_repel(
  data = ipa_2, aes(label = NAME),
  size = 3, fontface = 3,
  nudge_x = 2, nudge_y = 2) +ggtitle("IPA results of E0")+
  theme(plot.title = element_text(hjust = 0.5))  #也就加上这一行
dev.off()


setwd("C:/Users/yudonglin/Desktop/硕士论文写作/硕士论文/figure-整理/IPA分析")
library(ggplot2)

ipa<-read.csv("Upstream Regulators_E1.csv")
#ipa_2<-ipa[,c(1,4,5,7)]
ipa_2<-ipa[,c(3,6,7,11)]
head(ipa_2)
ipa_2[is.na(ipa_2)]<-0
colnames(ipa_2)
ipa_2$log_p<--log10(ipa_2$p.value.of.overlap)
ipa_2$Upstream.Regulator<-as.character(ipa_2$Upstream.Regulator)
ipa_2$NAME<-ifelse(abs(ipa_2$Activation.z.score)>=2 & ipa_2$p.value.of.overlap< 0.05,
                   ipa_2$Upstream.Regulator,"" )
ipa_2$NAME<-ifelse(abs(ipa_2$Activation.z.score)>=6 & ipa_2$p.value.of.overlap< 0.05,
                   ipa_2$Upstream.Regulator,"" )



#ipa_2$Predicted.Activation.State
ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(size=3.5)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("black","red","blue"))+
  geom_point()+geom_text(nudge_y = 1,size=3)+ylab("-log10(p-value)")

ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(size=.1)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("red","blue","black"))+
  geom_point()+ylab("-log10(p-value)")+geom_text(nudge_y = 1,size=3)





ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(aes(color = Predicted.Activation.State, size = log_p), show.legend = F)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("red","blue","black"))+
  geom_point()+ylab("-log10(p-value)")+geom_text(nudge_y = 1,size=3)

library(tidyverse)
library(ggrepel)
library(ggridges)
library(lubridate)
library(modelr)
library(patchwork)

ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(size=.1)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("red","blue","black"))+
  geom_point()+ylab("-log10(p-value)")+geom_text_repel(nudge_y = 1,size=3)
ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(aes(color = Predicted.Activation.State, size = log_p), show.legend = F)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("red","blue","black"))+
  geom_point()+ylab("-log10(p-value)")+geom_text_repel(nudge_y = 1,size=3)


pdf("E1_upstream.pdf",width = 6,height = 6)
ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(aes(color = Predicted.Activation.State, size = log_p), show.legend = F)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("red","blue","black"))+
  geom_point()+ylab("-log10(p-value)")+geom_text_repel(
    data = ipa_2, aes(label = NAME),
    size = 3, fontface = 3,
    nudge_x = 2, nudge_y = 2) +ggtitle("IPA results of E1")+
  theme(plot.title = element_text(hjust = 0.5))  #也就加上这一行
dev.off()


setwd("C:/Users/yudonglin/Desktop/硕士论文写作/硕士论文/figure-整理/IPA分析")
library(ggplot2)

ipa<-read.csv("Upstream Regulators_E2.csv")
#ipa_2<-ipa[,c(1,4,5,7)]
ipa_2<-ipa[,c(3,6,7,11)]
head(ipa_2)
ipa_2[is.na(ipa_2)]<-0
colnames(ipa_2)
ipa_2$log_p<--log10(ipa_2$p.value.of.overlap)
ipa_2$Upstream.Regulator<-as.character(ipa_2$Upstream.Regulator)
ipa_2$NAME<-ifelse(ipa_2$Activation.z.score>=4 |ipa_2$Activation.z.score<=-2 & ipa_2$p.value.of.overlap< 0.05,
                   ipa_2$Upstream.Regulator,"" )
ipa_2$NAME<-ifelse(abs(ipa_2$Activation.z.score)>=6 & ipa_2$p.value.of.overlap< 0.05,
                   ipa_2$Upstream.Regulator,"" )



#ipa_2$Predicted.Activation.State
ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(size=3.5)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("black","red","blue"))+
  geom_point()+geom_text(nudge_y = 1,size=3)+ylab("-log10(p-value)")

ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(size=.1)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("red","blue","black"))+
  geom_point()+ylab("-log10(p-value)")+geom_text(nudge_y = 1,size=3)





ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(aes(color = Predicted.Activation.State, size = log_p), show.legend = F)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("red","blue","black"))+
  geom_point()+ylab("-log10(p-value)")+geom_text(nudge_y = 1,size=3)

library(tidyverse)
library(ggrepel)
library(ggridges)
library(lubridate)
library(modelr)
library(patchwork)

ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(size=.1)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("red","blue","black"))+
  geom_point()+ylab("-log10(p-value)")+geom_text_repel(nudge_y = 1,size=3)
ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(aes(color = Predicted.Activation.State, size = log_p), show.legend = F)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("red","blue","black"))+
  geom_point()+ylab("-log10(p-value)")+geom_text_repel(nudge_y = 1,size=3)


pdf("E2_upstream.pdf",width = 6,height = 6)
ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(aes(color = Predicted.Activation.State, size = log_p), show.legend = F)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("red","blue","black"))+
  geom_point()+ylab("-log10(p-value)")+geom_text_repel(
    data = ipa_2, aes(label = NAME),
    size = 3, fontface = 3,
    nudge_x = 2, nudge_y = 2) +ggtitle("IPA results of E2")+
  theme(plot.title = element_text(hjust = 0.5))  #也就加上这一行
dev.off()


setwd("C:/Users/yudonglin/Desktop/硕士论文写作/硕士论文/figure-整理/IPA分析")
library(ggplot2)

ipa<-read.csv("Upstream Regulators_E3.csv")
#ipa_2<-ipa[,c(1,4,5,7)]
ipa_2<-ipa[,c(3,6,7,11)]
head(ipa_2)
ipa_2[is.na(ipa_2)]<-0
colnames(ipa_2)
ipa_2$log_p<--log10(ipa_2$p.value.of.overlap)
ipa_2$Upstream.Regulator<-as.character(ipa_2$Upstream.Regulator)
ipa_2$NAME<-ifelse(abs(ipa_2$Activation.z.score)>=2 & ipa_2$p.value.of.overlap< 0.05,
                   ipa_2$Upstream.Regulator,"" )
ipa_2$NAME<-ifelse(abs(ipa_2$Activation.z.score)>=6 & ipa_2$p.value.of.overlap< 0.05,
                   ipa_2$Upstream.Regulator,"" )



#ipa_2$Predicted.Activation.State
ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(size=3.5)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("black","red","blue"))+
  geom_point()+geom_text(nudge_y = 1,size=3)+ylab("-log10(p-value)")

ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(size=.1)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("red","blue","black"))+
  geom_point()+ylab("-log10(p-value)")+geom_text(nudge_y = 1,size=3)





ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(aes(color = Predicted.Activation.State, size = log_p), show.legend = F)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("red","blue","black"))+
  geom_point()+ylab("-log10(p-value)")+geom_text(nudge_y = 1,size=3)

library(tidyverse)
library(ggrepel)
library(ggridges)
library(lubridate)
library(modelr)
library(patchwork)

ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(size=.1)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("red","blue","black"))+
  geom_point()+ylab("-log10(p-value)")+geom_text_repel(nudge_y = 1,size=3)
ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(aes(color = Predicted.Activation.State, size = log_p), show.legend = F)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("red","blue","black"))+
  geom_point()+ylab("-log10(p-value)")+geom_text_repel(nudge_y = 1,size=3)


pdf("E3_upstream.pdf",width = 6,height = 6)
ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(aes(color = Predicted.Activation.State, size = log_p), show.legend = F)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("red","blue","black"))+
  geom_point()+ylab("-log10(p-value)")+geom_text_repel(
    data = ipa_2, aes(label = NAME),
    size = 3, fontface = 3,
    nudge_x = 1, nudge_y = 1) +ggtitle("IPA results of E3")+
  theme(plot.title = element_text(hjust = 0.5))  #也就加上这一行
dev.off()


setwd("C:/Users/yudonglin/Desktop/硕士论文写作/硕士论文/figure-整理/IPA分析")
library(ggplot2)

ipa<-read.csv("Upstream Regulators_E4.csv")
#ipa_2<-ipa[,c(1,4,5,7)]
ipa_2<-ipa[,c(3,6,7,11)]
head(ipa_2)
ipa_2[is.na(ipa_2)]<-0
colnames(ipa_2)
ipa_2$log_p<--log10(ipa_2$p.value.of.overlap)
ipa_2$Upstream.Regulator<-as.character(ipa_2$Upstream.Regulator)
ipa_2$NAME<-ifelse(abs(ipa_2$Activation.z.score)>=2 & ipa_2$p.value.of.overlap< 0.05,
                   ipa_2$Upstream.Regulator,"" )
ipa_2$NAME<-ifelse(abs(ipa_2$Activation.z.score)>=6 & ipa_2$p.value.of.overlap< 0.05,
                   ipa_2$Upstream.Regulator,"" )



#ipa_2$Predicted.Activation.State
ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(size=3.5)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("black","red","blue"))+
  geom_point()+geom_text(nudge_y = 1,size=3)+ylab("-log10(p-value)")

ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(size=.1)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("red","blue","black"))+
  geom_point()+ylab("-log10(p-value)")+geom_text(nudge_y = 1,size=3)





ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(aes(color = Predicted.Activation.State, size = log_p), show.legend = F)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("red","blue","black"))+
  geom_point()+ylab("-log10(p-value)")+geom_text(nudge_y = 1,size=3)

library(tidyverse)
library(ggrepel)
library(ggridges)
library(lubridate)
library(modelr)
library(patchwork)

ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(size=.1)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("red","blue","black"))+
  geom_point()+ylab("-log10(p-value)")+geom_text_repel(nudge_y = 1,size=3)
ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(aes(color = Predicted.Activation.State, size = log_p), show.legend = F)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("red","blue","black"))+
  geom_point()+ylab("-log10(p-value)")+geom_text_repel(nudge_y = 1,size=3)


pdf("E4_upstream.pdf",width = 6,height = 6)
ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(aes(color = Predicted.Activation.State, size = log_p), show.legend = F)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("red","blue","black"))+
  geom_point()+ylab("-log10(p-value)")+geom_text_repel(
    data = ipa_2, aes(label = NAME),
    size = 3, fontface = 3,
    nudge_x = 1, nudge_y = 1) +ggtitle("IPA results of E4")+
  theme(plot.title = element_text(hjust = 0.5))  #也就加上这一行
dev.off()


setwd("C:/Users/yudonglin/Desktop/硕士论文写作/硕士论文/figure-整理/IPA分析")
library(ggplot2)

ipa<-read.csv("Upstream Regulators_E5.csv")
#ipa_2<-ipa[,c(1,4,5,7)]
ipa_2<-ipa[,c(3,6,7,11)]
head(ipa_2)
ipa_2[is.na(ipa_2)]<-0
colnames(ipa_2)
ipa_2$log_p<--log10(ipa_2$p.value.of.overlap)
ipa_2$Upstream.Regulator<-as.character(ipa_2$Upstream.Regulator)
ipa_2$NAME<-ifelse(abs(ipa_2$Activation.z.score)>=2 & ipa_2$p.value.of.overlap< 0.05,
                   ipa_2$Upstream.Regulator,"" )
ipa_2$NAME<-ifelse(abs(ipa_2$Activation.z.score)>=6 & ipa_2$p.value.of.overlap< 0.05,
                   ipa_2$Upstream.Regulator,"" )



#ipa_2$Predicted.Activation.State
ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(size=3.5)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("black","red","blue"))+
  geom_point()+geom_text(nudge_y = 1,size=3)+ylab("-log10(p-value)")

ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(size=.1)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("red","blue","black"))+
  geom_point()+ylab("-log10(p-value)")+geom_text(nudge_y = 1,size=3)





ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(aes(color = Predicted.Activation.State, size = log_p), show.legend = F)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("red","blue","black"))+
  geom_point()+ylab("-log10(p-value)")+geom_text(nudge_y = 1,size=3)

library(tidyverse)
library(ggrepel)
library(ggridges)
library(lubridate)
library(modelr)
library(patchwork)

ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(size=.1)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("red","blue","black"))+
  geom_point()+ylab("-log10(p-value)")+geom_text_repel(nudge_y = 1,size=3)
ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(aes(color = Predicted.Activation.State, size = log_p), show.legend = F)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("red","blue","black"))+
  geom_point()+ylab("-log10(p-value)")+geom_text_repel(nudge_y = 1,size=3)


pdf("E5_upstream.pdf",width = 6,height = 6)
ggplot(ipa_2,aes(Activation.z.score,log_p,colour=Predicted.Activation.State,label=NAME))+
  geom_point(aes(color = Predicted.Activation.State, size = log_p), show.legend = F)+geom_vline(xintercept = c(-2,2),lty=2,color="grey")+
  geom_hline(yintercept = -log10(0.05),lty=4,color="grey")+theme_bw()+
  scale_color_manual(values = c("blue","black"))+
  geom_point()+ylab("-log10(p-value)")+geom_text_repel(
    data = ipa_2, aes(label = NAME),
    size = 3, fontface = 3,
    nudge_x =1, nudge_y = 1) +ggtitle("IPA results of E5")+
  theme(plot.title = element_text(hjust = 0.5))  #也就加上这一行
dev.off()


