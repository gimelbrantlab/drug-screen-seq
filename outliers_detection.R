library(tidyverse)
setwd("F:/MAEscreenDataAnalysis/Secondary_screen")

df <- read_csv("screen_seq_data_complete.csv")
df <- df %>% dplyr::select(readout_gene, target_gene, bias, num_reads) %>%
  filter(!is.na(bias)) %>%
  filter(num_reads >= 50)
colnames(df)[2] <- "drug"
df <- df %>% group_by(readout_gene) %>%
  mutate(lowerbound = quantile(bias, .15) - 3 * IQR(bias)) %>%
  mutate(upperbound = quantile(bias, .85) + 3 * IQR(bias)) %>%
  ungroup()
df <- df %>%
  dplyr::mutate(drug = factor(drug, levels = c("Control", "Garcinol", "Splitomicin", "BML-210",
                                                             "5-Aza-2-deoxycytidine", "Zebularine", "Tranylcypromine_hemisulfate",
                                                             "EX-527", "Nicotinamide", "BML-266", "Piceatannol", "Fluoro-SAHA",
                                                             "AGK2", "Salermide", "Anacardic_acid", "B2", "CTPB", "Sirtinol",
                                                             "Suramin6Na", "BML-278", "CI-994", "NSC-3852", "Aminoresveratrol_sulfate",
                                                             "BML-281", "Triacetylresveratrol")))
df$candidate <- "No"
df$candidate[df$bias < df$lowerbound | df$bias > df$upperbound ] <- "Yes"
ggplot(df, aes(drug, bias, col=candidate), alpha = .2) +
  geom_point(size = 1) +
  scale_color_manual(values=c("darkgrey", "red", "blue"))+
  facet_wrap(~readout_gene) +
  theme_bw()+
  theme(strip.background = element_blank(), strip.text = element_text(size = 10, margin = margin())) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.1, size = 4))+
  geom_vline(xintercept=c(5,10,15,20), color="grey",size=2,alpha=0.2)+
  scale_y_continuous(breaks=c(0,.5,1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

###############################################
df <- read_csv("screen_seq_data_complete.csv")
df <- df %>% dplyr::select(readout_gene, target_gene, bias, num_reads) %>%
  filter(!is.na(bias)) %>%
  filter(num_reads >= 50)
colnames(df)[2] <- "drug"

Nsigmas = 4

df <- df %>% group_by(readout_gene) %>%
  mutate(lowerbound = mean(bias) - Nsigmas * sd(bias)) %>%
  mutate(upperbound = mean(bias) + Nsigmas * sd(bias)) %>%
  # mutate(signal = (abs(bias - mean(bias)))/sd(bias)) %>% 
  mutate(signal = (abs(bias - mean(bias)))) %>% 
  ungroup()

df <- df %>%
  dplyr::mutate(drug = factor(drug, levels = c("Control", "Garcinol", "Splitomicin", "BML-210",
                                               "5-Aza-2-deoxycytidine", "Zebularine", "Tranylcypromine_hemisulfate",
                                               "EX-527", "Nicotinamide", "BML-266", "Piceatannol", "Fluoro-SAHA",
                                               "AGK2", "Salermide", "Anacardic_acid", "B2", "CTPB", "Sirtinol",
                                               "Suramin6Na", "BML-278", "CI-994", "NSC-3852", "Aminoresveratrol_sulfate",
                                               "BML-281", "Triacetylresveratrol")))
df$candidate <- "No"
df$candidate[df$bias < df$lowerbound | df$bias > df$upperbound ] <- "Yes"
ggplot(df, aes(drug, bias, col=candidate), alpha = .2) +
  geom_point(size = 1) +
  scale_color_manual(values=c("darkgrey", "red", "blue"))+
  facet_wrap(~readout_gene) +
  theme_bw()+
  theme(strip.background = element_blank(), strip.text = element_text(size = 10, margin = margin())) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.1, size = 4))+
  geom_vline(xintercept=c(5,10,15,20), color="grey",size=2,alpha=0.2)+
  # geom_hline(yintercept=df$upperbound[readout_gene], color="green",size=2,alpha=0.4)+
  # geom_hline(yintercept=c(lowerbound), color="green",size=2,alpha=0.4)+
  
  scale_y_continuous(breaks=c(0,.5,1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

df1 = df %>% arrange(signal)
p = ggplot(df1, aes(signal, bias, col=candidate))
p + geom_point()+
  scale_color_manual(values=c("darkgrey", "red", "blue"))

###############################################
xx = df %>% filter(readout_gene == "Adamtsl4")
Nsigmas = 3
mu = mean(xx$bias)
stdev = sd(xx$bias)
mu.plus=mu + Nsigmas*stdev
mu.minus=mu - Nsigmas*stdev

p = ggplot(xx, aes(drug, bias))
p +
  geom_point(size = 1) +
  scale_color_manual(values=c("darkgrey", "white", "blue"))+
  theme_bw()+
  theme(strip.background = element_blank(), strip.text = element_text(size = 10, margin = margin())) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.1))+
  geom_vline(xintercept=c(5,10,15,20), color="grey",size=2,alpha=0.2)+
  geom_hline(yintercept=c(mu), color="green",size=2,alpha=0.2)+
  geom_hline(yintercept=c(mu.plus), color="grey",size=2,alpha=0.4)+
  geom_hline(yintercept=c(mu.minus), color="grey",size=2,alpha=0.4)+
  geom_rect(inherit.aes=FALSE,
            aes(xmin=0,xmax=10,ymin=mu.plus,ymax=mu.minus, fill = c("grey"), alpha=0.2))+
  
  scale_y_continuous(breaks=c(0,.5,1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

