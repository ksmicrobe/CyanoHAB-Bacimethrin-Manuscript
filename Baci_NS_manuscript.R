###Packages to load 
library(tidyverse)
library(readxl)
library(dplyr)
library(ggplot2)
library(rstatix)
library(RColorBrewer)
library(gg.gap)
###

#Global note: enrichment plots are done first showing AFE then showing Xnet

###Necessary functions
calc_APE <- function(permil_samp, isotope, permil_bg){
  #Define natural abundance atom fractions from isotopes of interest
  natural_abund <- list("13C" = 0.011237, 
                        "15N" = 0.00367) #Can modify for other isotopes 
  if(!isotope %in% names(natural_abund)){
    stop("Unrecognized isotope; options include: '13C' and '15N'") #Can modify for your isotopes 
  }
  #Define standard ratio based on the chosen isotope
  std_ratio <- natural_abund[[isotope]]
  #Convert permil values to isotope fractions (R) then atom percents (f)
  R_samp <- std_ratio * (1 + permil_samp / 1000)
  f_samp <- R_samp / (1 + R_samp)
  R_bg <- std_ratio * (1 + permil_bg / 1000)
  f_bg <- R_bg / (1 + R_bg)
  #Calculate atom percent excess
  APE <- (f_samp - f_bg) * 100
  return(APE)
}
###

###F. LE24H monoculture C, N preference
le24h_sept_2025 <- read_xlsx(path = "~/OneDrive/Documents/Sept2025_LE24H_NS.xlsx") #Replace with your path to xlsx file
#Need to break this into individual treatments to calculate APE normalized that treatment's killed control
Glu_NH4_le24h_sept_2025 <- le24h_sept_2025 %>% filter(Treatment == "Glu_NH4")
Prot_noNitrgn_le24h_sept_2025 <- le24h_sept_2025 %>% filter(Treatment == "Prot_noNitrgn")
Glu_NH4_noNitrgn_le24h_sept_2025 <- le24h_sept_2025 %>% filter(Treatment == "Glu_NH4_noNitrgn")
Glu_NO3_noNitrgn_le24h_sept_2025 <- le24h_sept_2025 %>% filter(Treatment == "Glu_NO3_noNitrgn")

Glu_NH4_le24h_sept_2025_ape <- Glu_NH4_le24h_sept_2025 %>% mutate(APE_13C = 
                                                      calc_APE(permil_samp = C_permil, 
                                                               isotope = "13C", 
                                                               permil_bg = 89.9), #These will be unique to each killed control 
                                                    APE_15N = 
                                                      calc_APE(permil_samp = N_permil, 
                                                               isotope = "15N", 
                                                               permil_bg = 1755.9)) 

Prot_noNitrgn_le24h_sept_2025 <- Prot_noNitrgn_le24h_sept_2025 %>% mutate(APE_13C = 
                                                                    calc_APE(permil_samp = C_permil, 
                                                                             isotope = "13C", 
                                                                             permil_bg = 154.3), 
                                                                  APE_15N = 
                                                                    calc_APE(permil_samp = N_permil, 
                                                                             isotope = "15N", 
                                                                             permil_bg = 2020.3)) 

Glu_NH4_noNitrgn_le24h_sept_2025 <- Glu_NH4_noNitrgn_le24h_sept_2025 %>% mutate(APE_13C = 
                                                                    calc_APE(permil_samp = C_permil, 
                                                                             isotope = "13C", 
                                                                             permil_bg = 89.9), 
                                                                  APE_15N = 
                                                                    calc_APE(permil_samp = N_permil, 
                                                                             isotope = "15N", 
                                                                             permil_bg = 1755.9)) 

Glu_NO3_noNitrgn_le24h_sept_2025 <- Glu_NO3_noNitrgn_le24h_sept_2025 %>% mutate(APE_13C = 
                                                                    calc_APE(permil_samp = C_permil, 
                                                                             isotope = "13C", 
                                                                             permil_bg = 89.9), 
                                                                  APE_15N = 
                                                                    calc_APE(permil_samp = N_permil, 
                                                                             isotope = "15N", 
                                                                             permil_bg = 1001.4)) 

le24h_sept_2025_ape_1 <- rbind(Glu_NO3_noNitrgn_le24h_sept_2025, Glu_NH4_noNitrgn_le24h_sept_2025)
le24h_sept_2025_ape_2 <- rbind(le24h_sept_2025_ape_1, Prot_noNitrgn_le24h_sept_2025)
le24h_sept_2025_ape <- rbind(le24h_sept_2025_ape_2, Glu_NH4_le24h_sept_2025_ape)
unique(le24h_sept_2025_ape$Treatment) #Now we'll just use this since we've now corrected for killed controls

##Filter by numerators of isotope ratios 
totnum1_min <- quantile(le24h_sept_2025_ape$TOTNUM1, 0.30)
totnum2_min <- quantile(le24h_sept_2025_ape$TOTNUM2, 0.30)

dat_le24h_sept2025_filt <- le24h_sept_2025_ape %>% 
filter(TOTNUM1 > totnum1_min & TOTNUM2 > totnum2_min) #Keep only the 90th percentile of numerator values

dim(le24h_sept_2025_ape)
dim(dat_le24h_sept2025_filt) #See how many ROIs this dropped—should mainly be negative per mil values
##

##C and N boxplots by treatment 
jitter_pos <- position_jitter(width = 0.1, seed = 123)
dat_le24h_sept2025_filt_nokill <- dat_le24h_sept2025_filt
dat_le24h_sept2025_filt_nokill$Treatment <- factor(dat_le24h_sept2025_filt_nokill$Treatment, 
                                                   levels = c("Glu_NH4", "Glu_NH4_noNitrgn", "Glu_NO3_noNitrgn", "Prot_noNitrgn"))
#N AFE
dat_le24h_sept2025_filt_nokill %>% ggplot(aes(x = Treatment, y = APE_15N/100)) +
  geom_point(size = 1, position = jitter_pos, color = "gray35") + geom_boxplot(fill = NA) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
                                           axis.text.y = element_text(size = 12, color = "black"),
                                           axis.title = element_text(color = "black", size = 12)) +
  labs(y = expression(""^"15"*"N AFE"), x = "") +
  scale_x_discrete(labels = 
                     c('Prot_noNitrgn' = expression(""^"13"*"C, "^"15"*"N-Prot, "^"12"*"C-Glu"), 
                       'Glu_NH4' = expression(""^"13"*"C-Glu, "^"15"*"NH "[4]^"+"*", "^"14"*"NO"[3]^"-"), 
                       'Glu_NH4_noNitrgn' = expression(""^"13"*"C-Glu, "^"15"*"NH"[4]^"+"), 
                       'Glu_NO3_noNitrgn' = expression(""^"13"*"C-Glu, "^"15"*"NO"[3]^"-"))) +
  scale_y_continuous(breaks = c(0, .10, .20, .30, .40, .50))
#

#N net
Xnet_treatments <- c("Glu_NH4", "Glu_NH4_noNitrgn", "Prot_noNitrgn")
dat_le24h_sept2025_filt_nokill %>% filter(Treatment %in% Xnet_treatments) %>% 
  ggplot(aes(x = Treatment, y = Nnet*100)) +
  geom_point(size = 1, position = jitter_pos, color = "gray35") + geom_boxplot(fill = NA) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_text(color = "black", size = 12)) +
  labs(y = expression(N[net]*" (%)"), x = "") +
  scale_x_discrete(labels = 
                     c('Prot_noNitrgn' = expression(""^"13"*"C, "^"15"*"N-Prot, "^"12"*"C-Glu"), 
                       'Glu_NH4' = expression(""^"13"*"C-Glu, "^"15"*"NH "[4]^"+"*", "^"14"*"NO"[3]^"-"), 
                       'Glu_NH4_noNitrgn' = expression(""^"13"*"C-Glu, "^"15"*"NH"[4]^"+"))) +
  scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50))
#

#C AFE
dat_le24h_sept2025_filt_nokill %>% ggplot(aes(x = Treatment, y = APE_13C/100)) +
  geom_point(size = 1, position = jitter_pos, color = "gray35") + geom_boxplot(fill = NA) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
                                           axis.text.y = element_text(size = 12, color = "black"),
                                           axis.title = element_text(color = "black", size = 12)) +
  labs(y = expression(""^"13"*"C AFE"), x = "") +
  scale_x_discrete(labels = 
                     c('Prot_noNitrgn' = expression(""^"13"*"C, "^"15"*"N-Prot, "^"12"*"C-Glu"), 
                       'Glu_NH4' = expression(""^"13"*"C-Glu, "^"15"*"NH "[4]^"+"*", "^"14"*"NO"[3]^"-"), 
                       'Glu_NH4_noNitrgn' = expression(""^"13"*"C-Glu, "^"15"*"NH"[4]^"+"), 
                       'Glu_NO3_noNitrgn' = expression(""^"13"*"C-Glu, "^"15"*"NO"[3]^"-")))
#

#C net
dat_le24h_sept2025_filt_nokill %>% 
  ggplot(aes(x = Treatment, y = Cnet*100)) +
  geom_point(size = 1, position = jitter_pos, color = "gray35") + geom_boxplot(fill = NA) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_text(color = "black", size = 12)) +
  labs(y = expression(C[net]*" (%)"), x = "") +
  scale_x_discrete(labels = 
                     c('Prot_noNitrgn' = expression(""^"13"*"C, "^"15"*"N-Prot, "^"12"*"C-Glu"), 
                       'Glu_NH4' = expression(""^"13"*"C-Glu, "^"15"*"NH "[4]^"+"*", "^"14"*"NO"[3]^"-"), 
                       'Glu_NH4_noNitrgn' = expression(""^"13"*"C-Glu, "^"15"*"NH"[4]^"+"), 
                       'Glu_NO3_noNitrgn' = expression(""^"13"*"C-Glucose, "^"15"*"N-NO"[3]^"-"))) +
  scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50), limits = c(0, 50))
#
##

##Enrichment scatterplots and correlations 
comb_le24h <- dat_le24h_sept2025_filt_nokill %>% select(APE_13C, APE_15N, Treatment, Replicate, Nnet, Cnet)

#AFE
comb_le24h %>% 
  mutate(Treatment = factor(Treatment, levels = c("Glu_NH4_noNitrgn", 
                                                  "Glu_NH4", 
                                                  "Prot_noNitrgn", 
                                                  "Glu_NO3_noNitrgn"))) %>% 
  ggplot(aes(x = APE_15N/100, y = APE_13C/100, color = Treatment)) +
  geom_point(size = 1) + theme_classic() + 
  scale_color_brewer(palette = "Dark2", 
                     labels = c(expression("+"^"13"*"C-Glucose, +"^"15"*"NH"[4]^"+"), 
                                expression("+"^"13"*"C-Glucose, +"^"15"*"NH"[4]^"+"*", +"^"14"*"NO"[3]^"-"), 
                                expression("+"^"13"*"C, "^" 15"*"N-Protein, +"^"12"*"C-Glucose"), 
                                expression("+"^"13"*"C-Glucose, +"^"15"*"N-NO"[3]^"-"))) +
  theme(legend.position = "none", axis.text = element_text(color = "black", size = 12), 
        axis.title = element_text(color = "black", size = 12)) + geom_smooth(
          data = . %>% filter(Treatment %in% c("Glu_NH4", "Glu_NH4_noNitrgn", "Prot_noNitrgn")),
          se = FALSE, method = "lm", fullrange = TRUE
        ) + coord_cartesian(xlim = c(0, .60), ylim = c(0, NA)) + 
  labs(x = expression(""^"15"*"N AFE"), y = expression(""^"13"*"C AFE")) +
  scale_x_continuous(breaks = c(.0, .10, .20, .30, .40, .50, .60))

cor.test(Prot_noNitrgn_le24h_sept_2025$APE_13C/100, Prot_noNitrgn_le24h_sept_2025$APE_15N/100, 
         alternative = "two.sided", "spearman")
cor.test(Glu_NH4_le24h_sept_2025_ape$APE_13C/100, Glu_NH4_le24h_sept_2025_ape$APE_15N/100, 
         alternative = "two.sided", "spearman")
cor.test(Glu_NH4_noNitrgn_le24h_sept_2025$APE_13C/100, Glu_NH4_noNitrgn_le24h_sept_2025$APE_15N/100, 
         alternative = "two.sided", "spearman")
#AFE

#Xnet
comb_le24h %>% 
  mutate(Treatment = factor(Treatment, levels = c("Glu_NH4_noNitrgn", 
                                                  "Glu_NH4", 
                                                  "Prot_noNitrgn", 
                                                  "Glu_NO3_noNitrgn"))) %>% 
  filter(Treatment %in% Xnet_treatments) %>% 
  ggplot(aes(x = Nnet*100, y = Cnet*100, color = Treatment)) +
  geom_point(size = 1) + theme_classic() + 
  scale_color_brewer(palette = "Dark2", 
                     labels = c(expression("+"^"13"*"C-Glucose, +"^"15"*"NH"[4]^"+"), 
                                expression("+"^"13"*"C-Glucose, +"^"15"*"NH"[4]^"+"*", +"^"14"*"NO"[3]^"-"), 
                                expression("+"^"13"*"C, "^" 15"*"N-Protein, +"^"12"*"C-Glucose"))) +
  theme(legend.position = "none", axis.text = element_text(color = "black", size = 12), 
        axis.title = element_text(color = "black", size = 12)) + geom_smooth(
          data = . %>% filter(Treatment %in% c("Glu_NH4", "Glu_NH4_noNitrgn", "Prot_noNitrgn")),
          se = FALSE, method = "lm", fullrange = TRUE
        ) + ylim(0, NA) + labs(x = expression(N[net]*" (%)"), y = expression(C[net]*" (%)")) +
  scale_x_continuous(limits = c(0, 60), breaks = c(0, 10, 20, 30, 40, 50, 60))

cor.test(Prot_noNitrgn_le24h_sept_2025$Nnet, Prot_noNitrgn_le24h_sept_2025$Cnet, 
         alternative = "two.sided", "spearman")
cor.test(Glu_NH4_le24h_sept_2025_ape$Nnet, Glu_NH4_le24h_sept_2025_ape$Cnet, 
         alternative = "two.sided", "spearman")
cor.test(Glu_NH4_noNitrgn_le24h_sept_2025$Nnet, Glu_NH4_noNitrgn_le24h_sept_2025$Cnet, 
         alternative = "two.sided", "spearman")
#
##
###

###T. variabilis-F. LE24H C, N transfer 
A_le24h_sept_2025 <- read_xlsx(path = "~/OneDrive/Documents/Anabaena_data_Sept2025.xlsx") #Replace with your path to xlsx file
A_le24h_sept_2025 <- A_le24h_sept_2025 %>% mutate(
  treatment_merged = str_replace_all(Treatment, "_\\d+", "")) 
A_le24h_sept_2025_ape <- A_le24h_sept_2025 %>% 
  mutate(APE_13C = calc_APE(permil_13C, '13C', -46.3), 
         APE_15N = calc_APE(permil_15N, '15N', -411.9))
#These ROIs aren't filtered because all were drawn manually

##Boxplots 

#15N AFE
A_le24h_sept_2025_ape %>% filter(Identity == "Anabaena") %>% filter(Filter != "F5") %>% 
  mutate(APE_13C = APE_13C/100, APE_15N = APE_15N/100) %>% 
  ggplot(aes(y = APE_15N, x = treatment_merged)) + 
  geom_point(size = 2) + geom_boxplot(fill = NA) +
  theme_classic() + labs(x = "Culture", y = expression(""^"15"*"N AFE"), 
                         title = expression(italic("T. variabilis")*" ROIs")) +
  scale_x_discrete(labels = c(expression(italic("T. variabilis")), 
                              expression(italic("T. variabilis")*", "*italic("F")*". LE24H"))) +
  theme(axis.text = element_text(size = 12, color = "black"), 
        axis.title = element_text(size = 12, color = "black")) + 
  ggpubr::stat_compare_means(method = "wilcox.test", hide.ns = T, 
                             comparisons = list(c("A", "A_LE24H")))

A_le24h_sept_2025_ape %>% filter(Identity == "Heterotroph") %>% filter(Filter != "F5") %>% 
  mutate(APE_13C = APE_13C/100, APE_15N = APE_15N/100) %>% 
  ggplot(aes(y = APE_15N, x = treatment_merged)) + 
  geom_point(size = 2) + geom_boxplot(fill = NA) +
  theme_classic() + labs(x = "Culture", y = expression(""^"15"*"N AFE"), 
                         title = "Heterotroph ROIs") +
  scale_x_discrete(labels = c(expression(italic("T. variabilis")), 
                              expression(italic("T. variabilis")*", "*italic("F")*". LE24H"))) +
  theme(axis.text = element_text(size = 12, color = "black"), 
        axis.title = element_text(size = 12, color = "black")) + 
  ggpubr::stat_compare_means(method = "wilcox.test", hide.ns = T, 
                             comparisons = list(c("A", "A_LE24H"))) +
  scale_y_continuous(breaks = c(0, .1, .2, .3, .4), limits = c(0, .4))
#

#Nnet
A_le24h_sept_2025_ape %>% filter(Identity == "Anabaena") %>% filter(Filter != "F5") %>% 
  ggplot(aes(y = Nnet*100, x = treatment_merged)) + 
  geom_point(size = 2) + geom_boxplot(fill = NA) +
  theme_classic() + labs(x = "Culture", y = expression(N[net]*" (%)"), 
                         title = expression(italic("T. variabilis")*" ROIs")) +
  scale_x_discrete(labels = c(expression(italic("T. variabilis")), 
                              expression(italic("T. variabilis")*", "*italic("F")*". LE24H"))) +
  theme(axis.text = element_text(size = 12, color = "black"), 
        axis.title = element_text(size = 12, color = "black")) + 
  ggpubr::stat_compare_means(method = "wilcox.test", hide.ns = T, 
                             comparisons = list(c("A", "A_LE24H")))

A_le24h_sept_2025_ape %>% filter(Identity == "Heterotroph") %>% filter(Filter != "F5") %>% 
  ggplot(aes(y = Nnet*100, x = treatment_merged)) + 
  geom_point(size = 2) + geom_boxplot(fill = NA) +
  theme_classic() + labs(x = "Culture", y = expression(N[net]*" (%)"), 
                         title = "Heterotroph ROIs") +
  scale_x_discrete(labels = c(expression(italic("T. variabilis")), 
                              expression(italic("T. variabilis")*", "*italic("F")*". LE24H"))) +
  theme(axis.text = element_text(size = 12, color = "black"), 
        axis.title = element_text(size = 12, color = "black")) + 
  ggpubr::stat_compare_means(method = "wilcox.test", hide.ns = T, 
                             comparisons = list(c("A", "A_LE24H"))) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100), limits = c(0, 100))
#
##

##Enrichment scatterplots and correlations
#15N AFE
A_le24h_sept_2025_ape %>% filter(Identity == "Anabaena") %>% filter(Filter != "F5") %>%
  mutate(APE_13C = APE_13C/100, APE_15N = APE_15N/100) %>% 
  ggplot(aes(y = APE_13C, x = APE_15N, shape = treatment_merged)) + 
  geom_point(size = 2) + theme_classic() + 
  scale_shape_manual(values = c(19,1), 
                     labels = c(expression(italic("T. variabilis")), 
                                expression(italic("T. variabilis")*", "*italic("F")*". LE24H")), 
                     name = "Culture") +
  geom_smooth(aes(group = treatment_merged, linetype = treatment_merged), 
              method = lm, se = F, color = "black", linewidth = 0.7) +
  scale_linetype_manual(values = c("solid", "dashed"),
                        labels = c(expression(italic("T. variabilis")), 
                                   expression(italic("T. variabilis")*", "*italic("F")*". LE24H")),
                        name = "Culture") +
  labs(x = expression(""^"15"*"N AFE"), y = expression(""^"13"*"C AFE"), color = "Culture") + 
  theme(legend.position = "top", 
        axis.text = element_text(size = 12, color = "black"), 
        axis.title = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12, color = "black")) + 
  ggtitle(expression(italic("T. variabilis")*" ROIs")) +
  scale_y_continuous(breaks = c(.1, .2, .3), limits = c(0, .3))

Anabaena_A <- A_le24h_sept_2025_ape %>% filter(Identity == "Anabaena") %>% filter(Filter %in% c("F1", "F2"))
Anabaena_A_F <- A_le24h_sept_2025_ape %>% filter(Identity == "Anabaena") %>% filter(Filter %in% c("F3", "F4"))
cor.test(Anabaena_A$APE_13C/100, Anabaena_A$APE_15N/100, 
         alternative = "two.sided", "spearman")
cor.test(Anabaena_A_F$APE_13C/100, Anabaena_A_F$APE_15N/100, 
         alternative = "two.sided", "spearman")
#15N AFE

#Nnet
A_le24h_sept_2025_ape %>% filter(Identity == "Anabaena") %>% filter(Filter != "F5") %>%
  mutate(APE_13C = Cnet*100, APE_15N = Nnet*100) %>% 
  ggplot(aes(y = APE_13C, x = APE_15N, shape = treatment_merged)) + 
  geom_point(size = 2) + theme_classic() + 
  scale_shape_manual(values = c(19,1), 
                     labels = c(expression(italic("T. variabilis")), 
                                expression(italic("T. variabilis")*", "*italic("F")*". LE24H")), 
                     name = "Culture") +
  geom_smooth(aes(group = treatment_merged, linetype = treatment_merged), 
              method = lm, se = F, color = "black", linewidth = 0.7) +
  scale_linetype_manual(values = c("solid", "dashed"),
                        labels = c(expression(italic("T. variabilis")), 
                                   expression(italic("T. variabilis")*", "*italic("F")*". LE24H")),
                        name = "Culture") +
  labs(x = expression(N[net]*" (%)"), y = expression(C[net]*" (%)"), color = "Culture") + 
  theme(legend.position = "top", 
        axis.text = element_text(size = 12, color = "black"), 
        axis.title = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12, color = "black")) + 
  ggtitle(expression(italic("T. variabilis")*" ROIs")) +
  scale_y_continuous(limits = c(0, 30), breaks = c(0, 10, 20, 30))

cor.test(Anabaena_A$Cnet, Anabaena_A$Nnet, 
         alternative = "two.sided", "spearman")
cor.test(Anabaena_A_F$Cnet, Anabaena_A_F$Nnet, 
         alternative = "two.sided", "spearman")
#
##
###

###F. LE24H-T. variabilis bacimethrin additions
flavo_baci_oct2025 <- read_excel(path = "/Users/kellyshannon/Library/CloudStorage/OneDrive-Personal/Documents/flavo_NS_data_Oct2025.xlsx.xlsx")

flavo_baci_oct2025_ape <- flavo_baci_oct2025 %>% 
  mutate(APE_13C = calc_APE(R1, "13C", 0), 
         APE_15N = calc_APE(R2, "15N", 0))

flavo_baci_oct2025_ape$Treatment <- factor(flavo_baci_oct2025_ape$Treatment, 
                                           levels = c("Normal", "+Bacimethrin"))

totnum2_min_F_baci <- quantile(flavo_baci_oct2025_ape$TOTNUM2, 0.10)
flavo_baci_oct2025_ape_filt <- flavo_baci_oct2025_ape %>% 
  filter(TOTNUM2 > totnum2_min_F_baci)

flavo_baci_oct2025_ape_filt %>% ggplot(aes(x = Treatment, y = APE_15N/100)) +
  geom_point(size = 1, position = jitter_pos, color = "gray35") + geom_boxplot(fill = NA) +
  theme_classic() + theme(axis.text.x = 
                                             element_text(angle = 45, hjust = 1, colour = "black", size = 12), 
                                           axis.text.y = 
                                             element_text(colour = "black", size = 12), 
                                           axis.title = 
                                             element_text(colour = "black", size = 12)) +
  labs(y = expression(""^"15"*"N AFE")) +
  ggpubr::stat_compare_means(method = "wilcox.test", hide.ns = T, 
                             comparisons = list(c("Normal", "+Bacimethrin")))
###

###F. LE24H-cyanoHAB coculture experiments
NS_data_1 <- read_excel(path = "/Users/kellyshannon/Library/CloudStorage/OneDrive-Personal/Documents/Feb2026_NS_data_2.xlsx")
min_totnum2 <- quantile(NS_data_1$TOTNUM2, 0.3)
NS_data_filt <- NS_data_1 %>% filter(TOTNUM2 > min_totnum2)
NS_data_DIN <- NS_data_filt %>% filter(N_Substrate == "Ammonium")
NS_data_N2 <- NS_data_filt %>% filter(N_Substrate == "N2")

NS_data_DIN_ape <- NS_data_DIN %>% 
  mutate(APE_13C = calc_APE(permil_samp = C_per_mil, 
                            isotope = "13C", permil_bg = 0), 
         APE_15N = calc_APE(permil_samp = N_per_mil, 
                            isotope = "15N", permil_bg = 212.3)) %>% 
  mutate(AFE_13C = APE_13C/100, 
         AFE_15N = APE_15N/100)

NS_data_N2_ape <- NS_data_N2 %>% 
  mutate(APE_13C = calc_APE(permil_samp = C_per_mil, 
                            isotope = "13C", permil_bg = 0), 
         APE_15N = calc_APE(permil_samp = N_per_mil, 
                            isotope = "15N", permil_bg = 22.1)) %>% 
  mutate(AFE_13C = APE_13C/100, 
         AFE_15N = APE_15N/100)

NS_data_filt <- rbind(NS_data_DIN_ape, NS_data_N2_ape)

##Plotting 15N-enrichment of trichome-epibiont pairs
cyanos <- c("Afa", "Ma")
styles <- c("epi_Afa", "FL")
jitter_pos_new <- position_jitter(width = 0.1, seed = 123)

axes <-   theme(axis.text = element_text(colour = "black", size = 12), 
                axis.title = element_text(colour = "black", size = 12))

p <- NS_data_filt %>%
  filter(Taxon %in% cyanos) %>%
  filter(Lifestyle %in% styles) %>%
  filter(Epi_Pairs != "None") %>%
  mutate(color_var = ifelse(
    Epi_Pairs %in% unique(Epi_Pairs[Taxon == levels(factor(Taxon))[1] & AFE_15N > 0.04]),
    Epi_Pairs, NA_character_)) %>% 
  ggplot(aes(x = Taxon, y = AFE_15N, shape = Treatment)) +
  geom_point(size = 4, aes(color = color_var)) +
  geom_line(aes(group = Epi_Pairs, color = color_var), linetype = "dotted") +
  scale_shape_manual(values = c(16, 1)) +
  theme_classic() + ggtitle("15N2") + 
  scale_x_discrete(labels = c(expression(italic("D. flos-aquae")), 
                              expression(italic("M. aeruginosa")))) + 
  axes + theme(legend.position = "top") + 
  scale_color_brewer(palette = "Dark2", 
                     na.value = "black", guide = "none") +
  labs(y = expression(""^"15"*"N AFE")) + 
  scale_y_continuous(limits = c(0, 0.60))

gg.gap(plot = p,
       segments = c(0.125, 0.5),
       ylim = c(0, 0.60),
       tick_width = c(0.025, 0.1),
       rel_heights = c(5, 0.3, 1.3))


NS_data_N2_ape %>%
  filter(Taxon %in% cyanos) %>%
  filter(Lifestyle %in% styles) %>%
  filter(Epi_Pairs != "None") %>% 
  mutate(color_var = ifelse(
    Epi_Pairs %in% unique(Epi_Pairs[Taxon == levels(factor(Taxon))[1] & Nnet > 0.30]),
    Epi_Pairs, NA_character_)) %>% 
  ggplot(aes(x = Taxon, y = Nnet*100, shape = Treatment)) +
  geom_point(size = 4, aes(color = color_var)) +
  geom_line(aes(group = Epi_Pairs, color = color_var), linetype = "dotted") +
  scale_shape_manual(values = c(16, 1)) +
  theme_classic() +  
  scale_x_discrete(labels = c(expression(italic("D. flos-aquae")), 
                              expression(italic("M. aeruginosa")))) + 
  axes + theme(legend.position = "top") +
  scale_color_brewer(palette = "Dark2", 
                     na.value = "black", guide = "none") +
  labs(y = expression(N[net]*" (%)")) + 
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80), limits = c(0, 80))

epi_dat <- NS_data_N2_ape %>% filter(Taxon %in% cyanos) %>%
  filter(Lifestyle %in% styles) %>% filter(Epi_Pairs != "None")
##

##15N-enrichment by taxon 
NS_data_N2_ape$Taxon <- factor(NS_data_N2_ape$Taxon, levels = 
                                 c("F", "Ma", "Afa"))
NS_data_N2_ape <- NS_data_N2_ape %>%
  mutate(point_size = ifelse(Lifestyle %in% c("epi_Afa", "epi_Ma"), 2, 1))

#15N AFE N2 baci
NS_data_N2_ape %>% filter(Treatment == "Bacimethrin") %>% 
  ggplot(aes(x = Taxon, y = AFE_15N)) + scale_y_continuous(limits = c(0, .120)) +
  theme_classic() + 
  geom_point(position = jitter_pos, alpha = 0.8, 
             aes(color = Lifestyle, size = point_size)) + geom_boxplot(outlier.shape = NA, fill = NA) +
  scale_color_manual(values = c("darkgreen", "darkred", "darkorange"), 
                     labels = c(expression("Epiphytic, "*italic("D. flos-aquae")), 
                                expression("Epiphytic, "*italic("M. aeruginosa")), 
                                "Free-Living")) +
  labs(y = expression(""^"15"*"N AFE"), 
       title = expression(bold(""^"15"*"N"[2]*" AFE, +Bacimethrin")), 
       x = "", color = "Lifestyle") + 
  ggpubr::stat_compare_means(method = "kruskal.test") +
  scale_size_identity() + 
  stat_summary(fun = mean, geom = "point", shape = 15, size = 3, color = "black") +
  theme(axis.text = element_text(size = 12, colour = "black"), 
        axis.title = element_text(size = 12, colour = "black"), 
        legend.text = element_text(size = 12, colour = "black"), 
        legend.title = element_text(size = 12, colour = "black"), 
        legend.position = "top") +
  scale_x_discrete(labels = c(expression(italic("F")*". LE24H"), 
                              expression(italic("M. aeruginosa")), 
                              expression(italic("D. flos-aquae")))) +
  stat_summary(fun = mean, geom = "line", 
               aes(group = 1, shape = NULL), color = "black", linetype = "dashed")

#

#Nnet N2 baci
NS_data_N2_ape %>% filter(Treatment == "Bacimethrin") %>% 
  ggplot(aes(x = Taxon, y = Nnet*100)) + 
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80), limits = c(0, 80)) +
  theme_classic() + 
  geom_point(position = jitter_pos, alpha = 0.8, 
             aes(color = Lifestyle, size = point_size)) + geom_boxplot(outlier.shape = NA, fill = NA) +
  scale_color_manual(values = c("darkgreen", "darkred", "darkorange"), 
                     labels = c(expression("Epiphytic, "*italic("D. flos-aquae")), 
                                expression("Epiphytic, "*italic("M. aeruginosa")), 
                                "Free-Living")) +
  labs(y = expression(N[net]*" (%)"), 
       title = expression(""^bold("15")*bold("N")[bold("2")]*bold("-N"[net]*", +Bacimethrin")),
       x = "", color = "Lifestyle") + 
  ggpubr::stat_compare_means(method = "kruskal.test") +
  scale_size_identity() + 
  stat_summary(fun = mean, geom = "point", shape = 15, size = 3, color = "black") +
  theme(axis.text = element_text(size = 12, colour = "black"), 
        axis.title = element_text(size = 12, colour = "black"), 
        legend.text = element_text(size = 12, colour = "black"), 
        legend.title = element_text(size = 12, colour = "black"), 
        legend.position = "top") +
  scale_x_discrete(labels = c(expression(italic("F")*". LE24H"), 
                              expression(italic("M. aeruginosa")), 
                              expression(italic("D. flos-aquae")))) +
  stat_summary(fun = mean, geom = "line", 
               aes(group = 1, shape = NULL), color = "black", linetype = "dashed")
#

#15N AFE N2 
NS_data_N2_ape %>% filter(Treatment == "None") %>% 
  ggplot(aes(x = Taxon, y = AFE_15N)) +
  theme_classic() + 
  geom_point(position = jitter_pos, alpha = 0.8, 
             aes(color = Lifestyle, size = point_size)) + geom_boxplot(outlier.shape = NA, fill = NA) +
  scale_color_manual(values = c("darkgreen", "darkred", "darkorange"), 
                     labels = c(expression("Epiphytic, "*italic("D. flos-aquae")), 
                                expression("Epiphytic, "*italic("M. aeruginosa")), 
                                "Free-Living")) +
  labs(y = expression(""^"15"*"N AFE"), 
       title = expression(bold(""^"15"*"N"[2]*" AFE, -Bacimethrin")), 
       x = "", color = "Lifestyle") + scale_y_continuous(limits = c(0, .120)) +
  ggpubr::stat_compare_means(method = "kruskal.test",
                             aes(label = sprintf("Kruskal-Wallis, p = %.1e", 
                                                 ggplot2::after_stat(p)))) +
  scale_size_identity() +
  stat_summary(fun = mean, geom = "point", shape = 15, size = 3, color = "black") +
  theme(axis.text = element_text(size = 12, colour = "black"), 
        axis.title = element_text(size = 12, colour = "black"), 
        legend.text = element_text(size = 12, colour = "black"), 
        legend.title = element_text(size = 12, colour = "black"), 
        legend.position = "top") +
  scale_x_discrete(labels = c(expression(italic("F")*". LE24H"), 
                              expression(italic("M. aeruginosa")), 
                              expression(italic("D. flos-aquae")))) +
  stat_summary(fun = mean, geom = "line", 
               aes(group = 1, shape = NULL), color = "black", linetype = "dashed")
#

#Nnet N2
NS_data_N2_ape %>% filter(Treatment == "None") %>% 
  ggplot(aes(x = Taxon, y = Nnet*100)) + 
  theme_classic() + 
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80), limits = c(0, 80)) +
  geom_point(position = jitter_pos, alpha = 0.8, 
             aes(color = Lifestyle, size = point_size)) + geom_boxplot(outlier.shape = NA, fill = NA) +
  scale_color_manual(values = c("darkgreen", "darkred", "darkorange"), 
                     labels = c(expression("Epiphytic, "*italic("D. flos-aquae")), 
                                expression("Epiphytic, "*italic("M. aeruginosa")), 
                                "Free-Living")) +
  labs(y = expression(N[net]*" (%)"), 
       title = expression(""^bold("15")*bold("N")[bold("2")]*bold("-N"[net]*", -Bacimethrin")), 
       x = "", color = "Lifestyle") + 
  ggpubr::stat_compare_means(method = "kruskal.test",
                             aes(label = sprintf("Kruskal-Wallis, p = %.1e", 
                                                 ggplot2::after_stat(p)))) +
  scale_size_identity() + 
  stat_summary(fun = mean, geom = "point", shape = 15, size = 3, color = "black") +
  theme(axis.text = element_text(size = 12, colour = "black"), 
        axis.title = element_text(size = 12, colour = "black"), 
        legend.text = element_text(size = 12, colour = "black"), 
        legend.title = element_text(size = 12, colour = "black"), 
        legend.position = "top") +
  scale_x_discrete(labels = c(expression(italic("F")*". LE24H"), 
                              expression(italic("M. aeruginosa")), 
                              expression(italic("D. flos-aquae")))) +
  stat_summary(fun = mean, geom = "line", 
               aes(group = 1, shape = NULL), color = "black", linetype = "dashed")
#

NS_data_DIN_ape <- NS_data_DIN_ape %>%
  mutate(point_size = ifelse(Lifestyle %in% c("epi_Afa", "epi_Ma"), 2, 1))
NS_data_DIN_ape$Taxon <- factor(NS_data_DIN_ape$Taxon, levels = 
                                  c("F", "Ma", "Afa"))

#15N AFE DIN baci
NS_data_DIN_ape %>% filter(Treatment == "Bacimethrin") %>% 
  ggplot(aes(x = Taxon, y = AFE_15N)) + 
  theme_classic() + 
  geom_point(position = jitter_pos, alpha = 0.8, 
             aes(color = Lifestyle, size = point_size)) + geom_boxplot(outlier.shape = NA, fill = NA) +
  scale_color_manual(values = c("darkgreen", "darkred", "darkorange"), 
                     labels = c(expression("Epiphytic, "*italic("D. flos-aquae")), 
                                expression("Epiphytic, "*italic("M. aeruginosa")), 
                                "Free-Living")) +
  labs(y = expression(""^"15"*"N AFE"), 
       title = expression(bold(""^"15"*"NH"[4]^"+"*" AFE, +Bacimethrin")), 
       x = "", color = "Lifestyle") + 
  ggpubr::stat_compare_means(method = "kruskal.test",
                             aes(label = sprintf("Kruskal-Wallis, p = %.1e", 
                                                 ggplot2::after_stat(p)))) + 
  scale_size_identity() + 
  stat_summary(fun = mean, geom = "point", shape = 15, size = 3, color = "black") +
  theme(axis.text = element_text(size = 12, colour = "black"), 
        axis.title = element_text(size = 12, colour = "black"), 
        legend.text = element_text(size = 12, colour = "black"), 
        legend.title = element_text(size = 12, colour = "black"), 
        legend.position = "top") +
  scale_x_discrete(labels = c(expression(italic("F")*". LE24H"), 
                              expression(italic("M. aeruginosa")), 
                              expression(italic("D. flos-aquae")))) +
  stat_summary(fun = mean, geom = "line", 
               aes(group = 1, shape = NULL), color = "black", linetype = "dashed")
#

#Nnet DIN baci
NS_data_DIN_ape %>% filter(Treatment == "Bacimethrin") %>% 
  ggplot(aes(x = Taxon, y = Nnet*100)) + 
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80), limits = c(0, 80)) +
  theme_classic() + 
  geom_point(position = jitter_pos, alpha = 0.8, 
             aes(color = Lifestyle, size = point_size)) + geom_boxplot(outlier.shape = NA, fill = NA) +
  scale_color_manual(values = c("darkgreen", "darkred", "darkorange"), 
                     labels = c(expression("Epiphytic, "*italic("D. flos-aquae")), 
                                expression("Epiphytic, "*italic("M. aeruginosa")), 
                                "Free-Living")) +
  labs(y = expression(N[net]*" (%)"), 
       title = expression(""^bold("15")*bold("NH")[bold("4")]^bold("+")*bold("-N"[net]*", +Bacimethrin")), 
       x = "", color = "Lifestyle") + 
  ggpubr::stat_compare_means(method = "kruskal.test",
                             aes(label = sprintf("Kruskal-Wallis, p = %.1e", 
                                                 ggplot2::after_stat(p)))) + 
  scale_size_identity() + 
  stat_summary(fun = mean, geom = "point", shape = 15, size = 3, color = "black") +
  theme(axis.text = element_text(size = 12, colour = "black"), 
        axis.title = element_text(size = 12, colour = "black"), 
        legend.text = element_text(size = 12, colour = "black"), 
        legend.title = element_text(size = 12, colour = "black"), 
        legend.position = "top") +
  scale_x_discrete(labels = c(expression(italic("F")*". LE24H"), 
                              expression(italic("M. aeruginosa")), 
                              expression(italic("D. flos-aquae")))) +
  stat_summary(fun = mean, geom = "line", 
               aes(group = 1, shape = NULL), color = "black", linetype = "dashed")
#

#15N AFE DIN
NS_data_DIN_ape %>% filter(Treatment == "None") %>% ggplot(aes(x = Taxon, y = AFE_15N)) +
  theme_classic() +
  geom_point(position = jitter_pos, alpha = 0.8, 
             aes(color = Lifestyle, size = point_size)) + geom_boxplot(outlier.shape = NA, fill = NA) +
  scale_color_manual(values = c("darkgreen", "darkred", "darkorange"), 
                     labels = c(expression("Epiphytic, "*italic("D. flos-aquae")), 
                                expression("Epiphytic, "*italic("M. aeruginosa")), 
                                "Free-Living")) +
  labs(y = expression(""^"15"*"N AFE"), 
       title = expression(bold(""^"15"*"NH"[4]^"+"*" AFE, -Bacimethrin")), 
       x = "", color = "Lifestyle") + 
  ggpubr::stat_compare_means(method = "kruskal.test",
                             aes(label = sprintf("Kruskal-Wallis, p = %.1e", 
                                                 ggplot2::after_stat(p)))) +
  scale_size_identity() +
  stat_summary(fun = mean, geom = "point", shape = 15, 
               size = 3, color = "black") +
  theme(axis.text = element_text(size = 12, colour = "black"), 
        axis.title = element_text(size = 12, colour = "black"), 
        legend.text = element_text(size = 12, colour = "black"), 
        legend.title = element_text(size = 12, colour = "black"), 
        legend.position = "top") +
  scale_x_discrete(labels = c(expression(italic("F")*". LE24H"), 
                              expression(italic("M. aeruginosa")), 
                              expression(italic("D. flos-aquae")))) +
  stat_summary(fun = mean, geom = "line", 
               aes(group = 1, shape = NULL), color = "black", linetype = "dashed")
#

#Nnet DIN
NS_data_DIN_ape %>% filter(Treatment == "None") %>% 
  ggplot(aes(x = Taxon, y = Nnet*100)) + 
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80), limits = c(0, 80)) +
  theme_classic() + 
  geom_point(position = jitter_pos, alpha = 0.8, 
             aes(color = Lifestyle, size = point_size)) + geom_boxplot(outlier.shape = NA, fill = NA) +
  scale_color_manual(values = c("darkgreen", "darkred", "darkorange"), 
                     labels = c(expression("Epiphytic, "*italic("D. flos-aquae")), 
                                expression("Epiphytic, "*italic("M. aeruginosa")), 
                                "Free-Living")) +
  labs(y = expression(N[net]*" (%)"), 
       title = expression(""^bold("15")*bold("NH")[bold("4")]^bold("+")*bold("-N"[net]*", -Bacimethrin")),
       x = "", color = "Lifestyle") + 
  ggpubr::stat_compare_means(method = "kruskal.test",
                             aes(label = sprintf("Kruskal-Wallis, p = %.1e", 
                                                 ggplot2::after_stat(p)))) +
  scale_size_identity() + 
  stat_summary(fun = mean, geom = "point", shape = 15, size = 3, color = "black") +
  theme(axis.text = element_text(size = 12, colour = "black"), 
        axis.title = element_text(size = 12, colour = "black"), 
        legend.text = element_text(size = 12, colour = "black"), 
        legend.title = element_text(size = 12, colour = "black"), 
        legend.position = "top") +
  scale_x_discrete(labels = c(expression(italic("F")*". LE24H"), 
                              expression(italic("M. aeruginosa")), 
                              expression(italic("D. flos-aquae")))) +
  stat_summary(fun = mean, geom = "line", 
               aes(group = 1, shape = NULL), color = "black", linetype = "dashed")
#
##

##ROI count stacked boxplots
NS_data_filt$N_Substrate <- factor(NS_data_filt$N_Substrate, 
                                   levels = c("Ammonium", "N2"))
NS_data_filt$Taxon <- factor(NS_data_filt$Taxon, 
                             levels = c("Afa", "Ma", "F"))


Ma_dat <- NS_data_filt %>% filter(Taxon == "Ma")
Ma_dat$Lifestyle <- factor(Ma_dat$Lifestyle, 
                           levels = c("epi_Afa", "FL"))

F_dat <- NS_data_filt %>% filter(Taxon == "F")
F_dat$Lifestyle <- factor(F_dat$Lifestyle, 
                          levels = c("epi_Afa", "epi_Ma", "FL"))

Ma_dat %>%
  group_by(N_Substrate, Treatment, Lifestyle) %>%
  summarise(cell_count = n()) %>%
  ggplot(aes(x = interaction(N_Substrate, Treatment), y = cell_count, fill = Lifestyle)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5) +
  theme_classic() +
  scale_fill_manual(values = c("epi_Afa" = "darkgreen", "FL" = "darkorange"),
                    labels = c(expression("Epiphytic, "*italic("D. flos-aquae")), 
                               "Free-Living")) +
  labs(y = "ROI Count", x = "", title = expression(italic("M. aeruginosa")* " ROIs")) + 
  theme(axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"),
        legend.position = "top", 
        legend.text = element_text(size = 12, colour = "black")) + 
  scale_x_discrete(limits = c("Ammonium.None", "Ammonium.Bacimethrin", 
                              "N2.None", "N2.Bacimethrin"), 
                   labels = c(expression(NH[4]^"+"), 
                              expression(NH[4]^"+"*", +Bacimethrin"), 
                              expression(N[2]), 
                              expression(N[2]*", +Bacimethrin")))

F_dat %>%
  group_by(N_Substrate, Treatment, Lifestyle) %>%
  summarise(cell_count = n()) %>%
  ggplot(aes(x = interaction(N_Substrate, Treatment), y = cell_count, fill = Lifestyle)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5) +
  theme_classic() +
  scale_fill_manual(values = c("epi_Afa" = "darkgreen", "epi_Ma" = "darkred", "FL" = "darkorange"),
                    labels = c(expression("Epiphytic, "*italic("D. flos-aquae")), 
                               expression("Epiphytic, "*italic("M. aeruginosa")), 
                               "Free-Living")) +
  labs(y = "ROI Count", x = "", title = expression(italic("F")* ". LE24H ROIs")) +
  theme(axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"),
        legend.position = "top", 
        legend.text = element_text(size = 12, colour = "black")) + 
  scale_x_discrete(limits = c("Ammonium.None", "Ammonium.Bacimethrin", 
                              "N2.None", "N2.Bacimethrin"), 
                   labels = c(expression(NH[4]^"+"), 
                              expression(NH[4]^"+"*", +Bacimethrin"), 
                              expression(N[2]), 
                              expression(N[2]*", +Bacimethrin")))
##

##Stats
NS_data_N2_ape$Filter <- factor(NS_data_N2_ape$Filter, levels = c("F9", "F10", "F11", 
                                                                  "F12", "F13", "F14"))
NS_data_DIN_ape$Filter <- factor(NS_data_DIN_ape$Filter, levels = c("F9", "F10", "F11", 
                                                                    "F12", "F13", "F14"))

NS_data_N2_ape$Taxon <- factor(NS_data_N2_ape$Taxon, levels = 
                                 c("F", "Ma", "Afa"))
NS_data_DIN_ape$Taxon <- factor(NS_data_DIN_ape$Taxon, levels = 
                                  c("F", "Ma", "Afa"))

Flavo_NS_data_N2_ape <- NS_data_N2_ape %>% filter(Taxon == "F")
Ma_NS_data_N2_ape <- NS_data_N2_ape %>% filter(Taxon == "Ma")
Afa_NS_data_N2_ape <- NS_data_N2_ape %>% filter(Taxon == "Afa")

Flavo_NS_data_DIN_ape <- NS_data_DIN_ape %>% filter(Taxon == "F")
Ma_NS_data_DIN_ape <- NS_data_DIN_ape %>% filter(Taxon == "Ma")
Afa_NS_data_DIN_ape <- NS_data_DIN_ape %>% filter(Taxon == "Afa")

Flavo_NS_data_N2_ape_baci <- Flavo_NS_data_N2_ape %>% filter(Treatment == "Bacimethrin")
Flavo_NS_data_N2_ape_Nobaci <- Flavo_NS_data_N2_ape %>% filter(Treatment == "None")

mean(Flavo_NS_data_N2_ape_baci$Nnet)*100
mean(Flavo_NS_data_N2_ape_Nobaci$Nnet)*100

Ma_NS_data_N2_ape_baci <- Ma_NS_data_N2_ape %>% filter(Treatment == "Bacimethrin")
Ma_NS_data_N2_ape_Nobaci <- Ma_NS_data_N2_ape %>% filter(Treatment == "None")

mean(Ma_NS_data_N2_ape_baci$Nnet)*100
mean(Ma_NS_data_N2_ape_Nobaci$Nnet)*100


Afa_NS_data_N2_ape_baci <- Afa_NS_data_N2_ape %>% filter(Treatment == "Bacimethrin")
Afa_NS_data_N2_ape_Nobaci <- Afa_NS_data_N2_ape %>% filter(Treatment == "None")

mean(Afa_NS_data_N2_ape_baci$Nnet)*100
mean(Afa_NS_data_N2_ape_Nobaci$Nnet)*100


Afa_NS_data_N2_ape %>% wilcox_test(Nnet ~ Treatment, 
                                   comparisons = list("Bacimethrin", "None"), 
                                   ref.group = "None", paired = F)
Ma_NS_data_N2_ape %>% wilcox_test(Nnet ~ Treatment, 
                                  comparisons = list("Bacimethrin", "None"), 
                                  ref.group = "None", paired = F)
Flavo_NS_data_N2_ape %>% wilcox_test(Nnet ~ Treatment, 
                                     comparisons = list("Bacimethrin", "None"), 
                                     ref.group = "None", paired = F)



Flavo_NS_data_DIN_ape_baci <- Flavo_NS_data_DIN_ape %>% filter(Treatment == "Bacimethrin")
Flavo_NS_data_DIN_ape_Nobaci <- Flavo_NS_data_DIN_ape %>% filter(Treatment == "None")

mean(Flavo_NS_data_DIN_ape_baci$Nnet)*100
mean(Flavo_NS_data_DIN_ape_Nobaci$Nnet)*100


Ma_NS_data_DIN_ape_baci <- Ma_NS_data_DIN_ape %>% filter(Treatment == "Bacimethrin")
Ma_NS_data_DIN_ape_Nobaci <- Ma_NS_data_DIN_ape %>% filter(Treatment == "None")

mean(Ma_NS_data_DIN_ape_baci$Nnet)*100
mean(Ma_NS_data_DIN_ape_Nobaci$Nnet)*100


Afa_NS_data_DIN_ape_baci <- Afa_NS_data_DIN_ape %>% filter(Treatment == "Bacimethrin")
Afa_NS_data_DIN_ape_Nobaci <- Afa_NS_data_DIN_ape %>% filter(Treatment == "None")

mean(Afa_NS_data_DIN_ape_baci$Nnet)*100
mean(Afa_NS_data_DIN_ape_Nobaci$Nnet)*100


Afa_NS_data_DIN_ape %>% wilcox_test(Nnet ~ Treatment, 
                                    comparisons = list("Bacimethrin", "None"), 
                                    ref.group = "None", paired = F)
Ma_NS_data_DIN_ape %>% wilcox_test(Nnet ~ Treatment, 
                                   comparisons = list("Bacimethrin", "None"), 
                                   ref.group = "None", paired = F)
Flavo_NS_data_DIN_ape %>% wilcox_test(Nnet ~ Treatment, 
                                      comparisons = list("Bacimethrin", "None"), 
                                      ref.group = "None", paired = F)


F_epi_Afa_N2 <- Flavo_NS_data_N2_ape %>% filter(Lifestyle == "epi_Afa")
F_epi_Ma_N2 <- Flavo_NS_data_N2_ape %>% filter(Lifestyle == "epi_Ma")
F_FL_N2 <- Flavo_NS_data_N2_ape %>% filter(Lifestyle == "FL")

F_epi_Afa_DIN <- Flavo_NS_data_DIN_ape %>% filter(Lifestyle == "epi_Afa")
F_epi_Ma_DIN <- Flavo_NS_data_DIN_ape %>% filter(Lifestyle == "epi_Ma")
F_FL_DIN <- Flavo_NS_data_DIN_ape %>% filter(Lifestyle == "FL")


Ma_epi_Afa_N2 <- Ma_NS_data_N2_ape %>% filter(Lifestyle == "epi_Afa")
Ma_FL_N2 <- Ma_NS_data_N2_ape %>% filter(Lifestyle == "FL")

Ma_epi_Afa_DIN <- Ma_NS_data_DIN_ape %>% filter(Lifestyle == "epi_Afa")
Ma_FL_DIN <- Ma_NS_data_DIN_ape %>% filter(Lifestyle == "FL")

mean(F_epi_Afa_N2$Nnet)*100
mean(F_epi_Ma_N2$Nnet)*100
mean(F_FL_N2$Nnet)*100
mean(F_epi_Afa_DIN$Nnet)*100
mean(F_epi_Ma_DIN$Nnet)*100
mean(F_FL_DIN$Nnet)*100

mean(Ma_epi_Afa_N2$Nnet)*100
mean(Ma_FL_N2$Nnet)*100
mean(Ma_epi_Afa_DIN$Nnet)*100
mean(Ma_FL_DIN$Nnet)*100


NS_data_N2_ape %>% filter(Taxon == "F") %>% 
  wilcox_test(Nnet ~ Lifestyle, 
              comparisons = list(c("epi_Afa", "FL")), 
              ref.group = "FL", paired = F)

NS_data_N2_ape %>% filter(Taxon == "F") %>% 
  wilcox_test(Nnet ~ Lifestyle, 
              comparisons = list(c("epi_Ma", "FL")), 
              ref.group = "FL", paired = F)

NS_data_DIN_ape %>% filter(Taxon == "F") %>% 
  wilcox_test(Nnet ~ Lifestyle, 
              comparisons = list(c("epi_Afa", "FL")), 
              ref.group = "FL", paired = F)

NS_data_DIN_ape %>% filter(Taxon == "F") %>% 
  wilcox_test(Nnet ~ Lifestyle, 
              comparisons = list(c("epi_Ma", "FL")), 
              ref.group = "FL", paired = F)

NS_data_N2_ape %>% filter(Taxon == "Ma") %>% 
  wilcox_test(Nnet ~ Lifestyle, 
              comparisons = list(c("epi_Afa", "FL")), 
              ref.group = "FL", paired = F)
##







