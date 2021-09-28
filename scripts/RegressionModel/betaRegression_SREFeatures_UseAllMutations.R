rm(list = ls(all.names = TRUE))

setwd("~/Documents/LaederachLab/Papers/MAPT_splicing_structure/Figure4/")

library(betareg)
library(ggplot2)
library(cowplot)
library(moments)
library(reshape)

# Function to get standard normal or zscore
zscore <- function(x){
  return((x-mean(x))/sd(x))
}

muts_LowConfidence = c("Mut8","Mut10","Mut16","Mut31","Mut32","Mut36","Mut37","Mut40","Mut43","Mut45","Mut48","Mut51","Mut54","Mut56")

mutIDs_position = c("Mut1","Mut2","Mut4","Mut5","Mut6","Mut7","Mut9","Mut39","Mut11","Mut12","Mut13","Mut14","Mut15","Mut17","Mut41","Mut18","Mut19","Mut42","Mut44","Mut20",
                    "Mut21","Mut46","Mut47","Mut22","Mut23","Mut24","Mut25","Mut26","Mut49","Mut27","Mut50","Mut28","Mut52","Mut29","Mut53","Mut30","Mut33","Mut55","Mut57","Mut58","Mut59","Mut34","Mut35","Mut38")

allMuts_nonDoubtfulOnes <- c(mutIDs_position,"Mut60","Mut61","Mut63")

PSI_singleMuts <- read.table("data/MAPT_Model_SingleMutations_ExperimentalPSI.tsv",sep="\t",header=T)
PSI_doubleMuts <- read.table("data/MAPT_Model_DoubleMutations_ExperimentalPSI.tsv",sep="\t",header=T)
PSI_doubleMuts_NoWT <- PSI_doubleMuts[PSI_doubleMuts$MutID!="WT",]
PSI_testedSNPs <- read.table("data/MAPT_Model_TestedSNPs_ExperimentalPSI_AllWTPSI.tsv",sep="\t",header=T)


PSI_allMuts <- rbind(PSI_singleMuts,PSI_doubleMuts_NoWT,PSI_testedSNPs)
#allMuts_DoubtfulOnes <- PSI_allMuts$MutID[!(PSI_allMuts$MutID %in% c("WT",allMuts_nonDoubtfulOnes))]


# SREs
data_SS_SingleMuts <- read.table("data/MAPT_ModelFeature_SingleMutations_TrainingSet_WTvsMUT_DiffIn5pSpliceSitesScores.tsv",sep="\t",header=T)
data_SS_DoubleMuts <- read.table("data/MAPT_ModelFeature_DoubleMutations_TrainingSet_WTvsMUT_DiffIn5pSpliceSitesScores.tsv",sep="\t",header=T)
data_SS_DoubleMuts_NoWT <- data_SS_DoubleMuts[data_SS_DoubleMuts$MutID!="WT",]
data_SS_TestedSNPs <- read.table("data/MAPT_ModelFeature_TestedSNPs_TrainingSet_WTvsMUT_DiffIn5pSpliceSitesScores.tsv",sep="\t",header=T)

data_SS <- rbind(data_SS_SingleMuts,data_SS_DoubleMuts_NoWT,data_SS_TestedSNPs)

data_SRE_SingleMuts <- read.table("data/MAPT_ModelFeature_SingleMutations_TrainingSet_WTvsMUT_DiffStrengthOfSREs_MeanOfClustersPerMotifType.tsv",sep="\t",header = T)
data_SRE_DoubleMuts <- read.table("data/MAPT_ModelFeature_DoubleMutations_TrainingSet_WTvsMUT_DiffStrengthOfSREs_MeanOfClustersPerMotifType.tsv",sep="\t",header = T)
data_SRE_DoubleMuts_NoWT <- data_SRE_DoubleMuts[data_SRE_DoubleMuts$MutID!="WT",]
data_SRE_TestedSNPs <- read.table("data/MAPT_ModelFeature_TestedSNPs_TrainingSet_WTvsMUT_DiffStrengthOfSREs_MeanOfClustersPerMotifType.tsv",sep="\t",header = T)

data_SRE <- rbind(data_SRE_SingleMuts,data_SRE_DoubleMuts_NoWT,data_SRE_TestedSNPs)

data_toModelAndTest_SRE <- data_SRE
data_toModelAndTest_SRE$SS <- data_SS$Strength
data_toModelAndTest_SRE$Enhancer = data_SRE$ESE + data_SRE$ISE
data_toModelAndTest_SRE$Silencer = data_SRE$ESS + data_SRE$ISS
data_toModelAndTest_SRE$ExperimentalPSI <- PSI_allMuts$PSI
data_toModelAndTest_SRE$Label <- PSI_allMuts$Label
data_toModelAndTest_SRE$isSynonymous <- PSI_allMuts$isSynonymous
data_toModelAndTest_SRE$isNonSynonymous <- PSI_allMuts$isNonSynonymous
data_toModelAndTest_SRE$isExonic <- PSI_allMuts$isExonic
data_toModelAndTest_SRE$isIntronic <- PSI_allMuts$isIntronic
data_toModelAndTest_SRE$isDouble <- PSI_allMuts$isDouble


data_toModelAndTest <- data_toModelAndTest_SRE



data_toModelAndTest <- data_toModelAndTest[data_toModelAndTest$MutID %in% allMuts_nonDoubtfulOnes, ]
n=nrow(data_toModelAndTest)
data_toModelAndTest$PSI = (data_toModelAndTest$ExperimentalPSI*(n-1)+0.5)/n

beta_model <- betareg(PSI~SS, data = data_toModelAndTest)
beta_model <- betareg(PSI~SS+Enhancer+Silencer, data = data_toModelAndTest)
print(summary(beta_model))

data_toModelAndTest$PredictedPSI <- predict(beta_model,data_toModelAndTest)

p<-ggplot(data_toModelAndTest, aes(x=PSI, y=PredictedPSI)) +
  geom_point(aes(color=Label,shape=Label),size=5) +
  #scale_colour_brewer(palette = "Set1")+
  scale_color_manual(values=c('#A2D5C6','#3C53DD','#F52549','#A13EDB','#000000')) +
  scale_shape_manual(values=c(15, 17, 18,16,7))+
  #geom_abline(intercept = 0, slope=1, linetype = "dashed") +
  geom_smooth(method='lm',formula= y~x,linetype="dashed")+
  xlim(0,1) +
  ylim(0,1) +
  xlab("Experimental PSI Exon 10") +
  ylab("Predicted PSI Exon 10") + 
  theme_classic() +
  theme(text = element_text(size=25),
        axis.text.x = element_text(size=20),axis.text.y = element_text(size=20),
        panel.background = element_rect(colour = "black"),legend.position = "None")

p

#ggsave(p,file="figures/PredictedVsActualPSI_44Muts_JustSpliceSitestrength.svg",width = 10, height = 10, bg = "white")
#ggsave(p,file="figures/PredictedVsActualPSI_44Muts_SREstrength.svg",width = 10, height = 10, bg = "white")

cor.test(data_toModelAndTest$PSI,data_toModelAndTest$PredictedPSI)

for (label in c("NonSynonymous","Synonymous","Intronic","Double")){
  data_toModelAndTest_label = data_toModelAndTest[data_toModelAndTest$Label==label,]
  print(label)
  print(cor.test(data_toModelAndTest_label$PSI,data_toModelAndTest_label$PredictedPSI))
}

set.seed(20)

corr_coeff_all_1 = list()
corr_coeff_all_2 = list()


for (i in 1:10){
  
  #train_set = sample(allMuts_nonDoubtfulOnes,size=15,replace=FALSE)
  #test_set = sample(data_toModelAndTest$MutIDs,size = 36,replace = FALSE)
  train_set = sample(data_toModelAndTest$MutID,size = 35,replace = FALSE)
  
  if (!("WT" %in% train_set)){
    train_set = c(train_set,"WT")
  }
  
  data_toModel <- data_toModelAndTest[data_toModelAndTest$MutID %in% train_set,]
  
  beta_model_1 <- betareg(PSI~SS+Enhancer+Silencer, data = data_toModel)
  beta_model_2 <- betareg(PSI~SS, data = data_toModel)
  
  #data_toTest <- data_toModelAndTest[!(data_toModelAndTest$MutID %in% train_set),]
  data_toTest <- data_toModelAndTest[(data_toModelAndTest$MutID %in% train_set),]
  
  data_toTest$PredictedPSI_1 <- predict(beta_model_1,data_toTest)
  data_toTest$PredictedPSI_2 <- predict(beta_model_2,data_toTest)
  
  corr_coeff_1 <- c()
  corr_coeff_2 <- c()
  
  for (type_label in c("NonSynonymous","Synonymous","Intronic","Double")){
    #for (type_label in c("Exonic","Intronic","Double")){
    #for (type_label in c("Synonymous","Intronic")){
    #print(type_label)
    data_toTest_Label = data_toTest[data_toTest$Label==type_label,]
    if (nrow(data_toTest_Label) < 3){
      corr_coeff_1 <- c(corr_coeff_1,NA)
      corr_coeff_2 <- c(corr_coeff_2,NA)
    }
    else{
      cc_test = cor.test(data_toTest_Label$PSI,data_toTest_Label$PredictedPSI_1)  
      corr_coeff_1 <- c(corr_coeff_1,cc_test$estimate)
      cc_test = cor.test(data_toTest_Label$PSI,data_toTest_Label$PredictedPSI_2)  
      corr_coeff_2 <- c(corr_coeff_2,cc_test$estimate)
    }
  }
  
  corr_coeff_all_1[[i]] <- corr_coeff_1
  corr_coeff_all_2[[i]] <- corr_coeff_2
  
}

# Plot overall correlation plot for each complex for all mutations
corr_coeff_train_All_df <- data.frame(matrix(unlist(corr_coeff_all_2), nrow=length(corr_coeff_all_1), byrow=TRUE))
colnames(corr_coeff_train_All_df) <- c("Exonic: Non-Synonymous","Exonic: Synonymous","Intronic","Compensatory")
#summary_corr_coeff_train = data.frame(typeMuts=c("Non-Synonymous","Synonymous","Intronic","Compensatory"),meanVals=apply(corr_coeff_train_All_df,2,mean),sdVals=apply(corr_coeff_train_All_df,2,sd))
corr_coeff_train_All_df_melted = melt(corr_coeff_train_All_df)

p<-ggplot(corr_coeff_train_All_df_melted, aes(x=factor(variable,levels = c("Exonic: Non-Synonymous","Exonic: Synonymous","Intronic","Compensatory")), y=value,fill=variable)) + 
  geom_violin()+
  #geom_bar(position=position_dodge(), stat="identity") +
  #geom_errorbar(aes(ymin=meanVals-sdVals, ymax=meanVals+sdVals),
  #              width=.2,                    # Width of the error bars
  #              position=position_dodge(.9)) +
  #scale_fill_brewer(palette="Set1")+
  stat_summary(fun=median, geom="point", size=1.5, color="black")+
  scale_fill_manual(values=c('#F52549',"#A13EDB",'#3C53DD',"#a2d5c6"))+
  xlab("")+
  ylab("Correlation Coefficients")+
  #scale_x_discrete(labels=c("Pre-B complex","B complex","Pre-Bact complex","Bact complex"))+
  coord_cartesian(ylim = c(-0.4,1))+
  theme_classic() +
  theme(text = element_text(size=25),
        axis.text.x = element_text(size=15),axis.text.y = element_text(size=20),
        panel.background = element_rect(colour = "black"),legend.position = "None")
p

#ggsave(p,file="figures/ViolinPlots_SREstrength_47Mutations_TrainTest.svg",width = 10, height = 10, bg = "white")
#ggsave(p,file="figures/ViolinPlots_SplicesiteStrength_47Mutations_TrainTest.svg",width = 10, height = 10, bg = "white")