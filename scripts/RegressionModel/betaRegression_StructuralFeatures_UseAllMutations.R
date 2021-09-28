rm(list = ls(all.names = TRUE))

setwd("~/Documents/LaederachLab/Papers/MAPT_splicing_structure/Figure3/")

library(betareg)
library(ggplot2)
library(cowplot)
library(moments)
library(reshape)


# Function to get standard normal or zscore
zscore <- function(x){
  return((x-mean(x))/sd(x))
}

##### Mutations only in synonymous, intronic and double
singleMuts_level_order_nonDoubtfulOnes = c("Mut17","Mut22","Mut33","Mut34","Mut35","Mut38",
                                           "WT","Mut15","Mut23","Mut24",
                                           "Mut5","Mut9","Mut12","Mut20","Mut21","Mut25","Mut26","Mut27","Mut28","Mut29","Mut30"
)
doubleMuts_level_order_nonDoubtfulOnes <- c()



muts_LowConfidence = c("Mut8","Mut10","Mut16","Mut31","Mut32","Mut36","Mut37","Mut40","Mut43","Mut45","Mut48","Mut51","Mut54","Mut56")

allMuts_nonDoubtfulOnes = c(singleMuts_level_order_nonDoubtfulOnes,doubleMuts_level_order_nonDoubtfulOnes)

############################################################################################################################################3

PSI_singleMuts <- read.table("data/MAPT_Model_SingleMutations_ExperimentalPSI.tsv",sep="\t",header=T)
PSI_doubleMuts <- read.table("data/MAPT_Model_DoubleMutations_ExperimentalPSI.tsv",sep="\t",header=T)
PSI_doubleMuts_NoWT <- PSI_doubleMuts[PSI_doubleMuts$MutID!="WT",] 

PSI_allMuts <- rbind(PSI_singleMuts,PSI_doubleMuts_NoWT)

labelsForMuts = c("WT","Non-Synonymous","Non-Synonymous","Deletion","Non-Synonymous","Synonymous","Non-Synonymous","Non-Synonymous","Non-Synonymous","Synonymous","Deletion","Non-Synonymous","Synonymous","Non-Synonymous","Non-Synonymous","Synonymous","Non-Synonymous","Synonymous","Non-Synonymous","Non-Synonymous","Synonymous","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double")
#labelsForMuts = c("WT","Exonic","Exonic","Exonic","Exonic","Exonic","Exonic","Exonic","Exonic","Exonic","Exonic","Exonic","Exonic","Exonic","Exonic","Exonic","Exonic","Exonic","Exonic","Exonic","Exonic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double")

all_complexes = c("PreBComplex_6qx9","BComplex_5o9z","PreBactComplex_7abf","BactComplex_5z56")

corr_coeff <- c()

set.seed(20)

corr_coeff_train_df = data.frame(complex=all_complexes)
col_names_df = c("Complex")

corr_coeff_train <- list()
corr_coeff_train_JustComplex_list <- list()

# Bootstrapping 
for (i in 1:10){
  
  # exclude 30% of mutations
  test_set = sample(allMuts_nonDoubtfulOnes,size = 4,replace = FALSE)
  train_set = allMuts_nonDoubtfulOnes[!(allMuts_nonDoubtfulOnes %in% test_set)]
  #train_set = allMuts_nonDoubtfulOnes
  
  if (!("WT" %in% train_set)){
    train_set = c(train_set,"WT")
  }
  
  corr_coeff_complexes = list()
  corr_coeff_train_JustComplex = c()
  
  for (complex in all_complexes){
    dataFile = paste("data/AllMuts_Unfold_LengthRNAinSpliceosome_",complex,".tsv",sep="")
    data_ToRead <- read.table(dataFile,sep="\t",header = T)

    mean_Energy <- apply(data_ToRead,2,mean)
    sd_Energy <- apply(data_ToRead,2,sd)
    skew_Energy <- apply(data_ToRead,2,skewness)
    kurtosis_Energy <- apply(data_ToRead,2,kurtosis)


    data_toModelAndTest <- data.frame(MutIDs=PSI_allMuts$MutID,strucMean=mean_Energy,strucSD=sd_Energy,strucSkew=skew_Energy,strucKurt=kurtosis_Energy,ExperimentalPSI=PSI_allMuts$PSI,Labels=labelsForMuts)
    
    data_toModelAndTest <- data_toModelAndTest[data_toModelAndTest$MutIDs %in% allMuts_nonDoubtfulOnes, ]
    
    n=nrow(data_toModelAndTest)
    data_toModelAndTest$PSI = (data_toModelAndTest$ExperimentalPSI*(n-1)+0.5)/n

    data_toModel <- data_toModelAndTest[data_toModelAndTest$MutIDs %in% train_set,]

    beta_model <- betareg(PSI~strucMean+strucSD+strucSkew+strucKurt, data = data_toModel)

    #beta_model <- betareg(PSI~strucMean, data = data_toModel)

  
    #print(summary(beta_model))
    cc_train = cor.test(data_toModel$ExperimentalPSI,predict(beta_model,data_toModel))
  
    corr_coeff_train_JustComplex = c(corr_coeff_train_JustComplex,cc_train$estimate)

    data_toModel$PredictedPSI <- predict(beta_model,data_toModel)
  
    corr_coeff <- c()
    
    #for (type_label in c("Non-Synonymous","Synonymous","Intronic","Double")){
    #for (type_label in c("Exonic","Intronic","Double")){
    for (type_label in c("Synonymous","Intronic")){
      #print(type_label)
      data_toModel_type = data_toModel[data_toModel$Labels==type_label,]
      cc_test = cor.test(data_toModel_type$ExperimentalPSI,data_toModel_type$PredictedPSI)
      corr_coeff <- c(corr_coeff,cc_test$estimate)
    }
    
    corr_coeff_complexes[[complex]] <- corr_coeff
  
  }

  corr_coeff_train[[i]] <- corr_coeff_complexes
  corr_coeff_train_JustComplex_list[[i]] <- corr_coeff_train_JustComplex
}

# Plot overall correlation plot for each complex for all mutations
corr_coeff_train_JustComplex_df <- data.frame(matrix(unlist(corr_coeff_train_JustComplex_list), nrow=length(corr_coeff_train_JustComplex_list), byrow=TRUE))
colnames(corr_coeff_train_JustComplex_df) <- all_complexes
summary_corr_coeff_train = data.frame(Complex=all_complexes,meanVals=apply(corr_coeff_train_JustComplex_df,2,mean),sdVals=apply(corr_coeff_train_JustComplex_df,2,sd))



# Box plot of cc 
corr_coeff_train_JustComplex_df_melted = melt(corr_coeff_train_JustComplex_df)
p <- ggplot(corr_coeff_train_JustComplex_df_melted,aes(x=variable,y=value,fill=variable))+
  geom_boxplot()+
  scale_fill_manual(values=c("#edca82","#097770","#e0cdbe","#a9c0a6"))+
  xlab("")+
  ylab("Correlation Coefficients")+
  scale_x_discrete(labels=c("Pre-B complex","B complex","Pre-Bact complex","Bact complex"))+
  coord_cartesian(ylim = c(0.65,1))+
  theme_classic() +
  theme(text = element_text(size=25),
        axis.text.x = element_text(size=20),axis.text.y = element_text(size=20),
        panel.background = element_rect(colour = "black"),legend.position = "None")
p
#ggsave(p,file="figures/Boxplots_CorrelationCoefficients_AllMutations_JustIntronicAndSynonymous.svg",width = 10, height = 5, bg = "white")


wilcox.test(corr_coeff_train_JustComplex_df$PreBComplex_6qx9,corr_coeff_train_JustComplex_df$BactComplex_5z56)
wilcox.test(corr_coeff_train_JustComplex_df$BComplex_5o9z,corr_coeff_train_JustComplex_df$BactComplex_5z56)
wilcox.test(corr_coeff_train_JustComplex_df$PreBactComplex_7abf,corr_coeff_train_JustComplex_df$BactComplex_5z56)





##################################################################################################################
# Break down by intronic muts and synonymous muts

corr_coeff_train_df <- data.frame(matrix(unlist(corr_coeff_train[[1]]), nrow=length(corr_coeff_train[[1]]), byrow=TRUE))
#colnames(corr_coeff_train_df) <- c("Non-Synonymous","Synonymous","Intronic","Double")
#colnames(corr_coeff_train_df) <- c("Exonic","Intronic","Double")
colnames(corr_coeff_train_df) <- c("Synonymous","Intronic")
corr_coeff_train_df$Labels <- all_complexes
corr_coeff_train_df$SampleNum <- 1
corr_coeff_train_df_allSamples <- melt(corr_coeff_train_df,id.vars=c("Labels","SampleNum"))


for (samp in 2:10){
  corr_coeff_train_df <- data.frame(matrix(unlist(corr_coeff_train[[samp]]), nrow=length(corr_coeff_train[[samp]]), byrow=TRUE))
  #colnames(corr_coeff_train_df) <- c("Non-Synonymous","Synonymous","Intronic","Double")
  #colnames(corr_coeff_train_df) <- c("Exonic","Intronic","Double")
  colnames(corr_coeff_train_df) <- c("Synonymous","Intronic")
  corr_coeff_train_df$Labels <- all_complexes
  corr_coeff_train_df$SampleNum <- samp
  corr_coeff_train_df_melted <- melt(corr_coeff_train_df,id.vars=c("Labels","SampleNum"))
  corr_coeff_train_df_allSamples <- rbind(corr_coeff_train_df_allSamples,corr_coeff_train_df_melted)
}
  
  


corr_complex_Breakdown <- corr_coeff_train_df_allSamples[corr_coeff_train_df_allSamples$Labels=="BactComplex_5z56",]
p <- ggplot(corr_complex_Breakdown,aes(x=variable,y=value,fill=variable))+
  geom_violin()+
  stat_summary(fun=median, geom="point", size=1.5, color="black")+
  scale_fill_manual(values=c("#A13EDB",'#3C53DD'))+
  ylab("Correlation Coefficients")+
  #geom_bar(stat="identity", color="black", position=position_dodge())+
  #geom_violin()+
  coord_cartesian(y=c(0,1))+
  xlab("")+
  ylab("Correlation coefficient")+
  theme_classic() +
  theme(text = element_text(size=25),
        axis.text.x = element_text(size=20),axis.text.y = element_text(size=20),
        panel.background = element_rect(colour = "black"),axis.title.x=element_blank(),legend.position = "None")+
  scale_x_discrete(labels=c("Exonic: Synonymous","Intronic"))

p

#ggsave(p,file="figures/ViolinPlots_BactComplex_5z56_AllMutations_TrainTest_NonSynonymousAndDouble.svg",width = 10, height = 10, bg = "white")



##################################################################################################################
# Train on intronic and synonymous muts
# Test on everything else
dataFile = "data/AllMuts_Unfold_LengthRNAinSpliceosome_PreBComplex_6qx9.tsv"
dataFile = "data/AllMuts_Unfold_LengthRNAinSpliceosome_BComplex_5o9z.tsv"
dataFile = "data/AllMuts_Unfold_LengthRNAinSpliceosome_PreBactComplex_7abf.tsv"
dataFile = "data/AllMuts_Unfold_LengthRNAinSpliceosome_BactComplex_5z56.tsv"


data_ToRead <- read.table(dataFile,sep="\t",header = T)

mean_Energy <- apply(data_ToRead,2,mean)
sd_Energy <- apply(data_ToRead,2,sd)
skew_Energy <- apply(data_ToRead,2,skewness)
kurtosis_Energy <- apply(data_ToRead,2,kurtosis)

labelsForMuts = c("WT","Non-Synonymous","Non-Synonymous","Deletion","Non-Synonymous","Synonymous","Non-Synonymous","Non-Synonymous","Non-Synonymous","Synonymous","Deletion","Non-Synonymous","Synonymous","Non-Synonymous","Non-Synonymous","Synonymous","Non-Synonymous","Synonymous","Non-Synonymous","Non-Synonymous","Synonymous","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Intronic","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double","Double")

data_toModelAndTest <- data.frame(MutIDs=PSI_allMuts$MutID,strucMean=zscore(mean_Energy-mean_Energy[1]),strucSD=zscore(sd_Energy-sd_Energy[1]),strucSkew=zscore(skew_Energy-skew_Energy[1]),strucKurt=zscore(kurtosis_Energy-kurtosis_Energy[1]),ExperimentalPSI=PSI_allMuts$PSI,Labels=labelsForMuts)

data_toModelAndTest <- data_toModelAndTest[!(data_toModelAndTest$MutIDs %in% c("Mut3",muts_LowConfidence)), ]
#data_toModelAndTest <- data_toModelAndTest[data_toModelAndTest$MutIDs %in% allMuts_nonDoubtfulOnes, ]


n=nrow(data_toModelAndTest)
data_toModelAndTest$PSI = (data_toModelAndTest$ExperimentalPSI*(n-1)+0.5)/n
data_toModel <- data_toModelAndTest[data_toModelAndTest$MutIDs %in% allMuts_nonDoubtfulOnes, ]

beta_model <- betareg(PSI~strucMean+strucSD+strucSkew+strucKurt, data = data_toModel)
#beta_model <- betareg(PSI~strucMean, data = data_toModel)

summary(beta_model)
data_toPlot <- data_toModel
data_toPlot$PredictedPSI <- predict(beta_model,data_toModel)

p<-ggplot(data_toPlot, aes(x=PSI, y=PredictedPSI)) +
  geom_point(aes(color=Labels,shape=Labels),size=5) +
  #scale_colour_brewer(palette = "Set1")+
  scale_color_manual(values=c('#3C53DD',"#A13EDB","#000000")) +
  scale_shape_manual(values=c(17,16,10))+
  #geom_abline(intercept = 0, slope=1, linetype = "dashed") +
  geom_smooth(method='lm',formula= y~x,linetype="dashed")+
  xlim(0,1) +
  ylim(0,1) +
  xlab("Experimental PSI Exon 10") +
  ylab("Predicted PSI Exon 10") + 
  theme_classic() +
  #theme(text = element_text(size=25),panel.background = element_rect(colour = "black"),
  #      axis.text.x = element_text(size=20),axis.text.y = element_text(size=20),legend.position = "None")
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        legend.position = "None",panel.background = element_rect(colour = "black"))
p

ggsave(p,file="figures/PredictedVsActualPSI_AllMutations_NonDoubtfulOnes_BactComplex_5z56_JustIntronicAndSynonymous.svg",width = 10, height = 10, bg = "white")
#ggsave(p,file="figures/PredictedVsActualPSI_AllMutations_NonDoubtfulOnes_PreBComplex_6qx9_JustIntronicAndSynonymous.svg",width = 10, height = 10, bg = "white")
#ggsave(p,file="figures/PredictedVsActualPSI_AllMutations_NonDoubtfulOnes_BComplex_5o9z_JustIntronicAndSynonymous.svg",width = 10, height = 10, bg = "white")
#ggsave(p,file="figures/PredictedVsActualPSI_AllMutations_NonDoubtfulOnes_PreBactComplex_7abf_JustIntronicAndSynonymous.svg",width = 10, height = 10, bg = "white")


cor.test(data_toPlot$PSI,data_toPlot$PredictedPSI)


dataToTest <- data_toModelAndTest[!(data_toModelAndTest$MutIDs %in% c(allMuts_nonDoubtfulOnes)), ]
dataToTest$PredictedPSI <- predict(beta_model,dataToTest)

p<-ggplot(dataToTest, aes(x=PSI, y=PredictedPSI,color=Labels,shape=Labels)) +
  geom_point(size=5) +
  #scale_colour_brewer(palette = "Set1")+
  scale_color_manual(values=c("#a2d5c6",'#F52549')) +
  scale_shape_manual(values=c(15,18))+
  geom_abline(intercept = 0, slope=1, linetype = "dashed") +
  xlim(0,1) +
  ylim(0,1) +
  xlab("Experimental PSI Exon 10") +
  ylab("Predicted PSI Exon 10") + 
  theme_classic() +
  theme(text = element_text(size=25),
        axis.text.x = element_text(size=20),axis.text.y = element_text(size=20),
        panel.background = element_rect(colour = "black"), legend.position = "None")

p

#ggsave(p,file="figures/PredictedVsActualPSI_AllMutations_NonDoubtfulOnes_BactComplex_5z56_NonSynonymousAndDouble.svg",width = 10, height = 10, bg = "white")


cor.test(dataToTest$PSI,dataToTest$PredictedPSI)
