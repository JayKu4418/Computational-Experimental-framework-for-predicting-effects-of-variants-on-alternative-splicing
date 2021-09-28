rm(list = ls(all.names = TRUE))

setwd("~/Documents/LaederachLab/Papers/MAPT_splicing_structure/Figure5/")

library(betareg)
library(ggplot2)
library(moments)
library(reshape)
library(hash)


# Function to get standard normal or zscore
zscore <- function(x){
  return((x-mean(x))/sd(x))
}

muts_LowConfidence = c("Mut8","Mut10","Mut16","Mut31","Mut32","Mut36","Mut37","Mut40","Mut43","Mut45","Mut48","Mut51","Mut54","Mut56")

mutIDs_position = c("WT","Mut1","Mut2","Mut4","Mut5","Mut6","Mut7","Mut9","Mut39","Mut11","Mut12","Mut13","Mut14","Mut15","Mut17","Mut41","Mut18","Mut19","Mut42","Mut44","Mut20",
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

# Structure 
data_ToRead_Struc <- read.table("data/AllMuts_Unfold_LengthRNAinSpliceosome_BactComplex_5z56.tsv",sep="\t",header = T)

mean_Energy <- apply(data_ToRead_Struc,2,mean)

sd_Energy <- apply(data_ToRead_Struc,2,sd)

skew_Energy <- apply(data_ToRead_Struc,2,skewness)

kurtosis_Energy <- apply(data_ToRead_Struc,2,kurtosis)

data_toModelAndTest_Struc <- data.frame(strucMean=mean_Energy-mean_Energy[1],strucSD=sd_Energy-sd_Energy[1],strucSkew=skew_Energy-skew_Energy[1],strucKurt=kurtosis_Energy-kurtosis_Energy[1])

data_toModelAndTest <- cbind(data_toModelAndTest_SRE,data_toModelAndTest_Struc)



data_toModelAndTest <- data_toModelAndTest[data_toModelAndTest$MutID %in% allMuts_nonDoubtfulOnes, ]

data_toModelAndTest[,c("strucMean","strucSD","strucSkew","strucKurt","SS","Enhancer","Silencer")] <- apply(data_toModelAndTest[,c("strucMean","strucSD","strucSkew","strucKurt","SS","Enhancer","Silencer")],2,zscore)

#data_toModelAndTest$PSI <- data_toModelAndTest$ExperimentalPSI - 0.56


n=nrow(data_toModelAndTest)
data_toModelAndTest$PSI = (data_toModelAndTest$ExperimentalPSI*(n-1)+0.5)/n


#data_toModel <- data_toModelAndTest[data_toModelAndTest$MutIDs %in% train_set,]

#beta_model <- betareg(PSI~SS+Enhancer+Silencer+strucMean+strucSD+strucSkew+strucKurt, data = data_toModelAndTest)
#beta_model <- betareg(PSI~(SS+Enhancer+Silencer)*(strucMean+strucSD+strucSkew+strucKurt), data = data_toModelAndTest)
#beta_model <- betareg(PSI~strucMean+strucSD+strucSkew+strucKurt, data = data_toModelAndTest)
beta_model <- betareg(PSI~ ((Enhancer + Silencer+SS)*isNonSynonymous)+((strucMean+strucSD+strucSkew+strucKurt)*(isIntronic+isSynonymous)), data = data_toModelAndTest)



summary(beta_model)

data_toPlot <- data_toModelAndTest
data_toPlot$PredictedPSI <- predict(beta_model,data_toModelAndTest)



cor.test(data_toPlot$PSI,data_toPlot$PredictedPSI)

model_dict <- hash()
model_dict[["Structure"]] <-PSI~strucMean+strucSD+strucSkew+strucKurt
model_dict[["SRE"]] <-PSI~SS+Enhancer+Silencer
model_dict[["Combo"]] <-PSI~((Enhancer + Silencer+SS)*isNonSynonymous)+((strucMean+strucSD+strucSkew+strucKurt)*(isIntronic+isSynonymous))
#model_dict[["CombowithClass"]] <-PSI~SS*(strucMean+strucSD+strucSkew+strucKurt)+Enhancer+Silencer

set.seed(20)

corr_coeff_train <- list()
corr_coeff_overall <- list()

ntimes=10

for (i in 1:ntimes){
  
    
    train_set = sample(allMuts_nonDoubtfulOnes,size = 38,replace = FALSE)
    #train_set = allMuts_nonDoubtfulOnes[!(allMuts_nonDoubtfulOnes %in% test_set)]
    #train_set = allMuts_nonDoubtfulOnes
    
    if (!("WT" %in% train_set)){
      train_set = c(train_set,"WT")
    }
    
    data_toModel <- data_toModelAndTest[data_toModelAndTest$MutID %in% train_set,]
    
    corr_coeff_model = list()
    corr_coeff_model_overall = c()
    
    for (modelname in c("Structure","SRE","Combo")){
      tryCatch({
        model <- model_dict[[modelname]]
        print(modelname)
        beta_model <- betareg(model, data = data_toModel)
        
        data_toModel$PredictedPSI <- predict(beta_model,data_toModel)
        
        #data_toTest <- data_toModelAndTest[!(data_toModelAndTest$MutID %in% train_set),]
      
        #data_toTest$PredictedPSI <- predict(beta_model,data_toTest)
        cc_test_overall = cor.test(data_toModel$ExperimentalPSI,data_toModel$PredictedPSI)
        corr_coeff_model_overall <- c(corr_coeff_model_overall,cc_test_overall$estimate)
      
        corr_coeff <- c()
      
        for (type_label in c("NonSynonymous","Synonymous","Intronic","Double")){
        #for (type_label in c("Exonic","Intronic","Double")){
          #print(type_label)
          data_toModel_type = data_toModel[data_toModel$Label==type_label,]
          if (nrow(data_toModel_type) < 3){
            corr_coeff <- c(corr_coeff,NA)
          }
          else{
            cc_test = cor.test(data_toModel_type$ExperimentalPSI,data_toModel_type$PredictedPSI)
            corr_coeff <- c(corr_coeff,cc_test$estimate)
          }
        }
      
        corr_coeff_model[[modelname]] <- corr_coeff
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
    #col_names_df = c(col_names_df,paste("Run",i,sep=""))
    #corr_coeff_train_df = cbind(corr_coeff_train_df,corr_coeff_train)
    #orr_coeff_test_df = cbind(corr_coeff_test_df,corr_coeff_test)
    corr_coeff_train[[i]] <- corr_coeff_model
    corr_coeff_overall[[i]] <- corr_coeff_model_overall
  
}

corr_coeff_train_df <- data.frame(matrix(unlist(corr_coeff_train[[1]]), nrow=length(corr_coeff_train[[1]]), byrow=TRUE))
colnames(corr_coeff_train_df) <- c("Non-Synonymous","Synonymous","Intronic","Double")
#colnames(corr_coeff_train_df) <- c("Exonic","Intronic","Double")
corr_coeff_train_df$Labels <- c("Structure","SRE","Combo")
corr_coeff_train_df$SampleNum <- 1
corr_coeff_train_df_allSamples <- melt(corr_coeff_train_df,id.vars=c("Labels","SampleNum"))


for (samp in 2:ntimes){
  print(samp)
  if (length(corr_coeff_train[[samp]]) == 3) {
  corr_coeff_train_df <- data.frame(matrix(unlist(corr_coeff_train[[samp]]), nrow=length(corr_coeff_train[[samp]]), byrow=TRUE))
  colnames(corr_coeff_train_df) <- c("Non-Synonymous","Synonymous","Intronic","Double")
  #colnames(corr_coeff_train_df) <- c("Exonic","Intronic","Double")
  corr_coeff_train_df$Labels <- c("Structure","SRE","Combo")#,"CombowithClass")
  corr_coeff_train_df$SampleNum <- samp
  corr_coeff_train_df_melted <- melt(corr_coeff_train_df,id.vars=c("Labels","SampleNum"))
  corr_coeff_train_df_allSamples <- rbind(corr_coeff_train_df_allSamples,corr_coeff_train_df_melted)
  }
}



#corr_coeff_train_df_allSamples_Limited <- corr_coeff_train_df_allSamples[corr_coeff_train_df_allSamples$Labels %in% c("Combo","CombowithClass"),]
#cc_vals_mut_combo = corr_coeff_train_df_allSamples_Limited[(corr_coeff_train_df_allSamples_Limited$variable=="Non-Synonymous")&(corr_coeff_train_df_allSamples_Limited$Labels=="Combo"),"value"]
#cc_vals_mut_combowithClass = corr_coeff_train_df_allSamples_Limited[(corr_coeff_train_df_allSamples_Limited$variable=="Non-Synonymous")&(corr_coeff_train_df_allSamples_Limited$Labels=="CombowithClass"), "value"]

#wilcox.test(cc_vals_mut_combo,cc_vals_mut_combowithClass)

p <- ggplot(corr_coeff_train_df_allSamples,aes(x=factor(Labels,levels = c("Structure","SRE","Combo")),y=value,fill=variable))+
  #geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_violin(position=position_dodge(1),scale = "width") +
  geom_boxplot(position=position_dodge(1),outlier.shape = NA,width=0.35,lwd=0.25) + 
  #stat_summary(fun=median, geom="point", size=1.5, color="black")+
  scale_fill_manual(values=c('#F52549','#A13EDB','#3C53DD','#A2D5C6')) +
  #coord_cartesian(y=c(0,1))+
  xlab("")+
  ylab("Correlation coefficient")+
  #scale_x_discrete(labels=c("Structure + SRE strength","Structure + SRE strength \n interact with Mutation type"))+
theme_classic()+
  theme(text = element_text(size=25),
        axis.text.x = element_text(size=20),axis.text.y = element_text(size=20),
        panel.background = element_rect(colour = "black"),legend.position = "None")
p

#ggsave(p,file="figures/ViolinPlots_44Mutations_TrainTest_CombinedFeatures.svg",width = 15, height = 10, bg = "white")


for (i in 1:ntimes){
  if (length(corr_coeff_overall[[i]]) !=3){
    corr_coeff_overall[[i]] <- c(corr_coeff_overall[[i]],NA)
  }
}

corr_coeff_overall_df <- data.frame(matrix(unlist(corr_coeff_overall), nrow=length(corr_coeff_overall), byrow=TRUE))
colnames(corr_coeff_overall_df) <- c("Structure","SRE","Combo")
corr_coeff_overall_df_melted = melt(corr_coeff_overall_df)

colnames(corr_coeff_overall_df_melted) <- c("Labels","value")
corr_coeff_overall_df_melted$variable <- "Overall"


corr_coeff_overall_df_melted_ALL <- rbind(corr_coeff_overall_df_melted,corr_coeff_train_df_allSamples[,c("Labels","value","variable")])
corr_coeff_overall_df_melted_ALL$variable <- factor(corr_coeff_overall_df_melted_ALL$variable, levels = c('Non-Synonymous', 'Synonymous', 'Intronic',"Double","Overall"))

p <- ggplot(corr_coeff_overall_df_melted_ALL,aes(x=factor(Labels,levels = c("Structure","SRE","Combo")),y=value,fill=variable))+
  #geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_violin(position=position_dodge(1),scale = "width") +
  geom_boxplot(position=position_dodge(1),outlier.shape = NA,width=0.35,lwd=0.25) + 
  #stat_summary(fun=median, geom="point", size=1.5, color="black")+
  scale_fill_manual(values=c('#F52549','#A13EDB','#3C53DD','#A2D5C6',"#C95D00")) +
  #coord_cartesian(y=c(0,1))+
  xlab("")+
  ylab("Correlation Coefficient (R2)")+
  #scale_x_discrete(labels=c("Structure","SRE motif","Structure*SRE motif"))+
  theme_classic()+
  theme(text = element_text(size=25),
        axis.text.x = element_text(size=20),axis.text.y = element_text(size=20),
        panel.background = element_rect(colour = "black"),legend.position = "None")+
  facet_wrap(~ Labels,scales = "free_x")
p

