setwd("~/Documents/9semester/Master/R")

PheWeb = read.table(file = "~/Documents/9semester/Prosjektoppgave/PheWeb/snp-pheAll.txt", header = T, sep = "\t", dec = ".", stringsAsFactors = F) # reading the PheWas dataset

# Overview of LD blocks associated with each disease

library(sets)
a <- c(101, 102, 103, 114, 132, 142, 164, 171, 173, 174, 20, 205, 209, 29, 5, 58, 6, 8, 97) # Coronary atherosclerosis
b <- c(101, 102, 132, 142, 164, 171, 173, 174, 183, 186, 205, 209, 29, 5, 8, 97) # Myocardial infarction
c <- c(101, 102, 114, 132, 171, 173, 175, 186, 207, 209, 219, 33, 37, 43, 5, 58, 6, 8, 97) # Other chronic ischemic heart disease, unspecified
d <- c(101, 102, 114, 132, 142, 164, 169, 171, 173, 174, 175, 186, 20, 205, 209, 219, 29, 43, 5, 58, 6, 8, 95, 97) # Ischemic heart disease
e <- c(102, 103, 114, 132, 164, 171, 174, 175, 20, 207, 209, 33, 5, 58, 6, 95) # Angina pectoris
f <- c(102, 132) # Peripheral vascular disease 
common_values <- intersect(intersect(a,b),d)
common_values2 <- intersect(intersect(common_values, e), f)

# ------------------------------------------------------------------------------------------------------

# Overview of the SNPs and LD blocks connected to each disease above

CoronaryAtherosclerosis <- PheWeb2[which(PheWeb2$phenostring == "Coronary atherosclerosis"),5]
CA_LDblocks <- CoronaryAtherosclerosis[which(substr(CoronaryAtherosclerosis, 1, 2) == "LD")]
CA_snps <- CoronaryAtherosclerosis[which(substr(CoronaryAtherosclerosis, 1, 2) == "rs")]

MyocardialInfarction <- PheWeb2[which(PheWeb2$phenostring == "Myocardial infarction"),5]
MI_LDblocks <- MyocardialInfarction[which(substr(MyocardialInfarction, 1, 2) == "LD")]
MI_snps <- MyocardialInfarction[which(substr(MyocardialInfarction, 1, 2) == "rs")]

OtherIschemic <- PheWeb2[which(PheWeb2$phenostring == "Other chronic ischemic heart disease, unspecified"),5]
OI_LDblocks <- OtherIschemic[which(substr(OtherIschemic, 1, 2) == "LD")]
OI_snps <- OtherIschemic[which(substr(OtherIschemic, 1, 2) == "rs")]

IschemicHeartDisease <- PheWeb2[which(PheWeb2$phenostring == "Ischemic Heart Disease"),5]
IHD_LDblocks <- IschemicHeartDisease[which(substr(IschemicHeartDisease, 1, 2) == "LD")]
IHD_snps <- IschemicHeartDisease[which(substr(IschemicHeartDisease, 1, 2) == "rs")]

AnginaPectoris <- PheWeb2[which(PheWeb2$phenostring == "Angina pectoris"),5]
AP_LDblocks <- AnginaPectoris[which(substr(AnginaPectoris, 1, 2) == "LD")]
AP_snps <- AnginaPectoris[which(substr(AnginaPectoris, 1, 2) == "rs")]

PulmonaryHeartDisease <- PheWeb2[which(PheWeb2$phenostring == "Pulmonary heart disease"),5]
PHD_LDblocks <- PulmonaryHeartDisease[which(substr(PulmonaryHeartDisease, 1, 2) == "LD")]
PHD_snps <- PulmonaryHeartDisease[which(substr(PulmonaryHeartDisease, 1, 2) == "rs")]

# ------------------------------------------------------------------------------------------------------

# Find ICD9 codes for each disease

ICD9 <- read.csv(file = "icd9.csv", header = T, stringsAsFactors = F, dec = ".", sep = ",") # reading the file containing ICD9-codes
ICD9_CA <- ICD9[which(ICD9$PheCode=="411.4"),1] # Coronary atherosclerosis
ICD9_MI <- ICD9[which(ICD9$PheCode=="411.2"),1] # Myocardial infarction
ICD9_OIHD <- ICD9[which(ICD9$PheCode=="411.8"),1] # Other chronic ischemic heart disease, unspecified
ICD9_IHD <- ICD9[which(ICD9$PheCode=="411"),1] # Ischemic heart disease
ICD9_AP <- ICD9[which(ICD9$PheCode=="411.3"),1] # Angina Pectoris
ICD9_PHD <- ICD9[which(ICD9$PheCode=="415"),1] # Pulmonary heart disease

# Find ICD10 codes for each disease

ICD10 <- read.csv(file = "icd10.csv", header = T, stringsAsFactors = F, dec = ".", sep = ",") # reading the file containing ICD10-codes
ICD10_CA <- ICD10[which(ICD10$phecode=="411.4"),1]
ICD10_MI <- ICD10[which(ICD10$phecode=="411.2"),1]
ICD10_OIHD <- ICD10[which(ICD10$phecode=="411.8"),1]
ICD10_IHD <- ICD10[which(ICD10$phecode=="411"),1] # no ICD codes with phenocode 411.0
ICD10_AP <- ICD10[which(ICD10$phecode=="411.3"),1]
ICD10_PHD <- ICD10[which(ICD10$phecode=="415"),1] # Pulmonary heart disease (to replace ischemic heart disease)


SPN_LD_data <- read.csv(file = "SPN_LD_network_data.csv", header = T, stringsAsFactors = F, sep = ",", dec = ".") # importing network analysis data for SPN with LD
SPN_LD_data <- SPN_LD_data[,-c(3,4)]
CSD <- SPN_LD_data[which(SPN_LD_data$V2==1),] # Circulatory system diseases
CSD <- CSD[order(CSD$Degree),]

# ------------------------------------------------------------------------------------------------------

# Making lists of SNPs to be included in PRS analysis

# Coronary atherosclerosis

a <- list(CA_snps, MI_snps, OI_snps, IHD_snps, AP_snps, PHD_snps)

MakingDisease_snp <- function(snplist){
  cr <- c()
  for(i in 1:length(snplist)){
    cr <- append(cr, PheWeb[which(PheWeb$rsids == snplist[i])[1],1])
  }
  disease_snps <- cbind(snplist, cr)
  disease_snps <- data.frame(disease_snps)
  disease_snps$cr <- as.numeric(disease_snps$cr)
  return(disease_snps[order(disease_snps$cr),])
}

names <- c("CA", "MI", "OI", "IHD", "AP", "PHD")
for(i in 1:length(a)){
  nome <- paste(names[i], "snplist", sep = "_")
  dis_snps <- MakingDisease_snp(a[[i]])
  assign(nome, dis_snps)
}


cr <- c()
for(i in 1:length(CA_snps)){
  cr <- append(cr, PheWeb[which(PheWeb$rsids == CA_snps[i])[1],1])
}

CA_snps <- cbind(CA_snps, cr)
CA_snps <- data.frame(CA_snps)
CA_snps$cr <- as.numeric(CA_snps$cr)
CA_snps <- CA_snps[order(CA_snps$cr),]

ChrSNP_disease <- function(disease, snpList){
  SNPchrList <- list()
  unique_chrom <- unique(snpList[,2])
  for(i in 1:length(unique_chrom)){
    sn <- c()
    common <- which(snpList[,2] == unique_chrom[i])
    if(length(common)>1){
      for(j in 1:length(common)){
        sn <- append(sn, snpList[common[j],1])
      }
    } else{
      sn <- append(sn, snpList[common,1])
    }
    
    SNPchrList <- append(SNPchrList, list(sn))
    w <- which(snpList[,1]==SNPchrList[[i]][1])
    SNPchrList[[i]] <- append(SNPchrList[[i]], snpList[w,2])
  }
  return(SNPchrList)
}


CA_chrom <- ChrSNP_disease("CA",CA_snplist)
MI_chrom <- ChrSNP_disease("MI",MI_snplist)
OI_chrom <- ChrSNP_disease("OI", OI_snplist)
IHD_chrom <- ChrSNP_disease("IHD",IHD_snplist)
AP_chrom <- ChrSNP_disease("AP",AP_snplist)
PHD_chrom <- ChrSNP_disease("PHD",PHD_snplist)

b <- list(CA_chrom, MI_chrom, OI_chrom, IHD_chrom, AP_chrom, PHD_chrom)


for(i in 1:length(CA_chrom)){
  nr <- CA_chrom[[i]][length(CA_chrom[[i]])]
  nam <- paste("CA_chrom", nr, sep = "")
  listt <- CA_chrom[[i]][1:(length(CA_chrom[[i]])-1)]
  assign(nam, listt)
}

for(i in 1:length(MI_chrom)){
  nr <- MI_chrom[[i]][length(MI_chrom[[i]])]
  nam <- paste("MI_chrom", nr, sep = "")
  listt <- MI_chrom[[i]][1:(length(MI_chrom[[i]])-1)]
  assign(nam, listt)
}

for(i in 1:length(OI_chrom)){
  nr <- OI_chrom[[i]][length(OI_chrom[[i]])]
  nam <- paste("OI_chrom", nr, sep = "")
  listt <- OI_chrom[[i]][1:(length(OI_chrom[[i]])-1)]
  assign(nam, listt)
}

for(i in 1:length(IHD_chrom)){
  nr <- IHD_chrom[[i]][length(IHD_chrom[[i]])]
  nam <- paste("IHD_chrom", nr, sep = "")
  listt <- IHD_chrom[[i]][1:(length(IHD_chrom[[i]])-1)]
  assign(nam, listt)
}

for(i in 1:length(AP_chrom)){
  nr <- AP_chrom[[i]][length(AP_chrom[[i]])]
  nam <- paste("AP_chrom", nr, sep = "")
  listt <- AP_chrom[[i]][1:(length(AP_chrom[[i]])-1)]
  assign(nam, listt)
}

for(i in 1:length(PHD_chrom)){
  nr <- PHD_chrom[[i]][length(PHD_chrom[[i]])]
  nam <- paste("PHD_chrom", nr, sep = "")
  listt <- PHD_chrom[[i]][1:(length(PHD_chrom[[i]])-1)]
  assign(nam, listt)
}

# ------------------------------------------------------------------------------------------------------

# Must include SNPs in LD blocks - bare ta med en av SNPsene i LD?

xy <- list(LD1, LD2, LD3, LD4, LD5, LD7, LD8, LD9, LD10, LD11, LD12, LD13, LD14, LD15, LD16, LD17, 
           LD18, LD19, LD20, LD21, LD22)



