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

# ===========================================================================================================

# Overview of the SNPs and LD blocks connected to each disease above

CoronaryAtherosclerosis <- PheWeb[which(PheWeb$phenostring == "Coronary atherosclerosis"),5] # Finding SNPs associated with CA
CoronaryAtherosclerosiswLD <- PheWeb2[which(PheWeb2$phenostring == "Coronary atherosclerosis"),5] # Finding LD blocks associated with CA
CA_LDblocks <- CoronaryAtherosclerosiswLD[which(substr(CoronaryAtherosclerosiswLD, 1, 2) == "LD")]
CA_snps <- CoronaryAtherosclerosis[which(substr(CoronaryAtherosclerosis, 1, 2) == "rs")]

MyocardialInfarction <- PheWeb[which(PheWeb$phenostring == "Myocardial infarction"),5] # Finding SNPs associated with MI
MyocardialInfarctionwLD <- PheWeb2[which(PheWeb2$phenostring == "Myocardial infarction"),5] # Finding LD blocks associated with CA
MI_LDblocks <- MyocardialInfarctionwLD[which(substr(MyocardialInfarctionwLD, 1, 2) == "LD")]
MI_snps <- MyocardialInfarction[which(substr(MyocardialInfarction, 1, 2) == "rs")]

OtherIschemic <- PheWeb[which(PheWeb$phenostring == "Other chronic ischemic heart disease, unspecified"),5] # Finding SNPs associated with OI
OtherIschemicwLD <- PheWeb2[which(PheWeb2$phenostring == "Other chronic ischemic heart disease, unspecified"),5] # Finding LD blocks associated with CA
OI_LDblocks <- OtherIschemicwLD[which(substr(OtherIschemicwLD, 1, 2) == "LD")]
OI_snps <- OtherIschemic[which(substr(OtherIschemic, 1, 2) == "rs")]

IschemicHeartDisease <- PheWeb[which(PheWeb$phenostring == "Ischemic Heart Disease"),5] # Finding SNPs associated with IHD
IschemicHeartDiseasewLD <- PheWeb2[which(PheWeb2$phenostring == "Ischemic Heart Disease"),5] # Finding LD blocks associated with CA
IHD_LDblocks <- IschemicHeartDiseasewLD[which(substr(IschemicHeartDiseasewLD, 1, 2) == "LD")]
IHD_snps <- IschemicHeartDisease[which(substr(IschemicHeartDisease, 1, 2) == "rs")]

AnginaPectoris <- PheWeb[which(PheWeb$phenostring == "Angina pectoris"),5] # Finding SNPs associated with AP
AnginaPectoriswLD <- PheWeb2[which(PheWeb2$phenostring == "Angina pectoris"),5] # Finding LD blocks associated with CA
AP_LDblocks <- AnginaPectoriswLD[which(substr(AnginaPectoriswLD, 1, 2) == "LD")]
AP_snps <- AnginaPectoris[which(substr(AnginaPectoris, 1, 2) == "rs")]

PulmonaryHeartDisease <- PheWeb[which(PheWeb$phenostring == "Pulmonary heart disease"),5] # Finding SNPs associated with PHD
PulmonaryHeartDiseasewLD <- PheWeb2[which(PheWeb2$phenostring == "Pulmonary heart disease"),5] # Finding LD blocks associated with CA
PHD_LDblocks <- PulmonaryHeartDiseasewLD[which(substr(PulmonaryHeartDiseasewLD, 1, 2) == "LD")]
PHD_snps <- PulmonaryHeartDisease[which(substr(PulmonaryHeartDisease, 1, 2) == "rs")]

# ===========================================================================================================

# Find ICD9 codes for each disease

ICD9 <- read.csv(file = "icd9.csv", header = T, stringsAsFactors = F, dec = ".", sep = ",") # reading the file containing ICD9-codes
ICD9_CA <- ICD9[which(ICD9$PheCode=="411.4"),1] # Coronary atherosclerosis
ICD9_MI <- data.frame(ICD9[which(ICD9$PheCode=="411.2"),1]) # Myocardial infarction
ICD9_OIHD <- data.frame(ICD9[which(ICD9$PheCode=="411.8"),1]) # Other chronic ischemic heart disease, unspecified
ICD9_IHD <- ICD9[which(ICD9$PheCode=="411"),1] # Ischemic heart disease
ICD9_AP <- ICD9[which(ICD9$PheCode=="411.3"),1] # Angina Pectoris
ICD9_PHD <- data.frame(ICD9[which(ICD9$PheCode=="415"),1]) # Pulmonary heart disease
ICD9_EH <- data.frame(ICD9[which(ICD9$PheCode == "401.1"),1])
ICD9_PT <- data.frame(ICD9[which(ICD9$PheCode == "451.2"),1])
ICD9_HF <- data.frame(ICD9[which(ICD9$PheCode == "428.2"),1])

# Find ICD10 codes for each disease

ICD10 <- read.csv(file = "icd10.csv", header = T, stringsAsFactors = F, dec = ".", sep = ",") # reading the file containing ICD10-codes
ICD10_CA <- data.frame(ICD10[which(ICD10$phecode=="411.4"),1])
ICD10_MI <- data.frame(ICD10[which(ICD10$phecode=="411.2"),1])
ICD10_OIHD <- data.frame(ICD10[which(ICD10$phecode=="411.8"),1])
ICD10_IHD <- ICD10[which(ICD10$phecode=="411"),1] # no ICD codes with phenocode 411.0
ICD10_AP <- ICD10[which(ICD10$phecode=="411.3"),1]
ICD10_PHD <- data.frame(ICD10[which(ICD10$phecode=="415"),1])# Pulmonary heart disease (to replace ischemic heart disease)
ICD10_EH <- data.frame(ICD10[which(ICD10$phecode == "401.1"),1])
ICD10_PT <- data.frame(ICD10[which(ICD10$phecode == "451.2"),1])
ICD10_HF <- data.frame(ICD10[which(ICD10$phecode == "428.2"),1])

SPN_LD_data <- read.csv(file = "SPN_LD_network_data.csv", header = T, stringsAsFactors = F, sep = ",", dec = ".") # importing network analysis data for SPN with LD
SPN_LD_data <- SPN_LD_data[,-c(3,4)]
CSD <- SPN_LD_data[which(SPN_LD_data$V2==1),] # Circulatory system diseases
CSD <- CSD[order(CSD$Degree),]

# ===========================================================================================================

# Making lists of SNPs to be included in PRS analysis

a <- list(CA_snps, MI_snps, OI_snps, IHD_snps, AP_snps, PHD_snps) # Making a list of lists of the SNPs associated with each disease

MakingDisease_snp <- function(snplist){                           # Function that makes a matrix containing 
  cr <- c()                                                       # the SNPs associated with a particular   
  for(i in 1:length(snplist)){                                    # disease, and their chromosome number 
    cr <- append(cr, PheWeb[which(PheWeb$rsids == snplist[i])[1],1])
  }
  disease_snps <- cbind(snplist, cr)
  disease_snps <- data.frame(disease_snps)
  disease_snps$cr <- as.numeric(disease_snps$cr)
  return(disease_snps[order(disease_snps$cr),])
}

names <- c("CA", "MI", "OI", "IHD", "AP", "PHD")

for(i in 1:length(a)){                                           # Running the MakingDisease_snp function on
  nome <- paste(names[i], "snplist", sep = "_")                  # on each disease
  dis_snps <- MakingDisease_snp(a[[i]])
  assign(nome, dis_snps)
}


ChrSNP_disease <- function(disease, snpList){                    # A function that creates a list of lists, 
  SNPchrList <- list()                                           # where each list contains the SNPs on a 
  unique_chrom <- unique(snpList[,2])                            # particular chromosome
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

CAlist <- c()
for(i in 1:length(CA_chrom)){                                      # A for-loop which makes a SNP-vector for 
  nr <- CA_chrom[[i]][length(CA_chrom[[i]])]                       # each chromosome.
  nam <- paste("CA_chrom", nr, sep = "")                           # This is repeated for each disease.
  listt <- CA_chrom[[i]][1:(length(CA_chrom[[i]])-1)]
  listt <- cbind(rep(nr, length(listt)), listt)
  assign(nam, listt)
  if(i == 1){
    CAlist <- listt
  } else{
    CAlist <- rbind(CAlist, listt)
  }
}

MIlist <- list()
for(i in 1:length(MI_chrom)){
  nr <- MI_chrom[[i]][length(MI_chrom[[i]])]
  nam <- paste("MI_chrom", nr, sep = "")
  listt <- MI_chrom[[i]][1:(length(MI_chrom[[i]])-1)]
  listt <- cbind(rep(nr, length(listt)), listt)
  assign(nam, listt)
  if(i == 1){
    MIlist <- listt
  } else{
    MIlist <- rbind(MIlist, listt)
  }
}

OIlist <- list()
for(i in 1:length(OI_chrom)){
  nr <- OI_chrom[[i]][length(OI_chrom[[i]])]
  nam <- paste("OI_chrom", nr, sep = "")
  listt <- OI_chrom[[i]][1:(length(OI_chrom[[i]])-1)]
  listt <- cbind(rep(nr, length(listt)), listt)
  assign(nam, listt)
  if(i == 1){
    OIlist <- listt
  } else{
    OIlist <- rbind(OIlist, listt)
  }
}

IHDlist <- list()
for(i in 1:length(IHD_chrom)){
  nr <- IHD_chrom[[i]][length(IHD_chrom[[i]])]
  nam <- paste("IHD_chrom", nr, sep = "")
  listt <- IHD_chrom[[i]][1:(length(IHD_chrom[[i]])-1)]
  listt <- cbind(rep(nr, length(listt)), listt)
  assign(nam, listt)
  if(i == 1){
    IHDlist <- listt
  } else{
    IHDlist <- rbind(IHDlist, listt)
  }
}

APlist <- list()
for(i in 1:length(AP_chrom)){
  nr <- AP_chrom[[i]][length(AP_chrom[[i]])]
  nam <- paste("AP_chrom", nr, sep = "")
  listt <- AP_chrom[[i]][1:(length(AP_chrom[[i]])-1)]
  listt <- cbind(rep(nr, length(listt)), listt)
  assign(nam, listt)
  if(i == 1){
    APlist <- listt
  } else{
    APlist <- rbind(APlist, listt)
  }
}

PHDlist <- list()
for(i in 1:length(PHD_chrom)){
  nr <- PHD_chrom[[i]][length(PHD_chrom[[i]])]
  nam <- paste("PHD_chrom", nr, sep = "")
  listt <- PHD_chrom[[i]][1:(length(PHD_chrom[[i]])-1)]
  listt <- cbind(rep(nr, length(listt)), listt)
  assign(nam, listt)
  if(i == 1){
    PHDlist <- listt
  } else{
    PHDlist <- rbind(PHDlist, listt)
  }
}

# ===========================================================================================================


# Make one big list with all SNPs - include CA, MI, OI, AP, PHD,
# and one list for each disease 

CAsnp <- c()
for(i in 1:length(CAlist)){
  for(j in 1:length(CAlist[[i]])){
    AllSNPs <- append(AllSNPs, CAlist[[i]][j])
    CAsnp <- append(CAsnp, CAlist[[i]][j])
  }
}
MIsnp <- c()
for(i in 1:length(MIlist)){
  for(j in 1:length(MIlist[[i]])){
    AllSNPs <- append(AllSNPs, MIlist[[i]][j])
    MIsnp <- append(MIsnp, MIlist[[i]][j])
  }
}
OIsnp <- c()
for(i in 1:length(OIlist)){
  for(j in 1:length(OIlist[[i]])){
    AllSNPs <- append(AllSNPs, OIlist[[i]][j])
    OIsnp <- append(OIsnp, OIlist[[i]][j])
  }
}
APsnp <- c()
for(i in 1:length(APlist)){
  for(j in 1:length(APlist[[i]])){
    AllSNPs <- append(AllSNPs, APlist[[i]][j])
    APsnp <- append(APsnp, APlist[[i]][j])
  }
}
PHDsnp <- c()
for(i in 1:length(PHDlist)){
  for(j in 1:length(PHDlist[[i]])){
    AllSNPs <- append(AllSNPs, PHDlist[[i]][j])
    PHDsnp <- append(PHDsnp, PHDlist[[i]][j])
  }
}

# ===========================================================================================================

# Making files containing SNPs for each disease 

write.table(AllSNPs, file = "AllSNPs.txt", sep = "\t", row.names = F, col.names = F) # all SNPs from every disease

# All SNP lists for coronary atherosclerosis, for each chromosome

write.table(CAlist, file = "CAsnps.txt", sep = "\t", row.names = F, col.names = F) # all SNPs
write.table(CA_chrom1, file = "~/Documents/9semester/Master/R/CA/CA_chrom1.txt", sep = "\t", row.names = F, col.names = F)
write.table(CA_chrom2, file = "~/Documents/9semester/Master/R/CA/CA_chrom2.txt", sep = "\t", row.names = F, col.names = F)
write.table(CA_chrom3, file = "~/Documents/9semester/Master/R/CA/CA_chrom3.txt", sep = "\t", row.names = F, col.names = F)
write.table(CA_chrom4, file = "~/Documents/9semester/Master/R/CA/CA_chrom4.txt", sep = "\t", row.names = F, col.names = F)
write.table(CA_chrom6, file = "~/Documents/9semester/Master/R/CA/CA_chrom6.txt", sep = "\t", row.names = F, col.names = F)
write.table(CA_chrom7, file = "~/Documents/9semester/Master/R/CA/CA_chrom7.txt", sep = "\t", row.names = F, col.names = F)
write.table(CA_chrom8, file = "~/Documents/9semester/Master/R/CA/CA_chrom8.txt", sep = "\t", row.names = F, col.names = F)
write.table(CA_chrom9, file = "~/Documents/9semester/Master/R/CA/CA_chrom9.txt", sep = "\t", row.names = F, col.names = F)
write.table(CA_chrom10, file = "~/Documents/9semester/Master/R/CA/CA_chrom10.txt", sep = "\t", row.names = F, col.names = F)
write.table(CA_chrom11, file = "~/Documents/9semester/Master/R/CA/CA_chrom11.txt", sep = "\t", row.names = F, col.names = F)
write.table(CA_chrom12, file = "~/Documents/9semester/Master/R/CA/CA_chrom12.txt", sep = "\t", row.names = F, col.names = F)
write.table(CA_chrom13, file = "~/Documents/9semester/Master/R/CA/CA_chrom13.txt", sep = "\t", row.names = F, col.names = F)
write.table(CA_chrom15, file = "~/Documents/9semester/Master/R/CA/CA_chrom15.txt", sep = "\t", row.names = F, col.names = F)
write.table(CA_chrom16, file = "~/Documents/9semester/Master/R/CA/CA_chrom16.txt", sep = "\t", row.names = F, col.names = F)
write.table(CA_chrom17, file = "~/Documents/9semester/Master/R/CA/CA_chrom17.txt", sep = "\t", row.names = F, col.names = F)
write.table(CA_chrom19, file = "~/Documents/9semester/Master/R/CA/CA_chrom19.txt", sep = "\t", row.names = F, col.names = F)
write.table(CA_chrom21, file = "~/Documents/9semester/Master/R/CA/CA_chrom21.txt", sep = "\t", row.names = F, col.names = F)

# ===========================================================================================================

# All SNP lists for myocardial infarction, for each chromosome

write.table(MIlist, file = "MIsnps.txt", sep = "\t", row.names = F, col.names = F)
write.table(MI_chrom1, file = "~/Documents/9semester/Master/R/MI/MI_chrom1.txt", sep = "\t", row.names = F, col.names = F)
write.table(MI_chrom2, file = "~/Documents/9semester/Master/R/MI/MI_chrom2.txt", sep = "\t", row.names = F, col.names = F)
write.table(MI_chrom3, file = "~/Documents/9semester/Master/R/MI/MI_chrom3.txt", sep = "\t", row.names = F, col.names = F)
write.table(MI_chrom4, file = "~/Documents/9semester/Master/R/MI/MI_chrom4.txt", sep = "\t", row.names = F, col.names = F)
write.table(MI_chrom5, file = "~/Documents/9semester/Master/R/MI/MI_chrom5.txt", sep = "\t", row.names = F, col.names = F)
write.table(MI_chrom6, file = "~/Documents/9semester/Master/R/MI/MI_chrom6.txt", sep = "\t", row.names = F, col.names = F)
write.table(MI_chrom7, file = "~/Documents/9semester/Master/R/MI/MI_chrom7.txt", sep = "\t", row.names = F, col.names = F)
write.table(MI_chrom8, file = "~/Documents/9semester/Master/R/MI/MI_chrom8.txt", sep = "\t", row.names = F, col.names = F)
write.table(MI_chrom9, file = "~/Documents/9semester/Master/R/MI/MI_chrom9.txt", sep = "\t", row.names = F, col.names = F)
write.table(MI_chrom10, file = "~/Documents/9semester/Master/R/MI/MI_chrom10.txt", sep = "\t", row.names = F, col.names = F)
write.table(MI_chrom11, file = "~/Documents/9semester/Master/R/MI/MI_chrom11.txt", sep = "\t", row.names = F, col.names = F)
write.table(MI_chrom12, file = "~/Documents/9semester/Master/R/MI/MI_chrom12.txt", sep = "\t", row.names = F, col.names = F)
write.table(MI_chrom13, file = "~/Documents/9semester/Master/R/MI/MI_chrom13.txt", sep = "\t", row.names = F, col.names = F)
write.table(MI_chrom14, file = "~/Documents/9semester/Master/R/MI/MI_chrom14.txt", sep = "\t", row.names = F, col.names = F)
write.table(MI_chrom15, file = "~/Documents/9semester/Master/R/MI/MI_chrom15.txt", sep = "\t", row.names = F, col.names = F)
write.table(MI_chrom16, file = "~/Documents/9semester/Master/R/MI/MI_chrom16.txt", sep = "\t", row.names = F, col.names = F)
write.table(MI_chrom17, file = "~/Documents/9semester/Master/R/MI/MI_chrom17.txt", sep = "\t", row.names = F, col.names = F)
write.table(MI_chrom19, file = "~/Documents/9semester/Master/R/MI/MI_chrom19.txt", sep = "\t", row.names = F, col.names = F)
write.table(MI_chrom21, file = "~/Documents/9semester/Master/R/MI/MI_chrom21.txt", sep = "\t", row.names = F, col.names = F)
write.table(MI_chrom22, file = "~/Documents/9semester/Master/R/MI/MI_chrom22.txt", sep = "\t", row.names = F, col.names = F)

# ===========================================================================================================

# All SNP lists for Other chronic ischemic heart disease, for each chromosome

write.table(OIlist, file = "OIsnps.txt", sep = "\t", row.names = F, col.names = F)
write.table(OI_chrom1, file = "~/Documents/9semester/Master/R/OI/OI_chrom1.txt", sep = "\t", row.names = F, col.names = F)
write.table(OI_chrom2, file = "~/Documents/9semester/Master/R/OI/OI_chrom2.txt", sep = "\t", row.names = F, col.names = F)
write.table(OI_chrom3, file = "~/Documents/9semester/Master/R/OI/OI_chrom3.txt", sep = "\t", row.names = F, col.names = F)
write.table(OI_chrom4, file = "~/Documents/9semester/Master/R/OI/OI_chrom4.txt", sep = "\t", row.names = F, col.names = F)
write.table(OI_chrom6, file = "~/Documents/9semester/Master/R/OI/OI_chrom6.txt", sep = "\t", row.names = F, col.names = F)
write.table(OI_chrom7, file = "~/Documents/9semester/Master/R/OI/OI_chrom7.txt", sep = "\t", row.names = F, col.names = F)
write.table(OI_chrom8, file = "~/Documents/9semester/Master/R/OI/OI_chrom8.txt", sep = "\t", row.names = F, col.names = F)
write.table(OI_chrom9, file = "~/Documents/9semester/Master/R/OI/OI_chrom9.txt", sep = "\t", row.names = F, col.names = F)
write.table(OI_chrom10, file = "~/Documents/9semester/Master/R/OI/OI_chrom10.txt", sep = "\t", row.names = F, col.names = F)
write.table(OI_chrom11, file = "~/Documents/9semester/Master/R/OI/OI_chrom11.txt", sep = "\t", row.names = F, col.names = F)
write.table(OI_chrom12, file = "~/Documents/9semester/Master/R/OI/OI_chrom12.txt", sep = "\t", row.names = F, col.names = F)
write.table(OI_chrom13, file = "~/Documents/9semester/Master/R/OI/OI_chrom13.txt", sep = "\t", row.names = F, col.names = F)
write.table(OI_chrom14, file = "~/Documents/9semester/Master/R/OI/OI_chrom14.txt", sep = "\t", row.names = F, col.names = F)
write.table(OI_chrom15, file = "~/Documents/9semester/Master/R/OI/OI_chrom15.txt", sep = "\t", row.names = F, col.names = F)
write.table(OI_chrom16, file = "~/Documents/9semester/Master/R/OI/OI_chrom16.txt", sep = "\t", row.names = F, col.names = F)
write.table(OI_chrom17, file = "~/Documents/9semester/Master/R/OI/OI_chrom17.txt", sep = "\t", row.names = F, col.names = F)
write.table(OI_chrom18, file = "~/Documents/9semester/Master/R/OI/OI_chrom18.txt", sep = "\t", row.names = F, col.names = F)
write.table(OI_chrom19, file = "~/Documents/9semester/Master/R/OI/OI_chrom19.txt", sep = "\t", row.names = F, col.names = F)
write.table(OI_chrom20, file = "~/Documents/9semester/Master/R/OI/OI_chrom20.txt", sep = "\t", row.names = F, col.names = F)
write.table(OI_chrom21, file = "~/Documents/9semester/Master/R/OI/OI_chrom21.txt", sep = "\t", row.names = F, col.names = F)

# ===========================================================================================================

# All SNP lists for angina pectoris, for each chromosome

write.table(APlist, file = "APsnps.txt", sep = "\t", row.names = F, col.names = F)
write.table(AP_chrom1, file = "~/Documents/9semester/Master/R/AP/AP_chrom1.txt", sep = "\t", row.names = F, col.names = F)
write.table(AP_chrom2, file = "~/Documents/9semester/Master/R/AP/AP_chrom2.txt", sep = "\t", row.names = F, col.names = F)
write.table(AP_chrom3, file = "~/Documents/9semester/Master/R/AP/AP_chrom3.txt", sep = "\t", row.names = F, col.names = F)
write.table(AP_chrom4, file = "~/Documents/9semester/Master/R/AP/AP_chrom5.txt", sep = "\t", row.names = F, col.names = F)
write.table(AP_chrom5, file = "~/Documents/9semester/Master/R/AP/AP_chrom4.txt", sep = "\t", row.names = F, col.names = F)
write.table(AP_chrom6, file = "~/Documents/9semester/Master/R/AP/AP_chrom6.txt", sep = "\t", row.names = F, col.names = F)
write.table(AP_chrom7, file = "~/Documents/9semester/Master/R/AP/AP_chrom7.txt", sep = "\t", row.names = F, col.names = F)
write.table(AP_chrom8, file = "~/Documents/9semester/Master/R/AP/AP_chrom8.txt", sep = "\t", row.names = F, col.names = F)
write.table(AP_chrom9, file = "~/Documents/9semester/Master/R/AP/AP_chrom9.txt", sep = "\t", row.names = F, col.names = F)
write.table(AP_chrom10, file = "~/Documents/9semester/Master/R/AP/AP_chrom10.txt", sep = "\t", row.names = F, col.names = F)
write.table(AP_chrom11, file = "~/Documents/9semester/Master/R/AP/AP_chrom11.txt", sep = "\t", row.names = F, col.names = F)
write.table(AP_chrom12, file = "~/Documents/9semester/Master/R/AP/AP_chrom12.txt", sep = "\t", row.names = F, col.names = F)
write.table(AP_chrom13, file = "~/Documents/9semester/Master/R/AP/AP_chrom13.txt", sep = "\t", row.names = F, col.names = F)
write.table(AP_chrom14, file = "~/Documents/9semester/Master/R/AP/AP_chrom14.txt", sep = "\t", row.names = F, col.names = F)
write.table(AP_chrom15, file = "~/Documents/9semester/Master/R/AP/AP_chrom15.txt", sep = "\t", row.names = F, col.names = F)
write.table(AP_chrom16, file = "~/Documents/9semester/Master/R/AP/AP_chrom16.txt", sep = "\t", row.names = F, col.names = F)
write.table(AP_chrom17, file = "~/Documents/9semester/Master/R/AP/AP_chrom17.txt", sep = "\t", row.names = F, col.names = F)
write.table(AP_chrom18, file = "~/Documents/9semester/Master/R/AP/AP_chrom18.txt", sep = "\t", row.names = F, col.names = F)
write.table(AP_chrom19, file = "~/Documents/9semester/Master/R/AP/AP_chrom19.txt", sep = "\t", row.names = F, col.names = F)
write.table(AP_chrom21, file = "~/Documents/9semester/Master/R/AP/AP_chrom21.txt", sep = "\t", row.names = F, col.names = F)

# ===========================================================================================================

# All SNP lists for pulmonary heart disease, for each chromosome

write.table(PHDlist, file = "PHDsnps.txt", sep = "\t", row.names = F, col.names = F)
write.table(PHD_chrom1, file = "~/Documents/9semester/Master/R/PHD/PHD_chrom1.txt", sep = "\t", row.names = F, col.names = F)
write.table(PHD_chrom3, file = "~/Documents/9semester/Master/R/PHD/PHD_chrom3.txt", sep = "\t", row.names = F, col.names = F)
write.table(PHD_chrom4, file = "~/Documents/9semester/Master/R/PHD/PHD_chrom5.txt", sep = "\t", row.names = F, col.names = F)
write.table(PHD_chrom5, file = "~/Documents/9semester/Master/R/PHD/PHD_chrom4.txt", sep = "\t", row.names = F, col.names = F)
write.table(PHD_chrom7, file = "~/Documents/9semester/Master/R/PHD/PHD_chrom7.txt", sep = "\t", row.names = F, col.names = F)
write.table(PHD_chrom8, file = "~/Documents/9semester/Master/R/PHD/PHD_chrom8.txt", sep = "\t", row.names = F, col.names = F)
write.table(PHD_chrom9, file = "~/Documents/9semester/Master/R/PHD/PHD_chrom9.txt", sep = "\t", row.names = F, col.names = F)
write.table(PHD_chrom10, file = "~/Documents/9semester/Master/R/PHD/PHD_chrom10.txt", sep = "\t", row.names = F, col.names = F)
write.table(PHD_chrom11, file = "~/Documents/9semester/Master/R/PHD/PHD_chrom11.txt", sep = "\t", row.names = F, col.names = F)
write.table(PHD_chrom12, file = "~/Documents/9semester/Master/R/PHD/PHD_chrom12.txt", sep = "\t", row.names = F, col.names = F)
write.table(PHD_chrom14, file = "~/Documents/9semester/Master/R/PHD/PHD_chrom14.txt", sep = "\t", row.names = F, col.names = F)
write.table(PHD_chrom16, file = "~/Documents/9semester/Master/R/PHD/PHD_chrom16.txt", sep = "\t", row.names = F, col.names = F)
write.table(PHD_chrom19, file = "~/Documents/9semester/Master/R/PHD/PHD_chrom19.txt", sep = "\t", row.names = F, col.names = F)
write.table(PHD_chrom20, file = "~/Documents/9semester/Master/R/PHD/PHD_chrom20.txt", sep = "\t", row.names = F, col.names = F)


# ===========================================================================================================





