setwd("../EMBO Popgen/embo_popgen_2025_DCR/Me/Day6/EM ")
setwd("EMBO_Practical_Course_2024-main/")

AFR_EAS_Weir.fst<-read.table("AFR_EAS.weir.fst")
AFR_EUR_Weir.fst<-read.table("AFR_EUR.weir.fst")
EAS_EUR_Weir.fst<-read.table("EAS_EUR.weir.fst")

AFR_EAS_Weir.fst<-AFR_EAS_Weir.fst[-1,]

colnames(AFR_EAS_Weir.fst)<-c("Chrom","POS","fst")

AFR_EUR_Weir.fst<-AFR_EUR_Weir.fst[-1,]

colnames(AFR_EUR_Weir.fst)<-c("Chrom","POS","fst")

EAS_EUR_Weir.fst<-EAS_EUR_Weir.fst[-1,]

colnames(EAS_EUR_Weir.fst)<-c("Chrom","POS","fst")

EAS_EUR_Weir.fst<-EAS_EUR_Weir.fst %>% mutate(across(everything(), as.numeric))
AFR_EAS_Weir.fst<-AFR_EAS_Weir.fst %>% mutate(across(everything(), as.numeric))
AFR_EUR_Weir.fst<-AFR_EUR_Weir.fst %>% mutate(across(everything(), as.numeric))

library(dplyr)
length(AFR_EAS_Weir.fst$POS)
AFR_EAS_Weir.fst <- distinct(AFR_EAS_Weir.fst, POS,.keep_all = TRUE)
AFR_EUR_Weir.fst <- distinct(AFR_EUR_Weir.fst, POS, .keep_all = TRUE)
EAS_EUR_Weir.fst <- distinct(AFR_EAS_Weir.fst, POS, .keep_all = TRUE)

AFR_EAS_Weir.fst <-AFR_EAS_Weir.fst[!is.na(AFR_EAS_Weir.fst$fst),]
AFR_EUR_Weir.fst <-AFR_EUR_Weir.fst[!is.na(AFR_EUR_Weir.fst$fst),]
EAS_EUR_Weir.fst <-EAS_EUR_Weir.fst[!is.na(EAS_EUR_Weir.fst$fst),]
length(AFR_EAS_Weir.fst$POS)


#Overlapping 

overlapping_pos <- intersect(intersect(AFR_EAS_Weir.fst$POS, EAS_EUR_Weir.fst$POS), AFR_EUR_Weir.fst$POS)

AFR_EAS_Weir.fst <- AFR_EAS_Weir.fst %>% filter(POS %in% overlapping_pos)
AFR_EUR_Weir.fst <- AFR_EUR_Weir.fst %>% filter(POS %in% overlapping_pos)
EAS_EUR_Weir.fst <- EAS_EUR_Weir.fst %>% filter(POS %in% overlapping_pos)


head(AFR_EAS_Weir.fst)
length(AFR_EAS_Weir.fst$POS)


length(EAS_EUR_Weir.fst$POS)
length(AFR_EUR_Weir.fst$POS)
length(AFR_EAS_Weir.fst$POS)

EAS_EUR_Weir.fst$fst[EAS_EUR_Weir.fst$fst<0]<-0
AFR_EUR_Weir.fst$fst[AFR_EUR_Weir.fst$fst<0]<-0
AFR_EAS_Weir.fst$fst[AFR_EAS_Weir.fst$fst<0]<-0

#Check if the SNP rs3827760, located at position 109513601, 
#is an outlier in the FST distribution for any of the populatio n pair


POSs<-109513601
EAS_EUR_Weir.fst$fst[EAS_EUR_Weir.fst$POS==POS]
AFR_EUR_Weir.fst$fst[AFR_EUR_Weir.fst$POS==POS]
AFR_EAS_Weir.fst$fst[AFR_EAS_Weir.fst$POS==POS]
#Check which quartile percentile the rs3827760 
#distribution fall in each analyzed population pair?

quantile(as.numeric(EAS_EUR_Weir.fst$fst),probs = c(.05,.1,.5,.25,.5,.75,.95,.99,.995))
#5%          10%          50%          25%          50% 
#0.0000674337 0.0033766800 0.1117750000 0.0278812750 0.1117750000 
#75%          95%          99%        99.5% 
#  0.2225142500 0.4648131000 0.6490841100 0.7049496650 


X<- - log(1-EAS_EUR_Weir.fst$fst)
Y<- - log(1-AFR_EAS_Weir.fst$fst)
Y<- - log(1-AFR_EUR_Weir.fst$fst)

PBS<-(X+Y-Z)/2

PBS[PBS<0]<-0

EAS_EUR_Weir.fst$PBS<-PBS

EAS_EUR_Weir.fst$PBS[EAS_EUR_Weir.fst$POS==POSs]
quantile(as.numeric(EAS_EUR_Weir.fst$PBS),probs = c(.05,.1,.5,.25,.5,.75,.95,.99,.995))

plot(x=EAS_EUR_Weir.fst$POS[EAS_EUR_Weir.fst$POS>=(POSs-10000) & EAS_EUR_Weir.fst$POS<=(POSs+10000)],
     EAS_EUR_Weir.fst$fst[EAS_EUR_Weir.fst$POS>=(POSs-10000) & EAS_EUR_Weir.fst$POS<=(POSs+10000)])

plot(x=EAS_EUR_Weir.fst$POS[EAS_EUR_Weir.fst$POS>=(POSs-10000) & EAS_EUR_Weir.fst$POS<=(POSs+10000)],
     EAS_EUR_Weir.fst$PBS[EAS_EUR_Weir.fst$POS>=(POSs-10000) & EAS_EUR_Weir.fst$POS<=(POSs+10000)])
points(x=EAS_EUR_Weir.fst$POS[EAS_EUR_Weir.fst$POS>=(POSs-10000) & EAS_EUR_Weir.fst$POS<=(POSs+10000) & EAS_EUR_Weir.fst$PBS >0.51619501 ],
      y=EAS_EUR_Weir.fst$PBS[EAS_EUR_Weir.fst$POS>=(POSs-10000) & EAS_EUR_Weir.fst$POS<=(POSs+10000) & EAS_EUR_Weir.fst$PBS >0.51619501], col="red")




###### SOLUTION PART 1 #######

input="/Users/patriciopezo/Desktop/EMBO_course_2025/input/"
out="/Users/patriciopezo/Desktop/EMBO_course_2025/data_process/"

#I. Read the files with the Fst estimates (AFR_EUR.weir.fst, AFR_EAS.weir.fst and EAS_EUR.weir.fst)
names_header <- c("CHROM","POS","WEIR_AND_COCKERHAM_FST","NUM","DEN")

FST_AFR_EAS <- read.table(paste0(input,"AFR_EAS.weir.fst"), header=F, skip=1, col.names=names_header, fill = TRUE) #582,963 SNPs
FST_AFR_EUR <- read.table(paste0(input,"AFR_EUR.weir.fst"), header=F, skip=1, col.names=names_header, fill = TRUE) #582,964 SNPs
FST_EAS_EUR <- read.table(paste0(input,"EAS_EUR.weir.fst"), header=F, skip=1, col.names=names_header, fill = TRUE) #582,964 SNPs


#II. Remove duplicated positions
FST_AFR_EAS_filter <- FST_AFR_EAS[!duplicated(FST_AFR_EAS$POS),] #582,963 SNPs
FST_AFR_EUR_filter <- FST_AFR_EUR[!duplicated(FST_AFR_EUR$POS),] #582,964 SNPs
FST_EAS_EUR_filter <- FST_EAS_EUR[!duplicated(FST_EAS_EUR$POS),] #582,964 SNPs


#III. Take a look at the weir.fst file

head(FST_AFR_EAS_filter)
head(FST_AFR_EUR_filter)
head(FST_EAS_EUR_filter)


#IV. Exclude NAs position in Fst estimations
FST_AfrEas_data <- FST_AFR_EAS_filter[-which(is.na(FST_AFR_EAS_filter[,3])),] #582,694 SNPs
FST_AfrEur_data <- FST_AFR_EUR_filter[-which(is.na(FST_AFR_EUR_filter[,3])),] #582,363 SNPs
FST_EasEur_data <- FST_EAS_EUR_filter[-which(is.na(FST_EAS_EUR_filter[,3])),] #566,650 SNPs


#V. Overlaping SNPs
#565,779 SNPs

overlap_AfrEas_AfrEur <- FST_AfrEas_data[FST_AfrEas_data$POS %in% FST_AfrEur_data$POS,] 
overlap_AfrEasEur_EasEur <- overlap_AfrEas_AfrEur[overlap_AfrEas_AfrEur$POS %in% FST_EasEur_data$POS,]

FST_AfrEas_data_clean <- FST_AfrEas_data[FST_AfrEas_data$POS %in% overlap_AfrEasEur_EasEur$POS,]
FST_AfrEur_data_clean <- FST_AfrEur_data[FST_AfrEur_data$POS %in% overlap_AfrEasEur_EasEur$POS,]
FST_EasEur_data_clean <- FST_EasEur_data[FST_EasEur_data$POS %in% overlap_AfrEasEur_EasEur$POS,]

#VI. Convert negative values to zero
FST_AfrEas_data_clean[which(FST_AfrEas_data_clean[,3]<0),3] <- 0
FST_AfrEur_data_clean[which(FST_AfrEur_data_clean[,3]<0),3] <- 0
FST_EasEur_data_clean[which(FST_EasEur_data_clean[,3]<0),3] <- 0

#VII. Check if the SNP rs3827760 (pos 109513601) is a candidate for natural selection
#1. Check if the SNP rs3827760, located at position 109513601, is an outlier in the FST distribution for any of the population pairs

POS <- 109513601

FST_AfrEas_data_clean[FST_AfrEas_data_clean$POS==POS,]
##        CHROM       POS WEIR_AND_COCKERHAM_FST      NUM      DEN
## 266093     2 109513601               0.872881 0.762039 0.873016

FST_AfrEur_data_clean[FST_AfrEur_data_clean$POS==POS,]
##        CHROM       POS WEIR_AND_COCKERHAM_FST         NUM       DEN
## 266093     2 109513601             0.00997189 0.000108929 0.0109237

FST_EasEur_data_clean[FST_EasEur_data_clean$POS==POS,]
##        CHROM       POS WEIR_AND_COCKERHAM_FST      NUM      DEN
## 266093     2 109513601               0.859066 0.743056 0.864958


#2. Check which quartile percentile the rs3827760 distribution fall in each analyzed population pair?
FST_AfrEas_distr <- sort(FST_AfrEas_data_clean[,3])
FST_AfrEur_distr <- sort(FST_AfrEur_data_clean[,3])
FST_EasEur_distr <- sort(FST_EasEur_data_clean[,3])

FST_AfrEas_distrQT <- quantile(FST_AfrEas_distr, c(0.01, 0.05, 0.1, .25, .50,  .75, .90, 0.95, .99))

FST_AfrEas_distrQT
##           1%           5%          10%          25%          50%          75% 
## 0.0000000000 0.0000079206 0.0031241100 0.0262112250 0.1056270000 0.2222135000 
##          90%          95%          99% 
## 0.3675160000 0.4685402500 0.6521344200



FST_AfrEur_distrQT <- quantile(FST_AfrEur_distr, c(0.01, 0.05, 0.1, .25, .50,  .75, .90, 0.95, .99))
FST_AfrEur_distrQT
##          1%          5%         10%         25%         50%         75% 
## 0.000000000 0.000000000 0.002228373 0.018979525 0.081086300 0.186617500 
##         90%         95%         99% 
## 0.308151200 0.395242000 0.551998730



FST_EasEur_distrQT <- quantile(FST_EasEur_distr, c(0.01, 0.05, 0.1, .25, .50,  .75, .90, 0.95, .99))
FST_EasEur_distrQT
##           1%           5%          10%          25%          50%          75% 
## 0.0000000000 0.0000000000 0.0004382527 0.0092101950 0.0476573500 0.1294852500 
##          90%          95%          99% 
## 0.2363264000 0.3175715000 0.4790830000




#3. Ploting FST values in a 10,000 base pair region adjacent to the SNP at position 109513601. Highlight the SNPs that are outliers in the 95th percentile in each population pair.


###### Part 2 #######

setwd("../EMBO Popgen/embo_popgen_2025_DCR/Me/Day6/EMBO_Practical_Course_2024-main")

library(rehh)

Chr2_EDAR_CHS<-data2haplohh(hap_file = "Chr2_EDAR_CHS_500K.recode.vcf",polarize_vcf = FALSE)
Chr2_EDAR_LWK<-data2haplohh(hap_file = "Chr2_EDAR_LWK_500K.recode.vcf",polarize_vcf = FALSE)

EHH_LWK<-calc_ehh(Chr2_EDAR_LWK,mrk ="rs3827760" )
EHH_CHS<-calc_ehh(Chr2_EDAR_CHS,mrk ="rs3827760" )

plot(EHH_CHS)
plot(EHH_LWK)

plot(EHH_LWK[[3]]$POSITION,EHH_LWK[[3]]$EHH_A)
plot(EHH_CHS[[3]]$POSITION,EHH_CHS[[3]]$EHH_A)

FUR_LWK<-calc_furcation(Chr2_EDAR_LWK,mrk ="rs3827760" )
FUR_CHS<-calc_furcation(Chr2_EDAR_CHS,mrk ="rs3827760" )

plot(FUR_LWK)
plot(FUR_CHS)


EHH2_CHS<-scan_hh(Chr2_EDAR_CHS)
EHH2_LWK<-scan_hh(Chr2_EDAR_LWK)

EHH2_CHS["rs3827760",]
EHH2_LWK["rs3827760",]


IHS_CHS<-ihh2ihs(scan = EHH2_CHS,freqbin = 0.01,min_maf = 0.02)
IHS_LWK<-ihh2ihs(scan = EHH2_LWK,freqbin = 0.01,min_maf = 0.02)

IHS_CHS$ihs["rs3827760",]
IHS_LWK$ihs["rs3827760",]

plot(IHS_CHS$ihs$POSITION,IHS_CHS$ihs$IHS)


slideFunct <- function(data, window, step)
  { 
  total <- length(data) 
  spots <- seq(from=1, to = (total - window + 1), by = step)
  result<- vector(length = length(spots)) 
  for(i in 1:length(spots)){ 
   result[i] <- mean(abs(data[spots[i]:(spots[i] + window + 1)]),na.rm=TRUE )
   }
  return(result) 
  }


plot(slideFunct(IHS_CHS$ihs$IHS,100,5))
plot(slideFunct(IHS_LWK$ihs$IHS,100,5))
Slid<-slideFunct(IHS_CHS$ihs$IHS,50,40)

slidePos <- function(data, window, step)
{ 
  total <- length(data) 
  spots <- seq(from = 1, to = (total - window + 1), by = step) 
  result <- vector(length = length(spots)) 
  for(i in 1:length(spots))
  { 
    result[i] <- data[spots[i]] 
  } 
  return(result) 
  }

Pos<-slidePos(IHS_CHS$ihs$POSITION,50,40)

#pos (813 of 2137)
#109513601
Table<-as.data.frame(cbind(Slid,Pos))

target_pos <- 109513601  #
# Calculate the absolute differences

# Filter the data frame to include only Pos values smaller than the target
filtered_data <- Table[Table$Pos< target_pos, ]

# Calculate the absolute differences
differences <- abs(filtered_data$Pos - target_pos)

# Find the index of the minimum difference
closest_index <- which.min(differences)

# Retrieve the corresponding Slid value
closest_slid <- filtered_data$Slid[closest_index]

# Print the result
print(closest_slid)

quantile(na.omit(Table$Slid),c(.75,.95,.99))


Slid<-slideFunct(IHS_LWK$ihs$IHS,50,40)

Pos<-slidePos(IHS_LWK$ihs$POSITION,50,40)

#pos (813 of 2137)
#109513601
Table<-as.data.frame(cbind(Slid,Pos))

target_pos <- 109513601  #
# Calculate the absolute differences
differences <- abs(Table$Pos - target_pos)

# Find the index of the minimum difference
closest_index <- which.min(differences)

# Retrieve the corresponding Slid value
closest_slid <- Table$Slid[closest_index]

# Print the result
print(closest_slid)


manhattanplot(IHS_LWK,mrk="rs3827760")



