#------FEATURE PRE-FILTERING AND SELECTION FOR PREDICTION OF PATHOLOGICAL STAGE IN PROSTATE CANCER------#

#------Data Pre-Processing------#

#----Clinical Variables----#

#--Reading the file--#
prostate_data  <- read.csv("D:\\SEMESTER 3\\Data\\data.csv", header = TRUE)

#--Removing unwanted columns--#
prostate_data <- prostate_data[,-c(27:39)]

#--Removing row with NAs in all columns--#
prostate_data <- prostate_data[-c(1,2,3,4,5,6,7,8,9,10,11,12,13,232),]
prostate_data <- prostate_data[,-c(21,22)]

#--Recategorization of variables--#
prostate_data$NeoAdjRadTx[is.na(prostate_data$NeoAdjRadTx)] <- "No"
prostate_data$NeoAdjRadTx <- ifelse(prostate_data$NeoAdjRadTx == "EXTERNAL BEAM BRACHY"|prostate_data$NeoAdjRadTx == "EXTERNAL BEAM"|prostate_data$NeoAdjRadTx == "BRACHY"|prostate_data$NeoAdjRadTx == "CRYOTHERAPY EXTERNAL BEAM","Yes","No")


prostate_data$ChemoTx[is.na(prostate_data$ChemoTx)] <- "No"
prostate_data$ChemoTx <- ifelse(prostate_data$ChemoTx == "PostCHEMO"|prostate_data$ChemoTx == "PostCHEMO"|prostate_data$ChemoTx == "Neo_Adjuvant_CHEMO"|prostate_data$ChemoTx == "Neo_Adjuvant_CHEMO_PostCHEMO"|prostate_data$ChemoTx == "Chemo"|prostate_data$ChemoTx == "TAXOL"|prostate_data$ChemoTx == "Taxotere"|prostate_data$ChemoTx == "QS-21"|prostate_data$ChemoTx == "Docetaxel"|prostate_data$ChemoTx == "Navelbine, taxotere"|prostate_data$ChemoTx == "Docetaxel, MDX_010","Yes","No")


prostate_data$HormTx[is.na(prostate_data$HormTx)] <- "No"
prostate_data$HormTx <- ifelse(prostate_data$HormTx == "PostHORM"|prostate_data$HormTx == "Neoadjuvant HORM"|prostate_data$HormTx == "Neoadjuvant HORM_PostHORM","Yes","No")


prostate_data$RadTxType[is.na(prostate_data$RadTxType)] <- "No"
prostate_data$RadTxType <- ifelse(prostate_data$RadTxType == "EXTERNAL BEAM"|prostate_data$RadTxType == "Primary EXTERNAL BEAM"|prostate_data$RadTxType == "Primary BRACHY","Yes","No")


prostate_data$Race <- plyr::mapvalues(prostate_data$Race, from = c("Black Non Hispanic", "White Non Hispanic",
                                                                   "Black Hispanic", "White Hispanic", "Asian",
                                                                   "Unknown"),to =c("2Black","1White","2Black","1White","3Asian","4Unknown"))

prostate_data$ClinT_Stage <- plyr::mapvalues(prostate_data$ClinT_Stage, from = c("T1C","T2","T2A","T2B","T2C","T3","T3A","T3B","T3C","T4"),to =c("Early_clinical","Early_clinical","Early_clinical","Early_clinical","Early_clinical","Advanced_Clinical","Advanced_Clinical","Advanced_Clinical","Advanced_Clinical","Advanced_Clinical"))

prostate_data$PathStage <- plyr::mapvalues(prostate_data$PathStage, from = c("T2A","T2B","T2C","T3A","T3B","T3C","T4"),to =c("T2","T2","T2","T3","T3","T3","T4"))

prostate_data$MetSite <- ifelse(is.na(prostate_data$MetSite),'None', 'Yes')


#--Removing rows with NA values--#
prostate_data <- prostate_data[-c(which(is.na(prostate_data$PathStage))),]
prostate_data <- prostate_data[-c(which(is.na(prostate_data$PathGG1))),]
prostate_data <- prostate_data[-c(which(is.na(prostate_data$ClinT_Stage))),]
prostate_data <- prostate_data[-c(which(is.na(prostate_data$PreTxPSA))),]

#--Correlation between numerical clinical variables--#
pros_corr_data <- prostate_data[,-c(1)]
pros_corr_data <- pros_corr_data[,-c(1,2,3,6,7,8,10:23)]

pros_corr_data$PreDxBxPSA <- as.numeric(pros_corr_data$PreDxBxPSA)
pros_corr_data$DxAge <- as.numeric(pros_corr_data$DxAge)
pros_corr_data$PreTxPSA <- as.numeric(pros_corr_data$PreTxPSA)

corrplot(cor(pros_corr_data),method='number')

#--converting categorical features into factore--#
prostate_data$Type = as.factor(prostate_data$Type)
prostate_data$MetSite = as.factor(prostate_data$MetSite)
prostate_data$Race = as.factor(prostate_data$Race)
prostate_data$ClinT_Stage = as.factor(prostate_data$ClinT_Stage)
prostate_data$RP_Type = as.factor(prostate_data$RP_Type)
prostate_data$SMS = as.factor(prostate_data$SMS)
prostate_data$ECE = as.factor(prostate_data$ECE)
prostate_data$SVI = as.factor(prostate_data$SVI)
prostate_data$LNI = as.factor(prostate_data$LNI)
prostate_data$PathStage = as.factor(prostate_data$PathStage)
prostate_data$NeoAdjRadTx = as.factor(prostate_data$NeoAdjRadTx)
prostate_data$ChemoTx = as.factor(prostate_data$ChemoTx)
prostate_data$HormTx = as.factor(prostate_data$HormTx)
prostate_data$RadTxType = as.factor(prostate_data$RadTxType)
prostate_data$BxGG1 = as.factor(prostate_data$BxGG1)
prostate_data$BxGG2 = as.factor(prostate_data$BxGG2)
prostate_data$BxGGS = as.factor(prostate_data$BxGGS)

#--Converting pathological staging into two class--#
pros_new_dat = prostate_data
pros_new_dat$PathStage <- plyr::mapvalues(pros_new_dat$PathStage, from = c("T2","T3","T4"),to =c("Early_Stage","Advanced_Stage","Advanced_Stage"))

#--Converting numeric factors into numeric--#
pros_new_dat$Type = as.numeric(prostate_data$Type)
pros_new_dat$MetSite = as.numeric(prostate_data$MetSite)
pros_new_dat$Race = as.numeric(prostate_data$Race)
pros_new_dat$ClinT_Stage = as.numeric(prostate_data$ClinT_Stage)
pros_new_dat$RP_Type = as.numeric(prostate_data$RP_Type)
pros_new_dat$SMS = as.numeric(prostate_data$SMS)
pros_new_dat$ECE = as.numeric(prostate_data$ECE)
pros_new_dat$SVI = as.numeric(prostate_data$SVI)
pros_new_dat$LNI = as.numeric(prostate_data$LNI)
pros_new_dat$PathStage = as.factor(pros_new_dat$PathStage)
pros_new_dat$NeoAdjRadTx = as.numeric(pros_new_dat$NeoAdjRadTx)
pros_new_dat$ChemoTx = as.numeric(pros_new_dat$ChemoTx)
pros_new_dat$HormTx = as.numeric(pros_new_dat$HormTx)
pros_new_dat$RadTxType = as.numeric(pros_new_dat$RadTxType)
pros_new_dat$PathGG1 = as.numeric(prostate_data$PathGG1)
pros_new_dat$PathGG2 = as.numeric(prostate_data$PathGG2)
pros_new_dat$PathGGS = as.numeric(prostate_data$PathGGS)
pros_new_dat$BxGG1 = as.numeric(prostate_data$BxGG1)
pros_new_dat$BxGG2 = as.numeric(prostate_data$BxGG2)
pros_new_dat$BxGGS = as.numeric(prostate_data$BxGGS)
na.omit(pros_new_dat)

#--Removing post-RP and highly correlated features--#
pros_new_dat = pros_new_dat[,-c(16,17,18,19,20,22,23,24)]
rownames(pros_new_dat) = pros_new_dat[,1]
pros_new_dat = pros_new_dat[,-c(1)]
pros_new_dat = pros_new_dat[,-c(4)]
pros_new_dat = pros_new_dat[,-c(2)]

#----mRNA Variables----#
data_mRNA <- read.table(file = "D:\\SEMESTER 3\\Data\\MSKCC_PCa_mRNA_data.txt", header = TRUE)

#--taking transpose--#
transpose_mRNA <- t(data_mRNA)

#--Removing gene ID and making gene symbol as column names--#
transpose_mRNA <- transpose_mRNA[-1,]
colnames(transpose_mRNA) <- transpose_mRNA[1,]
transpose_mRNA <- transpose_mRNA[-1,]
transpose_mRNA[1,]
transpose_mRNA[,1]

#--Merging pathStage with mRNA data--#
pros_mRNA_copy <- pros_new_dat
pros_mRNA_stage <- data.frame(rownames(pros_mRNA_copy),pros_mRNA_copy[,c(13)])
colnames(pros_mRNA_stage)=c('sample_id','PathStage')
rownames(pros_mRNA_stage) <- pros_mRNA_stage[,1]
transpose_mRNA <- merge(transpose_mRNA,pros_mRNA_stage, by=0)

#--Removing duplicate sample id column and making sample id as rownames--#
rownames(transpose_mRNA) <- transpose_mRNA[,1]
transpose_mRNA <- transpose_mRNA[,-c(1,26449)]

#--Converting all mRNA columns into numeric--#
transpose_mRNA[,1:26447] <- lapply(transpose_mRNA[,1:26447], function(x) as.numeric(as.character(x)))

#--Finding correlation between mRNA variables and removing highly correlated mRNA variables--#
mrna_correlation <- cor(transpose_mRNA[,1:26447])
mrna_correlation[c(1:5),c(1:5)]
corrplot(cor(mrna_correlation[,highlyCor][c(1:50),c(1:50)]))
mrna_correlation[,highlyCor][c(1:50),c(1:50)]
memory.limit(size=56000)
highlyCor <- caret::findCorrelation(mrna_correlation, 0.9,exact = F)
transpose_mRNA <- transpose_mRNA[,-highlyCor]
na.omit(transpose_mRNA)

#----Combining mRNA and Clinical Features----#
transpose_mrna_copy <- transpose_mRNA
transpose_mrna_copy[,c(18562:18565)]


pros_clinical_mrna_dat <- merge(transpose_mrna_copy,pros_mRNA_copy, by=0)
nrow(pros_clinical_mrna_dat)
ncol(pros_clinical_mrna_dat)

pros_clinical_mrna_dat[,c(18566:18579)]
pros_clinical_mrna_dat[,c(1,18566)]

rownames(pros_clinical_mrna_dat) = pros_clinical_mrna_dat[,1]
pros_clinical_mrna_dat = pros_clinical_mrna_dat[,-c(1,18566)]
pros_clinical_mrna_dat[,c(18565:18577)]
colnames(pros_clinical_mrna_dat)[18577] = 'PathStage'

summary(pros_clinical_mrna_dat$PathStage)
pros_clin_mrna_numeric <- pros_clinical_mrna_dat[,c(1:18564,18567,18571,18577)]
pros_clin_mrna_categorical <- pros_clinical_mrna_dat[,c(18565,18566,18568:18570,18572:18577)]

#--Removing Special Characters in column names--#
pros_clinical_mrna_dat_copy <- pros_clinical_mrna_dat
colnames(pros_clinical_mrna_dat_copy) <- gsub("-", "", colnames(pros_clinical_mrna_dat_copy))
colnames(pros_clinical_mrna_dat_copy) <- gsub("/", "", colnames(pros_clinical_mrna_dat_copy))
colnames(pros_clinical_mrna_dat_copy) <- gsub("\\.", "", colnames(pros_clinical_mrna_dat_copy))


pros_clinical_mrna_dat_copy[,18565] = as.numeric(pros_clinical_mrna_dat_copy[,18565])
pros_clinical_mrna_dat_copy[,18566] = as.numeric(pros_clinical_mrna_dat_copy[,18566])
pros_clinical_mrna_dat_copy[,18572] = as.numeric(pros_clinical_mrna_dat_copy[,18572])
pros_clinical_mrna_dat_copy[,18573] = as.numeric(pros_clinical_mrna_dat_copy[,18573])
pros_clinical_mrna_dat_copy[,18574] = as.numeric(pros_clinical_mrna_dat_copy[,18574])
pros_clinical_mrna_dat_copy[,18575] = as.numeric(pros_clinical_mrna_dat_copy[,18575])
pros_clinical_mrna_dat_copy[,18576] = as.numeric(pros_clinical_mrna_dat_copy[,18576])

#--Combine variables--#
pros_clinical_mrna_dat_copy
