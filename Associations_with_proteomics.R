rm(list =ls())
library(data.table)
j<-1
p <- fread('ukb_protein.tsv.gz')
p <- p[p$ins_index==0,]
p_dim<- dim(p)[2]-3
p_data <- p[,4:1466]
p_data <- as.data.frame(lapply(p_data,as.numeric))
p_data_name <- colnames(p_data)
p_eid <- p[,2]

snp <- read.csv('gene.csv')
snp_list <- colnames(snp)
snp_list <- snp_list[-1]
cov <- read.csv('cov_Proteomics.csv')

result <- data.frame(array(NA, dim=c(p_dim,7)))
names(result) <- c('gene','category','phenotype','beta','se','p','sample')

count <- 1
for (i in 1:p_dim){
  print(count)
  snp_data <- snp[,c(1,j+1)]
  names(snp_data) <- c('eid','gene')
  pheno <- cbind(p_eid,p_data[i])
  names(pheno) <- c('eid','pheno')
  merge_data <- merge(snp_data,pheno,by='eid')
  merge_data <- na.omit(merge_data)
  merge_data <- merge(merge_data,cov,by='eid')
  merge_data <- na.omit(merge_data)
  
  model <- lm(pheno~gene +age+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
                          site10003+site11001+site11002+site11003+site11004+site11005+
                          site11006+site11007+site11008+site11009+site11010+site11011+
                          site11012+site11013+site11014+site11016+site11017+site11018+
                          site11020+site11021+site11022,
              data=merge_data)
  summary <- summary(model)
  x <- summary$coefficients
  result[count,1] <- snp_list[j]
  result[count,2] <- 'Proteomics'
  result[count,3] <- p_data_name[i]
  result[count,4] <- x[2,1]
  result[count,5] <- x[2,2]
  result[count,6] <- x[2,4]
  result[count,7] <- as.numeric(dim(merge_data)[1])
  count <- count+1
}


result1 <- result[!is.na(result$gene),]
output1 <- paste0('Proteomics_',snp_list[j],'.csv')
write.csv(result1,output1)

