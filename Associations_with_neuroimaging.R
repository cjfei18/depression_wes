rm(list =ls())

cov <- read.csv('cov.csv')
gene<- read.csv('gene.csv')
len1 <- dim(gene)[2]-1
men<-read.table('UKB_DTI.txt',header=T,sep=',')
len2 <- dim(men)[2]-1
men <- men[,c(len2+1,1:len2)]

result <- data.frame(array(NA, dim=c(len1*len2,7)))
names(result) <- c('gene','category','phenotype','beta','se','p','sample')
count <- 1

for (i in c(1:len1)){
  for(j in c(1:len2)){
    print(count)
    data_gene <- gene[,c(1,i+1)]
    data_men <- men[,c(1,j+1)]
    
    category_name <- 'DTI'
    gene_name <- colnames(data_gene)[2]
    men_name <- colnames(data_men)[2]
    
    names(data_gene)<-c('eid','gene')
    names(data_men)<-c('eid','pheno')
    merge_data <- merge(data_gene,data_men,by='eid')
    merge_data <- merge(merge_data,cov,by='eid')
    merge_data <- na.omit(merge_data)
    
    model <- lm(pheno~gene +age+sex+site1+site2+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data=merge_data)
    summary <- summary(model)
    x <- summary$coefficients
    result[count,1] <- gene_name
    result[count,2] <- category_name
    result[count,3] <- men_name
    result[count,4] <- x[2,1]
    result[count,5] <- x[2,2]
    result[count,6] <- x[2,4]
    result[count,7] <- as.numeric(dim(merge_data)[1])
    count <- count+1
    
  }
}
result$p_bon <- p.adjust(result$p, method = "bonferroni")
result$p_fdr <- p.adjust(result$p, method = "BH")

write.csv(result,'DTI.csv')
