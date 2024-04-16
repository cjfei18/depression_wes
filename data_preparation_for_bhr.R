library(bhr)
library(tidyr)
library(data.table)

pheno_single<-read.csv("/home1/Huashan1/lzy_data/WES_Depression/result/PHQ4_single_result.csv")
pheno_name <- 'PHQ4'
pheno_single <- pheno_single[,2:15]

anno<-as.data.frame(fread("/home1/Huashan1/wbs_data/wbs_wes/lof_missense_synonymous/snpeff/result/ukb_wes_chr_all_with_function_snpeff_chr_all_out_extracted_all_sorted.txt",header=F))
colnames(anno)<-c("MarkerID","symbol","ENSG","consequence")
pheno_single <- merge(pheno_single,anno,by="MarkerID")

###all function
BHR_input<-pheno_single[,c(16,2,3,14,9,7)]
BHR_input$phenotype_key<-pheno_name
colnames(BHR_input)<-c("gene","chromosome","gene_position","N","beta","AF","phenotype_key")
write_path<-paste0("/home1/Huashan1/lzy_data/WES_Depression/bhr/",pheno_name,".txt")
write.table(BHR_input,write_path,sep = "\t",row.names = F,quote = F)

###plof
pheno_single2_plof<-subset(pheno_single,pheno_single2$consequence=="1_pLOFs")
BHR_input_plof<-pheno_single2_plof[,c(16,2,3,14,9,7)]
BHR_input_plof$phenotype_key<-pheno_name
colnames(BHR_input_plof)<-c("gene","chromosome","gene_position","N","beta","AF","phenotype_key")
write_path<-paste0("/home1/Huashan1/lzy_data/WES_Depression/bhr/", pheno_name, "_plof.txt")
write.table(BHR_input_plof,write_path,sep = "\t",row.names = F,quote = F)

###missense
pheno_single2_missense<-subset(pheno_single,pheno_single2$consequence=="2_missense")
BHR_input_missense<-pheno_single2_missense[,c(16,2,3,14,9,7)]
BHR_input_missense$phenotype_key<-pheno_name
colnames(BHR_input_missense)<-c("gene","chromosome","gene_position","N","beta","AF","phenotype_key")
write_path<-paste0("/home1/Huashan1/lzy_data/WES_Depression/bhr/",pheno_name,"_missense.txt")
write.table(BHR_input_missense,write_path,sep = "\t",row.names = F,quote = F)

###synonymous
pheno_single2_synonymous<-subset(pheno_single,pheno_single2$consequence=="3_synonymous")
BHR_input_synonymous<-pheno_single2_synonymous[,c(16,2,3,14,9,7)]
BHR_input_synonymous$phenotype_key<-pheno_name
colnames(BHR_input_synonymous)<-c("gene","chromosome","gene_position","N","beta","AF","phenotype_key")
write_path<-paste0("/home1/Huashan1/lzy_data/WES_Depression/bhr/",pheno_name,"_synonymous.txt")
write.table(BHR_input_synonymous,write_path,sep = "\t",row.names = F,quote = F)
