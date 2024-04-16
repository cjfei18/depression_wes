library(bhr)
library(tidyr)
library(data.table)

pheno_name <- 'PHQ4'
pheno_list <- c(0)
BHR_input_path<-paste0("/home1/Huashan1/lzy_data/WES_Depression/bhr/",pheno_name,"_plof.txt")
BHR_input<-as.data.frame(fread(BHR_input_path))
baseline <- read.table("/home1/Huashan1/wbs_data/software/BHR/data/ms_baseline_oe5.txt")

BHR_input_rare<-subset(BHR_input,BHR_input$AF<1e-4 & BHR_input$AF>=1e-5)
pheno_name_bin1<-paste(pheno_name,"bin1_plof",sep="_")
assign(pheno_name_bin1,BHR_input_rare)
pheno_list<-c(pheno_list,pheno_name_bin1)

BHR_input_rare<-subset(BHR_input,BHR_input$AF<1e-3 & BHR_input$AF>=1e-4)
pheno_name_bin2<-paste(pheno_name,"bin2_plof",sep="_")
assign(pheno_name_bin2,BHR_input_rare)
pheno_list<-c(pheno_list,pheno_name_bin2)

BHR_input_rare<-subset(BHR_input,BHR_input$AF<1e-2 & BHR_input$AF>=1e-3)
pheno_name_bin3<-paste(pheno_name,"bin3_plof",sep="_")
assign(pheno_name_bin3,BHR_input_rare)
pheno_list<-c(pheno_list,pheno_name_bin3)

###missense
BHR_input_path<-paste0("/home1/Huashan1/lzy_data/WES_Depression/bhr/",pheno_name,"_missense.txt")
BHR_input<-as.data.frame(fread(BHR_input_path))
baseline <- read.table("/home1/Huashan1/wbs_data/software/BHR/data/ms_baseline_oe5.txt")

BHR_input_rare<-subset(BHR_input,BHR_input$AF<1e-4 & BHR_input$AF>=1e-5)
pheno_name_bin1<-paste(pheno_name,"bin1_missense",sep="_")
assign(pheno_name_bin1,BHR_input_rare)
pheno_list<-c(pheno_list,pheno_name_bin1)

BHR_input_rare<-subset(BHR_input,BHR_input$AF<1e-3 & BHR_input$AF>=1e-4)
pheno_name_bin2<-paste(pheno_name,"bin2_missense",sep="_")
assign(pheno_name_bin2,BHR_input_rare)
pheno_list<-c(pheno_list,pheno_name_bin2)

BHR_input_rare<-subset(BHR_input,BHR_input$AF<1e-2 & BHR_input$AF>=1e-3)
pheno_name_bin3<-paste(pheno_name,"bin3_missense",sep="_")
assign(pheno_name_bin3,BHR_input_rare)
pheno_list<-c(pheno_list,pheno_name_bin3)

pheno_list<-pheno_list[-1]


result_aggregate<-data.frame(pheno_name=as.character(0), aggregate_h2=as.numeric(0), aggregate_h2_se=as.numeric(0), maf_func_bins=as.character(0))     #创建一个数据框用于存储结果

##all
BHR_aggregate<-BHR(mode = "aggregate", 
ss_list = list(get(paste(pheno_name,"bin1_plof",sep = "_")),get(paste(pheno_name,"bin2_plof",sep = "_")),get(paste(pheno_name,"bin3_plof",sep = "_")),
get(paste(pheno_name,"bin1_missense",sep = "_")),get(paste(pheno_name,"bin2_missense",sep = "_")),get(paste(pheno_name,"bin3_missense",sep = "_"))),
trait_list = list(pheno_name),
annotations = list(baseline))
aggregate_h2<-BHR_aggregate$aggregated_mixed_model_h2
aggregate_h2_se<-BHR_aggregate$aggregated_mixed_model_h2se
result<-cbind(pheno_name,aggregate_h2,aggregate_h2_se,"plof+missense")
colnames(result)<-c("pheno_name","aggregate_h2","aggregate_h2_se","maf_func_bins")
result_aggregate<-rbind(result_aggregate,result)

##lof
BHR_aggregate<-BHR(mode = "aggregate", 
ss_list = list(get(paste(pheno_name,"bin1_plof",sep = "_")),get(paste(pheno_name,"bin2_plof",sep = "_")),get(paste(pheno_name,"bin3_plof",sep = "_"))),
trait_list = list(pheno_name),
annotations = list(baseline))
aggregate_h2<-BHR_aggregate$aggregated_mixed_model_h2
aggregate_h2_se<-BHR_aggregate$aggregated_mixed_model_h2se
result<-cbind(pheno_name,aggregate_h2,aggregate_h2_se,"plof")
colnames(result)<-c("pheno_name","aggregate_h2","aggregate_h2_se","maf_func_bins")
result_aggregate<-rbind(result_aggregate,result)

##missense
BHR_aggregate<-BHR(mode = "aggregate", 
ss_list = list(get(paste(pheno_name,"bin1_missense",sep = "_")),get(paste(pheno_name,"bin2_missense",sep = "_")),get(paste(pheno_name,"bin3_missense",sep = "_"))),
trait_list = list(pheno_name),
annotations = list(baseline))
aggregate_h2<-BHR_aggregate$aggregated_mixed_model_h2
aggregate_h2_se<-BHR_aggregate$aggregated_mixed_model_h2se
result<-cbind(pheno_name,aggregate_h2,aggregate_h2_se,"missense")
colnames(result)<-c("pheno_name","aggregate_h2","aggregate_h2_se","maf_func_bins")
result_aggregate<-rbind(result_aggregate,result)

result_aggregate<-result_aggregate[-1,]
write_path<-"/home1/Huashan1/lzy_data/WES_Depression/bhr/phq4_aggregate.txt"
write.table(result_aggregate,write_path,sep = "\t",row.names = F,quote = F)
