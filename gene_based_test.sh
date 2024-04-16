#!/bin/bash
#SBATCH -J depression_burden
#SBATCH -p DCU
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH -o %j.out
#SBATCH -e %j.err

module load miniconda3/RSAIGE

for i in {1..22};
do
  Rscript /public/home/zyli/UKB_WES_3rd/SAIGE_Rfile/step1_fitNULLGLMM.R \
    --sparseGRMFile=/public/home/zyli/UKB_WES_3rd/GRM/UKB_GRM_relatednessCutoff_0.05_5000_randomMarkersUsed_unrelated_3rd_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx \
    --sparseGRMSampleIDFile=/public/home/zyli/UKB_WES_3rd/GRM/UKB_GRM_relatednessCutoff_0.05_5000_randomMarkersUsed_unrelated_3rd_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
    --plinkFile=/public/home/zyli/UKB_WES_3rd/unrelated_0_0442/ukb_wes_chr${i}_sample_qc_final_unrelated \
    --useSparseGRMtoFitNULL=FALSE   \
    --useSparseGRMforVarRatio=TRUE \
    --phenoFile=/public/home/zyli/UKB_WES_3rd/PHQ4.csv \
    --phenoCol=PHQ4 \
    --covarColList=sex,age,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
    --qCovarColList=sex  \
    --sampleIDColinphenoFile=eid \
    --SampleIDIncludeFile=/public/home/zyli/UKB_WES_3rd/british_white_id.txt \
    --isCovariateOffset=FALSE \
    --traitType=quantitative  \
    --invNormalize=TRUE \
    --isCateVarianceRatio=TRUE \
    --outputPrefix=/public/home/zyli/Depression_WES/PHQ4/step1/PHQ4_3rd_10PC_chr${i} \
    --IsOverwriteVarianceRatioFile=TRUE

  Rscript /public/home/zyli/UKB_WES_3rd/SAIGE_Rfile/step2_SPAtests.R \
    --bedFile=/public/home/zyli/UKB_WES_3rd/unrelated_0_0442/ukb_wes_chr${i}_sample_qc_final_unrelated.bed \
    --bimFile=/public/home/zyli/UKB_WES_3rd/unrelated_0_0442/ukb_wes_chr${i}_sample_qc_final_unrelated.bim \
    --famFile=/public/home/zyli/UKB_WES_3rd/unrelated_0_0442/ukb_wes_chr${i}_sample_qc_final_unrelated.fam \
    --SAIGEOutputFile=/public/home/zyli/Depression_WES/PHQ4/step2_single/PHQ4_3rd_10PC_chr${i}.txt \
    --AlleleOrder=alt-first \
    --minMAF=0 \
    --minMAC=0.5 \
    --GMMATmodelFile=/public/home/zyli/Depression_WES/PHQ4/step1/PHQ4_3rd_10PC_chr${i}.rda \
    --varianceRatioFile=/public/home/zyli/Depression_WES/PHQ4/step1/PHQ4_3rd_10PC_chr${i}.varianceRatio.txt \
    --sparseGRMFile=/public/home/zyli/UKB_WES_3rd/GRM/UKB_GRM_relatednessCutoff_0.05_5000_randomMarkersUsed_unrelated_3rd_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx \
  	--sparseGRMSampleIDFile=/public/home/zyli/UKB_WES_3rd/GRM/UKB_GRM_relatednessCutoff_0.05_5000_randomMarkersUsed_unrelated_3rd_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
  	--groupFile=/public/home/zyli/UKB_WES_3rd/SNPEff/SnpEff_gene_group_chr${i}.txt \
    --annotation_in_groupTest="lof,missense:lof,missense" \
    --maxMAF_in_groupTest=0.00001,0.0001,0.001,0.01 \
    --is_output_markerList_in_groupTest=TRUE \
    --LOCO=FALSE \
    --is_fastTest=FALSE
done

