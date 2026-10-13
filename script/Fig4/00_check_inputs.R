# Run from the github/ repository root.
# Checks links without modifying or copying any input data.
required <- c(
  "class_converted_MF_chrX_260917.tsv.gz","clone_info_260917.tsv",
  "chrX_cgi_total_summary_pattern.tsv","ase_res_female_chrX_final_260919_v2.tsv",
  "allele_info_260918.tsv","final_rna_sample_list.tsv",
  "exp_RSEM_TPM_374clones_261013.tsv","exp_bc_TPM_374clones_261013.tsv",
  "Nature_GTEx_Tukiainene_2017.txt","unmethyl_x_cgis_pattern_manual_v3_260918.tsv",
  "HCs_final_res_260919_v20.tsv","HC50_RandomSample_100_final_res_260919_v20.tsv")
f <- file.path("input",required)
print(data.frame(file=required,exists=file.exists(f),is_link=nzchar(Sys.readlink(f)),
                 target=Sys.readlink(f)),row.names=FALSE)
cat("\nMissing required files:",sum(!file.exists(f)),"\n")
for(d in c("input/tree_F","input/blood_vaf_F")){
  cat(d,":",if(dir.exists(d)) length(list.files(d,recursive=TRUE)) else "MISSING","entries\n")
}
