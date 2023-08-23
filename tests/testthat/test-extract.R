piggyback::pb_download(repo = "soroushmdg/gwid",
                       tag = "v0.0.1",
                       dest = tempdir())
ibd_data_file <- paste0(tempdir(),"//chr3.ibd")
genome_data_file <- paste0(tempdir(),"//chr3.gds")
phase_data_file <- paste0(tempdir(),"//chr3.vcf")
case_control_data_file <- paste0(tempdir(),"//case-cont-RA.withmap.Rda")

case_control <- gwid::case_control(case_control_rda = case_control_data_file)
snp_data_gds <- gwid::build_gwas(gds_data = genome_data_file,caco = case_control,gwas_generator = TRUE)
haplotype_data <- gwid::build_phase(phased_vcf = phase_data_file,caco = case_control)
ibd_data <- gwid::build_gwid(ibd_data = ibd_data_file,gwas = snp_data_gds)


test_that("extract section works", {
  expect_s3_class(extract(ibd_data,snp_start = 118926294,snp_end = 122613594 ),"result_snps")

})
