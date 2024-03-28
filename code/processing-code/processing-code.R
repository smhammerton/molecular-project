# To install msa on a new computer: 
# install.packages("BiocManager")
# BiocManager::install("msa")


library(msa)

dat <-  readRDS(here::here("data/processed-data/sampled_dat.rds"))


table(dat$host_gender)

# Clean the data ---------------------------------------------------------------
# Wed Mar 13 18:01:34 2024 ------------------------------



clean_dat <- 
  dat |> 
  dplyr::mutate(
    host_gender = dplyr::case_when(
      host_gender == "female" | host_gender == "F" ~ "Female",
      host_gender == "male" | host_gender == "M" ~ "Male"
    )) |> 
  tidyr::separate(host_age, sep = " ", into = c("age", "unit")) |> 
  dplyr::mutate(age = as.numeric(age), 
                age = ifelse(unit == "months", floor(age/12), age)) |> 
  dplyr::select(!unit)


table(clean_dat$host_gender)


saveRDS(clean_dat, here::here("data/processed-data/clean_dat.rds"))


# Sequence alignment -----------------------------------------------------------
# Thu Mar 21 11:29:42 2024 ------------------------------

seqstring <- 
  Biostrings::readRNAStringSet(here::here("data/raw-data/BVBRC_genome_sequence.fasta"))
alignment <- 
  msa::msa(seqstring)
