# To install msa on a new computer: 
# install.packages("BiocManager")
# BiocManager::install("msa")


library(msa)
library(dplyr)

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



###
# Multiple Sequence Alignment for H1 strains
# Zane Billings
# 2024-03-30
# Perform a multiple sequence alignment on the sequences, and calculate
# multiple pairwise distance matrices
# H1 and H3 have to be done separately or for an unknown and unfindable
# reason the MSA package stops working.
###

# ---- Setup ----
# Declare package dependencies
box::use(
  readr,
  here,
  msa[...]
)



# separate H1 and H3
dat_h1 <- 
  clean_dat |>
  dplyr::filter(subtype == "H1N1") |> 
  dplyr::group_by(genome_id) |> 
  dplyr::mutate(strain_2 = ave(strain, strain, 
                               FUN = function(i) paste0(i, '_', seq_along(i))),
                date = ifelse(nchar(collection_date) == 4,
                              paste0(collection_date, "-01-01"),
                              collection_date))

dates_h1 <- setNames(lubridate::as_date(dat_h1$date), 
                     dat_h1$strain_2)
head(dates_h1)
saveRDS(dates_h1, here::here("results/dates_h1.rds"))

dat_h3 <- clean_dat |>
  dplyr::filter(subtype == "H3N2") |> 
  dplyr::mutate(strain_2 = ave(strain, strain, 
                               FUN = function(i) paste0(i, '_', seq_along(i))),
                date = ifelse(nchar(collection_date) == 4,
                              paste0(collection_date, "-01-01"),
                              collection_date))

dates_h3 <- setNames(lubridate::as_date(dat_h3$date), 
                     dat_h3$strain_2)
head(dates_h3)
saveRDS(dates_h3, here::here("results/dates_h3.rds"))


# ---- Nucleotide Alignment (H1) ----
h1_nuc_seqs <- with(dat_h1, 
                    rlang::set_names(
                      sequence, 
                      ave(strain, strain, 
                          FUN = function(i) paste0(i, '_', seq_along(i)))
                    )) |>
  # Replace T with U for MSA
  gsub(pattern = "t", replacement = "u")

h1_nuc_msa <-
  h1_nuc_seqs |>
  msa::msa(
    method = "Muscle",
    type = "rna",
    order = "input",
    verbose = TRUE
  )

 saveRDS(h1_nuc_msa, here::here("results/h1_nuc_msa.rds"))


h1_nuc_seqs_aligned <- alignment_to_character(h1_nuc_msa)

View(h1_nuc_seqs_aligned)

dat_seqs_h1 <-
  tibble::tibble(
    genome_id = dat_h1$genome_id,
    nuc_aligned = h1_nuc_seqs_aligned
  )


# Save the alignment to file
readr::write_rds(
  dat_seqs_h1,
  here::here("results", "h1-rna-alignment.Rds")
)


# ---- Nucleotide Alignment (H3) ----
h3_nuc_seqs <- with(dat_h3, 
                    rlang::set_names(
                      sequence, 
                      ave(strain, strain, 
                          FUN = function(i) paste0(i, '_', seq_along(i)))
                    )) |>
  # Replace T with U for MSA
  gsub(pattern = "t", replacement = "u")

h3_nuc_msa <-
  h3_nuc_seqs |>
  msa::msa(
    method = "Muscle",
    type = "rna",
    order = "input",
    verbose = TRUE
  )

saveRDS(h3_nuc_msa, here::here("results/h3_nuc_msa.rds"))


h3_nuc_seqs_aligned <- alignment_to_character(h3_nuc_msa)

View(h3_nuc_seqs_aligned)

dat_seqs_h3 <-
  tibble::tibble(
    genome_id = dat_h3$genome_id,
    nuc_aligned = h3_nuc_seqs_aligned
  )



# Save the alignment to file
readr::write_rds(
  dat_seqs_h3,
  here::here("results", "h3-rna-alignment.Rds")
)




