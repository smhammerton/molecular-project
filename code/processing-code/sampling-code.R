raw_dat <- readr::read_csv(here::here("data/raw-data/BVBRC_genome.csv")) |> 
  janitor::clean_names()

raw_seq <- 
  seqinr::read.fasta(file = 
                       here::here("data/raw-data/BVBRC_genome_sequence.fasta"),
                     as.string = TRUE) |> 
  do.call(what = c) |> 
  tibble::enframe(value = "sequence") |> 
  dplyr::mutate(name = stringr::str_sub(name, start = 6L))

complete_dat <- 
  dplyr::left_join(raw_dat, raw_seq, by = c("gen_bank_accessions" = "name"))

selected_dat <-
  complete_dat |> 
  dplyr::select(genome_id, genome_name, ncbi_taxon_id, taxon_lineage_i_ds,
                taxon_lineage_names, strain, segment, subtype, completion_date,
                bio_project_accession, gen_bank_accessions, sequencing_center,
                isolation_source, isolation_comments, collection_date, 
                collection_year, season, geographic_location, host_gender,
                host_age, passage, other_clinical, comments, sequence) |> 
  dplyr::filter(!grepl("VFFSP", strain))

selected_dat |> 
  dplyr::group_by(subtype, collection_year) |> 
  dplyr::count()

set.seed(42)
sampled_dat <- 
  selected_dat |> 
  dplyr::group_by(subtype, collection_year) |> 
  dplyr::slice_sample(n = 50) |> 
  dplyr::ungroup()

saveRDS(sampled_dat,
        here::here("data/processed-data/sampled_dat.rds"))
