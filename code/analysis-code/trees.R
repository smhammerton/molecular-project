# Setup-------------------------------------------------------------------------

library(phangorn)
library(msa)
library(dplyr)
library(ape)
library(ggtree)
library(deeptime)
library(treeio)
library(patchwork)


## Data loading-----------------------------------------------------------------

clean_dat <- readRDS(here::here("data/processed-data/clean_dat.rds"))

dates_h1 <- readRDS(here::here("results/dates_h1.rds"))
dates_h3 <- readRDS(here::here("results/dates_h3.rds"))

tree_dat_h1 <- 
  phangorn::as.phyDat(readRDS(here::here("results/h1_nuc_msa.rds")))
tree_dat_h3 <- 
  phangorn::as.phyDat(readRDS(here::here("results/h3_nuc_msa.rds")))


dat_h1 <- 
  clean_dat |>
  dplyr::filter(subtype == "H1N1") |> 
  dplyr::group_by(genome_id) |> 
  dplyr::mutate(strain_2 = ave(strain, strain, 
                               FUN = function(i) paste0(i, '_', seq_along(i))),
                date = ifelse(nchar(collection_date) == 4,
                              paste0(collection_date, "-01-01"),
                              collection_date)) |> 
  dplyr::ungroup()

dat_h3 <- clean_dat |>
  dplyr::filter(subtype == "H3N2") |> 
  dplyr::mutate(strain_2 = ave(strain, strain, 
                               FUN = function(i) paste0(i, '_', seq_along(i))),
                date = ifelse(nchar(collection_date) == 4,
                              paste0(collection_date, "-01-01"),
                              collection_date)) |> 
  dplyr::ungroup()

# Trees-------------------------------------------------------------------------
## H1---------------------------------------------------------------------------
dm  <- dist.ml(tree_dat_h1)
treeUPGMA  <- upgma(dm)
treeNJ  <- NJ(dm)



fun <- function(x) upgma(dist.ml(x))
bs_upgma <- bootstrap.phyDat(tree_dat_h1,  fun)

  



plotBS(treeUPGMA, bs_upgma, main="UPGMA", type = "unrooted")

plot(treeNJ, "unrooted")

h1_models <- 
  tree_dat_h1 |> 
  phangorn::modelTest(model=c("JC", "K80", "SYM"), 
                      multicore = TRUE, mc.cores = 2,
                      control = pml.control(
                        epsilon = 1e-08, maxit = 10, trace = 1))

h1_models

fit_h1 <- pml_bb(tree_dat_h1,
                 model = "SYM+G(4)+I",
                 method="tipdated", tip.dates=as.numeric(dates_h1),
                 rearrangement = "none")
fit_h1

plot(fit_h1, align.tip.label=TRUE)


tree_h1_tib <- 
  tibble::as_tibble(fit_h1$tree)  |> 
  dplyr::rename(strain_2 = label) |> 
  dplyr::left_join(dplyr::select(dat_h1, c(strain_2)),
                                 by = "strain_2") |> 
  dplyr::distinct()

states <- cbind(state.name, state.region, state.division) |> 
  as.data.frame() 

tree_tib <- 
  tree_h1_tib |> 
  dplyr::mutate(label = strain_2) |> 
  tidyr::separate_wider_delim(strain_2, delim = "/", names = c("A", "state.name", "Num", "Year")) |> 
  dplyr::mutate(Year = stringr::str_sub(Year, 1L, -3L)) |> 
  dplyr::mutate(
    state.name = dplyr::case_when(
      state.name == "Baltimore" ~ "Maryland",
      state.name == "Boston" ~ "Massachusetts",
      state.name == "Houston" ~ "Texas",
      state.name == "Santa Clara" ~ "California",
      state.name == "Gainesville" ~ "Florida",
      TRUE ~ state.name
    )) |> 
  dplyr::left_join(states, by = "state.name") |> 
  dplyr::mutate(
    Region = dplyr::case_when(
    state.region == 1 ~ "Northeast",
    state.region == 2 ~ "South",
    state.region == 3 ~ "North Central",
    state.region == 4 ~ "West",
    TRUE ~ NA
  )) |> 
  treeio::as.treedata() |> 
  root(1)



strat_tree_h1n1 <- 
  tree_tib |> 
  ggtree::ggtree(mrsd = "2017-12-18") +
  theme_tree2() +
  geom_tippoint(aes(color = Region)) +
  ggplot2::labs(title = "H1N1") +
  ggplot2::theme(legend.position = "right") +
  ggplot2::theme(legend.position = "bottom")


strat_tree_h1n1

plot(fit_h1$tree) 

ape::write.tree(fit_h1$tree, here::here("results/h1.tree"))

h1n1_tree <-
  fit_h1$tree |>
  ggtree() + 
  theme_tree() +
  geom_treescale(x=0, y=45, fontsize=4, linesize=2, offset=2, width=150) +
  geom_tippoint() +
  ggplot2::labs(title = "H1N1") 

h1n1_tree

h1tree <- here::here("results/h1.tree")

plot(fit_h1$tree)

## H3---------------------------------------------------------------------------
dm_h3  <- dist.ml(tree_dat_h3)
treeUPGMA_h3  <- upgma(dm_h3)
treeNJ_h3  <- NJ(dm_h3)


bs_upgma_h3 <- bootstrap.phyDat(tree_dat_h3,  fun)


h3_models <- 
  tree_dat_h3 |> 
  phangorn::modelTest(model=c("JC", "K80", "SYM"), 
                      multicore = TRUE, mc.cores = 2,
                      control = pml.control(
                        epsilon = 1e-08, maxit = 10, trace = 1))

h3_models

fit_h3 <- pml_bb(tree_dat_h3,
                 model = "SYM+G(4)+I",
                 method="tipdated", tip.dates=as.numeric(dates_h3),
                 rearrangement = "none")
fit_h3

saveRDS(fit_h3, here::here("results/fit_h3.rds"))

fit_h3 <- readRDS(here::here("results/fit_h3.rds"))

fit_h3$tree |> ggtree()

h3n2_tree <- 
  fit_h3$tree |>
  ggtree() +
  theme_tree() +
  geom_treescale(x=0, y=45, fontsize=4, linesize=2, offset=2, width=150) +
  geom_tippoint() +
  ggplot2::labs(title = "H3N2")

tree_h3_tib <- 
  tibble::as_tibble(fit_h3$tree)  |> 
  dplyr::rename(strain_2 = label) |> 
  dplyr::left_join(dplyr::select(dat_h3, c(strain_2)),
                   by = "strain_2") |> 
  dplyr::distinct()


tree_tib_h3 <- 
  tree_h3_tib |> 
  dplyr::mutate(label = strain_2) |> 
  tidyr::separate_wider_delim(strain_2, delim = "/", names = c("A", "state.name", "Num", "Year")) |> 
  dplyr::mutate(Year = stringr::str_sub(Year, 1L, -3L)) |> 
  dplyr::mutate(
    state.name = dplyr::case_when(
      state.name == "Baltimore" ~ "Maryland",
      state.name == "Boston" ~ "Massachusetts",
      state.name == "Houston" ~ "Texas",
      state.name == "Santa Clara" ~ "California",
      state.name == "Gainesville" ~ "Florida",
      TRUE ~ state.name
    )) |> 
  dplyr::left_join(states, by = "state.name") |> 
  dplyr::mutate(
    Region = dplyr::case_when(
      state.region == 1 ~ "Northeast",
      state.region == 2 ~ "South",
      state.region == 3 ~ "North Central",
      state.region == 4 ~ "West",
      TRUE ~ NA
    )) |> 
  treeio::as.treedata() |> 
  root(1)


strat_tree_h3n2 <- 
  tree_tib_h3 |> 
  ggtree::ggtree(mrsd = "2017-12-30") +
  theme_tree2() +
  geom_tippoint(aes(color = Region)) +
  ggplot2::labs(title = "H3N2") +
  ggplot2::theme(legend.position = "right") +
  ggplot2::theme(legend.position = "bottom")


strat_tree_h3n2
# BS and Cons-------------------------------------------------------------------

h1_sums <- 
  purrr::map(bs_upgma, "edge.length") |> 
  purrr::map_dbl(sum) 

h1_res <- 
  tibble::tribble(
    ~subtype, ~mean, ~lwr, ~upr, 
    "H1", mean(h1_sums), quantile(h1_sums, probs = 0.05), quantile(h1_sums, probs = 0.95))

h3_sums <- 
  purrr::map(bs_upgma_h3, "edge.length") |> 
  purrr::map_dbl(sum) 

h3_res <- 
  tibble::tribble(
    ~subtype, ~mean, ~lwr, ~upr,
    "H3", mean(h3_sums), quantile(h3_sums, probs = 0.05), quantile(h3_sums, probs = 0.95))

res <-
  rbind(h1_res, h3_res)


contrasts_df <- 
  cbind(h1_sums, h3_sums) |>
  data.frame() |> 
  dplyr::mutate(contrast = h1_sums - h3_sums) 

con_res <- 
  tibble::tribble(
    ~subtype, ~mean, ~lwr, ~upr,
    "Contrast", mean(contrasts_df$contrast), quantile(contrasts_df$contrast, probs = 0.05), quantile(contrasts_df$contrast, probs = 0.95))
res <-
  rbind(h1_res, h3_res, con_res)

res
saveRDS(res, here::here("results/tables/res.rds"))

# Figs--------------------------------------------------------------------------

dist_fig <- 
  contrasts_df |> 
  tidyr::pivot_longer(cols = c(h1_sums, h3_sums), names_to = "Subtype") |> 
  dplyr::mutate(Subtype = ifelse(Subtype == "h1_sums", "H1N1", "H3N2")) |> 
  ggplot2::ggplot() +
  ggplot2::aes(x = value, color = Subtype, group = Subtype, fill = Subtype) +
  ggplot2::geom_histogram(alpha = 0.75) +
  ggplot2::theme_minimal() +
  ggplot2::labs(title = "Bootstrap distributions of total distances",
                x = "Total distance", y = "Count")

ggplot2::ggsave(here::here("results/figures/dist_fig.png"),
                plot = dist_fig,
                height = 4, width = 4.5)

trees <- h1n1_tree + h3n2_tree

ggplot2::ggsave(here::here("results/figures/trees.png"),
                plot = trees,
                width = 7.5, height = 5)

strat_trees <- strat_tree_h1n1 + strat_tree_h3n2 + 
  patchwork::plot_layout(guides = "collect") & theme(legend.position = 'bottom')
strat_trees

ggplot2::ggsave(here::here("results/figures/strat_trees.png"),
                plot = strat_trees,
                width = 7.5, height = 6)

# Tabs--------------------------------------------------------------------------
h1_mods <- 
  h1_models |> 
  dplyr::mutate(Subtype = "H1N1")

h3_mods <- 
  h3_models |> 
  dplyr::mutate(Subtype = "H3N2")

mods_tab <- 
  rbind(h1_mods, h3_mods) |> 
  gt::gt(groupname_col = "Subtype")  

saveRDS(mods_tab,
        here::here("results/tables/mods_tab.rds"))
