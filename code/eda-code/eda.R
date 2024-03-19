library(patchwork)

# Set themes for ggplots
ggplot2::theme_set(ggplot2::theme_minimal())

dat <-  readRDS(here::here("data/processed-data/clean_dat.rds"))


# Strain occurrence fig
strain_nums <- 
  dat |> 
  dplyr::select(collection_year, strain, subtype) |> 
  dplyr::group_by(strain) |> 
  dplyr::add_count() |> 
  ggplot2::ggplot() +
  ggplot2::aes(x = n, label = ..count..) +
  ggplot2::geom_bar(alpha = 0.75, fill = "darkblue") +
  ggplot2::stat_bin(binwidth=1, geom='text', color='black', 
                    position=position_nudge(y=10)) +
  ggplot2::labs(y = "Number of strains",
                x = "Number of occurances over observation period")
strain_nums

ggplot2::ggsave(here::here("results/figures/strain_nums.png"),
                plot = strain_nums,
                height = 3, width = 4)

strain_nums_fct <- 
  strain_nums + 
  ggplot2::facet_wrap(~subtype)

strain_nums_fct

ggplot2::ggsave(here::here("results/figures/strain_nums_fct.png"),
                plot = strain_nums_fct,
                height = 3, width = 6)


dat |> 
  dplyr::group_by(subtype) |> 
  dplyr::select(geographic_location) |> 
  gtsummary::tbl_summary()

h1n1_loc <- 
  dat |> 
  dplyr::filter(subtype == "H1N1") |> 
  dplyr::select(geographic_location) |> 
  gtsummary::tbl_summary() |> 
  gtsummary::as_flex_table()

h1n1_loc

h3n2_loc <-
  dat |> 
  dplyr::filter(subtype == "H3N2") |> 
  dplyr::select(geographic_location) |> 
  gtsummary::tbl_summary()  

h3n2_loc


# H1N1 geographic locations plots 
h1n1_loc_fig <- 
  dat |> 
  dplyr::filter(subtype == "H1N1") |> 
  ggplot2::ggplot() +
  ggplot2::aes(x = geographic_location) +
  ggplot2::geom_bar(alpha = 0.75, fill = "darkblue", stat = "count") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                     vjust = 0.2)) +
  ggplot2::labs(x = NULL, y = "Number of occurances over observation period",
                title = "H1N1")
h1n1_loc_fig

# H3N2 geographic locations plots 
h3n2_loc_fig <- 
  dat |> 
  dplyr::filter(subtype == "H3N2") |> 
  ggplot2::ggplot() +
  ggplot2::aes(x = geographic_location) +
  ggplot2::geom_bar(alpha = 0.75, fill = "darkblue", stat = "count") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                     vjust = 0.2)) +
  ggplot2::labs(x = NULL, y = "Number of occurances over observation period",
                title = "H3N2")

h3n2_loc_fig

ylab <- h1n1_loc_fig$labels$y
h1n1_loc_fig$labels$y <- h3n2_loc_fig$labels$y <- " "

h1n1_loc_fig / h3n2_loc_fig 
grid::grid.draw(grid::textGrob(ylab, x = 0.02, rot = 90))
ggplot2::ggsave(filename = here::here("results/figures/loc_fig.png"),
                plot = ggplot2::last_plot(),
                height = 8, width = 6)


# Strain timing figs 
collection <- 
  dat |> 
  dplyr::mutate(
    collection_date = lubridate::as_date(collection_date)) |> 
  ggplot2::ggplot()+
  ggplot2::aes(x = collection_date) +
  ggplot2::geom_bar(stat = "count",
                    color = "black") +
  ggplot2::facet_wrap(~subtype, nrow = 2) +
  ggplot2::scale_x_date(date_minor_breaks = "14 days", date_breaks = "1 month") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                     size = 6)) +
  ggplot2::labs(x = "Collection Date", y = "Number of samples")

collection

ggplot2::ggsave(here::here("results/figures/collection_timing.png"),
                plot = collection,
                height = 3, width = 6)
