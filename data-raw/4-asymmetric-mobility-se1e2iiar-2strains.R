se1e2iiar_2strains_pop <- data.table::data.table(
  "location_code" = c("a", "b", "c"),
  "S" = c(1000, 1000, 2000),
  "E1" = c(0, 0, 0),
  "E2" = c(0, 0, 0),
  "I" = c(50, 0, 0),
  "Ia" = c(0, 0, 0),
  "R" = c(0, 0, 0),
  "E1_b" = c(0, 0, 0),
  "E2_b" = c(0, 0, 0),
  "I_b" = c(50, 0, 0),
  "Ia_b" = c(0, 0, 0)
)

temp <- data.table::data.table(
  from = c("a", "a", "b", "b", "c", "c"),
  to = c("b", "c", "a", "c", "a", "b"),
  n = c(50, 10, 10, 10, 10, 10)
)
mobility_matrix <- vector("list", length = 20)
for (i in seq_along(mobility_matrix)) {
  mobility_matrix[[i]] <- data.table::copy(temp)
  data.table::setnames(mobility_matrix[[i]], c("from", "to", "n"))
}

dayEach <- rep(1:5, each = 4)
timeEach <- rep(c(0, 6, 12, 18), 5)
location_codes <- c(rep("a", 5 * 4), rep("b", 5 * 4), rep("c", 5 * 4))
betas <- c(rep(0.6, 5 * 4), rep(0.4, 5 * 2), rep(0.2, 5 * 2), rep(0.7, 5 * 3), rep(0.4, 5))
days <- rep(dayEach, 3)
times <- rep(timeEach, 3)
betas <- data.table("location_code" = location_codes, "day" = days, "time" = times, "beta" = betas)

asymmetric_mobility_se1e2iiar_dummy_se1e2iiar_2strains_pop <- se1e2iiar_2strains_pop
usethis::use_data(asymmetric_mobility_se1e2iiar_dummy_se1e2iiar_2strains_pop, overwrite = TRUE, version = 3, compress = "xz")

asymmetric_mobility_se1e2iiar_dummy_mobility_matrix <- mobility_matrix
usethis::use_data(asymmetric_mobility_se1e2iiar_dummy_mobility_matrix, overwrite = TRUE, version = 3, compress = "xz")

asymmetric_mobility_se1e2iiar_dummy_betas <- betas
usethis::use_data(asymmetric_mobility_se1e2iiar_dummy_betas, overwrite = TRUE, version = 3, compress = "xz")

asymmetric_mobility_se1e2iiar_2strains_dummy_dynamic_seeds <- data.table::data.table(
  "location_code" = c("a"),
  "day" = c(3, 7, 8),
  "n" = c(3, 5, 0),
  "n_b" = c(1, 0, 5)
)
usethis::use_data(asymmetric_mobility_se1e2iiar_2strains_dummy_dynamic_seeds, overwrite = TRUE, version = 3, compress = "xz")
