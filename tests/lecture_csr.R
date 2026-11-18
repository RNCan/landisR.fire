ns <- asNamespace("landisR.fire")
decode <- get("basefire_csr_decode_states_255", ns)
expand_legacy <- get("basefire_expand_cohorts_251", ns)
expand <- get("basefire_expand_csr", ns)
maximum <- get("basefire_csr_max_age", ns)

store <- list(offset = c(0L, 3L, 3L, 5L), species_id = c(1L, 0L, 0L, 1L, 0L),
              age = c(40, 50, 20, 15, 50), biomass = c(2.125, 0, 3, 2, 1))
species <- c("PICE.MAR", "ABIE.BAL")
cells <- c(6L, 4L, 9L, 2L)
for (cell_state in list(c(3L, 2L, 1L, 3L), rep(2L, 4L), rep(1L, 4L))) {
  states <- decode(store, species, 0L)
  legacy <- expand_legacy(states, cell_state, cells, cells, cells, rep("eco", length(cells)))
  expected <- data.frame(cell_id = legacy$cell_id, species_code = legacy$species,
                         age = legacy$age, biomass = legacy$biomass)
  expected_max <- vapply(cells, function(id) max(c(0, legacy$age[legacy$cell_id == id])),
                         numeric(1L))
  stopifnot(identical(expand(store, cell_state, cells, species), expected),
            identical(maximum(store, cell_state), expected_max))
}
stopifnot(inherits(try(maximum(store, 4L), silent = TRUE), "try-error"))
