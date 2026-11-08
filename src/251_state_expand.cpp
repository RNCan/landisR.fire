// 251: expand landisR compact state to the exact public cohort data.frame.
// No change in stochastic behaviour, species, cohort order, or biomass values.
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::DataFrame basefire_expand_cohorts_251(
    Rcpp::List states, Rcpp::IntegerVector cell_state,
    Rcpp::IntegerVector cell_ids, Rcpp::IntegerVector rows,
    Rcpp::IntegerVector cols, Rcpp::CharacterVector ecoregions) {
  R_xlen_t n = cell_ids.size();
  if (cell_state.size()!=n || rows.size()!=n || cols.size()!=n || ecoregions.size()!=n)
    Rcpp::stop("251: metadata lengths do not match active ids");
  R_xlen_t ns = states.size();
  std::vector<R_xlen_t> sizes(static_cast<size_t>(ns));
  std::vector<Rcpp::CharacterVector> spp(static_cast<size_t>(ns));
  std::vector<Rcpp::NumericVector> ages(static_cast<size_t>(ns));
  std::vector<Rcpp::NumericVector> biomass(static_cast<size_t>(ns));
  for (R_xlen_t s=0; s<ns; ++s) {
    Rcpp::DataFrame d = states[s];
    spp[s] = Rcpp::as<Rcpp::CharacterVector>(d["species"]);
    ages[s] = Rcpp::as<Rcpp::NumericVector>(d["age"]);
    biomass[s] = Rcpp::as<Rcpp::NumericVector>(d["biomass"]);
    sizes[s] = spp[s].size();
    if (ages[s].size()!=sizes[s] || biomass[s].size()!=sizes[s])
      Rcpp::stop("251: malformed state cohort table");
  }
  R_xlen_t total=0;
  for (R_xlen_t i=0; i<n; ++i) {
    int s = cell_state[i];
    if (s==NA_INTEGER || s < 1 || static_cast<R_xlen_t>(s)>ns)
      Rcpp::stop("251: invalid cell_state index");
    total += sizes[static_cast<size_t>(s-1)];
  }
  Rcpp::IntegerVector out_cell(total), out_row(total), out_col(total);
  Rcpp::CharacterVector out_eco(total), out_species(total);
  Rcpp::NumericVector out_age(total), out_biomass(total);
  R_xlen_t k=0;
  for (R_xlen_t i=0; i<n; ++i) {
    size_t s=static_cast<size_t>(cell_state[i]-1);
    for (R_xlen_t j=0; j<sizes[s]; ++j,++k) {
      out_cell[k]=cell_ids[i];
      out_row[k]=rows[i];
      out_col[k]=cols[i];
      out_eco[k]=ecoregions[i];
      out_species[k]=spp[s][j];
      out_age[k]=ages[s][j];
      out_biomass[k]=biomass[s][j];
    }
  }
  return Rcpp::DataFrame::create(
    Rcpp::Named("cell_id")=out_cell,
    Rcpp::Named("row")=out_row,
    Rcpp::Named("col")=out_col,
    Rcpp::Named("ecoregion")=out_eco,
    Rcpp::Named("species")=out_species,
    Rcpp::Named("age")=out_age,
    Rcpp::Named("biomass")=out_biomass);
}
