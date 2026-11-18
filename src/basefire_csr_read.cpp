#include <Rcpp.h>
#include <algorithm>
#include <climits>

namespace {
class CohortStates {
 public:
  explicit CohortStates(const Rcpp::List& states)
    : offsets(states["offset"]), ages(states["age"]) {
    if (offsets.size() == 0 || offsets[0] != 0 ||
        offsets[offsets.size() - 1] != ages.size())
      Rcpp::stop("Invalid Base Fire CSR offsets.");
    for (R_xlen_t i = 1; i < offsets.size(); ++i) {
      if (offsets[i] == NA_INTEGER || offsets[i] < offsets[i - 1])
        Rcpp::stop("Base Fire CSR offsets must be nondecreasing.");
    }
  }

  int state_index(int state) const {
    if (state < 1 || state >= offsets.size())
      Rcpp::stop("Invalid Base Fire CSR cell state.");
    return state - 1;
  }

  int begin(R_xlen_t state) const { return offsets[state]; }
  int end(R_xlen_t state) const { return offsets[state + 1]; }
  R_xlen_t length() const { return offsets.size() - 1; }

  const Rcpp::IntegerVector offsets;
  const Rcpp::NumericVector ages;
};
}

// Only materialize the columns used by Base Fire, in active cell order and
// in the stored cohort order of each state.
// [[Rcpp::export(rng = false)]]
Rcpp::DataFrame basefire_expand_csr(
    const Rcpp::List& states,
    const Rcpp::IntegerVector& cell_state,
    const Rcpp::IntegerVector& cell_id,
    const Rcpp::CharacterVector& species) {
  const CohortStates view(states);
  const Rcpp::IntegerVector species_id = states["species_id"];
  const Rcpp::NumericVector biomass = states["biomass"];
  if (cell_state.size() != cell_id.size() || biomass.size() != view.ages.size() ||
      species_id.size() != view.ages.size())
    Rcpp::stop("Invalid Base Fire CSR cohort dimensions.");
  for (int sp : species_id) {
    if (sp < 0 || sp >= species.size() || species[sp] == NA_STRING)
      Rcpp::stop("Invalid Base Fire CSR species code.");
  }

  R_xlen_t total = 0;
  for (int state : cell_state) {
    const int index = view.state_index(state);
    total += view.end(index) - view.begin(index);
    if (total > INT_MAX) Rcpp::stop("Too many rows for the Base Fire cohort table.");
  }

  Rcpp::IntegerVector out_cell(Rcpp::no_init(total));
  Rcpp::CharacterVector out_species(total);
  Rcpp::NumericVector out_age(Rcpp::no_init(total)), out_biomass(Rcpp::no_init(total));
  R_xlen_t row = 0;
  for (R_xlen_t i = 0; i < cell_id.size(); ++i) {
    const int index = cell_state[i] - 1;
    for (int cohort = view.begin(index); cohort < view.end(index); ++cohort, ++row) {
      out_cell[row] = cell_id[i];
      out_species[row] = species[species_id[cohort]];
      out_age[row] = view.ages[cohort];
      out_biomass[row] = biomass[cohort];
    }
  }
  return Rcpp::DataFrame::create(
    Rcpp::Named("cell_id") = out_cell,
    Rcpp::Named("species_code") = out_species,
    Rcpp::Named("age") = out_age,
    Rcpp::Named("biomass") = out_biomass
  );
}

// Each unique state needs one scan, regardless of how many cells share it.
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector basefire_csr_max_age(
    const Rcpp::List& states,
    const Rcpp::IntegerVector& cell_state) {
  const CohortStates view(states);
  Rcpp::NumericVector max_by_state(view.length());
  for (R_xlen_t i = 0; i < view.length(); ++i) {
    const int start = view.begin(i);
    const int end = view.end(i);
    if (start == end) continue;
    double oldest = view.ages[start];
    for (int cohort = start + 1; cohort < end; ++cohort) {
      const double age = view.ages[cohort];
      if (R_IsNA(age)) {
        oldest = NA_REAL;
        break;
      }
      if (!R_IsNA(oldest))
        oldest = R_IsNaN(age) ? age : std::max(oldest, age);
    }
    max_by_state[i] = oldest;
  }
  Rcpp::NumericVector ages(Rcpp::no_init(cell_state.size()));
  for (R_xlen_t i = 0; i < cell_state.size(); ++i) {
    ages[i] = max_by_state[view.state_index(cell_state[i])];
  }
  return ages;
}
