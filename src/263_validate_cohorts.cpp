// 263: Fast, conservative validator for Base Fire cohorts.
// Returns TRUE only for a subset that is *provably* accepted by the original
// .basefire_validate_cohort_table. FALSE always falls back to the original R
// validator, thereby retaining its error behavior on all unusual inputs.
// No RNG, no mutation of package objects, no processing beyond the input.
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_set>
#include <cmath>
#include <climits>
#include <cstring>
using namespace Rcpp;

static bool id263(SEXP x, R_xlen_t i, int& value) {
  if (TYPEOF(x) == INTSXP) {
    int z=INTEGER(x)[i];
    if (z==NA_INTEGER) return false;
    value=z; return true;
  }
  if (TYPEOF(x) == REALSXP) {
    double z=REAL(x)[i];
    if (!std::isfinite(z) || z<static_cast<double>(INT_MIN) ||
        z>static_cast<double>(INT_MAX) || z!=std::trunc(z)) return false;
    value=static_cast<int>(z); return true;
  }
  return false;
}

// [[Rcpp::export(rng = false)]]
bool basefire_valid_cohorts_fast_263(SEXP cohorts, SEXP active_cell_ids) {
  if (TYPEOF(cohorts)!=VECSXP || !Rf_inherits(cohorts,"data.frame")) return false;
  SEXP nm=Rf_getAttrib(cohorts,R_NamesSymbol);
  if (TYPEOF(nm)!=STRSXP || Rf_xlength(nm)!=Rf_xlength(cohorts)) return false;
  int ic=-1,is=-1,ia=-1;
  for (R_xlen_t k=0;k<Rf_xlength(nm);++k) {
    if (STRING_ELT(nm,k)==NA_STRING) return false;
    const char* s=CHAR(STRING_ELT(nm,k));
    if (strcmp(s,"cell_id")==0) ic=static_cast<int>(k);
    if (strcmp(s,"species_code")==0) is=static_cast<int>(k);
    if (strcmp(s,"age")==0) ia=static_cast<int>(k);
  }
  if (ic<0 || is<0 || ia<0) return false;
  SEXP cell=VECTOR_ELT(cohorts,ic), sp=VECTOR_ELT(cohorts,is),
       age=VECTOR_ELT(cohorts,ia);
  const R_xlen_t n=Rf_xlength(cell);
  // A data.frame is a VECSXP (list of columns), normally without dim.
  // Rf_nrows() would return its NUMBER OF COLUMNS, not its row count.
  // Check row.names instead, accepting both R's compact c(NA, -n)
  // encoding and its expanded row-name representation.
  SEXP rn = Rf_getAttrib(cohorts, R_RowNamesSymbol);
  if (Rf_isNull(rn)) return false;
  bool row_count_ok = false;
  if (TYPEOF(rn)==INTSXP && Rf_xlength(rn)==2 &&
      INTEGER(rn)[0]==NA_INTEGER && INTEGER(rn)[1]<=0) {
    // Cast before negating, including the pathological INT_MIN case.
    const long long nr = -static_cast<long long>(INTEGER(rn)[1]);
    row_count_ok = (nr == static_cast<long long>(n));
  } else {
    row_count_ok = (Rf_xlength(rn)==n);
  }
  if (!row_count_ok) return false;
  if (Rf_xlength(sp)!=n || Rf_xlength(age)!=n ||
      !(TYPEOF(cell)==INTSXP || TYPEOF(cell)==REALSXP) ||
      TYPEOF(sp)!=STRSXP ||
      !(TYPEOF(age)==INTSXP || TYPEOF(age)==REALSXP) ||
      !(TYPEOF(active_cell_ids)==INTSXP || TYPEOF(active_cell_ids)==REALSXP))
    return false;
  // An exotic class may have coercion/validation semantics different from
  // the plain atomic vector types. In this case let the original R code run.
  if (!Rf_isNull(Rf_getAttrib(cell,R_ClassSymbol)) ||
      !Rf_isNull(Rf_getAttrib(sp,R_ClassSymbol)) ||
      !Rf_isNull(Rf_getAttrib(age,R_ClassSymbol)) ||
      !Rf_isNull(Rf_getAttrib(active_cell_ids,R_ClassSymbol))) return false;
  if (n==0) return true; // Original validator returns after column checks.
  std::unordered_set<int> active;
  active.reserve(static_cast<size_t>(Rf_xlength(active_cell_ids))*2u+1u);
  for (R_xlen_t j=0;j<Rf_xlength(active_cell_ids);++j) {
    int z;
    if(!id263(active_cell_ids,j,z)) return false;
    active.insert(z);
  }
  for (R_xlen_t j=0;j<n;++j) {
    int z;
    if (!id263(cell,j,z) || active.find(z)==active.end()) return false;
    SEXP sj=STRING_ELT(sp,j);
    if(sj==NA_STRING) return false;
    const char* p=CHAR(sj);
    if(!*p) return false;
    // Conservative: unusual/whitespace species labels use original validator.
    for (; *p; ++p) {
      unsigned char b=static_cast<unsigned char>(*p);
      if (b>127 || b==' ' || b=='\t' || b=='\n' || b=='\r' ||
          b=='\v' || b=='\f') return false;
    }
    const double a=(TYPEOF(age)==INTSXP) ?
      static_cast<double>(INTEGER(age)[j]):REAL(age)[j];
    if (!std::isfinite(a) || a<0.0) return false;
  }
  return true;
}
