// 255 - Decode compact CSR into exact R cohort-state data.frames without
// per-state Rcpp DataFrame constructors or R-language loops.
// Based on per-site cohort-state access in LANDIS-II Base Fire C#.
// This function ONLY reads input; no random draws and no mutation.
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <cstring>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP basefire_csr_decode_states_255(SEXP store_sexp, SEXP species_sexp,
                                    int species_base = 0) {
  if (TYPEOF(store_sexp) != VECSXP || TYPEOF(species_sexp) != STRSXP)
    Rcpp::stop("255: CSR store must be a list and species names a character vector");
  if (species_base != 0 && species_base != 1)
    Rcpp::stop("255: invalid species index base");
  Rcpp::List store(store_sexp);
  if (!store.containsElementNamed("offset") ||
      !store.containsElementNamed("species_id") ||
      !store.containsElementNamed("age") ||
      !store.containsElementNamed("biomass"))
    Rcpp::stop("255: required CSR fields absent");
  SEXP offset_sexp = store["offset"];
  SEXP species_codes_sexp = store["species_id"];
  SEXP age_sexp = store["age"];
  SEXP biomass_sexp = store["biomass"];
  if (TYPEOF(offset_sexp)!=INTSXP || TYPEOF(species_codes_sexp)!=INTSXP ||
      TYPEOF(age_sexp)!=REALSXP || TYPEOF(biomass_sexp)!=REALSXP)
    Rcpp::stop("255: unexpected CSR storage types; abort rather than coerce");
  const R_xlen_t noff = XLENGTH(offset_sexp);
  const R_xlen_t ncoh = XLENGTH(species_codes_sexp);
  const R_xlen_t nsp = XLENGTH(species_sexp);
  if (noff < 1 || XLENGTH(age_sexp)!=ncoh || XLENGTH(biomass_sexp)!=ncoh)
    Rcpp::stop("255: inconsistent CSR dimensions");
  const int* offsets = INTEGER(offset_sexp);
  const int* codes = INTEGER(species_codes_sexp);
  const double* ages = REAL(age_sexp);
  const double* biomass = REAL(biomass_sexp);
  if (offsets[0] != 0 || offsets[noff-1] != ncoh)
    Rcpp::stop("255: this kernel supports only zero-based CSR offsets");
  for (R_xlen_t s=0; s+1<noff; ++s) {
    if (offsets[s] < 0 || offsets[s+1] < offsets[s] || offsets[s+1] > ncoh)
      Rcpp::stop("255: nonmonotone CSR offsets");
  }
  for (R_xlen_t k=0; k<ncoh; ++k) {
    int sp = codes[k]-species_base;
    if (sp < 0 || sp >= nsp || STRING_ELT(species_sexp,sp)==NA_STRING)
      Rcpp::stop("255: invalid species code in cohort CSR");
  }
  const R_xlen_t nstate = noff-1;
  // Allocate and reuse the immutable, shared data.frame metadata.
  SEXP output = PROTECT(Rf_allocVector(VECSXP,nstate));
  SEXP colnames = PROTECT(Rf_allocVector(STRSXP,3));
  SET_STRING_ELT(colnames,0,Rf_mkChar("species"));
  SET_STRING_ELT(colnames,1,Rf_mkChar("age"));
  SET_STRING_ELT(colnames,2,Rf_mkChar("biomass"));
  SEXP klass = PROTECT(Rf_mkString("data.frame"));
  // A compact row.names attribute is exactly what base data.frame constructs.
  for (R_xlen_t s=0; s<nstate; ++s) {
    const int start = offsets[s], end = offsets[s+1];
    const int n = end-start;
    SEXP df = PROTECT(Rf_allocVector(VECSXP,3));
    SEXP spp = PROTECT(Rf_allocVector(STRSXP,n));
    SEXP age = PROTECT(Rf_allocVector(REALSXP,n));
    SEXP bio = PROTECT(Rf_allocVector(REALSXP,n));
    for (int k=0; k<n; ++k) {
      const int i=start+k;
      SET_STRING_ELT(spp,k,STRING_ELT(species_sexp,codes[i]-species_base));
      REAL(age)[k]=ages[i];
      REAL(bio)[k]=biomass[i];
    }
    SET_VECTOR_ELT(df,0,spp);
    SET_VECTOR_ELT(df,1,age);
    SET_VECTOR_ELT(df,2,bio);
    Rf_setAttrib(df,R_NamesSymbol,colnames);
    SEXP rn = PROTECT(Rf_allocVector(INTSXP,n>0 ? 2 : 0));
    if (n>0) { INTEGER(rn)[0]=NA_INTEGER; INTEGER(rn)[1]=-n; }
    Rf_setAttrib(df,R_RowNamesSymbol,rn);
    Rf_setAttrib(df,R_ClassSymbol,klass);
    SET_VECTOR_ELT(output,s,df);
    UNPROTECT(5);
  }
  UNPROTECT(3);
  return output;
}
