// 253 | LANDIS-R: update only changed cohorts directly in compact CSR form.
// EXPERIMENTAL, never installed in the landisR repository.
// Preserve raw CSR state layout while matching the canonical R state-key semantics.
// Stable state ordering, exact cohort biomass values, no RNG usage.
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <cstring>
using namespace Rcpp;

namespace {
struct Cohort253 {
  int species;
  double age;
  double biomass;
};
struct State253 {
  int eco;
  std::vector<Cohort253> cohorts;
};
template<typename T> void bytes253(std::string &s, const T &x) {
  const char *p=reinterpret_cast<const char *>(&x);
  s.append(p,sizeof(T));
}
// The legacy R key is computed by landisR_state_key(), which first calls
// landisR_normalize_cohort_state(): remove zero-biomass rows and order
// by canonical species index and then increasing age.  IMPORTANT: these
// transformations apply to the KEY only for EXISTING CSR states.  The raw
// CSR cohort row order must stay unchanged when those states are reused,
// because the reference keeps it until bs52_states_list_to_csr_cpp().
std::string key253(const State253 &state) {
  std::vector<Cohort253> normalized;
  normalized.reserve(state.cohorts.size());
  for (const auto &c : state.cohorts) {
    if (c.biomass > 0.0) normalized.push_back(c);
  }
  std::sort(normalized.begin(), normalized.end(),
    [](const Cohort253 &a,const Cohort253 &b) {
      if (a.species != b.species) return a.species < b.species;
      return a.age < b.age;
    });
  std::string key;
  key.reserve(sizeof(int)*2+normalized.size()*(sizeof(int)+2*sizeof(double)));
  bytes253(key,state.eco);
  int count=static_cast<int>(normalized.size());
  bytes253(key,count);
  for (const auto &c:normalized) {
    bytes253(key,c.species);
    bytes253(key,c.age);
    bytes253(key,c.biomass);
  }
  return key;
}
}

// [[Rcpp::export]]
Rcpp::List basefire_apply_csr_253(
    Rcpp::List old_store,
    Rcpp::IntegerVector old_state_eco,
    Rcpp::IntegerVector old_cell_state,
    Rcpp::IntegerVector active_ids,
    Rcpp::IntegerVector active_eco,
    Rcpp::DataFrame replacement,
    Rcpp::IntegerVector clear_cells,
    Rcpp::CharacterVector species_names,
    int species_base) {
  if(species_base!=0 && species_base!=1) stop("253: invalid species base");
  const IntegerVector offsets=as<IntegerVector>(old_store["offset"]);
  const IntegerVector old_species=as<IntegerVector>(old_store["species_id"]);
  const NumericVector old_age=as<NumericVector>(old_store["age"]);
  const NumericVector old_biomass=as<NumericVector>(old_store["biomass"]);
  const int n_state=old_state_eco.size();
  const int n_active=active_ids.size();
  if(offsets.size()!=n_state+1 || old_species.size()!=old_age.size() ||
     old_species.size()!=old_biomass.size() || old_cell_state.size()!=n_active ||
     active_eco.size()!=n_active) stop("253: invalid CSR/metadata dimensions");
  if (offsets.size()<2 || (offsets[0]!=0 && offsets[0]!=1) ||
      offsets[offsets.size()-1]-offsets[0]!=old_species.size())
    stop("253: unsupported CSR offset convention");
  const int offset_base=offsets[0];
  for(int i=1;i<offsets.size();++i)
    if(offsets[i]<offsets[i-1]) stop("253: CSR offsets non-monotone");
  std::vector<State253> states;
  states.reserve(static_cast<size_t>(n_state)+std::min(n_active,4096));
  std::unordered_map<std::string,int> state_map;
  state_map.reserve(static_cast<size_t>(n_state)*2+1024);
  for(int s=0;s<n_state;++s) {
    State253 state;
    state.eco=old_state_eco[s];
    const int from=offsets[s]-offset_base;
    const int to=offsets[s+1]-offset_base;
    state.cohorts.reserve(static_cast<size_t>(to-from));
    for(int i=from;i<to;++i) {
      if(old_species[i]==NA_INTEGER || old_species[i]-species_base<0 ||
         old_species[i]-species_base>=species_names.size())
        stop("253: invalid CSR species code for requested base");
      state.cohorts.push_back({old_species[i],old_age[i],old_biomass[i]});
    }
    std::string key=key253(state);
    states.push_back(std::move(state));
    state_map.emplace(std::move(key),s); // first match, same as match(key, state_keys)
  }
  std::unordered_map<int,int> active_pos;
  active_pos.reserve(static_cast<size_t>(n_active)*2);
  for(int i=0;i<n_active;++i) {
    if(old_cell_state[i]<1 || old_cell_state[i]>n_state)
      stop("253: invalid old cell_state");
    active_pos.emplace(active_ids[i],i);
  }
  std::unordered_map<std::string,int> species_codes;
  species_codes.reserve(static_cast<size_t>(species_names.size())*2);
  for(int i=0;i<species_names.size();++i) {
    if(species_names[i]==NA_STRING) stop("253: missing species name");
    species_codes.emplace(as<std::string>(species_names[i]),i+species_base);
  }
  const IntegerVector repl_cell=as<IntegerVector>(replacement["cell_id"]);
  const CharacterVector repl_species=as<CharacterVector>(replacement["species"]);
  const NumericVector repl_age=as<NumericVector>(replacement["age"]);
  const NumericVector repl_biomass=as<NumericVector>(replacement["biomass"]);
  const int n_replacement=repl_cell.size();
  if(repl_species.size()!=n_replacement || repl_age.size()!=n_replacement ||
     repl_biomass.size()!=n_replacement)
    stop("253: invalid replacement column lengths");

  // Preserve the R semantics of unique(c(clear_cells, cohorts$cell_id)).
  std::vector<int> affected;
  affected.reserve(static_cast<size_t>(n_replacement+clear_cells.size()));
  std::unordered_set<int> seen_cells;
  seen_cells.reserve(static_cast<size_t>(n_replacement+clear_cells.size())*2+1);
  auto add_cell=[&](int id) {
    if(id==NA_INTEGER || active_pos.find(id)==active_pos.end())
      stop("253: inactive or unknown replacement cell");
    if(seen_cells.insert(id).second) affected.push_back(id);
  };
  for(int i=0;i<clear_cells.size();++i) add_cell(clear_cells[i]);
  std::unordered_map<int,std::vector<Cohort253>> grouped;
  grouped.reserve(static_cast<size_t>(n_replacement/4+1));
  std::unordered_set<std::string> unique_triplets;
  unique_triplets.reserve(static_cast<size_t>(n_replacement)*2+1);
  for(int i=0;i<n_replacement;++i) {
    const int id=repl_cell[i];
    add_cell(id);
    if(STRING_ELT(repl_species, i) == NA_STRING) stop("253: missing species");
    auto sp=species_codes.find(as<std::string>(repl_species[i]));
    if(sp==species_codes.end()) stop("253: unknown species");
    const double age=repl_age[i],biomass=repl_biomass[i];
    if(!R_finite(age) || age<1 || age!=std::trunc(age) ||
       !R_finite(biomass) || biomass<0)
       stop("253: invalid age or biomass");
    std::string key;
    bytes253(key,id);bytes253(key,sp->second);bytes253(key,age);
    if(!unique_triplets.insert(key).second) stop("253: duplicated cohort");
    if(biomass>0) grouped[id].push_back({sp->second,age,biomass});
  }

  IntegerVector next_cell_state=clone(old_cell_state);
  for(int id:affected) {
    const int position=active_pos.at(id);
    State253 replacement_state;
    replacement_state.eco=active_eco[position];
    auto g=grouped.find(id);
    if(g!=grouped.end()) replacement_state.cohorts=std::move(g->second);
    std::sort(replacement_state.cohorts.begin(),replacement_state.cohorts.end(),
        [](const Cohort253&a,const Cohort253&b) {
          return a.species<b.species || (a.species==b.species && a.age<b.age);
        });
    auto key=key253(replacement_state);
    auto found=state_map.find(key);
    int state_id;
    if(found!=state_map.end()) state_id=found->second;
    else {
      state_id=static_cast<int>(states.size());
      states.push_back(std::move(replacement_state));
      state_map.emplace(std::move(key),state_id);
    }
    next_cell_state[position]=state_id+1;
  }

  // R's unique(cell_state) visits states in active-cell order. Preserve that.
  std::vector<int> used;
  used.reserve(states.size());
  std::vector<int> remap(states.size(),0);
  for(int i=0;i<n_active;++i) {
    int id=next_cell_state[i]-1;
    if(id<0 || id>=static_cast<int>(states.size()))
      stop("253: invalid updated state index");
    if(remap[id]==0) {
      used.push_back(id);
      remap[id]=static_cast<int>(used.size());
    }
    next_cell_state[i]=remap[id];
  }
  size_t ncoh=0;
  for(int id:used)ncoh+=states[id].cohorts.size();
  IntegerVector out_offsets(used.size()+1);
  IntegerVector out_species(ncoh),out_eco(used.size());
  NumericVector out_age(ncoh),out_biomass(ncoh);
  size_t j=0;
  out_offsets[0]=offset_base;
  for(size_t s=0;s<used.size();++s) {
    const State253 &state=states[used[s]];
    out_eco[s]=state.eco;
    for(const auto &c:state.cohorts) {
      out_species[j]=c.species;
      out_age[j]=c.age;
      out_biomass[j]=c.biomass;
      ++j;
    }
    out_offsets[s+1]=static_cast<int>(j)+offset_base;
  }
  Rcpp::List new_store=clone(old_store);
  new_store["offset"]=out_offsets;
  new_store["species_id"]=out_species;
  new_store["age"]=out_age;
  new_store["biomass"]=out_biomass;
  return Rcpp::List::create(
    Rcpp::Named("states")=new_store,
    Rcpp::Named("state_eco")=out_eco,
    Rcpp::Named("cell_state")=next_cell_state);
}
