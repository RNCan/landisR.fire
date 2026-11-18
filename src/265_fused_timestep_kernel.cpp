// 247: Experimental fused Base Fire event engine (Rcpp).
// Portable prototype for testing only; does NOT change landisR.fire source.
// Matches native R RNG call order, including sample.int() on selected neighbors.
// Directly computes fire severity, burn chance and cohort mortality on contiguous
// vectors, avoiding repeated R data.frame construction and validation.
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <landisR/mt19937.hpp>
#include <landisR/progress.hpp>
#include <memory>
#include <chrono>
#include <vector>
#include <algorithm>
#include <cmath>
#include <climits>
#include <R_ext/Random.h>
using namespace Rcpp;

// The compatibility stream lives for the entire timestep, including its
// ignition scan and the propagation of every accepted ignition.
struct FireV7Random {
  landis_rng::Stream stream;
  NumericVector mean_size, min_size, max_size;
  double cell_area;

  explicit FireV7Random(const List& config) :
    stream(0), mean_size(config["mean_size"]), min_size(config["min_size"]),
    max_size(config["max_size"]), cell_area(as<double>(config["cell_area"])) {
    stream.restore(as<std::string>(config["state"]));
    stream.draws = landis_rng::Stream::offset(as<double>(config["offset"]));
  }

  int size(int region) {
    const double mean = mean_size[region];
    if (mean <= 0) return 0;
    double u;
    do { u = stream.next_double(); } while (u < 1e-6);
    const double value = -std::log(u) / (1. / mean);
    return static_cast<int>((value < min_size[region] ? min_size[region] :
       value > max_size[region] ? max_size[region] : value) / cell_area);
  }
};

// Stable cohort ranges indexed by active cell position. Rows in the common
// cell-ordered table are their own permutation. Other tables retain input
// row order within each cell and keep the same row IDs for event journals.
struct FireCohortIndex {
  std::vector<int> offsets;
  std::vector<int> row_order;

  void build(int n_cells, const IntegerVector& sites) {
    if (sites.size() > INT_MAX) stop("Too many cohorts for Base Fire row indices");
    offsets.assign(static_cast<size_t>(n_cells) + 1, 0);
    bool in_cell_order = true;
    int previous = 0;
    for (int site : sites) {
      if (site < 0 || site > n_cells) stop("Invalid Base Fire cohort cell position");
      in_cell_order = in_cell_order && site > 0 && site >= previous;
      if (site > 0) ++offsets[site];
      previous = site;
    }
    for (size_t site = 1; site < offsets.size(); ++site) {
      offsets[site] += offsets[site - 1];
    }
    if (in_cell_order) return;

    row_order.resize(offsets.back());
    std::vector<int> next = offsets;
    for (int row = 0; row < sites.size(); ++row) {
      const int site = sites[row];
      if (site > 0) row_order[next[site - 1]++] = row;
    }
  }
};

// 265: One local workspace per fused timestep, reused by every event.
// No global cache. No per-event O(N) scratch initialization in fused mode.
// Legacy standalone event calls retain their historical scratch allocation.
struct FireWorkspace265 {
  int n;
  int generation;
  int events;
  size_t max_queue;
  size_t max_checked;
  std::vector<int> queue_mark;
  std::vector<double> site_severity;
  std::vector<int> queue;
  std::vector<int> checked_order;
  FireCohortIndex cohort_index;
  FireV7Random* v7 = nullptr;

  explicit FireWorkspace265(int n_cells) : n(n_cells), generation(0),
      events(0), max_queue(0), max_checked(0),
      queue_mark(n_cells,0), site_severity(n_cells,NA_REAL) {
    if (n_cells < 1) Rcpp::stop("265: empty workspace");
    queue.reserve(std::min(n_cells,128));
    checked_order.reserve(std::min(n_cells,128));
  }
  void start() {
    if (generation == INT_MAX) {
      // This branch is unreachable for ordinary simulated fire counts, but
      // prevents counter overflow and stale membership for very long runs.
      std::fill(queue_mark.begin(),queue_mark.end(),0);
      generation = 0;
    }
    ++generation;
    ++events;
    queue.clear();
    checked_order.clear();
  }
  void record() {
    max_queue = std::max(max_queue,queue.size());
    max_checked = std::max(max_checked,checked_order.size());
  }
};

// Event-level performance counters (not part of the event return contract).
static double g_event_uniform_draws_249 = 0.;
static double g_event_size_callbacks_249 = 0.;
static double g_event_shuffle_callbacks_249 = 0.;
static double g_event_count_249 = 0.;

// Explicit R RNG lease. This kernel is exported with rng=false.
// If an R callback is needed, the state MUST be published first.
// The destructor also commits on a C++ exception if the lease is held.
class RngLease249 {
 public:
  explicit RngLease249(bool enabled = true) : held_(enabled) { if (held_) GetRNGstate(); }
  ~RngLease249() noexcept { if (held_) PutRNGstate(); }
  RngLease249(const RngLease249&) = delete;
  RngLease249& operator=(const RngLease249&) = delete;

  double uniform() {
    if (!held_) Rcpp::stop("249: RNG lease not held");
    return R::runif(0.0, 1.0);
  }
  void publish() {
    if (!held_) Rcpp::stop("249: RNG lease already published");
    PutRNGstate();
    held_ = false;
  }
  void reload() {
    if (held_) Rcpp::stop("249: RNG lease already held");
    GetRNGstate();
    held_ = true;
  }
 private:
  bool held_;
};

// Small independent RNG interleaving regression before the simulation.
// Uses a callback for an exponential and a callback for sample.int,
// with uniform draws both before and after each callback.
// [[Rcpp::export(rng = false)]]
List basefire_rng_sequence_249(
    int before_n, int middle_n, int after_n,
    Function draw_exp, Function shuffle_fun) {
  if (before_n < 0 || middle_n < 0 || after_n < 0 ||
      before_n > 1000 || middle_n > 1000 || after_n > 1000)
    stop("249: invalid RNG test draw count");
  NumericVector before(before_n), middle(middle_n), after(after_n);
  RngLease249 lease;
  for (int i = 0; i < before_n; ++i) before[i] = lease.uniform();
  lease.publish();
  RObject exp_value(draw_exp());
  lease.reload();
  for (int i = 0; i < middle_n; ++i) middle[i] = lease.uniform();
  lease.publish();
  RObject shuffle_value(shuffle_fun());
  lease.reload();
  for (int i = 0; i < after_n; ++i) after[i] = lease.uniform();
  return List::create(_["before"] = before, _["exp"] = exp_value,
    _["middle"] = middle, _["shuffle"] = shuffle_value,
    _["after"] = after);
}

// [[Rcpp::export(rng = false)]]
List basefire_kernel_stats_249() {
  return List::create(_["events"] = g_event_count_249,
    _["uniform_draws"] = g_event_uniform_draws_249,
    _["size_callbacks"] = g_event_size_callbacks_249,
    _["shuffle_callbacks"] = g_event_shuffle_callbacks_249);
}

// [[Rcpp::export(rng = false)]]
List basefire_event_kernel_247(
    int initiation_cell_index,
    int current_time,
    const IntegerVector& cell_id,
    const IntegerVector& region_index,
    int n_regions,
    const IntegerMatrix& neighbors_full,
    const IntegerMatrix& second_neighbors_full,
    const IntegerVector& full_to_active,
    const LogicalVector& full_active,
    const NumericVector& fire_spread_age,
    const NumericMatrix& fuel_curve,
    const NumericMatrix& wind_curve,
    bool has_wind,
    const NumericVector& last_wind_year,
    const NumericVector& last_fire_year_input,
    const LogicalVector& event_occupied_input,
    const NumericVector& last_severity_input,
    const List& rows_by_site,
    Environment cohort_env,
    const NumericVector& cohort_age,
    const IntegerVector& cohort_species_index,
    const IntegerVector& species_fire_tolerance,
    const NumericVector& species_longevity,
    const NumericVector& damage_max_age,
    const IntegerVector& damage_severity_tolerance_diff,
    Function size_one,
    bool sparse_result_264 = false,
    SEXP workspace_ptr_265 = R_NilValue) {

  const int n = cell_id.size();
  if (n < 1 || initiation_cell_index < 1 || initiation_cell_index > n ||
      region_index.size() != n || fire_spread_age.size() != n ||
      last_fire_year_input.size() != n || event_occupied_input.size() != n ||
      neighbors_full.nrow() != n ||
      neighbors_full.ncol() != 8 || second_neighbors_full.nrow() != n ||
      second_neighbors_full.ncol() != 8 || n_regions < 1 ||
      fuel_curve.nrow() < n_regions || fuel_curve.ncol() != 5 ||
      wind_curve.nrow() < n_regions || wind_curve.ncol() != 5 ||
      (has_wind && last_wind_year.size() != n) ||
      cohort_species_index.size() != cohort_age.size() ||
      species_longevity.size() != species_fire_tolerance.size() ||
      damage_max_age.size() != damage_severity_tolerance_diff.size()) {
    stop("247: Invalid event kernel cache dimensions");
  }

  // The 265 fused path uses the single explicitly-owned workspace passed by
  // the timestep. All scratch arrays are created once for that timestep.
  // The legacy standalone path is intentionally kept functionally unchanged.
  FireWorkspace265* workspace265 = NULL;
  if (sparse_result_264 && workspace_ptr_265 != R_NilValue) {
    Rcpp::XPtr<FireWorkspace265> xp265(workspace_ptr_265);
    workspace265 = xp265.get();
    if (!workspace265 || workspace265->n != n)
      stop("265: workspace incompatible with landscape");
    workspace265->start();
  }
  const FireCohortIndex* cohort_index =
    workspace265 && !workspace265->cohort_index.offsets.empty() ?
      &workspace265->cohort_index : nullptr;
  if (!cohort_index && rows_by_site.size() != n)
    stop("247: Invalid event kernel cohort index dimensions");
  // 249: exact R RNG, but without a callback into stats::runif(1) per draw.
  // GetRNGstate / PutRNGstate bracket C++ draws. Before any callback that
  // itself draws (size_one/rexp and sample.int), flush .Random.seed, call R,
  // then reload its NEW state. Never hold an Rcpp RNGScope across a callback.
  FireV7Random* v7 = workspace265 ? workspace265->v7 : nullptr;
  RngLease249 synced_rng(!v7);
  auto draw_uniform = [&]() -> double {
    ++g_event_uniform_draws_249;
    return v7 ? v7->stream.next_double() : synced_rng.uniform();
  };
  auto draw_size = [&](int region) -> int {
    if (v7) return v7->size(region);
    synced_rng.publish();
    const int size = as<int>(size_one(region + 1));
    ++g_event_size_callbacks_249;
    synced_rng.reload();
    return size;
  };

  // 264: when called by the fused timestep, update its *owned* current-state
  // arrays directly.  Standalone event calls retain the historical clone.
  NumericVector last_fire_year = sparse_result_264 ?
      NumericVector(last_fire_year_input) : clone(last_fire_year_input);
  LogicalVector event_occupied = sparse_result_264 ?
      LogicalVector(event_occupied_input) : clone(event_occupied_input);
  LogicalVector alive = cohort_env["alive"];
  if (alive.size() != cohort_age.size()) {
    stop("247: Deferred cohort alive mask missing or incompatible");
  }

  // Same RNG calls and order as the original R event.
  const int wind_direction_raw = static_cast<int>(draw_uniform() * 8.);
  const int wind_direction_used = std::min(wind_direction_raw, 7);
  const double wind_speed = draw_uniform();
  const double base = (4. - wind_speed) / 8.;
  const double wplus = base * (1. + wind_speed);
  const double wminus = base * (1. - wind_speed);
  double wind_prob[8] = {wplus,wplus,base,wminus,wminus,wminus,base,wplus};

  IntegerVector size_targets(n_regions);
  IntegerVector sites_in_event(n_regions);
  const int source = initiation_cell_index - 1;
  int eco = region_index[source] - 1;
  if (eco < 0 || eco >= n_regions) stop("247: Invalid ignition region");
  size_targets[eco] = draw_size(eco);
  ++g_event_count_249;
  int max_fire_region_size = size_targets[eco];

  std::vector<int> queue_local;
  std::vector<int> checked_local;
  std::vector<int>& queue = workspace265 ? workspace265->queue : queue_local;
  std::vector<int>& checked_order = workspace265 ?
      workspace265->checked_order : checked_local;
  std::vector<unsigned char> in_queue;
  LogicalVector event_member, burned, disturbed;
  NumericVector cell_severity;
  if (!workspace265) {
    queue.reserve(n);
    checked_order.reserve(n);
    in_queue.assign(n,0);
    event_member = LogicalVector(n,false);
    burned = LogicalVector(n,false);
    disturbed = LogicalVector(n,false);
    cell_severity = NumericVector(n,NA_REAL);
  }
  queue.push_back(source);
  if (workspace265) workspace265->queue_mark[source] = workspace265->generation;
  else in_queue[source] = 1;
  size_t qhead = 0;
  int total_sites_in_event = 0;
  int total_sites_damaged = 0;
  int cohorts_killed = 0;
  std::vector<int> killed_cohort_rows_264;
  double total_site_severities = 0.;

  // Original R uses sample.int(length(selected)) when >1, not a C++ shuffle.
  Function sample_int = Environment::base_env()["sample.int"];
  while (qhead < queue.size() &&
         sites_in_event[eco] < size_targets[eco] &&
         total_sites_in_event < max_fire_region_size) {
    const int site = queue[qhead++];
    if (workspace265) workspace265->queue_mark[site] = 0;
    else in_queue[site] = 0;
    checked_order.push_back(site+1);
    int site_region = region_index[site] - 1;
    if (site_region < 0 || site_region >= n_regions) stop("247: Invalid site region");
    if (site_region != eco) {
      eco = site_region;
      if (size_targets[eco] < 1) {
        size_targets[eco] = draw_size(eco);
        if (size_targets[eco] > max_fire_region_size)
          max_fire_region_size = size_targets[eco];
      }
    }
    event_occupied[site] = true;
    if (!workspace265) event_member[site] = true;

    const double time_since_fire = current_time - last_fire_year[site];
    int severity = 0;
    for (int k=0;k<5;++k) {
      const double threshold = fuel_curve(site_region,k);
      if (threshold != -1. && time_since_fire >= threshold) severity = k+1;
    }
    if (has_wind) {
      const double since_wind = current_time - last_wind_year[site];
      for (int k=0;k<5;++k) {
        const double threshold = wind_curve(site_region,k);
        if (threshold != -1. && since_wind <= threshold && k+1 > severity)
          severity = k+1;
      }
    }
    if (workspace265) workspace265->site_severity[site] = severity;
    else cell_severity[site] = severity;
    if (severity <= 0) continue;

    ++sites_in_event[eco];
    ++total_sites_in_event;
    total_site_severities += severity;
    if (!workspace265) {
      burned[site] = true;
      disturbed[site] = true;
    }
    last_fire_year[site] = current_time;

    // Deferred mortality mask from 243. The table is compacted by R only once
    // after the timestep, not for each burned site.
    int killed_here = 0;
    auto apply_mortality = [&](int row) {
      if (row < 0 || row >= alive.size()) stop("247: Invalid cohort row index");
      if (!alive[row]) return;
      bool kill = (severity == 5);
      if (!kill) {
        const int sp = cohort_species_index[row] - 1;
        if (sp < 0 || sp >= species_fire_tolerance.size())
          stop("247: Invalid cohort species index");
        const int difference = severity - species_fire_tolerance[sp];
        const float age_fraction = static_cast<float>(
          static_cast<float>(cohort_age[row]) /
          static_cast<float>(species_longevity[sp]));
        for (int k = 0; k < damage_max_age.size(); ++k) {
          if (difference >= damage_severity_tolerance_diff[k] &&
              static_cast<float>(damage_max_age[k]) >= age_fraction) {
            kill = true;
            break;
          }
        }
      }
      if (kill) {
        alive[row] = false;
        ++killed_here;
        if (sparse_result_264) killed_cohort_rows_264.push_back(row + 1);
      }
    };
    if (cohort_index) {
      const int start = cohort_index->offsets[site];
      const int end = cohort_index->offsets[site + 1];
      if (cohort_index->row_order.empty()) {
        for (int row = start; row < end; ++row) apply_mortality(row);
      } else {
        for (int i = start; i < end; ++i) apply_mortality(cohort_index->row_order[i]);
      }
    } else {
      SEXP site_rows = rows_by_site[site];
      if (!Rf_isNull(site_rows)) {
        for (int row : IntegerVector(site_rows)) apply_mortality(row - 1);
      }
    }
    cohorts_killed += killed_here;
    if (killed_here > 0) ++total_sites_damaged;

    // Base Fire 4.0 consumes a draw even for out-of-bounds neighbors.
    // Its spatial API returns inactive sites for those positions, and any
    // selected inactive site still participates in the shuffle.
    int chosen_full[9];
    int n_chosen = 0;
    for (int j=0;j<8;++j) {
      const int full = neighbors_full(site,j);
      if (full == 0 && !v7) continue;
      const double p = wind_prob[(j+wind_direction_used) % 8];
      const double draw = draw_uniform();
      if (draw < p) chosen_full[n_chosen++] = full;
    }
    // v7 uses a positive row offset for the far neighbor at wind index 0.
    const int far_direction = v7 ? (12 - wind_direction_used) % 8 : wind_direction_used;
    const int second = second_neighbors_full(site, far_direction);
    if (v7 || second != 0) {
      const double draw = draw_uniform();
      if (draw < wind_speed) chosen_full[n_chosen++] = second;
    }
    if (n_chosen == 0) continue;

    IntegerVector permutation;
    if (v7) {
      // LANDIS-II samples from the remaining list, including its last item.
      for (int k = 0; k < n_chosen; ++k) {
        int chosen = k + v7->stream.next(n_chosen - k);
        const int full = chosen_full[chosen];
        for (; chosen > k; --chosen) chosen_full[chosen] = chosen_full[chosen - 1];
        chosen_full[k] = full;
      }
    } else if (n_chosen > 1) {
      synced_rng.publish();
      permutation = sample_int(n_chosen);
      ++g_event_shuffle_callbacks_249;
      synced_rng.reload();
    }
    for (int k=0;k<n_chosen;++k) {
      const int chosen = v7 ? k : (n_chosen > 1 ? permutation[k] - 1 : 0);
      const int full = chosen_full[chosen] - 1;
      if (v7 && full == -1) continue;
      if (full < 0 || full >= full_active.size()) stop("247: Invalid full neighbor");
      if (full_active[full] != TRUE) continue;
      const int neighbor = full_to_active[full] - 1;
      if (neighbor < 0 || neighbor >= n) stop("247: Invalid active neighbor");
      if (event_occupied[neighbor] == TRUE ||
          (workspace265 ? workspace265->queue_mark[neighbor] == workspace265->generation
                        : in_queue[neighbor])) continue;
      const double tsf = current_time - last_fire_year[neighbor];
      const double burn_probability = 1. - std::exp(tsf * (-1. / fire_spread_age[neighbor]));
      const double burn_draw = draw_uniform();
      if (burn_draw <= burn_probability) {
        queue.push_back(neighbor);
        if (workspace265) workspace265->queue_mark[neighbor] = workspace265->generation;
        else in_queue[neighbor] = 1;
        if (queue.size() > static_cast<size_t>(n))
          stop("247: Event queue overflow");
      }
    }
  }

  // Guarantee the R environment retains the cohort mask for later fires.
  cohort_env.assign("alive", alive);
  if (workspace265) workspace265->record();
  std::vector<int> member_idx;
  member_idx.reserve(checked_order.size());
  // 262: checked_order enumerates every event member exactly once.
  // Reproduce which(event_member) order without an O(n) full-landscape scan.
  for (int site_one : checked_order) member_idx.push_back(site_one - 1);
  std::sort(member_idx.begin(), member_idx.end());
  if (std::adjacent_find(member_idx.begin(), member_idx.end()) != member_idx.end())
    stop("262: duplicate event member; sparse enumeration rejected");
  const int nmember = static_cast<int>(member_idx.size());
  IntegerVector member_cell_index(nmember), member_cell_id(nmember),
    member_region_index(nmember);
  NumericVector member_severity(nmember);
  LogicalVector member_burned(nmember),member_disturbed(nmember);
  for (int i=0;i<nmember;++i) {
    const int k=member_idx[i];
    member_cell_index[i]=k+1;
    member_cell_id[i]=cell_id[k];
    member_region_index[i]=region_index[k];
    member_severity[i] = workspace265 ? workspace265->site_severity[k] : cell_severity[k];
    const bool is_burned = workspace265 ? member_severity[i] > 0. : burned[k] == TRUE;
    member_burned[i] = is_burned;
    member_disturbed[i] = workspace265 ? is_burned : disturbed[k];
  }
  DataFrame cells = DataFrame::create(
    _["cell_index"] = member_cell_index,
    _["cell_id"] = member_cell_id,
    _["region_index"] = member_region_index,
    _["severity"] = member_severity,
    _["burned"] = member_burned,
    _["disturbed"] = member_disturbed);
  const double event_severity = total_sites_damaged == 0 ? 0. :
    total_site_severities / total_sites_in_event;
  IntegerVector checked(checked_order.begin(), checked_order.end());
  List result = List::create(
    _["initiation_cell_index"] = initiation_cell_index,
    _["initiation_cell_id"] = cell_id[source],
    _["wind_direction_raw"] = wind_direction_raw,
    _["wind_direction"] = wind_direction_used,
    _["wind_speed"] = wind_speed,
    _["size_targets"] = size_targets,
    _["sites_in_event"] = sites_in_event,
    _["max_fire_region_size"] = max_fire_region_size,
    _["num_site_checked"] = static_cast<int>(checked_order.size()),
    _["checked_order"] = checked,
    _["total_sites_in_event"] = total_sites_in_event,
    _["total_sites_damaged"] = total_sites_damaged,
    _["cohorts_killed"] = cohorts_killed,
    _["severity"] = event_severity,
    _["cells"] = cells);
  if (sparse_result_264) {
    // Only touched cohorts, indexed into the pre-compaction cohort table.
    result["killed_cohort_rows"] = IntegerVector(
      killed_cohort_rows_264.begin(), killed_cohort_rows_264.end());
  } else {
    // Historical public event return contract, not used for timestep 264.
    result["event_member"] = event_member;
    result["event_occupied"] = event_occupied;
    result["last_fire_year"] = last_fire_year;
  }
  result.attr("class") = CharacterVector::create("landisR_fire_basefire_spread_event");
  return result;
}


// 260: A single Rcpp dispatch for the complete ignition/event loop of one
// Base Fire timestep. Existing 249 kernel remains the numerical ground truth.
// Internal R callbacks (fire-size and sample.int) preserve RNG exactly via
// RngLease249::publish()/reload(). This function does not reconstruct CSR or
// alter the biological rules. Experimental; no changes to installed packages.
// [[Rcpp::export(rng = false)]]
List basefire_timestep_kernel_260(
    int current_time,
    double timestep,
    const NumericVector& ignition_probability,
    const NumericVector& ignition_fire_spread_age,
    const IntegerVector& cell_id,
    const IntegerVector& region_index,
    int n_regions,
    const IntegerMatrix& neighbors_full,
    const IntegerMatrix& second_neighbors_full,
    const IntegerVector& full_to_active,
    const LogicalVector& full_active,
    const NumericVector& fire_spread_age,
    const NumericMatrix& fuel_curve,
    const NumericMatrix& wind_curve,
    bool has_wind,
    const NumericVector& last_wind_year,
    const NumericVector& last_fire_year_input,
    const LogicalVector& event_occupied_input,
    const NumericVector& last_severity_input,
    SEXP rows_by_site,
    Environment cohort_env,
    const NumericVector& cohort_age,
    const IntegerVector& cohort_species_index,
    const IntegerVector& species_fire_tolerance,
    const NumericVector& species_longevity,
    const NumericVector& damage_max_age,
    const IntegerVector& damage_severity_tolerance_diff,
    Function size_one, SEXP v7_config = R_NilValue) {
  using Clock = std::chrono::steady_clock;
  using Seconds = std::chrono::duration<double>;
  const bool timing = !Rf_isNull(v7_config) &&
    List(v7_config).containsElementNamed("timing") && as<bool>(List(v7_config)["timing"]);
  const auto start = timing ? Clock::now() : Clock::time_point{};
  const int n = cell_id.size();
  if (n < 1 || ignition_probability.size() != n ||
      ignition_fire_spread_age.size() != n ||
      last_fire_year_input.size() != n ||
      event_occupied_input.size() != n ||
      last_severity_input.size() != n ||
      region_index.size() != n || n_regions < 1 ||
      !R_finite(timestep) || timestep <= 0 ||
      fuel_curve.nrow() < n_regions || wind_curve.nrow() < n_regions ||
      fuel_curve.ncol() != 5 || wind_curve.ncol() != 5 ||
      (has_wind && last_wind_year.size() != n)) {
    stop("260: invalid input dimensions");
  }
  landis_progress::Reporter progress;
  progress.phase(0, "Preparing fire workspace");
  // 265: O(N) scratch initialized only ONCE per timestep, not per fire.
  Rcpp::XPtr<FireWorkspace265> ws265(new FireWorkspace265(n), true);
  List legacy_rows;
  const auto index_start = timing ? Clock::now() : Clock::time_point{};
  if (TYPEOF(rows_by_site) == INTSXP) {
    const IntegerVector cohort_sites(rows_by_site);
    if (cohort_sites.size() != cohort_age.size())
      stop("Invalid Base Fire cohort cell position count");
    ws265->cohort_index.build(n, cohort_sites);
  } else {
    legacy_rows = List(rows_by_site);
  }
  const auto index_end = timing ? Clock::now() : Clock::time_point{};
  std::unique_ptr<FireV7Random> v7;
  IntegerVector scan_order;
  if (!Rf_isNull(v7_config)) {
    const List config(v7_config);
    scan_order = config["scan_order"];
    if (scan_order.size() != n) stop("Invalid LANDIS-II v7 fire scan order");
    v7.reset(new FireV7Random(config));
    if (v7->mean_size.size() != n_regions || v7->min_size.size() != n_regions ||
        v7->max_size.size() != n_regions || !R_finite(v7->cell_area) || v7->cell_area <= 0)
      stop("Invalid LANDIS-II v7 fire size cache");
    ws265->v7 = v7.get();
  }
  NumericVector last_fire = clone(last_fire_year_input);
  LogicalVector occupied = clone(event_occupied_input);
  NumericVector severity(n);
  LogicalVector disturbed(n, false);
  IntegerVector event_owner(n);
  NumericVector last_severity = clone(last_severity_input);
  // 262: event output capacity scales with F instead of N. The List itself
  // remains R-protected across RNG callbacks, and individual event objects
  // remain identical to 260 (including full historical per-event vectors).
  List raw_events(std::min(n, 16));
  int scanned = 0;
  int events_run = 0;
  int ignition_draws = 0;
  // Keep a lease only across the ignition scan. The event kernel opens its
  // own lease and publishes before calling R's rexp/sample.int callbacks.
  // Never nest leases: publish this one before the event, then reacquire it.
  progress.start(0.08, 1, "Ignition, spread & severity", n);
  RngLease249 lease(!v7);
  auto draw_uniform = [&]() { return v7 ? v7->stream.next_double() : lease.uniform(); };
  const auto scan_start = timing ? Clock::now() : Clock::time_point{};
  double event_seconds = 0;
  for (int i = 0; i < n; ++i) {
    if (progress.due(i)) {
      if (!v7) lease.publish();
      progress.advance(i);
      if (!v7) lease.reload();
    }
    const int pos = v7 ? scan_order[i] - 1 : i;
    if (pos < 0 || pos >= n) stop("Invalid LANDIS-II v7 fire scan position");
    ++scanned;
    const double u0 = draw_uniform();
    ++ignition_draws;
    if (!(u0 <= ignition_probability[pos] * timestep)) continue;
    const double dt = current_time - last_fire[pos];
    const double p_fire = 1. - std::exp(dt * (-1. / ignition_fire_spread_age[pos]));
    const double u1 = draw_uniform();
    ++ignition_draws;
    if (!(u1 <= p_fire)) continue;
    int reg = region_index[pos] - 1;
    if (reg < 0 || reg >= n_regions) stop("260: invalid ignition region");
    const double time_since_fire = current_time - last_fire[pos];
    int s = 0;
    for (int k=0;k<5;++k) {
      const double threshold = fuel_curve(reg,k);
      if (threshold != -1. && time_since_fire >= threshold) s = k+1;
    }
    if (has_wind) {
      const double since_wind = current_time - last_wind_year[pos];
      for (int k=0;k<5;++k) {
        const double threshold = wind_curve(reg,k);
        if (threshold != -1. && since_wind <= threshold && k+1 > s)
          s = k+1;
      }
    }
    if (s <= 0) continue;

    const auto event_start = timing ? Clock::now() : Clock::time_point{};
    if (!v7) lease.publish();
    List ev = basefire_event_kernel_247(
      pos+1, current_time, cell_id, region_index, n_regions,
      neighbors_full, second_neighbors_full, full_to_active, full_active,
      fire_spread_age, fuel_curve, wind_curve, has_wind,
      last_wind_year, last_fire, occupied,
      last_severity_input, legacy_rows, cohort_env, cohort_age, cohort_species_index,
      species_fire_tolerance, species_longevity,
      damage_max_age, damage_severity_tolerance_diff, size_one,
      true, ws265); // 265: same sparse results, reusable scratch workspace.
    if (!v7) lease.reload();
    ++events_run;
    if (events_run > raw_events.size()) {
      const int cap = raw_events.size();
      const int next_cap = std::min(n, std::max(events_run, cap * 2));
      List grown(next_cap);
      for (int k=0;k<cap;++k) grown[k] = raw_events[k];
      raw_events = grown;
    }
    raw_events[events_run - 1] = ev;
    // Mirrors the R event loop: after each event update running state and
    // owner/severity arrays, BEFORE examining the next ignition candidate.
    // 264: last_fire and occupied were mutated in place.  Neither field is
    // copied into the event, or read back from it.
    // 262: event_cells already contains precisely the sorted sites marked
    // event_member in the 247 kernel. Eliminate the O(N) owner scan. Update
    // owner for ALL visited cells, including those with severity 0.
    DataFrame event_cells = ev["cells"];
    IntegerVector idx = event_cells["cell_index"];
    NumericVector sev = event_cells["severity"];
    LogicalVector burned = event_cells["burned"];
    for (int j=0;j<idx.size();++j) {
      const int ai = idx[j] - 1;
      if (ai < 0 || ai >= n) stop("262: invalid event cell index");
      event_owner[ai] = events_run;
      if (burned[j] != TRUE) continue;
      severity[ai] = sev[j];
      disturbed[ai] = true;
      last_severity[ai] = sev[j];
    }
    if (timing) event_seconds += Seconds(Clock::now() - event_start).count();
  }
  const auto scan_end = timing ? Clock::now() : Clock::time_point{};
  // The raw event objects are returned in identical order and format so
  // existing R event logs and metadata can be built unchanged, outside C++.
  List events(events_run);
  for (int i=0;i<events_run;++i) events[i] = raw_events[i];
  List result = List::create(
      _["raw_events"] = events,
      _["last_fire_year"] = last_fire,
      _["event_occupied"] = occupied,
      _["last_severity"] = last_severity,
      _["severity"] = severity,
      _["disturbed"] = disturbed,
      _["event_owner"] = event_owner,
      _["number_fires"] = events_run,
      _["scanned"] = scanned,
      _["ignition_rng_draws"] = ignition_draws,
      _["workspace_events_265"] = ws265->events,
      _["workspace_n_265"] = ws265->n,
      _["workspace_max_queue_265"] = static_cast<int>(ws265->max_queue),
      _["workspace_max_checked_265"] = static_cast<int>(ws265->max_checked),
      _["workspace_queue_capacity_265"] = static_cast<int>(ws265->queue.capacity()),
      _["workspace_checked_capacity_265"] = static_cast<int>(ws265->checked_order.capacity()));
  if (v7) {
    result["landis_rng"] = List::create(_["state"] = v7->stream.encode(),
        _["offset"] = static_cast<double>(v7->stream.draws));
  }
  if (timing) {
    result["timing"] = NumericVector::create(
        _["initialisation"] = Seconds(scan_start - start - (index_end - index_start)).count(),
        _["index_cohortes_natif"] = Seconds(index_end - index_start).count(),
        _["balayage_et_ignitions"] = Seconds(scan_end - scan_start).count() - event_seconds,
        _["propagation_et_mortalite"] = event_seconds,
        _["retour"] = Seconds(Clock::now() - scan_end).count());
  }
  return result;
}
