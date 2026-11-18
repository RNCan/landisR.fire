#include <Rcpp.h>
#include <algorithm>
#include <climits>
#include <cstdint>
#include <vector>

// Full-landscape indices match the strict Base Fire site table. Inactive
// sites retain an index and participate in the v7 neighbor selection order.
// [[Rcpp::export(rng = false)]]
Rcpp::List basefire_site_neighbors(
    const Rcpp::IntegerVector& rows,
    const Rcpp::IntegerVector& cols,
    const Rcpp::IntegerVector& active_to_full,
    const Rcpp::IntegerVector& drow,
    const Rcpp::IntegerVector& dcol,
    const Rcpp::CharacterVector& directions) {
  if (rows.size() == 0 || rows.size() != cols.size() || rows.size() > INT_MAX ||
      active_to_full.size() == 0 || active_to_full.size() > rows.size() ||
      drow.size() != 8 || dcol.size() != 8 || directions.size() != 8) {
    Rcpp::stop("Invalid Base Fire site geometry dimensions.");
  }

  int row_min = INT_MAX, col_min = INT_MAX;
  int row_max = INT_MIN, col_max = INT_MIN;
  const int n_full = static_cast<int>(rows.size());
  const int n_active = static_cast<int>(active_to_full.size());
  for (int i = 0; i < n_full; ++i) {
    const int row = rows[i], col = cols[i];
    if (row == NA_INTEGER || col == NA_INTEGER)
      Rcpp::stop("Base Fire landscape coordinates cannot contain NA.");
    row_min = std::min(row_min, row);
    row_max = std::max(row_max, row);
    col_min = std::min(col_min, col);
    col_max = std::max(col_max, col);
  }
  const std::int64_t height = static_cast<std::int64_t>(row_max) - row_min + 1;
  const std::int64_t width = static_cast<std::int64_t>(col_max) - col_min + 1;
  std::vector<int> lookup;
  if (static_cast<std::uint64_t>(width) > lookup.max_size() / height)
    Rcpp::stop("Base Fire landscape coordinates span too large a grid.");
  lookup.assign(static_cast<size_t>(height) * static_cast<size_t>(width), 0);
  auto raster_index = [height](std::int64_t row, std::int64_t col) {
    return static_cast<size_t>(row + col * height);
  };
  for (int i = 0; i < n_full; ++i) {
    const std::int64_t row = static_cast<std::int64_t>(rows[i]) - row_min;
    const std::int64_t col = static_cast<std::int64_t>(cols[i]) - col_min;
    int& slot = lookup[raster_index(row, col)];
    if (slot != 0)
      Rcpp::stop("Full Base Fire landscape row/col coordinates must be unique.");
    slot = i + 1;
  }

  Rcpp::IntegerMatrix first(n_active, 8), second(n_active, 8);
  for (int j = 0; j < 8; ++j) {
    for (int i = 0; i < n_active; ++i) {
      const int full = active_to_full[i];
      if (full < 1 || full > n_full)
        Rcpp::stop("Invalid active cell index in Base Fire site geometry.");
      const std::int64_t row = static_cast<std::int64_t>(rows[full - 1]) - row_min;
      const std::int64_t col = static_cast<std::int64_t>(cols[full - 1]) - col_min;
      const std::int64_t r1 = row + drow[j], c1 = col + dcol[j];
      const std::int64_t r2 = r1 + drow[j], c2 = c1 + dcol[j];
      if (r1 >= 0 && r1 < height && c1 >= 0 && c1 < width)
        first(i, j) = lookup[raster_index(r1, c1)];
      if (r2 >= 0 && r2 < height && c2 >= 0 && c2 < width)
        second(i, j) = lookup[raster_index(r2, c2)];
    }
  }
  Rcpp::colnames(first) = directions;
  Rcpp::colnames(second) = directions;
  return Rcpp::List::create(first, second);
}
