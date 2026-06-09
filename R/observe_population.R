# Observation processes kept separate from the biological generator.

sample_size_profiles <- function(
  census,
  detection = function(size, year) rep(0.6, length(size)),
  effort = NULL,
  measurement_sd = 0
) {
  years <- sort(unique(census$year))
  if (is.null(effort)) {
    effort <- rep(1, length(years))
  }
  stopifnot(length(effort) == length(years))

  observed <- lapply(seq_along(years), function(index) {
    year <- years[index]
    available <- census[census$year == year, ]
    probability <- pmin(
      1,
      pmax(0, effort[index] * detection(available$size, year))
    )
    selected <- rbinom(nrow(available), 1L, probability) == 1L
    result <- available[selected, c("year", "id", "size")]
    result$observed_size <- result$size + rnorm(nrow(result), sd = measurement_sd)
    result$detection_probability <- probability[selected]
    result$effort <- effort[index]
    result
  })

  do.call(rbind, observed)
}

sample_mark_recapture <- function(
  census,
  detection = function(size, year) rep(0.35, length(size)),
  measurement_sd = 0
) {
  probability <- pmin(
    1,
    pmax(0, detection(census$size, census$year))
  )
  captured <- rbinom(nrow(census), 1L, probability) == 1L
  result <- census[captured, c("year", "id", "age", "sex", "size")]
  result$observed_size <- result$size + rnorm(nrow(result), sd = measurement_sd)
  result$capture_probability <- probability[captured]
  result
}

bin_size_profiles <- function(profiles, breaks) {
  stopifnot(all(diff(breaks) > 0))
  years <- sort(unique(profiles$year))

  binned <- lapply(years, function(year) {
    sizes <- profiles$observed_size[profiles$year == year]
    counts <- hist(sizes, breaks = breaks, plot = FALSE)$counts
    data.frame(
      year = year,
      lower = head(breaks, -1L),
      upper = tail(breaks, -1L),
      midpoint = (head(breaks, -1L) + tail(breaks, -1L)) / 2,
      count = counts
    )
  })

  do.call(rbind, binned)
}
