# functions/data_helpers.R
#
# Loading, cleaning, event matching and curve factors for the BNR fixing.
# Maturities and curve weights are defined here, and every script takes them
# from this file.


# The Wald restrictions in 05 are built from this sequence, so the order is
# load-bearing.
MATURITIES <- c(dy_12m = "12M", dy_3y = "3Y", dy_5y = "5Y", dy_10y = "10Y")

# Fixed contrasts of the four quoted maturities, imposed by construction. 04
# compares them with an unrestricted principal component decomposition.
CURVE_WEIGHTS <- rbind(
  dy_level     = c( 0.25, 0.25, 0.25,  0.25),
  dy_slope     = c(-1.00, 0.00, 0.00,  1.00),
  dy_curvature = c(-1.00, 0.00, 2.00, -1.00)
)
colnames(CURVE_WEIGHTS) <- names(MATURITIES)


# BNR fixing

# The workbook opens with six pre-data rows - an accessibility notice, the
# title, a "Nota:" line, column headers, units and the BNR series codes - hence
# skip = 6 and column names supplied by hand.
load_bnr_raw <- function(path = "data/raw/titluri_de_stat_ro.xlsx") {
  raw <- readxl::read_excel(path, sheet = "Sheet1", col_names = FALSE, skip = 6)
  names(raw) <- c(
    "date",
    "bid_6m", "bid_12m", "bid_3y", "bid_5y", "bid_10y",
    "ask_6m", "ask_12m", "ask_3y", "ask_5y", "ask_10y"
  )
  raw
}

# In yield space the bid quote sits above the ask, since a bid price below the
# ask implies a bid yield above it. Mid is the simple average and the spread is
# bid - ask > 0. The output is sorted. Matching also requires one row per date,
# and 01 asserts that uniqueness immediately after cleaning, since the panel
# contract examines ECB event dates while the matching indexes BNR fixing
# dates.
clean_bnr <- function(raw) {
  raw |>
    dplyr::mutate(date = as.Date(date)) |>
    dplyr::mutate(dplyr::across(-date, as.numeric)) |>
    dplyr::filter(!dplyr::if_all(-date, is.na)) |>
    dplyr::arrange(date) |>
    dplyr::mutate(
      mid_6m  = (bid_6m  + ask_6m)  / 2,
      mid_12m = (bid_12m + ask_12m) / 2,
      mid_3y  = (bid_3y  + ask_3y)  / 2,
      mid_5y  = (bid_5y  + ask_5y)  / 2,
      mid_10y = (bid_10y + ask_10y) / 2,
      spread_10y_bps = 100 * (bid_10y - ask_10y)
    )
}


# Event matching

# Index of the first BNR fixing struck strictly after each ECB event.
#
# The fixing is struck at 12:00 Bucharest, i.e. 11:00 Frankfurt year-round.
# ECB announcements land at 13:45 Frankfurt through 20 July 2022 and at 14:15
# afterwards. The announcement-day fixing is therefore predetermined with
# respect to the surprise, and the response has to be read off the next one.
#
# findInterval counts the fixings at or before each event date, so +1 gives the
# first one struck strictly after. Indexing the fixings that actually exist
# carries weekends and Romanian public holidays correctly and keeps the one
# event that falls on a closed Romanian day, all of which a calendar lead would
# get wrong. Positions may exceed length(fixing_dates), and the caller trims.
first_available_fixing <- function(event_dates, fixing_dates) {
  findInterval(event_dates, fixing_dates) + 1L
}

# Response and placebo windows, all in basis points.
#
#   dy    = y[F(i)]     - y[F(i) - 1]   the response window
#   pre   = y[F(i) - 1] - y[F(i) - 2]   predetermined w.r.t. the surprise
#   post  = y[F(i) + 1] - y[F(i)]       the window following the response
#
# The placebos are built here, beside the response and by the same procedure,
# which supports the falsification exercise in 02. An event enters once
# it has two preceding and one following fixing, a condition every event in
# this sample meets.
build_event_panel <- function(bnr, events) {
  fix_dates <- bnr$date
  n_fix     <- length(fix_dates)
  short     <- sub("^dy_", "", names(MATURITIES))

  f      <- first_available_fixing(events$date, fix_dates)
  usable <- (f - 2L >= 1L) & (f + 1L <= n_fix)

  ev <- events[usable, , drop = FALSE]
  f  <- f[usable]

  panel <- ev |>
    dplyr::mutate(
      fixing_date     = fix_dates[f],
      days_to_fixing  = as.integer(fixing_date - date),
      on_fixing_day   = date %in% fix_dates
    )

  for (m in short) {
    col <- paste0("mid_", m)
    y   <- bnr[[col]]
    panel[[paste0("dy_",   m)]] <- 100 * (y[f]      - y[f - 1L])
    panel[[paste0("pre_",  m)]] <- 100 * (y[f - 1L] - y[f - 2L])
    panel[[paste0("post_", m)]] <- 100 * (y[f + 1L] - y[f])
  }

  panel
}


# CURVE_WEIGHTS applied to the event-aligned changes. Level, Slope and
# Curvature summarise the same four yields, and 01 asserts that the stored
# columns reproduce the contrasts exactly.
add_curve_factors <- function(panel) {
  X       <- as.matrix(panel[, colnames(CURVE_WEIGHTS)])
  factors <- X %*% t(CURVE_WEIGHTS)
  for (nm in rownames(CURVE_WEIGHTS)) panel[[nm]] <- as.vector(factors[, nm])
  panel
}
