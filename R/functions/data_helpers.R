# ==============================================================================
# functions/data_helpers.R
#
# BNR fixing: loading, cleaning, event matching, curve factors. Maturities and
# curve weights are defined here and nowhere else in the project.
# ==============================================================================


# Order matters: the Wald restrictions in 05 are built from this sequence.
MATURITIES <- c(dy_12m = "12M", dy_3y = "3Y", dy_5y = "5Y", dy_10y = "10Y")

# Fixed contrasts of the four quoted maturities, imposed rather than estimated.
# 04 checks them against an unrestricted principal component decomposition.
CURVE_WEIGHTS <- rbind(
  dy_level     = c( 0.25, 0.25, 0.25,  0.25),
  dy_slope     = c(-1.00, 0.00, 0.00,  1.00),
  dy_curvature = c(-1.00, 0.00, 2.00, -1.00)
)
colnames(CURVE_WEIGHTS) <- names(MATURITIES)


# ------------------------------------------------------------------------------
# BNR fixing
# ------------------------------------------------------------------------------

# Six pre-data rows in the workbook: accessibility notice, title, "Nota:",
# column headers, units, BNR series codes. Hence skip = 6 and names by hand.
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
# ask implies a bid yield above it. Mid is the simple average, spread is
# bid - ask > 0. Sorting is not cosmetic - findInterval() below assumes it, and
# so is uniqueness, which 01 asserts on this output: two rows sharing a fixing
# date would let findInterval() pick either one, and the panel contract cannot
# catch that because it checks the ECB event dates, not the BNR fixing dates.
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


# ------------------------------------------------------------------------------
# Event matching
# ------------------------------------------------------------------------------

# Index of the first BNR fixing struck strictly after each ECB event.
#
# The fixing is struck at 12:00 Bucharest, i.e. 11:00 Frankfurt year-round.
# ECB announcements land at 13:45 Frankfurt through 20 July 2022 and at 14:15
# afterwards. The announcement-day fixing is therefore predetermined with
# respect to the surprise, and the response has to be read off the next one.
#
# findInterval counts the fixings at or before each event date, so +1 gives the
# first strictly after. Not lead(date, 1): that mishandles weekends and
# Romanian public holidays, and drops the one event falling on a closed day.
# Positions may exceed length(fixing_dates); the caller trims.
first_available_fixing <- function(event_dates, fixing_dates) {
  findInterval(event_dates, fixing_dates) + 1L
}

# Response and placebo windows, all in basis points:
#
#   dy    = y[F(i)]     - y[F(i) - 1]   response
#   pre   = y[F(i) - 1] - y[F(i) - 2]   predetermined, should not react
#   post  = y[F(i) + 1] - y[F(i)]       should not persist
#
# The placebos are built here, next to the response and by the same procedure,
# which is what lets 02 read as a falsification design. Events without two
# preceding and one following fixing are dropped; none are on this sample.
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


# CURVE_WEIGHTS applied to the event-aligned changes. Summaries of the same
# four yields, not separately estimated factors; 01 asserts that the stored
# columns match the contrasts exactly.
add_curve_factors <- function(panel) {
  X       <- as.matrix(panel[, colnames(CURVE_WEIGHTS)])
  factors <- X %*% t(CURVE_WEIGHTS)
  for (nm in rownames(CURVE_WEIGHTS)) panel[[nm]] <- as.vector(factors[, nm])
  panel
}
