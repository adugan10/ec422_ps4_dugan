library(tidyverse)
library(fpp3)
library(BVAR)

dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
theme_set(theme_minimal(base_size = 12))

dm_row <- function(test, comparison) {
  tibble(Comparison = comparison, DM = unname(test$statistic),
         p_value = test$p.value)
}

liquor_raw <- read_csv(
  "data/MRTSSM4453USN.csv", na = ".",
  col_types = cols(observation_date = col_date(format = "%Y-%m-%d"),
                   MRTSSM4453USN = col_double())
)

liquor <- liquor_raw |>
  transmute(Month = yearmonth(observation_date),
            sales = MRTSSM4453USN) |>
  as_tsibble(index = Month) |>
  filter(Month >= yearmonth("1992 Jan"),
         Month <= yearmonth("2019 Dec"))

print(has_gaps(liquor))

stopifnot(
  nrow(liquor) == 336L,
  all(is.finite(liquor$sales)),
  all(liquor$sales > 0),
  !any(has_gaps(liquor)$.gaps)
)

sales_plot <- autoplot(liquor, sales) +
  labs(
    title = "U.S. Beer, Wine, and Liquor Store Sales",
    subtitle = "January 1992–December 2019",
    x = "Month",
    y = "Sales (millions of dollars)"
  )

season_plot <- gg_season(
  liquor, log(sales),
  pal = gray.colors(9, start = 0.68, end = 0.12)
) +
  labs(
    title = "Seasonality in Log Liquor Store Sales",
    x = "Month",
    y = "Log sales",
    colour = "Year"
  )

print(sales_plot)
print(season_plot)

ggsave("results/liquor_sales.png", sales_plot,
       width = 9, height = 4.8, dpi = 160)

ggsave("results/liquor_season.png", season_plot,
       width = 9, height = 5.2, dpi = 160)

liquor_fit <- liquor |>
  model(
    arima = ARIMA(log(sales)),
    ets = ETS(log(sales))
  )

liquor_fit |> select(arima) |> report()
liquor_fit |> select(ets) |> report()

saveRDS(liquor_fit, "results/liquor_fit.rds")

arima_dof <- liquor_fit |>
  select(arima) |>
  tidy() |>
  filter(str_detect(term, "^(ar|ma|sar|sma)[0-9]+$")) |>
  nrow()

residual_plot <- liquor_fit |>
  select(arima) |>
  gg_tsresiduals()

print(residual_plot)

ggsave("results/liquor_residuals.png", residual_plot,
       width = 9, height = 6.4, dpi = 160)

liquor_lb <- liquor_fit |>
  select(arima) |>
  augment() |>
  features(.innov, ljung_box, lag = 24, dof = arima_dof) |>
  mutate(
    lb_pvalue = pchisq(
      lb_stat,
      df = 24 - arima_dof,
      lower.tail = FALSE
    )
  )

print(liquor_lb)

saveRDS(
  list(dof = arima_dof, ljung_box = liquor_lb),
  "results/liquor_diagnostics.rds"
)

liquor_cv <- liquor |>
  stretch_tsibble(.init = 120, .step = 1) |>
  filter(.id < max(.id))

cat("Forecast origins:", n_distinct(liquor_cv$.id), "\n")

liquor_cv_fit <- liquor_cv |>
  model(
    snaive = SNAIVE(sales),
    ets = ETS(log(sales)),
    arima = ARIMA(log(sales))
  )

saveRDS(
  liquor_cv_fit,
  "results/liquor_cv_fit.rds",
  compress = "xz"
)

liquor_cv_fc <- liquor_cv_fit |>
  forecast(h = 1)

saveRDS(liquor_cv_fc, "results/liquor_cv_fc.rds")

liquor_accuracy <- liquor_cv_fc |>
  accuracy(liquor, measures = list(RMSE = RMSE, MAE = MAE))

print(range(liquor_cv_fc$Month))
print(liquor_accuracy)

saveRDS(liquor_accuracy, "results/liquor_accuracy.rds")

liquor_errs <- liquor_cv_fc |>
  as_tibble() |>
  left_join(
    as_tibble(liquor) |> select(Month, actual = sales),
    by = "Month"
  ) |>
  mutate(e = actual - .mean) |>
  select(Month, .model, e) |>
  pivot_wider(names_from = .model, values_from = e) |>
  arrange(Month)

stopifnot(nrow(liquor_errs) == 216L, !anyNA(liquor_errs))

dm_arima_ets <- forecast::dm.test(
  liquor_errs$arima, liquor_errs$ets,
  h = 1, power = 2
)

dm_arima_snaive <- forecast::dm.test(
  liquor_errs$arima, liquor_errs$snaive,
  h = 1, power = 2
)

liquor_dm <- bind_rows(
  dm_row(dm_arima_ets, "ARIMA versus ETS"),
  dm_row(dm_arima_snaive, "ARIMA versus seasonal naive")
)

print(liquor_dm)

saveRDS(liquor_errs, "results/liquor_errs.rds")
saveRDS(liquor_dm, "results/liquor_dm.rds")

