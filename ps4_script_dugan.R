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
  transmute(Month = yearmonth(observation_date), sales = MRTSSM4453USN) |>
  as_tsibble(index = Month) |>
  filter(Month >= yearmonth("1992 Jan"), Month <= yearmonth("2019 Dec"))
print(has_gaps(liquor))
stopifnot(nrow(liquor) == 336L, all(is.finite(liquor$sales)),
          all(liquor$sales > 0), !any(has_gaps(liquor)$.gaps))

sales_plot <- autoplot(liquor, sales) +
  labs(title = "U.S. Beer, Wine, and Liquor Store Sales",
       subtitle = "Monthly sales before seasonal adjustment, January 1992–December 2019",
       x = "Month", y = "Sales (millions of dollars)")
season_plot <- gg_season(liquor, log(sales), pal = gray.colors(9, start = 0.68, end = 0.12)) +
  labs(title = "Seasonality in Log Liquor Store Sales",
       x = "Month", y = "Log sales", colour = "Year")
print(sales_plot)
print(season_plot)
ggsave("results/liquor_sales.png", sales_plot, width = 9, height = 4.8, dpi = 160)
ggsave("results/liquor_season.png", season_plot, width = 9, height = 5.2, dpi = 160)

liquor_fit <- liquor |>
  model(arima = ARIMA(log(sales)), ets = ETS(log(sales)))
liquor_fit |> select(arima) |> report()
liquor_fit |> select(ets) |> report()
saveRDS(liquor_fit, "results/liquor_fit.rds")

arima_dof <- liquor_fit |> select(arima) |> tidy() |>
  filter(str_detect(term, "^(ar|ma|sar|sma)[0-9]+$")) |> nrow()
residual_plot <- liquor_fit |> select(arima) |> gg_tsresiduals()
print(residual_plot)
ggsave("results/liquor_residuals.png", residual_plot,
       width = 9, height = 6.4, dpi = 160)
liquor_lb <- liquor_fit |> select(arima) |> augment() |>
  features(.innov, ljung_box, lag = 24, dof = arima_dof) |>
  mutate(lb_pvalue = pchisq(lb_stat, df = 24 - arima_dof, lower.tail = FALSE))
print(liquor_lb)
saveRDS(list(dof = arima_dof, ljung_box = liquor_lb),
        "results/liquor_diagnostics.rds")

liquor_cv <- liquor |>
  stretch_tsibble(.init = 120, .step = 1) |>
  filter(.id < max(.id))
cat("Liquor forecast origins:", n_distinct(liquor_cv$.id), "\n")
message("Starting liquor cross-validation; this can take several minutes.")
liquor_cv_fit <- liquor_cv |>
  model(snaive = SNAIVE(sales), ets = ETS(log(sales)),
        arima = ARIMA(log(sales)))
saveRDS(liquor_cv_fit, "results/liquor_cv_fit.rds", compress = "xz")

liquor_cv_fc <- liquor_cv_fit |> forecast(h = 1)
saveRDS(liquor_cv_fc, "results/liquor_cv_fc.rds")
liquor_accuracy <- liquor_cv_fc |>
  accuracy(liquor, measures = list(RMSE = RMSE, MAE = MAE))
print(range(liquor_cv_fc$Month))
print(liquor_accuracy)
saveRDS(liquor_accuracy, "results/liquor_accuracy.rds")

liquor_errs <- liquor_cv_fc |> as_tibble() |>
  left_join(as_tibble(liquor) |> select(Month, actual = sales), by = "Month") |>
  mutate(e = actual - .mean) |> select(Month, .model, e) |>
  pivot_wider(names_from = .model, values_from = e) |> arrange(Month)
stopifnot(nrow(liquor_errs) == 216L, !anyNA(liquor_errs))
dm_arima_ets <- forecast::dm.test(liquor_errs$arima, liquor_errs$ets,
                                 h = 1, power = 2)
dm_arima_snaive <- forecast::dm.test(liquor_errs$arima, liquor_errs$snaive,
                                    h = 1, power = 2)
liquor_dm <- bind_rows(dm_row(dm_arima_ets, "ARIMA versus ETS"),
                      dm_row(dm_arima_snaive, "ARIMA versus seasonal naive"))
print(liquor_dm)
saveRDS(liquor_errs, "results/liquor_errs.rds")
saveRDS(liquor_dm, "results/liquor_dm.rds")

liquor_combo_fc <- liquor_cv_fit |>
  mutate(combo = (arima + ets) / 2) |> forecast(h = 1)
combo_accuracy <- liquor_combo_fc |>
  accuracy(liquor, measures = list(RMSE = RMSE, MAE = MAE))
error_correlation <- cor(liquor_errs$arima, liquor_errs$ets)
print(combo_accuracy)
print(error_correlation)
saveRDS(combo_accuracy, "results/combo_accuracy.rds")
saveRDS(error_correlation, "results/error_correlation.rds")
message("Parts 2 and 3 finished; saved model fits, forecasts, and accuracy.")

qd_raw <- read_csv("data/fredqd_2026-09.csv",
  col_types = cols(sasdate = col_character(), .default = col_double()),
  show_col_types = FALSE) |>
  slice(-(1:2))
growth <- function(x) 100 * difference(log(x))
qd <- qd_raw |>
  mutate(date = as.Date(sasdate, "%m/%d/%Y")) |>
  transmute(date, gdp = growth(GDPC1), cons = growth(PCECC96),
            inv = growth(GPDIC1), emp = growth(PAYEMS),
            infl = growth(GDPCTPI), unrate = difference(UNRATE),
            ffr = difference(FEDFUNDS), gs10 = difference(GS10)) |>
  filter(date >= as.Date("1960-03-01"), date <= as.Date("2019-12-01"))
y <- qd |> select(-date) |> as.matrix()
print(nrow(y))
print(head(y))
stopifnot(identical(dim(y), c(240L, 8L)), all(is.finite(y)))
N <- nrow(y)
origins <- 120:(N - 1)
actual <- y[121:N, ]

var_ols_fc <- function(z, p) {
  e <- embed(z, p + 1)
  M <- ncol(z)
  B <- coef(lm(e[, 1:M] ~ e[, -(1:M)]))
  x <- c(1, t(z[nrow(z):(nrow(z) - p + 1), ]))
  x %*% B
}
fc_var1 <- NULL
fc_var4 <- NULL
for (T in origins) {
  fc_var1 <- rbind(fc_var1, var_ols_fc(y[1:T, ], 1))
  fc_var4 <- rbind(fc_var4, var_ols_fc(y[1:T, ], 4))
}
rmse_by_var <- function(f) sqrt(colMeans((actual - f)^2))
print(rmse_by_var(fc_var1))
print(rmse_by_var(fc_var4))
saveRDS(list(actual = actual, fc_var1 = fc_var1, fc_var4 = fc_var4,
             dates = qd$date[121:N]), "results/ols_var_forecasts.rds")

bvar_cv <- function(lambda) {
  pr <- bv_priors(hyper = "lambda",
                  mn = bv_mn(lambda = bv_lambda(mode = lambda, sd = 0.4,
                                                min = 0.999 * lambda,
                                                max = 1.001 * lambda),
                             b = 0))
  fc <- NULL
  for (T in origins) {
    fit_T <- bvar(y[1:T, ], lags = 4, priors = pr,
                  n_draw = 2000, n_burn = 1000, verbose = FALSE)
    x <- c(1, t(y[T:(T - 3), ]))
    fc <- rbind(fc, x %*% coef(fit_T))
  }
  fc
}
message("Estimating BVAR(4), lambda = 0.3.")
set.seed(422)
fc_bvar_03 <- bvar_cv(0.3)
saveRDS(fc_bvar_03, "results/fc_bvar_03.rds")

message("Estimating BVAR(4), lambda = 0.05.")
set.seed(422)
fc_bvar_005 <- bvar_cv(0.05)
saveRDS(fc_bvar_005, "results/fc_bvar_005.rds")
var_results <- tibble(
  Variable = colnames(y), `VAR(1)` = rmse_by_var(fc_var1),
  `VAR(4)` = rmse_by_var(fc_var4),
  `BVAR(4), lambda = 0.3` = rmse_by_var(fc_bvar_03),
  `BVAR(4), lambda = 0.05` = rmse_by_var(fc_bvar_005)
)
print(var_results |> knitr::kable(digits = 3))
saveRDS(var_results, "results/var_results.rds")

e_var4 <- actual[, "gdp"] - fc_var4[, 1]
e_bvar <- actual[, "gdp"] - fc_bvar_03[, 1]
dm_gdp <- forecast::dm.test(e_var4, e_bvar, h = 1, power = 2)
e_var4_cons <- actual[, "cons"] - fc_var4[, 2]
e_bvar_cons <- actual[, "cons"] - fc_bvar_03[, 2]
dm_cons <- forecast::dm.test(e_var4_cons, e_bvar_cons, h = 1, power = 2)
var_dm <- bind_rows(dm_row(dm_gdp, "GDP: VAR(4) versus BVAR(4), lambda = 0.3"),
                   dm_row(dm_cons, "Consumption: VAR(4) versus BVAR(4), lambda = 0.3"))
print(var_dm)
saveRDS(var_dm, "results/var_dm.rds")
message("Part 4 finished.")

md_raw <- read_csv("data/fredmd_2026-09.csv",
  col_types = cols(sasdate = col_character(), .default = col_double()),
  show_col_types = FALSE)
tcode <- unlist(md_raw[1, -1])
md <- md_raw[-1, ] |> mutate(date = as.Date(sasdate, "%m/%d/%Y"))
transform_series <- function(x, code) {
  d <- function(v) c(NA, diff(v))
  switch(as.character(code),
         "1" = x, "2" = d(x), "3" = d(d(x)), "4" = log(x),
         "5" = 100 * d(log(x)), "6" = 100 * d(d(log(x))),
         "7" = d(x / dplyr::lag(x) - 1))
}
fm <- map2_dfc(md[names(tcode)], tcode, transform_series) |>
  mutate(date = md$date, .before = 1) |>
  filter(date >= as.Date("1960-03-01"), date <= as.Date("2019-12-01"))
fm <- fm[, colSums(is.na(fm)) == 0]
print(dim(fm))
print(tcode[c("INDPRO", "PAYEMS", "CPIAUCSL", "FEDFUNDS")])
stopifnot(identical(dim(fm), c(718L, 122L)))

Z <- fm |> select(-date) |> as.matrix()
pc <- prcomp(scale(Z))
factor_importance <- summary(pc)$importance[, 1:5]
print(factor_importance)
F1 <- pc$x[, 1] * sign(cor(pc$x[, 1], fm$INDPRO))
factor_plot <- tibble(date = fm$date, F1 = F1) |>
  ggplot(aes(date, F1)) + geom_line(colour = "grey20", linewidth = 0.4) +
  labs(title = "First Factor from FRED-MD",
       subtitle = "Sign chosen so that the factor rises with industrial production growth",
       x = NULL, y = "Factor score (arbitrary scale)") + theme_minimal()
print(factor_plot)
ggsave("results/first_factor.png", factor_plot, width = 9, height = 4.8, dpi = 160)
factor_top <- map_dfr(1:2, function(k) {
  vals <- sort(abs(cor(pc$x[, k], Z))[1, ], decreasing = TRUE)[1:8]
  tibble(Factor = k, Series = names(vals), Absolute_correlation = unname(vals))
})
print(factor_top)
saveRDS(list(fm = fm, tcode = tcode, importance = factor_importance,
             F1 = F1, top = factor_top), "results/factor_descriptive.rds")

message("Running the industrial production factor horse race.")
ip <- fm$INDPRO
n <- nrow(Z)
init <- which(fm$date == as.Date("1989-12-01"))
e_ar <- NULL
e_fa <- NULL
for (T in init:(n - 1)) {
  F_T <- prcomp(scale(Z[1:T, ]))$x[, 1:3]
  d <- data.frame(y_next = ip[2:T], y = ip[1:(T - 1)], F_T[1:(T - 1), ])
  new <- data.frame(y = ip[T], t(F_T[T, ]))
  e_ar <- c(e_ar, ip[T + 1] - predict(lm(y_next ~ y, d), new))
  e_fa <- c(e_fa, ip[T + 1] - predict(lm(y_next ~ ., d), new))
}
ip_results <- tibble(Model = c("AR(1)", "AR(1) + three factors"),
                     RMSE = c(sqrt(mean(e_ar^2)), sqrt(mean(e_fa^2))),
                     Forecasts = c(length(e_ar), length(e_fa)))
print(ip_results)
saveRDS(list(table = ip_results, e_ar = e_ar, e_fa = e_fa),
        "results/ip_results.rds")

direct_race <- function(series, h = 12, r = 3) {
  lv <- md |>
    mutate(log_x = log(.data[[series]]),
           g = 100 * (log_x - lag(log_x))) |>
    filter(date >= as.Date("1960-03-01"), date <= as.Date("2019-12-01"))
  stopifnot(identical(lv$date, fm$date))
  yh <- 100 * (lead(lv$log_x, h) - lv$log_x)
  g <- lv$g
  e_ar <- NULL
  e_fa <- NULL
  for (T in init:(n - h)) {
    F_T <- prcomp(scale(Z[1:T, ]))$x[, 1:r, drop = FALSE]
    i <- 1:(T - h)
    d <- data.frame(yh = yh[i], g = g[i], F_T[i, , drop = FALSE])
    new <- data.frame(g = g[T], F_T[T, , drop = FALSE])
    e_ar <- c(e_ar, yh[T] - predict(lm(yh ~ g, d), new))
    e_fa <- c(e_fa, yh[T] - predict(lm(yh ~ ., d), new))
  }
  test <- forecast::dm.test(e_ar, e_fa, h = h, power = 2,
                            varestimator = "bartlett")
  saveRDS(list(e_ar = e_ar, e_fa = e_fa, test = test,
               origin_dates = lv$date[init:(n - h)],
               endpoint_dates = lv$date[(init + h):n]),
          paste0("results/direct_errors_", series, ".rds"))
  tibble(series, forecasts = length(e_ar),
         RMSE_AR = sqrt(mean(e_ar^2)), RMSE_FA = sqrt(mean(e_fa^2)),
         MAE_AR = mean(abs(e_ar)), MAE_FA = mean(abs(e_fa)),
         DM_pvalue = test$p.value)
}
message("Running direct twelve-month CPI forecasts.")
cpi_results <- direct_race("CPIAUCSL")
print(cpi_results)
saveRDS(cpi_results, "results/cpi_results.rds")
print(100 * (1 - cpi_results$RMSE_FA / cpi_results$RMSE_AR))

writeLines(capture.output(sessionInfo()), "results/session_info.txt")
message("All estimation finished. Render ps4_writeup_dugan.qmd next.")
