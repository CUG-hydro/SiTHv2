pacman::p_load(
  Ipaper, dplyr, ggplot2, lubridate, data.table
)

na_approx <- function(x, y) {
  inds_bad <- which(is.na(y))
  inds_good <- which(!is.na(y))
  y[inds_bad] <- approx(x[inds_good], y[inds_good], xout = x[inds_bad])$y
  y
}

df <- fread("data/dat_栾城_ERA5L_1982-2019_raw.csv") |>
  rename(LAI_raw = LAI, VOD = VOD.x) |>
  mutate(
    date = ymd(date),
    Rn = Rn / 86400, # [J day-1 m-2] -> [W m-2]
    Tavg = Tavg - 273.15, # [K] -> [C]
    Prcp = Prcp * 1000, # [m] -> [mm]
    Pa = Pa / 1000, # [Pa] -> [kPa]
    LAI = na_approx(date, LAI_raw)
  ) |> 
  select(-LAI_raw, -VOD.y)
fwrite(df, "data/dat_栾城_ERA5L_1982-2019.csv")

dat <- melt(df |> select(-LULC), "date")

## 绘图检查输入数据
p <- ggplot(dat, aes(date, value)) +
  geom_line() +
  facet_wrap(~variable, scales = "free_y")
write_fig(p, "Figure1_Forcing_time-series.pdf", 10, 5)
