pacman::p_load(
    Ipaper, data.table, dplyr, ggplot2, gg.layers, lubridate
)

read_flux <- function(f) {
    vars <- fread(f, nrows = 1) |> names()
    fread(f, skip = 2) |> set_names(vars)
}

f = "data/OUTPUT_栾城_spin300y.csv"
d_obs <- read_flux("data/CRO_栾城_Day_Flux_200710-201809.csv") |> 
    mutate(date = date(date))
d_sim <- fread(f) |> rename(date = dates)
dat = merge(d_obs, d_sim, by = "date", suffixes = c("_obs", "_sim")) |> 
    filter(year(date) <= 2016)

# ET模拟结果
fmt_gof = "*KGE* = {str_num(KGE,2)}, *NSE* = {str_num(NSE,2)}, *R^2* = {str_num(R2, 2)} \n *RMSE* = {str_num(RMSE,2)}"
pdat = dat[, .(date, ET_obs, ET_sim)] |> melt("date")
p <- ggplot(pdat) + 
    geom_line(aes(date, value, color = variable)) + 
    geom_gof2(data = dat, aes(obs = ET_obs, sim = ET_sim), label.format = fmt_gof) + 
    labs(y = "Evapotranspiration (mm/d)", x = NULL)
write_fig(p, 'Figures/Figure1_栾城_ET.pdf', 10, 5)

# GW模拟结果
p <- ggplot(d_sim[GW>0]) + 
    geom_line(aes(date, GW)) 
write_fig(p, 'Figures/Figure1_栾城_GW.pdf', 10, 5)

# SM
pdat = d_sim[GW>0, .(date, SM1, SM2, SM3)] |> melt("date")
p <- ggplot(pdat) +
    geom_line(aes(date, value)) + 
    facet_wrap(~variable, scales = "free_y")
write_fig(p, "Figures/Figure1_栾城_SM.pdf", 10, 5)
