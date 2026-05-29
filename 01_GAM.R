# Data manipulation
library(dplyr)

# Function month()
library(lubridate)

# crossbasis(), crosspred(), equalknots()
library(dlnm)

# gam(), s(), family nb
library(mgcv)

# Plot
library(ggplot2)


# Create an randomized version of the dataset by randomly permuting the outcome n
# db_anon_MJJAS <- db_agg_capreg_may %>%
#   ungroup() %>% 
#   mutate(n = sample(n, size = n(), replace = FALSE)) 
# saveRDS(db_anon_MJJAS, file = "db_anon_MJJAS.rds")

# Load randomized version of the data
db_anon_MJJAS <- readRDS("db_anon_MJJAS.rds")

# Standardize mean temperature
db_anon_MJJAS$t_mean = as.numeric(scale(db_anon_MJJAS$t_mean))

# Keep only June-September observations for model fitting
db_anon_JJAS = db_anon_MJJAS %>% filter(month(dataric) %in% 6:9)

# Define temperature knots for the exposure-response function
temp_ord = sort(db_anon_JJAS$tmax)
lista_knots_var = list()
possibili_df_var = 3:10

for(df in seq_along(possibili_df_var)){
  lista_knots_var[[df]] = temp_ord[floor(length(temp_ord)/possibili_df_var[df]*(1:(possibili_df_var[df]-1)))]
}

names(lista_knots_var) = paste0(possibili_df_var, " df")

# Add an alternative knot specification based on selected quantiles
lista_knots_var[[9]] = as.vector(quantile(db_anon_JJAS$tmax, c(.5, .9)))
names(lista_knots_var)[9] = "3 df - quantile"

# lista_knots_var

# Define lag knots for the lag-response function
lista_knots_lag = list()
possibili_df_lag = 3:10

for(df in seq_along(possibili_df_lag)){
  lista_knots_lag[[df]] = equalknots(0:14, df = possibili_df_lag[df], intercept = T)
}

names(lista_knots_lag) = paste0(possibili_df_lag, " df")

# lista_knots_lag


# CROSSBASIS MATRIX 

# Build the full crossbasis matrix for temperature and lags up to 14 days
ns.basis_full <- crossbasis(db_anon_MJJAS$tmax,
                            group = list(db_anon_MJJAS$comres, db_anon_MJJAS$year),
                            lag = 14,
                            argvar = list(fun = "ns", knots = lista_knots_var[["5 df"]], intercept = F), #knots = temp_rif[c(length(temp_rif)/3, length(temp_rif)/3*2)]),
                            arglag = list(fun = "ns", knots = lista_knots_lag[["4 df"]], intercept = T))

# Select the rows corresponding to June-September
idx_JJASs <- month(db_anon_MJJAS$dataric) %in% 6:9
ns.basis = ns.basis_full[idx_JJASs, , drop = FALSE]

# Restore the attributes needed by crosspred after subsetting the crossbasis
att <- attributes(ns.basis_full)
att$dim <- dim(ns.basis)
att$dimnames <- dimnames(ns.basis)
attributes(ns.basis) <- att

# Fit a negative binomial GAM using the anonymized outcome
time <- Sys.time(); mod_nb = gam(formula = n ~ ns.basis + holiday1 + covid + year_ + weekend+
                                   s(seasonality) + offset(log(pop)) + t_mean + s(COD_RIP, bs = "re"),
                                 family = "nb",
                                 link = "log",
                                 data = db_anon_JJAS %>% 
                                   mutate(COD_RIP = factor(COD_RIP))); Sys.time()-time; #summary(mod_nb)

# Predict the overall cumulative exposure-response curve
pred.temp <- crosspred(ns.basis, mod_nb, from = 1, by=1, cen = 30, ci.level=.9)

# Plot the estimated relative risk curve
data.frame(temp = pred.temp$predvar,
           RR = pred.temp$allRRfit,
           RR_l = pred.temp$allRRlow,
           RR_u = pred.temp$allRRhigh) %>% 
  ggplot(aes(x = temp,
             y = RR,
             ymin = RR_l,
             ymax = RR_u))+
  geom_line()+
  geom_ribbon(alpha = .2)+
  theme_minimal()+
  labs(x = "Daily maximum temperature (°C)",
       caption = "(b)")+
  # scale_y_continuous(limits = c(0, 7))+
  geom_hline(aes(yintercept = 1), linetype = "dashed") 
