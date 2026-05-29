### WorkSpace

rm(list = ls())

library(dplyr)
library(mgcv)
library(dlnm)
library(lubridate)
library(sensobol)
library(doParallel)
library(foreach)
library(parallel)

# Load anonymized/randomized dataset
db_anon_MJJAS <- readRDS("db_anon_MJJAS.rds")

# Standardize mean temperature
db_anon_MJJAS$t_mean <- as.numeric(scale(db_anon_MJJAS$t_mean))

# Remove May, keeping June-September
db_anon_JJAS <- db_anon_MJJAS %>% 
  filter(month(dataric) != 5)

# Remove September, keeping June-August
db_anon_JJA <- db_anon_JJAS %>% 
  filter(month(dataric) != 9)

head(db_anon_JJAS)


# Definizione nodi --------------------------------------------------------

# Nodi temperatura con September

temp_ord <- sort(db_anon_JJAS$tmax)

lista_knots_var <- list()
possibili_df_var <- 3:10

for(df in seq_along(possibili_df_var)){
  lista_knots_var[[df]] <- temp_ord[
    floor(length(temp_ord) / possibili_df_var[df] * (1:(possibili_df_var[df] - 1)))
  ]
}

names(lista_knots_var) <- paste0(possibili_df_var, " df")

# Additional knot specification based on selected quantiles
lista_knots_var[[9]] <- as.vector(quantile(db_anon_JJAS$tmax, c(.5, .9)))
names(lista_knots_var)[9] <- "3 df - quantile"

lista_knots_var


# Nodi lag

lista_knots_lag <- list()
possibili_df_lag <- 3:10

for(df in seq_along(possibili_df_lag)){
  lista_knots_lag[[df]] <- equalknots(
    0:14, 
    df = possibili_df_lag[df], 
    intercept = TRUE
  )
}

names(lista_knots_lag) <- paste0(possibili_df_lag, " df")

lista_knots_lag


# GSA ---------------------------------------------------------------------

rm(db_anon_JJA, db_anon_JJAS)

N <- 2^9

params <- c("argvar", "arglag", "september")

mat_sob <- sobol_matrices(N = N, params = params)

mat_sob[, 1] <- ceiling(mat_sob[, 1] * 9)
mat_sob[, 2] <- ceiling(mat_sob[, 2] * 8)
mat_sob[, 3] <- round(mat_sob[, 3], 0)

mat_sob[, 1] <- pmax(mat_sob[, 1], 1)
mat_sob[, 2] <- pmax(mat_sob[, 2], 1)

head(mat_sob)


# Parallel setup ----------------------------------------------------------

n_cores <- ceiling(detectCores() * 0.75)

cl <- makeCluster(n_cores)
registerDoParallel(cl)


# Run GSA models ----------------------------------------------------------

Y <- foreach(
  i = 1:nrow(mat_sob),
  .combine = rbind,
  .packages = c("dlnm", "mgcv", "lubridate"),
  .export = c(
    "db_anon_MJJAS", 
    "lista_knots_var", 
    "lista_knots_lag", 
    "mat_sob"
  )
) %dopar% {
  
  data_full <- db_anon_MJJAS
  
  formulamod <- as.formula(
    n ~ cb + 
      factor(holiday1) + 
      factor(covid) + 
      year_ + 
      t_mean + 
      s(seasonality) + 
      offset(log(pop)) + 
      weekend + 
      s(COD_RIP, bs = "re")
  )
  
  if (mat_sob[i, 3] == 1) {
    idx_analysis <- month(data_full$dataric) %in% 6:8
  } else {
    idx_analysis <- month(data_full$dataric) %in% 6:9
  }
  
  data <- data_full[idx_analysis, , drop = FALSE]
  data$COD_RIP <- factor(data$COD_RIP)
  
  cb_full <- crossbasis(
    data_full$tmax,
    lag = 14,
    group = list(data_full$comres, data_full$year),
    argvar = list(
      fun = "ns",
      knots = lista_knots_var[[mat_sob[i, 1]]],
      intercept = FALSE
    ),
    arglag = list(
      fun = "ns",
      knots = lista_knots_lag[[mat_sob[i, 2]]],
      intercept = TRUE
    )
  )
  
  cb <- cb_full[idx_analysis, , drop = FALSE]
  
  att <- attributes(cb_full)
  att$dim <- dim(cb)
  att$dimnames <- dimnames(cb)
  attributes(cb) <- att
  
  time <- Sys.time()
  
  mod_nb <- gam(
    formula = formulamod,
    family = nb(),
    data = data
  )
  
  time_mod <- as.numeric(Sys.time() - time)
  
  pred.temp <- crosspred(
    basis = cb,
    model = mod_nb,
    from = 1,
    to = floor(max(data$tmax)),
    by = 1,
    cen = 30
  )
  
  c(
    pred.temp$allRRfit,
    aic = mod_nb$aic,
    dev = mod_nb$deviance,
    time = time_mod,
    pred.temp$allRRlow,
    pred.temp$allRRhigh
  )
}

stopCluster(cl)


# Plot GSA ----------------------------------------------------------------

sob <- function(Y, N, params) {
  
  sob <- sobol_indices(
    Y = Y, 
    N = N, 
    params = params
  )
  
  plot(sob)
}

apply(Y, 2, sob, N = N, params = params)
