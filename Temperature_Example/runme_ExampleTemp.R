# Example Sequential Consensus splitting data.frame ----

library(INLA)
library(inlabru)
library(ggplot2)
library(ggtext)
library(gridExtra)
library(dplyr)
library(sf)
library(leaflet)
library(tidyverse)

## Custom functions ----

## Loading data ----

DFsf_temp <- readRDS(file = "./Temperature_example/DF_temp.RDS")
mesh <- readRDS(file = "./Temperature_example/mesh.RDS")

## Exploring data ----

DF_temp <- data.frame(cbind(DFsf_temp %>% st_drop_geometry(.), DFsf_temp %>% st_coordinates(.)))

plot_func <- function(df, idx, title = NULL, name = "turbo") {
  ggplot(data = df[idx,], aes(x = X, y = Y, fill = max_sst)) +
    geom_tile() + labs(title = title) +
    scale_fill_viridis_c(option = name) + theme_minimal()
}

list_plots <- lapply(X = 1:36, FUN = function(i){plot_func(df = DF_temp, idx = which(DF_temp$month.id %in% i), title = paste("Month",i), name = "turbo")})

grid.arrange(arrangeGrob(grobs = list_plots, ncol = 6, nrow = 6))

DF_temp$month.id %>% max(.)
ntime <- 480
k_times <- 12
DFsimplified_temp <- DF_temp[which(DF_temp$month.id %in% 1:ntime),]

# Sequential consensus method 1 ----

latent_random_field <- list()
marginal.latent_random_field <- list()
mean.intercept <- c()
t_cpu <- c()
for(i in 1:k_times){
  # We start defining the
  spde <- inla.spde2.pcmatern(
    mesh = mesh, alpha = 2,
    prior.range = c(mesh$loc[c(mesh$segm$bnd$idx[,1], mesh$segm$bnd$idx[1,1]),1:2] %>% dist(x = .) %>% as.vector(.) %>% max(.)/5, 0.5),
    prior.sigma = c(1,0.5), constr = TRUE
  )
  idx <- DFsimplified_temp$month.id %in% ((ntime/k_times*(i-1)+1):(ntime/k_times*i))
  spde.idx <- inla.spde.make.index(name = "sp.idx", n.spde = spde$n.spde, n.group = ntime/k_times)
  ntimes_group <- DFsimplified_temp$month.id[idx] %>% table(.) %>% as.vector(.)
  A_spt <- inla.spde.make.A(mesh = mesh, loc = DFsf_temp[which(idx),] %>% st_coordinates(.), group = rep(1:(ntime/k_times), times = ntimes_group))
  
  stk_inf_i <- inla.stack(data = list(y = DFsimplified_temp$max_sst[idx]),
                          A = list(A_spt, 1),
                          effects = list(
                            spde.idx,
                            list(intercept = rep(1, times = sum(idx)))
                          ),
                          tag = paste0("stk_inf_i"),
                          compress = TRUE,
                          remove.unused = TRUE)
  
  if(i==1){
    t_cpu[i] <- Sys.time()
    formula_inla_i <- y ~ -1 + intercept + f(sp.idx, model = spde, group = sp.idx.group, control.group = list(model = "ar1"))
    model_inla_i <- inla(formula = formula_inla_i, family = "gaussian", 
                         data = inla.stack.data(stk_inf_i),
                         control.predictor = list(A = inla.stack.A(stk_inf_i)), 
                         control.compute = list(config = TRUE, waic = TRUE), 
                         verbose = FALSE)
  } else{
    mode.hyper <- model_inla_i$internal.summary.hyperpar$mode
    fixed.mean <- model_inla_i$summary.fixed$mean %>% list(.); fixed.prec <- model_inla_i$summary.fixed$sd**-2 %>% list(.)
    names(fixed.mean) <- names(fixed.prec) <- rownames(model_inla_i$summary.fixed)
    
    tab_internal.gaussian.prec <- paste0("table: ",
                                         paste(c(model_inla_i$internal.marginals.hyperpar[[1]] %>%
                                                   inla.smarginal(marginal = ., log = TRUE) %>% unlist(.) %>% as.numeric(.))),
                                         collapse = " ")
    
    tab_internal.sp.rho <- paste0("table: ",
                                  paste(c(model_inla_i$internal.marginals.hyperpar[[4]] %>%
                                            inla.smarginal(marginal = ., log = TRUE) %>% unlist(.) %>% as.numeric(.))),
                                  collapse = " ")
    
    alpha <- 2; d <- 2; nu <- alpha - d/2
    kappa0 <- log(8*nu)/2
    tau0 <- 0.5*(lgamma(nu) - lgamma(nu + d/2) - d/2*log(4*pi)) - nu*kappa0
    spde <- inla.spde2.matern(mesh = mesh, B.tau = cbind(tau0, nu, -1), 
                              B.kappa = cbind(kappa0, -1, 0), alpha = 2,
                              theta.prior.mean = model_inla_i$internal.summary.hyperpar[3:2,1], 
                              theta.prior.prec = model_inla_i$internal.summary.hyperpar[3:2,2]**-2, 
                              constr = TRUE)
    
    formula_inla_i <- y ~ -1 + intercept + f(sp.idx, model = spde, group = sp.idx.group, control.group = list(model = "ar1", hyper = list("logit correlation" = list(tab_internal.sp.rho))))
    model_inla_i <- inla(formula = formula_inla_i, family = "gaussian", data = inla.stack.data(stk_inf_i), 
                         control.predictor = list(A = inla.stack.A(stk_inf_i)), 
                         control.fixed = list(mean = fixed.mean, prec = fixed.prec),
                         control.family = list(hyper = list(theta1 = list(tab_internal.gaussian.prec))),
                         control.mode = list(theta = mode.hyper),
                         control.compute = list(config = TRUE, waic = TRUE), verbose = FALSE)
  }
  mean.intercept[i] <- model_inla_i$summary.fixed$mean[1]
  marginal.latent_random_field[[i]] <- model_inla_i$marginals.random
  latent_random_field[[i]] <- model_inla_i$summary.random
  t_cpu[i+1] <- Sys.time()
  print(paste("Model", i, "is done."))
}

# Sequential consensus method 2 ----

latent_random_field <- list()
marginal.latent_random_field <- list()
mean.intercept <- c()
latent_fixed_effects <- list()
marginal.latent_fixed_effects <- list()
t_cpu <- c()
for(i in 1:k_times){
  # We start defining the
  rho0 <- mesh$loc[c(mesh$segm$bnd$idx[,1], mesh$segm$bnd$idx[1,1]),1:2] %>% dist(x = .) %>% as.vector(.) %>% max(.)/5
  sigma0 <- 1
  alpha <- 2; d <- 2; nu <- alpha - d/2
  kappa0 <- log(8*nu)/2 - log(rho0)
  tau0 <- 0.5*(lgamma(nu) - lgamma(nu + d/2) - d/2*log(4*pi)) - log(sigma0) - nu*kappa0
  spde <- inla.spde2.matern(mesh = mesh, B.tau = cbind(tau0, nu, -1), 
                            B.kappa = cbind(kappa0, -1, 0), alpha = 2, constr = FALSE)
  
  idx <- DFsimplified_temp$month.id %in% ((ntime/k_times*(i-1)+1):(ntime/k_times*i))
  spde.idx <- inla.spde.make.index(name = "sp.idx", n.spde = spde$n.spde, n.group = ntime/k_times)
  ntimes_group <- DFsimplified_temp$month.id[idx] %>% table(.) %>% as.vector(.)
  A_spt <- inla.spde.make.A(mesh = mesh, loc = DFsf_temp[which(idx),] %>% st_coordinates(.), group = rep(1:(ntime/k_times), times = ntimes_group))
  
  stk_inf_i <- inla.stack(data = list(y = DFsimplified_temp$max_sst[idx]),
                          A = list(A_spt, 1),
                          effects = list(
                            spde.idx,
                            list(intercept = rep(1, times = sum(idx)))
                          ),
                          tag = paste0("stk_inf_i"),
                          compress = TRUE,
                          remove.unused = TRUE)
  
  if(i==1){
    t_cpu[i] <- Sys.time()
    formula_inla_i <- y ~ -1 + intercept + f(sp.idx, model = spde, group = sp.idx.group, control.group = list(model = "ar1"))
    model_inla_i <- inla(formula = formula_inla_i, family = "gaussian", 
                         data = inla.stack.data(stk_inf_i),
                         control.predictor = list(A = inla.stack.A(stk_inf_i)), 
                         control.compute = list(config = TRUE, waic = TRUE), 
                         verbose = FALSE)
  } else{
    mode.hyper <- model_inla_i$internal.summary.hyperpar$mode
    fixed.mean <- model_inla_i$summary.fixed$mean %>% list(.); fixed.prec <- model_inla_i$summary.fixed$sd**-2 %>% list(.)
    names(fixed.mean) <- names(fixed.prec) <- rownames(model_inla_i$summary.fixed)
    
    norm_internal.gaussian.prec <- list(prior = "normal", param = c(model_inla_i$internal.summary.hyperpar[1,1], model_inla_i$internal.summary.hyperpar[1,2]**-2), fixed = FALSE)
    norm_internal.sp.rho <- list(prior = "normal", param = c(model_inla_i$internal.summary.hyperpar[4,1], model_inla_i$internal.summary.hyperpar[4,2]**-2), fixed = FALSE)
    
    spde <- inla.spde2.matern(mesh = mesh, B.tau = cbind(tau0, nu, -1), 
                              B.kappa = cbind(kappa0, -1, 0), alpha = 2,
                              theta.prior.mean = model_inla_i$internal.summary.hyperpar[2:3,1], 
                              theta.prior.prec = model_inla_i$internal.summary.hyperpar[2:3,2]**-2)
    
    formula_inla_i <- y ~ -1 + intercept + f(sp.idx, model = spde, group = sp.idx.group, control.group = list(model = "ar1", hyper = list("logit correlation" = list(norm_internal.sp.rho))))
    model_inla_i <- inla(formula = formula_inla_i, family = "gaussian", data = inla.stack.data(stk_inf_i), 
                         control.predictor = list(A = inla.stack.A(stk_inf_i)), 
                         control.fixed = list(mean = fixed.mean, prec = fixed.prec),
                         control.family = list(hyper = list(theta1 = list(norm_internal.gaussian.prec))),
                         control.mode = list(theta = mode.hyper),
                         control.compute = list(config = TRUE, waic = TRUE), verbose = FALSE)
  }
  mean.intercept[i] <- model_inla_i$summary.fixed$mean[1]
  marginal.latent_random_field[[i]] <- model_inla_i$marginals.random
  latent_random_field[[i]] <- model_inla_i$summary.random
  marginal.latent_fixed_effects[[i]] <- model_inla_i$marginals.fixed
  latent_fixed_effects[[i]] <- model_inla_i$summary.fixed
  t_cpu[i+1] <- Sys.time()
  print(paste("Model", i, "is done."))
}

t_cpu %>% diff(.) %>% sum(.)/60

# Sequential consensus method 2 (re-fitting) ----
# In this case we will consider that pi(beta|y_{-i})=prod(pi(beta|y_{1:i}))prod(pi(y_{(i+1):n}|beta))

latent_random_field_rf <- list()
marginal.latent_random_field_rf <- list()
t_cpu_rf <- c()
for(i in 1:k_times){
  # We start defining the
  idx <- DFsimplified_temp$month.id %in% ((ntime/k_times*(i-1)+1):(ntime/k_times*i))
  spde.idx <- inla.spde.make.index(name = "sp.idx", n.spde = spde$n.spde, n.group = ntime/k_times)
  ntimes_group <- DFsimplified_temp$month.id[idx] %>% table(.) %>% as.vector(.)
  A_spt <- inla.spde.make.A(mesh = mesh, loc = DFsf_temp[which(idx),] %>% st_coordinates(.), group = rep(1:(ntime/k_times), times = ntimes_group))
  
  if(i==1){
    fixed.mean_corrected <- ((latent_fixed_effects[[k_times]]$sd**-2 - latent_fixed_effects[[i]]$sd**-2)**(-1)*(latent_fixed_effects[[k_times]]$sd**(-2)*latent_fixed_effects[[k_times]]$mean - latent_fixed_effects[[i]]$sd**(-2)*latent_fixed_effects[[i]]$mean)) %>% list(.)
    fixed.prec_corrected <- (latent_fixed_effects[[k_times]]$sd**-2 - latent_fixed_effects[[i]]$sd**-2) %>% list(.)
    names(fixed.mean_corrected) <- names(fixed.prec_corrected) <- rownames(model_inla_i$summary.fixed)
  } else if(i==k_times){
    fixed.mean_corrected <- latent_fixed_effects[[i-1]]$mean %>% list(.)
    fixed.prec_corrected <- latent_fixed_effects[[i-1]]$sd**-2 %>% list(.)
    names(fixed.mean_corrected) <- names(fixed.prec_corrected) <- rownames(model_inla_i$summary.fixed)
  } else{
    fixed.mean_corrected <- ((latent_fixed_effects[[i-1]]$sd**-2 + latent_fixed_effects[[k_times]]$sd**-2 - latent_fixed_effects[[i]]$sd**-2)**(-1)*(latent_fixed_effects[[i-1]]$sd**(-2)*latent_fixed_effects[[i-1]]$mean + latent_fixed_effects[[k_times]]$sd**(-2)*latent_fixed_effects[[k_times]]$mean - latent_fixed_effects[[i]]$sd**(-2)*latent_fixed_effects[[i]]$mean)) %>% list(.)
    fixed.prec_corrected <- (latent_fixed_effects[[i-1]]$sd**-2 + latent_fixed_effects[[k_times]]$sd**-2 - latent_fixed_effects[[i]]$sd**-2) %>% list(.)
    names(fixed.mean_corrected) <- names(fixed.prec_corrected) <- rownames(model_inla_i$summary.fixed)
  }
  
  stk_inf_i <- inla.stack(data = list(y = DFsimplified_temp$max_sst[idx]),
                          A = list(A_spt, 1),
                          effects = list(
                            spde.idx,
                            list(intercept = rep(1, times = sum(idx)),
                                 intercept_off = rep(fixed.mean$intercept, times = sum(idx)))
                          ),
                          tag = paste0("stk_inf_i"),
                          compress = TRUE,
                          remove.unused = TRUE)
  
  if(i==1){
    t_cpu_rf[i] <- Sys.time()
    formula_inla_i <- y ~ -1 + intercept + f(sp.idx, model = spde, group = sp.idx.group, control.group = list(model = "ar1"))
    # formula_inla_i <- y ~ -1 + offset(intercept_off) + f(sp.idx, model = spde, group = sp.idx.group, control.group = list(model = "ar1"))
    model_inla_i_rf <- inla(formula = formula_inla_i, family = "gaussian", 
                            data = inla.stack.data(stk_inf_i),
                            control.fixed = list(mean = fixed.mean_corrected, prec = fixed.prec_corrected),
                            # control.fixed = list(mean = fixed.mean, prec = fixed.prec),
                            control.predictor = list(A = inla.stack.A(stk_inf_i)), 
                            control.compute = list(config = TRUE, waic = TRUE),
                            control.inla = list(int.strategy = "user.expert", int.design = model_inla_i$joint.hyper),
                            inla.mode = "compact",
                            control.mode = list(theta = model_inla_i$misc$configs$config[[1]]$theta),
                            # control.mode = list(theta = model_inla_i$misc$configs$config[[1]]$theta, x = latent_random_field[[i]]$sp.idx$mode),
                            verbose = FALSE)
  } else{
    formula_inla_i <- y ~ -1 + intercept + f(sp.idx, model = spde, group = sp.idx.group, control.group = list(model = "ar1"))
    # formula_inla_i <- y ~ -1 + offset(intercept_off) + f(sp.idx, model = spde, group = sp.idx.group, control.group = list(model = "ar1", hyper = list("logit correlation" = list(tab_internal.sp.rho))))
    model_inla_i_rf <- inla(formula = formula_inla_i, family = "gaussian", data = inla.stack.data(stk_inf_i), 
                            control.predictor = list(A = inla.stack.A(stk_inf_i)), 
                            control.fixed = list(mean = fixed.mean_corrected, prec = fixed.prec_corrected),
                            # control.fixed = list(mean = fixed.mean, prec = fixed.prec),
                            control.compute = list(config = TRUE, waic = TRUE), 
                            control.inla = list(int.strategy = "user.expert", int.design = model_inla_i$joint.hyper),
                            inla.mode = "compact",
                            control.mode = list(theta = model_inla_i$misc$configs$config[[1]]$theta),
                            # control.mode = list(theta = model_inla_i$misc$configs$config[[1]]$theta, x = latent_random_field[[i]]$sp.idx$mode),
                            verbose = FALSE)
  }
  marginal.latent_random_field_rf[[i]] <- model_inla_i_rf$marginals.random
  latent_random_field_rf[[i]] <- model_inla_i_rf$summary.random
  t_cpu_rf[i+1] <- Sys.time()
  print(paste("Model", i, "is done."))
}

(t_cpu_rf %>% diff(.) %>% sum(.)/60)
(t_cpu %>% diff(.) %>% sum(.)/60) + (t_cpu_rf %>% diff(.) %>% sum(.)/60)

# Full data model ----

# model_inla1 <- inla(formula = formula_inla, family = "gaussian", data = inla.stack.data(stk_inf1), control.predictor = list(A = inla.stack.A(stk_inf1)), control.compute = list(config = TRUE, waic = TRUE), verbose = FALSE)
# model_inla2 <- inla(formula = formula_inla, family = "gaussian", data = inla.stack.data(stk_inf2), control.predictor = list(A = inla.stack.A(stk_inf2)), control.compute = list(config = TRUE, waic = TRUE), verbose = FALSE)

# spde <- inla.spde2.pcmatern(
#   mesh = mesh, alpha = 2,
#   prior.range = c(mesh$loc[c(mesh$segm$bnd$idx[,1], mesh$segm$bnd$idx[1,1]),1:2] %>% dist(x = .) %>% as.vector(.) %>% max(.)/5, 0.5),
#   prior.sigma = c(1,0.5)
# )

rho0 <- mesh$loc[c(mesh$segm$bnd$idx[,1], mesh$segm$bnd$idx[1,1]),1:2] %>% dist(x = .) %>% as.vector(.) %>% max(.)/5
sigma0 <- 1
alpha <- 2; d <- 2; nu <- alpha - d/2
kappa0 <- log(8*nu)/2 - log(rho0)
tau0 <- 0.5*(lgamma(nu) - lgamma(nu + d/2) - d/2*log(4*pi)) - log(sigma0) - nu*kappa0
spde <- inla.spde2.matern(mesh = mesh, B.tau = cbind(tau0, nu, -1), 
                          B.kappa = cbind(kappa0, -1, 0), alpha = 2, constr = FALSE)

spde.idx <- inla.spde.make.index(name = "sp.idx", n.spde = spde$n.spde, n.group = ntime)
A_spt <- inla.spde.make.A(mesh = mesh, loc = DFsf_temp[which(DFsf_temp$month.id %in% 1:ntime),] %>% st_coordinates(.), group = DFsimplified_temp$month.id)

stk_total <- inla.stack(data = list(y = DFsimplified_temp$max_sst),
                        A = list(A_spt, 1),
                        effects = list(
                          spde.idx,
                          list(intercept = rep(1, times = nrow(DFsimplified_temp)))
                        ),
                        tag = "stk_inf_total")

formula_inla <- y ~ -1 + intercept + f(sp.idx, model = spde, group = sp.idx.group, control.group = list(model = "ar1"))

model_inlat <- inla(formula = formula_inla, family = "gaussian", data = inla.stack.data(stk_total), control.predictor = list(A = inla.stack.A(stk_total)), control.compute = list(config = TRUE, waic = TRUE), verbose = FALSE)

j <- 1; k <- 1
jj <- j + ntime/k_times*(k-1)
mean_field <- inla.spde.make.A(mesh = mesh, loc = DFsf_temp[DFsf_temp$month.id==1,] %>% st_coordinates(.)) %*% cbind(model_inlat$summary.random$sp.idx$mean[(mesh$n*(jj-1)+1):(mesh$n*jj)], latent_random_field_rf[[k]]$sp.idx$mean[(mesh$n*(j-1)+1):(mesh$n*j)])
sd_field <- inla.spde.make.A(mesh = mesh, loc = DFsf_temp[DFsf_temp$month.id==1,] %>% st_coordinates(.)) %*% cbind(model_inlat$summary.random$sp.idx$sd[(mesh$n*(jj-1)+1):(mesh$n*jj)], latent_random_field_rf[[k]]$sp.idx$sd[(mesh$n*(j-1)+1):(mesh$n*j)])

gg_meanlf <- ggplot() + 
  # geom_tile(data = data.frame(DFsf_temp[DFsf_temp$month.id==1,] %>% st_coordinates(.), z = DFsimplified_temp$max_sst[DFsimplified_temp$month.id==1], group = "Data"), mapping = aes(x = X, y = Y, fill = z)) +
  # scale_fill_viridis_c(option = "turbo") + theme_minimal()
  geom_tile(data = data.frame(DFsf_temp[DFsf_temp$month.id==1,] %>% st_coordinates(.), z = mean_field[,1] + model_inlat$summary.fixed$mean, group = "Full model"), mapping = aes(x = X, y = Y, fill = z)) +
  geom_tile(data = data.frame(DFsf_temp[DFsf_temp$month.id==1,] %>% st_coordinates(.), z = mean_field[,2] + model_inla_i$summary.fixed$mean + 0*(mean.intercept[k] - last(mean.intercept)), group = "Sequential consensus"), mapping = aes(x = X, y = Y, fill = z)) +
  scale_fill_viridis_c(option = "turbo") + 
  facet_wrap(facets = ~group) + theme_minimal() +
  labs(title = "(A) Mean spatial latent field", fill = "Values") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

gg_stdevlf <- ggplot() + 
  geom_tile(data = data.frame(DFsf_temp[DFsf_temp$month.id==1,] %>% st_coordinates(.), z = sd_field[,1], group = "Full model"), mapping = aes(x = X, y = Y, fill = z)) +
  geom_tile(data = data.frame(DFsf_temp[DFsf_temp$month.id==1,] %>% st_coordinates(.), z = sd_field[,2], group = "Sequential consensus"), mapping = aes(x = X, y = Y, fill = z)) +
  scale_fill_viridis_c(option = "turbo") + 
  labs(title = "(B) Stdev. spatial latent field", fill = "Values") +
  facet_wrap(facets = ~group) + theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
grid.arrange(arrangeGrob(grobs = list(gg_meanlf, gg_stdevlf), ncol = 2))

j <- 1; k <- 1
jj <- j + mesh$n*ntime/k_times*(k-1)
ggplot() +
  geom_line(data = data.frame(model_inlat$marginals.random$sp.idx[[1+mesh$n*ntime/k_times*(1-1)]], group = "Full model", k = "Group 1"), mapping = aes(x = x, y = y, color = group)) +
  geom_line(data = data.frame(marginal.latent_random_field[[1]]$sp.idx[[1]], group = "Sequential consensus (A1)", k = "Group 1"), mapping = aes(x = x, y = y, color = group)) +
  geom_line(data = data.frame(marginal.latent_random_field_rf[[1]]$sp.idx[[1]], group = "Sequential consensus (A2)", k = "Group 1"), mapping = aes(x = x, y = y, color = group)) +
  geom_line(data = data.frame(model_inlat$marginals.random$sp.idx[[1+mesh$n*ntime/k_times*(2-1)]], group = "Full model", k = "Group 2"), mapping = aes(x = x, y = y, color = group)) +
  geom_line(data = data.frame(marginal.latent_random_field[[2]]$sp.idx[[1]], group = "Sequential consensus (A1)", k = "Group 2"), mapping = aes(x = x, y = y, color = group)) +
  geom_line(data = data.frame(marginal.latent_random_field_rf[[2]]$sp.idx[[1]], group = "Sequential consensus (A2)", k = "Group 2"), mapping = aes(x = x, y = y, color = group)) +
  geom_line(data = data.frame(model_inlat$marginals.random$sp.idx[[1+mesh$n*ntime/k_times*(3-1)]], group = "Full model", k = "Group 3"), mapping = aes(x = x, y = y, color = group)) +
  geom_line(data = data.frame(marginal.latent_random_field[[3]]$sp.idx[[1]], group = "Sequential consensus (A1)", k = "Group 3"), mapping = aes(x = x, y = y, color = group)) +
  geom_line(data = data.frame(marginal.latent_random_field_rf[[3]]$sp.idx[[1]], group = "Sequential consensus (A2)", k = "Group 3"), mapping = aes(x = x, y = y, color = group)) +
  geom_line(data = data.frame(model_inlat$marginals.random$sp.idx[[1+mesh$n*ntime/k_times*(4-1)]], group = "Full model", k = "Group 4"), mapping = aes(x = x, y = y, color = group)) +
  geom_line(data = data.frame(marginal.latent_random_field[[4]]$sp.idx[[1]], group = "Sequential consensus (A1)", k = "Group 4"), mapping = aes(x = x, y = y, color = group)) +
  geom_line(data = data.frame(marginal.latent_random_field_rf[[4]]$sp.idx[[1]], group = "Sequential consensus (A2)", k = "Group 4"), mapping = aes(x = x, y = y, color = group)) +
  geom_line(data = data.frame(model_inlat$marginals.random$sp.idx[[1+mesh$n*ntime/k_times*(5-1)]], group = "Full model", k = "Group 5"), mapping = aes(x = x, y = y, color = group)) +
  geom_line(data = data.frame(marginal.latent_random_field[[5]]$sp.idx[[1]], group = "Sequential consensus (A1)", k = "Group 5"), mapping = aes(x = x, y = y, color = group)) +
  geom_line(data = data.frame(marginal.latent_random_field_rf[[5]]$sp.idx[[1]], group = "Sequential consensus (A2)", k = "Group 5"), mapping = aes(x = x, y = y, color = group)) +
  geom_line(data = data.frame(model_inlat$marginals.random$sp.idx[[1+mesh$n*ntime/k_times*(6-1)]], group = "Full model", k = "Group 6"), mapping = aes(x = x, y = y, color = group)) +
  geom_line(data = data.frame(marginal.latent_random_field[[6]]$sp.idx[[1]], group = "Sequential consensus (A1)", k = "Group 6"), mapping = aes(x = x, y = y, color = group)) +
  geom_line(data = data.frame(marginal.latent_random_field_rf[[6]]$sp.idx[[1]], group = "Sequential consensus (A2)", k = "Group 6"), mapping = aes(x = x, y = y, color = group)) +
  labs(color = "Approach") + theme_minimal() + facet_wrap(facets = ~ k, scales = "free")

j <- 1; k <- 1
jj <- j + mesh$n*ntime/k_times*(k-1)
ggplot() + # model_inlat$marginals.random$sp.idx[[(mesh$n*(jj-1)+1):(mesh$n*jj)]] # marginal.latent_random_field[[k]]$sp.idx[[(mesh$n*(j-1)+1):(mesh$n*j)]]
  geom_line(data = data.frame(model_inlat$marginals.random$sp.idx[[jj]], group = "Full model"), mapping = aes(x = x, y = y, color = group)) +
  geom_line(data = data.frame(marginal.latent_random_field[[k]]$sp.idx[[j]], group = "Sequential consensus (A1)"), mapping = aes(x = x, y = y, color = group)) +
  geom_line(data = data.frame(marginal.latent_random_field_rf[[k]]$sp.idx[[j]], group = "Sequential consensus (A2)"), mapping = aes(x = x, y = y, color = group)) +
  labs(color = "Approach") + theme_minimal()

gg_intercept <- ggplot() + 
  geom_line(data = data.frame(model_inlat$marginals.fixed$intercept, group = "Full model"), mapping = aes(x = x, y = y, color = "Full model")) +
  geom_line(data = data.frame(model_inla_i$marginals.fixed$intercept, group = "Sequential consensus"), mapping = aes(x = x, y = y, color = "Sequential consensus")) +
  labs(title = "(A) \u03b2<sub>0</sub> posterior distribution", color = "Approach") + theme_minimal() +
  theme(plot.title = element_markdown(hjust = 0.5, face = "bold"), legend.position = "left")

gg_tau <- ggplot() + 
  geom_line(data = data.frame(model_inlat$internal.marginals.hyperpar$`Log precision for the Gaussian observations`, group = "Full model"), mapping = aes(x = x, y = y, color = "Full model")) +
  geom_line(data = data.frame(model_inla_i$internal.marginals.hyperpar$`Log precision for the Gaussian observations`, group = "Sequential consensus"), mapping = aes(x = x, y = y, color = "Sequential consensus")) +
  labs(title = "(B) log(\u03c4) posterior distribution", color = "Approach") + theme_minimal() + 
  theme(plot.title = element_markdown(hjust = 0.5, face = "bold"))  + guides(color = "none")

gg_rho <- ggplot() + 
  geom_line(data = data.frame(model_inlat$marginals.hyperpar$`GroupRho for sp.idx`, group = "Full model"), mapping = aes(x = x, y = y, color = "Full model")) +
  geom_line(data = data.frame(model_inla_i$marginals.hyperpar$`GroupRho for sp.idx`, group = "Sequential consensus"), mapping = aes(x = x, y = y, color = "Sequential consensus")) +
  labs(title = "(C) \u03c1 posterior distribution", color = "Approach") + theme_minimal() +
  theme(plot.title = element_markdown(hjust = 0.5, face = "bold")) + guides(color = "none")

gg_theta1 <- ggplot() + 
  geom_line(data = data.frame(model_inlat$marginals.hyperpar$`Theta1 for sp.idx`, group = "Full model"), mapping = aes(x = x, y = y, color = "Full model")) +
  geom_line(data = data.frame(model_inla_i$marginals.hyperpar$`Theta1 for sp.idx`, group = "Sequential consensus"), mapping = aes(x = x, y = y, color = "Sequential consensus")) +
  labs(title = "(D) \u03b8<sub>1</sub> posterior distribution", color = "Approach") + theme_minimal() +
  theme(plot.title = element_markdown(hjust = 0.5, face = "bold")) + guides(color = "none")

gg_theta2 <- ggplot() + 
  geom_line(data = data.frame(model_inlat$marginals.hyperpar$`Theta2 for sp.idx`, group = "Full model"), mapping = aes(x = x, y = y, color = "Full model")) +
  geom_line(data = data.frame(model_inla_i$marginals.hyperpar$`Theta2 for sp.idx`, group = "Sequential consensus"), mapping = aes(x = x, y = y, color = "Sequential consensus")) +
  labs(title = "(E) \u03b8<sub>2</sub> posterior distribution", color = "Approach") + theme_minimal() +
  theme(plot.title = element_markdown(hjust = 0.5, face = "bold")) + guides(color = "none")
grid.arrange(arrangeGrob(grobs = list(gg_intercept, gg_tau, gg_rho, gg_theta1, gg_theta2), layout_matrix = matrix(data = c(1,1,2,4,3,5), nrow = 2)))

l <- 4
ggplot() +
  geom_line(data = model_inla_i$internal.marginals.hyperpar[[l]] %>% inla.smarginal(marginal = .) %>% data.frame(.),
            mapping = aes(x = x, y = y), color = "red") +
  geom_line(data = 
              data.frame(
                x = seq(from = model_inla_i$internal.summary.hyperpar[l,1]-5*model_inla_i$internal.summary.hyperpar[l,2],
                        to = model_inla_i$internal.summary.hyperpar[l,1]+5*model_inla_i$internal.summary.hyperpar[l,2], length.out = 1E3),
                y = dnorm(x = seq(from = model_inla_i$internal.summary.hyperpar[l,1]-5*model_inla_i$internal.summary.hyperpar[l,2], 
                                  to = model_inla_i$internal.summary.hyperpar[l,1]+5*model_inla_i$internal.summary.hyperpar[l,2], length.out = 1E3),
                          mean = model_inla_i$internal.summary.hyperpar[l,1],
                          sd = model_inla_i$internal.summary.hyperpar[l,2])),
            mapping = aes(x = x, y = y), color = "blue"
  )


