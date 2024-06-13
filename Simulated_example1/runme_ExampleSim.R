
# Simulation and inference of an integrated model (stratified random sampling and preferential sampling)

## Libraries ----

library(INLA)
library(inlabru)
library(fmesher)
library(ggplot2)
library(gridExtra)
library(ggtext)
library(viridis)
library(dplyr)
library(sf)

## Custom functions ----

colsc <- function(...) { # Función para definir la escala de colores de los gráficos 
  scale_fill_gradientn(
    colours = turbo(n = 11),
    limits = range(..., na.rm = TRUE)
  )
}

RW2 <- function(N, x0 = 0, mu = 0, prec = 1, constr = TRUE, seed){
  # N: Number of nodes
  # x0: starting value
  # mu: shift in the normal error term
  # prec: precision value for the normal error term
  # constr: sum to zero constraint (by default is set to TRUE)
  # seed: seed for reproducibility of the results
  if(!missing(seed)){set.seed(seed)}
  x <- vector(mode="numeric", length=N)
  x[1] <- x0 + rnorm(1, mean=mu, sd=prec**(-1/2))
  x[2] <- 2*x[1]-x0+rnorm(1, mean=mu, sd=prec**(-1/2))
  for(i in 3:N){
    x[i] <- 2*x[i-1]+rnorm(1, mean=mu,sd=prec**(-1/2))-x[i-2]
  }
  if(constr){
    x <- x - mean(x)
  }
  return(x)
}

random_stratified_sampling <- function(seed.sample, length.xy=1, size.per.cell, n.cells.dim, ...){ # Completly random and stratified sampling design
  require(sp)
  require(sf)
  
  if(!missing(seed.sample)) set.seed(seed.sample)
  # Building the grid stratification for the sample process
  area_length <- length.xy; delta_length <- 0.001; numb_rec_by_axis <- n.cells.dim # Total number of rectangles n.cells.dim*n.cells.dim (numb_rec_by_axis). The delta_length is a delta extension to avoid Na's in the over operation between points and areas
  coord_agg_grid <- expand.grid(x=seq(-delta_length/2+area_length/numb_rec_by_axis/2, area_length+delta_length/2-area_length/numb_rec_by_axis/2, length.out=numb_rec_by_axis),
                                y=seq(-delta_length/2+area_length/numb_rec_by_axis/2, area_length+delta_length/2-area_length/numb_rec_by_axis/2, length.out=numb_rec_by_axis))
  spcoord <- sp::SpatialPoints(coord_agg_grid)
  sppixels <- sp::SpatialPixels(spcoord)
  PolygonsPixels <- as(sppixels, Class="SpatialPolygons")
  GridPolygons <- sf::st_as_sf(PolygonsPixels)
  
  # Locations simulated from each polygon of the grid 
  xy_sample <- do.call(rbind, lapply(X = 1:nrow(GridPolygons), FUN = function(x){
    sf::st_coordinates(sf::st_sample(x = GridPolygons[x,], size = size.per.cell, type = "random"))
  })
  )
  
  # Simulation of the response variable in the above defined locations
  # sample <- simulation_function(sim_coord = xy_sample, ...)
  return(list(xy_sample = xy_sample, GridPolygons = GridPolygons))
}

pref_sampling <- function(df, epoints, xstep, ystep){ # Preferential sampling design
  # df: it must be a gridded simulated data set, 
  # with the first column being the x-coordinate, the second the y-coordinate, 
  # it must contain a column called linpred, that is, the linear predictor of the model, 
  # and also have an id column for different temporal realizations
  # epoints: it is the expect number of points for the point process 
  # xstep: it is the the length of the step in the x-coord for the gridded simulated data (df)
  # ystep: it is the the length of the step in the y-coord for the gridded simulated data (df)
  require(spatstat)
  require(spatstat.random)
  require(spatstat.geom)
  
  a_optim <- log(epoints) - log(sum(xstep*ystep*exp(df$linpred)))
  
  list_lGCP <- lapply(X = 1:length(unique(df$id)), FUN = function(i){
    linpred_lgcp <- im(mat = exp(df$linpred[df$id==i] + a_optim),
                       xcol = seq(from = min(df[df$id==i,1]), to = max(df[df$id==i,1]), length.out = length(unique(df[df$id==i,1]))), 
                       yrow = seq(from = min(df[df$id==i,2]), to = max(df[df$id==i,2]), length.out = length(unique(df[df$id==i,2]))))
    xy_lGCP <- rpoispp(lambda = linpred_lgcp)
    return(xy_lGCP)})
  return(list_lGCP)
}

sfPolygonsBoundary <- function(mesh){
  BoundaryIndx <- sapply(1:(nrow(mesh$segm$bnd$idx)-1), function(i){
    return(mesh$segm$bnd$idx[i+1,1]==mesh$segm$bnd$idx[i,2])
  })
  BoundaryIndxLims <- which(!BoundaryIndx)
  BoundaryIndxLims <- c(0, BoundaryIndxLims, nrow(mesh$segm$bnd$idx))
  SpatialPolygonsBoundary <- st_polygon(sapply(2:(length(BoundaryIndxLims)), function(i){
    VertixBoundaryIndex <- c(mesh$segm$bnd$idx[(BoundaryIndxLims[i-1]+1),1], 
                             mesh$segm$bnd$idx[(BoundaryIndxLims[i-1]+1):BoundaryIndxLims[i],2])
    return(list(cbind(mesh$loc[VertixBoundaryIndex,1:2])))
  }))
  return(SpatialPolygonsBoundary)
}

sf_Triangulation <- function(mesh, crs_triangulation = NULL){
  polygon_Triangles <- lapply(1:nrow(mesh$graph$tv), function(i){
    return(st_polygon(x=list(as.matrix(mesh$loc[c(mesh$graph$tv[i,], mesh$graph$tv[i,1]),1:2]))))
  })
  sf_Triangulation <- do.call(c, polygon_Triangles)
  if(is.null(crs_triangulation)){
    coord <- try(st_sfc(sf_Triangulation, crs=st_crs(mesh)))
    if(class(coord)==class(try())){
      coord <- NA
    }
  } else{
    coord <- st_sfc(sf_Triangulation, crs=crs_triangulation)
  }
  sf_Triangulation <- coord
  return(sf_Triangulation)
}


mesh.dual <- function(mesh){
  if (mesh$manifold=='R2'){
    ce <- t(sapply(1:nrow(mesh$graph$tv), function(i)
      colMeans(mesh$loc[mesh$graph$tv[i, ], 1:2])))
    library(parallel)
    pls <- mclapply(1:mesh$n, function(i) {
      p <- unique(Reduce('rbind', lapply(1:3, function(k) {
        j <- which(mesh$graph$tv[,k]==i)
        if (length(j)>0)
          return(rbind(ce[j, , drop=FALSE],
                       cbind(mesh$loc[mesh$graph$tv[j, k], 1] +
                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 1],
                             mesh$loc[mesh$graph$tv[j, k], 2] +
                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 2])/2))
        else return(ce[j, , drop=FALSE])
      })))
      j1 <- which(mesh$segm$bnd$idx[,1]==i)
      j2 <- which(mesh$segm$bnd$idx[,2]==i)
      if ((length(j1)>0) | (length(j2)>0)) {
        p <- unique(rbind(mesh$loc[i, 1:2], p,
                          mesh$loc[mesh$segm$bnd$idx[j1, 1], 1:2]/2 +
                            mesh$loc[mesh$segm$bnd$idx[j1, 2], 1:2]/2,
                          mesh$loc[mesh$segm$bnd$idx[j2, 1], 1:2]/2 +
                            mesh$loc[mesh$segm$bnd$idx[j2, 2], 1:2]/2))
        yy <- p[,2]-mean(p[,2])/2-mesh$loc[i, 2]/2
        xx <- p[,1]-mean(p[,1])/2-mesh$loc[i, 1]/2
      }
      else {
        yy <- p[,2]-mesh$loc[i, 2]
        xx <- p[,1]-mesh$loc[i, 1]
      }
      p <- p[order(atan2(yy,xx)),]
      return(rbind(p,p[1,]))
    })
    return(
      lapply(1:mesh$n, function(i)
        st_polygon(x=list(pls[[i]]))
      )
    )
  }
  else stop("It only works for R2!")
}

marginal_consensus <- function(x, quantiles = c(0.025, 0.5, 0.975), constr = TRUE){
  # x: a list of the different summaries of the marginal distributions for the component over which perform the consensus
  # quantiles: the quantiles to compute in the posterior summary
  mat_prec <- lapply(X = x, FUN = function(i){i$sd**-2}) %>% unlist(.) %>% matrix(., ncol = length(x))
  prec_consensus <- apply(X = mat_prec, MARGIN = 1, FUN = sum) 
  mat_mean <- lapply(X = x, FUN = function(i){i$mean}) %>% unlist(.) %>% matrix(., ncol = length(x))
  mean_consensus <- sapply(X = 1:nrow(mat_mean), FUN = function(i){prec_consensus[i]**(-1)*sum(mat_prec[i,]*mat_mean[i,])})
  if(constr){mean_consensus <- mean_consensus - mean(mean_consensus)}
  
  if(length(quantiles)>0){
    vec_rename <- c(paste0("X", 1:length(quantiles))); names(vec_rename) <- c(paste0("q", quantiles))
    df_quantiles <- sapply(X = quantiles, FUN = function(q){qnorm(p = q, mean = mean_consensus, sd = prec_consensus**(-1/2))}) %>% data.frame(.) %>% 
      rename(., all_of(vec_rename))
    summary_consensus_marginals <- data.frame(ID = x[[1]]$ID, mean = mean_consensus, sd = prec_consensus**(-1/2), df_quantiles)
  } else{
    summary_consensus_marginals <- data.frame(ID = x[[1]]$ID, mean = mean_consensus, sd = prec_consensus**(-1/2))
  }
  return(summary_consensus_marginals)
}

fix.Q <- function(Q){
  d <- diag(Q)
  Q <- Q + t(Q)
  diag(Q) <- d
  return (Q)
}

## Defining the spatio-temporal geostatistical model for the abundance ----

seed <- 12345
set.seed(seed)

beta0 <- -1 # Intercept
beta1 <- 2 # Beta coefficient
prec <- 15 # Precision Gamma
prec_rw2 <- 15 # Precision random walk of second order (rw2)

# Building the mesh to simulate the spatial effect and a covariate with a smooth spatial structure
mesh2d_sim <- fm_mesh_2d_inla(loc.domain = matrix(data = c(0,0,10,0,10,10,0,10,0,0), ncol = 2), max.edge = c(0.25, 1), offset = c(0, -0.1))

u_nodes <- fm_matern_sample(x = mesh2d_sim, n = 1, rho = 3, sigma = 1) # Spatial effect
u_nodes <- u_nodes - mean(u_nodes) # Zero global mean constraint
cov_nodes <- fm_matern_sample(x = mesh2d_sim, n = 1, rho = 6, sigma = 0.25) %>% apply(X = ., MARGIN = 2, FUN = function(x){x-min(x)+0.2}) # Covariate with spatial structure
t_rw2 <- RW2(N = 10, prec = prec_rw2, seed = seed)
ggrw2 <- ggplot() +
  geom_line(data = data.frame(id = seq_along(t_rw2), rw2 = t_rw2), mapping = aes(x = id, y = rw2)) +
  labs(title = "(A) Temporal trend") + theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
# ggsave(filename = "./Example_simulated1/rw2_plot.png", plot = ggrw2, device = "png", width = 3000, height = 3000, units = "px", dpi = 600)

# Prediction over the xy_gridpred for the LGCP simulation
xy_gridpred <- expand.grid(x = seq(from = 10/2.5E2/2, to = 10-10/2.5E2/2, length.out = 2.5E2), y = seq(from = 10/2.5E2/2, to = 10-10/2.5E2/2, length.out = 2.5E2)) %>% as.matrix(.)
A_xy_pred <- fm_basis(x = mesh2d_sim, loc = xy_gridpred) # Projection from mesh2d_sim nodes to prediction grid locations
linpred_y_pred <- rep(x = beta0 + beta1*drop(A_xy_pred %*% cov_nodes) + drop(A_xy_pred %*% u_nodes), times = length(t_rw2)) + rep(t_rw2, each = nrow(A_xy_pred))

ggcov <- ggplot() +
  geom_tile(data = data.frame(x = xy_gridpred[,1], y = xy_gridpred[,2], cov = drop(A_xy_pred %*% cov_nodes)), mapping = aes(x = x, y = y, fill = cov)) +
  scale_fill_viridis_c(option = "turbo") + labs(fill = "Values") +
  # guides(fill = "none") + 
  theme_minimal() +
  labs(title = "(B) Covariate structure") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
# ggsave(filename = "./Example_simulated1/cov_plot.png", plot = ggcov, device = "png", width = 3000, height = 3000, units = "px", dpi = 600)

ggspde2 <- ggplot() +
  geom_tile(data = data.frame(x = xy_gridpred[,1], y = xy_gridpred[,2], spde2 = drop(A_xy_pred %*% u_nodes)), mapping = aes(x = x, y = y, fill = spde2)) +
  scale_fill_viridis_c(option = "turbo") + labs(fill = "Values") +
  # guides(fill = "none") + 
  theme_minimal() + 
  labs(title = "(C) Spatial effect") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
# ggsave(filename = "./Example_simulated1/spde2_plot.png", plot = ggspde2, device = "png", width = 3000, height = 3000, units = "px", dpi = 600)

beta_copy <- 1.5
linpred_lgcp_pred <- rep(x = beta0 + drop(A_xy_pred %*% u_nodes), times = length(t_rw2)) + rep(beta_copy*t_rw2, each = nrow(A_xy_pred))

alpha_pred <- exp(2*linpred_y_pred)*prec; beta_pred <- exp(linpred_y_pred)*prec
y_lgcp_pred <- rgamma(n = length(linpred_y_pred), shape = alpha_pred, rate = beta_pred)
# if(any(y_lgcp_pred<1E-8)){y_lgcp_pred[y_lgcp_pred<1E-8] <- 1E-8} # Simulated data can be 0 or too close to 0, we correct this by replacing pretty low values by 1E-8

dfpred_to_lgcp_Gamma <- st_sf(data.frame(id = rep(1:length(t_rw2), each = nrow(A_xy_pred)), y = y_lgcp_pred, linpred = linpred_y_pred), geometry = st_sfc(lapply(X=1:nrow(xy_gridpred), FUN = function(i){st_point(x = xy_gridpred[i,])})))
df_Gamma_lgcp <- data.frame(xy_gridpred, y_gamma = dfpred_to_lgcp_Gamma$y, linpred = linpred_lgcp_pred, id = dfpred_to_lgcp_Gamma$id) 
df_Gamma_lgcp %>% ggplot(.) + geom_tile(aes(x=x, y=y, fill=y_gamma)) + scale_fill_viridis_c(option = "turbo") + facet_wrap(facets = ~id, ncol = 2) + theme_minimal()

ggypred <- ggplot() +
  geom_tile(data = data.frame(x = df_Gamma_lgcp$x[df_Gamma_lgcp$id==1], y = df_Gamma_lgcp$y[df_Gamma_lgcp$id==1], ypred = df_Gamma_lgcp$y_gamma[df_Gamma_lgcp$id==1]), mapping = aes(x = x, y = y, fill = ypred)) +
  scale_fill_viridis_c(option = "turbo") + labs(fill = "Values") + 
  # guides(fill = "none") +
  theme_minimal() +
  labs(title = "(D) Abundance simulation") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
# ggsave(filename = "./Example_simulated1/ypred_plot.png", plot = ggypred, device = "png", width = 3000, height = 3000, units = "px", dpi = 600)

## Simulation of the two sampling designs ----

stratified_sampling <- lapply(X = 1:length(t_rw2), FUN = function(i){random_stratified_sampling(seed.sample = i, length.xy = 10, size.per.cell = 10, n.cells.dim = 5)})

i <- 1
ggstratified <- ggplot() +
  geom_tile(data = df_Gamma_lgcp[df_Gamma_lgcp$id==i,], mapping = aes(x=x, y=y, fill=y_gamma)) +
  geom_sf(data = stratified_sampling[[i]]$GridPolygons, mapping = aes(), alpha = 0.5, color = "black", linewidth = 1.25) +
  geom_point(data = data.frame(x = stratified_sampling[[i]]$xy_sample[,1], y = stratified_sampling[[i]]$xy_sample[,2]), mapping = aes(x=x, y=y)) +
  scale_fill_viridis_c(option = "turbo") + 
  guides(fill="none") +
  theme_minimal() +
  labs(title = "(E) Stratified random sampling") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
# ggsave(filename = "./Example_simulated1/stratified_plot.png", plot = ggstratified, width = 3000, height = 3000, units = "px", dpi = 600) 

A_list_stratified <- lapply(X = stratified_sampling, FUN = function(x){fm_basis(x = mesh2d_sim, loc = x$xy_sample)}) %>% do.call(what = rbind, .) # the projection matrix for the all the samples
linpred_rstratified_samples <- beta0 + beta1*drop(A_list_stratified %*% cov_nodes) + drop(A_list_stratified %*% u_nodes) + rep(t_rw2, times = lapply(X = stratified_sampling, FUN = function(x){x$xy_sample %>% nrow(.)}) %>% unlist(.))

alpha_pred <- exp(2*linpred_rstratified_samples)*prec; beta_pred <- exp(linpred_rstratified_samples)*prec
y_stratified_samples <- rgamma(n = length(linpred_rstratified_samples), shape = alpha_pred, rate = beta_pred)
# if(any(y_stratified_samples<1E-8)){y_stratified_samples[y_stratified_samples<1E-8] <- 1E-8} # Simulated data can be 0 or too close to 0, we correct this by replacing pretty low values by 1E-8

df_ystratified_samples <- data.frame(x = lapply(X = stratified_sampling, FUN = function(x){x$xy_sample[,1]}) %>% unlist(.), 
                                     y = lapply(X = stratified_sampling, FUN = function(x){x$xy_sample[,2]}) %>% unlist(.), 
                                     y_gamma = y_stratified_samples, 
                                     id = rep(x = 1:length(t_rw2), times = lapply(X = stratified_sampling, FUN = function(x){nrow(x$xy_sample)}) %>% unlist(.)))
ggplot() +
  geom_tile(data = df_Gamma_lgcp, mapping = aes(x=x, y=y, fill=y_gamma), alpha = 0.85) +
  geom_point(data = df_ystratified_samples, mapping = aes(x=x, y=y, fill = y_gamma), color = "black", shape=21, size=3, stroke = 1) +
  scale_color_viridis_c(option = "turbo") + scale_fill_viridis_c(option = "turbo") +
  facet_wrap(facets = ~id, ncol = 2) +
  theme_minimal()

# Prediction over the xy_gridpred for the LGCP simulation

xy_sim_lgcp <- pref_sampling(df=df_Gamma_lgcp, epoints = 25*10*10, xstep = 10/2.5E2, ystep = 10/2.5E2)
lapply(X = xy_sim_lgcp, FUN = function(i){i$n}) %>% do.call(what = sum, .) # the value is cole to the expoected umber of points (epoints)
a_optim <- log(2500) - log(sum(0.04*0.04*exp(df_Gamma_lgcp$linpred)))

A_list_lgcp <- lapply(X = xy_sim_lgcp, FUN = function(x){fm_basis(x = mesh2d_sim, loc = cbind(x$y,x$x))}) %>% do.call(what = rbind, .)
linpred_lgcp_samples <- beta0 + beta1*drop(A_list_lgcp %*% cov_nodes) + drop(A_list_lgcp %*% u_nodes)  + rep(t_rw2, times = lapply(X = xy_sim_lgcp, FUN = function(x){x$n}) %>% unlist(.))

alpha_pred <- exp(2*linpred_lgcp_samples)*prec; beta_pred <- exp(linpred_lgcp_samples)*prec
y_lgcp_samples <- rgamma(n = length(linpred_lgcp_samples), shape = alpha_pred, rate = beta_pred)
# if(any(y_lgcp_samples<1E-8)){y_lgcp_samples[y_lgcp_samples<1E-8] <- 1E-8} # Simulated data can be 0 or too close to 0, we correct this by replacing pretty low values by 1E-8

df_ylgcp_samples <- data.frame(x = lapply(X = xy_sim_lgcp, FUN = function(x){x$y}) %>% unlist(.), 
                               y = lapply(X = xy_sim_lgcp, FUN = function(x){x$x}) %>% unlist(.), 
                               y_gamma = y_lgcp_samples, 
                               id = rep(x = 1:length(t_rw2), times = lapply(X = xy_sim_lgcp, FUN = function(x){x$n}) %>% unlist(.)))

i <- 1
gglgcp <- ggplot() +
  geom_tile(data = df_Gamma_lgcp[df_Gamma_lgcp$id==i,], mapping = aes(x=x, y=y, fill=y_gamma), alpha = 0.85) +
  geom_point(data = df_ylgcp_samples[df_ylgcp_samples$id==i,], mapping = aes(x=x, y=y, fill = y_gamma), color = "black", shape=21, size=3, stroke = 1) +
  scale_color_viridis_c(option = "turbo") + scale_fill_viridis_c(option = "turbo") + labs(fill = "Values") +
  # guides(fill = "none") + 
  theme_minimal() +
  labs(title = "(F) Preferential sampling") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
# ggsave(filename = "./Example_simulated1/lgcp_plot.png", plot = gglgcp, width = 3000, height = 3000, units = "px", dpi = 600) 

ggplot() +
  geom_tile(data = df_Gamma_lgcp, mapping = aes(x=x, y=y, fill=y_gamma), alpha = 0.85) +
  geom_point(data = df_ylgcp_samples, mapping = aes(x=x, y=y, fill = y_gamma), color = "black", shape=21, size=3, stroke = 1) +
  scale_color_viridis_c(option = "turbo") + scale_fill_viridis_c(option = "turbo") +
  facet_wrap(facets = ~id, ncol = 2) +
  theme_minimal()

# csc <- colsc(df_Gamma_lgcp$y_gamma[df_Gamma_lgcp$id==1], df_ylgcp_samples[df_ylgcp_samples$id==1,], df_Gamma_lgcp$cov, drop(A_xy_pred %*% u_nodes))
ggsim <- grid.arrange(arrangeGrob(grobs = list(ggrw2, ggcov, ggspde2, ggypred, ggstratified, gglgcp), ncol = 3))
# ggsave(filename = "./Simulated_example1/Figures/simulation_plots.png", plot = ggsim, device = "png", width = 9000, height = 4500, units = "px", dpi = 600)

# Integrated model inference ----

stk_y_stratified <- inla.stack(data = list(y = cbind(df_ystratified_samples$y_gamma, NA)),
                               A = list(A_list_stratified, 1),
                               effects = list(
                                 list(sp.idx_gamma = 1:mesh2d_sim$n),
                                 list(
                                   intercept_gamma = rep(1, nrow(df_ystratified_samples)),
                                   cov_gamma = drop(A_list_stratified %*% cov_nodes),
                                   time_gamma = df_ystratified_samples$id
                                 )
                               ),
                               tag = "stk_y_stratified")

stk_y_lgcp <- inla.stack(data = list(y = cbind(df_ylgcp_samples$y_gamma, NA)),
                         A = list(A_list_lgcp, 1),
                         effects = list(
                           list(sp.idx_gamma = 1:mesh2d_sim$n),
                           list(
                             intercept_gamma = rep(1, nrow(df_ylgcp_samples)),
                             cov_gamma = drop(A_list_lgcp %*% cov_nodes),
                             time_gamma = df_ylgcp_samples$id
                           )
                         ),
                         tag = "stk_y_lgcp")

ldomain <- unique(mesh2d_sim$loc[mesh2d_sim$segm$int$idx,1:2])
ldomain <- rbind(ldomain, ldomain[1,])
dmesh <- mesh.dual(mesh = mesh2d_sim)
domainPolygon <- st_polygon(x=list(ldomain))
w <- sapply(1:length(dmesh), function(i) {
  if (sf::st_intersects(dmesh[[i]], domainPolygon, sparse=FALSE))
    return(sf::st_area(sf::st_intersection(dmesh[[i]], domainPolygon)))
  else return(0)
})

A_diag_lgcp <- rep(list(Diagonal(n = mesh2d_sim$n, x = 1)), times = length(unique(df_ylgcp_samples$id))) %>% do.call(what = rbind, .)
stk_lgcp <- inla.stack(data = list(y = cbind(NA, c(rep(0,mesh2d_sim$n*length(t_rw2)), rep(1,nrow(df_ylgcp_samples)))), 
                                   e = c(rep(w, times = length(unique(df_ylgcp_samples$id))), rep(0, times = nrow(df_ylgcp_samples)))),
                       A = list(rbind(A_diag_lgcp, A_list_lgcp), 1),
                       effects = list(
                         list(sp.idx_lgcp = 1:mesh2d_sim$n),
                         list(
                           intercept_lgcp = rep(1, times = A_diag_lgcp %>% nrow(.) + nrow(df_ylgcp_samples)),
                           # cov = c(drop(A_diag_lgcp %*% cov_nodes), drop(A_list_lgcp %*% cov_nodes)),
                           time_lgcp = c(rep(unique(df_ylgcp_samples$id), each = mesh2d_sim$n), df_ylgcp_samples$id)
                         )
                       ),
                       tag = "stk_lgcp")

stk_total <- inla.stack(stk_y_stratified, stk_y_lgcp, stk_lgcp)

rho0 <- mesh2d_sim$loc[c(mesh2d_sim$segm$bnd$idx[,1], mesh2d_sim$segm$bnd$idx[1,1]),1:2] %>% dist(x = .) %>% as.vector(.) %>% max(.)/5
sigma0 <- 1
alpha <- 2; d <- 2; nu <- alpha - d/2
kappa0 <- log(8*nu)/2 - log(rho0)
tau0 <- 0.5*(lgamma(nu) - lgamma(nu + d/2) - d/2*log(4*pi)) - log(sigma0) - nu*kappa0
spde <- inla.spde2.matern(mesh = mesh2d_sim, B.tau = cbind(tau0, nu, -1), 
                          B.kappa = cbind(kappa0, -1, 0), alpha = 2, constr = TRUE)

formula_integrated <- y ~ -1 + intercept_gamma + cov_gamma + f(time_gamma, model = "rw2", constr = TRUE) + f(sp.idx_gamma, model = spde) +
  intercept_lgcp + f(time_lgcp, model = "rw2", constr = TRUE) + f(sp.idx_lgcp, copy = "sp.idx_gamma", fixed = TRUE)

integrated_model <- inla(data = inla.stack.data(stk_total), E = inla.stack.data(stk_total)$e,
                         family = c("gamma", "poisson"), formula = formula_integrated,
                         control.predictor = list(A = inla.stack.A(stk_total)),
                         control.compute = list(config = TRUE, waic = TRUE),
                         verbose = FALSE)

ggplot() + 
  geom_line(data = data.frame(id = seq_along(t_rw2), rw2 = t_rw2), mapping = aes(x = id, y = rw2), color = "red") +
  geom_ribbon(data = integrated_model$summary.random$time_gamma %>% rename(., all_of(c(q1='0.025quant', q3='0.975quant'))),
              mapping = aes(x = ID, ymin = q1, ymax = q3), alpha=0.4) +
  geom_line(data = integrated_model$summary.random$time_gamma, mapping = aes(x = ID, y = mean)) +
  theme_minimal()

ggplot() + 
  geom_line(data = data.frame(id = seq_along(t_rw2), rw2 = beta_copy*t_rw2), mapping = aes(x = id, y = rw2), color = "red") +
  geom_ribbon(data = integrated_model$summary.random$time_lgcp %>% rename(., all_of(c(q1='0.025quant', q3='0.975quant'))),
              mapping = aes(x = ID, ymin = q1, ymax = q3), alpha=0.4) +
  geom_line(data = integrated_model$summary.random$time_lgcp, mapping = aes(x = ID, y = mean)) +
  theme_minimal()

csc_sp <- colsc(drop(A_xy_pred %*% u_nodes), drop(A_xy_pred %*% integrated_model$summary.random$sp.idx_gamma$mean))
ggspde2_sim <- ggplot() +
  geom_tile(data = data.frame(x = xy_gridpred[,1], y = xy_gridpred[,2], spde2 = drop(A_xy_pred %*% u_nodes)), mapping = aes(x = x, y = y, fill = spde2)) +
  labs(fill = "Values") + theme_minimal() + 
  labs(title = "(B) Simulated spatial effect") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsp_im <- 
  ggplot() + 
  geom_tile(data = data.frame(x = xy_gridpred[,1], y = xy_gridpred[,2], sp = drop(A_xy_pred %*% integrated_model$summary.random$sp.idx_gamma$mean)), 
            mapping = aes(x = x, y = y, fill = sp)) + labs(fill = "Values") +
  theme_minimal() + labs(title = "(A) Infered spatial effect") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))

grid.arrange(arrangeGrob(grobs = list(ggsp_im + csc_sp, ggspde2_sim + csc_sp), ncol = 2))

integrated_model$summary.fixed
integrated_model$summary.hyperpar

# Sequential consensus inference ----

stk_y_stratified_sc <- inla.stack(data = list(y = cbind(df_ystratified_samples$y_gamma)),
                                  A = list(A_list_stratified, 1),
                                  effects = list(
                                    list(sp.idx_gamma = 1:mesh2d_sim$n),
                                    list(
                                      intercept_gamma = rep(1, nrow(df_ystratified_samples)),
                                      cov_gamma = drop(A_list_stratified %*% cov_nodes),
                                      time_gamma = df_ystratified_samples$id
                                    )
                                  ),
                                  tag = "stk_y_stratified")

stk_y_lgcp_sc <- inla.stack(data = list(y = cbind(df_ylgcp_samples$y_gamma)),
                            A = list(A_list_lgcp, 1),
                            effects = list(
                              list(sp.idx_gamma = 1:mesh2d_sim$n),
                              list(
                                intercept_gamma = rep(1, nrow(df_ylgcp_samples)),
                                cov_gamma = drop(A_list_lgcp %*% cov_nodes),
                                time_gamma = df_ylgcp_samples$id
                              )
                            ),
                            tag = "stk_y_lgcp")

## Sequential consensus inference step 1 (stratified random sampling) ----

hyper.pc_rw2 <- list(prec = list(prior = "pc.prec", param = c(1,0.5)))
formula_integrated_sc_y1gamma <- y ~ -1 + intercept_gamma + cov_gamma + f(time_gamma, model = "rw2", constr = TRUE) + f(sp.idx_gamma, model = spde)

t1 <- Sys.time()
sc_y1_model <- inla(data = inla.stack.data(stk_y_stratified_sc),
                    family = c("gamma"), formula = formula_integrated_sc_y1gamma,
                    control.predictor = list(A = inla.stack.A(stk_y_stratified_sc)),
                    control.compute = list(config = TRUE, waic = TRUE),
                    verbose = FALSE)

fixed.mean <- sc_y1_model$summary.fixed$mean %>% lapply(X = ., FUN=function(x){x}); fixed.prec <- sc_y1_model$summary.fixed$sd**-2 %>% lapply(X = ., FUN=function(x){x})
names(fixed.mean) <- names(fixed.prec) <- rownames(sc_y1_model$summary.fixed)

norm_internal.gamma.prec <- list(prior = "normal", param = c(sc_y1_model$internal.summary.hyperpar[1,1], sc_y1_model$internal.summary.hyperpar[1,2]**-2), fixed = FALSE)
norm_internal.rw2 <- list(prior = "normal", param = c(sc_y1_model$internal.summary.hyperpar[2,1], sc_y1_model$internal.summary.hyperpar[2,2]**-2), fixed = FALSE)

spde_sc_y2gammma <- inla.spde2.matern(mesh = mesh2d_sim, B.tau = cbind(tau0, nu, -1), 
                          B.kappa = cbind(kappa0, -1, 0), alpha = 2, constr = TRUE,
                          theta.prior.mean = sc_y1_model$internal.summary.hyperpar[3:4,1], 
                          theta.prior.prec = sc_y1_model$internal.summary.hyperpar[3:4,2]**-2)

formula_integrated_sc_y2gamma <- y ~ -1 + intercept_gamma + cov_gamma + f(time_gamma, model = "rw2", constr = TRUE, hyper = list(theta = norm_internal.rw2)) + f(sp.idx_gamma, model = spde_sc_y2gammma)

sc_y2_model <- inla(data = inla.stack.data(stk_y_lgcp_sc),
                    family = c("gamma"), formula = formula_integrated_sc_y2gamma,
                    control.predictor = list(A = inla.stack.A(stk_y_lgcp_sc)),
                    control.fixed = list(mean = fixed.mean, prec = fixed.prec),
                    control.family = list(hyper = list(theta = list(norm_internal.gamma.prec))),
                    control.mode = list(theta = sc_y1_model$internal.summary.hyperpar$mode),
                    control.compute = list(config = TRUE, waic = TRUE),
                    verbose = FALSE)

stk_lgcp_sc <- inla.stack(data = list(y = cbind(c(rep(0,mesh2d_sim$n*length(t_rw2)), rep(1,nrow(df_ylgcp_samples)))), 
                                      e = c(rep(w, times = length(unique(df_ylgcp_samples$id))), rep(0, times = nrow(df_ylgcp_samples)))),
                          A = list(rbind(A_diag_lgcp, A_list_lgcp), 1),
                          effects = list(
                            list(sp.idx_lgcp = 1:mesh2d_sim$n),
                            list(
                              intercept_lgcp = rep(1, times = A_diag_lgcp %>% nrow(.) + nrow(df_ylgcp_samples)),
                              # cov = c(drop(A_diag_lgcp %*% cov_nodes), drop(A_list_lgcp %*% cov_nodes)),
                              time_lgcp = c(rep(unique(df_ylgcp_samples$id), each = mesh2d_sim$n), df_ylgcp_samples$id)
                            )
                          ),
                          tag = "stk_lgcp")


spde_sc_lgcp <- inla.spde2.matern(mesh = mesh2d_sim, B.tau = cbind(tau0, nu, -1), 
                                      B.kappa = cbind(kappa0, -1, 0), alpha = 2, constr = TRUE,
                                      theta.prior.mean = c(sc_y2_model$internal.summary.hyperpar[3:4,1]), 
                                      theta.prior.prec = c(sc_y2_model$internal.summary.hyperpar[3:4,2]**-2))

formula_integrated_sc_lgcp <- y ~ -1 + intercept_lgcp + f(time_lgcp, model = "rw2", constr = TRUE) + f(sp.idx_lgcp, model = spde_sc_lgcp)

sclgcp_model <- inla(data = inla.stack.data(stk_lgcp_sc), E = inla.stack.data(stk_lgcp_sc)$e,
                    family = c("poisson"), formula = formula_integrated_sc_lgcp,
                    control.predictor = list(A = inla.stack.A(stk_lgcp_sc)),
                    control.mode = list(theta = c(NA, sc_y2_model$internal.summary.hyperpar$mode[3:4])),
                    control.compute = list(config = TRUE, waic = TRUE),
                    verbose = FALSE)

t2 <- Sys.time()
difftime(t2,t1, units = "mins")
integrated_model$cpu.used/60

### Graphical results after sc y_stratified modelling ----

ggplot() + 
  geom_line(data = integrated_model$marginals.hyperpar$`Precision parameter for the Gamma observations`, mapping = aes(x=x, y=y), color = "red") +
  geom_line(data = sc_y1_model$marginals.hyperpar$`Precision parameter for the Gamma observations`, mapping = aes(x=x, y=y), color = "blue") +
  geom_vline(xintercept = prec) +
  theme_minimal()

ggplot() +
  geom_line(data = integrated_model$marginals.hyperpar$`Precision for time_gamma`, mapping = aes(x=x, y=y), color = "red") +
  geom_line(data = sc_y1_model$marginals.hyperpar$`Precision for time_gamma`, mapping = aes(x=x, y=y), color = "blue") +
  geom_vline(xintercept = prec_rw2) +
  theme_minimal()

ggplot() +
  geom_line(data = integrated_model$marginals.hyperpar$`Theta1 for sp.idx_gamma` %>% inla.tmarginal(fun = function(x){rho0*exp(x)}, marginal = .), mapping = aes(x=x, y=y), color = "red") +
  geom_line(data = sc_y1_model$marginals.hyperpar$`Theta1 for sp.idx_gamma` %>% inla.tmarginal(fun = function(x){rho0*exp(x)}, marginal = .), mapping = aes(x=x, y=y), color = "blue") +
  geom_vline(xintercept = 3) +
  theme_minimal()

ggplot() +
  geom_line(data = integrated_model$marginals.hyperpar$`Theta2 for sp.idx_gamma` %>% inla.tmarginal(fun = function(x){exp(x)}, marginal = .), mapping = aes(x=x, y=y), color = "red") +
  geom_line(data = sc_y1_model$marginals.hyperpar$`Theta2 for sp.idx_gamma` %>% inla.tmarginal(fun = function(x){exp(x)}, marginal = .), mapping = aes(x=x, y=y), color = "blue") +
  geom_vline(xintercept = 1) +
  theme_minimal()

ggplot() +
  geom_line(data = integrated_model$marginals.fixed$intercept_gamma, mapping = aes(x=x, y=y), color = "red") +
  geom_line(data = sc_y1_model$marginals.fixed$intercept_gamma, mapping = aes(x=x, y=y), color = "blue") +
  geom_vline(xintercept = -1) +
  theme_minimal()

ggplot() +
  geom_line(data = integrated_model$marginals.fixed$cov_gamma, mapping = aes(x=x, y=y), color = "red") +
  geom_line(data = sc_y1_model$marginals.fixed$cov_gamma, mapping = aes(x=x, y=y), color = "blue") +
  geom_vline(xintercept = 2) +
  theme_minimal()

### Graphical results after sc y_preferential modelling ----

ggplot() +
  geom_line(data = integrated_model$marginals.hyperpar$`Precision parameter for the Gamma observations`, mapping = aes(x=x, y=y), color = "red") +
  geom_line(data = sc_y2_model$marginals.hyperpar$`Precision parameter for the Gamma observations`, mapping = aes(x=x, y=y), color = "blue") +
  geom_vline(xintercept = prec) +
  theme_minimal()

ggplot() +
  geom_line(data = integrated_model$marginals.hyperpar$`Precision for time_gamma`, mapping = aes(x=x, y=y), color = "red") +
  geom_line(data = sc_y2_model$marginals.hyperpar$`Precision for time_gamma`, mapping = aes(x=x, y=y), color = "blue") +
  geom_vline(xintercept = prec_rw2) +
  theme_minimal()

ggplot() +
  geom_line(data = integrated_model$marginals.hyperpar$`Theta1 for sp.idx_gamma` %>% inla.tmarginal(fun = function(x){rho0*exp(x)}, marginal = .), mapping = aes(x=x, y=y), color = "red") +
  geom_line(data = sc_y2_model$marginals.hyperpar$`Theta1 for sp.idx_gamma` %>% inla.tmarginal(fun = function(x){rho0*exp(x)}, marginal = .), mapping = aes(x=x, y=y), color = "blue") +
  geom_vline(xintercept = 3) +
  theme_minimal()

ggplot() +
  geom_line(data = integrated_model$marginals.hyperpar$`Theta2 for sp.idx_gamma` %>% inla.tmarginal(fun = function(x){exp(x)}, marginal = .), mapping = aes(x=x, y=y), color = "red") +
  geom_line(data = sc_y2_model$marginals.hyperpar$`Theta2 for sp.idx_gamma` %>% inla.tmarginal(fun = function(x){exp(x)}, marginal = .), mapping = aes(x=x, y=y), color = "blue") +
  geom_vline(xintercept = 1) +
  theme_minimal()

ggplot() +
  geom_line(data = integrated_model$marginals.fixed$intercept_gamma, mapping = aes(x=x, y=y), color = "red") +
  geom_line(data = sc_y2_model$marginals.fixed$intercept_gamma, mapping = aes(x=x, y=y), color = "blue") +
  geom_vline(xintercept = -1) +
  theme_minimal()

ggplot() +
  geom_line(data = integrated_model$marginals.fixed$cov_gamma, mapping = aes(x=x, y=y), color = "red") +
  geom_line(data = sc_y2_model$marginals.fixed$cov_gamma, mapping = aes(x=x, y=y), color = "blue") +
  geom_vline(xintercept = 2) +
  theme_minimal()

### Graphical results after sc lgcp modelling ---- 

ggplot() +
  geom_line(data = integrated_model$marginals.hyperpar$`Precision for time_gamma`, mapping = aes(x=x, y=y), color = "red") +
  geom_line(data = sc_y2_model$marginals.hyperpar$`Precision for time_gamma`, mapping = aes(x=x, y=y), color = "blue") +
  geom_vline(xintercept = prec_rw2/1.5**2) +
  theme_minimal()

ggplot() +
  geom_line(data = integrated_model$marginals.hyperpar$`Theta1 for sp.idx_gamma` %>% inla.tmarginal(fun = function(x){rho0*exp(x)}, marginal = .), mapping = aes(x=x, y=y), color = "red") +
  geom_line(data = sclgcp_model$marginals.hyperpar$`Theta1 for sp.idx_lgcp` %>% inla.tmarginal(fun = function(x){rho0*exp(x)}, marginal = .), mapping = aes(x=x, y=y), color = "blue") +
  geom_vline(xintercept = 3) +
  theme_minimal()

ggplot() +
  geom_line(data = integrated_model$marginals.hyperpar$`Theta2 for sp.idx_gamma` %>% inla.tmarginal(fun = function(x){exp(x)}, marginal = .), mapping = aes(x=x, y=y), color = "red") +
  geom_line(data = sclgcp_model$marginals.hyperpar$`Theta2 for sp.idx_lgcp` %>% inla.tmarginal(fun = function(x){exp(x)}, marginal = .), mapping = aes(x=x, y=y), color = "blue") +
  geom_vline(xintercept = 1) +
  theme_minimal()

gg_intery <- ggplot() +
  geom_line(data = integrated_model$marginals.fixed$intercept_gamma, mapping = aes(x=x, y=y), color = "red") +
  geom_line(data = sc_y2_model$marginals.fixed$intercept_gamma, mapping = aes(x=x, y=y), color = "blue") +
  geom_vline(xintercept = -1) + theme_minimal() +
  labs(title = "(A) Posterior distribution \u03B2<sub>0</sub>") + theme(plot.title = element_markdown(hjust = 0.5, face = "bold"))

gg_covy <- ggplot() +
  geom_line(data = integrated_model$marginals.fixed$cov_gamma, mapping = aes(x=x, y=y), color = "red") +
  geom_line(data = sc_y2_model$marginals.fixed$cov_gamma, mapping = aes(x=x, y=y), color = "blue") +
  geom_vline(xintercept = 2) + theme_minimal() +
  labs(title = "(B) Posterior distribution \u03B2<sub>1</sub>") + theme(plot.title = element_markdown(hjust = 0.5, face = "bold"))

gg_interlgcp <- ggplot() +
  geom_line(data = integrated_model$marginals.fixed$intercept_lgcp, mapping = aes(x=x, y=y), color = "red") +
  geom_line(data = sclgcp_model$marginals.fixed$intercept_lgcp, mapping = aes(x=x, y=y), color = "blue") +
  geom_vline(xintercept = -1 + a_optim) + theme_minimal() + 
  labs(title = "(C) Posterior distribution \u03B2<sub>0</sub><sup>*</sup>") + theme(plot.title = element_markdown(hjust = 0.5, face = "bold"))

gg_rw2lgcp <- ggplot() + 
  geom_line(data = data.frame(id = seq_along(t_rw2), rw2 = beta_copy*t_rw2), mapping = aes(x = id, y = rw2), color = "black") +
  geom_line(data = integrated_model$summary.random$time_lgcp, mapping = aes(x = ID, y = mean), color  = "red") +
  geom_ribbon(data = integrated_model$summary.random$time_lgcp %>% rename(., all_of(c(q1='0.025quant', q3='0.975quant'))),
              mapping = aes(x = ID, ymin = q1, ymax = q3), fill = "red", alpha=0.4) +
  geom_line(data = sclgcp_model$summary.random$time_lgcp, mapping = aes(x = ID, y = mean), color = "blue") +
  geom_ribbon(data = sclgcp_model$summary.random$time_lgcp %>% rename(., all_of(c(q1='0.025quant', q3='0.975quant'))),
              mapping = aes(x = ID, ymin = q1, ymax = q3), alpha=0.4, fill = "blue") +
  theme_minimal() + 
  labs(title = "(E) Mean and CI 95% (LGCP-rw2)") + theme(plot.title = element_markdown(hjust = 0.5, face = "bold"))

## Consensus of the latent field ----

### Marginal consensus ----

list_spde_margsum <- list(sc_y1_model$summary.random$sp.idx_gamma, sc_y2_model$summary.random$sp.idx_gamma, sclgcp_model$summary.random$sp.idx_lgcp)
list_yrw2_margsum <- list(sc_y1_model$summary.random$time_gamma, sc_y2_model$summary.random$time_gamma)

consensus_yrw2 <- marginal_consensus(x = list_yrw2_margsum, constr = FALSE)

ggrw2_meanq <- ggplot() + 
  geom_line(data = data.frame(id = seq_along(t_rw2), rw2 = t_rw2), mapping = aes(x = id, y = rw2), color = "black") +
  geom_line(data = integrated_model$summary.random$time_gamma, mapping = aes(x = ID, y = mean), color = "red") +
  geom_ribbon(data = integrated_model$summary.random$time_gamma %>% rename(., all_of(c(q1='0.025quant', q3='0.975quant'))),
              mapping = aes(x = ID, ymin = q1, ymax = q3), fill = "red", alpha=0.4) +
  geom_line(data = consensus_yrw2, mapping = aes(x = ID, y = mean), color = "blue") +
  geom_ribbon(data = consensus_yrw2 %>% rename(., all_of(c(q1='q0.025', q3='q0.975'))),
              mapping = aes(x = ID, ymin = q1, ymax = q3), alpha=0.4, fill = "blue") +
  theme_minimal() +
  labs(title = "(D) Mean and CI 95% (Abundance-rw2)") + theme(plot.title = element_markdown(hjust = 0.5, face = "bold"))

grid.arrange(arrangeGrob(grobs = list(gg_intery, gg_covy, gg_interlgcp, ggrw2_meanq, gg_rw2lgcp), layout_matrix = matrix(c(1,1,2,2,3,3,4,4,4,5,5,5), nrow = 2, byrow = TRUE)))

# ggsave(filename = "./Simulated_example1/Figures/rw2_meanq.png", plot = ggrw2_meanq, device = "png", width = 4500, height = 3000, units = "px", dpi = 600)

consensus_yspde2 <- marginal_consensus(x = list_spde_margsum, constr = FALSE)
order_names <- c("Integrated model (mean)", "Integrated model (q1)", "Integrated model (q3)", "Simulated - Integrated model", 
                 "Sequential consensus (mean)", "Sequential consensus (q1)", "Sequential consensus (q3)", "Simulated - Sequential consensus")
ggspde_meanqres <- ggplot() + 
  geom_tile(data = data.frame(xy_gridpred, z = drop(A_xy_pred %*% integrated_model$summary.random$sp.idx_gamma$mean), approach = "Integrated model (mean)"), mapping = aes(x=x, y=y, fill=z)) +
  geom_tile(data = data.frame(xy_gridpred, z = drop(A_xy_pred %*% integrated_model$summary.random$sp.idx_gamma$`0.025quant`), approach = "Integrated model (q1)"), mapping = aes(x=x, y=y, fill=z)) +
  geom_tile(data = data.frame(xy_gridpred, z = drop(A_xy_pred %*% integrated_model$summary.random$sp.idx_gamma$`0.975quant`), approach = "Integrated model (q3)"), mapping = aes(x=x, y=y, fill=z)) +
  geom_tile(data = data.frame(xy_gridpred, z = drop(A_xy_pred %*% consensus_yspde2$mean), approach = "Sequential consensus (mean)"), mapping = aes(x=x, y=y, fill=z)) +
  geom_tile(data = data.frame(xy_gridpred, z = drop(A_xy_pred %*% consensus_yspde2$q0.025), approach = "Sequential consensus (q1)"), mapping = aes(x=x, y=y, fill=z)) +
  geom_tile(data = data.frame(xy_gridpred, z = drop(A_xy_pred %*% consensus_yspde2$q0.975), approach = "Sequential consensus (q3)"), mapping = aes(x=x, y=y, fill=z)) +
  geom_tile(data = data.frame(xy_gridpred, z = drop(A_xy_pred %*% (integrated_model$summary.random$sp.idx_gamma$mean - u_nodes)), approach = "Simulated - Integrated model"), mapping = aes(x=x, y=y, fill=z)) +
  geom_tile(data = data.frame(xy_gridpred, z = drop(A_xy_pred %*% (consensus_yspde2$mean - u_nodes)), approach = "Simulated - Sequential consensus"), mapping = aes(x=x, y=y, fill=z)) +
  scale_fill_viridis(option = "turbo") + facet_wrap(facets = ~factor(approach, order_names), nrow = 2, dir = "h") +
  # labs(title = "(A) Mean, quantiles of the spatial random effect") + 
  theme_minimal() + 
  # theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(fill = "Values") + theme(strip.text = element_text(size=14, face = "bold"))
  
# ggsave(filename = "./Simulated_example1/Figures/spde_meanqres.png", plot = ggspde_meanqres, device = "png", width = 6000, height = 3000, units = "px", dpi = 600)

drop(A_xy_pred %*% (integrated_model$summary.random$sp.idx_gamma$mean - u_nodes))**2 %>% mean(.)
drop(A_xy_pred %*% (consensus_yspde2$mean - u_nodes))**2 %>% mean(.)

### Product of multivariate gaussian distribution consensus ----

sc_y1_model$misc$configs$contents$tag
idx.Q1 <- 1:sc_y1_model$misc$configs$contents$length[3]
Q1 <- fix.Q(sc_y1_model$misc$configs$config[[1]]$Q[idx.Q1, idx.Q1])
mu1 <- sc_y1_model$misc$configs$config[[1]]$improved.mean[idx.Q1]
sc_y2_model$misc$configs$contents$tag
idx.Q2 <- 1:sc_y2_model$misc$configs$contents$length[3]
Q2 <- fix.Q(sc_y1_model$misc$configs$config[[1]]$Q[idx.Q2, idx.Q2])
mu2 <- sc_y2_model$misc$configs$config[[1]]$improved.mean[idx.Q2]

Qt <- Q1 + Q2
Qt_inv <- solve(Qt)
# Qt_inv <- inla.qinv(Qt)
mut <- Qt_inv %*% (Q1 %*% mu1 + Q2 %*% mu2) %>% drop(.); sdt <- diag(Qt)**(-1/2)
q1t <- qnorm(p = 0.025, mean = mut, sd = sdt) ;q3t <- qnorm(p = 0.975, mean = mut, sd = sdt)
df_prod_yrw2 <- data.frame(ID = seq_along(t_rw2), mean = mut, sd = sdt, q1 = q1t, q3 = q3t)

ggplot() + 
  geom_line(data = data.frame(id = seq_along(t_rw2), rw2 = t_rw2), mapping = aes(x = id, y = rw2), color = "black") +
  geom_line(data = integrated_model$summary.random$time_gamma, mapping = aes(x = ID, y = mean), color = "red") +
  geom_ribbon(data = integrated_model$summary.random$time_gamma %>% rename(., all_of(c(q1='0.025quant', q3='0.975quant'))),
              mapping = aes(x = ID, ymin = q1, ymax = q3), fill = "red", alpha=0.4) +
  geom_line(data = consensus_yrw2, mapping = aes(x = ID, y = mean), color = "blue") +
  geom_ribbon(data = consensus_yrw2 %>% rename(., all_of(c(q1='q0.025', q3='q0.975'))),
              mapping = aes(x = ID, ymin = q1, ymax = q3), alpha=0.4, fill = "blue") +
  geom_line(data = df_prod_yrw2, mapping = aes(x = ID, y = mean), color = "green") +
  geom_ribbon(data = df_prod_yrw2,
              mapping = aes(x = ID, ymin = q1, ymax = q3), alpha=0.4, fill = "green") +
  theme_minimal() +
  labs(title = "(D) Mean and CI 95% (Abundance-rw2)") + theme(plot.title = element_markdown(hjust = 0.5, face = "bold"))

# Integrated model v2 ----

formula_integrated_v2 <- y ~ -1 + intercept_gamma + cov_gamma + f(time_gamma, model = "rw2", constr = TRUE) + f(sp.idx_gamma, model = spde) +
  intercept_lgcp + f(time_lgcp, copy = "time_gamma", fixed = FALSE, hyper = list(beta = list(prior = "normal", params = c(0,0.01)))) + 
  f(sp.idx_lgcp, copy = "sp.idx_gamma", fixed = TRUE)

integrated_model_v2 <- inla(data = inla.stack.data(stk_total), E = inla.stack.data(stk_total)$e,
                            family = c("gamma", "poisson"), formula = formula_integrated_v2,
                            control.predictor = list(A = inla.stack.A(stk_total)),
                            control.compute = list(config = TRUE, waic = TRUE),
                            verbose = FALSE)

# Integrated model v3 ----

formula_integrated_v3 <- y ~ -1 + intercept_gamma + cov_gamma + f(time_gamma, model = "rw2", constr = TRUE) + f(sp.idx_gamma, model = spde) +
  intercept_lgcp + f(time_lgcp, copy = "time_gamma", fixed = FALSE, hyper = list(beta = list(prior = "normal", params = c(0,0.01)))) + 
  f(sp.idx_lgcp, copy = "sp.idx_gamma", fixed = FALSE)

integrated_model_v3 <- inla(data = inla.stack.data(stk_total), E = inla.stack.data(stk_total)$e,
                            family = c("gamma", "poisson"), formula = formula_integrated_v3,
                            control.predictor = list(A = inla.stack.A(stk_total)),
                            control.compute = list(config = TRUE, waic = TRUE),
                            verbose = FALSE)

# Consensus. We can work from the previous results of the models 

log_rw2.ytau <- inla.tmarginal(fun = function(x){log(x)}, marginal = sc_y2_model$marginals.hyperpar$`Precision for time_gamma`) %>% inla.zmarginal(.) %>% data.frame(.)
log_rw2.lgcptau <- inla.tmarginal(fun = function(x){log(x)}, marginal = sclgcp_model$marginals.hyperpar$`Precision for time_lgcp`) %>% inla.zmarginal(.) %>% data.frame(.)
inla.tmarginal(fun = function(x){log(x)}, marginal = sc_y2_model$marginals.hyperpar$`Precision for time_gamma`) %>% plot(. , type = "l")
inla.tmarginal(fun = function(x){log(x)}, marginal = sclgcp_model$marginals.hyperpar$`Precision for time_lgcp`) %>% plot(., type= "l")

alpha.mean_stdev <- c(log_rw2.lgcptau$mean - log_rw2.ytau$mean, sqrt(log_rw2.ytau$sd**2 + log_rw2.lgcptau$sd**2))

x_log.alpha <- seq(from = alpha.mean_stdev[1] - 5*alpha.mean_stdev[2], to = alpha.mean_stdev[1] + 5*alpha.mean_stdev[2], length.out = 1E3)
log_alpha <- data.frame(x = x_log.alpha, y= dnorm(x = x_log.alpha, mean = alpha.mean_stdev[1], sd = alpha.mean_stdev[2]))
log_alpha %>% inla.tmarginal(fun = exp, marginal = .) %>% plot(., type = "l")

## Beta-copy parameter for the temporal trend ----

### Gaussian approximation by using a second-Order Taylor approximation ----

mu_y <- sclgcp_model$summary.random$time_lgcp$mean; sd_y <- sclgcp_model$summary.random$time_lgcp$sd
mu_x <- consensus_yrw2$mean; sd_x <- consensus_yrw2$sd
cor_xy <- 0 # Compare 0 and 1
# cor_xy <- cor(x = sclgcp_model$summary.random$time_lgcp$mean, y = consensus_yrw2$mean)

mu_alpha <- mu_y/mu_x + sd_x**2*mu_y/mu_x**2 - cor_xy*sd_x*sd_y/mu_x**2
sd_alpha <- sqrt(sd_x**2*mu_y**2/mu_x**4 + sd_y**2/mu_x**2 - 2*cor_xy*sd_x*sd_y*mu_y/mu_x**3)

# mu_alpha_join <- mean(mu_alpha)
# sd_alpha_join <- 1/length(sd_alpha)*sqrt(sum(sd_alpha**2))

sd_alpha_join <- sum(sd_alpha**-2)**(-1/2)
mu_alpha_join <- sd_alpha_join**2*sum(sd_alpha**-2*mu_alpha)

x_alpha <- seq(from = mu_alpha_join - 5*sd_alpha_join, to = mu_alpha_join + 5*sd_alpha_join, length.out = 1E3)
df_alpha <- data.frame(x = x_alpha, y = dnorm(x_alpha, mean = mu_alpha_join, sd = sd_alpha_join))

gg_betacopy <- ggplot() + 
  geom_line(data = df_alpha, mapping = aes(x = x, y = y), color = "blue") +
  geom_line(data = integrated_model_v2$marginals.hyperpar$`Beta for time_lgcp`, mapping = aes(x = x, y = y), color = "red") +
  geom_vline(xintercept = (sclgcp_model$summary.random$time_lgcp$mean/consensus_yrw2$mean) %>% median(.), color = "blue") + 
  geom_vline(xintercept = 1.5) + theme_minimal() +
  labs(title = "(A) Posterior distribution  \u03b1-parameter") + theme(plot.title = element_markdown(size = 13, hjust = 0.5, face = "bold"))

### Point estimation of the beta-copy parameter 

(sclgcp_model$summary.random$time_lgcp$mean/consensus_yrw2$mean) %>% mean(.)
(sclgcp_model$summary.random$time_lgcp$mean/consensus_yrw2$mean) %>% median(.)

betacopy_pe <- median(sclgcp_model$summary.random$time_lgcp$mean/consensus_yrw2$mean)

sclgcp_rw2 <- sclgcp_model$summary.random$time_lgcp
sclgcp_rw2$mean <- sclgcp_rw2$mean/betacopy_pe; sclgcp_rw2$sd <- sclgcp_rw2$sd*betacopy_pe
list_rw2_v2 <- list(sc_y1_model$summary.random$time_gamma, sc_y2_model$summary.random$time_gamma, sclgcp_rw2)
consensus_yrw2_v2 <- marginal_consensus(x = list_rw2_v2)

gg_rw2_v2 <- ggplot() + 
  geom_line(data = data.frame(id = seq_along(t_rw2), rw2 = t_rw2), mapping = aes(x = id, y = rw2), color = "black") +
  geom_line(data = integrated_model_v2$summary.random$time_gamma, mapping = aes(x = ID, y = mean), color = "red") +
  geom_ribbon(data = integrated_model_v2$summary.random$time_gamma %>% rename(., all_of(c(q1='0.025quant', q3='0.975quant'))),
              mapping = aes(x = ID, ymin = q1, ymax = q3), fill = "red", alpha=0.4) +
  geom_line(data = consensus_yrw2_v2, mapping = aes(x = ID, y = mean), color = "blue") +
  geom_ribbon(data = consensus_yrw2_v2 %>% rename(., all_of(c(q1='q0.025', q3='q0.975'))),
              mapping = aes(x = ID, ymin = q1, ymax = q3), alpha=0.4, fill = "blue") +
  theme_minimal() +
  labs(title = "(B) Mean and CI 95% (Abundance-rw2)") + theme(plot.title = element_markdown(size = 13, hjust = 0.5, face = "bold"))

gg_rw2_alpha_v2 <- grid.arrange(arrangeGrob(grobs = list(gg_betacopy, gg_rw2_v2), ncol = 2))

ggsave(filename = "./Simulated_example1/Figures/rw2_alpha_v2.png", plot = gg_rw2_alpha_v2, device = "png", dpi = 400, width = 7500, height = 2500, units = "px")

## Beta-copy parameter for the spatial structured effect ----

### Gaussian approximation by using a second-Order Taylor approximation ----

inside <- sapply(1:mesh2d_sim$n, function(i) {
  return(sf::st_intersects(sf::st_point(x = mesh2d_sim$loc[i,1:2]), domainPolygon, sparse=FALSE))
}) %>% which(.)
loc_matrix <- mesh2d_sim$loc[inside,1:2]
# loc_matrix <- rbind(df_ystratified_samples[,1:2], df_ylgcp_samples[,1:2]) %>% as.matrix(.)

mu_y <- drop(fm_basis(x = mesh2d_sim, loc = loc_matrix) %*% sclgcp_model$summary.random$sp.idx_lgcp$mean) 
sd_y <- drop(fm_basis(x = mesh2d_sim, loc = loc_matrix) %*% sclgcp_model$summary.random$sp.idx_lgcp$sd)
mu_x <- drop(fm_basis(x = mesh2d_sim, loc = loc_matrix) %*%  marginal_consensus(x = list(sc_y1_model$summary.random$sp.idx_gamma, sc_y2_model$summary.random$sp.idx_gamma))$mean)
sd_x <-  drop(fm_basis(x = mesh2d_sim, loc = loc_matrix) %*% marginal_consensus(x = list(sc_y1_model$summary.random$sp.idx_gamma, sc_y2_model$summary.random$sp.idx_gamma))$sd)
cor_xy <- 0 # Compare 0 and 1
# cor_xy <- cor(x = sclgcp_model$summary.random$time_lgcp$mean, y = consensus_yrw2$mean)

mu_alpha <- mu_y/mu_x + sd_x**2*mu_y/mu_x**2 - cor_xy*sd_x*sd_y/mu_x**2
sd_alpha <- sqrt(sd_x**2*mu_y**2/mu_x**4 + sd_y**2/mu_x**2 - 2*cor_xy*sd_x*sd_y*mu_y/mu_x**3)

# mu_alpha_join <- median(mu_alpha)
# sd_alpha_join <- 1/length(sd_alpha)*sqrt(sum(sd_alpha**2))

sd_alpha_join <- sum(sd_alpha**-2)**(-1/2)
mu_alpha_join <- sd_alpha_join**2*sum(sd_alpha**-2*mu_alpha)

x_alpha <- seq(from = mu_alpha_join - 5*sd_alpha_join, to = mu_alpha_join + 5*sd_alpha_join, length.out = 1E3)
df_alpha <- data.frame(x = x_alpha, y = dnorm(x_alpha, mean = mu_alpha_join, sd = sd_alpha_join))

ggplot() + 
  geom_line(data = df_alpha, mapping = aes(x = x, y = y), color = "blue") +
  geom_line(data = integrated_model_v3$marginals.hyperpar$`Beta for sp.idx_lgcp`, mapping = aes(x = x, y = y), color = "red") +
  geom_vline(xintercept = 1) + theme_minimal()

### Point estimation of the beta-copy parameter ----

(mu_y/mu_x) %>% median(.)
(mu_y/mu_x)[(mu_y/mu_x) > ((mu_y/mu_x) %>% quantile(x = ., probs = c(0.05, 0.95)))[1] & (mu_y/mu_x) < ((mu_y/mu_x) %>% quantile(x = ., probs = c(0.025, 0.975)))[2]] %>% sd(.)
(mu_y/mu_x) %>% density(.) %>% plot(., type = "l")
