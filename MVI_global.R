require(magrittr)
require(missForest)
require(abind)
# Missing at Random Function ----------------------------------------------
MAR_generate <- function(data, mis_prop = 0.5) {
  all_idx <- which(data != Inf, arr.ind = T)
  rdm_idx <- sample(1:nrow(all_idx), round(nrow(all_idx)*mis_prop))
  slc_idx <- all_idx[rdm_idx, ]
  data_res <- data
  data_res[slc_idx] <- NA
  return(list(data_res = data_res, mis_idx = slc_idx))
}

# column-wise NRMSE calculation -------------------------------------------
nrmse_col <- function(data_imp, data_miss, data_true) {
  idx_miss <- data_miss %>% colSums %>% is.na %>% which
  nrmse_vec <- c()
  for (idx in idx_miss) {
    nrmse_vec <- c(nrmse_vec, nrmse(data_imp[, idx], data_miss[, idx], data_true[, idx]))
  }
  names(nrmse_vec) <- names(idx_miss)
  return(nrmse_vec)
}

# MNAR imputation compare -------------------------------------------------
MNAR_generate <- function (data_c, mis_var = 0.5, var_prop = seq(.3, .6, .1)) {
  data_mis <- data_c
  if (is.numeric(mis_var)) var_mis_list <- sample(1:ncol(data_c), round(ncol(data_c)*mis_var))
  else if (is.character(mis_var)) var_mis_list <- which(colnames(data_c) %in% mis_var)
  for (i in 1:length(var_mis_list)) {
    var_idx <- var_mis_list[i]
    cur_var <- data_mis[, var_idx]
    cutoff <- quantile(cur_var, sample(var_prop, 1))
    cur_var[cur_var < cutoff] <- NA
    data_mis[, var_idx] <- cur_var
  }
  mis_idx_df <- which(is.na(data_mis), arr.ind = T)
  return (list(data_mis = data_mis, mis_idx_df = mis_idx_df))
}

# Scale and recover -------------------------------------------------------
scale_recover <- function(data, method='scale', param_df = NULL) {
  results <- list()
  data_res <- data
  if (!is.null(param_df)) {
    if (method=='scale') {
      data_res[] <- scale(data, center=param_df$mean, scale=param_df$std)
    } else if (method=='recover') {
      data_res[] <- t(t(data)*param_df$std+param_df$mean)
    }
  } else {
    if (method=='scale') {
      param_df <- data.frame(mean=apply(data, 2, function(x) mean(x, na.rm=T)), 
                             std=apply(data, 2, function(x) sd(x, na.rm=T)))
      data_res[] <- scale(data, center=param_df$mean, scale=param_df$std)
    } else {stop('no param_df found for recover...')}
  }
  results[[1]] <- data_res
  results[[2]] <- param_df
  return(results)
}

# Multiplot 4 ggplot2 -----------------------------------------------------
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# Parallel combination ----------------------------------------------------
require(abind)
cbind_abind <- function(a, b) {
  res <- list()
  res$y_imp <- cbind(a$y_imp, b$y_imp)
  res$gibbs_res <- abind(a$gibbs_res, b$gibbs_res, along=2)
  return(res)
}
