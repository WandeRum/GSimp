require(vegan)
require(foreach)
require(doParallel)
require(reshape2)
require(ggplot2)
require(ropls)
source('MVI_global.R')

# MNAR generation and imputation ------------------------------------------
MNAR_gen_imp <- function(data_c, mis_var_prop=seq(.1, .8, .1), var_mis_prop=seq(.1, .8, .1), 
                         impute_list=c('kNN_wrapper', 'SVD_wrapper', 'HM_wrapper', 'QRILC_wrapper'), cores=5) {
  cl <- makeCluster(cores, type='FORK')
  registerDoParallel(cl)
  results <- foreach(prop=mis_var_prop) %dopar% {
    data_m_res <- MNAR_generate(data_c, mis_var=prop, var_prop=var_mis_prop)
    data_m <- data_m_res[[1]]
    mis_idx <- data_m_res[[2]]
    result_list <- list()
    for (i in seq_along(impute_list)) {
      method=impute_list[i]
      print(method)
      result_list[[i]] <- do.call(method, list(data_m))
    }
    result_list[[i+1]] <- mis_idx
    result_list
  }
  stopCluster(cl)
  return(list(results=results, data_c=data_c, impute_list=impute_list))
}

# NRMSE rank calculate and plot -------------------------------------------
NRMSE_rank_cal_plot <- function(impute_results, plot=T, x='Miss_Prop', colors=NULL, shapes=NULL) {
  name_list=gsub('(.*)_.*', '\\1', impute_results$impute_list)
  data_c <- impute_results$data_c
  impute_results_list <- impute_results$results
  nrmse_rank_df <- data.frame(matrix(0, length(impute_results_list), length(name_list) + 2))
  colnames(nrmse_rank_df) <- c(name_list, 'Miss_Prop', 'Miss_Num')
  for (i in seq_along(impute_results_list)) {
    # print(i)
    results_temp <- impute_results_list[[i]]
    mis_idx_temp <- results_temp[[length(results_temp)]]
    var_idx_temp <- unique(mis_idx_temp[, 2])
    data_m_temp <- data_c
    data_m_temp[results_temp[[length(results_temp)]]] <- NA
    temp_list <- rep(0, length(name_list))
    for (j in var_idx_temp) {
      for (k in seq_along(name_list)) {
        temp_list[k] <- nrmse(results_temp[[k]][, j], data_m_temp[, j], data_c[, j])
      }
      rank_temp <- rank(temp_list)
      nrmse_rank_df[i, seq_along(name_list)] <- nrmse_rank_df[i, seq_along(name_list)] + rank_temp
    }
    nrmse_rank_df$Miss_Prop[i] <- mean(is.na(data_m_temp))
    nrmse_rank_df$Miss_Num[i] <- length(var_idx_temp)
  }
  
  nrmse_rank_df_melt <- melt(nrmse_rank_df[, c(name_list, x)], id.vars=x)
  colnames(nrmse_rank_df_melt)[-1] <- c('Method', 'NRMSE_Rank')
  if (plot) {
    if (is.null(colors) & is.null(shapes)) {
      print(ggplot(aes_string(x=x, y='NRMSE_Rank', color='Method'), data=nrmse_rank_df_melt) + 
              geom_point(aes(shape=Method), size=4) + geom_line() + theme_bw())
    } else if (!is.null(colors) & is.null(shapes)) {
      print(ggplot(aes_string(x=x, y='NRMSE_Rank', color='Method'), data=nrmse_rank_df_melt) + 
              geom_point(aes(shape=Method), size=4) + geom_line() + theme_bw() + scale_color_manual(values=colors))
    } else if (is.null(colors) & !is.null(shapes)) {
      print(ggplot(aes_string(x=x, y='NRMSE_Rank', color='Method'), data=nrmse_rank_df_melt) + 
              geom_point(aes(shape=Method), size=4) + geom_line() + theme_bw() + scale_shape_manual(values=shapes))
    } else {
      print(ggplot(aes_string(x=x, y='NRMSE_Rank', color='Method'), data=nrmse_rank_df_melt) + 
              geom_point(aes(shape=Method), size=4) + geom_line() + theme_bw() + 
              scale_shape_manual(values=shapes) + scale_color_manual(values=colors))
    }
  }
  return(list(NRMSE_rank=nrmse_rank_df, NRMSE_rank_melt=nrmse_rank_df_melt))
}

# Procrustes and plot -----------------------------------------------------
Procrustes_cal_plot <- function(impute_results, DR='PCA', lg=F, nPCs=2, outcome=NULL, x='Miss_Prop', plot=T, colors=NULL, shapes=NULL) {
  name_list=gsub('(.*)_.*', '\\1', impute_results$impute_list)
  data_c <- impute_results$data_c
  if (lg) {data_c %<>% log1p}
  impute_results_list <- impute_results$result
  procruste_df <- data.frame(matrix(NA, length(impute_results_list), length(name_list) + 2))
  colnames(procruste_df) <- c(name_list, 'Miss_Prop', 'Miss_Num')
  if (DR == 'PCA') {
    data_c_pca <- prcomp(data_c, scale.=T, center=T)$x[, 1:nPCs]
    for (i in seq_along(impute_results_list)) {
      # print(i)
      results_temp <- impute_results_list[[i]]
      mis_idx_temp <- results_temp[[length(results_temp)]]
      var_idx_temp <- unique(mis_idx_temp[, 2])
      data_m_temp <- data_c
      data_m_temp[results_temp[[length(results_temp)]]] <- NA
      for (j in seq_along(name_list)) {
        # print(name_list[j])
        data_imp_temp <- results_temp[[j]]
        if (lg) {data_imp_temp %<>% log1p}
        procruste_df[i, j] <- procrustes(data_c_pca, prcomp(data_imp_temp, scale.=T, center=T)$x[, 1:nPCs], symmetric=T)$ss
      }
      procruste_df$Miss_Prop[i] <- mean(is.na(data_m_temp))
      procruste_df$Miss_Num[i] <- length(var_idx_temp)
    }
  }
  else if (DR == 'PLS' & !is.null(outcome)) {
    data_c_pls <- data.frame(opls(data_c, outcome, predI=nPCs, permI=0, plotL=F, printL=F)@scoreMN)
    for (i in seq_along(impute_results_list)) {
      # print(i)
      results_temp <- impute_results_list[[i]]
      mis_idx_temp <- results_temp[[length(results_temp)]]
      var_idx_temp <- unique(mis_idx_temp[, 2])
      data_m_temp <- data_c
      data_m_temp[results_temp[[length(results_temp)]]] <- NA
      for (j in seq_along(name_list)) {
        data_imp_temp <- results_temp[[j]]
        if (lg) {data_imp_temp %<>% log1p}
        procruste_df[i, j] <- procrustes(data_c_pls, 
                                         data.frame(opls(data_imp_temp, outcome, predI=nPCs, permI=0, plotL=F, printL=F)@scoreMN),
                                         symmetric=T)$ss
      }
      procruste_df$Miss_Prop[i] <- mean(is.na(data_m_temp))
      procruste_df$Miss_Num[i] <- length(var_idx_temp)
    }
  }
  procruste_df_melt <- melt(procruste_df[, c(name_list, x)], id.vars=x)
  colnames(procruste_df_melt)[-1] <- c('Method', 'Pro_SS')
  if (plot) {
    if (is.null(colors) & is.null(shapes)) {
      print(ggplot(aes_string(x=x, y='Pro_SS', color='Method'), data=procruste_df_melt) + geom_point(aes(shape=Method), size=4) + geom_line() + theme_bw())
    } else if (!is.null(colors) & is.null(shapes)) {
      print(ggplot(aes_string(x=x, y='Pro_SS', color='Method'), data=procruste_df_melt) + geom_point(aes(shape=Method), size=4) + geom_line() + theme_bw() + scale_color_manual(values=colors))
    } else if (is.null(colors) & !is.null(shapes)) {
      print(ggplot(aes_string(x=x, y='Pro_SS', color='Method'), data=procruste_df_melt) + geom_point(aes(shape=Method), size=4) + geom_line() + theme_bw() + scale_shape_manual(values=shapes))
    } else {
      print(ggplot(aes_string(x=x, y='Pro_SS', color='Method'), data=procruste_df_melt) + geom_point(aes(shape=Method), size=4) + geom_line() + theme_bw() + 
              scale_shape_manual(values=shapes) + scale_color_manual(values=colors))
    }
  }
  return(list(Pro_SS=procruste_df, Pro_SS_melt=procruste_df_melt))
}


# Univariate analysis and plot --------------------------------------------
Ttest_cor_cal_plot <- function(impute_results, group=NULL, plot=T, x='Miss_Prop', cor='P', colors=NULL, shapes=NULL) {
  name_list <- gsub('(.*)_.*', '\\1', impute_results$impute_list)
  data_c <- impute_results$data_c
  impute_results_list <- impute_results$results
  Ttest_p_c <- apply(data_c, 2, function(x) t.test(x ~ group)$p.value)
  
  pcor_df <- data.frame(matrix(NA, length(impute_results_list), length(name_list) + 2))
  scor_df <- data.frame(matrix(NA, length(impute_results_list), length(name_list) + 2))
  colnames(pcor_df) <- c(name_list, 'Miss_Prop', 'Miss_Num')
  colnames(scor_df) <- c(name_list, 'Miss_Prop', 'Miss_Num')
  
  for (i in seq_along(impute_results_list)) {
    # print(i)
    results_temp <- impute_results_list[[i]]
    mis_idx_temp <- results_temp[[length(results_temp)]]
    var_idx_temp <- unique(mis_idx_temp[, 2])
    data_m_temp <- data_c
    data_m_temp[results_temp[[length(results_temp)]]] <- NA
    for (j in seq_along(name_list)) {
      pcor_df[i, j] <- cor.test(log(Ttest_p_c[var_idx_temp]), 
                                log(apply(results_temp[[j]][, var_idx_temp], 2, 
                                          function(x) t.test(x ~ group)$p.value)))$estimate
      scor_df[i, j] <- cor.test(Ttest_p_c[var_idx_temp], 
                                apply(results_temp[[j]][, var_idx_temp], 2, 
                                          function(x) t.test(x ~ group)$p.value), method='spearman')$estimate
    }
    pcor_df$Miss_Prop[i] <- mean(is.na(data_m_temp))
    scor_df$Miss_Prop[i] <- mean(is.na(data_m_temp))
    pcor_df$Miss_Num[i] <- length(var_idx_temp)
    scor_df$Miss_Num[i] <- length(var_idx_temp)
  }
  pcor_df_melt <- melt(pcor_df[, c(name_list, x)], id.var=x)
  colnames(pcor_df_melt)[-1] <- c('Method', 'Cor')
  scor_df_melt <- melt(scor_df[, c(name_list, x)], id.var=x)
  colnames(scor_df_melt)[-1] <- c('Method', 'Cor')
  if (plot & cor == 'P') {
    if (is.null(colors) & is.null(shapes)) {
      print(ggplot(aes_string(x=x, y='Cor', color='Method'), data=pcor_df_melt) + geom_point(aes(shape=Method), size=4) + geom_line() + theme_bw() + geom_point(aes(shape=Method), size=4))
    } else if (!is.null(colors) & is.null(shapes)) {
      print(ggplot(aes_string(x=x, y='Cor', color='Method'), data=pcor_df_melt) + geom_point(aes(shape=Method), size=4) + geom_line() + theme_bw() + geom_point(aes(shape=Method), size=4) + scale_color_manual(values=colors))
    } else if (is.null(colors) & !is.null(shapes)) {
      print(ggplot(aes_string(x=x, y='Cor', color='Method'), data=pcor_df_melt) + geom_point(aes(shape=Method), size=4) + geom_line() + theme_bw() + geom_point(aes(shape=Method), size=4) + scale_shape_manual(values=shapes))
    } else {
      print(ggplot(aes_string(x=x, y='Cor', color='Method'), data=pcor_df_melt) + geom_point(aes(shape=Method), size=4) + geom_line() + theme_bw() + geom_point(aes(shape=Method), size=4) + scale_shape_manual(values=shapes) + scale_color_manual(values=colors))
    }
  } else if (plot & cor == 'S') {
    if (is.null(colors) & is.null(shapes)) {
      print(ggplot(aes_string(x=x, y='Cor', color='Method'), data=scor_df_melt) + geom_point(aes(shape=Method), size=4) + geom_line() + theme_bw() + geom_point(aes(shape=Method), size=4))
    } else if (!is.null(colors) & is.null(shapes)) {
      print(ggplot(aes_string(x=x, y='Cor', color='Method'), data=scor_df_melt) + geom_point(aes(shape=Method), size=4) + geom_line() + theme_bw() + geom_point(aes(shape=Method), size=4) + scale_color_manual(values=colors))
    } else if (is.null(colors) & !is.null(shapes)) {
      print(ggplot(aes_string(x=x, y='Cor', color='Method'), data=scor_df_melt) + geom_point(aes(shape=Method), size=4) + geom_line() + theme_bw() + geom_point(aes(shape=Method), size=4) + scale_shape_manual(values=shapes))
    } else {
      print(ggplot(aes_string(x=x, y='Cor', color='Method'), data=scor_df_melt) + geom_point(aes(shape=Method), size=4) + geom_line() + theme_bw() + geom_point(aes(shape=Method), size=4) + scale_shape_manual(values=shapes) + scale_color_manual(values=colors))
    }
  }
  return(list(P_cor=pcor_df, S_cor=scor_df, P_cor_melt=pcor_df_melt, S_cor_melt=scor_df_melt))
}

Ttest_P_cal_plot <- function(impute_results, group=NULL, plot=T, p_cut=.05, x='Miss_Prop', colors=NULL, shapes=NULL) {
  name_list <- gsub('(.*)_.*', '\\1', impute_results$impute_list)
  data_c <- impute_results$data_c
  impute_results_list <- impute_results$results
  Ttest_p_c <- apply(data_c, 2, function(x) t.test(x ~ group)$p.value)
  sig_idx <- which(Ttest_p_c < p_cut)
  
  P_df <- data.frame(matrix(NA, length(impute_results_list), length(name_list) + 2))
  colnames(P_df) <- c(name_list, 'Miss_Prop', 'Miss_Num')
  
  for (i in seq_along(impute_results_list)) {
    # print(i)
    results_temp <- impute_results_list[[i]]
    mis_idx_temp <- results_temp[[length(results_temp)]]
    var_idx_temp <- unique(mis_idx_temp[, 2])
    data_m_temp <- data_c
    data_m_temp[results_temp[[length(results_temp)]]] <- NA
    for (j in seq_along(name_list)) {
      sig_idx_temp <- which(apply(results_temp[[j]], 2, 
                                  function(x) t.test(x ~ group)$p.value) < p_cut)
      P_df[i, j] <- mean(sig_idx %in% sig_idx_temp)
    }
    P_df$Miss_Prop[i] <- mean(is.na(data_m_temp))
    P_df$Miss_Num[i] <- length(var_idx_temp)
  }
  P_df_melt <- melt(P_df[, c(name_list, x)], id.var=x)
  colnames(P_df_melt)[-1] <- c('Method', 'Po')
  if (plot) {
    if (is.null(colors) & is.null(shapes)) {
      print(ggplot(aes_string(x=x, y='Po', color='Method'), data=P_df_melt) + geom_point(aes(shape=Method), size=4) + geom_line() + theme_bw() + geom_point(aes(shape=Method), size=4))
    } else if (!is.null(colors) & is.null(shapes)) {
      print(ggplot(aes_string(x=x, y='Po', color='Method'), data=P_df_melt) + geom_point(aes(shape=Method), size=4) + geom_line() + theme_bw() + geom_point(aes(shape=Method), size=4) + scale_color_manual(values=colors))
    } else if (is.null(colors) & !is.null(shapes)) {
      print(ggplot(aes_string(x=x, y='Po', color='Method'), data=P_df_melt) + geom_point(aes(shape=Method), size=4) + geom_line() + theme_bw() + geom_point(aes(shape=Method), size=4) + scale_shape_manual(values=shapes))
    } else {
      print(ggplot(aes_string(x=x, y='Po', color='Method'), data=P_df_melt) + geom_point(aes(shape=Method), size=4) + geom_line() + theme_bw() + geom_point(aes(shape=Method), size=4) + scale_shape_manual(values=shapes) + scale_color_manual(values=colors))
    }
  }
  return(list(P_df=P_df, P_df_melt=P_df_melt))
}
