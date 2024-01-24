library(dplyr)
library(tidybayes)
library(brms)
library(tidyr)
library(rstan)
options(mc.cores = parallel::detectCores())
library(ggplot2)
library(bayesplot)
library(bayestestR)


#### forest plot ####

custom_forest_plot <- function(mod, variable) {
  post <- as_draws_array(mod) 
  ci.95 <- ci(post, method = 'ETI') %>% filter(Parameter %in% variable) %>% rowwise() %>% mutate(mu = (CI_low+CI_high)/2)
  ci.85 <- ci(post, method = 'ETI',ci = 0.85) %>% filter(Parameter %in%  variable) %>% rowwise() %>% mutate(mu = (CI_low+CI_high)/2)
  ci.50 <- ci(post, method = 'ETI', ci = 0.5) %>% filter(Parameter %in%  variable) %>% rowwise() %>% mutate(mu = (CI_low+CI_high)/2)  
  plot <- ggplot(ci.95, aes(x = mu, y = Parameter))+
    geom_vline(xintercept = 0, color = 'grey', linewidth = 0.8)+
    geom_linerange(aes(xmin = CI_low, xmax = CI_high), linewidth = 1, color = '#464645')+
    geom_linerange(data = ci.85, aes(xmin = CI_low, xmax = CI_high), linewidth = 2.6, color = '#a40b0b')+
    geom_point(shape = 21, size = 5, color = 'black', fill = '#dc3a3a')+
    theme_bw(base_size = 20, base_line_size = 20/44)+
    theme(axis.text = element_text(color = 'black'))+
    labs(x = NULL, y = NULL)
  return(plot)
  }



#### prediction regression ####

custom_reg_plot <- function(variable, mod, CI = 0.9, dpar = NULL) {
  plot <- conditional_effects(mod, effects = variable, prob = CI, dpar = dpar)[[1]] %>%
    mutate(Density = factor(Density, levels = c('0', '1'))) %>%
    ggplot(aes(x = .[,1], y = `estimate__`, fill = Density))+
    geom_line(aes(linetype = Density, color = Density), 
              linewidth = 1)+
    geom_ribbon(aes(ymin = `lower__`, ymax = `upper__`), alpha = 0.2)+
    theme_bw(base_size = 16)+
    theme(legend.position = 'none')+
    scale_fill_manual(values = c('#f1c40f', '#c9183a'))+
    scale_color_manual(values = c('#f1c40f', '#c9183a'))+
    xlab(variable)
  return(plot)
}

custom_reg_plot_hurdle <- function(variable, mod, CI = 0.9, dpar = NULL) {
  plot <- conditional_effects(mod, effects = variable, prob = CI, dpar = dpar)[[1]] %>%
    mutate(Density = factor(Density, levels = c('0', '1'))) %>%
    ggplot(aes(x = .[,1], y = log10(`estimate__`+ 1), fill = Density))+
    geom_line(aes(linetype = Density, color = Density), 
              linewidth = 1)+
    geom_ribbon(aes(ymin = log10(`lower__`+1), ymax = log10(`upper__`+1)), alpha = 0.2)+
    theme_bw(base_size = 22)+
    theme(legend.position = 'none')+
    scale_fill_manual(values = c('#f1c40f', '#c9183a'))+
    scale_color_manual(values = c('#f1c40f', '#c9183a'))+
    xlab(variable)
  return(plot)
}



#### prediction factor ####

custom_fac_plot <- function(variable, mod, CI = 0.89) {
  plot <- conditional_effects(mod, effects = variable, prob = CI)[[1]] %>%
    ggplot(aes(x = .[,1], y = `estimate__`, fill = as.factor(Density),  # x = .[,1] if Density is behind the continuous variable and .[,2] if it's before
               color = as.factor(Density)))+
    geom_errorbar(aes(ymin = lower__, ymax = upper__), width = 0.3, linewidth = 1)+
    geom_point(size = 3, color = 'black')+
    theme_bw(base_size = 22)+
    theme(legend.position = 'none')+
    scale_fill_manual(labels = c('0' = 'Low', '1' = 'High'), values = c('#f1c40f','#c9183a'))+
    scale_color_manual(labels = c('0' = 'Low', '1' = 'High'), values = c('#f1c40f','#c9183a'))+
    xlab(variable)
  return(plot)
}



### calculating proportion of posterior > or < 0

posterior_ratio <- function(mod, ndraws = 4000) {
  post <- as_draws_df(mod) %>% select(starts_with("b")) %>%
    mutate(across(starts_with("b"), ~ ifelse( . > 0, 1, 0 ), .names = "more0_{.col}")) %>%
    mutate(across(starts_with("b"), ~ ifelse( . < 0, 1, 0 ), .names = "less0_{.col}")) %>%
    select(!starts_with("b")) %>%
    summarise(across(everything(), list(sum))) %>%
    mutate(across(everything(), ~ ./ndraws))  %>%
    mutate(across(everything(), as.character))  %>%
    pivot_longer(cols = everything(), names_pattern = "(.*0_)(.*$)", names_to = c('sign', 'parameter')) %>%
    pivot_wider(names_from = sign, values_from = value)
  return(post)
}



posterior_density <- function(mod1, mod2, mod3, mod4, mod5, variable) {
  p.dwva <- as_draws_df(mod1) %>% select(dwva = variable)
  p.dwvb <- as_draws_df(mod2) %>% select(dwvb = variable) 
  p.bqcv <- as_draws_df(mod3) %>% select(bqcv = variable) 
  p.abpv <- as_draws_df(mod4) %>% select(abpv = variable) 
  p.sbv <- as_draws_df(mod5) %>% select(sbv = variable) 
  post.path <- cbind(p.dwva, p.dwvb, p.bqcv, p.abpv, p.sbv) %>% 
    pivot_longer(cols = everything(), names_to = 'Virus', values_to = 'Post') %>%
  ggplot(aes(x = Post, col = Virus, fill = Virus))+
    geom_density(alpha = 0.1, linewidth = 1.2)+
    theme_bw(base_size = 20)+
    geom_vline(xintercept = 0, linetype = 'dashed', linewidth = 1)
return(post.path)
}

posterior_density_inv <- function(mod1, mod2, mod3, mod4, mod5, variable) {
  p.dwva <- as_draws_df(mod1) %>% select(dwva = variable)
  p.dwvb <- as_draws_df(mod2) %>% select(dwvb = variable) 
  p.bqcv <- as_draws_df(mod3) %>% select(bqcv = variable) 
  p.abpv <- as_draws_df(mod4) %>% select(abpv = variable) 
  p.sbv <- as_draws_df(mod5) %>% select(sbv = variable) 
  post.path <- cbind(p.dwva, p.dwvb, p.bqcv, p.abpv, p.sbv) %>% 
    pivot_longer(cols = everything(), names_to = 'Virus', values_to = 'Post') %>%
    ggplot(aes(x = Post*-1, col = Virus, fill = Virus))+
    geom_density(alpha = 0.1, linewidth = 1.2)+
    theme_bw(base_size = 20)+
    geom_vline(xintercept = 0, linetype = 'dashed', linewidth = 1)
  return(post.path)
}

posterior_density2 <- function(mod1, mod2, mod3, mod4, mod5, mod6, mod7, variable) {
  p.dwva <- as_draws_df(mod1) %>% select(dwva = variable)
  p.dwvb <- as_draws_df(mod2) %>% select(dwvb = variable) 
  p.bqcv <- as_draws_df(mod3) %>% select(bqcv = variable) 
  p.abpv <- as_draws_df(mod4) %>% select(abpv = variable) 
  p.sbv <- as_draws_df(mod5) %>% select(sbv = variable) 
  p.tryp <- as_draws_df(mod6) %>% select(tryp = variable) 
  p.neo <- as_draws_df(mod7) %>% select(neo = variable) 
  post.path <- cbind(p.dwva, p.dwvb, p.bqcv, p.abpv, p.sbv, p.tryp, p.neo) %>% 
    pivot_longer(cols = everything(), names_to = 'Virus', values_to = 'Post') %>%
    ggplot(aes(x = Post, col = Virus, fill = Virus))+
    geom_density(alpha = 0.1, linewidth = 1.2)+
    theme_bw(base_size = 20)+
    geom_vline(xintercept = 0, linetype = 'dashed', linewidth = 1)
  return(post.path)
}

#### summary

custom_summary <- function(mod) {
  sum <- summary(mod)[['fixed']] %>% select(Estimate, Est.Error, Bulk_ESS, Tail_ESS) %>% cbind(posterior_ratio(mod)[-c(1)]) %>% 
    cbind(ci(as_draws_array(mod), ci = 0.9, method = 'ETI') %>% 
            filter(grepl('b_', Parameter)) %>% select(-c(CI, Parameter)) %>% rename(CI_low_90 = CI_low, CI_high_90 = CI_high)) %>% 
    mutate(pd = as.numeric(ifelse(more0_ > less0_, more0_, less0_))) %>% select(!c(more0_, less0_)) %>% mutate(across(everything(), ~ round(., 2))) %>%
    relocate(c(CI_low_90, CI_high_90), .after = Est.Error) %>% mutate(Response = mod[["formula"]][["resp"]], Predictor = rownames(.)) %>% relocate(c(Response, Predictor), .before = Estimate)
  rownames(sum) <- NULL
  return(sum)
}
