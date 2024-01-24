source('data_prep.R') # load all the data and packages

#micro.w2022 <- data.w2022 %>% left_join(micro %>% select(!c(Density)), by = 'Site') %>% drop_na(No_bees) %>%
#  mutate(Lat.s = scale(Lat), Long.s = scale(Long), Score.s = scale(Beekeeper_score), Wind.p = as.factor(Wind_protection))


library(dagitty)

dag <- dagitty("dag{H -> T -> V -> P -> B;
                    H -> T -> P; 
                    H -> B; 
                    D -> T;
                    D -> B;
                    V -> B
                    B [outcome]
                    T [unobserved]
                    H [exposure]
                    D [exposure]
                    V [observed]
                    P [observed]
               }")


plot(dag)
impliedConditionalIndependencies(dag)
backDoorGraph(dag)


### modelling longitudinal data:
# 1) with varying slope of month (time) and random intercept of colony
# 2) with gaussian processes of a distance between each month matrix (and random intercept of colony?)
# 3) with month as a covariate but with a spline and random intercept of colony

###### 1. COLONY GROWTH MODELS ######
set.seed(11)

### best performing - student family
bee.both <- brm(No_bees.s ~ Round + (fl.sum.s + H2) * Density  + (1+Round|Site/Sample), family = student(), 
                data = data.w2022, prior = c(prior(normal(0,10), class = 'b')), control = list(adapt_delta = 0.95))

#save(bee.both, file = 'Results/Bee_mod240116.RData')

#### underestimating, the zeros are biasing the result!
bee.both.nb <- brm(No_bees ~ Round + (fl.sum.s + H2) * Density  + (1 + Round|Site/Sample), family = negbinomial, 
               data = data.w2022, prior = c(prior(normal(0,10), class = 'b'),
                                            prior(student_t(2,3,2.5), class = 'Intercept')), control = list(adapt_delta = 0.95))

summary(bee.both)
pp_check(bee.both, ndraws = 100)

mod.both %>% mcmc_pairs(pars = c('b_Round', 'b_fl.sum.s', 'b_H2', 'b_Density1', 'b_fl.sum.s:Density1', 'b_H2:Density1'))
p_direction(mod.both)
conditional_effects(mod.both, effects = 'fl.sum.s:Density')

### checking out some contrast testing
hypothesis(mod.H, "H2 + H2:Density1 = 0")
emmeans(mod.H, specs = pairwise ~ H2*Density)
emmeans(mod.fl, specs = pairwise ~ fl.sum.s:Density)
fit.d0 <- fitted(mod.H, newdata = expand.grid(Round = 3,
                             H2 = seq(min(data.w2022$H2, max(data.w2022$H2), by = 0.1)),
                             Density = 0,
                             Site = 'Goe1214',
                             Sample = 'Goe1214H1'), summary = T)

df_post <- as_draws_df(mod.H)
hm <- data.frame(draws = df_post$b_Intercept + df_post$`b_H2:Density1`)

###### 2. VARROA MODELS ######

var.both <- brm(Varroa.p2 ~ (fl.sum.s + H2) * Density + Round + (1 + Round|Site/Sample),  chains = 4, family = Gamma(link = 'log'), 
            data.v2022, prior = c(prior(normal(0,10), class = 'b')), control = list(adapt_delta = 0.95))

var.both %>% mcmc_pairs(pars = c('b_Round', 'b_fl.sum.s', 'b_H2', 'b_Density1', 'b_fl.sum.s:Density1', 'b_H2:Density1'))

#save(var.both, file = 'Results/Varroa_mod240116.RData')
summary(var.both)
conditional_effects(var.both)
post.radio <- posterior_ratio(var.both)
p_direction(var.both)
post <- as_draws_array(var.both) 
ci.90 <- ci(post, method = 'ETI', ci = 0.9)
pp_check(var.both, ndraws = 100)+ coord_cartesian(xlim = c(-1, 6))

## other stuff

### remove 4th sampling because acid was used??? much bigger varroa fall

### SEM
  
sem1 <- brm(bf(Varroa.p2 ~ (H2 + fl.sum.s) * Density + Round + (1+Round|Site/Sample), family = Gamma(link = 'log'))+
              bf(No_bees.s ~ Varroa.p2 + (H2 + fl.sum.s) * Density + Round + (1+Round|Site/Sample), family = student)+
              set_rescor(F), data = data.v2022)
sem2 <- brm(bf(Varroa.p2 ~ (H2 + fl.sum.s) * Density + Round + (1+Round|Site/Sample), family = Gamma(link = 'log'))+
              bf(No_bees ~ Varroa.p2 + (H2 + fl.sum.s) * Density + Round + (1+Round|Site/Sample), family = negbinomial)+
              set_rescor(F), data = data.v2022)

summary(sem1)
summary(sem2)
conditional_effects(sem1)
conditional_effects(sem2)
p_direction(sem1)
p_direction(sem2)
pp_check(sem2, resp = 'Nobees', ndraws = 100)

sem2.dwvb <- brm(bf(Varroa.p2 ~ (H2 + fl.sum.s) * Density + (1|Site/Sample), family = Gamma(link = 'log'))+
              bf(dwvb ~ (H2 + fl.sum.s) * Density + Varroa.p2 + (1|Site/Sample), family = bernoulli)+
              bf(No_bees ~ Varroa.p2 + (H2 + fl.sum.s) * Density + dwvb + (1|Site/Sample), family = negbinomial)+
              set_rescor(F), data = data.path.w)

summary(sem2.dwvb)
p_direction(sem2.dwvb)
conditional_effects(sem2.dwvb)


#### 4. ABS QUANTIFICATION ####

#data.path.w2 <- data.w2022 %>% filter(Month %in% c('July', 'August')) %>% drop_na(dwva) %>%
#  mutate(Month = factor(Month, levels = c('July', 'August'))) %>%  mutate(Quant.m_DWVA = ifelse(dwva == 0, 0, Quant.m_DWVA),
#                                                                          Quant.m_DWVB = ifelse(dwvb == 0, 0, Quant.m_DWVB),
#                                                                          Quant.m_BQCV = ifelse(bqcv == 0, 0, Quant.m_BQCV),
#                                                                          Quant.m_ABPV = ifelse(abpv == 0, 0, Quant.m_ABPV),
#                                                                          Quant.m_SBV3 = ifelse(sbv == 0, 0, Quant.m_SBV3)) %>%
#  mutate(fl.sum.s = scale(fl.sum))

data.w2022 %>% select(Site, Month, fl.sum.s) %>% distinct()
data.path.w %>% select(Site, Month, fl.sum.s)

#### 4.1 PATHOGENS HURDLE #####

#hu <- rlogis(10000, 2.96, 1)
#plot(density(inv_logit_scaled(hu)))

data.path.w <- data.path.w %>% mutate(across(Quant.m_ABPV:Quant.m_SBV3, function(x) ifelse(x < 1, 0, x)))

## DWV-A

dwva2 <- glmer(abpv ~ (fl.sum.s + H2) * Density + Varroa.s + (1|Site/Sample), family = 'binomial',data.path.w)
summary(dwva2)
plot(ggpredict(dwva2, terms = c('H2', 'Density')))


dwva1.2 <- brm(dwva ~ (fl.sum.s + H2) * Density + Varroa.s + (1|Site/Sample), family = 'bernoulli',
             data.path.w, prior = c(prior(normal(0,5), class = 'b')), control = list(adapt_delta = 0.98, max_treedepth = 13))

custom_reg_plot(dwva1.2, variable = 'H2:Density')
summary(dwva1.2)

dwva1 <- brm(bf(Quant.m_DWVA ~ (fl.sum.s + H2) * Density + Varroa.s + (1|Site/Sample),
                hu ~ (fl.sum.s + H2) * Density + Varroa.s + (1|Site/Sample)), family = hurdle_lognormal(), 
             data.path.w, prior = c(prior(logistic(0,1), dpar = 'hu', class = 'Intercept'),
                                    prior(normal(0,5), class = 'b')), control = list(adapt_delta = 0.98, max_treedepth = 13))

save(dwva1, file = "Results/dwva_mod240116.RData")
custom_summary(dwva1)

pred <- posterior_predict(dwva1)
y1 <- dwva1[["data"]][["Quant.m_DWVA"]]
ppc_dens_overlay(y = log1p(y1), 
                 yrep = log1p(pred[1:length(y1),]))+
  coord_cartesian(ylim=c(0,4))

conditional_effects(dwva1, effects = 'H2:Density', dpar = 'hu')
conditional_effects(dwva1, effects = 'H2:Density', dpar = 'mu')
conditional_effects(dwva1, effects = 'fl.sum.s:Density', dpar = 'hu')
conditional_effects(dwva1, effects = 'fl.sum.s:Density', dpar = 'mu')



## DWV-B
#data.path.w$Density <- factor(data.path.w$Density, levels = c('1', '0'))

dwvb1 <- brm(bf(Quant.m_DWVB ~ (fl.sum.s + H2) * Density  + Varroa.s  + (1|Site/Sample),
                hu ~ (fl.sum.s + H2) * Density + Varroa.s + (1|Site/Sample)), family = hurdle_lognormal(),
             prior = c(prior(logistic(0,1), dpar = 'hu', class = 'Intercept'),
                       prior(normal(0,5), class = 'b')), control = list(adapt_delta = 0.98),
             data.path.w) 

save(dwvb1, file = "Results/dwvb_mod240116.RData")
hm <-  epred_draws(dwvb1, newdata = expand.grid(Density = c('0','1'),
                                          fl.sum.s = seq(min(data.path.w$fl.sum.s), max(data.path.w$fl.sum.s), by = 0.1),
                                          Site = 'Goe1214',
                                          Sample = 'Goe1214H1',
                                          Varroa.s = 0))


custom_summary(dwvb1)
conditional_effects(dwvb1, dpar = 'mu')
conditional_effects(dwvb1, dpar = 'hu')

pred <- posterior_predict(dwvb1)
y1 <- dwvb1[["data"]][["Quant.m_DWVB"]]
ppc_dens_overlay(y = log1p(y1), yrep = log1p(pred[1:length(y1),]))


conditional_effects(dwvb1, dpar = 'mu')
dwvb1 %>% emmeans( ~ Density, var = 'fl.sum.s', epred = T)
dwvb1 %>% emtrends( ~ Density, var = 'fl.sum.s', dpar = 'hu')


## BQCV 

bqcv1 <- brm(bf(Quant.m_BQCV ~ (fl.sum.s + H2) * Density + (1|Site/Sample),
                hu ~ (fl.sum.s + H2) * Density + (1|Site/Sample)), family = hurdle_lognormal(),
             prior = c(prior(logistic(0,1), dpar = 'hu', class = 'Intercept'),
                       prior(normal(0,5), class = 'b')), control = list(adapt_delta = 0.98),
             data.path.w) 

save(bqcv1, file = "Results/bqcv_mod240116.RData")
custom_summary(bqcv1)

conditional_effects(bqcv1, dpar = 'hu')
conditional_effects(bqcv1, dpar = 'mu')


## ABPV

abpv1 <- brm(bf(Quant.m_ABPV ~ (fl.sum.s + H2) * Density + Varroa.s + (1|Site/Sample),
                hu ~ (fl.sum.s + H2) * Density + Varroa.s + (1|Site/Sample)), family = hurdle_lognormal(), 
             prior = c(prior(logistic(0,1), dpar = 'hu', class = 'Intercept'),
                       prior(normal(0,5), class = 'b')), control = list(adapt_delta = 0.98), 
             data.path.w)
save(abpv1, file = "Results/abpv_mod240116.RData")

custom_summary(abpv1)
conditional_effects(abpv1, dpar = 'mu')
conditional_effects(abpv1, dpar = 'hu')

pred <- posterior_predict(abpv2)
y1 <- abpv2[["data"]][["Quant.m_ABPV"]]
ppc_dens_overlay(y = log1p(y1), 
                            yrep = log1p(pred[1:length(y1),]))

## sbv
sbv1 <- brm(bf(Quant.m_SBV3 ~ (fl.sum.s + H2) * Density + (1|Site/Sample),
               hu ~ (fl.sum.s + H2) * Density +  (1|Site/Sample)), family = hurdle_lognormal(), 
            prior = c(prior(logistic(0,1), dpar = 'hu', class = 'Intercept'),
                      prior(normal(0,5), class = 'b')), control = list(adapt_delta = 0.98), data.path.w)
save(sbv1, file = "Results/sbv_mod240116.RData")


custom_summary(sbv1)

y1 <- sbv1[["data"]][["Quant.m_SBV3"]]
pred <- posterior_predict(sbv1)
ppc_dens_overlay(y = log1p(y1), 
                            yrep = log1p(pred[1:124,]))

conditional_effects(sbv1, dpar = 'hu')
conditional_effects(sbv1, dpar = 'mu')

pathogen.summary <- lapply(list(dwva1, dwvb1, bqcv1, abpv1, sbv1), custom_summary) %>% bind_rows()


###### 5. SURVIVAL #######

surv1 <- brm(Status ~ dwvb.A * Density + Varroa.m.s + (1|Site), data.surv, 
             prior = c(prior(normal(0,2), class = 'b'),
                       prior(normal(0,2), class = 'Intercept')), family = bernoulli)

surv1 <- brm(Status ~ (fl.sum.s + H2) * Density + Varroa.m.s + (1|Site), data.surv, family = bernoulli)

#save(surv1, file = 'Results/Survival_mod240116.RData')


summary(surv1)
conditional_effects(surv1)
posterior_ratio(surv1)
p_direction(surv1)

custom_reg_plot(surv1, variable = 'fl.sum.s:Density')+
  labs(x = 'Flower covered area [ha]', y = 'Probability of survival')+
  geom_point(data = data.surv, aes(fl.sum.s, Status, col = Density))
custom_reg_plot(surv1, variable = 'Ann.fl:Density')
custom_reg_plot(surv1, variable = 'SNH:Density')

ggplot(data.surv, aes(fl.sum.s, Status, col = Density))+
  geom_point(position = position_dodge2(width = 0.2))

p <- posterior_ratio(surv1)


#### other stuff

surv <- data.full %>% select(Sample, Site, Status, Round, No_bees, Varroa_bottom, Varroa_p) %>% 
  group_by(Sample, Round) %>% summarise(No_bees = mean(No_bees), Varroa_bottom = mean(Varroa_bottom), Varroa_p = mean(Varroa_p)) %>%
  filter(Round != 0) %>% filter(Round != 4)

cummulative <- data.full %>% select(Sample, Site, Status, Round, No_bees, Varroa_bottom, Varroa_p) %>% 
  group_by(Sample, Round) %>% summarise(No_bees = mean(No_bees), Varroa_bottom = mean(Varroa_bottom), Varroa_p = mean(Varroa_p)) %>%
  filter(Round != 0) %>% filter(Round != 4) %>% 
  ## making cumulative sliding window
  group_by(Sample) %>% arrange(desc(Round)) %>% reframe(No_bees_c = runner(No_bees, f = sum, na_pad = T),
                                                        Varroa_p_c = runner(Varroa_p, f = sum, na_pad = T)) %>%
  group_by(Sample) %>% reframe(dup = sum(duplicated(Sample)), No_bees_c, Varroa_p_c) %>%
  filter(dup == 3) %>%
  mutate(C = rep(c('5','53', '532', '5321'), times = 120)) %>% mutate(C = as.factor(C)) %>% 
  pivot_wider(names_from = C, values_from = c(No_bees_c, Varroa_p_c)) %>% 
  left_join(data.full %>% distinct(Sample, .keep_all = T) %>% select(Sample, Site, Status, Density), by = 'Sample') %>%
  left_join(data.growth.w %>% distinct(Site, .keep_all = T) %>% select(Site, Org.farm.w, Ann.fl.w, SNH.w), by = 'Site') %>% 
  drop_na()

separate <- data.full %>% select(Sample, Site, Status, Round, Varroa_p) %>% 
  group_by(Sample, Round) %>% summarise(Varroa_p = mean(Varroa_p)) %>%
  filter(Round != 0) %>% filter(Round != 4) %>% arrange(Sample) %>%
  ## making cumulative sliding window
  pivot_wider(names_from = Round, values_from = Varroa_p, names_prefix = 'Varroa_p_') %>% 
  left_join(data.full %>% distinct(Sample, .keep_all = T) %>% select(Sample, Site, Status, Density), by = 'Sample')  %>%
  left_join(data.growth.w %>% distinct(Site, .keep_all = T) %>% select(Site, Org.farm.w, Ann.fl.w, SNH.w), by = 'Site') %>% 
  drop_na()

