#
### Generalized Additive Models, Logistic regression for RSV serology.
### Authors: Stijn Andeweg, Jan van de Kassteele and Michiel van Boven
#

# Load packages
library(mgcv)
library(ggplot2)
library(tidyverse)
library(lubridate)
library(gridExtra)
library(grid)

#set WD
folder <- "tensor_splines"
PATH <- "/home/andewegs/1_RSV_scripts/"
setwd(file.path(PATH, folder))
getwd()

# set figure save path
figure_folder <- '/home/andewegs/1_RSV_scripts/tensor_splines/figures/final_models/'

# set colors in figures
color_1 = '#0072B2' #Blue
color_2 = '#E69F00' #Orange

#
# Data ----
#

# Import data 
rsv.data <- read_csv(file = "data/Riskfactors2_csv.txt")

# Some modifications
rsv.data <- rsv.data %>%
  mutate(
    # Get day of the year of birthday (= number between 1 and 365)
    Birth_doy = birthday %>% yday()
        ) %>%
  mutate(
    household04_counts = case_when(
      age_days/365 > 5 ~ (household04 + 1),  #because the child itself is also included in household size
      age_days/365 <= 5 ~ household04       #the household04 of children of age 5 should also include the child itself
                            )
        ) %>%
  select(
    # Drop unused variables
    -household04,
    -household02,
    -household03,
    -n_household,
    -visitnursery_house,
    -pregnancytime,
    -contacttotal,
    -contact04,
    -contact59,
    -sex) %>%
  mutate(
    # Make number of siblings 0-4 years two level factor variable
    Resident04 = case_when(
      household04_counts <= 1 ~ 'False',
      household04_counts > 1 ~ 'True'
                          )
      ) %>%
  mutate(
    # Make number of siblings 5-9 years two level factor variable
    Resident59 = case_when(
      household59 <= 0 ~ 'False',
      household59 > 0 ~ 'True'
    )
  ) %>% 
  mutate(
    # Set nursery 0 1 to False True
    Nursery = case_when(
      visitnursery_child == 0 ~ 'False',
      visitnursery_child == 1 ~ 'True'
    )
  )
# Set the variables to factors
rsv.data$Resident04 <- factor(rsv.data$Resident04)
rsv.data$Resident59 <- factor(rsv.data$Resident59)
rsv.data$Nursery <- factor(rsv.data$Nursery)

# Only keep data with the variables of intrest
completeVec <- complete.cases(rsv.data[, c("Resident04", "Nursery")])
rsv.data <- rsv.data[completeVec,]

#
# Analysis ----
#

# Model 1 with age and birth doy
model1 <- gam(
  formula = infection ~ 
    ti(age_days, bs = "ps", k = 25) + 
    ti(Birth_doy, bs = "cp", k = 11),# +
   #ti(age_days, Birth_doy, bs = c("ps", "cp")),
  knots = list(Birth_doy = c(1, 365)),
  family = binomial,
  method = "REML",  #restricted maximum likelihood 
  data = rsv.data)
summary(model1)

# Generate some fake data for the model prediction
rsv.preddata <- expand.grid(
  age_days = as.matrix(c(30,60,90,120,150,180,270,365,547,730,900,1095)),
  Birth_doy = seq(from = 1, to = 365, by = 5))

# Make predictions for each combination using the GAM model
# Do this in a temporary tibble, because two columns are produced: fit and se.fit
tmp <- predict(
  object = model1,
  newdata = rsv.preddata,
  type = "link",
  se.fit = TRUE) %>%
  as_tibble

# Bind tmp to rsv.preddata and calculate the p_inf including 95% lower and upper bound
rsv.doy_age.preddata <- bind_cols(rsv.preddata, tmp) %>%
  mutate(
    fit_lwr = fit + qnorm(0.025)*se.fit,   
    fit_upr = fit + qnorm(0.975)*se.fit,
    p_inf = fit %>% plogis,                #create vector of probabilities from log odds
    p_inf_lwr = fit_lwr %>% plogis,
    p_inf_upr = fit_upr %>% plogis)

#
# Plot results ----
#

# Plot Figure 3

doy_age <- ggplot(
  data = rsv.doy_age.preddata,
  mapping = aes(
    x = Birth_doy, y = p_inf, ymin = p_inf_lwr, ymax = p_inf_upr)) +
  geom_ribbon(
    alpha = 0.25) +
  geom_line() +
  ylab("Probability of infection") +
  xlab("Birth day of year") +
  facet_wrap(
    facets = ~ age_days, labeller = labeller(age_days = 
                                               c("30" = "Age: 30 days",
                                                 "60" = "Age: 60 days",
                                                 "90" = "Age: 90 days",
                                                 "120" = "Age: 120 days",
                                                 "150" = "Age: 150 days",
                                                 "180" = "Age: 180 days",
                                                 "270" = "Age: 270 days",
                                                 "365" = "Age: 365 days",
                                                 "547" = "Age: 547 days",
                                                 "730" = "Age: 730 days",
                                                 "900" = "Age: 900 days",
                                                 "1095" = "Age: 1095 days"))
  ) + theme_bw() + 
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(colour = "grey65", size = 0.2),
        strip.background = element_rect(
          color="black", fill="white", linetype = "blank"),
        strip.text = element_text(size=8)
  )

doy_age
ggsave(file = file.path(figure_folder, 'age_doy.svg'), plot = doy_age)


# Model1 for age time steps of 73 days
rsv.preddata <- expand.grid(
  age_days = as.matrix(c(73,146,219,292,365,
                         73+365,146+365,219+365,292+365,365+365,
                         73+365+365,146+365+365,219+365+365,292+365+365,365+365+365)),
  Birth_doy = seq(from = 1, to = 365, by = 5))

# Make predictions for each combination using the GAM model
# Do this in a temporary tibble, because two columns are produced: fit and se.fit
tmp <- predict(
  object = model1,
  newdata = rsv.preddata,
  type = "link",
  se.fit = TRUE) %>%
  as_tibble

# Bind tmp to rsv.preddata and calculate the p_inf including 95% lower and upper bound
rsv.doy_age.preddata <- bind_cols(rsv.preddata, tmp) %>%
  mutate(
    fit_lwr = fit + qnorm(0.025)*se.fit,   
    fit_upr = fit + qnorm(0.975)*se.fit,
    p_inf = fit %>% plogis,                #create vector of probabilities from log odds
    p_inf_lwr = fit_lwr %>% plogis,
    p_inf_upr = fit_upr %>% plogis)

#
# Plot results ----
#

# Plot Supplement figure
doy_age_73 <- ggplot(
  data = rsv.doy_age.preddata,
  mapping = aes(
    x = Birth_doy, y = p_inf, ymin = p_inf_lwr, ymax = p_inf_upr)) +
  geom_ribbon(
    alpha = 0.25) +
  geom_line() +
  ylab("Probability of infection") +
  xlab("Birth day of year") +
  facet_wrap(facets = ~ age_days, ncol = 5, labeller = labeller(age_days = 
                                                                  c("73" = "Age: 73 days",
                                                                    "146" = "Age: 146 days",
                                                                    "219" = "Age: 219 days",
                                                                    "292" = "Age: 292 days",
                                                                    "365" = "Age: 365 days",
                                                                    "438" = "Age: 438 days",
                                                                    "511" = "Age: 511 days",
                                                                    "584" = "Age: 584 days",
                                                                    "657" = "Age: 657 days",
                                                                    "730" = "Age: 730 days",
                                                                    "803" = "Age: 803 days",
                                                                    "876" = "Age: 876 days",
                                                                    "949" = "Age: 949 days",
                                                                    "1022" = "Age: 1022 days",
                                                                    "1095" = "Age: 1095 days"))) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(colour = "grey65", size = 0.2),
        strip.background = element_rect(
          color="black", fill="white", linetype = "blank"),
        strip.text = element_text(size=8)
  )

doy_age_73
ggsave(file = file.path(figure_folder, 'age_doy_15grid.svg'), plot = doy_age_73)

## Model 2 with age, birth doy and resident04
model2 <- gam(
  formula = infection ~ 
    Resident04 + 
    ti(Birth_doy, bs = "cp", k = 11) + 
    ti(age_days, bs = "ps", k = 25),# +
    #ti(age_days, bs = "ps", k = 15, by = Resident04),
  knots = list(birth_doy = c(1, 365)),
  family = binomial,
  method = "REML",  #restricted maximum likelihood 
  data = rsv.data)
summary(model2)

# Generate some fake data for the model prediction
rsv.preddata <- expand.grid(
  age_days = c(30,90,180,270,365,547,730,900,1095),
  Birth_doy = seq(from = 1, to = 365, by = 5),
  Resident04 = factor(c('True', 'False')))

# Make predictions for each combination using the GAM model
# Do this in a temporary tibble, because two columns are produced: fit and se.fit
tmp <- predict(
  object = model2,
  newdata = rsv.preddata,
  type = "link",
  se.fit = TRUE) %>%
  as_tibble

# Bind tmp to rsv.preddata and calculate the p_inf including 95% lower and upper bound
rsv.doy_age_sib.preddata <- bind_cols(rsv.preddata, tmp) %>%
  mutate(
    fit_lwr = fit + qnorm(0.025)*se.fit,   
    fit_upr = fit + qnorm(0.975)*se.fit,
    p_inf = fit %>% plogis,                #create vector of probabilities from log odds
    p_inf_lwr = fit_lwr %>% plogis,
    p_inf_upr = fit_upr %>% plogis)

#
# Plot results ----
#

# Plot Supplement Figure
doy_age_sib <- ggplot(
  data = rsv.doy_age_sib.preddata,
  mapping = aes(
    x = Birth_doy, y = p_inf, ymin = p_inf_lwr, ymax = p_inf_upr,
    group = Resident04, fill = Resident04, color = Resident04
    )) +
  geom_ribbon(
    alpha = 0.25, colour = NA) +
  geom_line() +
  scale_fill_manual(values=c(color_1, color_2)) + 
  scale_color_manual(values=c(color_1, color_2)) + 
  ylab("Probability of infection") +
  xlab("Birth day of year") +
  facet_wrap(
    facets = ~ age_days, labeller = labeller(age_days = 
                                               c("30" = "Age: 30 days",
                                                 "90" = "Age: 90 days",
                                                 "180" = "Age: 180 days",
                                                 "270" = "Age: 270 days",
                                                 "365" = "Age: 365 days",
                                                 "547" = "Age: 547 days",
                                                 "730" = "Age: 730 days",
                                                 "900" = "Age: 900 days",
                                                 "1095" = "Age: 1095 days"))
  )  + theme_bw() + 
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(colour = "grey65", size = 0.2),
        strip.background = element_rect(
          color="black", fill="white", linetype = "blank"),
        strip.text = element_text(size=8)
  )

doy_age_sib
ggsave(file = file.path(figure_folder, 'age_doy_sib.svg'), plot = doy_age_sib)

# Additional model 2, with age, birth doy and Resident59
model2.1 <- gam(
  formula = infection ~ 
    Resident59 + 
    ti(Birth_doy, bs = "cp", k = 11) + 
    #ti(household04n, bs = "re", k = 10) +
    ti(age_days, bs = "ps", k = 25),# +
  #ti(age_days, bs = "ps", k = 15, by = siblings04),
  knots = list(birth_doy = c(1, 365)),
  family = binomial,
  method = "REML",  #restricted maximum likelihood 
  data = rsv.data)
summary(model2.1)

# Generate some fake data for the model prediction
rsv.preddata <- expand.grid(
  age_days = c(30,90,180,270,365,547,730,900,1095),
  Birth_doy = seq(from = 1, to = 365, by = 5),
  Resident59 = factor(c('True', 'False')))

# Make predictions for each combination using the GAM model
# Do this in a temporary tibble, because two columns are produced: fit and se.fit
tmp <- predict(
  object = model2.1,
  newdata = rsv.preddata,
  type = "link",
  se.fit = TRUE) %>%
  as_tibble

# Bind tmp to rsv.preddata and calculate the p_inf including 95% lower and upper bound
rsv.doy_age_sib.preddata <- bind_cols(rsv.preddata, tmp) %>%
  mutate(
    fit_lwr = fit + qnorm(0.025)*se.fit,   
    fit_upr = fit + qnorm(0.975)*se.fit,
    p_inf = fit %>% plogis,                #create vector of probabilities from log odds
    p_inf_lwr = fit_lwr %>% plogis,
    p_inf_upr = fit_upr %>% plogis)

#
# Plot results ----
#

# Plot Supplement Figure
doy_age_sib59 <- ggplot(
  data = rsv.doy_age_sib.preddata,
  mapping = aes(
    x = Birth_doy, y = p_inf, ymin = p_inf_lwr, ymax = p_inf_upr, 
    group= Resident59, col= Resident59, fill = Resident59)) +
  geom_ribbon(
    alpha = 0.25, colour = NA) +
  geom_line() +
  scale_fill_manual(values=c("green4", "grey25")) + 
  scale_color_manual(values=c("green4", "grey25")) + 
  ylab("Probability of infection") +
  xlab("Birth day of year") +
  facet_wrap(
    facets = ~ age_days, labeller = labeller(age_days = 
                                               c("30" = "Age: 30 days",
                                                 "90" = "Age: 90 days",
                                                 "180" = "Age: 180 days",
                                                 "270" = "Age: 270 days",
                                                 "365" = "Age: 365 days",
                                                 "547" = "Age: 547 days",
                                                 "730" = "Age: 730 days",
                                                 "900" = "Age: 900 days",
                                                 "1095" = "Age: 1095 days"))
  )  + theme_bw() + 
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(colour = "grey65", size = 0.2),
        strip.background = element_rect(
          color="black", fill="white", linetype = "blank"),
        strip.text = element_text(size=8)
  )

doy_age_sib59
ggsave(file = file.path(figure_folder, 'age_doy_sib59.svg'), plot = doy_age_sib59)

# Model3 with age, birth doy and nursery
model3 <- gam(
  formula = infection ~ 
    Nursery + 
    ti(age_days, bs = "ps", k = 25) + 
    ti(Birth_doy, bs = "cp", k = 11),# +
    #ti(contact04, bs = "re", k = 20) +
    #ti(age_days, bs = "ps", k = 25, by = visitnursery_child, m = 1) + 
    #ti(birth_doy, bs = "cp", k = 10, by = visitnursery_child, m = 1),
   # ti(contact04, bs = "re", k = 20, by = siblings04, m = 1), 

  knots = list(birth_doy = c(1, 365)),
  family = binomial,
  method = "REML",  #restricted maximum likelihood 
  data = rsv.data)
summary(model3)

# Generate some fake data for the model prediction
rsv.preddata <- expand.grid(
  age_days = c(30,90,180,270,365,547,730,900,1095),
  Birth_doy = seq(from = 1, to = 365, by = 5),
  Nursery = factor(c('False', 'True')))

# Make predictions for each combination using the GAM model
# Do this in a temporary tibble, because two columns are produced: fit and se.fit
tmp <- predict(
  object = model3,
  newdata = rsv.preddata,
  type = "link",
  se.fit = TRUE) %>%
  as_tibble

# Bind tmp to rsv.preddata and calculate the p_inf including 95% lower and upper bound
rsv.doy_age_nur.preddata <- bind_cols(rsv.preddata, tmp) %>%
  mutate(
    fit_lwr = fit + qnorm(0.025)*se.fit,   
    fit_upr = fit + qnorm(0.975)*se.fit,
    p_inf = fit %>% plogis,                #create vector of probabilities from log odds
    p_inf_lwr = fit_lwr %>% plogis,
    p_inf_upr = fit_upr %>% plogis)

#
# Plot results ----
#

# Plot Supplement Figure
doy_age_nur <- ggplot(
  data = rsv.doy_age_nur.preddata,
  mapping = aes(
    x = Birth_doy, y = p_inf, ymin = p_inf_lwr, ymax = p_inf_upr, 
    group= Nursery, col= Nursery, fill = Nursery)) +
  geom_ribbon(
    alpha = 0.25, colour = NA) +
  geom_line() +
  scale_fill_manual(values=c(color_1, color_2)) + 
  scale_color_manual(values=c(color_1, color_2)) + 
  ylab("Probability of infection") +
  xlab("Birth day of year") +
  facet_wrap(
    facets = ~ age_days, labeller = labeller(age_days = 
                                               c("30" = "Age: 30 days",
                                                 "90" = "Age: 90 days",
                                                 "180" = "Age: 180 days",
                                                 "270" = "Age: 270 days",
                                                 "365" = "Age: 365 days",
                                                 "547" = "Age: 547 days",
                                                 "730" = "Age: 730 days",
                                                 "900" = "Age: 900 days",
                                                 "1095" = "Age: 1095 days"))
  ) + theme_bw() + 
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(colour = "grey65", size = 0.2),
        strip.background = element_rect(
          color="black", fill="white", linetype = "blank"),
        strip.text = element_text(size=8)
  )

doy_age_nur   
ggsave(file = file.path(figure_folder, 'age_doy_nur.svg'), plot = doy_age_nur)

# Model4 with age, birth doy, resident04 and nursery
model4 <- gam(
  formula = infection ~ 
    Nursery + 
    Resident04 +
    ti(age_days, bs = "ps", k = 25) + 
    ti(Birth_doy, bs = "cp", k = 11),# +
    #ti(contact04, bs = "re", k = 10) +
    #ti(age_days, bs = "ps", k = 25, by = siblings04, m = 1) + 
    #ti(birth_doy, bs = "cp", k = 10, by = siblings04, m = 1) +
   # ti(contact04, bs = "re", k = 10, by = siblings04, m = 1) + 
    #ti(age_days, bs = "ps", k = 25, by = visitnursery_child, m = 1) + 
    #ti(birth_doy, bs = "cp", k = 10, by = visitnursery_child, m = 1),# +
   # ti(contact04, bs = "re", k = 10, by = visitnursery_child, m = 1),# 
  
  knots = list(Birth_doy = c(1, 365)),
  family = binomial,
  method = "REML",  #restricted maximum likelihood 
  data = rsv.data)
summary(model4)

# Generate some fake data for the model prediction
rsv.preddata <- expand.grid(
  age_days = c(30,180,365,547,730,1095),
  Birth_doy = seq(from = 1, to = 365, by = 5),
  Resident04 = factor(c('True', 'False')),
  Nursery = factor(c('True', 'False')))

# Make predictions for each combination using the GAM model
# Do this in a temporary tibble, because two columns are produced: fit and se.fit
tmp <- predict(
  object = model4,
  newdata = rsv.preddata,
  type = "link",
  se.fit = TRUE) %>%
  as_tibble

# Bind tmp to rsv.preddata and calculate the p_inf including 95% lower and upper bound
rsv.doy_age_sib_nur.preddata <- bind_cols(rsv.preddata, tmp) %>%
  mutate(
    fit_lwr = fit + qnorm(0.025)*se.fit,   
    fit_upr = fit + qnorm(0.975)*se.fit,
    p_inf = fit %>% plogis,                #create vector of probabilities from log odds
    p_inf_lwr = fit_lwr %>% plogis,
    p_inf_upr = fit_upr %>% plogis)

#
# Plot results ----
#

# Plot Figure 4
doy_age_sib_nur <- ggplot(
  data = rsv.doy_age_sib_nur.preddata,
  mapping = aes(
    x = Birth_doy, y = p_inf, ymin = p_inf_lwr, ymax = p_inf_upr, 
    group=Nursery , col= Nursery, fill = Nursery)) +
  geom_ribbon(
    alpha = 0.25, colour = NA) +
  geom_line() +
  scale_fill_manual(values=c(color_1, color_2)) + 
  scale_color_manual(values=c(color_1, color_2)) + 
  ylab("Probability of infection") +
  xlab("Birth day of year") +
  facet_wrap(
    facets = ~ age_days + Resident04, labeller = labeller(age_days = 
                                                            c("30" = "Age: 30 days",
                                                              "90" = "Age: 90 days",
                                                              "180" = "Age: 180 days",
                                                              "270" = "Age: 270 days",
                                                              "365" = "Age: 365 days",
                                                              "547" = "Age: 547 days",
                                                              "730" = "Age: 730 days",
                                                              "900" = "Age: 900 days",
                                                              "1095" = "Age: 1095 days"),
                                                          Resident04 = c("False" = "Res04: False",
                                                                         "True" = "Res04: True"), .multi_line = FALSE)
    ) + theme_bw() + 
  theme(panel.grid.minor = element_blank(), #remove grid lines between ticks
        panel.grid.major = element_line(colour = "grey65", size = 0.2), #set size and color of grid lines
        strip.background = element_rect(color="black", fill="white", linetype = "blank"), #set the background
        legend.position = "none", #remove legend
        strip.text = element_text(size=8) #set text size for text above the plots
  )

doy_age_sib_nur
ggsave(file = file.path(figure_folder, 'age_doy_sib_nur.svg'), plot = doy_age_sib_nur)

#
### Model fit measurements ----
#
AIC.data <- data.frame(
  model = c('model1', 'model2', 'model3', 'model4'),
  score = c(model1$aic, model2$aic, model3$aic, model4$aic)
)

BIC.data <- data.frame(
  model = c('model1', 'model2','model3', 'model4'),
  score = c(BIC(model1), BIC(model2), BIC(model3),BIC(model4))
)

r.sq.data <- data.frame(
  model = c('model1', 'model2', 'model3', 'model4'),
  score = c(summary(model1)$r.sq, summary(model2)$r.sq, summary(model3)$r.sq, summary(model4)$r.sq)
)
AIC.data$method <- "AIC"
BIC.data$method <- "BIC"
r.sq.data$method <- "r.sq"
model_scores <- rbind(AIC.data, BIC.data, r.sq.data)