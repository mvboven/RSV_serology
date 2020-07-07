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
folder <- "RSV_serology"
PATH <- "/home/andewegs/" #"/home/andewegs/1_RSV_scripts/"
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
rsv.data <- read_csv(file = "data/rsv_serology_csv.txt")#"data/Riskfactors2_csv.txt")
season_border = "10-01" #MM-DD
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
  mutate(
    # Make number of siblings 0-4 years two level factor variable
    Siblings04 = case_when(
      household04_counts <= 1 ~ 'False',
      household04_counts > 1 ~ 'True'
                          )
      ) %>%
  mutate(
    # Make number of siblings 5-9 years two level factor variable
    Siblings59 = case_when(
      household59 <= 0 ~ 'False',
      household59 > 0 ~ 'True'
    )
  ) %>% 
  mutate(
    # Set nursery 0 1 to False True
    Nursery = case_when(
      visitnursery_child == 0 ~ 'False',
      visitnursery_child == 1 ~ 'True'
    ) ) %>% 
#   mutate(
     # Set variable for Pienter2 or Pienter3
 #     PIENTER = case_when(
  #      startsWith(ID, "P2") ~ "P2",
   #     startsWith(ID, "P3") ~ "P3"
  #    )
 # ) %>%
  # Set variable for the different seasons.
  mutate(seasons = case_when(
    consultdate < paste("2006-", season_border, sep = "")  ~ "2005/2006", 
    (consultdate >= paste("2006-", season_border, sep = "") &  consultdate < "2010-01-01") ~ "2006/2007",
    (consultdate >= "2010-01-01" &  consultdate < paste("2016-", season_border, sep = "")) ~ "2015/2016",
    consultdate >= paste("2016-", season_border, sep = "")  ~ "2016/2017"))
# Set the variables to factors
rsv.data$Siblings04 <- factor(rsv.data$Siblings04)
rsv.data$Siblings59 <- factor(rsv.data$Siblings59)
rsv.data$Nursery <- factor(rsv.data$Nursery)
rsv.data$visitnursery_house <- factor(rsv.data$visitnursery_house)
rsv.data$seasons <- factor(rsv.data$seasons)

# Only keep data with the variables of intrest n = 616
completeVec <- complete.cases(rsv.data[, c("Siblings04", "Nursery")])
rsv.data <- rsv.data[completeVec,]

#
# Analysis ----
#

# Model 1 with age and birth doy
model1 <- gam(
  formula = infection ~ 
    breastfeeding +
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
  Birth_doy = seq(from = 1, to = 365, by = 1))

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

# Generate some fake data for the model prediction
rsv.preddata <- expand.grid(
  age_days = seq(from = 1, to = 1095, by = 1),#as.matrix(c(30,60,90,120,150,180,270,365,547,730,900,1095)),
  Birth_doy = seq(from = 1, to = 365, by = 1))

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


ratio_mean = rsv.doy_age.preddata[rsv.doy_age.preddata$Birth_doy == 182,]$p_inf / rsv.doy_age.preddata[rsv.doy_age.preddata$Birth_doy == 1,]$p_inf
ratio_lwr = rsv.doy_age.preddata[rsv.doy_age.preddata$Birth_doy == 182,]$p_inf_lwr / rsv.doy_age.preddata[rsv.doy_age.preddata$Birth_doy == 1,]$p_inf_lwr
ratio_upr = rsv.doy_age.preddata[rsv.doy_age.preddata$Birth_doy == 182,]$p_inf_upr / rsv.doy_age.preddata[rsv.doy_age.preddata$Birth_doy == 1,]$p_inf_upr
age = rsv.doy_age.preddata[rsv.doy_age.preddata$Birth_doy == 182,]$age_days
ratio.data <- data.frame(age_days = age, ratio_mean = ratio_mean, ratio_lwr = ratio_lwr, ratio_upr = ratio_upr)

difference_mean = rsv.doy_age.preddata[rsv.doy_age.preddata$Birth_doy == 182,]$p_inf - rsv.doy_age.preddata[rsv.doy_age.preddata$Birth_doy == 1,]$p_inf
diff.data <- data.frame(age_days = age, diff_mean = difference_mean)
age_ratio <- ggplot(
  data = ratio.data,
  mapping = aes(
    x = age_days, y = ratio_mean)) +
  #geom_ribbon(
  #  alpha = 0.25) +
  geom_line() +
  ylab("Ratio of Pinf (July 1/ Januari 1)") +
  xlab("Age (days)") +
  scale_x_continuous(limits = c(0, 1095), expand = c(0, 0)) +
  scale_y_continuous(limits = c(1, 3),  expand = c(0, 0)) + 
  #coord_cartesian(xlim = c(3, plotmaxage), ylim = c(0, y_max)) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(colour = "grey65", size = 0.2),
        strip.background = element_rect(
          color="black", fill="white", linetype = "blank"),
        strip.text = element_text(size=8)
  )

age_diff <- ggplot(
  data = diff.data,
  mapping = aes(
    x = age_days, y = diff_mean)) +
  #geom_ribbon(
  #  alpha = 0.25) +
  geom_line() +
  ylab("Difference of Pinf (July 1 - Januari 1)") +
  xlab("Age (days)") +
  scale_x_continuous(limits = c(0, 1095), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.5),  expand = c(0, 0)) + 
  #coord_cartesian(xlim = c(3, plotmaxage), ylim = c(0, y_max)) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(colour = "grey65", size = 0.2),
        strip.background = element_rect(
          color="black", fill="white", linetype = "blank"),
        strip.text = element_text(size=8)
  )

grid.arrange(age_ratio, age_diff, ncol=2)
ggsave(file = file.path(figure_folder, "diff_ratio_pinf.svg"), arrangeGrob(age_ratio, age_diff, ncol=2))


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
rsv.doy_age.preddata_eqstep <- bind_cols(rsv.preddata, tmp) %>%
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
  data = rsv.doy_age.preddata_eqstep,
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

## Model 2 with age, birth doy and Siblings04
model2 <- gam(
  formula = infection ~ 
    Siblings04 + 
    ti(Birth_doy, bs = "cp", k = 11) + 
    ti(age_days, bs = "ps", k = 25),# +
    #ti(age_days, bs = "ps", k = 15, by = Siblings04),
  knots = list(birth_doy = c(1, 365)),
  family = binomial,
  method = "REML",  #restricted maximum likelihood 
  data = rsv.data)
summary(model2)

# Generate some fake data for the model prediction
rsv.preddata <- expand.grid(
  age_days = c(30,90,180,270,365,547,730,900,1095),
  Birth_doy = seq(from = 1, to = 365, by = 5),
  Siblings04 = factor(c('True', 'False')))

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
    group = Siblings04, fill = Siblings04, color = Siblings04
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

# Additional model 2, with age, birth doy and Siblings59
model2.1 <- gam(
  formula = infection ~ 
    Siblings59 + 
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
  Siblings59 = factor(c('True', 'False')))

# Make predictions for each combination using the GAM model
# Do this in a temporary tibble, because two columns are produced: fit and se.fit
tmp <- predict(
  object = model2.1,
  newdata = rsv.preddata,
  type = "link",
  se.fit = TRUE) %>%
  as_tibble

# Bind tmp to rsv.preddata and calculate the p_inf including 95% lower and upper bound
rsv.doy_age_sib59.preddata <- bind_cols(rsv.preddata, tmp) %>%
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
  data = rsv.doy_age_sib59.preddata,
  mapping = aes(
    x = Birth_doy, y = p_inf, ymin = p_inf_lwr, ymax = p_inf_upr, 
    group= Siblings59, col= Siblings59, fill = Siblings59)) +
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
  ) + guides(fill=guide_legend(title="Day-care"), col=guide_legend(title="Day-care"))

doy_age_nur   
ggsave(file = file.path(figure_folder, 'age_doy_nur.svg'), plot = doy_age_nur)

# Model4 with age, birth doy, Siblings04 and nursery
model4 <- gam(
  formula = infection ~ 
    Nursery + 
    Siblings04 +
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
  Siblings04 = factor(c('True', 'False')),
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
    facets = ~ age_days + Siblings04, labeller = labeller(age_days = 
                                                            c("30" = "Age: 30 days",
                                                              "90" = "Age: 90 days",
                                                              "180" = "Age: 180 days",
                                                              "270" = "Age: 270 days",
                                                              "365" = "Age: 365 days",
                                                              "547" = "Age: 547 days",
                                                              "730" = "Age: 730 days",
                                                              "900" = "Age: 900 days",
                                                              "1095" = "Age: 1095 days"),
                                                          Siblings04 = c("False" = "Siblings04: False",
                                                                         "True" = "Siblings04: True"), .multi_line = FALSE)
    ) + theme_bw() + 
  theme(panel.grid.minor = element_blank(), #remove grid lines between ticks
        panel.grid.major = element_line(colour = "grey65", size = 0.2), #set size and color of grid lines
        strip.background = element_rect(color="black", fill="white", linetype = "blank"), #set the background
        legend.position = "bottom", #remove legend
        strip.text = element_text(size=7.4) #set text size for text above the plots
  ) + guides(fill=guide_legend(title="Day-care"), col=guide_legend(title="Day-care"))

doy_age_sib_nur
ggsave(file = file.path(figure_folder, 'age_doy_sib_nur.svg'), plot = doy_age_sib_nur)

#The examples/numbers given in the manuscript
# Abstract summer vs winter 
rsv.doy_age.preddata[rsv.doy_age.preddata$age_days == 365 & rsv.doy_age.preddata$Birth_doy == 1,c('p_inf', 'p_inf_lwr', 'p_inf_upr')]
rsv.doy_age.preddata[rsv.doy_age.preddata$age_days == 365 & rsv.doy_age.preddata$Birth_doy == 182,c('p_inf', 'p_inf_lwr', 'p_inf_upr')]
# Abstract siblings04 
rsv.doy_age_sib.preddata[rsv.doy_age_sib.preddata$age_days == 365 & 
                           rsv.doy_age_sib.preddata$Birth_doy == 181 & 
                           rsv.doy_age_sib.preddata$Siblings04 == 'True',c('p_inf', 'p_inf_lwr', 'p_inf_upr')]
rsv.doy_age_sib.preddata[rsv.doy_age_sib.preddata$age_days == 365 & 
                           rsv.doy_age_sib.preddata$Birth_doy == 181 & 
                           rsv.doy_age_sib.preddata$Siblings04 == 'False',c('p_inf', 'p_inf_lwr', 'p_inf_upr')]
# Abstract nursery
rsv.doy_age_nur.preddata[rsv.doy_age_nur.preddata$age_days == 365 & 
                           rsv.doy_age_nur.preddata$Birth_doy == 181 & 
                           rsv.doy_age_nur.preddata$Nursery == 'True',]
rsv.doy_age_nur.preddata[rsv.doy_age_nur.preddata$age_days == 365 & 
                           rsv.doy_age_nur.preddata$Birth_doy == 181 & 
                           rsv.doy_age_nur.preddata$Nursery == 'False',c('p_inf', 'p_inf_lwr', 'p_inf_upr')]
# Results July vs Januari
#half a year
rsv.doy_age.preddata[rsv.doy_age.preddata$age_days == 182 & rsv.doy_age.preddata$Birth_doy == 15,c('p_inf', 'p_inf_lwr', 'p_inf_upr')] #jan. 15
rsv.doy_age.preddata[rsv.doy_age.preddata$age_days == 182 & rsv.doy_age.preddata$Birth_doy == 196,c('p_inf', 'p_inf_lwr', 'p_inf_upr')]# juli 15
# one year
rsv.doy_age.preddata[rsv.doy_age.preddata$age_days == 365 & rsv.doy_age.preddata$Birth_doy == 15,c('p_inf', 'p_inf_lwr', 'p_inf_upr')] #jan. 15
rsv.doy_age.preddata[rsv.doy_age.preddata$age_days == 365 & rsv.doy_age.preddata$Birth_doy == 196,c('p_inf', 'p_inf_lwr', 'p_inf_upr')]# juli 15

# Results siblings0-4
# age half a year
rsv.doy_age_sib.preddata[rsv.doy_age_sib.preddata$age_days == 180 & 
                           rsv.doy_age_sib.preddata$Birth_doy == 196 & 
                           rsv.doy_age_sib.preddata$Siblings04 == 'True',c('p_inf', 'p_inf_lwr', 'p_inf_upr')]
rsv.doy_age_sib.preddata[rsv.doy_age_sib.preddata$age_days == 180 & 
                           rsv.doy_age_sib.preddata$Birth_doy == 196 & 
                           rsv.doy_age_sib.preddata$Siblings04 == 'False',c('p_inf', 'p_inf_lwr', 'p_inf_upr')]
# age 1 year
rsv.doy_age_sib.preddata[rsv.doy_age_sib.preddata$age_days == 365 & 
                           rsv.doy_age_sib.preddata$Birth_doy == 196 & 
                           rsv.doy_age_sib.preddata$Siblings04 == 'True',c('p_inf', 'p_inf_lwr', 'p_inf_upr')]
rsv.doy_age_sib.preddata[rsv.doy_age_sib.preddata$age_days == 365 & 
                           rsv.doy_age_sib.preddata$Birth_doy == 196 & 
                           rsv.doy_age_sib.preddata$Siblings04 == 'False',c('p_inf', 'p_inf_lwr', 'p_inf_upr')]

# Results day-care
# age half a year
rsv.doy_age_nur.preddata[rsv.doy_age_nur.preddata$age_days == 180 & 
                           rsv.doy_age_nur.preddata$Birth_doy == 196 & 
                           rsv.doy_age_nur.preddata$Nursery == 'True',c('p_inf', 'p_inf_lwr', 'p_inf_upr')]
rsv.doy_age_nur.preddata[rsv.doy_age_nur.preddata$age_days == 180 & 
                           rsv.doy_age_nur.preddata$Birth_doy == 196 & 
                           rsv.doy_age_nur.preddata$Nursery == 'False',c('p_inf', 'p_inf_lwr', 'p_inf_upr')]
# age a year
rsv.doy_age_nur.preddata[rsv.doy_age_nur.preddata$age_days == 365 & 
                           rsv.doy_age_nur.preddata$Birth_doy == 196 & 
                           rsv.doy_age_nur.preddata$Nursery == 'True',c('p_inf', 'p_inf_lwr', 'p_inf_upr')]
rsv.doy_age_nur.preddata[rsv.doy_age_nur.preddata$age_days == 365 & 
                           rsv.doy_age_nur.preddata$Birth_doy == 196 & 
                           rsv.doy_age_nur.preddata$Nursery == 'False',c('p_inf', 'p_inf_lwr', 'p_inf_upr')]

#HR vs LR
# at half a year
rsv.doy_age_sib_nur.preddata[rsv.doy_age_sib_nur.preddata$age_days == 180 & 
                               rsv.doy_age_sib_nur.preddata$Birth_doy == 196 & 
                               rsv.doy_age_sib_nur.preddata$Siblings04 == 'True' &
                               rsv.doy_age_sib_nur.preddata$Nursery == 'True',c('p_inf', 'p_inf_lwr', 'p_inf_upr')]

rsv.doy_age_sib_nur.preddata[rsv.doy_age_sib_nur.preddata$age_days == 180 & 
                               rsv.doy_age_sib_nur.preddata$Birth_doy == 196 & 
                               rsv.doy_age_sib_nur.preddata$Siblings04 == 'False' &
                               rsv.doy_age_sib_nur.preddata$Nursery == 'False',c('p_inf', 'p_inf_lwr', 'p_inf_upr')]

# at a year
rsv.doy_age_sib_nur.preddata[rsv.doy_age_sib_nur.preddata$age_days == 365 & 
                               rsv.doy_age_sib_nur.preddata$Birth_doy == 196 & 
                               rsv.doy_age_sib_nur.preddata$Siblings04 == 'True' &
                               rsv.doy_age_sib_nur.preddata$Nursery == 'True',c('p_inf', 'p_inf_lwr', 'p_inf_upr')]

rsv.doy_age_sib_nur.preddata[rsv.doy_age_sib_nur.preddata$age_days == 365 & 
                               rsv.doy_age_sib_nur.preddata$Birth_doy == 196 & 
                               rsv.doy_age_sib_nur.preddata$Siblings04 == 'False' &
                               rsv.doy_age_sib_nur.preddata$Nursery == 'False',c('p_inf', 'p_inf_lwr', 'p_inf_upr')]

#Discussion
#month one sumer
rsv.doy_age.preddata[rsv.doy_age.preddata$age_days == 30 & rsv.doy_age.preddata$Birth_doy == 196, c('p_inf', 'p_inf_lwr', 'p_inf_upr')]
#month 24 sumer
rsv.doy_age.preddata[rsv.doy_age.preddata$age_days == 395 & rsv.doy_age.preddata$Birth_doy == 196, c('p_inf', 'p_inf_lwr', 'p_inf_upr')]
# month one winter
rsv.doy_age.preddata[rsv.doy_age.preddata$age_days == 30 & rsv.doy_age.preddata$Birth_doy == 15, c('p_inf', 'p_inf_lwr', 'p_inf_upr')]
# month 24 winter
rsv.doy_age.preddata[rsv.doy_age.preddata$age_days == 395 & rsv.doy_age.preddata$Birth_doy == 15, c('p_inf', 'p_inf_lwr', 'p_inf_upr')]



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