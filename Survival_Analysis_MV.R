library(dplyr)
library(survival)
library(survminer)
mouseMV <- read.csv(file="Mouse_OS_data_MVirus.csv")
s <- Surv(mouseMV$time, mouseMV$status)
survfit(s~1)
survfit(Surv(time, status)~1, data=mouseMV)
sfit <- survfit(Surv(time, status)~1, data=mouseMV)
sfit
summary(sfit)
sfit <- survfit(Surv(time, status)~Group_Responses, data=mouseMV)
sfit
summary(sfit)
plot(sfit)
library(survminer)
ggsurvplot(sfit)

survplots <- list()
survplots[[1]] <- ggsurvplot(sfit, data = mouseMV, conf.int=TRUE, pval=TRUE, risk.table=TRUE, risk.table.col = "strata",
          legend.title="MV Response",  
           palette=c("green", "red", "#E7B800", "#2E9FDF"), 
           title="Kaplan-Meier Curve for Mouse MV Medullablastoma Survival", 
           risk.table.height = 0.25, 
           ggtheme = theme_bw())



survplots[[2]] <- ggsurvplot(
  sfit,                     # survfit object with calculated statistics.
  data = mouseMV,             # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = FALSE,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = c("green", "red", "#E7B800", "#2E9FDF"),
  title="Kaplan-Meier Curve for Mouse MV Medullablastoma Survival", 
  xlim = c(0,150),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in days",   # customize X axis label.
  break.time.by = 10,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
    # change legend labels.
)
# Arrange multiple ggsurvplots and print the output
MV_survPlots <- arrange_ggsurvplots(survplots, print = TRUE,
                    ncol = 2, nrow = 1, risk.table.height = 0.4)

MV_survPlots <- arrange_ggsurvplots(survplots, print = FALSE)
ggsave("MV_survPlots.pdf", width=20, height=10, MV_survPlots)

##############################################################################################

mouseToca <- read.csv(file="Mouse_OS_data_Toca.csv")
s1 <- Surv(mouseToca$time, mouseToca$status)
survfit(s1~1)
survfit(Surv(time, status)~1, data=mouseToca)
sfit1 <- survfit(Surv(time, status)~1, data=mouseToca)
sfit1
summary(sfit1)
sfit1 <- survfit(Surv(time, status)~Group_Responses, data=mouseToca)
sfit1
summary(sfit1)
plot(sfit1)
library(survminer)
ggsurvplot(sfit1)

survplots <- list()
survplots[[1]] <- ggsurvplot(sfit1, data = mouseToca, conf.int=TRUE, pval=TRUE, risk.table=TRUE, risk.table.col = "strata",
                          legend.title="Toca Response",  
                          palette=c("green", "red", "#E7B800", "#2E9FDF"), 
                          title="Kaplan-Meier Curve for Mouse Toca Medullablastoma Survival", 
                          risk.table.height = 0.25, 
                          ggtheme = theme_bw())



survplots[[2]] <- ggsurvplot(
  sfit1,                     # survfit object with calculated statistics.
  data = mouseToca,             # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = FALSE,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = c("green", "red", "#E7B800", "#2E9FDF"),
  title="Kaplan-Meier Curve for Mouse Toca Medullablastoma Survival", 
  xlim = c(0,150),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in days",   # customize X axis label.
  break.time.by = 10,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
      # change legend labels.
)



# Arrange multiple ggsurvplots and print the output
Toca_survPlots <- arrange_ggsurvplots(survplots, print = TRUE,
                                    ncol = 2, nrow = 1, risk.table.height = 0.4)

Toca_survPlots <- arrange_ggsurvplots(survplots, print = FALSE)
ggsave("Toca_survPlots.pdf", width=20, height=10, Toca_survPlots)


library(tidyverse)
library(tidytidbits)
library(survivalAnalysis)
library(gridExtra)
library(ggplot2)
library(meta)
mouseMV <- read.csv(file="Mouse_OS_data_MVirus.csv")
## Plot forest plot

require("survival")
model <- coxph( Surv(time, status) ~ Group_Responses + Tumour_weight,
                data = mouseMV )
ggforest(model)

library("forestmodel")
library("survival")
library("dplyr")
mouseForest <- mouseMV %>%
  transmute(time,
            status,
            Tumour_weight = Tumour_weight,
            Group_Responses = factor(Group_Responses))

plot1 <- print(forest_model(coxph(Surv(time, status) ~ ., mouseForest))) +
  ggtitle("Mouse Measles Virus Survival Forest Plot")
