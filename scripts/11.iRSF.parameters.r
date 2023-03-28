
library(metafor)
library(ggplot2)

rm(list = ls())
gc()

data <- read.csv('./data/rsf/all_periods_RSF_metafor_230301.csv')

periods <- unique(data$period)
params <- c("slope", "elev", "heli.ten", "proj.age")

period_res <- list()

#loop over the periods (indexed by i)
for(i in 1:length(periods)){
  
  subset <- data[which(data$period == periods[i]),]
  
  param_results <- list()
  # loop over the coefficients (indexed by j)
  for(j in 1:length(params)){
    
    param_subset <- data.frame(yi = subset[,params[j]],
                               vi = subset[,paste(params[j],".cov",sep = "")])
    
    
    res <- rma(yi, vi,
               mods = ~ 1,
               data=param_subset,
               method = "REML",
               control=list(maxiter=10000, stepadj=0.5))
    res
    
    preds <- predict(res)
    preds[c("pred", "ci.lb", "ci.lb")]
    
    results <- data.frame(period = periods[i],
                          param = params[j],
                          beta_hat = preds$pred,
                          beta_hat_min = preds$ci.lb,
                          beta_hat_max = preds$ci.ub)
    
    param_results[[j]] <- results
  }
  
  period_res[[i]] <- do.call(rbind, param_results)
  
}

results <- do.call(rbind, period_res)

blue1 = "#4774A0"
blue2 = "#7BA2C9"
blueA = "#BCD5EE"
red = "#FF0000"

pd = position_dodge(width = 0.5)
ggplot(results, position = pd) +
  geom_vline(xintercept = 0, col = "grey70", linetype = "dashed") +
  geom_pointrange(aes(x = beta_hat, y = param, xmin = beta_hat_min, xmax = beta_hat_max, color = period), size = 0.5, position = pd) +
  theme_bw() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30", size = 0.5, alpha = 0.7) +
  scale_y_discrete(name = '', 
                   labels = c('Elevation',
                              'Heli Tenure', 
                              'Forest Age', 
                              'Slope')) +
  scale_x_continuous(name = 'Parameter Estimate') +
  guides(color = guide_legend(reverse = TRUE)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="top",
        axis.title.y = element_text(size=12, family = "sans"),
        axis.title.x = element_text(size=12, family = "sans"),
        axis.text.y = element_text(size=10, family = "sans", vjust = 0.4),
        axis.text.x  = element_text(size=10, family = "sans")) +
  scale_colour_manual(values = c(blueA, red, blue1, blue2), labels = c('2022', '2021', '2020', '2019')) +
  coord_cartesian(xlim = c(-.5, 0.5))
ggsave(file = "./UBC/thesis/figures/2022/rsf/PARAMETER_ESTIMATES_221221.jpg", units = 'cm', width = 17, height = 9)


