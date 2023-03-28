#script to run glmm on caribou_results_no_error.csv which has been assigned spatial covariates
#and renamed caribou_results_with_covariates2.csv
#to get here, only the 95% home range estimates for each period were selected (exclude upper and lower CI).
#the script 7.extract.raster.95hr.estimate.r was used to extract covariates for each HR estimate.

rm(list = ls())
gc()

#required libraries
library(lme4)
#for p values:
library(lmtest)
library(dplyr)
library(MuMIn)
library(sjPlot)
library(ggplot2)
library(gridExtra)

getwd()

#read in the dataset
data = read.csv('./data/home_range/Caribou_results_with_error_covariates_230214.csv')

plot(data$HR)
plot(data$diffusion)
#there are some NAs in this file, but do not need to remove all, some animals can be
#used for different analyses.

# Some carpentry to force the during period to be the baseline against
# which Before 1 and Before 2 are compared
# Because factors go by alphabetical order in the model
#change 'after' to 'return' to keep the temporal period alphabetical:
data[which(data$period == "after"),"period"] <- "return"

#Convert the other variables to factors
data$period <-  factor(data$period)
data$ID <- as.factor(data$ID)
data$herd <- as.factor(data$herd)

#Convert to HR size to km^2
data$HR <- data$HR/1e+6
data$HR_min <- data$HR_min/1e+6
data$HR_max <- data$HR_max/1e+6

#subset for records used for HR analysis to get rid of NA and remove
#diffusion columns so we can use as much data as possible:
dat.hr = droplevels(data[!(data$ID == 44095),-c(10:15, 22:24)])
dat.hr = na.omit(dat.hr)
plot(dat.hr$HR)
#for looking at skier days:
#period = c('prior2', 'prior1', 'during', 'return')
#skier.days = c(120000,90000,12500,105000)

#dat.heli = data.frame(period, skier.days, stringsAsFactors = TRUE)
#rm(period)

#sample size:
hr.ss = dat.hr %>%
  dplyr::group_by(herd) %>%
  dplyr::summarise(n_animals = n_distinct(ID))

#join skier days to dat.hr
#dat.hr = merge(dat.hr, dat.heli, by = 'period')

#corrplot to determine variables to include in full model
library(corrplot)
dat.cor = dat.hr[,c("HR","elev_2023","slope_2023","roughness_2023","proj_age_numeric_2023")]
names(dat.cor) = c("Home Range","Elevation","Slope","Roughness","Forest Age")
corrplot(cor(dat.cor), method = 'number')

#rescale variables
#dat.hr$elev.scl = scale(dat.hr$elev_2023, center = TRUE, scale = TRUE)
##dat.hr$slope.scl = scale(dat.hr$slope_2023, center = TRUE, scale = TRUE)
#dat.hr$age.scl = scale(dat.hr$proj_age_numeric_2023, center = TRUE, scale = TRUE)


#Build the full model with uncorrelated covariates
m1 <- glmer(HR ~ period + heli_ten_220809 + elev_220809 + slope_220809 + proj_age_220809 + 
              (1|herd/ID),
            family = Gamma('log'),
            data = dat.hr,
            na.action = "na.fail")

summary(m1)

sink("./data/home_range/m1_HR.txt")
print(summary(m1))
sink() 

#AICc based model selection
sel <- dredge(m1); sel
write.csv(sel, "./data/home_range/selmod_HR_230217.csv", row.names = FALSE)


#Selected model (same as full model):
sel_mod <- glmer(HR ~ period + heli_ten_220809 + elev_220809 + slope_220809 + proj_age_220809 + 
                   (1|herd/ID),
                 family = Gamma('log'),
                 data = dat.hr,
                 na.action = "na.fail")

summary(sel_mod)


plot_model(sel_mod,
           type = "pred",
           show.data = T,
           terms = c("heli_ten_220809", "period"))


ggplot(data = data) +
  geom_boxplot(aes(x = period, y = heli_ten_220809))

#Generate a figure of the modelling results
blue1 = "#4774A0"
blue2 = "#7BA2C9"
blueA = "#BCD5EE"
red = "#FF0000"

ggplot() +
  geom_point(data = dat.hr, aes(x = herd, y = HR, color = period)) +
  scale_colour_manual(values = c(red, blue1, blue2, blueA),
                      labels = c("During Lockdown", "1 Year Prior", "2 Years Prior", "After"),
                      name = "period")


#get predicted values for panel e to order the axis:
per.table = get_model_data(sel_mod, type = c('pred'))[[1]]
per.table$period = as.factor(ifelse(per.table$x == 1, '2021',
                                    ifelse(per.table$x == 2, '2020',
                                           ifelse(per.table$x == 3, '2019','2022'))))
per.table = per.table[,c(8,2,4,5)]

#generate plots:
a <- 
  plot_model(sel_mod,
             show.intercept = FALSE,
             sort.est = TRUE,
             transform = NULL,
             value.offset = .4,
             dot.size = 1,
             line.size = 0.2,
             value.size = 2.8,
             colors = "black",
             show.values = TRUE,
             axis.labels = c("Elevation",
                             "Slope",
                             "2019",
                             "2022",
                             "2020",
                             "Forest age",
                             "Heli tenure")) +
  ggtitle("a)") +
  geom_hline(yintercept = 0, linetype = "dashed", col = "grey70", alpha = 0.8) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        #plot.margin = margin(2, 2, 2, 2, "cm"),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=10, family = "sans"),
        axis.title.x = element_text(size=10, family = "sans"),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=8, family = "sans"),
        legend.text  = element_text(size=6, family = "sans"),
        legend.title  = element_text(size=6, family = "sans"),
        legend.key.size = unit(0.35, 'cm'),
        legend.key.height = unit(0.35, 'cm'),
        plot.title = element_text(hjust = -0.05, size = 10, family = "sans"),
        legend.position = c(0.8,0.85)) +
  scale_y_continuous(limits = c(-3,3), expand = c(0,0.01), breaks = c(-3,-2,-1,0,1,2,3)) + 
  #scale_x_discrete(limits = 0,0.1) +
  ylab(expression(beta))

b <- 
  ggplot(per.table, aes(x=period, y=predicted)) + 
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width = 0.07) +
  geom_line() +
  geom_point() +
  ggtitle("b)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=10, family = "sans"),
        axis.title.x = element_text(size=10, family = "sans"),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=8, family = "sans"),
        legend.text  = element_text(size=6, family = "sans"),
        legend.title  = element_text(size=6, family = "sans"),
        legend.key.size = unit(0.35, 'cm'),
        legend.key.height = unit(0.35, 'cm'),
        plot.title = element_text(hjust = -0.05, size = 10, family = "sans"),
        legend.position = c(0.5,0.85)) +
  xlab("Time Period (year)") +
  ylab(expression(paste("Predicted home range area (",km^2,")")))

#for presentation:
#plot_model(sel_mod,
#           show.intercept = TRUE,
#           sort.est = TRUE,
#           transform = NULL,
#           value.offset = .3,
#           dot.size = 3,
#           line.size = 1,
#           title = 'Home Range',
#           value.size = 4,
#           colors = "black",
#           show.values = TRUE,
#           axis.labels = c("2019",
#                           "2022",
#                           "2020",
#                           "Forest age",
#                           "Slope",
#                           "Heli tenure",
#                           "Elevation",
#                           "Intercept")) +
#  geom_hline(yintercept = 0, linetype = "dashed", col = "grey70", alpha = 0.8) +
#  theme_bw() +
#  theme(panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        axis.title.y = element_text(size=16, family = "sans"),
#        axis.title.x = element_text(size=16, family = "sans"),
#        axis.text.y = element_text(size=14, family = "sans", face = 'bold'),
#        axis.text.x  = element_text(size=14, family = "sans", face = 'bold'),
#        #plot.title = element_text(hjust = -0.05, size = 12, family = "sans"),
#        legend.position = c(0.8,0.85)) +
#  scale_y_continuous(limits = c(-3,9), expand = c(0,0.001), breaks = c(-2,0,2,4,6,8,10)) + 
#  scale_x_discrete(expand = c(0,1)) +
#  ylab(expression(beta))


c <- 
  plot_model(sel_mod,
             type = "pred",
             line.size = 0.2,
             terms = c("slope_220809")) +
  geom_point(data = dat.hr, aes(x = slope_220809, y = HR, col = period), alpha = 0.6, pch = 16, size = 0.5) +
  scale_colour_manual(values = c(red, blue1, blue2, blueA),
                      labels = c("2021", "2020", "2019", "2022"),
                      name = "period") + 
  ggtitle("c)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=8, family = "sans"),
        axis.title.x = element_text(size=8, family = "sans"),
        axis.text.y = element_text(size=6, family = "sans"),
        axis.text.x  = element_text(size=6, family = "sans"),
        legend.text  = element_text(size=6, family = "sans"),
        legend.title  = element_text(size=6, family = "sans"),
        legend.key.size = unit(0.35, 'cm'),
        legend.key.height = unit(0.35, 'cm'),
        plot.title = element_text(hjust = -0.05, size = 10, family = "sans"),
        legend.position = "none") + 
  #scale_x_continuous(limits = c(0,1), expand = c(0,0.02)) + 
  #scale_y_continuous(limits = c(0,2100), expand = c(0,0.1)) + 
  ylab(expression(paste("Predicted home range area (",km^2,")"))) +
  xlab("Standardised slope") +
  coord_cartesian(ylim =  c(0,2100))

d <- plot_model(sel_mod,
                type = "pred",
                line.size = 0.2,
                terms = c("elev_220809")) +
  geom_point(data = dat.hr, aes(x = elev_220809, y = HR, col = period), alpha = 0.6, pch = 16, size = 0.5) +
  #scale_x_discrete(limits = c('2021','2020','2019','2022')) +
  scale_colour_manual(values = c(red, blue1, blue2, blueA),
                      labels = c("2021", "2020", "2019", "2022"),
                      name = "period") + 
  ggtitle("d)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=8, family = "sans"),
        axis.title.x = element_text(size=8, family = "sans"),
        axis.text.y = element_text(size=6, family = "sans"),
        axis.text.x  = element_text(size=6, family = "sans"),
        legend.text  = element_text(size=6, family = "sans"),
        legend.title  = element_text(size=6, family = "sans"),
        legend.key.size = unit(0.35, 'cm'),
        legend.key.height = unit(0.35, 'cm'),
        plot.title = element_text(hjust = -0.05, size = 10, family = "sans"),
        legend.position = "none") + 
  #scale_x_continuous(limits = c(790,2115), expand = c(0,0.02)) + 
  # scale_y_continuous(limits = c(300,200), expand = c(0,0.1)) + 
  ylab(expression(paste("Predicted home range area (",km^2,")"))) +
  xlab("Standardised elevation") +
  coord_cartesian(ylim =  c(0,2100))

e <-
  plot_model(sel_mod,
             type = "pred",
             terms = c("proj_age_220809"),
             pred.type = "re",
             title = "",
             line.size = 0.2,
             axis.title = c("Standardised forest age", "Home range area (km^2)")) +
  ggtitle("e)") +
  geom_point(data = dat.hr, aes(x = proj_age_220809, y = HR, col = period), alpha = 0.6, pch = 16, size = 0.5) +
  scale_colour_manual(values = c(red, blue1, blue2, blueA),
                      labels = c("During Lockdown", "1 Year Prior", "2 Years Prior", "After"),
                      name = "period") + 
  scale_fill_manual(values = c(red, blue1, blue2, blueA),
                    labels = c("During Lockdown", "1 Year Prior", "2 Years Prior", "After")) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=8, family = "sans"),
        axis.title.x = element_text(size=8, family = "sans"),
        axis.text.y = element_text(size=6, family = "sans"),
        axis.text.x  = element_text(size=6, family = "sans"),
        legend.text  = element_text(size=6, family = "sans"),
        legend.title  = element_text(size=6, family = "sans"),
        legend.key.size = unit(0.35, 'cm'),
        legend.key.height = unit(0.35, 'cm'),
        plot.title = element_text(hjust = -0.05, size = 10, family = "sans"),
        legend.position = "none") +
  #scale_x_continuous(limits = c(-3.01,4.1), expand = c(0,0.001)) + 
  scale_y_continuous(expand = c(0,0.1)) + 
  ylab(expression(paste("Predicted home range area (",km^2,")"))) +
    coord_cartesian(ylim =  c(0,2100))

f <-
  plot_model(sel_mod,
             type = "pred",
             line.size = 0.2,
             terms = c("heli_ten_220809")) +
  geom_point(data = dat.hr, aes(x = heli_ten_220809, y = HR, col = period), alpha = 0.6, pch = 16, size = 0.5) +
  #scale_x_discrete(limits = c('2021','2020','2019','2022')) +
  scale_colour_manual(values = c(red, blue1, blue2, blueA),
                      labels = c("2021", "2020", "2019", "2022"),
                      name = "period") + 
  ggtitle("f)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=8, family = "sans"),
        axis.title.x = element_text(size=8, family = "sans"),
        axis.text.y = element_text(size=6, family = "sans"),
        axis.text.x  = element_text(size=6, family = "sans"),
        legend.text  = element_text(size=6, family = "sans"),
        legend.title  = element_text(size=6, family = "sans"),
        legend.key.size = unit(0.35, 'cm'),
        legend.key.height = unit(0.35, 'cm'),
        plot.title = element_text(hjust = -0.05, size = 10, family = "sans"),
        legend.position = "none") + 
  scale_x_continuous(limits = c(0,1), expand = c(0,0.02)) + 
 # scale_y_continuous(limits = c(0,200), expand = c(0,0.1)) + 
  ylab(expression(paste("Predicted home range area (",km^2,")"))) +
  xlab("Mean Heli-ski tenure score") +
  coord_cartesian(ylim =  c(0,2100))

top = grid.arrange(a,b, ncol = 2)
bot <- grid.arrange(c,d,e,f, ncol = 4)
Fig <- grid.arrange(top,bot, ncol = 1, heights = c(3,2))

#using px for better scaling:
ggsave(Fig,
       width = 2200, height = 1400, units = "px", dpi = 300,
       file="./figures/HR_regression_ALL_230217_PX.png")

#-------------------------------------------------------------------------------
# skier days:
set_theme(geom.alpha = 0) #I just want to plot the smoothed line, so make the plot_model line transparent
plot_model(sel_mod,
           type = "pred",
           ci.lvl = 0.95,
           terms = c("skier.days")) +
  stat_smooth(method = "loess", 
              span = 1.4, 
              se = FALSE,
              color = 'black') +
  scale_x_discrete(limits = c('2021','2020','2019','2022')) +
  scale_colour_manual(values = c(red, blue1, blue2, blueA),
                      labels = c("2021", "2020", "2019", "2022"),
                      name = "period") + 
  ggtitle("g)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=8, family = "sans"),
        axis.title.x = element_text(size=8, family = "sans"),
        axis.text.y = element_text(size=6, family = "sans"),
        axis.text.x  = element_text(size=6, family = "sans"),
        legend.text  = element_text(size=6, family = "sans"),
        legend.title  = element_text(size=6, family = "sans"),
        legend.key.size = unit(0.35, 'cm'),
        legend.key.height = unit(0.35, 'cm'),
        plot.title = element_text(hjust = -0.05, size = 10, family = "sans"),
        legend.position = "none") + 
  scale_x_continuous(limits = c(12000,122000), expand = c(0,0.02)) + 
  # scale_y_continuous(limits = c(0,200), expand = c(0,0.1)) + 
  ylab(expression(paste("Predicted home range area (",km^2,")"))) +
  xlab("Helicat Canada Skier Days")
ggsave(width = 6, height = 4, units = "in", dpi = 600,
       file="./UBC/thesis/figures/2022/HR_skier_days.jpg")

#-------------------------------------------------------------------------------

#  p1 individual boxplot proportion ordered by proportion
dat.hr$ID = as.factor(dat.hr$ID) 

#how many caribou with >90% inside tenure?
n.bou = dat.hr %>%
  dplyr::group_by(period) %>%
  dplyr::summarise(total_bou = n_distinct(ID), 
                   gt90 = sum(heli_ten_220809 >= 0.95),
                   p_gt90 = round(gt90 / total_bou,2))


j = ggplot() +
  geom_boxplot(data = dat.hr[dat.hr$period == 'prior2',], aes(x = reorder(ID, heli_ten_220809, sum), y = heli_ten_220809)) +
  #geom_hline(yintercept = 0.95, linetype = 'dashed', colour = 'red') +
  ggtitle(paste('2019',n.bou[n.bou$period == 'prior2',c(4)], sep = ' - ')) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=8, family = "sans"),
        axis.title.x = element_text(size=8, family = "sans"),
        axis.text.y = element_text(size=6, family = "sans"),
        axis.text.x  = element_blank(), #element_text(size=6, family = "sans", angle = 90),
        legend.text  = element_text(size=6, family = "sans"),
        legend.title  = element_text(size=6, family = "sans"),
        legend.key.size = unit(0.35, 'cm'),
        legend.key.height = unit(0.35, 'cm'),
        plot.title = element_text(hjust = -0.05, size = 10, family = "sans"),
        legend.position = "none") + 
  ylab(expression(paste("Proportion of HR in heli-tenure (",km^2,")"))) +
  xlab("Individual caribou")

k = ggplot() +
  geom_boxplot(data = dat.hr[dat.hr$period == 'prior1',], aes(x = reorder(ID, heli_ten_220809, sum), y = heli_ten_220809)) +
  #geom_hline(yintercept = 0.95, linetype = 'dashed', colour = 'red') +
  ggtitle(paste('2020',n.bou[n.bou$period == 'prior1',c(4)], sep = ' - ')) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=8, family = "sans"),
        axis.title.x = element_text(size=8, family = "sans"),
        axis.text.y = element_text(size=6, family = "sans"),
        axis.text.x  = element_blank(), #element_text(size=6, family = "sans", angle = 90),
        legend.text  = element_text(size=6, family = "sans"),
        legend.title  = element_text(size=6, family = "sans"),
        legend.key.size = unit(0.35, 'cm'),
        legend.key.height = unit(0.35, 'cm'),
        plot.title = element_text(hjust = -0.05, size = 10, family = "sans"),
        legend.position = "none") + 
  ylab(expression(paste("Proportion of HR in heli-tenure (",km^2,")"))) +
  xlab("Individual caribou")

l = ggplot() +
  geom_boxplot(data = dat.hr[dat.hr$period == 'during',], aes(x = reorder(ID, heli_ten_220809, sum), y = heli_ten_220809)) +
  #geom_hline(yintercept = 0.95, linetype = 'dashed', colour = 'red') +
  ggtitle(paste('2021',n.bou[n.bou$period == 'during',c(4)], sep = ' - ')) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=8, family = "sans"),
        axis.title.x = element_text(size=8, family = "sans"),
        axis.text.y = element_text(size=6, family = "sans"),
        axis.text.x  = element_blank(), #element_text(size=6, family = "sans", angle = 90),
        legend.text  = element_text(size=6, family = "sans"),
        legend.title  = element_text(size=6, family = "sans"),
        legend.key.size = unit(0.35, 'cm'),
        legend.key.height = unit(0.35, 'cm'),
        plot.title = element_text(hjust = -0.05, size = 10, family = "sans"),
        legend.position = "none") + 
  ylab(expression(paste("Proportion of HR in heli-tenure (",km^2,")"))) +
  xlab("Individual caribou")

m = ggplot() +
  geom_boxplot(data = dat.hr[dat.hr$period == 'return',], aes(x = reorder(ID, heli_ten_220809, sum), y = heli_ten_220809)) +
  #geom_hline(yintercept = 0.95, linetype = 'dashed', colour = 'red') +
  ggtitle(paste('2022',n.bou[n.bou$period == 'return',c(4)], sep = ' - ')) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=8, family = "sans"),
        axis.title.x = element_text(size=8, family = "sans"),
        axis.text.y = element_text(size=6, family = "sans"),
        axis.text.x  = element_blank(), #element_text(size=6, family = "sans", angle = 90),
        legend.text  = element_text(size=6, family = "sans"),
        legend.title  = element_text(size=6, family = "sans"),
        legend.key.size = unit(0.35, 'cm'),
        legend.key.height = unit(0.35, 'cm'),
        plot.title = element_text(hjust = -0.05, size = 10, family = "sans"),
        legend.position = "none") + 
  ylab(expression(paste("Proportion of HR in heli-tenure (",km^2,")"))) +
  xlab("Individual caribou")

Figa <- grid.arrange(j,k,l,m, ncol = 2, heights = c(2,2))
ggsave(Figa,
       width = 2000, height = 1400, units = "px", dpi = 300,
       file="./figures/HR_proportion_tenure_2300217_PX.png")
