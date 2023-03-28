library(sf)
library(terra)
library(ggplot2)
library(lme4)
library(MuMIn)

rm(list = ls())
gc()

data = read.csv('./data/occurrence/OD_data_rsf_230323.csv')
herd = read.csv('./data/input_data/herd.csv')
dat <- merge(data, herd, by = "ID")

dat$covid = ifelse(dat$period == 'during',1,0)

#Does habitat use differ between years?
#rsf
rsf_mod <- lmer(rsf ~ period +  (1|herd/ID),
                data = dat,
                na.action = "na.fail",
                REML = F)
rsf_null <- lmer(rsf ~ 1 +  (1|herd/ID),
                 data = dat,
                 na.action = "na.fail",
                 REML = F)

rsf_res <- anova(rsf_mod, rsf_null)
summary(rsf_mod)
dredge(rsf_mod)

paste("p = ",round(rsf_res$`Pr(>Chisq)`[2],2), sep = "")

blue1 = "#4774A0"
blue2 = "#7BA2C9"
blueA = "#BCD5EE"
red = "#FF0000"

#reorder the period levels so they plot in chronological order:
lvl.ordr.x = c('prior2','prior1','during','after')

ggplot(dat) +
  geom_boxplot(aes(x = factor(period, level = lvl.ordr.x), y = rsf, fill = period),  size = 0.2, outlier.size = 0.2) +
  ggtitle("e)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y  = element_text(size=10, family = "sans"),
        plot.title = element_text(size=10, hjust = 0, family = "sans"),
        axis.text.y  = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=8, family = "sans"),
        legend.position="none",
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=8, family = "sans"),
        legend.key.size = unit(0.2, "cm"),
        legend.key = element_blank()) +
  ylab("RSF Score") +
  scale_fill_manual(labels = c("2019", "2020", "2021", "2022"), values = c(blueA, red, blue1,blue2)) +
  scale_x_discrete(labels = c("2019", "2020", "2021", "2022"))


ggsave(file="./figures/RSF_boxplot_230323.png",
       width = 1700,
       height=900,
       units = "px",
       dpi = 300)


#EOF
