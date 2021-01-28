#######################
### Needed packages ###
#######################
library(tidyverse)
library(mgcv)
library(mgcViz)
library(tsibble)
library(feasts)
library(lubridate)
library(forecast)
library(modelr)
library(reshape2)
library(gam)

###########################################
### download data of interest from NEON ###
###########################################
# 1 - water quality product
if(file.exists("waq_arik.rda")) {
  waq_arik <- readRDS("waq_arik.rda")
} else {
  waq_arik <- neonUtilities::loadByProduct(
    dpID = "DP1.20288.001",
    site = "ARIK",
    startdate = "2018-01",
    enddate = "2019-12",
    package = "expanded",
    token = Sys.getenv("NEON_TOKEN"),
    check.size = FALSE
  )
  saveRDS(waq_arik, "waq_arik.rda")
}

# Take only water quality at the down station
water_quality_arik <- waq_arik$waq_instantaneous %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  filter(horizontal_position == "102") %>%
  select(c(5,19,31,33,45,47,59,61,73,75,87,89,101,103,118)) %>%
  mutate(
    start_date_time = ymd_hms(start_date_time) - second(start_date_time)
  ) %>%
  rename(
    spec_cond = specific_conductance,
    label_spec_cond = specific_cond_final_qf,
    oxygen = dissolved_oxygen,
    label_oxygen = dissolved_oxygen_final_qf,
    oxygen_sat = dissolved_oxygen_saturation,
    label_oxygen_sat = dissolved_oxygen_sat_final_qf,
    ph = p_h,
    label_ph = p_h_final_qf,
    chloro = chlorophyll,
    label_chloro = chlorophyll_final_qf,
    label_turbidity = turbidity_final_qf,
    label_f_dom = f_dom_final_qf
  )

# 2 -  Nitrate in Suface Water product
if(file.exists("nsw_arik.rda")) {
  nsw_arik <- readRDS("nsw_arik.rda")
} else {
  nsw_arik <-  neonUtilities::loadByProduct(
    dpID="DP1.20033.001",
    site="ARIK",
    startdate = "2018-01",
    enddate = "2019-12",
    package="expanded",
    token = Sys.getenv("NEON_TOKEN"),
    check.size = FALSE
  )
  saveRDS(nsw_arik, "nsw_arik.rda")
}

# Extract nitrate at 15 minute intervals
nitrate_arik <- nsw_arik$NSW_15_minute %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  mutate(
    start_date_time = ymd_hms(start_date_time) - second(start_date_time)
  ) %>%
  select(c(start_date_time, surf_water_nitrate_mean, final_qf)) %>%
  rename(
    nitrate_mean = surf_water_nitrate_mean,
    label_nitrate_mean = final_qf
  )

# 3 - Temperature Suface Water product
if(file.exists("temp_arik.rda")) {
  temp_arik <- readRDS("temp_arik.rda")
} else {
  temp_arik <-  neonUtilities::loadByProduct(
    dpID="DP1.20053.001",
    site="ARIK",
    startdate = "2018-01",
    enddate = "2019-12",
    package="expanded",
    token = Sys.getenv("NEON_TOKEN"),
    check.size = FALSE
  )
  saveRDS(temp_arik, "temp_arik.rda")
}

# Extract temperature at 15 minute intervals
temp_arik <- temp_arik$TSW_5min %>%
  as_tibble() %>%
  janitor::clean_names()  %>%
  filter(horizontal_position == "101") %>%
  mutate(start_date_time = ymd_hms(start_date_time) - second(start_date_time)) %>%
  select(c(start_date_time, surf_water_temp_mean, final_qf)) %>%
  rename(
    temp_mean = surf_water_temp_mean,
    label_temp_mean = final_qf
  )


# Join water_quality and nitrate, and clean up
data_arik <- left_join(
  nitrate_arik,
  temp_arik,
  by = "start_date_time")

data_arik <- left_join(
  data_arik,
  water_quality_arik,
  by = "start_date_time") 



data_arik <- data_arik %>%
  mutate(
    spec_cond = if_else(label_spec_cond == 1, NA_real_, spec_cond),
    oxygen = if_else(label_oxygen == 1, NA_real_, oxygen),
    oxygen_sat = if_else(label_oxygen_sat == 1, NA_real_, oxygen_sat),
    turbidity = if_else(label_turbidity == 1, NA_real_, turbidity),
    ph = if_else(label_ph == 1, NA_real_, ph),
    chloro = if_else(label_chloro == 1, NA_real_, chloro),
    f_dom = if_else(label_f_dom == 1, NA_real_, f_dom),
    temp_mean = if_else(label_temp_mean == 1, NA_real_, temp_mean),
    nitrate_mean = if_else(label_nitrate_mean == 1, NA_real_ , nitrate_mean)
  )  %>%
  distinct(start_date_time, .keep_all=TRUE)

#clear impossible data points
data_arik$turbidity = if_else(data_arik$turbidity < 0, NA_real_ , data_arik$turbidity)
data_arik$spec_cond = if_else(data_arik$spec_cond < 100, NA_real_ , data_arik$spec_cond)
data_arik$nitrate_mean = if_else(data_arik$nitrate_mean  > 25, NA_real_ , data_arik$nitrate_mean)
data_arik$nitrate_mean = if_else(data_arik$nitrate_mean <= 0, NA_real_ , data_arik$nitrate_mean)


# Create a time_recod column
data_arik$time_recod <- seq(from = 1, to = nrow(data_arik))

###############################
### SI - data visualisation ###
###############################

data<-data_arik
## Plot data
data %>%
  dplyr::select(start_date_time, nitrate_mean, turbidity, temp_mean, spec_cond, oxygen ) %>%
  pivot_longer(-start_date_time, names_to = "variable") %>%
  ggplot(aes(x = start_date_time , y = value, col = variable)) +
  geom_point() +
  facet_grid(rows = vars(variable),  scales = "free", space="free_x") +
  labs(y = "", x = "Time") +
  scale_x_datetime(date_labels ="%Y-%m") +
  theme(legend.position = "none")

#########################
### Figure 1 - 5 days ###
#########################

data<-data_arik[c(52313:52793),] # select 5 days
g2<-data %>%
  dplyr::select(start_date_time, nitrate_mean, turbidity, temp_mean, spec_cond, oxygen ) %>%
  pivot_longer(-start_date_time, names_to = "variable") %>%
  ggplot(aes(x = start_date_time , y = value, col = variable)) +
  geom_point() +
  facet_grid(rows = vars(variable),  scales = "free", space="free_x") +
  labs(y = "", x = "") +
  scale_x_datetime(date_labels ="%D") +
  theme(legend.position = "none",
        strip.text.x = element_text(size = 20))

g2



################
### Figure 2 ### 
################


data_plot_arik<-cbind(data_arik[,c(2,4,6,8)], log(data_arik$turbidity +1))
colnames(data_plot_arik) <- c("Nitrate","Temperature","Spec. Cond.", "Oxygen C.", "log_Turbidity")
data_plot_arik<-melt(data_plot_arik)
data_plot_arik$site<-rep("Arikaree", length = length(data_plot_arik[,1]))

gg<-ggplot(data = data_plot_arik[which(data_plot_arik$variable == "Nitrate"),], aes(factor(1), value, color = site)) +
  geom_boxplot() +  facet_wrap(~variable,scales="free",ncol=3)+ 
  scale_color_manual(values=c("#440154FF", "#277F8EFF", "#9FDA3AFF"))+
  theme(axis.text.x=element_blank(),
        axis.text.y= element_text(size = 11),
        legend.text=element_text(size=15),
        axis.title.x=element_blank(), 
        axis.title.y = element_blank(),
        legend.position="bottom",
        legend.direction="vertical",
        legend.title = element_blank(),
        strip.text.x = element_text(size = 15))

gg1<-ggplot(data = data_plot_arik[which(data_plot_arik$variable == "Temperature"),], aes(factor(1), value, color = site)) +
  geom_boxplot() +  facet_wrap(~variable,scales="free",ncol=3)+ 
  scale_color_manual(values=c("#440154FF", "#277F8EFF", "#9FDA3AFF"))+
  theme(axis.text.x=element_blank(),
        axis.text.y= element_text(size = 11),
        axis.title.x=element_blank(), 
        axis.title.y = element_blank(),
        legend.position= "none",
        legend.title = element_blank(),
        strip.text.x = element_text(size = 15))

gg2<-ggplot(data = data_plot_arik[which(data_plot_arik$variable == "Spec. Cond."),], aes(factor(1), value, color = site)) +
  geom_boxplot() +  facet_wrap(~variable,scales="free",ncol=3)+ 
  scale_color_manual(values=c("#440154FF", "#277F8EFF", "#9FDA3AFF"))+
  theme(axis.text.x=element_blank(),
        axis.text.y= element_text(size = 11),
        axis.title.x=element_blank(), 
        axis.title.y = element_blank(),
        legend.position= "none",
        legend.title = element_blank(),
        strip.text.x = element_text(size = 15))
gg3<-ggplot(data = data_plot_arik[which(data_plot_arik$variable == "Oxygen C."),], aes(factor(1), value, color = site)) +
  geom_boxplot() +  facet_wrap(~variable,scales="free",ncol=3)+ 
  scale_color_manual(values=c("#440154FF", "#277F8EFF", "#9FDA3AFF"))+
  theme(axis.text.x=element_blank(),
        axis.text.y= element_text(size = 11),
        axis.title.x=element_blank(), 
        axis.title.y = element_blank(),
        legend.position= "none",
        legend.title = element_blank(),
        strip.text.x = element_text(size = 15))
gg4<-ggplot(data = data_plot_arik[which(data_plot_arik$variable == "log_Turbidity"),], aes(factor(1), value, color = site)) +
  geom_boxplot() +  facet_wrap(~variable,scales="free",ncol=3)+ 
  scale_color_manual(values=c("#440154FF", "#277F8EFF", "#9FDA3AFF"))+
  theme(axis.text.x=element_blank(),
        axis.text.y= element_text(size = 11),
        axis.title.x=element_blank(), 
        axis.title.y = element_blank(),
        legend.position= "none",
        legend.title = element_blank(),
        strip.text.x = element_text(size = 15))

gridExtra::grid.arrange(gg, gg1, gg2,gg3,gg4, layout_matrix = rbind(c(1,2,3),c(1,4,5)), widths = c(1.2, 1, 1))


###################
### GAM at ARIK ###
###################

#base model with gam package
model_gam_arik<-gam::gam(nitrate_mean ~ spec_cond + oxygen + turbidity + temp_mean, data = data_arik)

# Stepwise selection of variables
step_model_gam<-gam::step.Gam(model_gam_arik, scope=list("spec_cond"=~1+spec_cond+s(spec_cond,4)+s(spec_cond,6)+s(spec_cond,12),
                                                         "oxygen"=~1+oxygen+s(oxygen,4)+s(oxygen,5)+s(oxygen,6)+s(oxygen,12),
                                                         "turbidity"=~1+turbidity+s(turbidity,4)+s(turbidity,5)+s(turbidity,6)+s(turbidity,12),
                                                         "temp_mean"=~1+temp_mean+s(temp_mean,4)+s(temp_mean,6)+s(temp_mean,12)))
# best model selected 
best_model_gam_arik<-mgcv::gam(nitrate_mean ~  s(spec_cond, k = 12) + s(oxygen, k = 12) +  s(temp_mean, k = 12) + s(log(turbidity + 1 ), k = 12)+ s(time_recod, k = 12),  data = data_arik, na.action=na.exclude)

#fill where we deleted na
data_arik <- data_arik %>%
  mutate(gam_res = residuals(best_model_gam_arik)) %>%
  as_tsibble(index=start_date_time) %>%
  fill_gaps()

gam_prediction<-predict(best_model_gam_arik, data_arik)

data_rmse <- cbind(data_arik$nitrate_mean,gam_prediction)
data_rmse <- data_rmse[complete.cases(data_rmse),]
rmse_arik <- Metrics::rmse(data_rmse[,1], data_rmse[,2])

arik_res <- c(rmse_arik, best_model_gam_arik$gcv.ubre)


# Same model for visualization
best_model_gam_arik<-mgcv::gam(nitrate_mean ~  s(spec_cond, k = 12) + s(oxygen, k = 12) + 
                                 +  s(temp_mean, k = 12) +
                                 s(log(turbidity+1) ,  k = 12)+ s(time_recod, k = 12), data = data_arik)



## RMSE and GCV

res<-round(arik_res,2)
colnames(res)<- c("RMSE", "GCV")
res

################
### Figure 3 ###
################

a <- getViz(best_model_gam_arik)

arik1 <- plot(sm(a, 1)) + l_fitLine(colour = "red") + l_ciLine(mul = 5, colour = "blue", linetype = 2) + xlab("Spec. Cond.") + ylab("")  +
  theme_classic()
arik2 <- plot(sm(a, 2)) + l_fitLine(colour = "red") + l_ciLine(mul = 5, colour = "blue", linetype = 2) + xlab("Oxygen") +ylab("") + ylim(-5,5) + theme_classic()
arik3 <- plot(sm(a, 3)) + l_fitLine(colour = "red") +  l_ciLine(mul = 5, colour = "blue", linetype = 2) +  xlab("Temperature") +ylab("")+ ylim(-5,5) + 
  theme_classic()
arik4 <- plot(sm(a, 4)) + l_fitLine(colour = "red") + l_ciLine(mul = 5, colour = "blue", linetype = 2) + xlab("log(Turbidity)") + ylab("") +
  theme_classic()
arik5 <- plot(sm(a, 5)) + l_fitLine(colour = "red") +  l_ciLine(mul = 5, colour = "blue", linetype = 2) + xlab("Time") + ylab("") +
  theme_classic()


gridPrint(arik1, arik2, arik3, arik4,  arik5, ncol = 1)


############################################
### Calculating importance of variables ###
###########################################

#This is long to run, aproximately 30 minutes !!

unloadNamespace("mgcViz")
unloadNamespace("qgam")
unloadNamespace("gamm4")
unloadNamespace("mgcv")
library("vimp")
library("SuperLearner")
library("gam")

#Create a GAM SuperLearner
sl_gam<-create.Learner("SL.gam", params =list(deg.gam = 12))


dfa <- data_arik[, c("nitrate_mean","temp_mean", "spec_cond", "oxygen", "turbidity","time_recod")]
dfa$log_turbidity <- log(dfa$turbidity +1)
dfa <- dfa[, c("nitrate_mean","temp_mean", "spec_cond", "oxygen", "log_turbidity","time_recod")]
dfa<-na.omit(dfa)

## estimate variable importance with learner gam
imp_gam_tempa <- vimp_rsquared(Y = dfa$nitrate_mean, X = dfa[,-1],
                               indx = 1, run_regression = TRUE, SL.library = "SL.gam_1")
imp_gam_conda <- vimp_rsquared(Y = dfa$nitrate_mean, X = dfa[,-1],
                               indx = 2, run_regression = TRUE, SL.library = "SL.gam_1")
imp_gam_oxygena <- vimp_rsquared(Y = dfa$nitrate_mean, X = dfa[,-1],
                                 indx = 3, run_regression = TRUE, SL.library = "SL.gam_1" )
imp_gam_turbiditya <- vimp_rsquared(Y = dfa$nitrate_mean, X = dfa[,-1],
                                    indx = 4, run_regression = TRUE, SL.library = "SL.gam_1")
imp_gam_timea <- vimp_rsquared(Y = dfa$nitrate_mean, X = dfa[,-1],
                               indx = 5, run_regression = TRUE, SL.library = "SL.gam_1")

imp_gam_totala <- vimp_rsquared(Y = dfa$nitrate_mean, X = dfa[,-1],
                                indx = c(1:5), run_regression = TRUE, SL.library = "SL.gam_1")


################
### Figure 4 ###
################

estsa <- merge_vim(imp_gam_tempa, imp_gam_conda, imp_gam_oxygena, imp_gam_turbiditya, imp_gam_timea)
arik_imp<-cbind(estsa$mat[1],(estsa$mat[2:5]/imp_gam_totala$est)*100)
arik_imp$Site <- "Arikaree"


get_nm <- function(s) {
  if (s == 1) {
    return("Temperature")
  } else if (s == 2) {
    return("Conductance")
  } else if (s == 3) {
    return("Oxygen")
  } else if (s == 4) {
    return("Log Turbidity")
  } else if (s == 5) {
    return("Time")
  } 
}
get_nms <- function(ests) {
  return(apply(matrix(as.numeric(ests$s)), 1, get_nm))
}



total_imp %>%
  mutate(ord_group = get_nms(total_imp)) %>%
  ggplot(aes(x = est, y = ord_group, color = Site)) +
  geom_point(position = position_dodge(0.5)) +
  geom_errorbarh(aes(xmin = cil, xmax = ciu),position = position_dodge(0.5)) +
  xlab(expression(paste("Variable importance (%)")))+
  ylab("") +
  ggtitle("")+ 
  scale_color_manual(values=c("#440154FF", "#277F8EFF", "#9FDA3AFF"))+
  theme(legend.text = element_text(size = 14 ),
        legend.title = element_blank(),
        legend.position = c(0.8, 0.4),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        axis.line = element_line(color = "grey") )

