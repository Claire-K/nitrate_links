#' @title Data acquisition, processing, and modeling to understand links between
#'  water quality and nitrate 
#' @author Claire Kermorvant, Guy Litt
#' @description Companion code to "Understanding links between water-quality 
#' variables and nitrate concentration in freshwater streams using 
#' high-frequency sensor data" https://arxiv.org/abs/2106.01719
#' @details Requires using neonUtilities >= version 2.2.1
# Changelog/contributions
# 2021-05-21 Originally created, Claire Kermorvant
# 2023-04-02 Update data acquisition & results display, generalized to multiple sites,
#    change temp horizontal position from 101 (upstream) to 102 (downstream) Guy Litt
#     Added temporal aggregation method (mean/max)

#TODO substitute end_date_time for start_date_time
#TODO describe SWE as surface water elevation
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
library(neonUtilities) # Must be >= version 2.2.1 
library(dplyr)
library(lubridate)
library(padr)
library(svglite)
# library(devtools)
# 
# install package in github used for custom plotting
devtools::install_github(repo = "https://github.com/remichel/rmTools.git")
###########################################
### download data of interest from NEON ###
###########################################
neon_sites <-  c("ARIK","LEWI","CARI")

# Define the file saving directory
user_base_dir <- Sys.getenv("USERPROFILE")
save_dir <- file.path(user_base_dir,"Documents","waq_nitrate_analysis")
plot_dir <- file.path(save_dir, "plots")

# Create the save directory
if(!dir.exists(plot_dir)){
  print(paste0("Creating the following directory path: ", plot_dir))
  dir.create(plot_dir,recursive=TRUE)
}

neonPkgVer <- utils::packageVersion("neonUtilities")
if(neonPkgVer < base::package_version('2.2.1')){
  warning("The neonUtilities package version should at least be >= 2.2.1")
}

# Define the instantaneous water quality data columns of interest
waq_cols_sel <- c("startDateTime","horizontalPosition",  "specificConductance","specificCondFinalQF",
                  "dissolvedOxygen","dissolvedOxygenFinalQF","dissolvedOxygenSaturation",
                  "dissolvedOxygenSatFinalQF","pH","pHFinalQF",               
                  "chlorophyll","chlorophyllFinalQF","turbidity",                
                  "turbidityFinalQF","fDOM","fDOMFinalQF")


# ---------------------------------------------------------------------------- #
#   Download data
# ---------------------------------------------------------------------------- #
# 1 - water quality product
if(file.exists(file.path(save_dir,"waq.rda"))) {
  waq <- readRDS(file.path(save_dir,"waq.rda"))
} else {
  print("Downloading water quality data (may take a while)")
  waq <- neonUtilities::loadByProduct(
    dpID = "DP1.20288.001",
    site = neon_sites,
    startdate = "2018-01",
    enddate = "2019-12",
    package = "expanded",
    token = Sys.getenv("NEON_TOKEN"),
    check.size = FALSE,
    release="RELEASE-2021"
  )
  saveRDS(waq, file.path(save_dir,"waq.rda"))
}


# 2 -  Nitrate in Suface Water data product
if(file.exists(file.path(save_dir,"nsw.rda"))) {
  nsw <- readRDS(file.path(save_dir,"nsw.rda"))
} else {
  print("Downloading nitrate in surface water data")
  nsw <-  neonUtilities::loadByProduct(
    dpID="DP1.20033.001",
    site=neon_sites,
    startdate = "2018-01",
    enddate = "2019-12",
    package="expanded",
    token = Sys.getenv("NEON_TOKEN"),
    check.size = FALSE
  )
  saveRDS(nsw, file.path(save_dir,"nsw.rda"))
}

# 3 - Temperature Suface Water product
if(file.exists(file.path(save_dir,"temp_all.rda"))) {
  temp_all <- readRDS(file.path(save_dir,"temp_all.rda"))
} else {
  print("Downloading surface water temperature data")
  temp_all <-  neonUtilities::loadByProduct(
    dpID="DP1.20053.001",
    site=neon_sites,
    startdate = "2018-01",
    enddate = "2019-12",
    package="expanded",
    token = Sys.getenv("NEON_TOKEN"),
    check.size = FALSE
  )
  saveRDS(temp_all, file.path(save_dir,"temp_all.rda"))
}


# 4 - Elevation of surface water
if(file.exists(file.path(save_dir,"swe_all.rda"))) {
  swe_all <- readRDS(file.path(save_dir,"swe_all.rda"))
} else {
  print("Downloading elevation of surface water data")
  swe_all <-  neonUtilities::loadByProduct(
    dpID="DP1.20016.001",
    site=neon_sites,
    startdate = "2018-01",
    enddate = "2019-12",
    package="expanded",
    token = Sys.getenv("NEON_TOKEN"),
    check.size = FALSE
  )
  saveRDS(swe_all, file.path(save_dir,"swe_all.rda"))
}

print("Finished data download")
#  End data download
# ---------------------------------------------------------------------------- #
#     Conduct analysis by site
# ---------------------------------------------------------------------------- #

# Extract water quality units from variables file:
origUnits <- unlist(lapply(waq_cols_sel, function(x) waq$variables_20288$units[which(waq$variables_20288$fieldName==x)] ))
dfUnits <- data.frame(orig_cols = waq_cols_sel, units = origUnits)
print(dfUnits)



for(site in neon_sites){
  print(paste0("Conducting analysis for NEON siteID = ",site))
  
  # Filter for water quality data at NEON's downstream stations
  water_quality_site <- waq$waq_instantaneous %>% base::subset(siteID == site) %>%
    dplyr::as_tibble() %>%
    dplyr::filter(horizontalPosition == "102") %>% # select the downstream station
    dplyr::select(dplyr::all_of(waq_cols_sel)) %>%
    janitor::clean_names() %>%
    dplyr::mutate(
      start_date_time = lubridate::ymd_hms(start_date_time) - 
        lubridate::second(start_date_time)
          ) %>%
    dplyr::rename(
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
      label_fdom = f_dom_final_qf
    )
  # TODO add 15-min aggregation interval here??

  # Extract nitrate at 15 minute intervals
  nitrate_site <- nsw$NSW_15_minute %>% subset(siteID==site) %>%
    as_tibble() %>%
    janitor::clean_names() %>%
    # mutate( start_date_time = lubridate::ymd_hms(start_date_time) - second(start_date_time)) %>%
    select(c(start_date_time, surf_water_nitrate_mean, final_qf)) %>%
    rename(
      nitrate_mean = surf_water_nitrate_mean,
      label_nitrate_mean = final_qf
    ) %>% padr::thicken(interval = '15 min') %>% # rounding ensures 15 minute intervals
    group_by(start_date_time_15_min) %>% 
    summarise(nitrate_mean = mean(nitrate_mean,na.rm=TRUE), label_nitrate_mean = max(label_nitrate_mean,na.rm = TRUE)) %>% 
    rename(start_date_time = start_date_time_15_min)
  
  
  # Extract temperature at 15 minute intervals for site of interest
  temp_site <- temp_all$TSW_5min %>% subset(siteID==site) %>%
    as_tibble() %>%
    janitor::clean_names()  %>%
    # filter(horizontal_position == "101") %>%
    filter(horizontal_position == "102") %>%
    # mutate(start_date_time = ymd_hms(start_date_time) - second(start_date_time)) %>%
    select(c(start_date_time, surf_water_temp_mean, final_qf)) %>%
    rename(
      temp_mean = surf_water_temp_mean,
      label_temp_mean = final_qf
    ) %>% padr::thicken(interval = '15 min') %>% 
    group_by(start_date_time_15_min) %>% 
    summarise(temp_mean = mean(temp_mean,na.rm=TRUE), label_temp_mean = max(label_temp_mean,na.rm=TRUE)) %>% 
    rename(start_date_time = start_date_time_15_min)
  
  # Extract surface water elevation at 15 minute intervals for site of interest
  elev_site <- swe_all$EOS_5_min %>% subset(siteID == site) %>%
    as_tibble() %>%
    janitor::clean_names() %>%
    filter(horizontal_position=="102") %>%
    # mutate(start_date_time= ymd_hms(start_date_time) - second(start_date_time)) %>%
    select(c(start_date_time, surfacewater_elev_mean, s_wat_elev_final_qf)) %>%
    rename(
      elev = surfacewater_elev_mean,
      label_elev_mean = s_wat_elev_final_qf
    ) %>% padr::thicken(interval = '15 min') %>% 
    group_by(start_date_time_15_min) %>% 
    summarise(elev = mean(elev,na.rm=TRUE), label_elev_mean = max(label_elev_mean,na.rm = TRUE)) %>% 
    rename(start_date_time = start_date_time_15_min)
  
  # Join water_quality and nitrate, and clean up
  data_site <- dplyr::left_join(
    nitrate_site,
    temp_site,
    by = "start_date_time")
  
  data_site <- dplyr::left_join(
    data_site,
    water_quality_site,
    by = "start_date_time") 
  
  data_site <- dplyr::left_join(
    data_site,
    elev_site,
    by = "start_date_time"
  )
  
  data_site <- data_site %>%
    dplyr::mutate(
      spec_cond = if_else(label_spec_cond == 1, NA_real_, spec_cond),
      oxygen = if_else(label_oxygen == 1, NA_real_, oxygen),
      oxygen_sat = if_else(label_oxygen_sat == 1, NA_real_, oxygen_sat),
      turbidity = if_else(label_turbidity == 1, NA_real_, turbidity),
      ph = if_else(label_ph == 1, NA_real_, ph),
      chloro = if_else(label_chloro == 1, NA_real_, chloro),
      fdom = if_else(label_fdom == 1, NA_real_, fdom),
      temp_mean = if_else(label_temp_mean == 1, NA_real_, temp_mean),
      nitrate_mean = if_else(label_nitrate_mean == 1, NA_real_ , nitrate_mean),
      elev = if_else(label_elev_mean == 1, NA_real_, elev)
    )  %>%
    dplyr::distinct(start_date_time, .keep_all=TRUE)
  
  #clear impossible data points
  data_site$turbidity = if_else(data_site$turbidity < 0, NA_real_ , data_site$turbidity)
  # data_site$spec_cond = if_else(data_site$spec_cond < 100, NA_real_ , data_site$spec_cond) #TODO Edit GL - low SpC can be valid
  # data_site$nitrate_mean = if_else(data_site$nitrate_mean  > 25, NA_real_ , data_site$nitrate_mean) #TODO Edit GL - high NO3 can be valid
  data_site$nitrate_mean = if_else(data_site$nitrate_mean <= 0, NA_real_ , data_site$nitrate_mean)
  
  
  # Create a time_recod column
  #TODO Problem here, failed to gap-fill first!
  # data_site$time_recod <- seq(from = 1, to = nrow(data_site))
  
  # Edit GL:
  # First regularize timeseries to 15-min intervals
  
  data_site <- data_site %>%
        padr::pad(interval = "15 mins")
  data_site$time_recod <- seq(from = 1, to = nrow(data_site))
  

  ###############################
  ### SI 1 - Original timeseries ###
  ###############################
  
  data<-data_site
  dataPlot <- data %>% dplyr::select(start_date_time, nitrate_mean, turbidity, 
                                     temp_mean, spec_cond, oxygen ) %>% 
    dplyr::rename(time= start_date_time, Nitrate=nitrate_mean,
           Temp=temp_mean, SpC=spec_cond, DO= oxygen )

  
  # Assign units: nitrate : micromoles/L; fDOM : QSU; turbidity : FNU, DO : mg/L; SpC microsiemens per Centimeter; 
  dfUnits <-  data.frame(var = names(dataPlot)[-1],
             units = c("\u03BCMol/L","FNU","°C", "\u03BCS/cm", "mg/L"))
  dfUnits <- rbind(dfUnits,data.frame(var="SWE",units="m")) # add in elevation

  longDf <- dataPlot %>% pivot_longer(-time, names_to = "variable") 
  longDf <- merge(longDf,dfUnits, by.x = "variable", by.y = "var")
  longDf$lablUnit <- paste0(longDf$variable, " [",longDf$units,"]")
  
  ## Plot data
  plotTs <- longDf %>% #dataPlot %>%
    # dplyr::select(start_date_time, nitrate_mean, turbidity, temp_mean, spec_cond, oxygen ) %>%
    # pivot_longer(-time, names_to = "variable") %>%
    ggplot(aes(x = time , y = value, col = variable)) +
    geom_line(alpha=0.6) +
    facet_grid(rows = vars(lablUnit),  scales = "free", space="free_x") +
    ggplot2::xlab("Time") +
    ggplot2::ylab("") +
    scale_x_datetime(date_labels ="%Y-%m") +
    theme_minimal() +
    theme(legend.position = "none") + 
    ggtitle(paste0(site, " timeseries data"))
    # theme_light()
  ggsave(plot = plotTs,
         filename = file.path(plot_dir,paste0("timeseries_",site,".png")),
         height = 7,width = 5, units = "in")
  ggsave(plot = plotTs,
         filename = file.path(plot_dir,paste0("timeseries_",site,".svg")),
         height = 7,width = 5, units = "in")
  
  #################
  ### Figure X  ###
  #################
  # CustomFacetWrap <- ggproto(
  #   "CustomFacetWrap", FacetWrap,
  #   init_scales = function(self, layout, x_scale = NULL, y_scale = NULL, params) {
  #     # make the initial x, y scales list
  #     scales <- ggproto_parent(FacetWrap, self)$init_scales(layout, x_scale, y_scale, params)
  #     
  #     if(is.null(params$scale_overrides)) return(scales)
  #     
  #     max_scale_x <- length(scales$x)
  #     max_scale_y <- length(scales$y)
  #     
  #     # ... do some modification of the scales$x and scales$y here based on params$scale_overrides
  #     for(scale_override in params$scale_overrides) {
  #       which <- scale_override$which
  #       scale <- scale_override$scale
  #       
  #       if("x" %in% scale$aesthetics) {
  #         if(!is.null(scales$x)) {
  #           if(which < 0 || which > max_scale_x) stop("Invalid index of x scale: ", which)
  #           scales$x[[which]] <- scale$clone()
  #         }
  #       } else if("y" %in% scale$aesthetics) {
  #         if(!is.null(scales$y)) {
  #           if(which < 0 || which > max_scale_y) stop("Invalid index of y scale: ", which)
  #           scales$y[[which]] <- scale$clone()
  #         }
  #       } else {
  #         stop("Invalid scale")
  #       }
  #     }
  #     
  #     # return scales
  #     scales
  #   }
  # )
  # 
  # # Custom function from https://github.com/remichel/rmTools/blob/master/R/facet_wrap_costum.R
  # facet_wrap_custom <- function(..., scale_overrides = NULL) {
  #   # take advantage of the sanitizing that happens in facet_wrap
  #   facet_super <- facet_wrap(...)
  #   
  #   # sanitize scale overrides
  #   if(inherits(scale_overrides, "scale_override")) {
  #     scale_overrides <- list(scale_overrides)
  #   } else if(!is.list(scale_overrides) ||
  #             !all(vapply(scale_overrides, inherits, "scale_override", FUN.VALUE = logical(1)))) {
  #     stop("scale_overrides must be a scale_override object or a list of scale_override objects")
  #   }
  #   
  #   facet_super$params$scale_overrides <- scale_overrides
  #   
  #   ggproto(NULL, CustomFacetWrap,
  #           shrink = facet_super$shrink,
  #           params = facet_super$params
  #   )
  # }
  # 
  # scale_override <- function(which, scale) {
  #   # Custom function from https://github.com/remichel/rmTools/blob/master/R/scale_override.R
  #   if(!is.numeric(which) || (length(which) != 1) || (which %% 1 != 0)) {
  #     stop("which must be an integer of length 1")
  #   }
  #   
  #   if(is.null(scale$aesthetics) || !any(c("x", "y") %in% scale$aesthetics)) {
  #     stop("scale must be an x or y position scale")
  #   }
  #   
  #   structure(list(which = which, scale = scale), class = "scale_override")
  # }
  
  # Add surface water elevation to units

  
  if(site == "ARIK" || site == "CARI"){
    # arik 5 days
    start <- which(data_site$start_date_time == "2019-07-01 10:00:00 UTC")
    end <- which(data_site$start_date_time == "2019-07-06 10:00:00 UTC")
  } else if (site == "LEWI"){
    # TODO determine if this is the appropriate date range/year
    start <- which(data_site$start_date_time == "2019-03-22 10:00:00 UTC")
    end <- which(data_site$start_date_time == "2019-03-27 10:00:00 UTC")
  }
  
  # TODO add water level
  
  data_piv <-data_site[c(start:end),] %>%  dplyr::select(start_date_time, nitrate_mean, turbidity, temp_mean, spec_cond, oxygen, elev) %>%
    dplyr::rename(start_date_time = start_date_time, "Nitrate" = nitrate_mean, DO = oxygen, SpC =  spec_cond, SWE = elev, Temp = temp_mean, turbidity = turbidity) %>%
    pivot_longer(-start_date_time, names_to = "variable")
  data_piv <- merge(data_piv,dfUnits, by.x = "variable", by.y = "var")
  data_piv$lablUnit <- paste0(data_piv$variable, " [",data_piv$units,"]")
  
  
  plotFig2_top <- data_piv %>%
    ggplot(aes(x = start_date_time , y = value, col = variable)) + ggtitle("Arikaree") +
    geom_line() +
    facet_grid(rows = vars(lablUnit),  scales = "free", space="free_x") +
  
    # facet_wrap_custom(~ variable,  ncol = 1, scales = "free_y", strip.position="right",  scale_overrides = list(
    #   scale_override(1, scale_y_continuous(breaks = c(535,540))),
    #   scale_override(2, scale_y_continuous(breaks = c(1179.75,1179.79))),
    #   scale_override(3, scale_y_continuous(breaks = c(4,6))),
    #   scale_override(4, scale_y_continuous(breaks = c(2.5,7.5))),
    #   scale_override(5, scale_y_continuous(breaks = c(21,25))),
    #   scale_override(6, scale_y_continuous(breaks = c(3,7)))
    # )) +
    labs(y = "", x = "time") +
    scale_x_datetime(date_labels ="%b:%d", breaks=scales::date_breaks("1 day")) +
    theme_minimal() +
    theme(legend.position = "none",
          strip.text.y = element_text(size = 6),
          axis.text.y= element_text(size = 8),
          axis.text.x= element_text(size = 8))
  
  if(site == "ARIK" || site == "LEWI"){
    # arik with peack
    start2 <- which(data_site$start_date_time == "2019-07-20 19:15:00 UTC")
    end2 <- which(data_site$start_date_time == "2019-07-23 19:15:00 UTC")
  } else if (site == "CARI"){
    # TODO determine if this is the appropriate date range/ year
    start2 <- which(data_site$start_date_time == "2019-08-01 22:00:00 UTC")
    end2 <- which(data_site$start_date_time == "2019-08-06 10:00:00 UTC")
  }
  
  data_piv2<-data_site[c(start:end),] 
  data %>%
    dplyr::select(start_date_time, nitrate_mean, turbidity, temp_mean, spec_cond, oxygen, elev ) %>%
    dplyr::rename(start_date_time = start_date_time, Elev. = elev, "N03-" = nitrate_mean, O2 = oxygen, Cond =  spec_cond, Temp. = temp_mean, Turb. = turbidity) %>%
    pivot_longer(-start_date_time, names_to = "variable") %>%
    ggplot(aes(x = start_date_time , y = value, col = variable)) +
    geom_point() +
    facet_wrap_custom(~ variable,  ncol = 1, scales = "free_y", strip.position="right",  scale_overrides = list(
      scale_override(1, scale_y_continuous(breaks = c(400,550))),
      scale_override(2, scale_y_continuous(breaks = c(1179.52,1179.60))),
      scale_override(3, scale_y_continuous(breaks = c(0,15))),
      scale_override(4, scale_y_continuous(breaks = c(2,6))),
      scale_override(5, scale_y_continuous(breaks = c(20,27.5))),
      scale_override(6, scale_y_continuous(breaks = c(10,40)))
    )) +
    labs(y = "", x = "") +
    scale_x_datetime(date_labels ="%b %d", breaks=date_breaks("1 day")) +
    
    theme(legend.position = "none",
          strip.text.y = element_text(size = 6),
          axis.text.y= element_text(size = 8),
          axis.text.x= element_text(size = 8))
  
  
  ################
  ### Figure 1 (boxplots) ### 
  ################
  
  
  data_plot_site<-cbind(data_site[,c(2,4,6,8)], log(data_site$turbidity +1))
  colnames(data_plot_site) <- c("Nitrate","Temperature","Spec. Cond.", "Oxygen C.", "log_Turbidity")
  data_plot_site<-melt(data_plot_site)
  data_plot_site$site<-rep("Arikaree", length = length(data_plot_site[,1]))
  
  gg<-ggplot(data = data_plot_site[which(data_plot_site$variable == "Nitrate"),], aes(factor(1), value, color = site)) +
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
  
  gg1<-ggplot(data = data_plot_site[which(data_plot_site$variable == "Temperature"),], aes(factor(1), value, color = site)) +
    geom_boxplot() +  facet_wrap(~variable,scales="free",ncol=3)+ 
    scale_color_manual(values=c("#440154FF", "#277F8EFF", "#9FDA3AFF"))+
    theme(axis.text.x=element_blank(),
          axis.text.y= element_text(size = 11),
          axis.title.x=element_blank(), 
          axis.title.y = element_blank(),
          legend.position= "none",
          legend.title = element_blank(),
          strip.text.x = element_text(size = 15))
  
  gg2<-ggplot(data = data_plot_site[which(data_plot_site$variable == "Spec. Cond."),], aes(factor(1), value, color = site)) +
    geom_boxplot() +  facet_wrap(~variable,scales="free",ncol=3)+ 
    scale_color_manual(values=c("#440154FF", "#277F8EFF", "#9FDA3AFF"))+
    theme(axis.text.x=element_blank(),
          axis.text.y= element_text(size = 11),
          axis.title.x=element_blank(), 
          axis.title.y = element_blank(),
          legend.position= "none",
          legend.title = element_blank(),
          strip.text.x = element_text(size = 15))
  gg3<-ggplot(data = data_plot_site[which(data_plot_site$variable == "Oxygen C."),], aes(factor(1), value, color = site)) +
    geom_boxplot() +  facet_wrap(~variable,scales="free",ncol=3)+ 
    scale_color_manual(values=c("#440154FF", "#277F8EFF", "#9FDA3AFF"))+
    theme(axis.text.x=element_blank(),
          axis.text.y= element_text(size = 11),
          axis.title.x=element_blank(), 
          axis.title.y = element_blank(),
          legend.position= "none",
          legend.title = element_blank(),
          strip.text.x = element_text(size = 15))
  gg4<-ggplot(data = data_plot_site[which(data_plot_site$variable == "log_Turbidity"),], aes(factor(1), value, color = site)) +
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
  model_gam_site<-gam::gam(nitrate_mean ~ spec_cond + oxygen + turbidity + temp_mean, data = data_site)
  
  # Stepwise selection of variables
  step_model_gam<-gam::step.Gam(model_gam_site, scope=list("spec_cond"=~1+spec_cond+s(spec_cond,4)+s(spec_cond,6)+s(spec_cond,12),
                                                           "oxygen"=~1+oxygen+s(oxygen,4)+s(oxygen,5)+s(oxygen,6)+s(oxygen,12),
                                                           "turbidity"=~1+turbidity+s(turbidity,4)+s(turbidity,5)+s(turbidity,6)+s(turbidity,12),
                                                           "temp_mean"=~1+temp_mean+s(temp_mean,4)+s(temp_mean,6)+s(temp_mean,12)))
  # best model selected 
  best_model_gam_site<-mgcv::gam(nitrate_mean ~  s(spec_cond, k = 12) + s(oxygen, k = 12) +  s(temp_mean, k = 12) + s(log(turbidity + 1 ), k = 12)+ s(time_recod, k = 12),  data = data_site, na.action=na.exclude)
  
  #fill where we deleted na
  data_site <- data_site %>%
    mutate(gam_res = residuals(best_model_gam_site)) %>%
    as_tsibble(index=start_date_time) %>%
    fill_gaps()
  
  
 
  
  # Same model for visualization
  best_model_gam_site<-mgcv::gam(nitrate_mean ~  s(spec_cond, k = 12) + s(oxygen, k = 12) + 
                                   +  s(temp_mean, k = 12) +
                                   s(log(turbidity+1) ,  k = 12)+ s(time_recod, k = 12), data = data_site)
  
  
  ## RMSE and GCV
  site_res <- c(site, round(best_model_gam_site$gcv.ubre,2))
  site_res
  
  ################
  ### Figure 3 ###
  ################
  
  a <- getViz(best_model_gam_site)
  
  site1 <- plot(sm(a, 1)) + l_fitLine(colour = "red") + l_ciLine(mul = 5, colour = "blue", linetype = 2) + xlab("Spec. Cond.") + ylab("")  +
    theme_classic()
  site2 <- plot(sm(a, 2)) + l_fitLine(colour = "red") + l_ciLine(mul = 5, colour = "blue", linetype = 2) + xlab("Oxygen") +ylab("") + ylim(-5,5) + theme_classic()
  site3 <- plot(sm(a, 3)) + l_fitLine(colour = "red") +  l_ciLine(mul = 5, colour = "blue", linetype = 2) +  xlab("Temperature") +ylab("")+ ylim(-5,5) + 
    theme_classic()
  site4 <- plot(sm(a, 4)) + l_fitLine(colour = "red") + l_ciLine(mul = 5, colour = "blue", linetype = 2) + xlab("log(Turbidity)") + ylab("") +
    theme_classic()
  site5 <- plot(sm(a, 5)) + l_fitLine(colour = "red") +  l_ciLine(mul = 5, colour = "blue", linetype = 2) + xlab("Time") + ylab("") +
    theme_classic()
  
  
  gridPrint(site1, site2, site3, site4,  site5, ncol = 1)
  
  
  ############################################
  ### Calculating importance of variables ###
  ###########################################
  
  #Deviance of the GAMM 
  d_null<-deviance(gam(nitrate_mean ~  1 , data = data_site))
  ar_model<-auto.arima(residuals(best_model_gam_site))  # selecting best ARIMA model
  D_GAMM_total <- deviance(best_model_gam_site)
  D_RES <- sum(residuals(ar_model)^2, na.rm=T)
  Deviance_GAMM_total <- (d_null - D_RES) / d_null
  
  
  # Deviances of GAMMs minus one by one covariates
  GAM_elev<-mgcv::gam(nitrate_mean ~  s(spec_cond, k = 12) + s(oxygen, k = 12) +  s(temp_mean, k = 12) +
                        s(log(turbidity+1),  k = 12)+ s(time_recod, k = 12), data = data_site, sp=best_model_gam_site$sp[c(1:5)])
  D_GAM <- deviance(GAM_elev)
  D_RES <- sum(residuals(arima(residuals(GAM_elev), order = c(2,0,5)))^2, na.rm=T)
  Deviance_elev <- (D_GAMM_total - D_RES) / D_GAMM_total
  Deviance_elev
  
  GAM_time<-mgcv::gam(nitrate_mean ~  s(spec_cond, k = 12) + s(oxygen, k = 12) +  s(temp_mean, k = 12) +
                        s(log(turbidity+1),  k = 12)+ s(elev, k = 12), data = data_site, sp=best_model_gam_site$sp[c(1:4,6)])
  D_GAM <- deviance(GAM_time)
  D_RES <- sum(residuals(arima(residuals(GAM_time), order = c(2,0,5)))^2, na.rm=T)
  Deviance_time <- (D_GAMM_total - D_RES) / D_GAMM_total
  Deviance_time
  
  
  GAM_turbi<-mgcv::gam(nitrate_mean ~  s(spec_cond, k = 12) + s(oxygen, k = 12) +  s(temp_mean, k = 12) + s(time_recod, k = 12)+
                         data = data_site, sp=best_model_gam_site$sp[c(1:3,5,6)])
  D_GAM <- deviance(GAM_turbi)
  D_RES <- sum(residuals(arima(residuals(GAM_turbi), order = c(2,0,5)))^2, na.rm=T)
  Deviance_turbi <- (D_GAMM_total - D_RES) / D_GAMM_total
  Deviance_turbi
  
  
  GAM_temp<-mgcv::gam(nitrate_mean ~   s(spec_cond, k = 12) +  s(oxygen, k = 12) + s(log(turbidity+1),  k = 12)+s(time_recod, k = 12)+
                        s(elev, k = 12), data = data_site, sp=best_model_gam_site$sp[c(1,3:6)])
  D_GAM <- deviance(GAM_temp)
  D_RES <- sum(residuals(arima(residuals(GAM_temp), order = c(2,0,5)))^2, na.rm=T)
  Deviance_temp <- (D_GAMM_total - D_RES) / D_GAMM_total
  Deviance_temp
  
  GAM_oxygen<-mgcv::gam(nitrate_mean ~   s(spec_cond, k = 12) +  s(temp_mean, k = 12) + s(log(turbidity+1),  k = 12)+s(time_recod, k = 12)+
                          s(elev, k = 12), data = data_site, sp=best_model_gam_site$sp[c(1:2,4:6)])
  D_GAM <- deviance(GAM_oxygen)
  D_RES <- sum(residuals(arima(residuals(GAM_oxygen), order = c(2,0,5)))^2, na.rm=T)
  Deviance_oxygen <- (D_GAMM_total - D_RES) / D_GAMM_total
  Deviance_oxygen
  
  
  GAM_cond<-mgcv::gam(nitrate_mean ~   s(oxygen, k = 12) +  s(temp_mean, k = 12) + s(log(turbidity+1),  k = 12)+s(time_recod, k = 12)+
                        s(elev, k = 12), data = data_site, sp=best_model_gam_site$sp[c(2:6)])
  D_GAM <- deviance(GAM_cond)
  D_RES <- sum(residuals(arima(residuals(GAM_cond), order = c(2,0,5)))^2, na.rm=T)
  Deviance_cond <- (D_GAMM_total - D_RES) / D_GAMM_total
  Deviance_cond
  
  # deviances at arikaree
  site_dev<-c(Deviance_elev, Deviance_time,Deviance_turbi, Deviance_temp, Deviance_oxygen, Deviance_cond )
  site_imp <- data.frame("Importance" = (1-site_dev)*100, "Var" = c("Elevation", "Time", "Log(turbidity)", "Temperature", "Dissolved Oxygen", "Special Conductance"))
  site_imp
  
  ### aAIC calculation ###
  
  aAIC_site <- 14973*log(best_model_gam_site$sig2) + (2*sum(best_model_gam_site$edf)) # for the GAM
  GAMM_site<-forecast::auto.arima(best_model_gam_site$residuals) # 5 df
  aAIC_GAMM_site <- (14973*log(GAMM_site$sigma2)) + (sum(best_model_gam_lewi$edf) + 5) # for the GAMM
  
  
  ################
  ### FIGURE 4 ###
  ################
  
  total_imp<- site_imp
  total_imp$site<-c( rep("Arikaree", 6))
  
  total_imp %>%
    ggplot(aes(x = Importance, y = Var, color = site, shape = site)) +
    geom_point(size = 4) +
    xlab(expression(paste("Variable importance (% of total deviance)")))+
    ylab("") +
    ggtitle("")+ 
    scale_color_manual(values=c("#440154FF", "#277F8EFF", "#9FDA3AFF"))+
    scale_shape_manual(values=c(15,16,17 ))+
    theme(legend.text = element_text(size = 13 ),
          legend.title = element_blank(),
          legend.position = c(0.85, 0.12),
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),
          axis.title.x = element_text(size=14),
          panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          axis.line = element_line(color = "grey") )
  
}

