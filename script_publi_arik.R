#' @title Data acquisition, processing, and modeling to understand links between
#'  water quality and nitrate 
#' @author Claire Kermorvant, Guy Litt
#' @description Companion code to "Understanding links between water-quality 
#' variables and nitrate concentration in freshwater streams using 
#' high-frequency sensor data" https://arxiv.org/abs/2106.01719
#' @details Requires using neonUtilities >= version 2.2.1
# Changelog/contributions
# 2021-05-21 Originally created, Claire Kermorvant
# 2023-04-05 Update data acquisition & results display, generalized to multiple sites,
#    Added temporal aggregation method (mean/max)
#    Manually remove anomalous CARI turbidity data, Guy Litt
# 2023-04-09 Generate waq characteristics .csv ,
#            Condensed Fig2 panel labels & reset 10pm GMT bounds, 
#            Update GAM with surface water elevation 
#            Subset model dataset to user-defined time ranges, Guy Litt
# 2023-04-13 Add option to regularize timeseries
# 
# TODO change n to be length of timeseries rather than hard-coded
# TODO address hard-coding of x-axis date in Fig 3
#TODO substitute end_date_time for start_date_time (low priority)
# TODO re-upload all figures to overleaf
# TODO ARIK period of data modeled in Table 2 mis-represents model range??

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
# install package in github used for custom plotting
# devtools::install_github(repo = "https://github.com/remichel/rmTools.git")
###########################################
### download data of interest from NEON ###
###########################################
neon_sites <-  c("ARIK","LEWI","CARI")

# Define the file saving directory
user_base_dir <- Sys.getenv("USERPROFILE")
save_dir <- file.path(user_base_dir,"Documents","waq_nitrate_analysis")
plot_dir_base <- file.path(save_dir, "plots")


neonPkgVer <- utils::packageVersion("neonUtilities")
if(neonPkgVer < base::package_version('2.2.1')){
  warning("The neonUtilities package version should at least be >= 2.2.1")
}

# Define the instantaneous water quality data columns of interest
waq_cols_sel <- c("startDateTime","horizontalPosition",  "specificConductance","specificCondFinalQF",
                  "dissolvedOxygen","dissolvedOxygenFinalQF","dissolvedOxygenSaturation",
                  "dissolvedOxygenSatFinalQF","pH","pHFinalQF",               
                  "chlorophyll","chlorophyllFinalQF","turbidity",                
                  "turbidityFinalQF","fDOM","fDOMFinalQFSciRvw")
# The dimension of the bases used to represent the smooth term. See ?mgcv::s()
smoothDimK <- 6
stepGam <- FALSE # Boolean, perform step-wise check on best GAM? Helps decide smoothDimK
regularizeTs <- c(TRUE,FALSE)[1] # Should timeseries be regularized to 15-min intervals?
# Define the time range to apply the model
dataSetup <- base::data.frame(site = c("ARIK","CARI","LEWI"),
                                   startDate = as.POSIXct(c("2018-09-01","2018-06-01","2018-01-01"),tz="GMT"),
                                   endDate = as.POSIXct(c("2020-01-01","2019-12-31","2019-12-31"),tz="GMT"),
                              modlStart = as.POSIXct(c("2018-09-01","2018-06-01","2018-01-01"), tz="GMT"),
                              modlEnd = as.POSIXct(c("2019-12-31","2019-10-31","2019-12-31"),tz="GMT"),
                              elevHor = c("101","102","102"),
                              tempHor = c("101","102","102"),
                              defaultHor = "102")



if(regularizeTs){regStr <- "Reg_"}else{regStr <- "NonReg_"}     
sitesElev102 <- dataSetup$site[which(dataSetup$elevHor=="102")]
sitesTemp102 <- dataSetup$site[which(dataSetup$tempHor=="102")]
uniqId <- paste0("smoothDimK",smoothDimK, "_",regStr,"TempDownStrm_",paste0(sitesTemp102,collapse=""),
                  "_ElevDownStrm_",paste0(sitesElev102,collapse=""))

plot_dir <- file.path(plot_dir_base,uniqId)
# Create the save directory
if(!dir.exists(plot_dir)){
  print(paste0("Creating the following directory path: ", plot_dir))
  dir.create(plot_dir,recursive=TRUE)
}
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

siteNameDf <- base::list(siteID = c("ARIK","CARI","LEWI"),
                         siteName = c("Arikaree River","Caribou Creek","Lewis Run"))

# Extract water quality units from variables file:
origUnits <- unlist(lapply(waq_cols_sel, function(x) waq$variables_20288$units[which(waq$variables_20288$fieldName==x)] ))
dfUnits <- data.frame(orig_cols = waq_cols_sel, units = origUnits)
print(dfUnits)

# Create lists that will be populated with various objects for compilation across sites
lsRsltsWaqChars <- base::list()
lsPlotTs <- base::list()
lsFig2 <- base::list()
lsBoxFig1 <- base::list()
lsSmoothFig_3 <- base::list()
lsImportFig4 <- base::list()
lsSiteImp <- list()
lsAic <- list()
lsARIMAmodls <- list()
lsGamGammARMASumm <- list()
lsDevCalc <- list()
lsImport <- list()
for(site in neon_sites){
  siteName <- siteNameDf$siteName[siteNameDf$siteID == site]
  print(paste0("Conducting analysis for NEON siteID = ",site))
  
  # Filter for water quality data at NEON's downstream stations
  water_quality_site <- waq$waq_instantaneous %>% base::subset(siteID == site) %>%
    dplyr::as_tibble() %>%
    dplyr::filter(horizontalPosition == dataSetup$defaultHor) %>% # select the downstream station
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
      fdom = f_dom,
      label_fdom = f_dom_final_qf_sci_rvw
    )

  # NEON Issue 36904 water quality algorithm error caused additional timestamps
  # infrequently. This was fixed in the 2023 data release.
  idxsDupTime <- which(duplicated(water_quality_site$start_date_time))
  if(length(idxsDupTime)>0){
    water_quality_site <- water_quality_site[-idxsDupTime,]
  }
  
  # Extract nitrate at 15 minute intervals
  nitrate_site <- nsw$NSW_15_minute %>% subset(siteID==site) %>%
    as_tibble() %>%
    janitor::clean_names() %>%
    mutate( start_date_time = lubridate::ymd_hms(start_date_time) - second(start_date_time)) %>%
    select(c(start_date_time, surf_water_nitrate_mean, final_qf)) %>%
    rename(
      nitrate_mean = surf_water_nitrate_mean,
      label_nitrate_mean = final_qf)  
  
  # Extract temperature at 15 minute intervals for site of interest
  temp_site <- temp_all$TSW_5min %>% 
    subset(siteID==site) %>%
    as_tibble() %>%
    janitor::clean_names()  %>%
    filter(horizontal_position == dataSetup$tempHor[dataSetup$site==site]) %>%
    mutate(start_date_time = ymd_hms(start_date_time) - second(start_date_time)) %>%
    select(c(start_date_time, surf_water_temp_mean, final_qf)) %>%
    rename(
      temp_mean = surf_water_temp_mean,
      label_temp_mean = final_qf)
  
  # Extract surface water elevation at 15 minute intervals for site of interest
  elev_site <- swe_all$EOS_5_min %>% subset(siteID == site) %>%
    as_tibble() %>%
    janitor::clean_names() %>%
    filter(horizontal_position==dataSetup$elevHor[dataSetup$site==site]) %>%
    mutate(start_date_time= ymd_hms(start_date_time) - second(start_date_time)) %>%
    select(c(start_date_time, surfacewater_elev_mean, s_wat_elev_final_qf)) %>%
    rename(
      elev = surfacewater_elev_mean,
      label_elev_mean = s_wat_elev_final_qf) #%>%

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
  # Remove poor quality CARI turbidity data:
  if(length(unique(data_site$start_date_time)) != nrow(data_site)){
    stop("Problem with duplicate timestamps")
  }
  #clear impossible data points
  if(site == "CARI"){
    data_site$turbidity[data_site$turbidity > 500] <- NA_real_
    data_site$nitrate_mean <- if_else(data_site$nitrate_mean > 50, NA_real_ , data_site$nitrate_mean)
    data_site$chloro = if_else(data_site$chloro > 200, NA_real_ , data_site$chloro)
    data_site$spec_cond = if_else(data_site$spec_cond > 300, NA_real_ , data_site$spec_cond)
  } else if (site == "ARIK"){
    data_site$turbidity = if_else(data_site$turbidity < 0, NA_real_ , data_site$turbidity)
    data_site$spec_cond = if_else(data_site$spec_cond < 100, NA_real_ , data_site$spec_cond)
    data_site$nitrate_mean = if_else(data_site$nitrate_mean  > 25, NA_real_ , data_site$nitrate_mean)
    data_site$nitrate_mean = if_else(data_site$nitrate_mean <= 0, NA_real_ , data_site$nitrate_mean)
    data_site$fdom = if_else(data_site$fdom  > 140, NA_real_ , data_site$fdom)
  } else if (site == "LEWI"){

    data_site$nitrate_mean = if_else(data_site$nitrate_mean <= 5, NA_real_ , data_site$nitrate_mean)
    data_site$turbidity = if_else(data_site$turbidity < 0, NA_real_ , data_site$turbidity)
    data_site$temp_mean = if_else(data_site$temp_mean < 0, NA_real_ , data_site$temp_mean)
  }
  
  # First regularize timeseries to 15-min intervals
  if(regularizeTs){
    data_site <- data_site %>%
      padr::pad(interval = "15 mins")
  }
  
  # Subset data to modlStart and modlEnd timeframes:
  idxsModl <- intersect(which(data_site$start_date_time >= dataSetup$modlStart[dataSetup$site==site]),
                        which(data_site$start_date_time <= dataSetup$modlEnd[dataSetup$site==site]))
  data_site <- data_site[idxsModl,]
  
  data_site$time_recod <- seq(from = 1, to = nrow(data_site))
  
  # Check the total number of NAs by column
  naCountCols <- lapply(ncol(data_site), function(j) length(which(is.na(data_site[,j]))))
  dfNaCount <- data.frame(colname = names(data_site), naCount = unlist(naCountCols))
  ###############################
  ### SI 1 - Original timeseries ###
  ###############################
  
  data<-data_site
  dataPlot <- data %>% dplyr::select(start_date_time, nitrate_mean, turbidity, 
                                     temp_mean, elev, spec_cond, oxygen ) %>% 
    dplyr::rename(time= start_date_time, Nitrate=nitrate_mean,
           Temp=temp_mean,SWE=elev, SpC=spec_cond, DO= oxygen )

  
  # Assign units: nitrate : micromoles/L; fDOM : QSU; turbidity : FNU, DO : mg/L; SpC microsiemens per Centimeter; 
  dfUnits <-  data.frame(var = names(dataPlot)[-1],
             units = c("\u03BCMol/L","FNU","°C", "m", "\u03BCS/cm", "mg/L"))
  # dfUnits <- rbind(dfUnits,data.frame(var="SWE",units="m")) # add in elevation

  longDf <- dataPlot %>% pivot_longer(-time, names_to = "variable") 
  longDf <- merge(longDf,dfUnits, by.x = "variable", by.y = "var")
  longDf$lablUnit <- paste0(longDf$variable, "\n[",longDf$units,"]")
  
  ## Plot data
  plotTs <- longDf %>% #dataPlot %>%
    # dplyr::select(start_date_time, nitrate_mean, turbidity, temp_mean, spec_cond, oxygen ) %>%
    # pivot_longer(-time, names_to = "variable") %>%
    ggplot(aes(x = time , y = value, col = variable)) +
    geom_point(size=0.01,alpha=0.3) +
    facet_grid(rows = vars(lablUnit),  scales = "free", space="free_x") +
    ggplot2::xlab("Time") +
    ggplot2::ylab("") +
    scale_x_datetime(date_labels ="%Y-%m") +
    theme_minimal() +
    theme(legend.position = "none") + 
    ggtitle(paste0(siteName, " timeseries data"))
    # theme_light()
  ggsave(plot = plotTs,
         filename = file.path(plot_dir,paste0("timeseries_",site,".png")),
         height = 7,width = 5, units = "in")
  ggsave(plot = plotTs,
         filename = file.path(plot_dir,paste0("timeseries_",site,".svg")),
         height = 7,width = 5, units = "in")
  lsPlotTs[[site]] <- plotTs
  
  
  # ----------------------
  #   Results description:
  # -----------------------

  lsRsltsWaqChars[[site]] <- longDf %>% dplyr::group_by(variable) %>% summarise(mean = mean(value,na.rm=TRUE),
                                                     median = median(value,na.rm=TRUE),
                                                     min = min(value,na.rm = TRUE),
                                                     max = max(value,na.rm = TRUE),
                                                     site = site)
  
  #################
  ### Figure 2  ###
  #################
  #TODO ensure SWE is described in figure caption
  if(site == "ARIK" || site == "CARI"){
    # arik 5 days
    start <- which(data_site$start_date_time == "2019-07-01 10:00:00 UTC")
    end <- which(data_site$start_date_time == "2019-07-06 10:00:00 UTC")
  } else if (site == "LEWI"){
    # TODO determine if this is the appropriate date range 
    start <- which(data_site$start_date_time == "2018-03-22 10:00:00 UTC")
    end <- which(data_site$start_date_time == "2018-03-27 10:00:00 UTC")
  }

  sub_data <-data_site[c(start:end),] %>%  dplyr::select(start_date_time, nitrate_mean, turbidity, temp_mean, spec_cond, oxygen, elev) %>%
    dplyr::rename(start_date_time = start_date_time, Nitrate = nitrate_mean, DO = oxygen, SpC =  spec_cond, SWE = elev, Temp = temp_mean, turbidity = turbidity)
  data_piv <- sub_data %>% pivot_longer(-start_date_time, names_to = "variable")
  
  data_piv <- merge(data_piv,dfUnits, by.x = "variable", by.y = "var")
  data_piv$lablUnit <- paste0(data_piv$variable, " [",data_piv$units,"]")
  
  
  plotFig2_top <- data_piv %>%
    ggplot(aes(x = start_date_time , y = value, col = variable)) + 
    ggtitle(siteName) +
    geom_point(size=0.01, alpha=0.5) +
    facet_grid(rows = vars(lablUnit),  scales = "free", space="free_x") +
    labs(y = "", x = "") +
    scale_x_datetime(date_labels ="%Y%m%d", breaks=scales::date_breaks("1 day")) +
    theme_minimal() +
    theme(legend.position = "none",
          strip.text.y = element_text(size = 6),
          axis.text.y= element_text(size = 8),
          axis.text.x= element_text(size = 8, angle=30))
  
  
  
  if(site == "ARIK" || site == "LEWI"){
    # arik with peack
    start2 <- which(data_site$start_date_time == "2019-07-20 22:00:00 UTC")
    end2 <- which(data_site$start_date_time == "2019-07-23 22:00:00 UTC")
  } else if (site == "LEWI"){
    start2 <- which(data_site$start_date_time == "2019-07-21 22:00:00 UTC")
    end2 <- which(data_site$start_date_time == "2019-07-24 22:00:00 UTC")
  } else if (site == "CARI"){
    # TODO determine if this is the appropriate date range/ year
    start2 <- which(data_site$start_date_time == "2019-08-01 22:00:00 UTC")
    end2 <- which(data_site$start_date_time == "2019-08-06 22:00:00 UTC")
  }
  
  data_piv2<-data_site[c(start2:end2),] %>%  dplyr::select(start_date_time, nitrate_mean, turbidity, temp_mean, spec_cond, oxygen, elev) %>%
    dplyr::rename(start_date_time = start_date_time, Nitrate = nitrate_mean, DO = oxygen, SpC =  spec_cond, SWE = elev, Temp = temp_mean, turbidity = turbidity) %>% 
    tidyr::pivot_longer(-start_date_time, names_to = "variable")
  data_piv2 <- merge(data_piv2,dfUnits, by.x = "variable", by.y = "var")
  data_piv2$lablUnit <- paste0(data_piv2$variable, "\n [",data_piv2$units,"]")
  
  
  plotFig2_bttm <- data_piv2 %>%
    # dplyr::select(start_date_time, nitrate_mean, turbidity, temp_mean, spec_cond, oxygen, elev ) %>%
    # dplyr::rename(start_date_time = start_date_time, Elev. = elev, "N03-" = nitrate_mean, O2 = oxygen, Cond =  spec_cond, Temp. = temp_mean, Turb. = turbidity) %>%
    # pivot_longer(-start_date_time, names_to = "variable") %>%
    ggplot(aes(x = start_date_time , y = value, col = variable)) +
    geom_point(size=0.02, alpha=0.5) +
    facet_grid(rows = vars(lablUnit),  scales = "free", space="free_x") +
    labs(y = "", x = "") +
    scale_x_datetime(date_labels ="%Y%m%d", breaks=scales::date_breaks("1 day")) +
    theme_minimal() +
    theme(legend.position = "none",
          strip.text.y = element_text(size = 6),
          axis.text.y= element_text(size = 8),
          axis.text.x= element_text(size = 8,angle=30))
  
  
  # Combine top and bottom plots into single panel:
  
  plotZoomTsFig2 <- gridPrint(plotFig2_top,plotFig2_bttm, ncol = 1)
  lsFig2[[site]] <- plotZoomTsFig2
  
  ggsave(plot=plotZoomTsFig2,
         filename=file.path(plot_dir,paste0("Figure_2_ZoomTS_",site,".png")))
  ggsave(plot=plotZoomTsFig2,
         filename=file.path(plot_dir,paste0("Figure_2_ZoomTS_",site,".png")))

  ################
  ### Figure 1 (boxplots) ### 
  ################
  sel_cols <- c("nitrate_mean","temp_mean","spec_cond","oxygen")
  
  data_plot_site<-cbind(data_site[,sel_cols], log(data_site$turbidity +1))
  colnames(data_plot_site) <- c("Nitrate","Temperature","Spec. Cond.", "Oxygen C.", "log_Turbidity")
  data_plot_site<-melt(data_plot_site)
  data_plot_site$site<-rep(siteName, length = length(data_plot_site[,1]))
  
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
  
  lsBoxFig1[[site]] <- list(gg, gg1, gg2,gg3,gg4)
  gridExtra::grid.arrange(gg, gg1, gg2,gg3,gg4, layout_matrix = rbind(c(1,2,3),c(1,4,5)), widths = c(1.2, 1, 1))
  
  
  ###################
  ### GAM at ARIK ###
  ###################
  
  # Create subset of data for modeling
  bgnTime <- dataSetup$startDate[dataSetup$site == site]
  endTime <- dataSetup$endDate[dataSetup$site == site]

  modl_data <- data_site %>% filter(start_date_time >= bgnTime & start_date_time <= endTime )
  
  #base model with gam package
  model_gam_site<-gam::gam(nitrate_mean ~ spec_cond + oxygen + turbidity + temp_mean, data = modl_data)
  
  if(stepGam){
    # Stepwise selection of variables
    step_model_gam<-gam::step.Gam(model_gam_site, scope=list("spec_cond"=~1+spec_cond+s(spec_cond,4)+s(spec_cond,6)+s(spec_cond,12),
                                                             "oxygen"=~1+oxygen+s(oxygen,4)+s(oxygen,5)+s(oxygen,6)+s(oxygen,12),
                                                             "turbidity"=~1+turbidity+s(turbidity,4)+s(turbidity,5)+s(turbidity,6)+s(turbidity,12),
                                                             "temp_mean"=~1+temp_mean+s(temp_mean,4)+s(temp_mean,6)+s(temp_mean,12)))
    
  }
   # best model selected 
  best_model_gam_site<-mgcv::gam(nitrate_mean ~  s(spec_cond, k = smoothDimK) +
                                   s(oxygen, k = smoothDimK) +  s(temp_mean, k = smoothDimK) +
                                   s(log(turbidity + 1 ), k = smoothDimK) +
                                   s(time_recod, k = smoothDimK) + s(elev, k=smoothDimK),
                                 data = modl_data, na.action=na.exclude)
  
  #fill where we deleted na
  modl_data <- modl_data %>%
    mutate(gam_res = residuals(best_model_gam_site)) %>%
    as_tsibble(index=start_date_time) %>%
    fill_gaps()
  
  
  # Same model for visualization
  best_model_gam_site<-mgcv::gam(nitrate_mean ~  s(spec_cond, k = smoothDimK) + s(oxygen, k = smoothDimK) + 
                                   +  s(temp_mean, k = smoothDimK)+ s(log(turbidity+1),  k = smoothDimK)+ 
                                   s(time_recod, k = smoothDimK) + s(elev,k=smoothDimK), 
                                 data = modl_data)
  
  
  ## RMSE and GCV
  site_res <- c(site, round(best_model_gam_site$gcv.ubre,2))
  site_res
  
  ################
  ### Figure 3 ###
  ################
  a <- getViz(best_model_gam_site)
  
  if(site == "ARIK"){
    # TODO edit the breaks/labels to be exact rather than approximate
    ylabFig_3 <- "s(x)"
    labelsX <- c("Nov-18","Jun-19", "Dec-19")
    breaksX <- c(30000,50000,70000)
  } else if(site == "CARI"){
    ylabFig_3 <- ""
    labelsX <- c("Sept-18","Jul-19")
    breaksX <- c(10000,20000)
  } else if (site == "LEWI"){
    ylabFig_3 <- ""
    labelsX <- c("Jul-18","Feb-19", "Sep-19")
    breaksX <- c(20000,40000,60000)
  }
  
  var1 <- plot(sm(a, 1)) + l_fitLine(colour = "red") + l_ciLine(mul = 5, colour = "blue", linetype = 2) +
    xlab(paste0("Spec. Cond. [",dfUnits$units[dfUnits$var == "SpC"],"]")) + ylab(ylabFig_3)  +
    theme_classic() + 
    ggtitle(siteName)
  var2 <- plot(sm(a, 2)) + l_fitLine(colour = "red") + l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
    xlab(paste0("Diss. Oxygen [",dfUnits$units[dfUnits$var == "DO"],"]")) +ylab(ylabFig_3) + ylim(-5,5) + theme_classic()
  var3 <- plot(sm(a, 3)) + l_fitLine(colour = "red") +  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
    xlab(paste0("Temp [",dfUnits$units[dfUnits$var == "Temp"],"]")) +ylab(ylabFig_3)+ ylim(-5,5) + 
    theme_classic()
  var4 <- plot(sm(a, 4)) + l_fitLine(colour = "red") + l_ciLine(mul = 5, colour = "blue", linetype = 2) +
    xlab(paste0("Log(Turbidity) [",dfUnits$units[dfUnits$var == "turbidity"],"]")) + ylab(ylabFig_3) +
    theme_classic()
  # TODO change to date
  # DateBreaks <- modl_data$start_date_time[1]
  var5 <- plot(sm(a, 5)) + l_fitLine(colour = "red") +  l_ciLine(mul = 5, colour = "blue", linetype = 2) +
    xlab("Time") + ylab(ylabFig_3) +
    theme_classic() + scale_x_continuous(labels= labelsX , breaks = breaksX )
  
  
  var6 <- plot(sm(a,6)) + l_fitLine(colour = "red") + l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
    xlab("Surface Water\nElevation [m]") + ylab(ylabFig_3) +
    theme_classic()
  
  
  
  lsSmoothFig_3[[site]] <- gridPrint(var1, var2, var3, var4,  var5, var6, ncol = 1) 
  
  ############################################
  ### Calculating importance of variables ###
  ###########################################
  
  #Deviance of the GAMM 
  d_null<-deviance(gam(nitrate_mean ~  1 , data = modl_data))
  ar_model<-auto.arima(residuals(best_model_gam_site))  # selecting best ARIMA model
  D_GAMM_total <- deviance(best_model_gam_site)
  D_RES <- sum(residuals(ar_model)^2, na.rm=T)
  Deviance_GAMM_total <- (d_null - D_RES) / d_null
  
  arimaOrderAuto <- forecast::arimaorder(ar_model)
  print(paste0("ARIMA ORDER: ", arimaOrderAuto))
  if(site == "ARIK"){
    
    armaOrder <- c(2,0,5) # 2,0,3 -> manuscript
  } else if (site=="CARI"){
    armaOrder <- c(3,0,2)
    # armaOrder <- c(2,0,5)
  } else if (site == "LEWI"){
    armaOrder <- c(2,0,1)
    # armaOrder <- c(2,0,5)
  }
  df_dev <- data.frame(type = "arma",
             dev_residuals = D_RES,
             dev_null = d_null,
             dev_gamm_totl = Deviance_GAMM_total)
  
  # Deviances of GAMMs minus one by one covariates
  GAM_elev<-mgcv::gam(nitrate_mean ~  s(spec_cond, k = smoothDimK) + s(oxygen, k = smoothDimK) +  s(temp_mean, k = smoothDimK) +
                        s(log(turbidity+1),  k = smoothDimK)+ s(time_recod, k = smoothDimK),
                      data = modl_data, sp=best_model_gam_site$sp[-grep("elev",names(best_model_gam_site$sp))])
  D_GAM <- deviance(GAM_elev)
  D_RES <- sum(residuals(arima(residuals(GAM_elev), order = armaOrder))^2, na.rm=T)
  Deviance_elev <- (D_GAMM_total - D_RES) / D_GAMM_total
  Deviance_elev
  df_dev <- rbind(df_dev,data.frame(type = "elev",
                       dev_residuals = D_RES,
                       dev_null = D_GAMM_total,
                       dev_gamm_totl = Deviance_elev))
  
  GAM_time<-mgcv::gam(nitrate_mean ~  s(spec_cond, k = smoothDimK) + s(oxygen, k = smoothDimK) +  s(temp_mean, k = smoothDimK) +
                        s(log(turbidity+1),  k = smoothDimK)+ s(elev, k = smoothDimK),
                      data = modl_data, sp=best_model_gam_site$sp[-grep("time",names(best_model_gam_site$sp))])
  D_GAM <- deviance(GAM_time)
  D_RES <- sum(residuals(arima(residuals(GAM_time), order = armaOrder))^2, na.rm=T)
  Deviance_time <- (D_GAMM_total - D_RES) / D_GAMM_total
  Deviance_time
  df_dev <- rbind(df_dev,data.frame(type = "time",
                                    dev_residuals = D_RES,
                                    dev_null = D_GAMM_total,
                                    dev_gamm_totl = Deviance_time))
  
  GAM_turbi<-mgcv::gam(nitrate_mean ~  s(spec_cond, k = smoothDimK) + s(oxygen, k = smoothDimK) +  s(temp_mean, k = smoothDimK) +
                         s(time_recod, k = smoothDimK)+ s(elev, k=smoothDimK),# Edit GL, add elev
                         data = modl_data, sp=best_model_gam_site$sp[-grep("turb",names(best_model_gam_site$sp))])
  D_GAM <- deviance(GAM_turbi)
  D_RES <- sum(residuals(arima(residuals(GAM_turbi), order = armaOrder))^2, na.rm=T)
  Deviance_turbi <- (D_GAMM_total - D_RES) / D_GAMM_total
  Deviance_turbi
  df_dev <- rbind(df_dev,data.frame(type = "turbidity",
                                    dev_residuals = D_RES,
                                    dev_null = D_GAMM_total,
                                    dev_gamm_totl = Deviance_turbi))
  
  GAM_temp<-mgcv::gam(nitrate_mean ~   s(spec_cond, k = smoothDimK) +  s(oxygen, k = smoothDimK) + s(log(turbidity+1),  k = smoothDimK)+
                        s(time_recod, k = smoothDimK)+ s(elev, k = smoothDimK),
                      data = modl_data, sp=best_model_gam_site$sp[-grep("temp",names(best_model_gam_site$sp))])
  D_GAM <- deviance(GAM_temp)
  D_RES <- sum(residuals(arima(residuals(GAM_temp), order = armaOrder))^2, na.rm=T)
  Deviance_temp <- (D_GAMM_total - D_RES) / D_GAMM_total
  Deviance_temp
  df_dev <- rbind(df_dev,data.frame(type = "temp",
                                    dev_residuals = D_RES,
                                    dev_null = D_GAMM_total,
                                    dev_gamm_totl = Deviance_temp))
  
  
  GAM_oxygen<-mgcv::gam(nitrate_mean ~   s(spec_cond, k = smoothDimK) +  s(temp_mean, k = smoothDimK) + s(log(turbidity+1),  k = smoothDimK)+s(time_recod, k = smoothDimK)+
                          s(elev, k = smoothDimK), data = modl_data, sp=best_model_gam_site$sp[-grep("oxygen",names(best_model_gam_site$sp))])
  D_GAM <- deviance(GAM_oxygen)
  D_RES <- sum(residuals(arima(residuals(GAM_oxygen), order = armaOrder))^2, na.rm=T)
  Deviance_oxygen <- (D_GAMM_total - D_RES) / D_GAMM_total
  Deviance_oxygen
  df_dev <- rbind(df_dev,data.frame(type = "oxygen",
                                    dev_residuals = D_RES,
                                    dev_null = D_GAMM_total,
                                    dev_gamm_totl = Deviance_oxygen))
  
  
  GAM_cond<-mgcv::gam(nitrate_mean ~   s(oxygen, k = smoothDimK) +  s(temp_mean, k = smoothDimK) + s(log(turbidity+1),  k = smoothDimK)+s(time_recod, k = smoothDimK)+
                        s(elev, k = smoothDimK), data = modl_data, sp=best_model_gam_site$sp[-grep("cond",names(best_model_gam_site$sp))])
  D_GAM <- deviance(GAM_cond)
  D_RES <- sum(residuals(arima(residuals(GAM_cond), order = armaOrder))^2, na.rm=T)
  Deviance_cond <- (D_GAMM_total - D_RES) / D_GAMM_total
  Deviance_cond
  df_dev <- rbind(df_dev,data.frame(type = "cond",
                                    dev_residuals = D_RES,
                                    dev_null = D_GAMM_total,
                                    dev_gamm_totl = Deviance_cond))
  df_dev$site <- site
  lsDevCalc[[site]] <- df_dev
  # deviances at arikaree
  site_dev<-c(Deviance_elev, Deviance_time,Deviance_turbi, Deviance_temp, Deviance_oxygen, Deviance_cond )
  site_imp <- data.frame("Importance" = (1-site_dev)*100, "Var" = c("Elevation", "Time", "Log(turbidity)", "Temperature", "Dissolved Oxygen", "Specific Conductance"))
  lsSiteImp[[site]] <- site_imp
  
  ### aAIC calculation ###
  print(paste0("Length of timeseries: ", nrow(modl_data)))
  aAIC_site <- nrow(modl_data)*log(best_model_gam_site$sig2) + (2*sum(best_model_gam_site$edf)) # for the GAM
  GAMM_site<-forecast::auto.arima(best_model_gam_site$residuals) # 5 df
  aAIC_GAMM_site <- (nrow(modl_data)*log(GAMM_site$sigma2)) + (sum(best_model_gam_site$edf) + 5) # for the GAMM
  
  # Extract autoarima results:
  sumGam <- summary(ar_model)
  sumGamm <- summary(GAMM_site)
  lsArmaRslt <- base::list(site = site,GAM = sumGam, GAMM = sumGamm)
  lsARIMAmodls[[site]] <- base::list(ar_modl = ar_model, GAMM = GAMM_site)
  
  lsGamGammARMASumm[[site]] <- lsArmaRslt
  dfAic <- base::data.frame(site = site, aAIC = aAIC_site, aAIC_GAMM = aAIC_GAMM_site)
  # --------- methods in script_paper_estimates.Rmd
  # summary(best_model_gam_site)

  n <- length(best_model_gam_site$residuals)
  print(paste0("Total residuals at ",site,": ", n))
  if(site == "LEWI"){ #16421
    #n=13341 
    dfAdd <- 3 # TODO what does this mean?
  } else if (site == "ARIK") { #14894
    #n=14973
    dfAdd <- 5
  } else if (site == "CARI") { #8650
   # n=8040 # Number of residuals in model
    dfAdd <- 5
  }
  
  aAIC_site_rmd <- n*log(best_model_gam_site$sig2) + (2*sum(best_model_gam_site$edf))
  
  aAIC_GAMM_site_rmd <-(n*log(GAMM_site$sigma2)) + (sum(best_model_gam_site$edf) + dfAdd)
  
  dfAic_rmd <- base::data.frame(aAICrmd = aAIC_site_rmd, aAIC_GAMM = aAIC_GAMM_site_rmd)
  dfAic <- cbind(dfAic, dfAic_rmd)

  # --------- END methods in script_paper_estimates.Rmd
  # Compile Gam Results

  
  lsAic[[site]] <- dfAic
  
  ################
  ### FIGURE 4 ###
  ################
  
  total_imp<- site_imp
  total_imp$site<-siteName
  lsImport[[site]] <- total_imp
  
  lsImportFig4[[site]] <- total_imp %>%
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
  
  
  ###############
  #### SI 2 ####
  ##############
  # gridExtra::grid.arrange(g1, g2,g3, layout_matrix = rbind(c(1,2,3)), widths = c(1, 1, 1))
}

# Combine all results of water quality characteristics:
dtRsltWaq <- data.table::rbindlist(lsRsltsWaqChars)
write.csv(dtRsltWaq,file.path(plot_dir, "WaterQualityCharacteristicsAmongSites.csv"))

# Combine and write AIC results:
dtAic <- data.table::rbindlist(lsAic)
write.csv(dtAic, file.path(plot_dir,"aAIC_results_GAM_GAMM.csv"))

# Combine site results for calculations for variable importance
dtDevCalc <- data.table::rbindlist(lsDevCalc)
write.csv(dtDevCalc,file.path(plot_dir, "deviance_calcs_GAMM.csv"))

# Combine panels to create Figure 2:
plotAllFig2 <- mgcViz::gridPrint(lsFig2[["ARIK"]],lsFig2[["CARI"]],lsFig2[["LEWI"]], ncol=3)
ggplot2::ggsave(plot=plotAllFig2,
                filename=file.path(plot_dir,"Fig_2_all_sites.png"),
                width = 10, height = 8)
ggplot2::ggsave(plot=plotAllFig2,
                filename=file.path(plot_dir,"Fig_2_all_sites.svg"),
                width = 10, height = 8)
# Combine panels to create Fig 3
plotAllFig_3 <- mgcViz::gridPrint(lsSmoothFig_3[["ARIK"]],lsSmoothFig_3[["CARI"]],lsSmoothFig_3[["LEWI"]],ncol=3)
ggplot2::ggsave(plot=plotAllFig_3,
                filename=file.path(plot_dir,"Fig_3_all_sites.png"),
                width = 10, height = 8)
ggplot2::ggsave(plot=plotAllFig_3,
                filename=file.path(plot_dir,"Fig_3_all_sites.svg"),
                width = 10, height = 8)

# Combine data to create Fig 4:
lsAllSiteImp <- lapply(names(lsSiteImp), function(n) mutate(lsSiteImp[[n]], site = n) )
dtAllSiteImp <- data.table::rbindlist(lsAllSiteImp)
dtAllSiteImp1 <- merge(dtAllSiteImp,siteNameDf,by.x = "site", by.y = "siteID") 

plotImportFig4 <- dtAllSiteImp1 %>%
  ggplot(aes(x = Importance, y = Var, color = siteName, shape = siteName)) +
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
ggplot2::ggsave(plot=plotImportFig4,
                filename=file.path(plot_dir, "Fig_4_all_sites_importance.png"),
                width=6,height = 6, units = "in")
ggplot2::ggsave(plot=plotImportFig4,
                filename=file.path(plot_dir, "Fig_4_all_sites_importance.svg"),
                width=6,height=6, units = "in")


# Print results and save to file:
sink(file.path(plot_dir,"ARMA_model_summaries.txt"))
for(sn in names(lsGamGammARMASumm)){
  print(sn)
  print("ARMA MODEL, GAM residuals:\n")
  summary(lsARIMAmodls[[sn]]$ar_modl)
  print("ARMA MODEL, GAMM residuals: \n")
  summary(lsARIMAmodls[[sn]]$GAMM)
}
sink(file=NULL)