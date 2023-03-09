OS = Sys.info()['sysname']
PROGNAME =	"metals_imputation_final_model_20iter_01262022.R"
if(tolower(OS)=="windows"){
  sink(paste0("M:/WELLWISE/logs/", PROGNAME, "out"), split=TRUE)
  BASE = "M:/WELLWISE/"
}

if(tolower(OS)=="darwin"){
  if(dir.exists("/Users/akeil/temp/mdrive/m/WELLWISE/share/")){
    BASE = "/Users/akeil/temp/mdrive/m/WELLWISE/share/"
  } else{
    BASE = "/Volumes/share/"
  }
  sink(paste0(BASE, "logs/", PROGNAME, "out"), split=TRUE)
}

print(Sys.time())
###############################################################################-
#
# Well water metals imputation
# Purpose: Impute metals via model chosen in metals_imputation_final_model_08192024.R
# Reference: mice.impute.leftcenslognorm in 'qgcomp' documentation (Author: Alex Keil)
#
# Author: alex Keil
# Date: 8/24/2021
# changes from previous version
#  implement parallel processing
#  keep census tract with data
#  fixed bug in logxform matrix
#
###############################################################################-
library("qgcomp")
library("ggplot2")
library("mice")
library("haven")
library("dplyr")
library("readr")
library("survival")
#library("tidyverse")
#library("knitr")
options(scipen=999)


# sets directory
if(Sys.info()[["sysname"]]=="Linux"){
  # running on longleaf
  mdrive = "/proj/keil_lab/projects/imputation/"
  indir = paste0(mdrive, "/data/")
  outdir = paste0(mdrive, "/output/")
}
if(Sys.info()[["sysname"]]=="Darwin"){
  mdrive = "/Users/akeil/temp/mdrive/m/WELLWISE/share"
  indir = paste0(mdrive, "/projects/multiple_imputation/data/")
  outdir = paste0(mdrive, "/projects/multiple_imputation/output/")
  #outdir = paste0(mdrive, "/user_Butts/output/final/300 it/")
  #indir = paste0(mdrive, "/projects/multiple_imputation/data/")
}
if(Sys.info()[["sysname"]]=="Windows"){
  mdrive = "M:/WELLWISE/"
  indir = paste0(mdrive, "/projects/multiple_imputation/data/")
  outdir = paste0(mdrive, "/projects/multiple_imputation/output/")
}


setwd(indir)
source(paste0(indir,"../code/logxform_micefun2.R"))
#setwd("M:/projects/multiple_imputation/data")

compoundlist <- c(
  "Chromium", "Copper","Mercury", # added 1/22
  "Manganese", "Arsenic", "Lead", "Cadmium", "Calcium", "Chloride", "Iron", "Magnesium", "Total_Hardness",
  "Total_Alkalinity", "Sodium", "Zinc", "Sulfate"
)
lodlist <- paste0(compoundlist, "_lod")


#variables for imputation
keepcols = c("sample_id", "MMDDYY", "county", "census_tract", compoundlist, lodlist
)
#"Total_Alkalinity", "Sodium", "Zinc", "Sulfate"

# current well water data (will get replaced by permanent version)
# this is hard coded to speed up the read, but will need to be replaced if file structure changes
metalcols <- readr::cols(
  sample_id = col_character(),
  sample_date = col_character(),
  MMDDYY = col_character(),
  longitude = col_double(),
  latitude = col_double(),
  census_tract = col_double(),
  county = col_character(),
  Arsenic = col_double(),
  Arsenic_lod = col_double(),
  Lead = col_double(),
  Lead_lod = col_double(),
  Manganese = col_double(),
  Manganese_lod = col_double(),
  Cadmium = col_double(),
  Cadmium_lod = col_double(),
  Acidity = col_double(),
  Aluminum = col_double(),
  Antimony = col_double(),
  Barium = col_double(),
  Beryllium = col_double(),
  Boron = col_double(),
  Calcium = col_double(),
  Chloride = col_double(),
  Chromium = col_double(),
  Cobalt = col_double(),
  Color = col_double(),
  Conductivity = col_double(),
  Copper = col_double(),
  Cyanide = col_double(),
  Fluoride = col_double(),
  Gold = col_double(),
  Hexavalent_Chromium = col_double(),
  Insoluble_Iron = col_double(),
  Insoluble_Manganese = col_double(),
  Iron = col_double(),
  Lithium = col_double(),
  Magnesium = col_double(),
  Mercury = col_double(),
  Molybdenum = col_double(),
  Nickel = col_double(),
  Nitrate = col_double(),
  Nitrite = col_double(),
  Orthophosphate = col_double(),
  pH = col_double(),
  Potassium = col_double(),
  Selenium = col_double(),
  Settleable_Solids = col_double(),
  Silica = col_double(),
  Silver = col_double(),
  Sodium = col_double(),
  Soluble_Iron = col_double(),
  Soluble_Manganese = col_double(),
  Strontium = col_double(),
  Sulfate = col_double(),
  Thallium = col_double(),
  Tin = col_double(),
  Total_Alkalinity = col_double(),
  Total_Dissolved_Solids = col_double(),
  Total_Hardness = col_double(),
  Total_Phosphate = col_double(),
  Total_Suspended_Solids = col_double(),
  Turbidity = col_double(),
  Uranium = col_double(),
  Vanadium = col_double(),
  Zinc = col_double(),
  Acidity_lod = col_double(),
  Aluminum_lod = col_double(),
  Antimony_lod = col_double(),
  Barium_lod = col_double(),
  Beryllium_lod = col_double(),
  Boron_lod = col_double(),
  Calcium_lod = col_double(),
  Chloride_lod = col_double(),
  Chromium_lod = col_double(),
  Cobalt_lod = col_double(),
  Color_lod = col_double(),
  Conductivity_lod = col_double(),
  Copper_lod = col_double(),
  Cyanide_lod = col_double(),
  Fluoride_lod = col_double(),
  Gold_lod = col_double(),
  Hexavalent_Chromium_lod = col_double(),
  Insoluble_Iron_lod = col_double(),
  Insoluble_Manganese_lod = col_double(),
  Iron_lod = col_double(),
  Lithium_lod = col_double(),
  Magnesium_lod = col_double(),
  Mercury_lod = col_double(),
  Molybdenum_lod = col_double(),
  Nickel_lod = col_double(),
  Nitrate_lod = col_double(),
  Nitrite_lod = col_double(),
  Orthophosphate_lod = col_double(),
  pH_lod = col_double(),
  Potassium_lod = col_double(),
  Selenium_lod = col_double(),
  Settleable_Solids_lod = col_double(),
  Silica_lod = col_double(),
  Silver_lod = col_double(),
  Sodium_lod = col_double(),
  Soluble_Iron_lod = col_double(),
  Soluble_Manganese_lod = col_double(),
  Strontium_lod = col_double(),
  Sulfate_lod = col_double(),
  Thallium_lod = col_double(),
  Tin_lod = col_double(),
  Total_Alkalinity_lod = col_double(),
  Total_Dissolved_Solids_lod = col_double(),
  Total_Hardness_lod = col_double(),
  Total_Phosphate_lod = col_double(),
  Total_Suspended_Solids_lod = col_double(),
  Turbidity_lod = col_double(),
  Uranium_lod = col_double(),
  Vanadium_lod = col_double(),
  Zinc_lod = col_double()
)

#cdat_1 = readr::read_csv("wellwise_NCwellmetals_1998_2019_20200825.csv", guess_max = 2,
#                         col_types = metalcols)       #reads in contaminant data
#cdat_1 = readr::read_csv("wellwise_NCwellmetals_1998_2019_20210824.csv", guess_max = 2,
#                         col_types = metalcols)       #reads in contaminant data
cdat_1 = readr::read_csv("wellwise_NCwellmetals_1998_2019_20220127.csv", guess_max = 2,
                         col_types = metalcols)       #reads in contaminant data
setwd(outdir)
cdat_2 = cdat_1[, keepcols]                                                                       #keep contaminant data for imputation


cdat_2$Year <- substr(cdat_2$MMDDYY,7,8)                                                          #creates year variable (for imputing missing Mn LOD)
cdat_2$date <- as.Date(cdat_2$MMDDYY, format = "%m/%d/%y")                                                          #creates year variable (for imputing missing Mn LOD)
cdat_2$month_year = format(cdat_2$date, "%Y%m")

# a new variable that is equal to the LOD for manganese
# manganese is used because it has lots of non-missing >LOD values and is useful for judging accuracy of imputations
# NOTE: CORRECTED A BUG IN  VALUE FOR TOTAL HARDNESS
# NOTE: CORRECTED SORTING BUG 1/13 - ak
# NOTE: set infinite lods to 0 1/29 - ak
real_lod_value <- function(metal, lodvar){
  ym = cdat_2$month_year
  cdat_2$ord = 1:nrow(cdat_2)
  ymvals = sort(unique(ym))
  suppressWarnings(lods <- tapply(metal, ym, function(x) min(x[x>0], na.rm=TRUE)))
  #lods[is.infinite(lods)] = NA
  lods[is.infinite(lods)] = 0
  #lods = lods[ & !is.na(lods)]
  df1 = data.frame(lod=lods, month_year=names(lods))
  anymiss = sum(is.na(df1$lod)) + length(setdiff(df1$month_year, ymvals)) # finds elements of first set that are not in second set
  if(anymiss>0){
    # carry forward imputation if there are any months with no LOD
    #df1 = merge(data.frame(month_year=ymvals), df1, all.x=TRUE, all.y = FALSE)
    df1 = merge(data.frame(month_year=ymvals), df1, all.x=TRUE, all.y = FALSE, sort=TRUE)
    isna = is.na(df1$lod)
    isinf = is.infinite(df1$lod)
    if(isna[1] | isinf[1]) df1$lod[1] = df1$lod[!isna & !isinf ][1]  # carry backward if first value is missing
    for(row in 2:nrow(df1)){
      if(isna[row] | isinf[row]) df1$lod[row] = df1$lod[row-1]
    }
  }
  #df2 = merge(cdat_2[,c("sample_id","month_year"),drop=FALSE], df1, all.x=TRUE, all.y=FALSE,by="month_year")
  df2 = merge(cdat_2[,c("ord","month_year"),drop=FALSE], df1, all.x=TRUE, all.y=FALSE,by="month_year", sort=FALSE)
  df2 = df2[ order(df2$ord), ]# 1/13 - ak
  df2$lod
}

#apply to all exposures
lodvaluelist <- paste0(lodlist, "_value")
for(i in 1:length(compoundlist)){
  cdat_2[,lodvaluelist[i]] <- real_lod_value(metal = cdat_2[,compoundlist[i], drop=TRUE], lodvar=cdat_2[,lodlist[i], drop=TRUE])
}


#cdat_3 <- subset(cdat_2, county == "Orange" | county == "Mecklenburg" | county == "Union")       #restricts dataset to Orange, Mecklenburg, and Union counties
cdat_4 <- cdat_2 %>%                                                                              #excludes observations with missing values for all key contaminants
  #dplyr::filter(!(is.na(Arsenic) & is.na(Cadmium) & is.na(Lead) & is.na(Manganese)))
  mutate()

# need indicator to tell MICE to impute missing values (in these data we'll have at least one non-missing metal per obs)
#create indicator LOD variable
#creates new LOD indicator variable for contaminants with 1: missing values or values below LOD, 0: non-missing and above LOD
pseudo_lod_indicator <- function(lodindicator){
  print(paste(sum(is.na(lodindicator)), "missing values"))
  ifelse(is.na(lodindicator), 1, lodindicator)
}

#apply to all exposures
pseudolodlist <- paste0(compoundlist, "_pseudolod")
for(i in 1:length(compoundlist)){
  cdat_4[,pseudolodlist[i]] <- pseudo_lod_indicator(lodindicator=cdat_4[,lodlist[i], drop=TRUE])
}


#resets LOD values                                                                                #replaces LOD indicator variable with interval scale measures; missing LOD values are set to the maximum observed value
# coding trick: mice will impute any value <LOD, so for true missing values, set LOD to high value and mice will just impute a nice value
#   E.G.:
#     cdat_4$Arsenic_lodvalue <- ifelse(cdat_4$Arsenic_lod == 1, cdat_4$Arsenic, cdat_4$Arsenic_lod)
#     cdat_4$Arsenic_lodvalue <- ifelse(is.na(cdat_4$Arsenic_lod), max(cdat_4$Arsenic, na.rm=TRUE), cdat_4$Arsenic_lod) # keep imputed values below the max observed value
pseudo_lod_value <- function(metal, lodvalue, bigfakelodvalue=NULL){
  #lodvalue = cdat_4$Manganese_lod_value
  #metal = cdat_4$Manganese
  # if metal is missing, use really high value
  # if metal is non-missing (either <LOD or >=LOD), use actual lod
  if(is.null(bigfakelodvalue)) bigfakelodvalue <- max(metal, na.rm=TRUE)
  ifelse(!is.na(metal), lodvalue, bigfakelodvalue)
}


#apply to all exposures
pseudolodvaluelist <- paste0(pseudolodlist, "_value")
for(i in 1:length(compoundlist)){
  cdat_4[,pseudolodvaluelist[i]] <- pseudo_lod_value(cdat_4[,compoundlist[i], drop=TRUE], cdat_4[,lodvaluelist[i], drop=TRUE])
}


#sets measurements below LOD to missing for imputation
# this tells 'mice' to impute new values
setmissing_belowlod <- function(metal, lod){
  #ifelse(cdat_4$Arsenic > cdat_4$Arsenic_lod, cdat_4$Arsenic, NA)
  whichlod  = (metal <= lod)
  metal[whichlod] <- NA
  metal
}


nalist <- paste0(compoundlist, "_na")
natlist <- paste0(nalist, "t")
for(i in 1:length(compoundlist)){
  # na variables are missing below lod
  # nat variables are missing below lod and always missing in test set (if applicable)
  cdat_4[,nalist[i]] <-
    cdat_4[,natlist[i]] <-
    setmissing_belowlod(cdat_4[,compoundlist[i], drop=TRUE]       , cdat_4[,lodvaluelist[i], drop=TRUE])
}


#create new county variable "county2" for use in models when county is grouped; counties with 0 Cadmium values above LOD are grouped
# NOTE: we want to use "county" to impute values that should be similar due to proximity, but using all counties is problematic
Cd_lod <- cdat_4[, c("Cadmium_lod", "county")]
Cd_lod$Cadmium_gtLOD <- ifelse(!is.na(Cd_lod$Cadmium_lod), 1-Cd_lod$Cadmium_lod, 0)

Hg_lod <- cdat_4[, c("Mercury_lod", "county")]
Hg_lod$Mercury_gtLOD <- ifelse(!is.na(Hg_lod$Mercury_lod), 1-Hg_lod$Mercury_lod, 0)


Cd_sum <- aggregate(Cd_lod$Cadmium_gtLOD, by =list(Category =Cd_lod$county), FUN = sum)
Cd_sum$grouping_ind <- ifelse(Cd_sum$x <= 1 , 1, 0) #55/100 counties
Cd_sum_grouping <- Cd_sum[,c(1,3)]


cdat_4 <- merge(cdat_4, Cd_sum_grouping, by.x = "county", by.y = "Category")
cdat_4$countygrp <- ifelse(cdat_4$grouping_ind == 1, "grouped", cdat_4$county)


Hg_sum <- aggregate(Hg_lod$Mercury_gtLOD, by =list(Category2 =Hg_lod$county), FUN = sum)
Hg_sum$grouping_ind2 <- ifelse(Hg_sum$x <= 1 , 1, 0) #76/100 counties
Hg_sum_grouping <- Hg_sum[,c(1,3)]
cdat_4 <- merge(cdat_4, Hg_sum_grouping, by.x = "county", by.y = "Category2")
cdat_4$countygrp2 <- ifelse(cdat_4$grouping_ind2 == 1, "grouped", cdat_4$county)

cdat_4$county1 = as.factor(c(cdat_4[,"county", drop=TRUE]))
cdat_4$county2 = as.factor(c(cdat_4[,"countygrp", drop=TRUE]))  # Valid for Cadmium and all others except mercury
cdat_4$county3 = as.factor(c(cdat_4[,"countygrp2", drop=TRUE])) # valid for Mercury


## Imputation data ----
#mdat <- data.frame(Manganese = cdat_4$Manganese_na, Manganese_pseudolod = cdat_4$Manganese_pseudolod, Manganese_pseudolod_value = cdat_4$Manganese_pseudolod_value, county = county)
mdat3 <- cdat_4


init_vals <- function(data, nalist, lodlist, lodvaluelist){
  nm = length(nalist)
  for(i in 1:nm){
    idx = which(data[,lodlist[i]] == 1)
    data[idx,nalist[i]] = data[idx,lodvaluelist[i]] / sqrt(2)
    idx2 = which(is.na(data[,nalist[i]]))
    data[idx2,nalist[i]] = sample(data[-idx2,nalist[i]], size = length(data[idx2,nalist[i]]), replace = TRUE)
  }
  data
}


###########-
## Imputation model creation function ----
###########-
impsetup_final <- function(
  metals,
  extravars,
  dat=mdat3,
  metals2 = NULL,
  extravars2 = NULL,
  idvars = c("sample_id", "census_tract")
){
  #idvars = c("sample_id")
  allmetals = c(metals, metals2)
  allextravars = c(extravars, extravars2)
  metalnames = paste0(metals, "_na")
  metalnames2 =  paste0(metals2, "_na")
  allmetalnames = paste0(allmetals, "_na")
  #
  datainit = init_vals(dat, nalist, lodlist, lodvaluelist)
  lodindlist = list()
  for(idvar in idvars){
    lodindlist[[idvar]] <- rep(NA,nrow(dat))
  }
  for(m in allmetals){
    lodindlist[[m]] <- dat[,paste0(m, "_pseudolod_value")]
  }
  for(v in allextravars){
    lodindlist[[v]] <- rep(NA,nrow(dat))
  }
  nvars = length(c(allmetals, allextravars))
  #logxformmat = do.call("rbind", lapply(1:nvars, function(x) c(rep(1, length(allmetals)), rep(0, length(allextravars)))))
  logxformmat = do.call("rbind", lapply(1:(nvars+length(idvars)), function(x) c(rep(0, length(idvars)), rep(1, length(allmetals)), rep(0, length(allextravars)))))
  predmat =     do.call("rbind", lapply(1:(nvars+length(idvars)), function(x) c(rep(0, length(idvars)), rep(1, nvars))))
  colnames(predmat) <- rownames(predmat) <- colnames(logxformmat) <- rownames(logxformmat) <- c(idvars, allmetalnames, allextravars)
  # dont use self to predict
  diag(predmat) <- 0
  # dont impute id variables
  predmat[which(rownames(predmat) %in% c(idvars)),] <- 0
  # dont use extravars 2 to predict metals 1
  predmat[which(rownames(predmat) %in% metalnames), which(colnames(predmat) %in% extravars2)] <- 0
  # dont use extravars to predict metals 2
  predmat[which(rownames(predmat) %in% metalnames2), which(colnames(predmat) %in% extravars)] <- 0
  #dimnames(predmat) <- NULL # errors without this step due to name conflicts
  # dropping rows wont work unless specifying blocs
  #predmat = predmat[-which(rownames(predmat) %in% c(idvars,extravars,extravars2)),]


  cat(paste("\n Using", paste(allmetalnames, collapse=",")))
  list(lodindlist=lodindlist,
       midat = dat[,c(idvars, allmetalnames, allextravars)],
       predmat = predmat,
       logxformmat = logxformmat,
       methodstring = c(rep("", length(idvars)), rep("leftcenslognorm", length(allmetals)), rep("", length(allextravars))),
       inits = datainit
       )
}



############-
## Imputation model 11b: final model ----
print("Model 11b")
#source(paste0(outdir,"../code/logxform_micefun.R"))
library("qgcomp")
mlst = impsetup_final(
  metals=c(                                 # use only other metals and extravars to predict
    "Calcium", "Total_Alkalinity","Sodium",
    "Total_Hardness", "Magnesium",
    "Manganese", "Chloride", "Iron", "Zinc", "Sulfate",
    "Copper", # added 1/22
    "Chromium",
    "Arsenic",  "Lead" #, "Mercury"  # dropped due to collinearity
  ),
  extravars = c("county1"),
  dat=mdat3,
  metals2=c("Cadmium"),                       # use only other metals and extravars2 to predict
  extravars2=c("county2"),
  idvars = c("sample_id", "census_tract")
)
#set.seed(1231)
set.seed(1231123) # actually does nothing with parallel implementation
mlst$methodstring = gsub("leftcenslognorm", "leftcenslognormlogX", mlst$methodstring)
#mlst$methodstring[7] <- "pmm"
#mlst$methodstring = gsub("leftcenslognorm", "pmm", mlst$methodstring)
seedval <- 1231123
numsets <- 20 # 20, total number of imputations
warmup <- 300 # 300, total number of imputations
#impdat_final = mice(data = mlst$midat,
impdat_final = parlmice(data = mlst$midat,         # changed 1/22 parallel processing
                    method = mlst$methodstring,
                    predictorMatrix = mlst$predmat,
                    logxformMatrix = mlst$logxformmat,
                    lod=mlst$lodlist,
                    debug=FALSE,
                    data.init = mlst$inits,
                    #m=numsets, maxit = warmup,
                    # changed 1/22 parallel processing
                    cl.type=ifelse(Sys.info()[["sysname"]]=="Windows","PSOCK", "FORK"),
                    n.core = numsets, n.imp.core = warmup, seed = ifelse(numsets>1, NA, seedval)
)
save(impdat_final, file = paste0(outdir, "mi_final_20_iter.RData"))


############-
## Imputation model performance (TODO) ----
# TODO: check stability across the 20 imputed data sets
#  - do this on inividual basis (take 30 random wells, look at SD of imputations across 20 imputed datasets)
#  - also on census tract averages: estimate SD for census tract means across 20 imputed datasets
# Tasks: 1) use 'complete' function for all 20 datasets: e.g. complete(impdat_final, 12) is the 12th imputed dataset
#        2) get census tract means for for Mn, As, Cd, Pb for all 20 datasets (one of these will be merged with defects data: #1)
#        3) look at performance metrics above
#        4) move final program and data into M:/Projects/multiple_imputation folder

complete(impdat_final, 1)
complete(impdat_final, 2)
complete(impdat_final, 3)
complete(impdat_final, 4)
complete(impdat_final, 5)
complete(impdat_final, 6)
complete(impdat_final, 7)
complete(impdat_final, 8)
complete(impdat_final, 9)
complete(impdat_final, 10)
complete(impdat_final, 11)
complete(impdat_final, 12)
complete(impdat_final, 13)
complete(impdat_final, 14)
complete(impdat_final, 15)
complete(impdat_final, 16)
complete(impdat_final, 17)
complete(impdat_final, 18)
complete(impdat_final, 19)
complete(impdat_final, 20)
dat2 <- complete(impdat_final)

#save imputed data as csv file
write.csv(dat2, paste0(outdir, "impdataset_m1.csv"), row.names = TRUE)

#plot mean and sd
plot(impdat_final)

plot(impdat_final, layout = c(2, 1))

#obtain mean and sd
t(as.data.frame(impdat_final$chainMean))
means <- t(as.data.frame(impdat_final$chainMean))

t(as.data.frame(impdat_final$chainVar))
vars <- t(as.data.frame(impdat_final$chainVar))

#save means and sd as csv files
write.csv(means,paste0(outdir, "means_m1.csv"), row.names = TRUE)
write.csv(vars,paste0(outdir, "var_m1.csv"), row.names = TRUE)

sink()
