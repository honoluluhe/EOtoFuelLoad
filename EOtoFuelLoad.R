#title: "EOtoFuelLoad"
#Author: Lulu He (lulu.he@adelaide.edu.au)
#Date: 08/08/2023




##########################################################################
###Step 1: IdentiFy and obtain relevant earth observation products
#The EO product used to identify burned areas (BAs) is “MCD64A1 Version 6 global BA” (called MODIS hereafter), a burn mask product at 500 m × 500 m resolution (Giglio et al., 2020). 
#The burn mask is generated through compositing the daily surface reflectance and daily active fire masks (Giglio et al., 2018; Vermote et al., 2002). 
#This product is selected because it has a daily temporal resolution (Giglio et al., 2018), low occurrence of false positives (Hantson et al., 2013), high accuracy of BA mapping 
#compared with other products, such as MCD45, Geoland 2 and Fire_cci (Padilla et al., 2015; Vermote et al., 2002), and covers a relatively long temporal archive from 2000 to 2022 (Giglio et al., 2020), 
#resulting in the 21 year analysis period mentioned above. Although this product has been discontinued, the product and its underlying algorithm provide the baseline for its replacements, 
#such as the Visible Infrared Imaging Radiometer Suite (VIIRS) VNP64A1 BA product (Giglio et al., 2018). MODIS contains information on the burn date, information on the  uncertainty of the burn date, 
#as well as an indicator of quality assurance (Giglio et al., 2015). Images covering the study area are downloaded in WGS84 projection from the USGS website (USGS, 2022), with a total cumulative BA of 
#approximately 3.6 million hectares for the 21-year study period.




##########################################################################
###Step 2: Generate fuel disturbance event maps using MODIS###
#MODIS burnt date file path
AppEEARS.MODIS.MCD64.data.dir <- file.path("C:/Users/a1677880/Earth Observation data pre-processing/Time since last disturbance/WA MCD64 raw data/WA")

#do fire year analysis, call in images from July in a given year to June (inclusive) in the following year. MODIS image has the first day of the Julian month in its file name, which can be used to extract the image. First, define the date vector for leap and non-leap year separately.
#define the date vector for non-leap year
dates.vec.nonleapYear.p1 <- c("182", "213", "244", "274", "305", "335") #first half of the financial year
dates.vec.nonleapYear.p2 <- c("001", "032", "060", "091", "121", "152") #second half of the financial year

#repeat the above step to leap year: add one day
dates.vec.leapYear.p1 <- c("183", "214", "245", "275", "306", "336") #first half of the financial year
dates.vec.leapYear.p2 <- c("001", "032", "061", "092", "122", "153") #second half of the financial year

#classify the year types and define the vector of years. Each type 
years.before.leapyear <- c(2003, 2007, 2011, 2015, 2019) #the years before leap year
leap.years <- c(2004, 2008, 2012, 2016, 2020) #leap years
years.after.leapyear <- c(2005, 2009, 2013, 2017, 2021, 2001, 2002, 2006, 2010,2014,2018, 2022) #all other years
#OtherYears <- c(2001, 2002, 2006, 2010,2014,2018, 2022)

#setting year for testing purpose, remove before running for loop
#year <- 2007

start = Sys.time() 
for (year in 2001:2021) {
  
  #read in the MODIS images for a fire year based on the year type
  if (year %in% years.before.leapyear) {
    
    rast1 <- list.files(file.path(AppEEARS.MODIS.MCD64.data.dir, year), pattern = paste0("^WA_MCD.+", dates.vec.nonleapYear.p1, collapse = "|"), full.names = T)
    rast2 <- list.files(file.path(AppEEARS.MODIS.MCD64.data.dir, (year + 1)), pattern = paste0("^WA_MCD.+", dates.vec.leapYear.p2, collapse = "|"), full.names = T)
    
  } else if (year %in% leap.years){
    rast1 <- list.files(file.path(AppEEARS.MODIS.MCD64.data.dir, year), pattern = paste0("^WA_MCD.+", dates.vec.leapYear.p1, collapse = "|"), full.names = T)
    rast2 <- list.files(file.path(AppEEARS.MODIS.MCD64.data.dir, (year + 1)), pattern = paste0("^WA_MCD.+", dates.vec.nonleapYear.p2, collapse = "|"), full.names = T)
    
  } else if (year %in% years.after.leapyear) {
    rast1 <- list.files(file.path (AppEEARS.MODIS.MCD64.data.dir, year), pattern = paste0("^WA_MCD.+", dates.vec.nonleapYear.p1, collapse = "|"), full.names = T)
    rast2 <- list.files(file.path (AppEEARS.MODIS.MCD64.data.dir, (year + 1)), pattern = paste0("^WA_MCD.+", dates.vec.nonleapYear.p2, collapse = "|"), full.names = T)
  } 
  
  #combine the two parts of raster images into one (= a given financial year)
  MODIS.financialyear.files <- c(rast1, rast2)
  MODIS.financialyear.files <- lapply(MODIS.financialyear.files, rast)
  
  #read in the 12 raster images in a fire year
  MODIS.financialyear.rast <- rast(MODIS.financialyear.files)
  
  #replace the cell values that are not zero to 1 (= the cell has burnt)
  MODIS.financialyear.rast <- subst(MODIS.financialyear.rast, from = 1:366, to = 1)
  
  #replace the cell values that zero and below zero (= the cell has NOT burnt)
  MODIS.financialyear.rast <- subst(MODIS.financialyear.rast, from = -2:0, to = 0)
  
  #sum all the layers together to create a single raster layer dataset
  MODIS.financialyear.rast <- app(MODIS.financialyear.rast, fun = "sum")
  
  #reproject to match the projection of the initial Fuel.Age layer
  MODIS.financialyear.rast <- terra::project(MODIS.financialyear.rast, Base.CRS, method = "near") #this line doesn't work, fixing it by re-open R
  
  #resample to match the resolution of the template
  MODIS.financialyear.rast <- terra::resample(MODIS.financialyear.rast, Template, method = "near")
  
  # Reclassify the cell values according to the plan below
  # If, cell value ≥ 1 => 1 (= fuel age for the year is reset to zero)
  # If, cell value = 0 => NA 
  MODIS.financialyear.rast <- ifel(MODIS.financialyear.rast >= 1, 1, 0)
  
  #mask out water body using LU mask
  crs(MODIS.financialyear.rast) <- crs(LU.mask) # set the two objects as the same crs
  MODIS.financialyear.rast <- mask(MODIS.financialyear.rast, LU.mask, updatevalue = NA) #as the LU mask is already use the template, so skip the crop step below
  
  #view the result images and count if there are 21 images (should be 21 images + a blank image)
  plot(MODIS.financialyear.rast)
  
  #export for later use
  writeRaster(MODIS.financialyear.rast, file.path(Output.file, glue("MODIS_burnt_area_financial_{year}.tif")), filetype = "GTiff", datatype = 'FLT4S', overwrite = TRUE)
  cat(paste("Done with year ", year))
  cat("\n")
  
}
Sys.time() - start 


##########################################################################
###Step 3: Obtain current (i.e. 2001) fuel age map using MCFH#############
#To initialise the calculation, a baseline fuel age map for the year 2000 is required, as the annual data for MODIS is only available from 2001 onwards, as mentioned previously. 
#This baseline fuel age map is obtained with the aid of a manually compiled fire history dataset (MCFH) (see Section 3.3 for details), from which we compute baseline fuel age from 
#1975 to align with the calculation of the Fire Behaviour Index (FBI) in the Australian Fire Danger Rating System (Matthews, 2022). Consequently, fuel age starts accumulating from 
#1975, resulting in a maximum fuel age of 25 years on the baseline fuel age map for areas unburned until 2000.

# Load Fire History Layers
Fire.History.raw <- Fire.scar # takes ?min for WA

# crop to the extent of the study area
Fire.History <- terra::crop(Fire.History.raw, terra::ext(terra::project(Template, terra::crs(Fire.History.raw)))+0.25) %>% terra::project(Base.CRS)

#terra::plot(Fire.History)

#fire scar temporal extent is from 1937 to present, but we only need till 2001 as MODIS EO data start in 2001

####################

year <- 2000

# Use maximum Fuel age starting at 10 years to give some visual progression of how increasing fuel age affects risk.
# This is an assumption which will affect how it "looks". Not necessarily made to be accurate, made to demonstrate value of temporal dynamics.
# Background.TSLF <-  10.0 + step_year
Background.TSLF <-  year - 1975 #Lulu:the SA fire scar data start in 1931, WA fire scar data start in 1937. # set 1975 as the starting because AFDRS calculation is using 25 years as background value (source: email contact with Stuart AFDRS)

# subset only the fire polygons which occurred before 2001:
Fire.History <- terra::subset(Fire.History, Fire.History$fih_date1 < paste0(year, "0101") & Fire.History$fih_date1 > paste0("19750101")) 

# Add a field for number of years between historical fires and the fire of interest
Fire.History$FuelAge <- lubridate::time_length(interval(paste0(year, "0101"), start = Fire.History$fih_date1), unit = "year")

# round the Fuel age to the nearest digit
Fire.History$FuelAge <- round(Fire.History$FuelAge, 0)

# Rasterize the fuel age layer, using the Background value as that defined above. 
Fuel.Age <- terra::rasterize(Fire.History, Template, field = "FuelAge", background = Background.TSLF)

# plot the results
terra::plot(Fuel.Age)

# export the fuel age layer for later use
writeRaster(Fuel.Age, file.path(MODIS.fuelage.output, glue("MODIS_FuelAge_Accumulation_{year}.tif")), filetype = "GTiff", datatype = 'INT2U', overwrite = TRUE)

fuelage2000 <- rast(file.path(MODIS.fuelage.accumulation.output, "Fuel_Age_2000.tif"))



##########################################################################
#### Step 4: generate fuel age maps using MODIS burned area maps#####
start = Sys.time()

#year <- 2001

for (year in 2001:2021) {
  #Load the initial Fuel Age map
  FuelAge <- rast(file.path(MODIS.fuelage.output, glue("MODIS_FuelAge_Accumulation_{year - 1}.tif")))
  #plot(FuelAge)
  
  #List the MODIS files for the year selected above
  MODIS.files.list <- list.files(file.path(MODIS.rawdata.file, year), pattern = "[.]tif$", full.names = TRUE)
  
  #read the data in memory as a list of rasters
  MODIS.yearly.Burn <- lapply(MODIS.files.list, rast) 
  #MODIS.bunt.date <- lapply(MODIS.files.list, function(x){rast(x, subds = 1)}) # subds = 1 means that we want to load the first band of each dataset
  
  #stack all the raster layers together to create a single dataset with the same number of layers as there are elements in the MODIS.files.list
  MODIS.yearly.Burn <- rast(MODIS.yearly.Burn)
  
  #replace the cell values that are not zero to 1 (= the cell has burnt)
  MODIS.yearly.Burn <- subst(MODIS.yearly.Burn, from = 1:366, to = 1)
  #replace the cell values that zero and below zero (= the cell has NOT burnt)
  MODIS.yearly.Burn <- subst(MODIS.yearly.Burn, from = -2:0, to = 0)
  
  #sum all the layers together to create a single raster layer dataset
  MODIS.yearly.Burn <- app(MODIS.yearly.Burn, fun = "sum")
  #plot(MODIS.yearly.Burn)
  
  #reproject to match the projection of the initial Fuel.Age layer
  MODIS.yearly.Burn <- project(MODIS.yearly.Burn, Base.CRS, method = "near")
  #plot(MODIS.yearly.Burn)
  
  #resample to match the resolution of the initial fuel age map
  MODIS.yearly.Burn <- resample(MODIS.yearly.Burn, LU.mask, method = "near")
  #plot(MODIS.yearly.Burn)
  #Reclassify the cell values according to the plan below
  # If, cell value ≥ 1 => 1 (= fuel age for the year is reset to zero)
  # If, cell value = 0 => NA 
  #MODIS.yearly.Burn.test <- ifel(MODIS.yearly.Burn >= 1, 1, NA) #28/03/2023: debugging, this step produced an incorrect map, so use replace function to replace 0 values
  #plot(MODIS.yearly.Burn)
  
  #replace value of 0 with NA(=remove zero values)
  MODIS.yearly.Burn[MODIS.yearly.Burn == 0] <- NA
  
  #update the fuel age layer from the MODIS.bunt.date raster
  FuelAge <- ifel(MODIS.yearly.Burn == 1, 0, (FuelAge + 1))
  #plot(FuelAge)
  
  #mask out water body using LU mask
  #FuelAge.present <- crop(FuelAge.present, LU.mask)
  MODIS.Fuel.Age.accumulation <- mask(FuelAge, LU.mask, updatevalue = NA)
  #plot(MODIS.Fuel.Age.accumulation)
  
  #visualise the results
  terra::plot(MODIS.Fuel.Age.accumulation, colNA = "navy", main = glue("MODIS fuel age with accumulation - {year}"))
  
  #export as tiff for later use
  writeRaster(MODIS.Fuel.Age.accumulation, file.path(MODIS.fuelage.output, glue("MODIS_FuelAge_Accumulation_{year}.tif")), filetype = "GTiff", datatype = 'INT2S', overwrite = TRUE)
  
  #export in rst format for use in MCK
  #writeRaster(MODIS.Fuel.Age.accumulation, file.path(MODIS.fuelage.output, glue("RST_MODIS_FuelAge_Accumulation_{year}.rst")), filetype = "RST", datatype = 'INT2U', overwrite = TRUE)
  
  cat(paste("Done with year ", year))
  cat("\n")
}

Sys.time() - start

#####################################################################
#### Step 5: generate fuel load maps using MODIS fuel age maps#####

#### vegetation map####
#file path to original vegetation type raster
VegeType.file <- file.path("C:/Users/a1677880/UNHaRMED WA SA/UNHaRMEDWesternAustralia/Data/Bushfire/Vegetation")
VegeTypeOriginal.raster <- rast(file.path(VegeType.file, "Fuel_Model_reclassified_100m_GDA94_LU_clip.tif")) #in this raster, there are 10 vegetation types

#Fuel age raster as a template for extent
FuelAge.template <- rast(file.path("C:/Users/a1677880/Earth Observation data pre-processing/Time since last disturbance/MODIS_FuelAge_Accumu_March28/MODIS_FuelAge_Accumulation_2000.tif"))

#crop the fuel type raster to the extent of 4236 x 4349 using the template
VegeTypeOriginal.raster <- crop(VegeTypeOriginal.raster, FuelAge.template)

####Olson curve#####
#create fuel load function based on Olson curve
FuelLoad.function <- function(VegType, FuelAge){ifelse(VegType == 50 & FuelAge == 1, 1, 
                                                       ifelse(VegType == 50 & FuelAge == 2, 1.5, 
                                                              ifelse(VegType == 50 & FuelAge >= 3, 2, 
                                                                     ifelse(VegType == 51, ((1 - exp(-0.252*FuelAge))*13.31), 
                                                                            ifelse(VegType == 52, 4.5, 
                                                                                   ifelse(VegType == 53, 1.5,
                                                                                          ifelse(VegType == 54 & FuelAge >= 0 & FuelAge <= 2, 2, 
                                                                                                 ifelse(VegType == 54 & FuelAge >= 3 & FuelAge <= 4, 10,
                                                                                                        ifelse(VegType == 54 & FuelAge >= 5 & FuelAge <= 9, 15,
                                                                                                               ifelse(VegType == 54 & FuelAge >= 10 & FuelAge <= 20, 20,
                                                                                                                      ifelse(VegType == 54 & FuelAge >= 21, 22,
                                                                                                                             ifelse(VegType == 55 & FuelAge >= 0 & FuelAge <= 2, 5, 
                                                                                                                                    ifelse(VegType == 55 & FuelAge >= 3 & FuelAge <= 4, 12,
                                                                                                                                           ifelse(VegType == 55 & FuelAge >= 5 & FuelAge <= 9, 18,
                                                                                                                                                  ifelse(VegType == 55 & FuelAge >= 10 & FuelAge <= 20, 22.5,
                                                                                                                                                         ifelse(VegType == 55 & FuelAge >= 21, 24.5,
                                                                                                                                                                ifelse(VegType == 56 & FuelAge >= 0 & FuelAge <= 2, 5.5, 
                                                                                                                                                                       ifelse(VegType == 56 & FuelAge >= 3 & FuelAge <= 4, 14.5,
                                                                                                                                                                              ifelse(VegType == 56 & FuelAge >= 5 & FuelAge <= 9, 19.5,
                                                                                                                                                                                     ifelse(VegType == 56 & FuelAge >= 10 & FuelAge <= 20, 24.5,
                                                                                                                                                                                            ifelse(VegType == 56 & FuelAge >= 21, 26, 
                                                                                                                                                                                                   ifelse(VegType == 57 & FuelAge >= 0 & FuelAge <= 2, 5.5, 
                                                                                                                                                                                                          ifelse(VegType == 57 & FuelAge >= 3 & FuelAge <= 4, 15.5,
                                                                                                                                                                                                                 ifelse(VegType == 57 & FuelAge >= 5 & FuelAge <= 9, 19.5,
                                                                                                                                                                                                                        ifelse(VegType == 57 & FuelAge >= 10 & FuelAge <= 20, 24.5,
                                                                                                                                                                                                                               ifelse(VegType == 57 & FuelAge >= 21, 32.5,
                                                                                                                                                                                                                                      ifelse(VegType == 58 & FuelAge == 1, 1,
                                                                                                                                                                                                                                             ifelse(VegType == 58 & FuelAge == 2, 1.6,
                                                                                                                                                                                                                                                    ifelse(VegType == 58 & FuelAge == 3, 2.5, 
                                                                                                                                                                                                                                                           ifelse(VegType == 58 & FuelAge == 4, 3.4, 
                                                                                                                                                                                                                                                                  ifelse(VegType == 58 & FuelAge == 5, 4.2, 
                                                                                                                                                                                                                                                                         ifelse(VegType == 58 & FuelAge == 6, 5, 
                                                                                                                                                                                                                                                                                ifelse(VegType == 58 & FuelAge == 7, 5.8, 
                                                                                                                                                                                                                                                                                       ifelse(VegType == 58 & FuelAge >= 8 & FuelAge <= 9, 6.5, 
                                                                                                                                                                                                                                                                                              ifelse(VegType == 58 & FuelAge >= 10 & FuelAge <= 11, 7.7, 
                                                                                                                                                                                                                                                                                                     ifelse(VegType == 58 & FuelAge >= 12 & FuelAge <= 14, 8.8, 
                                                                                                                                                                                                                                                                                                            ifelse(VegType == 58 & FuelAge >= 15 & FuelAge <= 19, 10.5, 
                                                                                                                                                                                                                                                                                                                   ifelse(VegType == 58 & FuelAge >= 20 & FuelAge <= 24, 12.7, 
                                                                                                                                                                                                                                                                                                                          ifelse(VegType == 58 & FuelAge >= 25, 14.8, 
                                                                                                                                                                                                                                                                                                                                 ifelse(VegType == 59 & FuelAge >= 1 & FuelAge <= 3, 1.1,
                                                                                                                                                                                                                                                                                                                                        ifelse(VegType == 59 & FuelAge >= 4 & FuelAge <= 8, 13,
                                                                                                                                                                                                                                                                                                                                               ifelse(VegType == 59 & FuelAge >= 9 & FuelAge <= 13, 16, 
                                                                                                                                                                                                                                                                                                                                                      ifelse(VegType == 59 & FuelAge >= 14 & FuelAge <= 20, 17, 
                                                                                                                                                                                                                                                                                                                                                             ifelse(VegType == 59 & FuelAge >= 21, 18, 0))))))))))))))))))))))))))))))))))))))))))))
}

# file path to recalculated fuel age 
MODIS.FuelAge <- file.path("C:/Users/a1677880/Earth Observation data pre-processing/Time since last disturbance/MODIS_FuelAge_Accumu")

#read in MODIS fuel age
EOFuelage.list <- list.files(MODIS.FuelAge, pattern = "MODIS_FuelAge_Accumulation.+[.]tif$", full.names = TRUE)
EOFuelAge.raster <- lapply(EOFuelage.list, rast)
EOFuelAge.raster <- rast(EOFuelAge.raster)

# combine the vegetype raster with fuel age raster, and convert to a data frame
EOAge.VegeType.raster <- c(VegeTypeOriginal.raster, EOFuelAge.raster)

#rename the raster layer in an automatic way
year.sequence <- seq(2000, 2021, 1)
file.names <- paste0("FuelAge", year.sequence)
names(EOAge.VegeType.raster) <- c("VegeType", file.names)

#convert to a data frame
EOAge.VegeType.table <- as.data.frame(EOAge.VegeType.raster, xy = T, cells = F) #note we don't need cells as this will make the following process slower

#export the table to save time if re-open R. Use saverds rather than csv to save time
saveRDS(EOAge.VegeType.table, file = file.path(MODIS.fuelload, "MODIS.FuelAge.VegeType.table.rds"))


#read in the table. !!!!Re-run this line before running the for loop 
EOAge.VegeType.table <- readRDS(file.path(MODIS.fuelload, "MODIS.FuelAge.VegeType.table.rds"))

#year <- 2021
start = Sys.time()
for (year in 2000:2021){
  #build an index column number to select the dynamic columns. This is to select from the 2000fuelage column
  indx = (year - 2000) + 4 
  
  #subset the x, y, vegetype, and dynamic fuelage column
  table.subset <- EOAge.VegeType.table[, c(1, 2, 3, indx)]
  
  #apply fuel load function to the selected columns, the glue function in this line is to keep the column name dynamic
  table.subset[, glue("FuelLoad{year}")] <- FuelLoad.function(table.subset[,3], table.subset[,4])  #this line works
  
  #select the x, y, and recalculated fuel load column
  FuelLoad.column <- table.subset[,c(1,2,5)]
  
  #merge the fuel load column to the original table. Not sure why this step?
  #EOAge.VegeType.table <- merge(EOAge.VegeType.table, FuelLoad.column, by = c("x", "y"))
  
  #create a raster layer of the fuel load
  MODIS.FuelLoad <- rast(FuelLoad.column, type = "xyz", crs = crs(EOAge.VegeType.raster))
  
  plot(MODIS.FuelLoad, main = glue("MODIS fuel load {year}"), colNA = "blue")
  #export the 2021 fuel load in RST for use in MCK
  #writeRaster(MODIS.FuelLoad, file.path(MODIS.fuelload, glue("RST_MODIS_FuelLoad_{year}.rst")), filetype = "RST", datatype = 'INT2U', overwrite = TRUE)
  
  
  #export the fuel load rasters for GIS analysis
  writeRaster(MODIS.FuelLoad, file.path(MODIS.fuelload, glue("MODIS_FuelLoad_{year}.tif")), filetype = "GTiff", datatype = 'FLT4S', overwrite = TRUE)
  cat(paste("Done with year ", year))
  cat("\n")
}

Sys.time() - start #Time difference of 4.04738 hours