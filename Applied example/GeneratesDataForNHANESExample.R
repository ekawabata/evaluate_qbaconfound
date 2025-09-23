################################################################################
# GENERATES DATA FOR THE NHANES EXAMPLE

# A. Preliminaries
# loads necessary packages and functions
# defines RemoveSeqn function
# takes sample_old, a list of data frames with a seqn column for individual participant ID
# takes Seqn2Remove, a vector of seqn numbers to remove from each data frame in sample_old
# returns a list of data frames without Seqn2Remove
# defines AlcoholLabel function
# takes data, a data frame with a drxfcld column which contains descriptions for individual food items
# adds a new column Alcoholic to the data
# sets Alcoholic to "alcoholic" if words such as beer, cocktail, gin, margarita, rum, vodka, whiskey, wine or alcoholic (but not nonalcoholic) appear in drxfcld
# sets "nonalcoholic" otherwise
# returns a data frame with Alcoholic column

# B. Imports data from the CDC website (requires an Internet connection)
# imports body measures data
# imports alcohol use data
# imports demographic variables
# imports smoking data
# imports physical activity data
# imports total nutrient intake data
# imports individual foods data
# imports dietary interview technical support
# imports reproductive health data

# C. Cleans data based on Butler et al 2016
# combines all types of data
# removes anyone who is not between 20 and 79 years old
# removes anyone with missing alcohol use data
# removes anyone with missing total nutrient intake data - first day
# removes anyone with missing total nutrient intake data - second day
# removes anyone with missing individual foods data - first day
# removes anyone with missing individual foods data - second day
# removes anyone with missing body measures
# removes anyone who is currently pregnant or breastfeeding
# removes anyone who is currently on a medical or weight-gain/loss diet
# removes anyone with missing demographic variables
# removes anyone with missing physical activity data
# removes anyone with missing smoking data
# removes anyone who is female
# creates complete data for 7,477 men

# D. Defines analysis variables based on Butler et al 2016
# defines continuous outcome: waistcircumferenceY (waist circumference (cm))
# define 5-level categorical exposure: alcoholX (1to2 (baseline), never, former, 3to4, atleast5)
# defines measured confounders:
# 3-level categorical: ageC (20to39, 40to59, 60to79)
# 5-level categorical: raceC (MexicanAmerican, Hispanic, White, Black, Other)
# 4-level categorical: educationC (<HighSchool, HighSchoolGraduate, College, CollegeGraduate)
# 3-level categorical: physicalactivitylevelC (low, moderate, high)
# 5-level categorical: yearC (2003to2004, 2005to2006, 2007to2008, 2009to2010, 2011to2012)
# binary: recall1C (weekend, weekday)
# binary: recall2C (weekend, weekday)
# 3-level categorical: smokingC (never, current, former)
# 3-level categorical: incomeC (0to130, 131to299, atleast300)
# 3-level categorical: misreporingC (accurate, under, over)
# defines unmeasured confounders:
# continuous: energyU (daily non-alcoholic energy intake (kcal))
# binary: maritalstatusU (0 currently in a relationship, 1 currently single)
# writes the data to a csv file
################################################################################

################################################################################
# A. Preliminaries
################################################################################

# REMOVE ALL OBJECTS (IF ANY)
rm(list=ls())
# CLOSE OPEN GRAPHICS DEVICES (IF ANY)
graphics.off()

# INSTALL Hmisc IF NECESSARY
if (!("Hmisc" %in% rownames(installed.packages())))
{
  install.packages("Hmisc")
}
# LOAD Hmisc
library(Hmisc)

# DEFINE RemoveSeqn FUNCTION
# take sample_old, a list of data frames such as nhanes, paqiaf, paq2003, paq2007, dr1iff and dr2iff, each with a seqn column for individual participant ID
# take Seqn2Remove, a vector of seqn numbers to remove from each data frame in sample_old
RemoveSeqn <- function(sample_old, Seqn2Remove)
{
  sample_new <- sample_old
  for (listitem in 1:length(sample_old))
  {
    sample_new[[listitem]] <- sample_old[[listitem]][!is.element(sample_old[[listitem]]$seqn, Seqn2Remove),]
  }
  # return a list of data frames without Seqn2Remove
  return(sample_new)
}

# DEFINE AlcoholLabel FUNCTION
# take data, which contains a drxfcld column which contains descriptions for individual food items
AlcoholLabel <- function(data)
{
  # change all letters in drxfcld to lowercase
  data$drxfcld <- tolower(data$drxfcld)
  # adds a new column Alcoholic to the data
  data <- data.frame(data, "Alcoholic"="nonalcoholic")
  # sets Alcoholic to "alcoholic" if words such as beer, cocktail, gin, margarita, rum, vodka, whiskey, wine or alcoholic (but not nonalcoholic) appear in drxfcld
  # sets "nonalcoholic" otherwise
  data[grep("beer", data$drxfcld), "Alcoholic"] <- "alcoholic"
  data[grep("cocktail", data$drxfcld), "Alcoholic"] <- "alcoholic"
  data[grep("gin", data$drxfcld), "Alcoholic"] <- "alcoholic"
  data[grep("ginger", data$drxfcld), "Alcoholic"] <- "nonalcoholic"
  data[grep("original", data$drxfcld), "Alcoholic"] <- "nonalcoholic"
  data[grep("margarita", data$drxfcld), "Alcoholic"] <- "alcoholic"
  data[grep("rum", data$drxfcld), "Alcoholic"] <- "alcoholic"
  data[grep("crum", data$drxfcld), "Alcoholic"] <- "nonalcoholic"
  data[grep("drumstick", data$drxfcld), "Alcoholic"] <- "nonalcoholic"
  data[grep("vodka", data$drxfcld), "Alcoholic"] <- "alcoholic"
  data[grep("whiskey", data$drxfcld), "Alcoholic"] <- "alcoholic"
  data[grep("wine", data$drxfcld), "Alcoholic"] <- "alcoholic"
  data[grep("alcoholic", data$drxfcld), "Alcoholic"] <- "alcoholic"
  data[grep("nonalcoholic", data$drxfcld), "Alcoholic"] <- "nonalcoholic"
  # return a data frame with Alcoholic column
  return(data)
}

################################################################################
# B. Import data from the CDC website
################################################################################

# IMPORT BODY MEASURES DATA
BmxVariables <- c("seqn", "bmxwaist", "bmxht", "bmxwt")

# import body measures from 2003-2004 cycle
bmx_c_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/BMX_C.xpt")
# remove anyone whose bmxwaist measure is outside the possible range of values
bmx_c_temp2 <- bmx_c_temp1[bmx_c_temp1[,"bmxwaist"]>=32 & bmx_c_temp1[,"bmxwaist"]<=170.7 | is.na(bmx_c_temp1[,"bmxwaist"]),]
# remove anyone whose bmxht measure is outside the possible range of values
bmx_c_temp3 <- bmx_c_temp2[bmx_c_temp2[,"bmxht"]>=79.7 & bmx_c_temp2[,"bmxht"]<=204.4 | is.na(bmx_c_temp2[,"bmxht"]),]
# remove anyone whose bmxwt measure is outside the possible range of values
bmx_c_temp4 <- bmx_c_temp3[bmx_c_temp3[,"bmxwt"]>=2.4 & bmx_c_temp3[,"bmxwt"]<=209.1 | is.na(bmx_c_temp3[,"bmxwt"]),]

# import body measures from 2005-2006 cycle
bmx_d_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/BMX_D.xpt")
# remove anyone whose bmxwaist measure is outside the possible range of values
bmx_d_temp2 <- bmx_d_temp1[bmx_d_temp1[,"bmxwaist"]>=40.2 & bmx_d_temp1[,"bmxwaist"]<=175 | is.na(bmx_d_temp1[,"bmxwaist"]),]
# remove anyone whose bmxht measure is outside the possible range of values
bmx_d_temp3 <- bmx_d_temp2[bmx_d_temp2[,"bmxht"]>=79.8 & bmx_d_temp2[,"bmxht"]<=204.1 | is.na(bmx_d_temp2[,"bmxht"]),]
# remove anyone whose bmxwt measure is outside the possible range of values
bmx_d_temp4 <- bmx_d_temp3[bmx_d_temp3[,"bmxwt"]>=2.6 & bmx_d_temp3[,"bmxwt"]<=371 | is.na(bmx_d_temp3[,"bmxwt"]),]

# import body measures from 2007-2008 cycle
bmx_e_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2007/DataFiles/BMX_E.xpt")
# remove anyone whose bmxwaist measure is outside the possible range of values
bmx_e_temp2 <- bmx_e_temp1[bmx_e_temp1[,"bmxwaist"]>=37.8 & bmx_e_temp1[,"bmxwaist"]<=178.2 | is.na(bmx_e_temp1[,"bmxwaist"]),]
# remove anyone whose bmxht measure is outside the possible range of values
bmx_e_temp3 <- bmx_e_temp2[bmx_e_temp2[,"bmxht"]>=81.5 & bmx_e_temp2[,"bmxht"]<=203.8 | is.na(bmx_e_temp2[,"bmxht"]),]
# remove anyone whose bmxwt measure is outside the possible range of values
bmx_e_temp4 <- bmx_e_temp3[bmx_e_temp3[,"bmxwt"]>=3.1 & bmx_e_temp3[,"bmxwt"]<=218.2 | is.na(bmx_e_temp3[,"bmxwt"]),]

# import body measures from 2009-2010 cycle
bmx_f_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/BMX_F.xpt")
# remove anyone whose bmxwaist measure is outside the possible range of values
bmx_f_temp2 <- bmx_f_temp1[bmx_f_temp1[,"bmxwaist"]>=40.7 & bmx_f_temp1[,"bmxwaist"]<=179 | is.na(bmx_f_temp1[,"bmxwaist"]),]
# remove anyone whose bmxht measure is outside the possible range of values
bmx_f_temp3 <- bmx_f_temp2[bmx_f_temp2[,"bmxht"]>=79.1 & bmx_f_temp2[,"bmxht"]<=202.7 | is.na(bmx_f_temp2[,"bmxht"]),]
# remove anyone whose bmxwt measure is outside the possible range of values
bmx_f_temp4 <- bmx_f_temp3[bmx_f_temp3[,"bmxwt"]>=2.7 & bmx_f_temp3[,"bmxwt"]<=239.4 | is.na(bmx_f_temp3[,"bmxwt"]),]

# import body measures from 2011-2012 cycle
bmx_g_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/BMX_G.xpt")
# remove anyone whose bmxwaist measure is outside the possible range of values
bmx_g_temp2 <- bmx_g_temp1[bmx_g_temp1[,"bmxwaist"]>=0 & bmx_g_temp1[,"bmxwaist"]<=200 | is.na(bmx_g_temp1[,"bmxwaist"]),]
# remove anyone whose bmxht measure is outside the possible range of values
bmx_g_temp3 <- bmx_g_temp2[bmx_g_temp2[,"bmxht"]>=70.0 & bmx_g_temp2[,"bmxht"]<=208.4 | is.na(bmx_g_temp2[,"bmxht"]),]
# remove anyone whose bmxwt measure is outside the possible range of values
bmx_g_temp4 <- bmx_g_temp3[bmx_g_temp3[,"bmxwt"]>=0 & bmx_g_temp3[,"bmxwt"]<=440 | is.na(bmx_g_temp3[,"bmxwt"]),]

# combine body measures from all survey cycles
bmx_temp1 <- rbind(bmx_c_temp4[,BmxVariables], bmx_d_temp4[,BmxVariables], bmx_e_temp4[,BmxVariables], bmx_f_temp4[,BmxVariables], bmx_g_temp4[,BmxVariables])
# set bmx_temp1 to bmx which may still include missing (NA) measurements
bmx <- bmx_temp1

# IMPORT ALCOHOL USE DATA
AlqVariables <- c("seqn", "alq101", "alq110", "alq120q", "alq130")

# import alcohol use data from 2003-2004 cycle
alq_c_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/ALQ_C.xpt")
alq_c_temp2 <- alq_c_temp1
# set alq120q measure to NA if recorded 777 or 999
alq_c_temp2[which(alq_c_temp2[,"alq120q"]==777 | alq_c_temp2[,"alq120q"]==999),"alq120q"] <- NA
# remove anyone whose alq120q measure is outside the possible range of values
alq_c_temp3 <- alq_c_temp2[alq_c_temp2[,"alq120q"]>=0 & alq_c_temp2[,"alq120q"]<=365 | is.na(alq_c_temp2[,"alq120q"]),]
alq_c_temp4 <- alq_c_temp3
# set alq130 measure to NA if recorded 77 or 99
alq_c_temp4[which(alq_c_temp4[,"alq130"]==77 | alq_c_temp4[,"alq130"]==99),"alq130"] <- NA
# remove anyone whose alq130 measure is outside the possible range of values
alq_c_temp5 <- alq_c_temp4[alq_c_temp4[,"alq130"]>=1 & alq_c_temp4[,"alq130"]<=50 | is.na(alq_c_temp4[,"alq130"]),]

# import alcohol use data from 2005-2006 cycle
alq_d_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/ALQ_D.xpt")
alq_d_temp2 <- alq_d_temp1
# set alq120q measure to NA if recorded 777 or 999
alq_d_temp2[which(alq_d_temp2[,"alq120q"]==777 | alq_d_temp2[,"alq120q"]==999),"alq120q"] <- NA
# remove anyone whose alq120q measure is outside the possible range of values
alq_d_temp3 <- alq_d_temp2[alq_d_temp2[,"alq120q"]>=0 & alq_d_temp2[,"alq120q"]<=365 | is.na(alq_d_temp2[,"alq120q"]),]
alq_d_temp4 <- alq_d_temp3
# set alq130 measure to NA if recorded 777 or 999
alq_d_temp4[which(alq_d_temp4[,"alq130"]==777 | alq_d_temp4[,"alq130"]==999),"alq130"] <- NA
# remove anyone whose alq130 measure is outside the possible range of values
alq_d_temp5 <- alq_d_temp4[alq_d_temp4[,"alq130"]>=1 & alq_d_temp4[,"alq130"]<=32 | is.na(alq_d_temp4[,"alq130"]),]

# import alcohol use data from 2007-2008 cycle
alq_e_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2007/DataFiles/ALQ_E.xpt")
alq_e_temp2 <- alq_e_temp1
# set alq120q measure to NA if recorded 777 or 999
alq_e_temp2[which(alq_e_temp2[,"alq120q"]==777 | alq_e_temp2[,"alq120q"]==999),"alq120q"] <- NA
# remove anyone whose alq120q measure is outside the possible range of values
alq_e_temp3 <- alq_e_temp2[alq_e_temp2[,"alq120q"]>=0 & alq_e_temp2[,"alq120q"]<=365 | is.na(alq_e_temp2[,"alq120q"]),]
alq_e_temp4 <- alq_e_temp3
# set alq130 measure to NA if recorded 777 or 999
alq_e_temp4[which(alq_e_temp4[,"alq130"]==777 | alq_e_temp4[,"alq130"]==999),"alq130"] <- NA
# remove anyone whose alq130 measure is outside the possible range of values
alq_e_temp5 <- alq_e_temp4[alq_e_temp4[,"alq130"]>=1 & alq_e_temp4[,"alq130"]<=83 | is.na(alq_e_temp4[,"alq130"]),]

# import alcohol use data from 2009-2010 cycle
alq_f_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/ALQ_F.xpt")
alq_f_temp2 <- alq_f_temp1
# set alq120q measure to NA if recorded 777 or 999
alq_f_temp2[which(alq_f_temp2[,"alq120q"]==777 | alq_f_temp2[,"alq120q"]==999),"alq120q"] <- NA
# remove anyone whose alq120q measure is outside the possible range of values
alq_f_temp3 <- alq_f_temp2[alq_f_temp2[,"alq120q"]>=0 & alq_f_temp2[,"alq120q"]<=360 | is.na(alq_f_temp2[,"alq120q"]),]
alq_f_temp4 <- alq_f_temp3
# set alq130 measure to NA if recorded 777 or 999
alq_f_temp4[which(alq_f_temp4[,"alq130"]==777 | alq_f_temp4[,"alq130"]==999),"alq130"] <- NA
# remove anyone whose alq130 measure is outside the possible range of values
alq_f_temp5 <- alq_f_temp4[alq_f_temp4[,"alq130"]>=1 & alq_f_temp4[,"alq130"]<=36 | is.na(alq_f_temp4[,"alq130"]),]

# import alcohol use data from 2011-2012 cycle
alq_g_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/ALQ_G.xpt")
alq_g_temp2 <- alq_g_temp1
# set alq120q measure to NA if recorded 777 or 999
alq_g_temp2[which(alq_g_temp2[,"alq120q"]==777 | alq_g_temp2[,"alq120q"]==999),"alq120q"] <- NA
# remove anyone whose alq120q measure is outside the possible range of values
alq_g_temp3 <- alq_g_temp2[alq_g_temp2[,"alq120q"]>=0 & alq_g_temp2[,"alq120q"]<=350 | is.na(alq_g_temp2[,"alq120q"]),]
alq_g_temp4 <- alq_g_temp3
# set alq130 measure to NA if recorded 777 or 999
alq_g_temp4[which(alq_g_temp4[,"alq130"]==777 | alq_g_temp4[,"alq130"]==999),"alq130"] <- NA
# remove anyone whose alq130 measure is outside the possible range of values
alq_g_temp5 <- alq_g_temp4[alq_g_temp4[,"alq130"]>=1 & alq_g_temp4[,"alq130"]<=82 | is.na(alq_g_temp4[,"alq130"]),]

# combine alcohol use data from all survey cycles
alq_temp1 <- rbind(alq_c_temp5[,AlqVariables], alq_d_temp5[,AlqVariables], alq_e_temp5[,AlqVariables], alq_f_temp5[,AlqVariables], alq_g_temp5[,AlqVariables])
alq_temp2 <- alq_temp1
# set alq101 measure to NA if recorded 7 or 9
alq_temp2[which(alq_temp2[,"alq101"]==7 | alq_temp2[,"alq101"]==9),"alq101"] <- NA
# remove anyone whose alq101 measure is outside the possible range of values
alq_temp3 <- alq_temp2[alq_temp2[,"alq101"]==1 | alq_temp2[,"alq101"]==2 | is.na(alq_temp2[,"alq101"]),]
alq_temp4 <- alq_temp3
# set alq110 measure to NA if recorded 7 or 9
alq_temp4[which(alq_temp4[,"alq110"]==7 | alq_temp4[,"alq110"]==9),"alq110"] <- NA
# remove anyone whose alq110 measure is outside the possible range of values
alq_temp5 <- alq_temp4[alq_temp4[,"alq110"]==1 | alq_temp4[,"alq110"]==2 | is.na(alq_temp4[,"alq110"]),]
# set alq_temp5 to alq which may still include missing (NA) measurements
alq <- alq_temp5

# IMPORT DEMOGRAPHIC VARIABLES
DemoVariables <- c("seqn", "riagendr", "ridageyr", "ridreth1", "dmdeduc2", "dmdmartl", "indfmpir", "cycle")

# import demographic variables from 2003-2004 cycle
demo_c_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/DEMO_C.xpt")
# remove anyone with ridageyr measure outside the possible range of values
demo_c_temp2 <- demo_c_temp1[demo_c_temp1["ridageyr"]>=0 & demo_c_temp1["ridageyr"]<=85 | is.na(demo_c_temp1["ridageyr"]),]

# import demographic variables from 2005-2006 cycle
demo_d_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/DEMO_D.xpt")
demo_d_temp2 <- demo_d_temp1[demo_d_temp1["ridageyr"]>=0 & demo_d_temp1["ridageyr"]<=85 | is.na(demo_d_temp1["ridageyr"]),]

# import demographic variables from 2007-2008 cycle
demo_e_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2007/DataFiles/DEMO_E.xpt")
demo_e_temp2 <- demo_e_temp1[demo_e_temp1["ridageyr"]>=0 & demo_e_temp1["ridageyr"]<=80 | is.na(demo_e_temp1["ridageyr"]),]

# import demographic variables from 2009-2010 cycle
demo_f_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/DEMO_F.xpt")
demo_f_temp2 <- demo_f_temp1[demo_f_temp1["ridageyr"]>=0 & demo_f_temp1["ridageyr"]<=80 | is.na(demo_f_temp1["ridageyr"]),]

# import demographic variables from 2011-2012 cycle
demo_g_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/DEMO_G.xpt")
demo_g_temp2 <- demo_g_temp1[demo_g_temp1["ridageyr"]>=0 & demo_g_temp1["ridageyr"]<=80 | is.na(demo_g_temp1["ridageyr"]),]

# add a new column called cycle for Nhanes survey cycle
demo_c_temp3 <- data.frame(demo_c_temp2, "cycle"="2003to2004")
demo_d_temp3 <- data.frame(demo_d_temp2, "cycle"="2005to2006")
demo_e_temp3 <- data.frame(demo_e_temp2, "cycle"="2007to2008")
demo_f_temp3 <- data.frame(demo_f_temp2, "cycle"="2009to2010")
demo_g_temp3 <- data.frame(demo_g_temp2, "cycle"="2011to2012")

# combine demographic variables from all survey cycles
demo_temp1 <- rbind(demo_c_temp3[,DemoVariables], demo_d_temp3[,DemoVariables], demo_e_temp3[,DemoVariables], demo_f_temp3[,DemoVariables], demo_g_temp3[,DemoVariables])
# remove anyone with ridreth1 measure is outside the possible range of values
demo_temp2 <- demo_temp1[demo_temp1["ridreth1"]==1 | demo_temp1["ridreth1"]==2 | demo_temp1["ridreth1"]==3 | demo_temp1["ridreth1"]==4 | demo_temp1["ridreth1"]==5 | is.na(demo_temp1["ridreth1"]),]
demo_temp3 <- demo_temp2
# set dmdeduc2 measure to NA if recorded 7 or 9
demo_temp3[which(demo_temp3[,"dmdeduc2"]==7 | demo_temp3[,"dmdeduc2"]==9),"dmdeduc2"] <- NA
# remove anyone with dmdeduc2 measure is outside the possible range of values
demo_temp4 <- demo_temp3[demo_temp3["dmdeduc2"]==1 | demo_temp3["dmdeduc2"]==2 | demo_temp3["dmdeduc2"]==3 | demo_temp3["dmdeduc2"]==4 | demo_temp3["dmdeduc2"]==5 | is.na(demo_temp3["dmdeduc2"]),]
demo_temp5 <- demo_temp4
# set dmdmartl measure to NA if recorded 77 or 99
demo_temp5[which(demo_temp5[,"dmdmartl"]==77 | demo_temp5[,"dmdmartl"]==99),"dmdmartl"] <- NA
# remove anyone with dmdmartl measure is outside the possible range of values
demo_temp6 <- demo_temp5[demo_temp5[,"dmdmartl"]==1 | demo_temp5[,"dmdmartl"]==2 | demo_temp5[,"dmdmartl"]==3 | demo_temp5[,"dmdmartl"]==4 | demo_temp5[,"dmdmartl"]==5 | demo_temp5[,"dmdmartl"]==6 | is.na(demo_temp5[,"dmdmartl"]),]
# remove anyone with indfmpir measure is outside the possible range of values
demo_temp7 <- demo_temp6[demo_temp6[,"indfmpir"]>=0 & demo_temp6[,"indfmpir"]<=5 | is.na(demo_temp6[,"indfmpir"]),]
# set demo_temp7 to demo which may still include missing (NA) measurements
demo <- demo_temp7

# IMPORT SMOKING DATA
SmqVariables <- c("seqn", "smq020", "smq040", "smq050q")

# import smoking data from 2003-2004 cycle
smq_c_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/SMQ_C.xpt")
smq_c_temp2 <- smq_c_temp1
# set smq050q measure to NA if recorded 77777 or 99999
smq_c_temp2[which(smq_c_temp2[,"smq050q"]==77777 | smq_c_temp2[,"smq050q"]==99999),"smq050q"] <- NA
# remove anyone whose smq050q measure is outside the possible range of values
smq_c_temp3 <- smq_c_temp2[smq_c_temp2[,"smq050q"]>=1 & smq_c_temp2[,"smq050q"]<=75 | is.na(smq_c_temp2[,"smq050q"]),]

# import smoking data from 2005-2006 cycle
smq_d_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/SMQ_D.xpt")
smq_d_temp2 <- smq_d_temp1
# set smq050q measure to NA if recorded 77777 or 99999
smq_d_temp2[which(smq_d_temp2[,"smq050q"]==77777 | smq_d_temp2[,"smq050q"]==99999),"smq050q"] <- NA
# remove anyone whose smq050q measure is outside the possible range of values
smq_d_temp3 <- smq_d_temp2[smq_d_temp2[,"smq050q"]>=1 & smq_d_temp2[,"smq050q"]<=72 | is.na(smq_d_temp2[,"smq050q"]),]

# import smoking data from 2007-2008 cycle
smq_e_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2007/DataFiles/SMQ_E.xpt")
smq_e_temp2 <- smq_e_temp1
# set smq050q measure to NA if recorded 77777 or 99999
smq_e_temp2[which(smq_e_temp2[,"smq050q"]==77777 | smq_e_temp2[,"smq050q"]==99999),"smq050q"] <- NA
# remove anyone whose smq050q measure is outside the possible range of values
smq_e_temp3 <- smq_e_temp2[smq_e_temp2[,"smq050q"]>=1 & smq_e_temp2[,"smq050q"]<=66 | is.na(smq_e_temp2[,"smq050q"]),]

# import smoking data from 2009-2010 cycle
smq_f_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/SMQ_F.xpt")
smq_f_temp2 <- smq_f_temp1
# set smq050q measure to NA if recorded 77777 or 99999
smq_f_temp2[which(smq_f_temp2[,"smq050q"]==77777 | smq_f_temp2[,"smq050q"]==99999),"smq050q"] <- NA
# remove anyone whose smq050q measure is outside the possible range of values
smq_f_temp3 <- smq_f_temp2[smq_f_temp2[,"smq050q"]>=1 & smq_f_temp2[,"smq050q"]<=67 | is.na(smq_f_temp2[,"smq050q"]),]

# import smoking data from 2011-2012 cycle
smq_g_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/SMQ_G.xpt")
smq_g_temp2 <- smq_g_temp1
# set smq050q measure to NA if recorded 77777 or 99999
smq_g_temp2[which(smq_g_temp2[,"smq050q"]==77777 | smq_g_temp2[,"smq050q"]==99999),"smq050q"] <- NA
# remove anyone whose smq050q measure is outside the possible range of values
smq_g_temp3 <- smq_g_temp2[smq_g_temp2[,"smq050q"]>=1 & smq_g_temp2[,"smq050q"]<=193 | is.na(smq_g_temp2[,"smq050q"]),]

# combine smoking data from all survey cycles
smq_temp1 <- rbind(smq_c_temp3[,SmqVariables], smq_d_temp3[,SmqVariables], smq_e_temp3[,SmqVariables], smq_f_temp3[,SmqVariables], smq_g_temp3[,SmqVariables])
smq_temp2 <- smq_temp1
# set smq020 measure to NA if recorded 7 or 9
smq_temp2[which(smq_temp2[,"smq020"]==7 | smq_temp2[,"smq020"]==9),"smq020"] <- NA
# remove anyone with smq020 measure outside the possible range of values
smq_temp3 <- smq_temp2[smq_temp2[,"smq020"]==1 | smq_temp2[,"smq020"]==2 | is.na(smq_temp2[,"smq020"]),]
smq_temp4 <- smq_temp3
# set smq040 measure to NA if recorded 7 or 9
smq_temp4[which(smq_temp4[,"smq040"]==7 | smq_temp4[,"smq040"]==9),"smq040"] <- NA
# remove anyone with smq040 measure outside the possible range of values
smq_temp5 <- smq_temp4[smq_temp4[,"smq040"]==1 | smq_temp4[,"smq040"]==2 | smq_temp4[,"smq040"]==3 | is.na(smq_temp4[,"smq040"]),]
# set smq_temp5 to smq which may still include missing (NA) measurements
smq <- smq_temp5

# IMPORT PHYSICAL ACTIVITY DATA (paqiaf) FROM 2003-2004 AND 2005-2006 CYCLES
PaqiafVariables <- c("seqn", "padtimes", "paddurat", "padmets")

# import physical activity data (paqiaf) from 2003-2004 and 2005-2006 cycles
paqiaf_c_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/PAQIAF_C.xpt")
paqiaf_d_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/PAQIAF_D.xpt")

# combine physical activity data (paqiaf) from 2003-2004 and 2005-2006 cycles
paqiaf_temp1 <- rbind(paqiaf_c_temp1[,PaqiafVariables], paqiaf_d_temp1[,PaqiafVariables])
# remove anyone with padtimes measure is outside the possible range of values
paqiaf_temp2 <- paqiaf_temp1[paqiaf_temp1[,"padtimes"]>=0 & paqiaf_temp1[,"padtimes"]<=300 | is.na(paqiaf_temp1[,"padtimes"]),]
# remove anyone with paddurat measure is outside the possible range of values
paqiaf_temp3 <- paqiaf_temp2[paqiaf_temp2[,"paddurat"]>=10 & paqiaf_temp2[,"paddurat"]<=600 | is.na(paqiaf_temp2[,"paddurat"]),]
# remove anyone with padmets measure is outside the possible range of values
paqiaf_temp4 <- paqiaf_temp3[paqiaf_temp3[,"padmets"]>=2.5 & paqiaf_temp3[,"padmets"]<=10 | is.na(paqiaf_temp3[,"padmets"]),]
# set paqiaf_temp4 to paqiaf which may still include missing (NA) measurements
paqiaf <- paqiaf_temp4

# IMPORT PHYSICAL ACTIVITY DATA (paq) FROM 2003-2004 AND 2005-2006 CYCLES
Paq2003Variables <- c("seqn", "pad020", "paq050q", "paq050u", "pad080", "paq100", "pad120", "pad160")

# import physical activity data (paq) from 2003-2004 cycle
paq_c_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/PAQ_C.xpt")
paq_c_temp2 <- paq_c_temp1
# set paq050q measure to NA if recorded 77777 or 99999
paq_c_temp2[which(paq_c_temp2[,"paq050q"]==77777 | paq_c_temp2[,"paq050q"]==99999),"paq050q"] <- NA
# remove anyone whose paq050q measure is outside the possible range of values
paq_c_temp3 <- paq_c_temp2[paq_c_temp2[,"paq050q"]>=1 & paq_c_temp2[,"paq050q"]<=124 | is.na(paq_c_temp2[,"paq050q"]),]
paq_c_temp4 <- paq_c_temp3
# set pad120 measure to NA if recorded 77777 or 99999
paq_c_temp4[which(paq_c_temp4[,"pad120"]==77777 | paq_c_temp4[,"pad120"]==99999),"pad120"] <- NA
# remove anyone whose pad120 measure is outside the possible range of values
paq_c_temp5 <- paq_c_temp4[paq_c_temp4[,"pad120"]>=1 & paq_c_temp4[,"pad120"]<=300 | is.na(paq_c_temp4[,"pad120"]),]
paq_c_temp6 <- paq_c_temp5
# set pad160 measure to NA if recorded 77777 or 99999
paq_c_temp6[which(paq_c_temp6[,"pad160"]==77777 | paq_c_temp6[,"pad160"]==99999),"pad160"] <- NA
# remove anyone whose pad160 measure is outside the possible range of values
paq_c_temp7 <- paq_c_temp6[paq_c_temp6[,"pad160"]>=1 & paq_c_temp6[,"pad160"]<=600 | is.na(paq_c_temp6[,"pad160"]),]

# import physical activity data (paq) from 2005-2006 cycle
paq_d_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/PAQ_D.xpt")
paq_d_temp2 <- paq_d_temp1
# set paq050q measure to NA if recorded 77777 or 99999
paq_d_temp2[which(paq_d_temp2[,"paq050q"]==77777 | paq_d_temp2[,"paq050q"]==99999),"paq050q"] <- NA
# remove anyone whose paq050q measure is outside the possible range of values
paq_d_temp3 <- paq_d_temp2[paq_d_temp2[,"paq050q"]>=1 & paq_d_temp2[,"paq050q"]<=100 | is.na(paq_d_temp2[,"paq050q"]),]
paq_d_temp4 <- paq_d_temp3
# set pad120 measure to NA if recorded 77777 or 99999
paq_d_temp4[which(paq_d_temp4[,"pad120"]==77777 | paq_d_temp4[,"pad120"]==99999),"pad120"] <- NA
# remove anyone whose pad120 measure is outside the possible range of values
paq_d_temp5 <- paq_d_temp4[paq_d_temp4[,"pad120"]>=1 & paq_d_temp4[,"pad120"]<=120 | is.na(paq_d_temp4[,"pad120"]),]
paq_d_temp6 <- paq_d_temp5
# set pad160 measure to NA if recorded 77777 or 99999
paq_d_temp6[which(paq_d_temp6[,"pad160"]==77777 | paq_d_temp6[,"pad160"]==99999),"pad160"] <- NA
# remove anyone whose pad160 measure is outside the possible range of values
paq_d_temp7 <- paq_d_temp6[paq_d_temp6[,"pad160"]>=10 & paq_d_temp6[,"pad160"]<=600 | is.na(paq_d_temp6[,"pad160"]),]

# combine physical activity data (paq) from 2003-2004 and 2005-2006 cycles
paq2003_temp1 <- rbind(paq_c_temp7[,Paq2003Variables], paq_d_temp7[,Paq2003Variables])
paq2003_temp2 <- paq2003_temp1
# set pad020 measure to NA if recorded 7 or 9
paq2003_temp2[which(paq2003_temp2[,"pad020"]==7 | paq2003_temp2[,"pad020"]==9),"pad020"] <- NA
# remove anyone whose pad020 measure is outside the possible range of values
paq2003_temp3 <- paq2003_temp2[paq2003_temp2[,"pad020"]==1 | paq2003_temp2[,"pad020"]==2 | paq2003_temp2[,"pad020"]==3 | is.na(paq2003_temp2[,"pad020"]),]
paq2003_temp4 <- paq2003_temp3
# set paq050u measure to NA if recorded 7 or 9
paq2003_temp4[which(paq2003_temp4[,"paq050u"]==7 | paq2003_temp4[,"paq050u"]==9),"paq050u"] <- NA
# remove anyone whose paq050u measure is outside the possible range of values
paq2003_temp5 <- paq2003_temp4[paq2003_temp4[,"paq050u"]==1 | paq2003_temp4[,"paq050u"]==2 | paq2003_temp4[,"paq050u"]==3 | is.na(paq2003_temp4[,"paq050u"]),]
paq2003_temp6 <- paq2003_temp5
# set pad080 measure to NA if recorded 77777 or 99999
paq2003_temp6[which(paq2003_temp6[,"pad080"]==77777 | paq2003_temp6[,"pad080"]==99999),"pad080"] <- NA
# remove anyone whose pad080 measure is outside the possible range of values
paq2003_temp7 <- paq2003_temp6[paq2003_temp6[,"pad080"]>=1 & paq2003_temp6[,"pad080"]<=600 | is.na(paq2003_temp6[,"pad080"]),]
paq2003_temp8 <- paq2003_temp7
# set paq100 measure to NA if recorded 7 or 9
paq2003_temp8[which(paq2003_temp8[,"paq100"]==7 | paq2003_temp8[,"paq100"]==9),"paq100"] <- NA
# remove anyone whose paq100 measure is outside the possible range of values
paq2003_temp9 <- paq2003_temp8[paq2003_temp8[,"paq100"]==1 | paq2003_temp8[,"paq100"]==2 | paq2003_temp8[,"paq100"]==3 | is.na(paq2003_temp8[,"paq100"]),]
# set paq2003_temp9 to paq2003 which may still include missing (NA) measurements
paq2003 <- paq2003_temp9

# IMPORT PHYSICAL ACTIVITY DATA (paq) FROM 2007-2008, 2009-2010 AND 2011-2012 CYCLES
Paq2007Variables <- c("seqn", "paq605", "paq610", "pad615", "paq620", "paq625", "pad630", "paq635", "paq640", "pad645", "paq650", "paq655", "pad660", "paq665", "paq670", "pad675")

# import physical activity data (paq) from 2007-2008 cycle
paq_e_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2007/DataFiles/PAQ_E.xpt")
paq_e_temp2 <- paq_e_temp1
# set pad615 measure to NA if recorded 7777 or 9999
paq_e_temp2[which(paq_e_temp2[,"pad615"]==7777 | paq_e_temp2[,"pad615"]==9999),"pad615"] <- NA
# remove anyone whose pad615 measure is outside the possible range of values
paq_e_temp3 <- paq_e_temp2[paq_e_temp2[,"pad615"]>=10 & paq_e_temp2[,"pad615"]<=960 | is.na(paq_e_temp2[,"pad615"]),]
paq_e_temp4 <- paq_e_temp3
# set pad630 measure to NA if recorded 7777 or 9999
paq_e_temp4[which(paq_e_temp4[,"pad630"]==7777 | paq_e_temp4[,"pad630"]==9999),"pad630"] <- NA
# remove anyone whose pad630 measure is outside the possible range of values
paq_e_temp5 <- paq_e_temp4[paq_e_temp4[,"pad630"]>=10 & paq_e_temp4[,"pad630"]<=1440 | is.na(paq_e_temp4[,"pad630"]),]
paq_e_temp6 <- paq_e_temp5
# set pad645 measure to NA if recorded 7777 or 9999
paq_e_temp6[which(paq_e_temp6[,"pad645"]==7777 | paq_e_temp6[,"pad645"]==9999),"pad645"] <- NA
# remove anyone whose pad645 measure is outside the possible range of values
paq_e_temp7 <- paq_e_temp6[paq_e_temp6[,"pad645"]>=10 & paq_e_temp6[,"pad645"]<=600 | is.na(paq_e_temp6[,"pad645"]),]
paq_e_temp8 <- paq_e_temp7
# set paq655 measure to NA if recorded 77 or 99
paq_e_temp8[which(paq_e_temp8[,"paq655"]==77 | paq_e_temp8[,"paq655"]==99),"paq655"] <- NA
# remove anyone whose paq655 measure is outside the possible range of values
paq_e_temp9 <- paq_e_temp8[paq_e_temp8[,"paq655"]>=1 & paq_e_temp8[,"paq655"]<=7 | is.na(paq_e_temp8[,"paq655"]),]
paq_e_temp10 <- paq_e_temp9
# set pad660 measure to NA if recorded 7777 or 9999
paq_e_temp10[which(paq_e_temp10[,"pad660"]==7777 | paq_e_temp10[,"pad660"]==9999),"pad660"] <- NA
# remove anyone whose pad660 measure is outside the possible range of values
paq_e_temp11 <- paq_e_temp10[paq_e_temp10[,"pad660"]>=10 & paq_e_temp10[,"pad660"]<=990 | is.na(paq_e_temp10[,"pad660"]),]
paq_e_temp12 <- paq_e_temp11
# set pad675 measure to NA if recorded 7777 or 9999
paq_e_temp12[which(paq_e_temp12[,"pad675"]==7777 | paq_e_temp12[,"pad675"]==9999),"pad675"] <- NA
# remove anyone whose pad675 measure is outside the possible range of values
paq_e_temp13 <- paq_e_temp12[paq_e_temp12[,"pad675"]>=10 & paq_e_temp12[,"pad675"]<=720 | is.na(paq_e_temp12[,"pad675"]),]

# import physical activity data (paq) from 2009-2010 cycle
paq_f_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/PAQ_F.xpt")
paq_f_temp2 <- paq_f_temp1
# set pad615 measure to NA if recorded 7777 or 9999
paq_f_temp2[which(paq_f_temp2[,"pad615"]==7777 | paq_f_temp2[,"pad615"]==9999),"pad615"] <- NA
# remove anyone whose pad615 measure is outside the possible range of values
paq_f_temp3 <- paq_f_temp2[paq_f_temp2[,"pad615"]>=10 & paq_f_temp2[,"pad615"]<=960 | is.na(paq_f_temp2[,"pad615"]),]
paq_f_temp4 <- paq_f_temp3
# set pad630 measure to NA if recorded 7777 or 9999
paq_f_temp4[which(paq_f_temp4[,"pad630"]==7777 | paq_f_temp4[,"pad630"]==9999),"pad630"] <- NA
# remove anyone whose pad630 measure is outside the possible range of values
paq_f_temp5 <- paq_f_temp4[paq_f_temp4[,"pad630"]>=10 & paq_f_temp4[,"pad630"]<=1080 | is.na(paq_f_temp4[,"pad630"]),]
paq_f_temp6 <- paq_f_temp5
# set pad645 measure to NA if recorded 7777 or 9999
paq_f_temp6[which(paq_f_temp6[,"pad645"]==7777 | paq_f_temp6[,"pad645"]==9999),"pad645"] <- NA
# remove anyone whose pad645 measure is outside the possible range of values
paq_f_temp7 <- paq_f_temp6[paq_f_temp6[,"pad645"]>=10 & paq_f_temp6[,"pad645"]<=960 | is.na(paq_f_temp6[,"pad645"]),]
paq_f_temp8 <- paq_f_temp7
# set paq655 measure to NA if recorded 77 or 99
paq_f_temp8[which(paq_f_temp8[,"paq655"]==77 | paq_f_temp8[,"paq655"]==99),"paq655"] <- NA
# remove anyone whose paq655 measure is outside the possible range of values
paq_f_temp9 <- paq_f_temp8[paq_f_temp8[,"paq655"]>=1 & paq_f_temp8[,"paq655"]<=7 | is.na(paq_f_temp8[,"paq655"]),]
paq_f_temp10 <- paq_f_temp9
# set pad660 measure to NA if recorded 7777 or 9999
paq_f_temp10[which(paq_f_temp10[,"pad660"]==7777 | paq_f_temp10[,"pad660"]==9999),"pad660"] <- NA
# remove anyone whose pad660 measure is outside the possible range of values
paq_f_temp11 <- paq_f_temp10[paq_f_temp10[,"pad660"]>=10 & paq_f_temp10[,"pad660"]<=480 | is.na(paq_f_temp10[,"pad660"]),]
paq_f_temp12 <- paq_f_temp11
# set pad675 measure to NA if recorded 7777 or 9999
paq_f_temp12[which(paq_f_temp12[,"pad675"]==7777 | paq_f_temp12[,"pad675"]==9999),"pad675"] <- NA
# remove anyone whose pad675 measure is outside the possible range of values
paq_f_temp13 <- paq_f_temp12[paq_f_temp12[,"pad675"]>=10 & paq_f_temp12[,"pad675"]<=600 | is.na(paq_f_temp12[,"pad675"]),]

# import physical activity data (paq) from 2011-2012 cycle
paq_g_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/PAQ_G.xpt")
paq_g_temp2 <- paq_g_temp1
# set pad615 measure to NA if recorded 7777 or 9999
paq_g_temp2[which(paq_g_temp2[,"pad615"]==7777 | paq_g_temp2[,"pad615"]==9999),"pad615"] <- NA
# remove anyone whose pad615 measure is outside the possible range of values
paq_g_temp3 <- paq_g_temp2[paq_g_temp2[,"pad615"]>=10 & paq_g_temp2[,"pad615"]<=780 | is.na(paq_g_temp2[,"pad615"]),]
paq_g_temp4 <- paq_g_temp3
# set pad630 measure to NA if recorded 7777 or 9999
paq_g_temp4[which(paq_g_temp4[,"pad630"]==7777 | paq_g_temp4[,"pad630"]==9999),"pad630"] <- NA
# remove anyone whose pad630 measure is outside the possible range of values
paq_g_temp5 <- paq_g_temp4[paq_g_temp4[,"pad630"]>=10 & paq_g_temp4[,"pad630"]<=960 | is.na(paq_g_temp4[,"pad630"]),]
paq_g_temp6 <- paq_g_temp5
# set pad645 measure to NA if recorded 7777 or 9999
paq_g_temp6[which(paq_g_temp6[,"pad645"]==7777 | paq_g_temp6[,"pad645"]==9999),"pad645"] <- NA
# remove anyone whose pad645 measure is outside the possible range of values
paq_g_temp7 <- paq_g_temp6[paq_g_temp6[,"pad645"]>=10 & paq_g_temp6[,"pad645"]<=720 | is.na(paq_g_temp6[,"pad645"]),]
paq_g_temp8 <- paq_g_temp7
# set paq655 measure to NA if recorded 77 or 99
paq_g_temp8[which(paq_g_temp8[,"paq655"]==77 | paq_g_temp8[,"paq655"]==99),"paq655"] <- NA
# remove anyone whose paq655 measure is outside the possible range of values
paq_g_temp9 <- paq_g_temp8[paq_g_temp8[,"paq655"]>=1 & paq_g_temp8[,"paq655"]<=7 | is.na(paq_g_temp8[,"paq655"]),]
paq_g_temp10 <- paq_g_temp9
# set pad660 measure to NA if recorded 7777 or 9999
paq_g_temp10[which(paq_g_temp10[,"pad660"]==7777 | paq_g_temp10[,"pad660"]==9999),"pad660"] <- NA
# remove anyone whose pad660 measure is outside the possible range of values
paq_g_temp11 <- paq_g_temp10[paq_g_temp10[,"pad660"]>=10 & paq_g_temp10[,"pad660"]<=720 | is.na(paq_g_temp10[,"pad660"]),]
paq_g_temp12 <- paq_g_temp11
# set pad675 measure to NA if recorded 77 or 99
paq_g_temp12[which(paq_g_temp12[,"pad675"]==77 | paq_g_temp12[,"pad675"]==99),"pad675"] <- NA
# remove anyone whose pad675 measure is outside the possible range of values
paq_g_temp13 <- paq_g_temp12[paq_g_temp12[,"pad675"]>=10 & paq_g_temp12[,"pad675"]<=600 | is.na(paq_g_temp12[,"pad675"]),]

# combine physical activity data (paq) from 2007-2008, 2009-2010 and 2011-2012 cycles
paq2007_temp1 <- rbind(paq_e_temp13[,Paq2007Variables], paq_f_temp13[,Paq2007Variables], paq_g_temp13[,Paq2007Variables])
paq2007_temp2 <- paq2007_temp1
# set paq605 measure to NA if recorded 7 or 9
paq2007_temp2[which(paq2007_temp2[,"paq605"]==7 | paq2007_temp2[,"paq605"]==9),"paq605"] <- NA
paq2007_temp3 <- paq2007_temp2[paq2007_temp2[,"paq605"]==1 | paq2007_temp2[,"paq605"]==2 | is.na(paq2007_temp2[,"paq605"]),]
paq2007_temp4 <- paq2007_temp3
# set paq610 measure to NA if recorded 77 or 99
paq2007_temp4[which(paq2007_temp4[,"paq610"]==77 | paq2007_temp4[,"paq610"]==99),"paq610"] <- NA
paq2007_temp5 <- paq2007_temp4[paq2007_temp4[,"paq610"]>=1 & paq2007_temp4[,"paq610"]<=7 | is.na(paq2007_temp4[,"paq610"]),]
paq2007_temp6 <- paq2007_temp5
# set paq620 measure to NA if recorded 7 or 9
paq2007_temp6[which(paq2007_temp6[,"paq620"]==7 | paq2007_temp6[,"paq620"]==9),"paq620"] <- NA
paq2007_temp7 <- paq2007_temp6[paq2007_temp6[,"paq620"]==1 | paq2007_temp6[,"paq620"]==2 | is.na(paq2007_temp6[,"paq620"]),]
paq2007_temp8 <- paq2007_temp7
# set paq625 measure to NA if recorded 77 or 99
paq2007_temp8[which(paq2007_temp8[,"paq625"]==77 | paq2007_temp8[,"paq625"]==99),"paq625"] <- NA
# remove anyone whose paq625 measure is outside the possible range of values
paq2007_temp9 <- paq2007_temp8[paq2007_temp8[,"paq625"]>=1 & paq2007_temp8[,"paq625"]<=7 | is.na(paq2007_temp8[,"paq625"]),]
paq2007_temp10 <- paq2007_temp9
# set paq635 measure to NA if recorded 7 or 9
paq2007_temp10[which(paq2007_temp10[,"paq635"]==7 | paq2007_temp10[,"paq635"]==9),"paq635"] <- NA
paq2007_temp11 <- paq2007_temp10[paq2007_temp10[,"paq635"]==1 | paq2007_temp10[,"paq635"]==2 | is.na(paq2007_temp10[,"paq635"]),]
paq2007_temp12 <- paq2007_temp11
# set paq640 measure to NA if recorded 77 or 99
paq2007_temp12[which(paq2007_temp12[,"paq640"]==77 | paq2007_temp12[,"paq640"]==99),"paq640"] <- NA
# remove anyone whose paq640 measure is outside the possible range of values
paq2007_temp13 <- paq2007_temp12[paq2007_temp12[,"paq640"]>=1 & paq2007_temp12[,"paq640"]<=7 | is.na(paq2007_temp12[,"paq640"]),]
paq2007_temp14 <- paq2007_temp13
# set paq650 measure to NA if recorded 7 or 9
paq2007_temp14[which(paq2007_temp14[,"paq650"]==7 | paq2007_temp14[,"paq650"]==9),"paq650"] <- NA
paq2007_temp15 <- paq2007_temp14[paq2007_temp14[,"paq650"]==1 | paq2007_temp14[,"paq650"]==2 | is.na(paq2007_temp14[,"paq650"]),]
paq2007_temp16 <- paq2007_temp15
# set paq665 measure to NA if recorded 7 or 9
paq2007_temp16[which(paq2007_temp16[,"paq665"]==7 | paq2007_temp16[,"paq665"]==9),"paq665"] <- NA
paq2007_temp17 <- paq2007_temp16[paq2007_temp16[,"paq665"]==1 | paq2007_temp16[,"paq665"]==2 | is.na(paq2007_temp16[,"paq665"]),]
paq2007_temp18 <- paq2007_temp17
# set paq670 measure to NA if recorded 77 or 99
paq2007_temp18[which(paq2007_temp18[,"paq670"]==77 | paq2007_temp18[,"paq670"]==99),"paq670"] <- NA
# remove anyone whose paq670 measure is outside the possible range of values
paq2007_temp19 <- paq2007_temp18[paq2007_temp18[,"paq670"]>=1 & paq2007_temp18[,"paq670"]<=7 | is.na(paq2007_temp18[,"paq670"]),]
# set paq2007_temp19 to paq2007 which may still include missing (NA) measurements
paq2007 <- paq2007_temp19

# IMPORT TOTAL NUTRIENT INTAKE DATA - FIRST DAY
Dr1totVariables <- c("seqn", "drqsdiet", "drqsdt1", "drqsdt8", "dr1day", "dr1tkcal")

# import total nutrient intake data - first day from 2003-2004 cycle
dr1tot_c_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/DR1TOT_C.xpt")
# remove anyone whose dr1tkcal measure is outside the possible range of values
dr1tot_c_temp2 <- dr1tot_c_temp1[dr1tot_c_temp1[,"dr1tkcal"]>=0 & dr1tot_c_temp1[,"dr1tkcal"]<=9724 | is.na(dr1tot_c_temp1[,"dr1tkcal"]),]

# import total nutrient intake data - first day from 2005-2006 cycle
dr1tot_d_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/DR1TOT_D.xpt")
# remove anyone whose dr1tkcal measure is outside the possible range of values
dr1tot_d_temp2 <- dr1tot_d_temp1[dr1tot_d_temp1[,"dr1tkcal"]>=0 & dr1tot_d_temp1[,"dr1tkcal"]<=10569 | is.na(dr1tot_d_temp1[,"dr1tkcal"]),]

# import total nutrient intake data - first day from 2007-2008 cycle
dr1tot_e_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2007/DataFiles/DR1TOT_E.xpt")
# remove anyone whose dr1tkcal measure is outside the possible range of values
dr1tot_e_temp2 <- dr1tot_e_temp1[dr1tot_e_temp1[,"dr1tkcal"]>=0 & dr1tot_e_temp1[,"dr1tkcal"]<=13509 | is.na(dr1tot_e_temp1[,"dr1tkcal"]),]

# import total nutrient intake data - first day from 2009-2010 cycle
dr1tot_f_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/DR1TOT_F.xpt")
# remove anyone whose dr1tkcal measure is outside the possible range of values
dr1tot_f_temp2 <- dr1tot_f_temp1[dr1tot_f_temp1[,"dr1tkcal"]>=55 & dr1tot_f_temp1[,"dr1tkcal"]<=10463 | is.na(dr1tot_f_temp1[,"dr1tkcal"]),]

# import total nutrient intake data - first day from 2011-2012 cycle
dr1tot_g_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/DR1TOT_G.xpt")
# remove anyone whose dr1tkcal measure is outside the possible range of values
dr1tot_g_temp2 <- dr1tot_g_temp1[dr1tot_g_temp1[,"dr1tkcal"]>=0 & dr1tot_g_temp1[,"dr1tkcal"]<=13687 | is.na(dr1tot_g_temp1[,"dr1tkcal"]),]

# combine total nutrient intake data - first day from all survey cycles
dr1tot_temp1 <- rbind(dr1tot_c_temp2[,Dr1totVariables], dr1tot_d_temp2[,Dr1totVariables], dr1tot_e_temp2[,Dr1totVariables], dr1tot_f_temp2[,Dr1totVariables], dr1tot_g_temp2[,Dr1totVariables])
# remove anyone with dr1day measure is outside the possible range of values
dr1tot_temp2 <- dr1tot_temp1[dr1tot_temp1[,"dr1day"]==1 | dr1tot_temp1[,"dr1day"]==2 | dr1tot_temp1[,"dr1day"]==3 | dr1tot_temp1[,"dr1day"]==4 | dr1tot_temp1[,"dr1day"]==5 | dr1tot_temp1[,"dr1day"]==6 | dr1tot_temp1[,"dr1day"]==7 | is.na(dr1tot_temp1[,"dr1day"]),]
# set dr1tot_temp2 to dr1tot which may still include missing (NA) measurements
dr1tot <- dr1tot_temp2

# IMPORT TOTAL NUTRIENT INTAKE DATA - SECOND DAY
Dr2totVariables <- c("seqn", "dr2day", "dr2tkcal")

# import total nutrient intake data - second day from 2003-2004 cycle
dr2tot_c_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/DR2TOT_C.xpt")
# remove anyone whose dr2tkcal measure is outside the possible range of values
dr2tot_c_temp2 <- dr2tot_c_temp1[dr2tot_c_temp1[,"dr2tkcal"]>=0 & dr2tot_c_temp1[,"dr2tkcal"]<=9265 | is.na(dr2tot_c_temp1[,"dr2tkcal"]),]

# import total nutrient intake data - second day from 2005-2006 cycle
dr2tot_d_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/DR2TOT_D.xpt")
# remove anyone whose dr2tkcal measure is outside the possible range of values
dr2tot_d_temp2 <- dr2tot_d_temp1[dr2tot_d_temp1[,"dr2tkcal"]>=0 & dr2tot_d_temp1[,"dr2tkcal"]<=9119 | is.na(dr2tot_d_temp1[,"dr2tkcal"]),]

# import total nutrient intake data - second day from 2007-2008 cycle
dr2tot_e_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2007/DataFiles/DR2TOT_E.xpt")
# remove anyone whose dr2tkcal measure is outside the possible range of values
dr2tot_e_temp2 <- dr2tot_e_temp1[dr2tot_e_temp1[,"dr2tkcal"]>=0 & dr2tot_e_temp1[,"dr2tkcal"]<=9608 | is.na(dr2tot_e_temp1[,"dr2tkcal"]),]

# import total nutrient intake data - second day from 2009-2010 cycle
dr2tot_f_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/DR2TOT_F.xpt")
# remove anyone whose dr2tkcal measure is outside the possible range of values
dr2tot_f_temp2 <- dr2tot_f_temp1[dr2tot_f_temp1[,"dr2tkcal"]>=0 & dr2tot_f_temp1[,"dr2tkcal"]<=8976 | is.na(dr2tot_f_temp1[,"dr2tkcal"]),]

# import total nutrient intake data - second day from 2011-2012 cycle
dr2tot_g_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/DR2TOT_G.xpt")
# remove anyone whose dr2tkcal measure is outside the possible range of values
dr2tot_g_temp2 <- dr2tot_g_temp1[dr2tot_g_temp1[,"dr2tkcal"]>=0 & dr2tot_g_temp1[,"dr2tkcal"]<=8213 | is.na(dr2tot_g_temp1[,"dr2tkcal"]),]

# combine total nutrient intake data - second day from all survey cycles
dr2tot_temp1 <- rbind(dr2tot_c_temp2[,Dr2totVariables], dr2tot_d_temp2[,Dr2totVariables], dr2tot_e_temp2[,Dr2totVariables], dr2tot_f_temp2[,Dr2totVariables], dr2tot_g_temp2[,Dr2totVariables])
# remove anyone with dr2day measure is outside the possible range of values
dr2tot_temp2 <- dr2tot_temp1[dr2tot_temp1[,"dr2day"]==1 | dr2tot_temp1[,"dr2day"]==2 | dr2tot_temp1[,"dr2day"]==3 | dr2tot_temp1[,"dr2day"]==4 | dr2tot_temp1[,"dr2day"]==5 | dr2tot_temp1[,"dr2day"]==6 | dr2tot_temp1[,"dr2day"]==7 | is.na(dr2tot_temp1[,"dr2day"]),]
# set dr2tot_temp2 to dr2tot which may still include missing (NA) measurements
dr2tot <- dr2tot_temp2

# IMPORT INDIVIDUAL FOODS DATA - FIRST DAY
Dr1iffVariables <- c("seqn", "dr1ifdcd", "dr1ikcal")

# import individual foods data - first day from 2003-2004 cycle
dr1iff_c_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/DR1IFF_C.xpt")
# remove anyone whose dr1ifdcd measure is outside the possible range of values
dr1iff_c_temp2 <- dr1iff_c_temp1[dr1iff_c_temp1[,"dr1ifdcd"]>=11000000 & dr1iff_c_temp1[,"dr1ifdcd"]<=94210100 | is.na(dr1iff_c_temp1[,"dr1ifdcd"]),]
# remove anyone whose dr1ikcal measure is outside the possible range of values
dr1iff_c_temp3 <- dr1iff_c_temp2[dr1iff_c_temp2[,"dr1ikcal"]>=0 & dr1iff_c_temp2[,"dr1ikcal"]<=4114 | is.na(dr1iff_c_temp2[,"dr1ikcal"]),]

# import individual foods data - first day from 2005-2006 cycle
dr1iff_d_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/DR1IFF_D.xpt")
# remove anyone whose dr1ifdcd measure is outside the possible range of values
dr1iff_d_temp2 <- dr1iff_d_temp1[dr1iff_d_temp1[,"dr1ifdcd"]>=11000000 & dr1iff_d_temp1[,"dr1ifdcd"]<=94210200 | is.na(dr1iff_d_temp1[,"dr1ifdcd"]),]
# remove anyone whose dr1ikcal measure is outside the possible range of values
dr1iff_d_temp3 <- dr1iff_d_temp2[dr1iff_d_temp2[,"dr1ikcal"]>=0 & dr1iff_d_temp2[,"dr1ikcal"]<=5559 | is.na(dr1iff_d_temp2[,"dr1ikcal"]),]

# import individual foods data - first day from 2007-2008 cycle
dr1iff_e_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2007/DataFiles/DR1IFF_E.xpt")
# remove anyone whose dr1ifdcd measure is outside the possible range of values
dr1iff_e_temp2 <- dr1iff_e_temp1[dr1iff_e_temp1[,"dr1ifdcd"]>=11000000 & dr1iff_e_temp1[,"dr1ifdcd"]<=94210200 | is.na(dr1iff_e_temp1[,"dr1ifdcd"]),]
# remove anyone whose dr1ikcal measure is outside the possible range of values
dr1iff_e_temp3 <- dr1iff_e_temp2[dr1iff_e_temp2[,"dr1ikcal"]>=0 & dr1iff_e_temp2[,"dr1ikcal"]<=5662 | is.na(dr1iff_e_temp2[,"dr1ikcal"]),]

# import individual foods data - first day from 2009-2010 cycle
dr1iff_f_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/DR1IFF_F.xpt")
# remove anyone whose dr1ifdcd measure is outside the possible range of values
dr1iff_f_temp2 <- dr1iff_f_temp1[dr1iff_f_temp1[,"dr1ifdcd"]>=11000000 & dr1iff_f_temp1[,"dr1ifdcd"]<=94300100 | is.na(dr1iff_f_temp1[,"dr1ifdcd"]),]
# remove anyone whose dr1ikcal measure is outside the possible range of values
dr1iff_f_temp3 <- dr1iff_f_temp2[dr1iff_f_temp2[,"dr1ikcal"]>=0 & dr1iff_f_temp2[,"dr1ikcal"]<=4632 | is.na(dr1iff_f_temp2[,"dr1ikcal"]),]

# import individual foods data - first day from 2011-2012 cycle
dr1iff_g_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/DR1IFF_G.xpt")
# remove anyone whose dr1ifdcd measure is outside the possible range of values
dr1iff_g_temp2 <- dr1iff_g_temp1[dr1iff_g_temp1[,"dr1ifdcd"]>=11000000 & dr1iff_g_temp1[,"dr1ifdcd"]<=95330500 | is.na(dr1iff_g_temp1[,"dr1ifdcd"]),]
# remove anyone whose dr1ikcal measure is outside the possible range of values
dr1iff_g_temp3 <- dr1iff_g_temp2[dr1iff_g_temp2[,"dr1ikcal"]>=0 & dr1iff_g_temp2[,"dr1ikcal"]<=4128 | is.na(dr1iff_g_temp2[,"dr1ikcal"]),]

# combine individual foods data - first day from all survey cycles
dr1iff_temp1 <- rbind(dr1iff_c_temp3[,Dr1iffVariables], dr1iff_d_temp3[,Dr1iffVariables], dr1iff_e_temp3[,Dr1iffVariables], dr1iff_f_temp3[,Dr1iffVariables], dr1iff_g_temp3[,Dr1iffVariables])
# set dr1iff_temp1 to dr1iff which may still include missing (NA) measurements
dr1iff <- dr1iff_temp1

# IMPORT INDIVIDUAL FOODS DATA - SECOND DAY
Dr2iffVariables <- c("seqn", "dr2ifdcd", "dr2ikcal")

# import individual foods data - second day from 2003-2004 cycle
dr2iff_c_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/DR2IFF_C.xpt")
# remove anyone whose dr2ifdcd measure is outside the possible range of values
dr2iff_c_temp2 <- dr2iff_c_temp1[dr2iff_c_temp1[,"dr2ifdcd"]>=11000000 & dr2iff_c_temp1[,"dr2ifdcd"]<=94210100 | is.na(dr2iff_c_temp1[,"dr2ifdcd"]),]
# remove anyone whose dr2ikcal measure is outside the possible range of values
dr2iff_c_temp3 <- dr2iff_c_temp2[dr2iff_c_temp2[,"dr2ikcal"]>=0 & dr2iff_c_temp2[,"dr2ikcal"]<=4480 | is.na(dr2iff_c_temp2[,"dr2ikcal"]),]

# import individual foods data - second day from 2005-2006 cycle
dr2iff_d_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/DR2IFF_D.xpt")
# remove anyone whose dr2ifdcd measure is outside the possible range of values
dr2iff_d_temp2 <- dr2iff_d_temp1[dr2iff_d_temp1[,"dr2ifdcd"]>=11000000 & dr2iff_d_temp1[,"dr2ifdcd"]<=94210200 | is.na(dr2iff_d_temp1[,"dr2ifdcd"]),]
# remove anyone whose dr2ikcal measure is outside the possible range of values
dr2iff_d_temp3 <- dr2iff_d_temp2[dr2iff_d_temp2[,"dr2ikcal"]>=0 & dr2iff_d_temp2[,"dr2ikcal"]<=3738 | is.na(dr2iff_d_temp2[,"dr2ikcal"]),]

# import individual foods data - second day from 2007-2008 cycle
dr2iff_e_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2007/DataFiles/DR2IFF_E.xpt")
# remove anyone whose dr2ifdcd measure is outside the possible range of values
dr2iff_e_temp2 <- dr2iff_e_temp1[dr2iff_e_temp1[,"dr2ifdcd"]>=11000000 & dr2iff_e_temp1[,"dr2ifdcd"]<=94210200 | is.na(dr2iff_e_temp1[,"dr2ifdcd"]),]
# remove anyone whose dr2ikcal measure is outside the possible range of values
dr2iff_e_temp3 <- dr2iff_e_temp2[dr2iff_e_temp2[,"dr2ikcal"]>=0 & dr2iff_e_temp2[,"dr2ikcal"]<=5056 | is.na(dr2iff_e_temp2[,"dr2ikcal"]),]

# import individual foods data - second day from 2009-2010 cycle
dr2iff_f_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/DR2IFF_F.xpt")
# remove anyone whose dr2ifdcd measure is outside the possible range of values
dr2iff_f_temp2 <- dr2iff_f_temp1[dr2iff_f_temp1[,"dr2ifdcd"]>=11000000 & dr2iff_f_temp1[,"dr2ifdcd"]<=94300100 | is.na(dr2iff_f_temp1[,"dr2ifdcd"]),]
# remove anyone whose dr2ikcal measure is outside the possible range of values
dr2iff_f_temp3 <- dr2iff_f_temp2[dr2iff_f_temp2[,"dr2ikcal"]>=0 & dr2iff_f_temp2[,"dr2ikcal"]<=5081 | is.na(dr2iff_f_temp2[,"dr2ikcal"]),]

# import individual foods data - second day from 2011-2012 cycle
dr2iff_g_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/DR2IFF_G.xpt")
# remove anyone whose dr2ifdcd measure is outside the possible range of values
dr2iff_g_temp2 <- dr2iff_g_temp1[dr2iff_g_temp1[,"dr2ifdcd"]>=11000000 & dr2iff_g_temp1[,"dr2ifdcd"]<=95342000 | is.na(dr2iff_g_temp1[,"dr2ifdcd"]),]
# remove anyone whose dr2ikcal measure is outside the possible range of values
dr2iff_g_temp3 <- dr2iff_g_temp2[dr2iff_g_temp2[,"dr2ikcal"]>=0 & dr2iff_g_temp2[,"dr2ikcal"]<=5443 | is.na(dr2iff_g_temp2[,"dr2ikcal"]),]

# combine individual foods data - second day from all survey cycles
dr2iff_temp1 <- rbind(dr2iff_c_temp2[,Dr2iffVariables], dr2iff_d_temp2[,Dr2iffVariables], dr2iff_e_temp2[,Dr2iffVariables], dr2iff_f_temp2[,Dr2iffVariables], dr2iff_g_temp2[,Dr2iffVariables])
# set dr2iff_temp1 to dr2iff which may still include missing (NA) measurements
dr2iff <- dr2iff_temp1

# IMPORT DIETARY INTERVIEW TECHNICAL SUPPORT
# import dietary interview technical support from 2003-2004 cycle
drxfcd_c_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/DRXFCD_C.xpt")
# remove anyone whose drxfdcd measure is outside the possible range of values
drxfcd_c_temp2 <- drxfcd_c_temp1[drxfcd_c_temp1[,"drxfdcd"]>=11000000 & drxfcd_c_temp1[,"drxfdcd"]<=94210100 | is.na(drxfcd_c_temp1[,"drxfdcd"]),]

# import dietary interview technical support from 2005-2006 cycle
drxfcd_d_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/DRXFCD_D.xpt")
# remove anyone whose drxfdcd measure is outside the possible range of values
drxfcd_d_temp2 <- drxfcd_d_temp1[drxfcd_d_temp1[,"drxfdcd"]>=11000000 & drxfcd_d_temp1[,"drxfdcd"]<=94210200 | is.na(drxfcd_d_temp1[,"drxfdcd"]),]

# import dietary interview technical support from 2007-2008 cycle
drxfcd_e_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2007/DataFiles/DRXFCD_E.xpt")
# remove anyone whose drxfdcd measure is outside the possible range of values
drxfcd_e_temp2 <- drxfcd_e_temp1[drxfcd_e_temp1[,"drxfdcd"]>=11000000 & drxfcd_e_temp1[,"drxfdcd"]<=94210200 | is.na(drxfcd_e_temp1[,"drxfdcd"]),]

# import dietary interview technical support from 2009-2010 cycle
drxfcd_f_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/DRXFCD_F.xpt")
# remove anyone whose drxfdcd measure is outside the possible range of values
drxfcd_f_temp2 <- drxfcd_f_temp1[drxfcd_f_temp1[,"drxfdcd"]>=11000000 & drxfcd_f_temp1[,"drxfdcd"]<=94300100 | is.na(drxfcd_f_temp1[,"drxfdcd"]),]

# import dietary interview technical support from 2011-2012 cycle
drxfcd_g_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/DRXFCD_G.xpt")
# remove anyone whose drxfdcd measure is outside the possible range of values
drxfcd_g_temp2 <- drxfcd_g_temp1[drxfcd_g_temp1[,"drxfdcd"]>=11000000 & drxfcd_g_temp1[,"drxfdcd"]<=95342000 | is.na(drxfcd_g_temp1[,"drxfdcd"]),]

# IMPORT REPRODUCTIVE HEALTH DATA
RhqVariables <- c("seqn", "rhd143", "rhq200")

# import reproductive health data from 2003-2004 cycle
rhq_c_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/RHQ_C.xpt")
# import reproductive health data from 2005-2006 cycle
rhq_d_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/RHQ_D.xpt")
# import reproductive health data from 2007-2008 cycle
rhq_e_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2007/DataFiles/RHQ_E.xpt")
# import reproductive health data from 2009-2010 cycle
rhq_f_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/RHQ_F.xpt")
# import reproductive health data from 2011-2012 cycle
rhq_g_temp1 <- sasxport.get(file="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/RHQ_G.xpt")

# combine body measures from all survey cycles
rhq_temp1 <- rbind(rhq_c_temp1[,RhqVariables], rhq_d_temp1[,RhqVariables], rhq_e_temp1[,RhqVariables], rhq_f_temp1[,RhqVariables], rhq_g_temp1[,RhqVariables])
# set rhq_temp1 to rhq
rhq <- rhq_temp1

################################################################################
# C. Cleans data based on Butler et al 2016
################################################################################

# COMBINE ALL TYPES OF DATA
# paqiaf, paq2003 and paq2007 are subsamples
# seqn is not unique in paqiaf as each seqn can perform more than one type of physical activity
# seqn is not unique in dr1iff and dr2iff as each seqn can consume more than one food item in a day
# seqn is unique and measured for all individuals in body measures, alcohol use data, demographic variables, smoking data, total nutrient intake data - first day and total nutrient intake data - second day
# create nhanes by combining them
# rhq is not included in nhanes because it is only used for exclusions
nhanes <- merge(merge(merge(merge(merge(bmx, alq, by="seqn"), demo, by="seqn"), smq, by="seqn"), dr1tot, by="seqn"), dr2tot, by="seqn")
# produce a vector of seqn that does not appear in all data frames
Seqn_temp1 <- setdiff(c(nhanes$seqn, paqiaf$seqn, paq2003$seqn, paq2007$seqn, dr1iff$seqn, dr2iff$seqn), intersect(intersect(intersect(nhanes$seqn, c(paqiaf$seqn, paq2003$seqn, paq2007$seqn)), dr1iff$seqn), dr2iff$seqn))
# remove participants that do not appear in all data frames
sample_temp1 <- RemoveSeqn(list("nhanes"=nhanes, "paqiaf"=paqiaf, "paq2003"=paq2003, "paq2007"=paq2007, "dr1iff"=dr1iff, "dr2iff"=dr2iff), Seqn_temp1)

# REMOVE ANYONE WHO IS NOT BETWEEN 20 AND 79 YEARS OLD
# remove anyone not ridageyr>=20 and ridageyr<=79
seqn_temp2 <- sample_temp1$nhanes$seqn[is.na(sample_temp1$nhanes$ridageyr) | sample_temp1$nhanes$ridageyr<20 | sample_temp1$nhanes$ridageyr>79]
sample_temp2 <- RemoveSeqn(sample_temp1, seqn_temp2)

# REMOVE ANYONE WITH MISSING ALCOHOLIC USE DATA
# remove anyone with NA to all alq101, alq110, alq120q and alq130
seqn_temp3a <- sample_temp2$nhanes$seqn[rowSums(!is.na(sample_temp2$nhanes[,c("alq101", "alq110", "alq120q", "alq130")]))==0]
sample_temp3a <- RemoveSeqn(sample_temp2, seqn_temp3a)

# remove anyone with NA to all alq110, alq120q and alq130 when alq101 is 1
seqn_temp3b <- sample_temp3a$nhanes$seqn[sample_temp3a$nhanes[,"alq101"]==1 & rowSums(!is.na(sample_temp3a$nhanes[,c("alq110", "alq120q", "alq130")]))==0]
sample_temp3b <- RemoveSeqn(sample_temp3a, seqn_temp3b)

# remove anyone with NA to all alq110, alq120q and alq130 when alq101 is 2
seqn_temp3c <- sample_temp3b$nhanes$seqn[sample_temp3b$nhanes[,"alq101"]==2 & rowSums(!is.na(sample_temp3b$nhanes[,c("alq110", "alq120q", "alq130")]))==0]
sample_temp3c <- RemoveSeqn(sample_temp3b, seqn_temp3c)

# remove anyone with NA to alq130 when alq101 or alq110 is 1 and alq120q is at least 1
seqn_temp3d <- sample_temp3c$nhanes$seqn[which((sample_temp3c$nhanes[,"alq101"]==1 | sample_temp3c$nhanes[,"alq110"]==1) & sample_temp3c$nhanes[,"alq120q"]>=1 & is.na(sample_temp3c$nhanes[,"alq130"]))]
sample_temp3d <- RemoveSeqn(sample_temp3c, seqn_temp3d)

# remove anyone with NA to all alq120q and alq130 when alq101 is 2 and alq110 is 1
seqn_temp3 <- sample_temp3d$nhanes$seqn[which(sample_temp3d$nhanes[,"alq101"]==2 & sample_temp3d$nhanes[,"alq110"]==1 & is.na(sample_temp3d$nhanes[,"alq120q"]) & is.na(sample_temp3d$nhanes[,"alq130"]))]
sample_temp3 <- RemoveSeqn(sample_temp3d, seqn_temp3)

# REMOVE ANYONE WITH MISSING TOTAL NUTRIENT INTAKE DATA - FIRST DAY
# remove anyone with NA to dr1tkcal
seqn_temp4 <- sample_temp3$nhanes$seqn[is.na(sample_temp3$nhanes$dr1tkcal)]
sample_temp4 <- RemoveSeqn(sample_temp3, seqn_temp4)

# REMOVE ANYONE WITH MISSING TOTAL NUTRIENT INTAKE DATA - SECOND DAY
# remove anyone with NA to dr2tkcal
seqn_temp5 <- sample_temp4$nhanes$seqn[is.na(sample_temp4$nhanes$dr2tkcal)]
sample_temp5 <- RemoveSeqn(sample_temp4, seqn_temp5)

# REMOVE ANYONE WITH MISSING INDIVIDUAL FOODS DATA - FIRST DAY
# remove anyone with NA to dr1ikcal
seqn_temp6 <- sample_temp5$dr1iff$seqn[is.na(sample_temp5$dr1iff$dr1ikcal)]
sample_temp6 <- RemoveSeqn(sample_temp5, seqn_temp6)

# REMOVE ANYONE WITH MISSING INDIVIDUAL FOODS DATA - SECOND DAY
# remove anyone with NA to dr2ikcal
seqn_temp7 <- sample_temp6$dr2iff$seqn[is.na(sample_temp6$dr2iff$dr2ikcal)]
sample_temp7 <- RemoveSeqn(sample_temp6, seqn_temp7)

# REMOVE ANYONE WITH MISSING BODY MEASURES
# remove anyone with NA to bmxwaist
seqn_temp8 <- sample_temp7$nhanes$seqn[is.na(sample_temp7$nhanes[,"bmxwaist"])]
sample_temp8 <- RemoveSeqn(sample_temp7, seqn_temp8)

# REMOVE ANYONE WHO IS CURRENTLY PREGNANT OR BREASTFEEDING
# remove anyone if rhd143 is 1
# note we don't remove anyone whose rhd143 is missing or responded "don't know" (implicitly assuming they are not pregnant)
seqn_temp9a <- rhq$seqn[which(rhq[,"rhd143"]==1)]
sample_temp9a <- RemoveSeqn(sample_temp8, seqn_temp9a)

# remove anyone if rhq200 is 1
# note we don't remove anyone whose rhq200 is missing or responded "don't know" (implicitly assuming they are not breastfeeding)
seqn_temp9 <- rhq$seqn[which(rhq[,"rhq200"]==1)]
sample_temp9 <- RemoveSeqn(sample_temp9a, seqn_temp9)

# REMOVE ANYONE WHO IS CURRENTLY ON A MEDICAL DIET
# remove anyone if drqsdiet is 1 but not drqsdt1 is 1 nor drqsdt8 is 8
# note we don't remove anyone whose drqsdiet is missing or responded "don't know" (implicitly assuming they are not on a medical diet)
seqn_temp10 <- sample_temp9$nhanes$seqn[setdiff(which(sample_temp9$nhanes[,"drqsdiet"]==1), c(which(sample_temp9$nhanes[,"drqsdt1"]==1), which(sample_temp9$nhanes[,"drqsdt8"]==8)))]
sample_temp10 <- RemoveSeqn(sample_temp9, seqn_temp10)

# REMOVE ANYONE WHO IS CURRENTLY ON A WEIGHT-LOSS DIET
# remove anyone if drqsdt1 is 1
# note we don't remove anyone whose drqsdt1 is missing (implicitly assuming they are not on a weight-loss diet)
seqn_temp11a <- sample_temp10$nhanes$seqn[which(sample_temp10$nhanes[,"drqsdt1"]==1)]
sample_temp11a <- RemoveSeqn(sample_temp10, seqn_temp11a)

# REMOVE ANYONE WHO IS CURRENTLY ON A WEIGHT-GAIN DIET
# remove anyone if drqsdt8 is 8
# note we don't remove anyone whose drqsdt8 is missing (implicitly assuming they are not on a weight-gain diet)
seqn_temp11 <- sample_temp11a$nhanes$seqn[which(sample_temp11a$nhanes[,"drqsdt8"]==8)]
sample_temp11 <- RemoveSeqn(sample_temp11a, seqn_temp11)

# REMOVE ANYONE WITH MISSING DEMOGRAPHIC VARIABLES
# remove anyone with NA to ridageyr
seqn_temp12 <- sample_temp11$nhanes$seqn[is.na(sample_temp11$nhanes$ridageyr)]
sample_temp12 <- RemoveSeqn(sample_temp11, seqn_temp12)

# remove anyone with NA to ridreth1
seqn_temp13 <- sample_temp12$nhanes$seqn[is.na(sample_temp12$nhanes$ridreth1)]
sample_temp13 <- RemoveSeqn(sample_temp12, seqn_temp13)

# remove anyone with NA to dmdeduc2
seqn_temp14 <- sample_temp13$nhanes$seqn[is.na(sample_temp13$nhanes$dmdeduc2)]
sample_temp14 <- RemoveSeqn(sample_temp13, seqn_temp14)

# REMOVE ANYONE WITH MISSING PHYSICAL ACTIVITY DATA
# remove anyone with NA to padtimes, paddurat or padmets
# note that seqn is not necessary unique in paqiaf but we assume if one recorded activity is unreliable then total activity will be unreliable
# so participants with at least one unreliable activity should be removed from the sample altogether
seqn_temp15a <- sample_temp14$paqiaf$seqn[is.na(sample_temp14$paqiaf$padtimes) | is.na(sample_temp14$paqiaf$paddurat) | is.na(sample_temp14$paqiaf$padmets)]
sample_temp15a <- RemoveSeqn(sample_temp14, seqn_temp15a)

# remove anyone with NA to pad020
seqn_temp15b <- sample_temp15a$paq2003$seqn[is.na(sample_temp15a$paq2003$pad020)]
sample_temp15b <- RemoveSeqn(sample_temp15a, seqn_temp15b)

# remove anyone with NA to at least one of paq050q, paq050u and pad080 when pad020 is 1
seqn_temp15c <- sample_temp15b$paq2003$seqn[sample_temp15b$paq2003$pad020==1 & (is.na(sample_temp15b$paq2003$paq050q) | is.na(sample_temp15b$paq2003$paq050u) | is.na(sample_temp15b$paq2003$pad080))]
sample_temp15c <- RemoveSeqn(sample_temp15b, seqn_temp15c)

# remove anyone with NA to paq100
seqn_temp15d <- sample_temp15c$paq2003$seqn[is.na(sample_temp15c$paq2003$paq100)]
sample_temp15d <- RemoveSeqn(sample_temp15c, seqn_temp15d)

# remove anyone with NA to at least one of pad120 and pad160 when paq100 is 1
seqn_temp15e <- sample_temp15d$paq2003$seqn[sample_temp15d$paq2003$paq100==1 & (is.na(sample_temp15d$paq2003$pad120) | is.na(sample_temp15d$paq2003$pad160))]
sample_temp15e <- RemoveSeqn(sample_temp15d, seqn_temp15e)

# remove anyone with NA to paq605
seqn_temp15f <- sample_temp15e$paq2007$seqn[is.na(sample_temp15e$paq2007$paq605)]
sample_temp15f <- RemoveSeqn(sample_temp15e, seqn_temp15f)

# remove anyone with NA to at least one of paq610 and pad615 when paq605 is 1
seqn_temp15g <- sample_temp15f$paq2007$seqn[sample_temp15f$paq2007$paq605==1 & (is.na(sample_temp15f$paq2007$paq610) | is.na(sample_temp15f$paq2007$pad615))]
sample_temp15g <- RemoveSeqn(sample_temp15f, seqn_temp15g)

# remove anyone with NA to paq620
seqn_temp15h <- sample_temp15g$paq2007$seqn[is.na(sample_temp15g$paq2007$paq620)]
sample_temp15h <- RemoveSeqn(sample_temp15g, seqn_temp15h)

# remove anyone with NA to at least one of paq625 and pad630 when paq620 is 1
seqn_temp15i <- sample_temp15h$paq2007$seqn[sample_temp15h$paq2007$paq620==1 & (is.na(sample_temp15h$paq2007$paq625) | is.na(sample_temp15h$paq2007$pad630))]
sample_temp15i <- RemoveSeqn(sample_temp15h, seqn_temp15i)

# remove anyone with NA to paq635
seqn_temp15j <- sample_temp15i$paq2007$seqn[is.na(sample_temp15i$paq2007$paq635)]
sample_temp15j <- RemoveSeqn(sample_temp15i, seqn_temp15j)

# remove anyone with NA to at least one of paq640 and pad645 when paq635 is 1
seqn_temp15k <- sample_temp15j$paq2007$seqn[sample_temp15j$paq2007$paq635==1 & (is.na(sample_temp15j$paq2007$paq640) | is.na(sample_temp15j$paq2007$pad645))]
sample_temp15k <- RemoveSeqn(sample_temp15j, seqn_temp15k)

# remove anyone with NA to paq650
seqn_temp15l <- sample_temp15k$paq2007$seqn[is.na(sample_temp15k$paq2007$paq650)]
sample_temp15l <- RemoveSeqn(sample_temp15k, seqn_temp15l)

# remove anyone with NA to at least one of paq655 and pad660 when paq650 is 1
seqn_temp15m <- sample_temp15l$paq2007$seqn[sample_temp15l$paq2007$paq650==1 & (is.na(sample_temp15l$paq2007$paq655) | is.na(sample_temp15l$paq2007$pad660))]
sample_temp15m <- RemoveSeqn(sample_temp15l, seqn_temp15m)

# remove anyone with NA to paq665
seqn_temp15n <- sample_temp15m$paq2007$seqn[is.na(sample_temp15m$paq2007$paq665)]
sample_temp15n <- RemoveSeqn(sample_temp15m, seqn_temp15n)

# remove anyone with NA to at least one of paq670 and pad675 when paq665 is 1
seqn_temp15 <- sample_temp15n$paq2007$seqn[sample_temp15n$paq2007$paq665==1 & (is.na(sample_temp15n$paq2007$paq670) | is.na(sample_temp15n$paq2007$pad675))]
sample_temp15 <- RemoveSeqn(sample_temp15n, seqn_temp15)

# remove anyone with NA to day of recall 1
seqn_temp16a <- sample_temp15$nhanes$seqn[is.na(sample_temp15$nhanes$dr1day)]
sample_temp16a <- RemoveSeqn(sample_temp15, seqn_temp16a)

# remove anyone with NA to day of recall 2
seqn_temp16 <- sample_temp16a$nhanes$seqn[is.na(sample_temp16a$nhanes$dr2day)]
sample_temp16 <- RemoveSeqn(sample_temp16a, seqn_temp16)

# REMOVE ANYONE WITH MISSING SMOKING DATA
# remove anyone with NA to all smq020, smq040 and smq050q
seqn_temp17a <- sample_temp16$nhanes$seqn[rowSums(!is.na(sample_temp16$nhanes[,c("smq020", "smq040", "smq050q")]))==0]
sample_temp17a <- RemoveSeqn(sample_temp16, seqn_temp17a)

# remove anyone with NA to all smq040 and smq050q when smq020 is 1
seqn_temp17b <- sample_temp17a$nhanes$seqn[sample_temp17a$nhanes$smq020==1 & is.na(sample_temp17a$nhanes$smq040) & is.na(sample_temp17a$nhanes$smq050q)]
sample_temp17b <- RemoveSeqn(sample_temp17a, seqn_temp17b)

# remove anyone with NA to smq050q when smq020 is 1 and smq040 is 3
seqn_temp17 <- sample_temp17b$nhanes$seqn[sample_temp17b$nhanes$smq020==1 & sample_temp17b$nhanes$smq040==3 & is.na(sample_temp17b$nhanes$smq050q)]
sample_temp17 <- RemoveSeqn(sample_temp17b, seqn_temp17)

# remove anyone with NA to indfmpir
seqn_temp18 <- sample_temp17$nhanes$seqn[is.na(sample_temp17$nhanes$indfmpir)]
sample_temp18 <- RemoveSeqn(sample_temp17, seqn_temp18)

# remove anyone with NA to bmxwt
seqn_temp19a <- sample_temp18$nhanes$seqn[is.na(sample_temp18$nhanes$bmxwt)]
sample_temp19a <- RemoveSeqn(sample_temp18, seqn_temp19a)

# remove anyone with NA to bmxht
seqn_temp19 <- sample_temp19a$nhanes$seqn[is.na(sample_temp19a$nhanes$bmxht)]
sample_temp19 <- RemoveSeqn(sample_temp19a, seqn_temp19)

# remove anyone with NA to dr1ifdcd
seqn_temp20a <- sample_temp19$dr1iff$seqn[is.na(sample_temp19$dr1iff$dr1ifdcd)]
sample_temp20a <- RemoveSeqn(sample_temp19, seqn_temp20a)

# remove anyone with NA to dr2ifdcd
seqn_temp20 <- sample_temp20a$dr2iff$seqn[is.na(sample_temp20a$dr2iff$dr2ifdcd)]
sample_temp20 <- RemoveSeqn(sample_temp20a, seqn_temp20)

# remove anyone with NA to dmdmartl
seqn_temp21 <- sample_temp20$nhanes$seqn[is.na(sample_temp20$nhanes$dmdmartl)]
sample_temp21 <- RemoveSeqn(sample_temp20, seqn_temp21)

# remove anyone with NA to riagendr
seqn_temp22a <- sample_temp21$nhanes$seqn[is.na(sample_temp21$nhanes$riagendr)]
sample_temp22a <- RemoveSeqn(sample_temp21, seqn_temp22a)

# REMOVE ANYONE WHO IS FEMALE
sex <- 1
seqn_temp22 <- sample_temp22a$nhanes$seqn[sample_temp22a$nhanes$riagendr!=sex]
sample_temp22 <- RemoveSeqn(sample_temp22a, seqn_temp22)
finalmalesample <- sample_temp22

# CREATE COMPLETE DATA FOR 7,477 MEN
# produce a vector unique seqn for the remaining male participants
finalmaleseqn <- finalmalesample$nhanes$seqn

################################################################################
# D. Defines analysis variables based on Butler et al 2016
################################################################################

# DEFINE CONTINUOUS OUTCOME: waistcircumferenceY (WAIST CIRCUMFERENCE (cm))
waistcircumferenceY <- finalmalesample$nhanes$bmxwaist

# DEFINE 5-LEVEL CATEGORICAL EXPOSURE: alcoholX (1to2 (baseline), never, former, 3to4, atleast5)
alcoholX <- vector(length=length(finalmaleseqn))
# define as defined by Butler et al 2016
alcoholX[finalmalesample$nhanes$alq101==2 & finalmalesample$nhanes$alq110==2] <- "never"
alcoholX[finalmalesample$nhanes$alq110==1 & finalmalesample$nhanes$alq120q==0] <- "former"
alcoholX[(finalmalesample$nhanes$alq101==1 | finalmalesample$nhanes$alq110==1) & finalmalesample$nhanes$alq120q>=1 & (finalmalesample$nhanes$alq130==1 | finalmalesample$nhanes$alq130==2)] <- "1to2"
alcoholX[(finalmalesample$nhanes$alq101==1 | finalmalesample$nhanes$alq110==1) & finalmalesample$nhanes$alq120q>=1 & (finalmalesample$nhanes$alq130==3 | finalmalesample$nhanes$alq130==4)] <- "3to4"
alcoholX[(finalmalesample$nhanes$alq101==1 | finalmalesample$nhanes$alq110==1) & finalmalesample$nhanes$alq120q>=1 & (finalmalesample$nhanes$alq130>=5)] <- "atleast5"
alcoholXOld <- alcoholX
# create additional definitions
# assign “former” if alq101 is 1 and alq120q is 0
alcoholX[which(finalmalesample$nhanes$alq101==1 & finalmalesample$nhanes$alq120q==0)] <- "former"
# assign “1to2” if alq101 is 1 or alq110 is 1 AND alq102 is NA AND alq130 is 1 or 2
alcoholXOld <- alcoholX
alcoholX[(finalmalesample$nhanes$alq101==1 | finalmalesample$nhanes$alq110==1) & is.na(finalmalesample$nhanes$alq120q) & (finalmalesample$nhanes$alq130==1 | finalmalesample$nhanes$alq130==2)] <- "1to2"
# assign “3to4” if alq101 is 1 or alq110 is 1 AND alq102 is NA AND alq130 is 3 or 4
alcoholXOld <- alcoholX
alcoholX[(finalmalesample$nhanes$alq101==1 | finalmalesample$nhanes$alq110==1) & is.na(finalmalesample$nhanes$alq120q) & (finalmalesample$nhanes$alq130==3 | finalmalesample$nhanes$alq130==4)] <- "3to4"
# assign “atleast” if alq101 is 1 or alq110 is 1 AND alq102 is NA AND alq130 is 5 or more
alcoholXOld <- alcoholX
alcoholX[(finalmalesample$nhanes$alq101==1 | finalmalesample$nhanes$alq110==1) & is.na(finalmalesample$nhanes$alq120q) & finalmalesample$nhanes$alq130>=5] <- "atleast5"

# DEFINE MEASURED CONFOUNDERS:
# 3-LEVEL CATEGORICAL: ageC (20to39, 40to59, 60to79)
ageC <- vector(length=length(finalmaleseqn))
ageC[finalmalesample$nhanes$ridageyr>=20 & finalmalesample$nhanes$ridageyr<40] <- "20to39"
ageC[finalmalesample$nhanes$ridageyr>=40 & finalmalesample$nhanes$ridageyr<60] <- "40to59"
ageC[finalmalesample$nhanes$ridageyr>=60 & finalmalesample$nhanes$ridageyr<=79] <- "60to79"

# 5-LEVEL CATEGORICAL: raceC (MexicanAmerican, Hispanic, White, Black, Other)
raceC <- vector(length=length(finalmaleseqn))
raceC[finalmalesample$nhanes$ridreth1==1] <- "MexicanAmerican"
raceC[finalmalesample$nhanes$ridreth1==2] <- "Hispanic"
raceC[finalmalesample$nhanes$ridreth1==3] <- "White"
raceC[finalmalesample$nhanes$ridreth1==4] <- "Black"
raceC[finalmalesample$nhanes$ridreth1==5] <- "Other"

# 4-LEVEL CATEGORICAL: educationC (<HighSchool, HighSchoolGraduate, College, CollegeGraduate)
educationC <- vector(length=length(finalmaleseqn))
educationC[finalmalesample$nhanes$dmdeduc2==1 | finalmalesample$nhanes$dmdeduc2==2] <- "<HighSchool"
educationC[finalmalesample$nhanes$dmdeduc2==3] <- "HighSchoolGraduate"
educationC[finalmalesample$nhanes$dmdeduc2==4] <- "College"
educationC[finalmalesample$nhanes$dmdeduc2==5] <- "CollegeGraduate"

# 3-LEVEL CATEGORICAL: physicalactivitylevelC (low, moderate, high)
# start by setting every participant's physical activity level to low
physicalactivitylevelC <- rep("low", length(finalmaleseqn))
# create TotalInstancesPerWeek to represent the number of instances per week where the participant has done (moderate or vigorous) activity
TotalInstancesPerWeek <- rep(0, length(finalmaleseqn))
# create VigorousDaysPerWeekAbove20Minutes to represent the number of days per week where the participant has done vigorous activity for at least 20 minutes at each instance
VigorousDaysPerWeekAbove20Minutes <- rep(0, length(finalmaleseqn))
# create ModerateDaysPerWeekAbove30Minutes to represent the number of days per week where the participant has done (moderate or vigorous) activity for at least 30 minutes at each instance
ModerateDaysPerWeekAbove30Minutes <- rep(0, length(finalmaleseqn))
# create VigorousDaysPerWeek to represent the number of days per week where the participant has done vigorous activity
VigorousDaysPerWeek <- rep(0, length(finalmaleseqn))
# create TotalMetMinutesPerWeek to represent the MET minutes of physical activity the participant has done in a week
MetMinutesPerWeek <- rep(0, length(finalmaleseqn))
for (seqn in finalmaleseqn)
{
  # check which cycle the participant comes from
  if (finalmalesample$nhanes$cycle[finalmalesample$nhanes$seqn==seqn]=="2003to2004" | finalmalesample$nhanes$cycle[finalmalesample$nhanes$seqn==seqn]=="2005to2006")
  {
    # define physical activity level for participants from 2003-2004 and 2005-2006 cycles
    # check if seqn appears in paqiaf
    if (sum(finalmalesample$paqiaf$seqn==seqn)>0)
    {
      # one participant can have multiple sets of padtimes, paddurat and padmets across multiple rows
      IndividualPaqiaf <- finalmalesample$paqiaf[finalmalesample$paqiaf$seqn==seqn,]
      # process one row at a time for simplicity
      for (activity in 1:nrow(IndividualPaqiaf))
      {
        # adjust padtimes from instances from past 30 days to per week
        TotalInstancesPerWeek[finalmaleseqn==seqn] <- TotalInstancesPerWeek[finalmaleseqn==seqn] + IndividualPaqiaf$padtimes[activity]/30*7
        MetMinutesPerWeek[finalmaleseqn==seqn] <- MetMinutesPerWeek[finalmaleseqn==seqn] + prod(IndividualPaqiaf[activity,c("padtimes", "paddurat", "padmets")])/30*7
      }
    }
    
    # check if seqn appears in paq2003
    if (sum(finalmalesample$paq2003$seqn==seqn)>0)
    {
      # each participant has one set (row) of pad020, paq050q, paq050u, pad080, paq100, pad120, pad160
      # first deal with pad020, paq050q, paq050u and pad080
      # if pad020 is 1 then add to TotalInstancesPerWeek and MetMinutesPerWeek
      # this activity has a suggested MET score of 4 https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/paq_c.htm#Appendix_1._Suggested_MET_Scores
      if (finalmalesample$paq2003$pad020[finalmalesample$paq2003$seqn==seqn]==1)
      {
        # if paq050u is per day
        if (finalmalesample$paq2003$paq050u[finalmalesample$paq2003$seqn==seqn]==1)
        {
          TotalInstancesPerWeek[finalmaleseqn==seqn] <- TotalInstancesPerWeek[finalmaleseqn==seqn] + finalmalesample$paq2003$paq050q[finalmalesample$paq2003$seqn==seqn]*7
          MetMinutesPerWeek[finalmaleseqn==seqn] <- MetMinutesPerWeek[finalmaleseqn==seqn] + 7 * finalmalesample$paq2003$pad080[finalmalesample$paq2003$seqn==seqn] * 4
          # if paq050u is per week
        } else if (finalmalesample$paq2003$paq050u[finalmalesample$paq2003$seqn==seqn]==2)
        {
          TotalInstancesPerWeek[finalmaleseqn==seqn] <- TotalInstancesPerWeek[finalmaleseqn==seqn] + finalmalesample$paq2003$paq050q[finalmalesample$paq2003$seqn==seqn]
          # if paq050u is per week and paq050q is at least 7 instances (assuming equivalent to 7 days per week)
          if (finalmalesample$paq2003$paq050q[finalmalesample$paq2003$seqn==seqn]>=7)
          {
            MetMinutesPerWeek[finalmaleseqn==seqn] <- MetMinutesPerWeek[finalmaleseqn==seqn] + 7 * finalmalesample$paq2003$pad080[finalmalesample$paq2003$seqn==seqn] * 4
            # if paq050u is per week and paq050q is less 7 instances (assuming equivalent to paq050q days per week)
          } else {
            MetMinutesPerWeek[finalmaleseqn==seqn] <- MetMinutesPerWeek[finalmaleseqn==seqn] + finalmalesample$paq2003$paq050q[finalmalesample$paq2003$seqn==seqn] * finalmalesample$paq2003$pad080[finalmalesample$paq2003$seqn==seqn] * 4
          }
          # if paq050u is per month
        } else if (finalmalesample$paq2003$paq050u[finalmalesample$paq2003$seqn==seqn]==3)
        {
          TotalInstancesPerWeek[finalmaleseqn==seqn] <- TotalInstancesPerWeek[finalmaleseqn==seqn] + finalmalesample$paq2003$paq050q[finalmalesample$paq2003$seqn==seqn]/30*7
          # if paq050u is per month and paq050q is at least 30 instances (assuming equivalent to 7 days per week)
          if (finalmalesample$paq2003$paq050q[finalmalesample$paq2003$seqn==seqn]>=30)
          {
            MetMinutesPerWeek[finalmaleseqn==seqn] <- MetMinutesPerWeek[finalmaleseqn==seqn] + 7 * finalmalesample$paq2003$pad080[finalmalesample$paq2003$seqn==seqn] * 4
            # if paq050u is per month and paq050q is less than 30 instances (assuming equivalent to paq050q/30*7 days per week)
          } else {
            MetMinutesPerWeek[finalmaleseqn==seqn] <- MetMinutesPerWeek[finalmaleseqn==seqn] + finalmalesample$paq2003$paq050q[finalmalesample$paq2003$seqn==seqn]/30*7 * finalmalesample$paq2003$pad080[finalmalesample$paq2003$seqn==seqn] * 4
          }
        }
      }
      # second deal with paq100, pad120 and pad160
      # if paq100 is 1 then add to TotalInstancesPerWeek and MetMinutesPerWeek
      # this activity has a suggested MET score of 4.5 https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/PAQ_C.htm#Appendix_1._Suggested_MET_Scores
      if (finalmalesample$paq2003$paq100[finalmalesample$paq2003$seqn==seqn]==1)
      {
        TotalInstancesPerWeek[finalmaleseqn==seqn] <- TotalInstancesPerWeek[finalmaleseqn==seqn] + finalmalesample$paq2003$pad120[finalmalesample$paq2003$seqn==seqn]/30*7
        MetMinutesPerWeek[finalmaleseqn==seqn] <- MetMinutesPerWeek[finalmaleseqn==seqn] + prod(finalmalesample$paq2003[finalmalesample$paq2003$seqn==seqn,c("pad120", "pad160")])/30*7 * 4.5
      }
    }
    # define physical activity level for participants from 2007-2008, 2009-2010 and 2011-2012 cycles
  } else {
    # each participant has one set (row) of paq605, paq610, pad615, paq620, paq625, pad630, paq635, paq640, pad645, paq650, paq655, pad660, paq665, paq670 and pad675
    # first deal with paq605, paq610 and pad615
    # if paq605 is 1 then add to VigorousDaysPerWeek, VigorousDaysPerWeekAbove20Minutes and MetMinutesPerWeek
    # this activity has a suggested MET score of 8 https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2007/DataFiles/paq_e.htm#Appendix_1.__Suggested_MET_Scores
    if (finalmalesample$paq2007$paq605[finalmalesample$paq2007$seqn==seqn]==1)
    {
      VigorousDaysPerWeek[finalmaleseqn==seqn] <- VigorousDaysPerWeek[finalmaleseqn==seqn] + finalmalesample$paq2007$paq610[finalmalesample$paq2007$seqn==seqn]
      # check if duration was at least 20 minutes at each day
      if (finalmalesample$paq2007$pad615[finalmalesample$paq2007$seqn==seqn]>=20)
      {
        VigorousDaysPerWeekAbove20Minutes[finalmaleseqn==seqn] <- VigorousDaysPerWeekAbove20Minutes[finalmaleseqn==seqn] + finalmalesample$paq2007$paq610[finalmalesample$paq2007$seqn==seqn]
      }
      MetMinutesPerWeek[finalmaleseqn==seqn] <- MetMinutesPerWeek[finalmaleseqn==seqn] + prod(finalmalesample$paq2007[finalmalesample$paq2007$seqn==seqn,c("paq610", "pad615")]) * 8
    }
    
    # second deal with paq620, paq625 and pad630
    # if paq620 is 1 then add to ModerateDaysPerWeekAbove30Minutes and MetMinutesPerWeek
    # this activity has a suggested MET score of 4 https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2007/DataFiles/paq_e.htm#Appendix_1.__Suggested_MET_Scores
    if (finalmalesample$paq2007$paq620[finalmalesample$paq2007$seqn==seqn]==1)
    {
      # check if duration was at least 30 minutes at each day
      if (finalmalesample$paq2007$pad630[finalmalesample$paq2007$seqn==seqn]>=30)
      {
        ModerateDaysPerWeekAbove30Minutes[finalmaleseqn==seqn] <- ModerateDaysPerWeekAbove30Minutes[finalmaleseqn==seqn] + finalmalesample$paq2007$paq625[finalmalesample$paq2007$seqn==seqn]
      }
      MetMinutesPerWeek[finalmaleseqn==seqn] <- MetMinutesPerWeek[finalmaleseqn==seqn] + prod(finalmalesample$paq2007[finalmalesample$paq2007$seqn==seqn,c("paq625", "pad630")]) * 4
    }
    
    # third deal with paq635, paq640 and pad645
    # if paq635 is 1 then add to ModerateDaysPerWeekAbove30Minutes and MetMinutesPerWeek
    # this activity has a suggested MET score of 4 https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2007/DataFiles/paq_e.htm#Appendix_1.__Suggested_MET_Scores
    if (finalmalesample$paq2007$paq635[finalmalesample$paq2007$seqn==seqn]==1)
    {
      # check if duration was at least 30 minutes at each day
      if (finalmalesample$paq2007$pad645[finalmalesample$paq2007$seqn==seqn]>=30)
      {
        # add to ModerateDaysPerWeekAbove30Minutes if TRUE
        ModerateDaysPerWeekAbove30Minutes[finalmaleseqn==seqn] <- ModerateDaysPerWeekAbove30Minutes[finalmaleseqn==seqn] + finalmalesample$paq2007$paq640[finalmalesample$paq2007$seqn==seqn]
      }
      # add to MET minutes from this activity TotalMetMinutesPerWeek
      MetMinutesPerWeek[finalmaleseqn==seqn] <- MetMinutesPerWeek[finalmaleseqn==seqn] + prod(finalmalesample$paq2007[finalmalesample$paq2007$seqn==seqn,c("paq640", "pad645")]) * 4
    }
    
    # fourth deal with paq650, paq655 and pad660
    # if paq635 is 1 then add to VigorousDaysPerWeek, VigorousDaysPerWeekAbove20Minutes and MetMinutesPerWeek
    # this activity has a suggested MET score of 8 https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2007/DataFiles/paq_e.htm#Appendix_1.__Suggested_MET_Scores
    if (finalmalesample$paq2007$paq650[finalmalesample$paq2007$seqn==seqn]==1)
    {
      VigorousDaysPerWeek[finalmaleseqn==seqn] <- VigorousDaysPerWeek[finalmaleseqn==seqn] + finalmalesample$paq2007$paq655[finalmalesample$paq2007$seqn==seqn]
      # check if duration was at least 20 minutes at each day
      if (finalmalesample$paq2007$pad660[finalmalesample$paq2007$seqn==seqn]>=20)
      {
        VigorousDaysPerWeekAbove20Minutes[finalmaleseqn==seqn] <- VigorousDaysPerWeekAbove20Minutes[finalmaleseqn==seqn] + finalmalesample$paq2007$paq655[finalmalesample$paq2007$seqn==seqn]
      }
      MetMinutesPerWeek[finalmaleseqn==seqn] <- MetMinutesPerWeek[finalmaleseqn==seqn] + prod(finalmalesample$paq2007[finalmalesample$paq2007$seqn==seqn,c("paq655", "pad660")]) * 8
    }
    
    # finally deal with paq665, paq670 and pad675
    # if paq665 is 1 then add to ModerateDaysPerWeekAbove30Minutes and MetMinutesPerWeek
    # this activity has a suggested MET score of 4 https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2007/DataFiles/paq_e.htm#Appendix_1.__Suggested_MET_Scores
    if (finalmalesample$paq2007$paq665[finalmalesample$paq2007$seqn==seqn]==1)
    {
      # check if duration was at least 30 minutes at each day
      if (finalmalesample$paq2007$pad675[finalmalesample$paq2007$seqn==seqn]>=30)
      {
        ModerateDaysPerWeekAbove30Minutes[finalmaleseqn==seqn] <- ModerateDaysPerWeekAbove30Minutes[finalmaleseqn==seqn] + finalmalesample$paq2007$paq670[finalmalesample$paq2007$seqn==seqn]
      }
      MetMinutesPerWeek[finalmaleseqn==seqn] <- MetMinutesPerWeek[finalmaleseqn==seqn] + prod(finalmalesample$paq2007[finalmalesample$paq2007$seqn==seqn,c("paq670", "pad675")]) * 4
    }
  }
  
  # now define physicalactivitylevelC using TotalInstancesPerWeek, VigorousDaysPerWeekAbove20Minutes, ModerateDaysPerWeekAbove30Minutes, VigorousDaysPerWeek and MetMinutesPerWeek
  # moderate: >=3 days of vigorous PA for >=20 minutes
  if (VigorousDaysPerWeekAbove20Minutes[finalmaleseqn==seqn]>=3)
  {
    physicalactivitylevelC[finalmaleseqn==seqn] <- "moderate"
  }
  # moderate: >=5 days of self-reported moderate PA [eg, biking and walking] for >=30 minutes
  if (ModerateDaysPerWeekAbove30Minutes[finalmaleseqn==seqn]>=5)
  {
    physicalactivitylevelC[finalmaleseqn==seqn] <- "moderate"
  }
  # moderate: >=5 instances/wk of total PA and >=600 MET-min/wk of PA
  if (TotalInstancesPerWeek[finalmaleseqn==seqn]>=5 & MetMinutesPerWeek[finalmaleseqn==seqn]>=600)
  {
    physicalactivitylevelC[finalmaleseqn==seqn] <- "moderate"
  }
  # high: >=3 days of vigorous PA and >=1,500 MET-min/wk of PA
  if (VigorousDaysPerWeek[finalmaleseqn==seqn]>=3 & MetMinutesPerWeek[finalmaleseqn==seqn]>=1500)
  {
    physicalactivitylevelC[finalmaleseqn==seqn] <- "high"
  }
  # high: >=7 instances/wk of total PA and >=3,000 MET-min/wk of PA
  if (TotalInstancesPerWeek[finalmaleseqn==seqn]>=7 & MetMinutesPerWeek[finalmaleseqn==seqn]>=3000)
  {
    physicalactivitylevelC[finalmaleseqn==seqn] <- "high"
  }
}

# 5-LEVEL CATEGORICAL: yearC (2003to2004, 2005to2006, 2007to2008, 2009to2010, 2011to2012)
yearC <- finalmalesample$nhanes$cycle

# BINARY: recall1C (weekend, weekday)
recall1C <- vector(length=length(finalmaleseqn))
recall1C[finalmalesample$nhanes$dr1day==1 | finalmalesample$nhanes$dr1day==6 | finalmalesample$nhanes$dr1day==7] <- "weekend"
recall1C[finalmalesample$nhanes$dr1day==2 | finalmalesample$nhanes$dr1day==3 | finalmalesample$nhanes$dr1day==4 | finalmalesample$nhanes$dr1day==5] <- "weekday"

# BINARY: recall2C (weekend, weekday)
recall2C <- vector(length=length(finalmaleseqn))
recall2C[finalmalesample$nhanes$dr2day==1 | finalmalesample$nhanes$dr2day==6 | finalmalesample$nhanes$dr2day==7] <- "weekend"
recall2C[finalmalesample$nhanes$dr2day==2 | finalmalesample$nhanes$dr2day==3 | finalmalesample$nhanes$dr2day==4 | finalmalesample$nhanes$dr2day==5] <- "weekday"

# 3-LEVEL CATEGORICAL: smokingC (never, current, former)
smokingC <- vector(length=length(finalmaleseqn))
# logically if one responded 2 (no) to smq020 (Have you smoked at least 100 cigarettes in your entire life?), then their response to smq040 (Do you now smoke cigarettes) should be 3 (not at all) but some participants who responded 2 to smq020 did not respond to smq040.
# so only use smq020 to define "never" smoking status
smokingC[which(finalmalesample$nhanes$smq020==2)] <- "never"
smokingC[which(finalmalesample$nhanes$smq040==1 | finalmalesample$nhanes$smq040==2)] <- "current"
smokingC[which(finalmalesample$nhanes$smq040==3 & finalmalesample$nhanes$smq050q>=1)] <- "former"

# 3-LEVEL CATEGORICAL: incomeC (0to130, 131to299, atleast300)
incomeC <- vector(length=length(finalmaleseqn))
incomeC[finalmalesample$nhanes$indfmpir<=1.3] <- "0to130"
incomeC[finalmalesample$nhanes$indfmpir>1.3 & finalmalesample$nhanes$indfmpir<3] <- "131to299"
incomeC[finalmalesample$nhanes$indfmpir>=3] <- "atleast300"

# 3-LEVEL CATEGORICAL: misreporingC (accurate, under, over)
# bmr in kcal using the goldberg revised method as explained in Mendez et al 2011 and based on the equation proposed by Mifflin et al 1990
BMR <- finalmalesample$nhanes$bmxwt*9.99+finalmalesample$nhanes$bmxht*6.25-finalmalesample$nhanes$ridageyr*4.92+5
# PAL is the mean PAL for the population under study as defined in Black 2000
PAL <- c(1.4, 1.55, 1.75)[1*(physicalactivitylevelC=="low")+2*(physicalactivitylevelC=="moderate")+3*(physicalactivitylevelC=="high")] # set by Butler et al 2017
# sd for the width of confidence limit for energy intake (EI_rep:BMR) (as defined in Black 2000)
SD <- 1.5 # set by Butler et al 2017
# within-subject coefficient of variation in energy intake as defined in Black 2000
CV_wEI <- 23 # extracted from Black 2000 in Clinical nutrition
# number of days of diet assessment
d <- 2
# coefficient of variation of repeated BMR measurements or the precision of estimated compared with measured BMR (as defined in Black 2000)
CV_wB <- 8.5 # extracted from Black 2000 in Clinical nutrition
# total variation in PAL (as defined in Black 2000)
CV_tP <- 15 # extracted from Black 2000 in Clinical nutrition
S <- sqrt(CV_wEI^2/d+CV_wB^2+CV_tP^2)
# number of subjects under study set to 1 for individual level (as defined in Black 2000 and https://www.efsa.europa.eu/sites/default/files/efsa_rep/blobserver_assets/3944A-8-2-1.pdf)
n <- 1
# the limits from Black 2000
UpperLimit <- PAL*exp(SD*S/100/sqrt(n))
LowerLimit <- PAL*exp(-SD*S/100/sqrt(n))
# energy intake (EI) measured in kcal is total daily energy intake from food and beverages
EI <- vector(length=length(finalmaleseqn))
for (seqn in finalmaleseqn)
{
  EI[finalmaleseqn==seqn] <- mean(c(finalmalesample$nhanes$dr1tkcal[finalmalesample$nhanes$seqn==seqn], finalmalesample$nhanes$dr2tkcal[finalmalesample$nhanes$seqn==seqn]))
}
# rEI:BMR = ratio between EI and BMR as defined in Black 2000
EI2BMR <- EI/BMR
misreportingC <- vector(length=length(finalmaleseqn))
misreportingC[LowerLimit<=EI2BMR & EI2BMR<=UpperLimit] <- "accurate"
misreportingC[EI2BMR<LowerLimit] <- "under"
misreportingC[EI2BMR>UpperLimit] <- "over"

# DEFINE UNMEASURED CONFOUNDERS:
# CONTINUOUS: energyU (daily non-alcoholic energy intake (kcal))
# define Alcoholic for each food item
drxfcd_c_temp3 <- AlcoholLabel(drxfcd_c_temp2)
drxfcd_d_temp3 <- AlcoholLabel(drxfcd_d_temp2)
drxfcd_e_temp3 <- AlcoholLabel(drxfcd_e_temp2)
drxfcd_f_temp3 <- AlcoholLabel(drxfcd_f_temp2)
drxfcd_g_temp3 <- AlcoholLabel(drxfcd_g_temp2)

# set drxfcd_c_temp3 to drxfcd_c
drxfcd_c <- drxfcd_c_temp3
# set drxfcd_d_temp3 to drxfcd_d
drxfcd_d <- drxfcd_d_temp3
# set drxfcd_e_temp3 to drxfcd_e
drxfcd_e <- drxfcd_e_temp3
# set drxfcd_f_temp3 to drxfcd_f
drxfcd_f <- drxfcd_f_temp3
# set drxfcd_g_temp3 to drxfcd_g
drxfcd_g <- drxfcd_g_temp3

# calculate average total alcoholic beverage energy intake (kcal)
alcoholicEI <- vector(length=length(finalmaleseqn))
for (seqn in finalmaleseqn)
{
  # get the individual EI data for seqn
  dr1iff4seqn <- finalmalesample$dr1iff[is.element(finalmalesample$dr1iff$seqn, seqn),]
  # set daily alcoholic beverage energy to 0
  TotalDay1AlcholicBeverageIntake4seqn <- 0
  # for each food item the participant consumed,
  for (item in 1:nrow(dr1iff4seqn))
  {
    # look up whether the item is alcoholic or not
    if (finalmalesample$nhanes[finalmalesample$nhanes$seqn==seqn,"cycle"]=="2003to2004")
    {
      AlcoholicStatus4seqn4item <- drxfcd_c[drxfcd_c$drxfdcd==dr1iff4seqn[item,"dr1ifdcd"],"Alcoholic"]
    } else if (finalmalesample$nhanes[finalmalesample$nhanes$seqn==seqn,"cycle"]=="2005to2006")
    {
      AlcoholicStatus4seqn4item <- drxfcd_d[drxfcd_d$drxfdcd==dr1iff4seqn[item,"dr1ifdcd"],"Alcoholic"]
    } else if (finalmalesample$nhanes[finalmalesample$nhanes$seqn==seqn,"cycle"]=="2007to2008")
    {
      AlcoholicStatus4seqn4item <- drxfcd_e[drxfcd_e$drxfdcd==dr1iff4seqn[item,"dr1ifdcd"],"Alcoholic"]
    } else if (finalmalesample$nhanes[finalmalesample$nhanes$seqn==seqn,"cycle"]=="2009to2010")
    {
      AlcoholicStatus4seqn4item <- drxfcd_f[drxfcd_f$drxfdcd==dr1iff4seqn[item,"dr1ifdcd"],"Alcoholic"]
    } else if (finalmalesample$nhanes[finalmalesample$nhanes$seqn==seqn,"cycle"]=="2011to2012")
    {
      AlcoholicStatus4seqn4item <- drxfcd_g[drxfcd_g$drxfdcd==dr1iff4seqn[item,"dr1ifdcd"],"Alcoholic"]
    }
    # if the item is alcoholic,
    if (AlcoholicStatus4seqn4item=="alcoholic")
    {
      # add the item's calorie to daily alcoholic beverage energy
      TotalDay1AlcholicBeverageIntake4seqn <- TotalDay1AlcholicBeverageIntake4seqn+dr1iff4seqn[item,"dr1ikcal"]
    }
  }
  # repeat for day 2
  dr2iff4seqn <- finalmalesample$dr2iff[is.element(finalmalesample$dr2iff$seqn, seqn),]
  TotalDay2AlcholicBeverageIntake4seqn <- 0
  for (item in 1:nrow(dr2iff4seqn))
  {
    if (finalmalesample$nhanes[finalmalesample$nhanes$seqn==seqn,"cycle"]=="2003to2004")
    {
      AlcoholicStatus4seqn4item <- drxfcd_c[drxfcd_c$drxfdcd==dr2iff4seqn[item,"dr2ifdcd"],"Alcoholic"]
    } else if (finalmalesample$nhanes[finalmalesample$nhanes$seqn==seqn,"cycle"]=="2005to2006")
    {
      AlcoholicStatus4seqn4item <- drxfcd_d[drxfcd_d$drxfdcd==dr2iff4seqn[item,"dr2ifdcd"],"Alcoholic"]
    } else if (finalmalesample$nhanes[finalmalesample$nhanes$seqn==seqn,"cycle"]=="2007to2008")
    {
      AlcoholicStatus4seqn4item <- drxfcd_e[drxfcd_e$drxfdcd==dr2iff4seqn[item,"dr2ifdcd"],"Alcoholic"]
    } else if (finalmalesample$nhanes[finalmalesample$nhanes$seqn==seqn,"cycle"]=="2009to2010")
    {
      AlcoholicStatus4seqn4item <- drxfcd_f[drxfcd_f$drxfdcd==dr2iff4seqn[item,"dr2ifdcd"],"Alcoholic"]
    } else if (finalmalesample$nhanes[finalmalesample$nhanes$seqn==seqn,"cycle"]=="2011to2012")
    {
      AlcoholicStatus4seqn4item <- drxfcd_g[drxfcd_g$drxfdcd==dr2iff4seqn[item,"dr2ifdcd"],"Alcoholic"]
    }
    if (AlcoholicStatus4seqn4item=="alcoholic")
    {
      TotalDay2AlcholicBeverageIntake4seqn <- TotalDay2AlcholicBeverageIntake4seqn+dr2iff4seqn[item,"dr2ikcal"]
    }
  }
  # define average total alcoholic beverage energy intake as the mean of the total alcoholic beverage energy intakes from days 1 and 2 
  alcoholicEI[finalmaleseqn==seqn] <- mean(c(TotalDay1AlcholicBeverageIntake4seqn, TotalDay2AlcholicBeverageIntake4seqn))
}
# define daily non-alcoholic energy intake (energyU) as daily energy intake minus daily alcoholic energy intake
energyU <- EI-alcoholicEI

# BINARY: maritalstatusU (0 currently in a relationship, 1 currently single)
maritalstatusU <- vector(length=length(finalmaleseqn))
maritalstatusU[finalmalesample$nhanes$dmdmartl==2 | finalmalesample$nhanes$dmdmartl==3 | finalmalesample$nhanes$dmdmartl==4 | finalmalesample$nhanes$dmdmartl==5] <- 1
maritalstatusU[finalmalesample$nhanes$dmdmartl==1 | finalmalesample$nhanes$dmdmartl==6] <- 0

# create a data frame by combining all analysis variables 
Data <- data.frame(waistcircumferenceY, alcoholX, ageC, raceC, educationC, physicalactivitylevelC, yearC, recall1C, recall2C, smokingC, incomeC, misreportingC, energyU, maritalstatusU)

# WRITE THE DATA TO A csv FILE
write.csv(Data, file="Data/Data.csv", row.names=F)
