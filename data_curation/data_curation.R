# this script gathers and prepares all data for the analyses

# load required packages
library(readxl)
library(tidyverse)
library(naniar)
library(lme4)
library(kit)

hints <- function(df, col, tops = 5, tp = 1){
  # a brief overview intended to give hints at potential issues in the data
  # tops can be used to determine the number of max/min values to print
  # if tp is anything but one, no separate results fro each tp will be produced 
  if(tp == 1){ # additionally separate for each timepoint
    vars = c(df[,col], df[df$tp == "bl", col], df[df$tp == "fu", col]) 
  } else{ # only for bl and fu together
    vars = c(df$col)
  } 
  intros = list("Results for all rows", "Results for all baseline rows", "Results for all follow-up rows")
  histogram_titles = list("all data", "baseline data", "follow-up data")
  n = 1
  for (var in vars) {
    print(intros[n])
    print(paste0("n of observations: ", length(which(!is.na(var))), "; n of NAs: ", length(which(is.na(var))), "; % of NAs: ", (length(which(is.na(var)))/length(var))*100))
    print("20 values that occur in the column: ")
    print(head(round(unique(var), digits = 3), n = 20))
    hist(var, main = histogram_titles[n])
    if(length(unique(var)) > 4){
      print(paste0("top ", tops, " max: "))
      print(topn(var, tops, decreasing = T, index = F))
      print(paste0("rown of max: ", which.max(var)))
      print(paste0("top ", tops, " min: "))
      print(topn(var, tops, decreasing = F, index = F))
      print(paste0("rown of min: ", which.min(var)))
      print(paste0("median: ", median(var, na.rm = T)))
      print(paste0("mean: ", mean(var, na.rm = T)))
      print(paste0("sd: ", sd(var, na.rm = T)))
    }
    else{ # for dichotomous variables
      print("n of measures == 1")
      print(length(which(var == 1)))
      print("n of measures == 0")
      print(length(which(var == 0)))
      print("percentage of measures == 1")
      print(paste0(length(which(var == 1))/length(which(!is.na(var))), " %"))
    }
    n = n + 1
  }
}

# firstly, we read in and curate all data provided by LIFE
# most data are in one large file for baseline and follow-up
bl_datajoin <- read_excel("/data/pt_life/ResearchProjects/LLammer/gamms/Data/PV772_Lammer/Adult Basis/data/PV0772_datajoin.xlsx")
fu_datajoin <- read_excel("/data/pt_life/ResearchProjects/LLammer/gamms/Data/PV772_Lammer/Adult Follow-up/data/PV0772_datajoin_ohne_Anamnesen.xlsx")
# read in all MRI data
radiological_assessment <- read_csv(file = "/data/gh_gr_agingandobesity_share/life_shared/Data/Preprocessed/derivatives/radiological_assessment/radio_assessment_long.csv")
colnames(radiological_assessment)[2] <- "mpi_id"
colnames(radiological_assessment)[20] <- "tp"
radiological_assessment$Messung.Datum_TI <- ifelse(is.na(radiological_assessment$Messung.Datum_TI), NA, paste0(substr(radiological_assessment$Messung.Datum_TI, 7, 10), "-", substr(radiological_assessment$Messung.Datum_TI, 4, 5), "-", substr(radiological_assessment$Messung.Datum_TI, 1, 2)))
samseg_data <- read_csv("/data/pt_life_freesurfer/samseg/results_summary.csv", col_select = 2:51)
colnames(samseg_data)[1] <- "mrt_pseudonym"
colnames(samseg_data)[2] <- "tp"
colnames(samseg_data)[50] <- "LESIONS"
samseg_data$tp <- ifelse(samseg_data$tp == "BL", "bl", "fu") # turn into lowercase for compatibility
LCL_QA <- read_csv("/data/pt_life_whm/Results/Tables/Datatable_for_LIFE/qa_info_LIFE_WMH_LST.csv")
colnames(LCL_QA)[1] <- "mrt_pseudonym"
################# enter ventricle expansion data here  #################
# load the ID-translation table to connect the MRI data with the other data
translator <- read_excel("/data/pt_life/ResearchProjects/LLammer/gamms/Data/PV772_Lammer/PV772_PV-MRT-Pseudonymliste.xlsx")
colnames(translator)[1] <- "ID"
samseg_data <- merge(samseg_data, translator, by = "mrt_pseudonym")
radiological_assessment <- merge(radiological_assessment, translator, by = "mrt_pseudonym")
LCL_QA <- merge(LCL_QA, translator, by = "mrt_pseudonym")
# some assessments are duplicated, keep only one, keep the one without NA in the evaluation if two exist
duplicated <- radiological_assessment %>%
  group_by(ID, tp) %>%
  filter(n() > 1) %>%
  filter(if_any(med_Befund.Bewertung_Bef, ~ !is.na(.)) | n() == 1) %>%
  slice_head(n = 1) %>%
  ungroup()
unduplicated <- radiological_assessment %>%
  group_by(ID, tp) %>%
  filter(n() == 1) %>%
  ungroup()
radiological_assessment <- rbind(duplicated, unduplicated)

# our first larger task is to handle duplicates
# the FU cardiological anamnesis is provided in a separate file 
# read in the FU cardiological anamnesis 
fu_card_anam <- read_excel("/data/pt_life/ResearchProjects/LLammer/gamms/Data/PV772_Lammer/Adult Follow-up/data/PV0772_T01226_NODUP.xlsx")
# for some participants the anamnesis was performed twice (in person and as a paper follow-up)
# if both are available, the in person anamnesis will be preferred
df_fu_card_anam_2 <- fu_card_anam %>%
  group_by(SIC) %>%
  filter(n() == 2) %>% # select all participants with 2 anamneses
  ungroup() %>%
  filter(GRUPPE == "A1_FU1") # choose the row with the in-person anamnesis
df_fu_card_anam_1 <- fu_card_anam %>%
  group_by(SIC) %>%
  filter(n() == 1) %>%
  ungroup() 
fu_card_anam <- rbind(df_fu_card_anam_1, df_fu_card_anam_2) # reconnect particpants with 1 and 2 anamneses into one df
fu_card_anam <- fu_card_anam %>% # rename cols to merge card anamnesis with other follow-up data
  rename(
    TEILNEHMER_SIC = SIC,
    SGROUP_card_anamn = GRUPPE
  )
fu_datajoin <- merge(fu_datajoin, fu_card_anam, by = "TEILNEHMER_SIC", all.y = T) # merge the dfs and keep all participants in the datajoin with no card anamnesis
# the same has to be reiterated for the general medical anamnesis of which some participants also have two versions
fu_med_anam <- read_excel("/data/pt_life/ResearchProjects/LLammer/gamms/Data/PV772_Lammer/Adult Follow-up/data/PV0772_T01228_NODUP.xlsx")
fu_med_anam_2 <- fu_med_anam %>%
  group_by(SIC) %>%
  filter(n() == 2) %>% # select all participants with 2 anamneses
  ungroup() %>%
  filter(GRUPPE == "A1_FU1") # choose the row with the in person anamnesis
fu_med_anam_1 <- fu_med_anam %>%
  group_by(SIC) %>%
  filter(n() == 1) %>%
  ungroup() 
fu_med_anam <- rbind(fu_med_anam_1, fu_med_anam_2) # reconnect particpants with 1 and 2 anamneses into one df
fu_med_anam <- fu_med_anam %>% # rename cols to merge card anamnesis with other follow-up data
  rename(
    TEILNEHMER_SIC = SIC,
    SGROUP_med_anamn = GRUPPE
  )
fu_datajoin <- merge(fu_datajoin, fu_med_anam, by = "TEILNEHMER_SIC", all.y = T) # merge the dfs and keep all participants in the datajoin with no card anamnesis
# and it has to be reiterated for the CES-D depression questionnaire
fu_cesd <- read_excel("/data/pt_life/ResearchProjects/LLammer/gamms/Data/PV772_Lammer/Adult Follow-up/data/PV0772_T00013_NODUP.xlsx")
fu_cesd_2 <- fu_cesd %>%
  group_by(SIC) %>%
  filter(n() == 2) %>% # select all participants with 2 anamneses
  ungroup() %>%
  filter(GRUPPE == "A1_FU1") # choose the row with the in person anamnesis
fu_cesd_1 <- fu_cesd %>%
  group_by(SIC) %>%
  filter(n() == 1) %>%
  ungroup() 
fu_cesd <- rbind(fu_cesd_1, fu_cesd_2) # reconnect particpants with 1 and 2 anamneses into one df
fu_cesd <- fu_cesd %>% # rename cols to merge card anamnesis with other follow-up data
  rename(
    TEILNEHMER_SIC = SIC,
    SGROUP_med_anamn = GRUPPE
  )
fu_datajoin <- merge(fu_datajoin, fu_cesd, by = "TEILNEHMER_SIC", all.y = T) # merge the dfs and keep all participants in the datajoin with no card anamnesis

# there are also duplicated rows in the bl datajoin
# as in previous work, we will only keep one row per participant preferring the main testing phase over pilots
bl_datajoin <- arrange(.data = bl_datajoin, "TEILNEHMER_SIC", factor(SGROUP, levels = c("A1_HAUPT01", "A1_PILOT_2", "A1_PILOT", "A1_FEASI001", "A1_TEST")))
bl_datajoin <- distinct(.data = bl_datajoin, TEILNEHMER_SIC, .keep_all = T)
# the fu datajoin contains even more duplicates because most participants that took part in the in-person FU also filled in the questionnaires for the papaer follow-up
# as some questionnaires were only part of the paper follow-up we can not simply choose the observations from the in person follow-up
# firstly, we turn the df into the wide format so that paper and in-person fu are in one row
fu_datajoin_wide <- reshape(fu_datajoin, timevar = "SGROUP", idvar = "TEILNEHMER_SIC", direction = "wide", sep = "_") # transform table into wide format
# we know that LSNS was either answered in person or in the paper FU
# we will thus fill in empty values in the in person columns with those from the paper FU
for (col in grep("^LSNS.*A1_FU1$", names(fu_datajoin_wide), value = TRUE)) {
  # Find the corresponding column for each "LSNS...A1_FU1" column
  ans_column <- gsub("A1_FU1$", "A1_FU1_ANSCH", col)
  # Replace NAs in the current column with values from the corresponding "_ANSCH" column
  fu_datajoin_wide[[col]] <- ifelse(is.na(fu_datajoin_wide[[col]]), fu_datajoin_wide[[ans_column]], fu_datajoin_wide[[col]])
}
fu_datajoin_wide$LSNS_DATUM_A1_FU1 <- as.Date(as.POSIXct(fu_datajoin_wide$LSNS_DATUM_A1_FU1)) # undo transformation
# reiterate the same for the CESD questionnaire
for (col in grep("^CES_D.*A1_FU1$", names(fu_datajoin_wide), value = TRUE)) {
  # Find the corresponding column 
  ans_column <- gsub("A1_FU1$", "A1_FU1_ANSCH", col)
  # Replace NAs in the current column with values from the corresponding "_ANSCH" column
  fu_datajoin_wide[[col]] <- ifelse(is.na(fu_datajoin_wide[[col]]), fu_datajoin_wide[[ans_column]], fu_datajoin_wide[[col]])
}
fu_datajoin_wide$CES_D_DATUM_A1_FU1 <- as.Date(as.POSIXct(fu_datajoin_wide$CES_D_DATUM_A1_FU1))
# reiterate the same for the GAD7 questionnaire
for (col in grep("^GAD7.*A1_FU1$", names(fu_datajoin_wide), value = TRUE)) {
  # Find the corresponding column 
  ans_column <- gsub("A1_FU1$", "A1_FU1_ANSCH", col)
  # Replace NAs in the current column with values from the corresponding "_ANSCH" column
  fu_datajoin_wide[[col]] <- ifelse(is.na(fu_datajoin_wide[[col]]), fu_datajoin_wide[[ans_column]], fu_datajoin_wide[[col]])
}
fu_datajoin_wide$GAD7_DATUM_A1_FU1 <- as.Date(as.POSIXct(fu_datajoin_wide$GAD7_DATUM_A1_FU1))
fu_datajoin_wide$GAD7_STARTZEIT_A1_FU1 <- as.Date(as.POSIXct(fu_datajoin_wide$GAD7_STARTZEIT_A1_FU1))
# reiterate the same for gender
for (col in grep("^TEILNEHMER_GESCHLECHT.*A1_FU1$", names(fu_datajoin_wide), value = TRUE)) {
  # Find the corresponding column 
  ans_column <- gsub("A1_FU1$", "A1_FU1_ANSCH", col)
  # Replace NAs in the current column with values from the corresponding "_ANSCH" column
  fu_datajoin_wide[[col]] <- ifelse(is.na(fu_datajoin_wide[[col]]), fu_datajoin_wide[[ans_column]], fu_datajoin_wide[[col]])
}
# and the month of birth
# reiterate the same for gender
for (col in grep("^TEILNEHMER_GEB_JJJJMM.*A1_FU1$", names(fu_datajoin_wide), value = TRUE)) {
  # Find the corresponding column 
  ans_column <- gsub("A1_FU1$", "A1_FU1_ANSCH", col)
  # Replace NAs in the current column with values from the corresponding "_ANSCH" column
  fu_datajoin_wide[[col]] <- ifelse(is.na(fu_datajoin_wide[[col]]), fu_datajoin_wide[[ans_column]], fu_datajoin_wide[[col]])
}
# for sociodemographic data the cols from the paper fu start with "SOZDEM" while those from the in person fu start with "SOZIO"
sozdem_columns <- grep("^SOZDEM", colnames(fu_datajoin_wide), value = TRUE)
sozio_columns <- c(grep("^SOZIO", colnames(fu_datajoin_wide), value = TRUE), rep("NA",15))
# before combining paper and in-person FU we have to homogenise one coding inconsistency
fu_datajoin_wide$SOZDEM_ERWERB_A1_FU1_ANSCH <- ifelse(fu_datajoin_wide$SOZDEM_ERWERB_A1_FU1_ANSCH == 4, 0, fu_datajoin_wide$SOZDEM_ERWERB_A1_FU1_ANSCH)
# all relevant columns have been asked in paper and in person
extract_suffix <- function(col_name, prefix) {
  sub(prefix, "", col_name)
}
sozdem_suffixes <- sapply(sozdem_columns, extract_suffix, prefix = "SOZDEM")
sozio_suffixes <- sapply(sozio_columns, extract_suffix, prefix = "SOZIO")
soz_match <- sozio_columns[sozio_suffixes %in% sozdem_suffixes] # this gets us all matching colnames (those asked in paper and in person) 
soz_match_person <- grep("A1_FU1$", soz_match, value = TRUE)
# now do the same as for CESD etc. for the sociodemographic data (insert paper into empty in person cells)
for (col in soz_match_person) {
  # Find the corresponding column (for each "SOZIO...A1_FU1" the matching "SOZDEM...A1_FU1_ANSCH")
  ans_column <- gsub("^SOZIO", "SOZDEM", col)
  ans_column <- gsub("A1_FU1$", "A1_FU1_ANSCH", ans_column)
  # Replace NAs in the current column with values from the corresponding "_ANSCH" column
  fu_datajoin_wide[[col]] <- ifelse(is.na(fu_datajoin_wide[[col]]), fu_datajoin_wide[[ans_column]], fu_datajoin_wide[[col]])
}
# the ATC drug codes from paper fu are in a col called IDOM_AMB_ATC_A1_FU1_ANSCH, those for in person in a col called IDOM_UNT_ATC_A1_FU1
fu_datajoin_wide[["IDOM_UNT_ATC_A1_FU1"]] <- ifelse(is.na(fu_datajoin_wide[["IDOM_UNT_ATC_A1_FU1"]]), fu_datajoin_wide[["IDOM_AMB_ATC_A1_FU1_ANSCH"]], fu_datajoin_wide[["IDOM_UNT_ATC_A1_FU1"]])

# after handling the duplicates and issues resulting from paper/in person fu, we will now calculate the education scores
# the educational attainment was only assessed at BL
# the score is derived from the highest school leaving qualification and the highest further education

# define some freetext equivalents
fhs <- c("Fernstudium stationäre Altenpflege", "Ingenieurfachschulabschluß", "Nachdiplomiert und anerkannt worden", "Ingenieur", "Waldorf-Lehrer",
"Nachdiplomiert, wurde nach der Wende als Diplomabschluss anerkannt", "Soziologie- und Arbeitswissenschaften", 
"Abschluss in Italien erworben an dortiger Universität.", "Damals Fachschule, heute gilt es als Diplom", "Hochschulpädagogik", 
"Fachschulstudium Staatswissenschaften zu DDR-Zeiten (wird heute nicht anerkannt)", "Ingenieurabschluss (Diplomingenieur)", 
"Fachhochschulabschluss mit Nachdiplomierung", "Ingenieurschule", "EDV-Abschluss (Diplom)", "Burufsschullehrer", 
"Fachhochschulabschl. wurde nach der Wende als Dipl.-Ing. anerkannt", "Hochschulabschluss", "Fachschulstudium Sozialarbeiterin", 
"Studium in Peru und Chile", "Sozialarbeiter", "Anerkennung Dipl. Wirtschaftsingenieur", "Konzertexamen", "Miltärakademie", 
"Ingenieurschule = Fachschule zu DDR- Zeiten", "Werkzeugmacher / Fachschule wurde überschrieben auf Diplom", 
"Sonderstudium während Berufstätigkeit, Diplomwirtschaftsingenieur; Brückenkurs zur Betriebswirtin", "Dipl. Ingeneur für Informationstechnik", 
"fachspezifischer Hochabschluss", "Ingenieurschulabschluss zu DDR-Zeiten (\"F\")", "Universitätsabschluss in der USA", "Ingeneurschule")
magprom <- c("Promoviert", "Promotion", "DR.", "Diplomlehrerin", "Promotion Tiermedizin", "Fachpsychologe,pädagogische Psychologei", "Lehrerin", "postgraduales Studium")

bl_datajoin <- bl_datajoin %>%
  mutate(education = case_when(
    SOZIO_F0054 %in% magprom ~ 7.0, # Highest priority: Assign 7.0 if condition is met
    (SOZIO_F0041 == 6 | SOZIO_F0041 == 7) & (SOZIO_F0050 == 1 | SOZIO_F0051 == 1 |  SOZIO_F0054 %in% fhs) ~ 6.1, # Abitur/Fachhoschulreife + (Fach)hochschulabschluss
    (SOZIO_F0041 == 4 | SOZIO_F0041 == 5) & (SOZIO_F0051 == 1 | SOZIO_F0054 %in% magprom) ~ 5.3, # Realschul- or POS-Abschluss and Hochschulabschluss
    SOZIO_F0041 == 3 & (SOZIO_F0051 == 1 | SOZIO_F0054 %in% magprom) ~ 5.0, # Hauptschul- and Hochschulabschluss
    (SOZIO_F0041 == 6 | SOZIO_F0041 == 7) & (SOZIO_F0047 == 1 | SOZIO_F0048 == 1 | SOZIO_F0049 == 1) ~ 4.8, # Abitur/Fachhochschulreife + apprenticeship
    (SOZIO_F0041 == 4 | SOZIO_F0041 == 5) & (SOZIO_F0050 == 1 | SOZIO_F0054 %in% fhs) ~ 4.85, # Realschul- or POS-Abschluss and Fachhochschulabschluss
    SOZIO_F0041 == 3 & (SOZIO_F0050 == 1 | SOZIO_F0054 %in% fhs) ~ 4.55, # Hauptschulabschluss and Fachhochschulabschluss
    (SOZIO_F0041 == 6 | SOZIO_F0041 == 7) & SOZIO_F0046 == 1 ~ 3.7, # Abitur/Fachhochschulreife but no apprenticeship
    (SOZIO_F0041 == 4 | SOZIO_F0041 == 5) & (SOZIO_F0047 == 1 | SOZIO_F0048 == 1 | SOZIO_F0049 == 1) ~ 3.6, # people that fulfill both categories above
    (SOZIO_F0041 == 1 | SOZIO_F0041 == 2 | SOZIO_F0041 == 3) & (SOZIO_F0047 == 1 | SOZIO_F0048 == 1 | SOZIO_F0049 == 1) ~ 3.0, # people without school degrees above but some form of completed apprenticeship
    (SOZIO_F0041 == 4 | SOZIO_F0041 == 5) & SOZIO_F0046 == 1 ~ 2.8, # people with a Realschul- or POS-Abschluss
    SOZIO_F0041 == 3 & SOZIO_F0046 == 1 ~ 1.7, # people with a Hauptschulabschluss
    (SOZIO_F0041 == 1 | SOZIO_F0041 == 2) & SOZIO_F0046 == 1 ~ 1.0 # people without any degree are assigned 1.0
  ))

# professional status was also only assessed at baseline
# different professional groups are assigned to different scores
bl_datajoin <- bl_datajoin %>%
  mutate(profstat = case_when(
    SOZIO_F0072 == "B3" ~ 7.0, # academic with 5+ employees
    SOZIO_F0072 == "B2" ~ 6.8, # academic with 1-4 employees
    SOZIO_F0072 == "D4" ~ 6.4, # higher civil servant
    SOZIO_F0072 == "B1" ~ 5.8, # academic without employees
    SOZIO_F0072 == "D3" ~ 5.2, # elevated civil servant
    SOZIO_F0072 == "E4" ~ 4.7, # employee with managerial tasks
    SOZIO_F0072 == "E3" | SOZIO_F0072 == "C3" | SOZIO_F0072 == "C4" ~ 4.2, # employee with leading role /self-employed with 5+ employees
    SOZIO_F0072 == "D2" ~ 4.1, # middle civil servant
    SOZIO_F0072 == "C2" | SOZIO_F0072 == "E2" ~ 3.6, # self-employed with 1-4 employees / employee with qualified tasks
    SOZIO_F0072 == "C1"  ~ 3.5, # self-employed without employees
    SOZIO_F0072 == "G1" | SOZIO_F0072 == "G2" | SOZIO_F0072 == "G3" | SOZIO_F0072 == "H" | SOZIO_F0072 == "D1" ~ 2.9, # apprentice, simple civil servant, employee in family business
    SOZIO_F0072 == "E1" | SOZIO_F0072 == "F5" ~ 2.4, # master craftsperson / other employees
    SOZIO_F0072 == "F3" ~ 2.1, # craftsperson
    SOZIO_F0072 == "F4" ~ 2.0, # foreperson
    SOZIO_F0072 == "F2" ~ 1.8, # semi-skilled worker
    SOZIO_F0072 == "F1" ~ 1.3, # unskilled worker
    SOZIO_F0072 == "A2" ~ 1.1, # farmer with > 10ha land
    SOZIO_F0072 == "A1" | SOZIO_F0072 == "A3" ~ 1.0, # farmer with < 10ha land / farmer in cooperative
  ))

# the income is turned into an equivalised income based on the OECD method
# net household income / 1 + 0.5 * additional household members >= 14a, + 0.3 * additional household members < 14a
min(bl_datajoin$SOZIO_F0079-bl_datajoin$SOZIO_F0080, na.rm = T) # is 1, there are no misunderstandings that have to be accounted for 
bl_datajoin$equiv_income <- bl_datajoin$SOZIO_F0083 / (1 + 0.5 * (bl_datajoin$SOZIO_F0079 - 1 - bl_datajoin$SOZIO_F0080) + 0.3 * bl_datajoin$SOZIO_F0080)
# many people living alone only stated their income as an individual income but not as a household income
bl_datajoin$equiv_income <- ifelse(bl_datajoin$SOZIO_F0079 == 1, bl_datajoin$SOZIO_F0086, bl_datajoin$equiv_income)
# this equivalised income is again converted into a score ranging from 1.0 to 7.0
bl_datajoin <- bl_datajoin %>%
  mutate(income_score = case_when(
    equiv_income >= 3000 ~ 7.0, 
    equiv_income < 3000 & equiv_income >= 2333.33 ~ 6.5,
    equiv_income < 2333.33 & equiv_income > 2000 ~ 6.0,
    equiv_income <= 2000 & equiv_income >= 1866.66 ~ 5.5,
    equiv_income < 1866.66 & equiv_income >= 1666.66 ~ 5.0,
    equiv_income < 1666.66 & equiv_income >= 1533.33 ~ 4.5,
    equiv_income < 1533.33 & equiv_income >= 1400 ~ 4.0,
    equiv_income < 1400 & equiv_income >= 1333.33 ~ 3.5,
    equiv_income < 1333.33 & equiv_income > 1200 ~ 3.0,
    equiv_income <= 1200 & equiv_income >= 1100 ~ 2.5,
    equiv_income < 1100 & equiv_income >= 980 ~ 2.0,
    equiv_income < 980 & equiv_income >= 800 ~ 1.5,
    equiv_income < 800 ~ 1.0
  ))
# the SES is calculated as the mean of the professional status, the education and the income score
bl_datajoin$SES <- bl_datajoin$profstat + bl_datajoin$education + bl_datajoin$income_score


# to merge baseline and follow-up data and to obtain better readibility of this code, we will now rename all relevant columns
bl_datajoin <- bl_datajoin %>%
  rename(
    ID = TEILNEHMER_SIC,
    LSNS_DATE = LSNS_DATUM,
    CES_D_DATE = CES_D_DATUM,
    GAD7_DATE = GAD7_DATUM,
    GENDER = TEILNEHMER_GESCHLECHT,
    MARITAL_STATUS = SOZIO_F0026,
    EDUCATION = education,
    INCOME_SCORE = income_score,
    PROFESSIONAL_STATUS = profstat,
    N_INHABITANTS = SOZIO_F0079,
    OCCUPATION = SOZIO_F0055,
    NO_OCCUPATION = SOZIO_F0056,
    MEDICATION = ADULT_MEDA_H_ATC,
    BP_SYST_1 = BLUTDRUCKMESS_F0014,
    BP_DIAST_1 = BLUTDRUCKMESS_F0015,
    BP_SYST_2 = BLUTDRUCKMESS_F0019,
    BP_DIAST_2 = BLUTDRUCKMESS_F0020,
    BP_SYST_3 = BLUTDRUCKMESS_F0024,
    BP_DIAST_3 = BLUTDRUCKMESS_F0025,
    HBA1C = HBA1C_E_NUM_VALUE,
    HEIGHT = ANTHRO_F0010,
    WEIGHT = ANTHRO_F0011,
    MOB = TEILNEHMER_GEB_JJJJMM,
    DIABETES_ANAMNESIS = MEDANAM_F0040,
    HYPERTENSION_ANAMNESIS = KARDANAM_F0109
  ) %>%
  rename_with(~ str_replace(., "^GAD7_GAD7", "GAD7"))

fu_datajoin_wide <- fu_datajoin_wide %>%
  rename_with(~ sub("_A1_FU1$", "", .), ends_with("_A1_FU1")) %>%
  rename(
    ID = TEILNEHMER_SIC,
    LSNS_DATE = LSNS_DATUM,
    CES_D_DATE = CES_D_DATUM,
    GAD7_DATE = GAD7_DATUM,
    GENDER = TEILNEHMER_GESCHLECHT,
    MARITAL_STATUS = SOZIO_FAMSTAND,
    N_INHABITANTS = SOZIO_HAUSH_GROE,
    OCCUPATION = SOZIO_ERWERB,
    NO_OCCUPATION = SOZIO_NERWERB_GR,
    MEDICATION = IDOM_UNT_ATC,
    BP_SYST_1 = BP_F0019,
    BP_DIAST_1 = BP_F0020,
    BP_SYST_2 = BP_F0024,
    BP_DIAST_2 = BP_F0025,
    BP_SYST_3 = BP_F0029,
    BP_DIAST_3 = BP_F0030,
    HEIGHT = ANTHRO_GROESSE,
    WEIGHT = ANTHRO_GEWICHT,
    MOB = TEILNEHMER_GEB_JJJJMM,
    DIABETES_ANAMNESIS = MEDIZ_AN_F5,
    HYPERTENSION_ANAMNESIS = KARD_AN_F10
  ) %>%
  rename_with(~ str_replace(., "^GAD7_GAD7", "GAD7"))
  

# create a tp variable to identify baseline and follow-up rows before binding the two dfs together
bl_datajoin$tp <- "bl"
fu_datajoin_wide$tp <- "fu"
df <- bind_rows(bl_datajoin, fu_datajoin_wide)
df <- merge(df, radiological_assessment, all.x = T, by = c("ID","tp"))
df <- merge(df, samseg_data, all.x = T, by = c("ID","tp"))
df <- merge(df, LCL_QA, all.x = T, by = c("ID"))

# SES and education are the same for baseline and FU
# no hypertension/diabetes at FU if it was present at BL is considered to be a mistake and corrected accordingly
df <- df %>%
  # Create a copy of the baseline values (for each participant)
  group_by(ID) %>%
  mutate(bl_ses = SES[tp == "bl"]) %>%
  mutate(bl_ed = EDUCATION[tp == "bl"]) %>%
  ungroup() %>%
  # Now, for the follow-up rows (fu), assign the baseline SES value
  mutate(SES = ifelse(tp == "fu", bl_ses, SES)) %>%
  mutate(EDUCATION = ifelse(tp == "fu", bl_ed, EDUCATION)) %>%
  select(-bl_ses, -bl_ed)  # Drop the temporary columns `bl_ses` and `bl_ed` after they are used

# from the SES scores, we can derive to which SES-quintile of the German population our participants belong
# we use the quintile breaks provided by Lampert et al.
breaks <- c(3.0, 7.7, 9.6, 11.7, 14.1, 21.0)
df$SES_QUINTILE <- cut(df$SES, breaks = breaks, labels = 1:5, right = TRUE)
# turn quintiles into weights
# Remove NA values for weight calculation
QUINTILES_no_NAs <- df[!is.na(df$SES_QUINTILE), ]
# Calculate the frequency of each level (excluding NAs)
level_counts <- table(QUINTILES_no_NAs$SES_QUINTILE)
# Calculate weights
level_weights <- sapply(level_counts, function(x) 0.2/(x/sum(level_counts)))
df$SES_WEIGHT <- level_weights[as.character(df$SES_QUINTILE)]
# give observations where SES is NA a weight of 1
df$SES_WEIGHT <- ifelse(is.na(df$SES_WEIGHT), 1,  df$SES_WEIGHT)


# as baseline and fu are now in one df, we can now calculate some variables 
# let's start with BMI
df$BMI <- df$WEIGHT / ((df$HEIGHT/100)^2)
df$BMI <- ifelse(df$BMI > 50 | df$BMI < 15, NA, df$BMI) # turn very low/high values into NA
# high degree of missingness is confronted with last observation carried forward
# Replace BMI for "fu" rows where BMI is NA with the corresponding "bl" BMI value
df <- df %>%
  left_join(
    filter(df, tp == "bl") %>% select(ID, BMI_bl = BMI),
    by = "ID"
  ) %>%
  mutate(BMI = ifelse(tp == "fu" & is.na(BMI), BMI_bl, BMI)) %>%
  select(-BMI_bl)  # Drop the BMI_bl column


# let's continue with the sum scores of LSNS, CES-D and GAD7
df$LSNS_SUM <- rowSums(df[,c(paste0("LSNS_", c(1:6)))])
df$CES_D_SUM <- rowSums(df[,c(paste0("CES_D_", c(1:20)))])
df$GAD7_SUM <- rowSums(df[,c(paste0("GAD7_", c(1:7)))])

# we will assign 0 / 1 to employed / unemployed participants 
df$EMPLOYED <- ifelse(df$OCCUPATION == 0, 1, 0) # assign 1 (i.e. unemployed) to people who stated that they are not gainfully employed
df$EMPLOYED <- ifelse(df$OCCUPATION == 98 | df$OCCUPATION == 99, NA, df$EMPLOYED) # turn codes into NAs
df$EMPLOYED <- ifelse(df$EMPLOYED == 1 & (df$NO_OCCUPATION == 1 | df$NO_OCCUPATION == 6), 0, df$EMPLOYED) # unless they are not gainfully employed because of studying/voluntary service   

# we will assign 0 / 1 to unmarried / married participants
df$MARRIED <- ifelse(df$MARITAL_STATUS == 1, 1, 0)
df$MARRIED <- ifelse(df$MARITAL_STATUS == 98 | df$MARITAL_STATUS == 99, NA, df$MARRIED) # turn codes into NAs

# we will assign 0 / 1 to participants living alone / cohabitating in some form
df$N_INHABITANTS <- ifelse(df$N_INHABITANTS == 98 | df$N_INHABITANTS == 99, NA, df$N_INHABITANTS) # turn codes into NA
df <- df %>%
  mutate(LIVE_ALONE = case_when(
    is.na(df$N_INHABITANTS) ~ NA_real_,
    df$N_INHABITANTS == 1 ~ 1, # those without cohabitees are coded as 1
    TRUE ~ 0 # everyone living with at least one person is coded as 0
  ))

# we may now calculate the participants' age
# first try LSNS date
df$AGE <- as.numeric(difftime(as.Date(df$LSNS_DATE), as.Date(paste0(substr(df$MOB, 1, 4), "-", substr(df$MOB, 5, 6), "-01")), units = "days")) / 365.25
# backup is MRI date
df$AGE <- ifelse(is.na(df$AGE), 
                 as.numeric(difftime(as.Date(df$Messung.Datum_TI), as.Date(paste0(substr(df$MOB, 1, 4), "-", substr(df$MOB, 5, 6), "-01")), units = "days")) / 365.25,
                 df$AGE)
# 2nd backup is date of CES-D
df$AGE <- ifelse(is.na(df$AGE), difftime(as.Date(df$CES_D_DATE), as.Date(paste0(substr(df$MOB, 1, 4), "-", substr(df$MOB, 5, 6), "-01")) , units = "days") / 365.25, df$AGE)
# 3rd backup is date of GAD7
df$AGE <- ifelse(is.na(df$AGE), difftime(as.Date(df$GAD7_DATE), as.Date(paste0(substr(df$MOB, 1, 4), "-", substr(df$MOB, 5, 6), "-01")) , units = "days") / 365.25, df$AGE)

# we will modify the gender variable so that females are our reference value (0) while males keep being 1
df$GENDER <- ifelse(df$GENDER == 2, 0, df$GENDER)
# we will also fill in two NAs from FU with their baseline values
df <- df %>%
  left_join(df %>% filter(tp == "bl") %>% select(ID, GENDER), by = "ID", suffix = c("", "_bl")) %>%
  mutate(GENDER = ifelse(tp == "fu" & is.na(GENDER), GENDER_bl, GENDER)) %>%
  select(-GENDER_bl) 
# for further variables we have to check whether participants take antidiabetic/antihypertensive medication
# remove hashtags and spaces around drug names
df$MEDICATION <- str_replace_all(df$MEDICATION, c("#" = "", " " = "")) 
# in the ATC nomenclature all general antihypertensive drugs begin with "C02", all diuretics begin with "C03", 
# all beta-blockers with "C07", all Ca-chanel blockers with "C08" and all RAAS affecting drugs with "C09"  
# all antidiabetic drugs begin with "A10" in the ATC nomenclature
# we want to except SGLT2-Inhibitors ( starting with "A10BK") as they have a broader set of indications and their use as an antidiabetic monotherapy is unlikely
# Wegovy (semaglutide for weight loss) actually has a different ATC-code
df <- df %>%
  mutate(BP_MED = case_when(
    grepl("^C02", MEDICATION) | grepl(",C02", MEDICATION) ~ 1,
    grepl("^C03", MEDICATION) | grepl(",C03", MEDICATION) ~ 1,
    grepl("^C07", MEDICATION) | grepl(",C07", MEDICATION) ~ 1,
    grepl("^C08", MEDICATION) | grepl(",C08", MEDICATION) ~ 1,
    grepl("^C09", MEDICATION) | grepl(",C09", MEDICATION) ~ 1,
    is.na(MEDICATION) ~ NA_real_,
    TRUE ~ 0
  )) %>%
  mutate(DIABETES_MED = case_when(
    grepl("^A10(?!BK)", MEDICATION, perl = TRUE) | grepl(",A10(?!BK)", MEDICATION, perl = TRUE) ~ 1,
    is.na(MEDICATION) ~ NA_real_,
    TRUE ~ 0
  ))

# to be categorised as suffering from diabetes, participants must fulfill one of three criteria
df <- df %>%
  mutate(DIABETES = case_when(
    DIABETES_MED == 1 ~ 1,# they must take antidiabetic medication
    HBA1C >= 6 ~ 1,# or have an HBa1c >= 6%
    DIABETES_ANAMNESIS == 1 ~ 1,# or have stated that they received a diabetes diagnosis in the anamnesis
    is.na(DIABETES_MED) & is.na(HBA1C) & is.na(DIABETES_ANAMNESIS) ~ NA_real_,
    TRUE ~ 0
  ))

# first calculate average systolic and diastolic blood pressure
df$BP_SYST <- rowSums(df[,c(paste0("BP_SYST_", c(1:3)))])/3 
df$BP_SYST <- ifelse(df$BP_SYST > 200, NA, df$BP_SYST) # turn very high values to NA as they are most likely to result from faulty measurements
df$BP_SYST <- ifelse(df$BP_SYST < 80, NA, df$BP_SYST) # same for very low values
df$BP_DIAST <- rowSums(df[,c(paste0("BP_DIAST_", c(1:3)))])/3
# to be categorised as suffering from hypertension, participants must fulfill one of three criteria
df <- df %>%
  mutate(HYPERTENSION = case_when(
    BP_MED == 1 ~ 1,# they must take antihypertensive medication
    df$BP_SYST > 140 | df$BP_DIAST > 90 ~ 1,# an average BPsyst > 140mmHg or BPdisat > 90mmHg
    HYPERTENSION_ANAMNESIS == 1 ~ 1,# or have stated that they received a hypertension diagnosis in the anamnesis
    is.na(BP_MED) & is.na(HBA1C) & is.na(HYPERTENSION_ANAMNESIS) ~ NA_real_,
    TRUE ~ 0
  ))
# no hypertension/diabetes at FU if it was present at BL is considered to be a mistake and corrected accordingly
df <- df %>%
  # Create a copy of the baseline values (for each participant)
  group_by(ID) %>%
  mutate(bl_hypertension = HYPERTENSION[tp == "bl"]) %>%
  mutate(bl_diabetes = DIABETES[tp == "bl"]) %>%
  ungroup() %>%
  # Now, for the follow-up rows (fu), assign the baseline values if FU is 0
  mutate(HYPERTENSION = ifelse(tp == "fu" & HYPERTENSION == 0, bl_hypertension, HYPERTENSION)) %>%
  mutate(DIABETES = ifelse(tp == "fu" & DIABETES == 0, bl_diabetes, DIABETES)) %>%
  select(-bl_hypertension, -bl_diabetes)  # Drop the temporary columns after they are used

# determine groups of columns belonging to the different cognitive tests 
b <- c("BUTTER", "ARM", "STRAND", "BRIEF", "KONIGI", "HUTTE", "STANGE", "KARTE", "GRAS", "MOTOR")
c <- c("KIRCHE", "KAFFEE", "BUTTER", "DOLLAR", "ARM", "STRAND", "FUNF", "BRIEF", "HOTEL", "BERG", "KONIGI", "HUTTE", 
       "PANTOF", "STANGE", "DORF", "BAND", "KARTE", "HEER", "GRAS", "MOTOR")
learning_cols_1 <- paste0("CERAD_WL1_", b)
learning_cols_2 <- paste0("CERAD_WL2_", b)
learning_cols_3 <- paste0("CERAD_WL3_", b)
recall_cols <- paste0("CERAD_WL4_", b)
recognition_cols <- paste0("CERAD_WLW_", c)
phon_flu_cols <- paste0("CERAD_S_", c(15, 30, 45, 60), "S_CORR")
sem_flu_cols <- paste0("CERAD_VF_", c(15, 30, 45, 60), "S_CORR")
# replace various codes indicating that values are not usable with NA
df <- replace_with_na_at(df, .vars = c("TMT_TIMEA", "TMT_TIMEB", sem_flu_cols), condition = ~.x %in% c(997, 998, 999))
df <- replace_with_na_at(df, .vars = c(recognition_cols, phon_flu_cols), condition =  ~.x %in% c(97, 98, 99))
# calculate scores for different iterations of the learning test
df$learning1 <- rowSums(df[,learning_cols_1]) 
df$learning2 <- rowSums(df[,learning_cols_2]) 
df$learning3 <- rowSums(df[,learning_cols_3]) 
df$recall <- rowSums(df[, recall_cols])
# see if the iterations have to be replaced by NA
df <- df %>%
  mutate(learning1 = if_else(CERAD_WL1_97 == 1 | CERAD_WL1_98 == 1 | CERAD_WL1_99 == 1, NA_real_, learning1)) %>%
  mutate(learning2 = if_else(CERAD_WL2_97 == 1 | CERAD_WL2_98 == 1 | CERAD_WL2_99 == 1, NA_real_, learning2)) %>%
  mutate(learning3 = if_else(CERAD_WL3_97 == 1 | CERAD_WL3_98 == 1 | CERAD_WL3_99 == 1, NA_real_, learning3)) %>%
  mutate(recall = if_else(CERAD_WL4_97 == 1 | CERAD_WL4_98 == 1 | CERAD_WL4_99 == 1, NA_real_, recall))
# calculate scores for different cognitive tests
df$learning <- df$learning1 + df$learning2 + df$learning3
df$recognition <- rowSums(df[, recognition_cols])
df$sem_flu <- rowSums(df[, sem_flu_cols])
df$phon_flu <- rowSums(df[,phon_flu_cols])
# calculate TMT composite score and negate it to facilitate later calculations
df$TMT <- -((df$TMT_TIMEB - df$TMT_TIMEA)/df$TMT_TIMEA)
# standardise the cognitive test results using the baseline mean and SD
df$learning_z <- (df$learning - mean(df[df$tp == "bl",]$learning, na.rm = T))/sd(df[df$tp == "bl",]$learning, na.rm = T)  
df$recognition_z <- (df$recognition - mean(df[df$tp == "bl",]$recognition, na.rm = T))/sd(df[df$tp == "bl",]$recognition, na.rm = T)  
df$recall_z <- (df$recall - mean(df[df$tp == "bl",]$recall, na.rm = T))/sd(df[df$tp == "bl",]$recall, na.rm = T)  
df$sem_flu_z <- (df$sem_flu - mean(df[df$tp == "bl",]$sem_flu, na.rm = T))/sd(df[df$tp == "bl",]$sem_flu, na.rm = T)  
df$phon_flu_z <- (df$phon_flu - mean(df[df$tp == "bl",]$phon_flu, na.rm = T))/sd(df[df$tp == "bl",]$phon_flu, na.rm = T)  
df$TMT_z <- (df$TMT - mean(df[df$tp == "bl",]$TMT, na.rm = T))/sd(df[df$tp == "bl",]$TMT, na.rm = T)  
df$TMT_TIMEA_z <- (df$TMT_TIMEA - mean(df[df$tp == "bl",]$TMT_TIMEA, na.rm = T))/sd(df[df$tp == "bl",]$TMT_TIMEA, na.rm = T)  

# calculate composite scores for cognitive functions
df$EF <- (df$phon_flu_z + df$sem_flu_z + df$TMT_z)/3
df$MEMO <- (df$learning_z + df$recall_z + df$recognition_z)/3
df$PS <- -df$TMT_TIMEA_z
df$CERAD <- df$learning_z + df$recall_z + df$recognition_z + df$phon_flu_z + df$sem_flu_z + df$TMT_z

# z-transform the cognitive scores for better interpretability of estimates
df <- df %>%
  mutate(across(c(EF, MEMO, PS, CERAD), ~ as.vector(scale(.))))

# let's now calculate the average hippocampus volume corrected by ICV
df$MEAN_HCV <- (df$`Left-Hippocampus` + df$`Right-Hippocampus`)/2 # first average HCV over both hemispheres
MEAN_ICV <- mean(df$sbTIV, na.rm = T) # calculate average ICV over all participants
res <- lmer(MEAN_HCV ~ sbTIV + (1|ID), data = df, REML = F, na.action = na.omit) # determine regression coefficient of ICV on HCV 
sres <- summary(res)
beta <- sres$coefficients[2]
df$HCV_ADJ <- df$MEAN_HCV - beta * (df$sbTIV - MEAN_ICV) # adjust HCV using beta from LME

# create a variable that indicates whether ventricle expansion was coded in the QA
df <- df %>%
  mutate(
    VENTRICLE_EXPANSION = case_when(
      tp == "bl" & qa_LST_bl == 2 ~ 1,  # For baseline with qa_LST_bl = 2
      tp == "fu" & qa_LST_fu == 2 ~ 1,  # For follow-up with qa_LST_fu = 2
      TRUE ~ 0  # Otherwise, set to 0
    )
  )

# create a variable that indicates whether white matter measure quality was coded as unusable in the QA
df <- df %>%
  mutate(
    WMHV_PROBLEMS = case_when(
      tp == "bl" & (qa_LST_bl == 1 | qa_LST_bl == 3) ~ 1,  # For baseline with qa_LST_bl = 1/3
      tp == "fu" & (qa_LST_bl == 1 | qa_LST_bl == 3) ~ 1,  # For follow-up with qa_LST_fu = 1/3
      TRUE ~ 0  # Otherwise, set to 0
    )
  )
# create a variable that indicates whether general MRI quality was coded as unusable in the radiological QA
df <- df %>%
  mutate(
    RADIO_PROBLEMS = case_when(
      is.na(med_Befund.Bewertung_Bef) ~ NA_real_, # NAs are NAs
      med_Befund.Bewertung_Bef == "verwendungsf\xe4hig nein" ~ 1,  # those deemed unusable are coded as 1
      med_Befund.Bewertung_Bef == "verwendungsf\xe4hig voll" ~ 0,  # those deemed usable are coded as 0
      med_Befund.Bewertung_Bef == "verwendungsf\xe4hig eingeschr\xe4nkt" &
      (is.na(med_Befund.Nebenbefund) | (
      med_Befund.Nebenbefund != "(post)isch\xe4misch" &
      med_Befund.Nebenbefund != "(post)h\xe4morrhagisch" &
      med_Befund.Nebenbefund != "(post)isch\xe4misch (post)h\xe4morrhagisch" &
      med_Befund.Nebenbefund != "(post)traumatisch" &
      med_Befund.Nebenbefund != "(post)h\xe4morrhagisch (post)traumatisch" &
      med_Befund.Nebenbefund != "DVA (post)h\xe4morrhagisch" &
      med_Befund.Nebenbefund != "DVA (post)isch\xe4misch" &
      med_Befund.Nebenbefund != "Cavernom (post)h\xe4morrhagisch" &
      med_Befund.Nebenbefund != "DVA kontrollbed\xfcrftig (post)h\xe4morrhagisch")) ~ 0,  # those deemed to have limited usability are coded as 0 if this is due to something else than a tumor or ischaemic, hemorrhagic or traumatic lesion 
      TRUE ~ 1  # Otherwise, set to 1
    )
  )
# visually test if variables are normally distributed
ggplot(gather(df[, c("LSNS_SUM", "BMI", "CES_D_SUM", "GAD7_SUM", "MEMO", "SES",
                     "PS", "EF", "AGE", "LESIONS", "HCV_ADJ")]), aes(value)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~key, scales = 'free') 
# CES_D can resemble a normal distribution more closely after an asinh transformation
# this will be done in the analysis script after imputation

# create a variable that indicated whether a participant was in the sub-cohort with extensive cognitive testing/MRIs
# at baseline everybody with an MRI measurement or some cognitive test result other than the TMT / semantic fluency are in the cognitive subcohort
# at FU the TMT and semantic fluency was also only performed by participants in this subcohort
df <- df %>%
  mutate(
    COG_COHORT = case_when(
      tp == "bl" & (if_any(all_of(learning_cols_1), ~ !is.na(.)) | if_any(all_of(learning_cols_2), ~ !is.na(.)) | 
                      if_any(all_of(learning_cols_3), ~ !is.na(.)) | if_any(all_of(recall_cols), ~ !is.na(.)) |
                    if_any(all_of(recognition_cols), ~ !is.na(.)) | if_any(all_of(phon_flu_cols), ~ !is.na(.)) |
                      !is.na(HCV_ADJ)) ~ 1,
      tp == "fu" & (if_any(all_of(learning_cols_1), ~ !is.na(.)) | if_any(all_of(learning_cols_2), ~ !is.na(.)) | 
                      if_any(all_of(learning_cols_3), ~ !is.na(.)) | if_any(all_of(recall_cols), ~ !is.na(.)) |
                      if_any(all_of(recognition_cols), ~ !is.na(.)) | if_any(all_of(phon_flu_cols), ~ !is.na(.)) |
                      if_any(all_of(sem_flu_cols), ~ !is.na(.)) | !is.na(TMT_TIMEA) | !is.na(HCV_ADJ)) ~ 1,
      TRUE ~ 0))

# every participant whose ID appears in the samseg_data file at least once is in the MRI cohort
df <- df %>%
  mutate(MRI_COHORT = if_else(ID %in% samseg_data$ID, 1, 0))

# replace string IDs with unique numbers
df <- df %>%
  mutate(ID = dense_rank(ID))

# write uncentered/non-z-scored data to csv
write_csv(as.data.frame(df), "/data/pt_life/ResearchProjects/LLammer/gamms/Data/assembled_unscaled_data.csv")

# center for better interpretability of estimates
df$AGE <- df$AGE - mean(df$AGE, na.rm = T)
# center for better interpretability of estimates
df$BMI <- df$BMI - mean(df$BMI, na.rm = T)
# scale to improve interpretability of estimates
df$SES <- (df$SES - mean(df$SES, na.rm = T))/sd(df$SES, na.rm = T)

write_csv(as.data.frame(df), "/data/pt_life/ResearchProjects/LLammer/gamms/Data/assembled_data.csv")
