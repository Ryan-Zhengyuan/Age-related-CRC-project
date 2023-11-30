library(dplyr)
library(tibble)
library(stringr)
library(DescTools)
'%!in%' <- function(x,y)!('%in%'(x,y))

### Metastudy - Nicola Segata 
## This Study contains both Ccancer and COntrols

Cancer_1 <-
   read.csv("~/Metadata/merged_metadata_CRC.tsv",
            sep = "\t") %>% 
   dplyr::select(sample_id,
                 subject_id,
                 disease,
                 study_condition,
                 age,
                 gender,
                 BMI,
                 country,
                 NCBI_accession) %>% # column select
   dplyr::filter(str_detect(disease, 'CRC')) %>% # select CRC disease
   dplyr::filter(is.na(age) == FALSE) %>%  # remove age with no data
   dplyr::filter(age != "unknown") %>% 
   dplyr::filter(age != "adult") %>% 
   dplyr::filter(is.na(BMI) == FALSE) %>% # remove BMI with no data
   dplyr::filter(BMI != "unknown") %>% 
   dplyr::mutate(country = gsub("AUT", "AUS",country)) %>% # replace country name from 'AUT' to 'AUS'
   dplyr::select(-disease) %>% # remove column "disease"
   tidyr::drop_na() %>% # remove any rows containing blank
   dplyr::mutate(sample_id = as.character(sample_id)) %>% 
   dplyr::mutate(subject_id = as.character(subject_id)) %>% # transferring id from string to character
   dplyr::mutate(country = gsub("DEU", "GER",country)) %>% 
   dplyr::mutate(sample_id = gsub("sub_", "",sample_id)) %>% 
   dplyr::select(-subject_id) %>% # remove column "subject_id"
   dplyr::rename(subject_id = "sample_id")  %>% # replace sample_id to subject_id
   mutate(NCBI_accession = gsub(";.*", "", NCBI_accession)) %>% 
   mutate(subject_id = ifelse(is.na(NCBI_accession) == T,as.character(subject_id), NCBI_accession )) %>% # replace subject_id to NCBI_accession
   dplyr::select(-NCBI_accession)



Control_1 <-
   read.csv("~/Metadata/merged_metadata_CRC.tsv",
            sep = "\t") %>% 
   dplyr::select(sample_id,
                 subject_id,
                 study_condition,
                 age,
                 gender,
                 BMI,
                 country,
                 NCBI_accession) %>% 
   dplyr::filter(str_detect(study_condition, 'control')) %>% # select control samples
   dplyr::filter(is.na(age) == FALSE) %>% 
   dplyr::filter(age != "unknown") %>% 
   dplyr::filter(is.na(BMI) == FALSE) %>% 
   dplyr::mutate(country = gsub("AUT", "AUS",country)) %>% 
   dplyr::mutate(country = gsub("DEU", "GER",country)) %>% 
   dplyr::filter(BMI != "unknown") %>% 
   dplyr::mutate(sample_id = as.character(sample_id)) %>% 
   dplyr::mutate(subject_id = as.character(subject_id)) %>% 
   dplyr::mutate(sample_id = gsub("sub_", "",sample_id)) %>% 
   dplyr::select(-subject_id) %>% 
   dplyr::rename(subject_id = "sample_id") %>% 
   mutate(NCBI_accession = gsub(";.*", "", NCBI_accession)) %>% 
   mutate(subject_id = ifelse(is.na(NCBI_accession) == T,as.character(subject_id), NCBI_accession )) %>% 
   dplyr::select(-NCBI_accession)





### Metagenomic and metabolomic analyses reveal distinct stage-specific phenotypes of the gut microbiota in colorectal cancer

Cancer_2 <-
   read.csv("~/Metadata/Metadata_Cancer_2.tsv",
            sep = "\t") %>% 
   dplyr::select(-c(Alcohol,Brinkman.Index,Tumor.location,Stage)) %>% 
   dplyr::rename(subject_id = "Subject_ID",
                 study_condition = "Group",
                 gender = "Gender",
                 age = "Age") %>%
   dplyr::mutate(gender = ifelse(gender == "M","male","female")) %>% # replace 'M' to male and no 'm' to female
   dplyr::filter(study_condition != "Healthy"& study_condition != "MP" & study_condition != "HS") %>% # filter samples with CRC
   dplyr::rename(Cancer_Stage = "study_condition") %>% 
   dplyr::mutate(study_condition = rep('CRC',nrow(.))) %>% # add 'study_condition' column filled with 'CRC'
   dplyr::mutate(country = rep('JPN',nrow(.))) %>% # add 'country' column filled with 'JPN' 
   tidyr::drop_na() %>% 
   dplyr::mutate(Cancer_Stage = as.character(Cancer_Stage)) %>% 
   dplyr::mutate(Cancer_Stage = ifelse(str_detect("Stage_I_II",Cancer_Stage),"stages I and II",
                                       ifelse(str_detect("Stage_III_IV",Cancer_Stage), "stages III and IV", "stage 0"))) %>% 
   dplyr::select(-Cancer_Stage)


Control_2 <-
   read.csv("~/Metadata/Metadata_Cancer_2.tsv",
            sep = "\t") %>% 
   dplyr::select(-c(Alcohol,Brinkman.Index,Tumor.location,Stage)) %>% 
   dplyr::rename(subject_id = "Subject_ID",
                 study_condition = "Group",
                 gender = "Gender",
                 age = "Age") %>% 
   dplyr::mutate(gender = ifelse(gender == "M","male","female")) %>% 
   dplyr::filter(BMI != "unknown") %>% 
   dplyr::filter(study_condition == "Healthy") %>% 
   dplyr::rename(Cancer_Stage = "study_condition") %>% 
   dplyr::mutate(study_condition = rep('control',nrow(.))) %>% 
   dplyr::mutate(country = rep('JPN',nrow(.))) %>% 
   tidyr::drop_na() %>% 
   dplyr::mutate(Cancer_Stage = gsub("Healthy", "control",Cancer_Stage)) 



### Multi-kingdom microbiota analyses identify bacterial–fungal interactions and biomarkers of colorectal cancer across cohorts

DB_1 <-
   read.csv("~/Metadata/merged_metadata_CRC.tsv",
            sep = "\t") %>% 
dplyr::select(sample_id,NCBI_accession) %>% 
   mutate(NCBI_accession = gsub(";.*", "", NCBI_accession)) %>% 
   dplyr::rename(subject_id = "sample_id")


Cancer_3 <-
   read.csv("~/Metadata/Metadata.tsv",
            sep = "\t") %>% 
   dplyr::select(SampleID,Group,
                 Gender,Age,BMI,Disease.Stage,country) %>% 
   dplyr::rename(gender = "Gender", age = "Age", Cancer_Stage = "Disease.Stage", 
                 study_condition = "Group", subject_id = "SampleID") %>% 
   dplyr::mutate(gender = ifelse(gender == "M","male","female")) %>% 
   dplyr::filter(study_condition == "CRC") %>% 
   dplyr::mutate(country = gsub("JAP", "JPN",country)) %>% 
   dplyr::mutate(country = gsub("ITA1", "ITA",country)) %>% 
   dplyr::mutate(country = gsub("China", "CHN",country)) %>% 
   dplyr::mutate(country = gsub("CHI", "CHN",country)) %>% 
   tidyr::drop_na() %>% 
   dplyr::select(-Cancer_Stage) %>% 
   left_join(DB_1, by = "subject_id") %>% # add NCBI ID on the left based on subject_id
   mutate(subject_id = ifelse(is.na(NCBI_accession) == T,as.character(subject_id), NCBI_accession )) %>% 
   dplyr::select(-NCBI_accession)





Control_3 <-
   read.csv("~/Metadata/Metadata.tsv",
            sep = "\t") %>% 
   dplyr::select(SampleID,Group,
                 Gender,Age,BMI,Disease.Stage,country) %>% 
   dplyr::rename(gender = "Gender", age = "Age", Cancer_Stage = "Disease.Stage", 
                 study_condition = "Group", subject_id = "SampleID") %>% 
   dplyr::mutate(gender = ifelse(gender == "M","male","female")) %>% 
   dplyr::filter(study_condition == "CTR") %>% 
   dplyr::mutate(country = gsub("JAP", "JPN",country)) %>% 
   dplyr::mutate(country = gsub("ITA1", "ITA",country)) %>% 
   dplyr::mutate(country = gsub("China", "CHN",country)) %>% 
   dplyr::mutate(country = gsub("CHI", "CHN",country)) %>% 
   tidyr::drop_na() %>% 
   dplyr::mutate(study_condition = gsub("CTR", "control",Cancer_Stage)) %>% 
   dplyr::mutate(Cancer_Stage = gsub("Early", "stages I and II",Cancer_Stage)) %>% 
   dplyr::mutate(Cancer_Stage = gsub("Advanced", "stages III and IV",Cancer_Stage)) %>% 
   dplyr::mutate(Cancer_Stage = gsub("CTR", "control",Cancer_Stage)) %>% 
   left_join(DB_1, by = "subject_id") %>% 
   mutate(subject_id = ifelse(is.na(NCBI_accession) == T,as.character(subject_id), NCBI_accession)) %>% 
   dplyr::select(-c(NCBI_accession,Cancer_Stage))


Cancer_1_t<-Cancer_1 %>% 
   dplyr::mutate(study_condition = as.character(study_condition)) %>% 
   dplyr::mutate(age = as.numeric(age)) %>% 
   dplyr::mutate(BMI = as.numeric(BMI)) %>% 
   dplyr::mutate(gender = as.character(gender)) 

Cancer_2_t <-Cancer_2 %>%
   dplyr::mutate(subject_id = as.character(subject_id)) %>%
   dplyr::mutate(study_condition = as.character(study_condition)) %>%
   dplyr::mutate(age = as.numeric(age)) %>%
   dplyr::mutate(BMI = as.numeric(as.character(BMI))) 


Cancer_3_t <-Cancer_3 %>% 
   dplyr::mutate(subject_id = as.character(subject_id)) %>% 
   dplyr::mutate(study_condition = as.character(study_condition)) %>% 
   dplyr::mutate(BMI = as.numeric(BMI)) %>%  
   dplyr::mutate(age = as.numeric(age)) 



### Keeping the unique samples from each study , focusing on study 3 which as meta-study contains information from the other 2

All_Cancer <- 
   dplyr::bind_rows(Cancer_3_t,
                    Cancer_1_t
   ) %>% 
   dplyr::mutate(Group = c(rep("Cancer 3",nrow(Cancer_3_t)),
                           rep("Cancer 1",nrow(Cancer_1_t)))) %>% # add the research that the sample came from
   dplyr::distinct(subject_id,
                   study_condition,
                   gender,
                   country,
                   .keep_all = T) %>% # de-replicating the data (only dereplicate data on the basis of subject_id is ok)    
   dplyr::mutate(continent = ifelse(country == "AUS"|country == "GER"|country == "FRA"| country == "ITA","Europe",
                                    ifelse(country == "CAN"|country == "USA","America","Αsia")) %>% as.factor()) %>% 
   dplyr::mutate(real_age = age) %>% 
   dplyr::mutate(age = ifelse(age <= 50,0,1)) %>% 
   dplyr::mutate(gender = ifelse(gender == "male",1,2) %>% as.factor()) %>% 
   dplyr::mutate(country = as.factor(country))

library(MatchIt)
library(tableone)
set.seed(123)

# testing the influence of other factors on the data analyses through calculating standardized mean difference (SMD)

##################################################### CRC_age_matching##############################################################
confounders <- c('gender', 'BMI', 'country', 'continent') 
table1 <- CreateTableOne(var=confounders, strata='age', data=All_Cancer, test=TRUE)
print(table1, smd=TRUE) # p<0.05 and SMD>0.1 together may represent the big influence of these factors on the experiment

matched_3 <- matchit(age ~ gender + BMI + continent, data = All_Cancer, method = "optimal", replace = FALSE, distance = "mahalanobis",
                     within = list(BMI=3), ratio=5)
# Assessing the result of MatchIt
summary(matched_3, un = T)

the_best_match_age_country <- function(matched) {
  
  #matched = matched_3
  matched_C <-  match.data(matched) %>% dplyr::filter(age == 0) %>% as.data.frame()
  matched_T <-  match.data(matched) %>% dplyr::filter(age == 1) %>% as.data.frame()
  
  
  true_matched_T = data.frame()
  for (i in 1:nrow(matched_C)) {
   
    
     print(i)
    current_matched <-
      matched_T %>% 
      dplyr::filter(subclass == i) %>% 
      filter(gender == matched_C %>% dplyr::filter(subclass == i) %>% pull(gender)) %>% 
      filter(country == matched_C %>% dplyr::filter(subclass == i) %>% pull(country) %>% as.character()) %>% 
      filter(BMI == DescTools::Closest(BMI, matched_C %>% dplyr::filter(subclass == i) %>% pull(BMI))) %>% 
      as.data.frame()
    
    
    true_matched_T <- rbind(true_matched_T,
                            current_matched)
  }
  
  final_data <- true_matched_T %>% 
    left_join(matched_C, by = "subclass") %>% 
    dplyr::select(subclass,
                  age.x,age.y,BMI.x,BMI.y,country.x,country.y,
                  continent.x,continent.y,
                  gender.x,gender.y,
                  real_age.x,real_age.y,
                  subject_id.x,subject_id.y)

  return(final_data)
}
the_best_match_age_continent <- function(matched) {
  
  #matched = matched_2
  matched_C <-  match.data(matched) %>% dplyr::filter(age == 0) %>% as.data.frame()
  matched_T <-  match.data(matched) %>% dplyr::filter(age == 1) %>% as.data.frame()
  
  
  true_matched_T = data.frame()
  for (i in 1:nrow(matched_C)) {
    
    
    print(i)
    current_matched <-
      matched_T %>% 
      dplyr::filter(subclass == i) %>% 
      filter(gender == matched_C %>% dplyr::filter(subclass == i) %>% pull(gender)) %>% 
      filter(continent == matched_C %>% dplyr::filter(subclass == i) %>% pull(continent) %>% as.character()) %>% 
      filter(BMI == DescTools::Closest(BMI, matched_C %>% dplyr::filter(subclass == i) %>% pull(BMI))) %>% 
      as.data.frame()
    
    
    true_matched_T <- rbind(true_matched_T,
                            current_matched)
  }
  
  final_data <- true_matched_T %>% 
    left_join(matched_C, by = "subclass") %>% 
    dplyr::select(subclass,
                  age.x,age.y,BMI.x,BMI.y,country.x,country.y,
                  continent.x,continent.y,
                  gender.x,gender.y,
                  real_age.x,real_age.y,
                  subject_id.x,subject_id.y)
  
  return(final_data)
}

############################## Extracting CRC-age matching ###################################

matched_3_res<- 
  the_best_match_age_continent(matched_3) %>% 
  filter(gender.x == gender.y) %>% 
  filter(continent.x == continent.y) %>% 
  filter(abs(BMI.x-BMI.y) <=3) 



####3 Importing the data rfom LU Controls


Controls_Lu_Continents <-
   read.csv("~/Controls/combine_metadata_sel.tsv", sep = '\t') %>% 
   dplyr::select(internal_sample_id,
                 study_name,
                 sample_id,
                 subject_id,
                 study_condition,
                 disease,
                 age,
                 gender,
                 country,
                 BMI) %>% 
   dplyr::filter(is.na(BMI) == FALSE) %>% 
   dplyr::filter(is.na(gender) == FALSE) %>% 
   dplyr::filter(is.na(age) == FALSE) %>% 
   dplyr::mutate(country = gsub("AUT", "AUS",country)) %>% 
   dplyr::mutate(country = gsub("DEU", "GER",country)) %>% 
   filter(country == "FRA" | country == "CHN" | country == "AUS" |
             country == "GER" | country == "JPN" | country == "USA" | country == "ITA" | country == "CAN") %>% 
   dplyr::mutate(continent = ifelse(country == "CHN"|country == "JPN","Αsia",
                                    ifelse(country == "CAN"|country == "USA","America","Europe")) %>% as.factor()) %>% 
   dplyr::select(sample_id,age,gender,country,BMI,continent) %>% 
   dplyr::mutate(age = as.numeric(age)) %>% 
   dplyr::mutate(gender = as.character(gender)) %>% 
   dplyr::rename(subject_id = "sample_id") 




Controls_Lu_Countries <-
   read.csv("~/Controls/combine_metadata_sel.tsv", sep = '\t') %>% 
   dplyr::select(internal_sample_id,
                 study_name,
                 sample_id,
                 subject_id,
                 study_condition,
                 disease,
                 age,
                 gender,
                 country,
                 BMI) %>% 
   dplyr::filter(is.na(BMI) == FALSE) %>% 
   dplyr::filter(is.na(gender) == FALSE) %>% 
   dplyr::filter(is.na(age) == FALSE) %>% 
   dplyr::mutate(country = gsub("AUT", "AUS",country)) %>% 
   dplyr::mutate(country = gsub("DEU", "GER",country)) %>% 
   dplyr::filter(country !=  "CMR" & country !=  "PER" & country !=  "MDG") %>% 
   dplyr::mutate(continent = ifelse(country == "CHN"|country == "JPN","Αsia",
                                    ifelse(country == "CAN"|country == "USA","America","Europe")) %>% as.factor()) %>% 
   dplyr::select(sample_id,age,gender,BMI,country,continent) %>% 
   dplyr::mutate(age = as.numeric(age)) %>% 
   dplyr::mutate(gender = as.character(gender)) %>% 
   dplyr::rename(subject_id = "sample_id") 



##### COmbining all the controls
Control_1_t <- Control_1 %>% 
   dplyr::mutate(study_condition = as.character(study_condition)) %>% 
   dplyr::mutate(age = as.numeric(age)) %>% 
   dplyr::mutate(gender = as.character(gender)) 


Control_3_t <-Control_3 %>% 
   dplyr::mutate(subject_id = as.character(subject_id)) %>% 
   dplyr::mutate(study_condition = as.character(study_condition)) %>% 
   dplyr::mutate(age = as.numeric(age)) 


### Keeping the uniqque samples from each study , focusing on study 3 which as meta-study contains information from the other 2

All_Control_Cancer <- 
   dplyr::bind_rows(Control_3_t,
                    Control_1_t
   ) %>% 
   dplyr::mutate(Group = c(rep("Control 3",nrow(Control_3_t)),
                           rep("Control 1",nrow(Control_1_t)))) %>% 
   dplyr::distinct(subject_id,
                   study_condition,
                   gender,
                   country,
                   .keep_all = T) %>%    
   dplyr::mutate(continent = ifelse(country == "AUS"|country == "GER"|country == "FRA"| country == "ITA","Europe",
                                    ifelse(country == "CAN"|country == "USA","America","Αsia")) %>% as.factor()) %>% 
   dplyr::select(-c(Group,study_condition))



All_Controls_Countries <- dplyr::bind_rows(All_Control_Cancer,
                                           Controls_Lu_Countries) %>% 
   dplyr::mutate(Group = c(rep("Caner_Controls",nrow(All_Control_Cancer)),
                           rep("Lu_Countries",nrow(Controls_Lu_Countries)))) %>% 
   dplyr::distinct(subject_id,
                   gender,
                   country,
                   .keep_all = T) %>% 
   dplyr::mutate(gender = ifelse(gender == "male",1,2) %>% as.factor()) %>% 
   dplyr::mutate(real_age = age) %>% 
   dplyr::mutate(age = ifelse(age <= 50,0,1)) 


All_Controls_Continents <- dplyr::bind_rows(All_Control_Cancer,
                                            Controls_Lu_Continents) %>% 
   dplyr::mutate(Group = c(rep("Caner_Controls",nrow(All_Control_Cancer)),
                           rep("Lu_Continents",nrow(Controls_Lu_Continents)))) %>% 
   dplyr::distinct(subject_id,
                   gender,
                   country,
                   .keep_all = T) %>% 
   dplyr::mutate(gender = ifelse(gender == "male",1,2) %>% as.factor()) %>% 
   dplyr::mutate(real_age = age) %>% 
   dplyr::mutate(age = ifelse(age <= 50,0,1)) %>%
   dplyr::mutate(country = as.factor(country))

All_Controls_Continents_under_50 <- All_Controls_Continents %>% dplyr::filter(real_age <= 50) %>% dplyr::select(-c(Group))
All_Controls_Continents_above_50 <- All_Controls_Continents %>% dplyr::filter(real_age > 50) %>% dplyr::select(-c(Group))

Cancer_Continents_under_50 <- matched_3_res %>%dplyr::select(age.y,gender.y,BMI.y,country.y,continent.y,real_age.y,subject_id.y) %>% 
   dplyr::rename(age = "age.y",BMI = "BMI.y", country = "country.y", gender = "gender.y",continent = "continent.y", real_age = "real_age.y",
                 subject_id = "subject_id.y")

Cancer_Continents_above_50 <- matched_3_res %>% dplyr::select(age.x,gender.x,BMI.x,country.x,continent.x,real_age.x,subject_id.x)%>% 
   dplyr::rename(age = "age.x",BMI = "BMI.x", country = "country.x", gender = "gender.x",continent = "continent.x", real_age = "real_age.x",
                 subject_id = "subject_id.x")

Continents_under_50 <- dplyr::bind_rows(Cancer_Continents_under_50,All_Controls_Continents_under_50) %>% 
   dplyr::mutate(Group = c(rep("CRC",nrow(Cancer_Continents_under_50)),
                           rep("Control",nrow(All_Controls_Continents_under_50))) %>% as.factor()) %>% 
   dplyr::mutate(country = as.factor(country)) %>% 
   dplyr::mutate(continent = as.factor(continent))

Continents_above_50 <- dplyr::bind_rows(Cancer_Continents_above_50,All_Controls_Continents_above_50) %>% 
   dplyr::mutate(Group = c(rep("CRC",nrow(Cancer_Continents_above_50)),
                           rep("Control",nrow(All_Controls_Continents_above_50))) %>% as.factor()) %>% 
   dplyr::mutate(country = as.factor(country))%>% 
   dplyr::mutate(continent = as.factor(continent))



############################################## CTR-CRC matching ###############################################################

set.seed(123)
CRC_Ctrl_Countries_u50 <- matchit(Group ~   real_age + BMI + country + gender, data = Countries_under_50, 
                                  replace = FALSE,
                                  distance = "mahalanobis", method = "optimal",  within = list(BMI=3,real_age = 5),
                                  ratio = 5)
CRC_Ctrl_Continents_u50 <- matchit(Group ~  real_age + BMI + continent  + gender, data = Continents_under_50, method = "optimal",
                                   replace = FALSE,
                                   distance = "mahalanobis", within = list(BMI=3,real_age = 5), ratio=5)

CRC_Ctrl_Countries_a50 <- matchit(Group ~    real_age + BMI +  country +  gender, data = Countries_above_50, method = "optimal",
                                  replace = FALSE,
                                  distance = "mahalanobis", within = list(BMI=3,real_age = 5), ratio=5)
CRC_Ctrl_Continents_a50 <- matchit(Group ~  real_age + BMI +  continent +  gender, data = Continents_above_50, method = "optimal",
                                   replace = FALSE,
                                   distance = "mahalanobis", within = list(BMI=3,real_age = 5), ratio=5)

summary(CRC_Ctrl_Continents_u50)
summary(CRC_Ctrl_Continents_a50)


the_best_match_CRC_country <- function(matched,dataset) {
   
   # dataset = Countries_under_50
   # matched = CRC_Ctrl_Countries_u50
   matched_C <-  match.data(matched) %>% dplyr::filter(Group == "CRC") %>% as.data.frame()
   matched_T <-  match.data(matched) %>% dplyr::filter(Group == "Control") %>% as.data.frame()
   
   
   true_matched_T = data.frame()
   for (i in 1:nrow(matched_C)) {
      print(i)
      current_matched <-
         matched_T %>% 
         dplyr::filter(subclass == i) %>% 
         filter(gender == matched_C %>% dplyr::filter(subclass == i) %>% pull(gender)) %>% 
         filter(country == matched_C %>% dplyr::filter(subclass == i) %>% pull(country) %>% as.character()) %>% 
         filter(BMI == DescTools::Closest(BMI, matched_C %>% dplyr::filter(subclass == i) %>% pull(BMI))) %>% 
         as.data.frame() %>% 
         .[1,]
      
      
      true_matched_T <- rbind(true_matched_T,
                              current_matched)
   }
   
   final_data <- true_matched_T %>% 
      left_join(matched_C, by = "subclass") %>% 
      dplyr::select(subclass,
                    age.x,age.y,BMI.x,BMI.y,country.x,country.y,
                    continent.x,continent.y,
                    gender.x,gender.y,
                    real_age.x,real_age.y,
                    subject_id.x,subject_id.y) 
   
   dim(final_data)
   unique(final_data$subject_id.x)
   unique(final_data$subject_id.y)
   
   Correct_Countries_1 <-
      final_data %>% 
      filter(country.x == country.y) %>% 
      filter(abs(BMI.x-BMI.y) <=3)  %>% 
      filter(abs(real_age.x-real_age.y) <=5) 
   
   #### Second Round
   
   set.seed(123)
   matched_2 <- matchit(Group ~   real_age + BMI + country + gender, 
                                     data = dataset %>% 
                                        filter(subject_id %!in% Correct_Countries_1$subject_id.y) %>% 
                                        filter(subject_id %!in% Correct_Countries_1$subject_id.x), 
                                     distance = "mahalanobis", method = "optimal",  within = list(BMI=3,real_age = 5),
                                     ratio = 5)
   
   matched_C_2 <-  match.data(matched_2) %>% dplyr::filter(Group == "CRC") %>% as.data.frame()
   matched_T_2 <-  match.data(matched_2) %>% dplyr::filter(Group == "Control") %>% as.data.frame()
   
   
   
   true_matched_T_2 = data.frame()
   for (i in 1:nrow(matched_C_2)) {
      print(i)
      current_matched <-
         matched_T_2 %>% 
         dplyr::filter(subclass == i) %>% 
         filter(gender == matched_C_2 %>% dplyr::filter(subclass == i) %>% pull(gender)) %>% 
         filter(country == matched_C_2 %>% dplyr::filter(subclass == i) %>% pull(country) %>% as.character()) %>% 
         filter(BMI == DescTools::Closest(BMI, matched_C_2 %>% dplyr::filter(subclass == i) %>% pull(BMI))) %>% 
         as.data.frame() %>% 
         .[1,]
      
      
      true_matched_T_2 <- rbind(true_matched_T_2,
                              current_matched)
   }
   
   final_data_2 <- true_matched_T_2 %>% 
      left_join(matched_C_2, by = "subclass") %>% 
      dplyr::select(subclass,
                    age.x,age.y,BMI.x,BMI.y,country.x,country.y,
                    continent.x,continent.y,
                    gender.x,gender.y,
                    real_age.x,real_age.y,
                    subject_id.x,subject_id.y) 
   
   
   Correct_Countries_2 <-
      final_data_2 %>% 
      filter(country.x == country.y) %>% 
      filter(abs(BMI.x-BMI.y) <=3)  %>% 
      filter(abs(real_age.x-real_age.y) <=5) 
   
   
   
   ### Third Round
   
   
   set.seed(123)
   matched_3 <- matchit(Group ~   real_age + BMI + country + gender, 
                        data = dataset %>% 
                           filter(subject_id %!in% Correct_Countries_1$subject_id.y) %>% 
                           filter(subject_id %!in% Correct_Countries_1$subject_id.x) %>% 
                           filter(subject_id %!in% Correct_Countries_2$subject_id.y) %>% 
                           filter(subject_id %!in% Correct_Countries_2$subject_id.x), 
                        distance = "mahalanobis", method = "optimal",  within = list(BMI=3,real_age = 5),
                        ratio = 5)
   
   matched_C_3 <-  match.data(matched_3) %>% dplyr::filter(Group == "CRC") %>% as.data.frame()
   matched_T_3 <-  match.data(matched_3) %>% dplyr::filter(Group == "Control") %>% as.data.frame()
   
   
   
   true_matched_T_3 = data.frame()
   for (i in 1:nrow(matched_C_3)) {
      print(i)
      current_matched <-
         matched_T_3 %>% 
         dplyr::filter(subclass == i) %>% 
         filter(gender == matched_C_3 %>% dplyr::filter(subclass == i) %>% pull(gender)) %>% 
         filter(country == matched_C_3 %>% dplyr::filter(subclass == i) %>% pull(country) %>% as.character()) %>% 
         filter(BMI == DescTools::Closest(BMI, matched_C_3 %>% dplyr::filter(subclass == i) %>% pull(BMI))) %>% 
         as.data.frame() %>% 
         .[1,]
      
      
      true_matched_T_3 <- rbind(true_matched_T_3,
                                current_matched)
   }
   
   final_data_3 <- true_matched_T_3 %>% 
      left_join(matched_C_3, by = "subclass") %>% 
      dplyr::select(subclass,
                    age.x,age.y,BMI.x,BMI.y,country.x,country.y,
                    continent.x,continent.y,
                    gender.x,gender.y,
                    real_age.x,real_age.y,
                    subject_id.x,subject_id.y) 
   
   
   Correct_Countries_3 <-
      final_data_3 %>% 
      filter(country.x == country.y) %>% 
      filter(abs(BMI.x-BMI.y) <=3)  %>% 
      filter(abs(real_age.x-real_age.y) <=5) 

   
   
   
   
   ##### 4th round
   
   set.seed(123)
   matched_4 <- matchit(Group ~   real_age + BMI + country + gender, 
                        data = dataset %>% 
                           filter(subject_id %!in% Correct_Countries_1$subject_id.y) %>% 
                           filter(subject_id %!in% Correct_Countries_1$subject_id.x) %>% 
                           filter(subject_id %!in% Correct_Countries_2$subject_id.y) %>% 
                           filter(subject_id %!in% Correct_Countries_2$subject_id.x) %>% 
                           filter(subject_id %!in% Correct_Countries_3$subject_id.y) %>% 
                           filter(subject_id %!in% Correct_Countries_3$subject_id.x), 
                        distance = "mahalanobis", method = "optimal",  within = list(BMI=3,real_age = 5),
                        ratio = 5)
   
   matched_C_4 <-  match.data(matched_4) %>% dplyr::filter(Group == "CRC") %>% as.data.frame()
   matched_T_4 <-  match.data(matched_4) %>% dplyr::filter(Group == "Control") %>% as.data.frame()
   
   
   
   true_matched_T_4 = data.frame()
   for (i in 1:nrow(matched_C_4)) {
      print(i)
      current_matched <-
         matched_T_4 %>% 
         dplyr::filter(subclass == i) %>% 
         filter(gender == matched_C_4 %>% dplyr::filter(subclass == i) %>% pull(gender)) %>% 
         filter(country == matched_C_4 %>% dplyr::filter(subclass == i) %>% pull(country) %>% as.character()) %>% 
         filter(BMI == DescTools::Closest(BMI, matched_C_4 %>% dplyr::filter(subclass == i) %>% pull(BMI))) %>% 
         as.data.frame() %>% 
         .[1,]
      
      
      true_matched_T_4 <- rbind(true_matched_T_4,
                                current_matched)
   }
   
   final_data_4 <- true_matched_T_4 %>% 
      left_join(matched_C_4, by = "subclass") %>% 
      dplyr::select(subclass,
                    age.x,age.y,BMI.x,BMI.y,country.x,country.y,
                    continent.x,continent.y,
                    gender.x,gender.y,
                    real_age.x,real_age.y,
                    subject_id.x,subject_id.y) 
   
   
   Correct_Countries_4 <-
      final_data_4 %>% 
      filter(country.x == country.y) %>% 
      filter(abs(BMI.x-BMI.y) <=3)  %>% 
      filter(abs(real_age.x-real_age.y) <=5) 
   
   
 
  
### Fifth round
   set.seed(123)
   matched_5 <- matchit(Group ~   real_age + BMI + country + gender, 
                        data = dataset %>% 
                          filter(subject_id %!in% Correct_Countries_1$subject_id.y) %>% 
                          filter(subject_id %!in% Correct_Countries_1$subject_id.x), 
                        distance = "mahalanobis", method = "optimal",  within = list(BMI=3,real_age = 5),
                        ratio = 5)
   
   matched_C_5 <-  match.data(matched_5) %>% dplyr::filter(Group == "CRC") %>% as.data.frame()
   matched_T_5 <-  match.data(matched_5) %>% dplyr::filter(Group == "Control") %>% as.data.frame()
   
   
   
   true_matched_T_5 = data.frame()
   for (i in 1:nrow(matched_C_5)) {
     print(i)
     current_matched <-
       matched_T_5 %>% 
       dplyr::filter(subclass == i) %>% 
       filter(gender == matched_C_5 %>% dplyr::filter(subclass == i) %>% pull(gender)) %>% 
       filter(country == matched_C_5 %>% dplyr::filter(subclass == i) %>% pull(country) %>% as.character()) %>% 
       filter(BMI == DescTools::Closest(BMI, matched_C_5 %>% dplyr::filter(subclass == i) %>% pull(BMI))) %>% 
       as.data.frame() %>% 
       .[1,]
     
     
     true_matched_T_5 <- rbind(true_matched_T_5,
                               current_matched)
   }
   
   final_data_5 <- true_matched_T_5 %>% 
     left_join(matched_C_5, by = "subclass") %>% 
     dplyr::select(subclass,
                   age.x,age.y,BMI.x,BMI.y,country.x,country.y,
                   continent.x,continent.y,
                   gender.x,gender.y,
                   real_age.x,real_age.y,
                   subject_id.x,subject_id.y) 
   
   
   Correct_Countries_5 <-
     final_data_5 %>% 
     filter(country.x == country.y) %>% 
     filter(abs(BMI.x-BMI.y) <=3)  %>% 
     filter(abs(real_age.x-real_age.y) <=5)
  
 

##################################################### continent #################################################################

the_best_match_CRC_continent <- function(matched,dataset) {
  
  # dataset = Continents_under_50
  # matched = CRC_Ctrl_Continents_u50
  matched_C <-  match.data(matched) %>% dplyr::filter(Group == "CRC") %>% as.data.frame()
  matched_T <-  match.data(matched) %>% dplyr::filter(Group == "Control") %>% as.data.frame()
  
  
  
  
  
  true_matched_T = data.frame()
  for (i in 1:nrow(matched_C)) {
    print(i)
    current_matched <-
      matched_T %>% 
      dplyr::filter(subclass == i) %>% 
      filter(gender == matched_C %>% dplyr::filter(subclass == i) %>% pull(gender)) %>% 
      filter(continent == matched_C %>% dplyr::filter(subclass == i) %>% pull(continent) %>% as.character()) %>% 
      filter(BMI == DescTools::Closest(BMI, matched_C %>% dplyr::filter(subclass == i) %>% pull(BMI))) %>% 
      as.data.frame() %>% 
       .[1,]
    
    
    true_matched_T <- rbind(true_matched_T,
                            current_matched)
  }
  
  final_data <- true_matched_T %>% 
    left_join(matched_C, by = "subclass") %>% 
    dplyr::select(subclass,
                  age.x,age.y,BMI.x,BMI.y,country.x,country.y,
                  continent.x,continent.y,
                  gender.x,gender.y,
                  real_age.x,real_age.y,
                  subject_id.x,subject_id.y) 
  
  Correct_Continent_1 <-
     final_data %>% 
     filter(continent.x == continent.y) %>% 
     filter(abs(BMI.x-BMI.y) <=3)  %>% 
     filter(abs(real_age.x-real_age.y) <=5) 
  
  #### Second Round
  
  set.seed(123)
  matched_2 <- matchit(Group ~   real_age + BMI + continent + gender, 
                       data = dataset %>% 
                          filter(subject_id %!in% Correct_Continent_1$subject_id.y) %>% 
                          filter(subject_id %!in% Correct_Continent_1$subject_id.x), 
                       distance = "mahalanobis", method = "optimal",  within = list(BMI=3,real_age = 5),
                       ratio = 5)
  
  matched_C_2 <-  match.data(matched_2) %>% dplyr::filter(Group == "CRC") %>% as.data.frame()
  matched_T_2 <-  match.data(matched_2) %>% dplyr::filter(Group == "Control") %>% as.data.frame()
  
  
  
  true_matched_T_2 = data.frame()
  for (i in 1:nrow(matched_C_2)) {
     print(i)
     current_matched <-
        matched_T_2 %>% 
        dplyr::filter(subclass == i) %>% 
        filter(gender == matched_C_2 %>% dplyr::filter(subclass == i) %>% pull(gender)) %>% 
        filter(continent == matched_C_2 %>% dplyr::filter(subclass == i) %>% pull(continent) %>% as.character()) %>% 
        filter(BMI == DescTools::Closest(BMI, matched_C_2 %>% dplyr::filter(subclass == i) %>% pull(BMI))) %>% 
        as.data.frame() %>% 
        .[1,]
     
     
     true_matched_T_2 <- rbind(true_matched_T_2,
                               current_matched)
  }
  
  final_data_2 <- true_matched_T_2 %>% 
     left_join(matched_C_2, by = "subclass") %>% 
     dplyr::select(subclass,
                   age.x,age.y,BMI.x,BMI.y,country.x,country.y,
                   continent.x,continent.y,
                   gender.x,gender.y,
                   real_age.x,real_age.y,
                   subject_id.x,subject_id.y) 
  
  
  Correct_Continent_2 <-
     final_data_2 %>% 
     filter(continent.x == continent.y) %>% 
     filter(abs(BMI.x-BMI.y) <=3)  %>% 
     filter(abs(real_age.x-real_age.y) <=5) 
  
  
  
  ### Third Round
  
  
  set.seed(123)
  matched_3 <- matchit(Group ~   real_age + BMI + continent + gender, 
                       data = dataset %>% 
                          filter(subject_id %!in% Correct_Continent_1$subject_id.y) %>% 
                          filter(subject_id %!in% Correct_Continent_1$subject_id.x) %>% 
                          filter(subject_id %!in% Correct_Continent_2$subject_id.y) %>% 
                          filter(subject_id %!in% Correct_Continent_2$subject_id.x), 
                       distance = "mahalanobis", method = "optimal",  within = list(BMI=3,real_age = 5),
                       ratio = 5)
  
  matched_C_3 <-  match.data(matched_3) %>% dplyr::filter(Group == "CRC") %>% as.data.frame()
  matched_T_3 <-  match.data(matched_3) %>% dplyr::filter(Group == "Control") %>% as.data.frame()
  
  
  
  true_matched_T_3 = data.frame()
  for (i in 1:nrow(matched_C_3)) {
     print(i)
     current_matched <-
        matched_T_3 %>% 
        dplyr::filter(subclass == i) %>% 
        filter(gender == matched_C_3 %>% dplyr::filter(subclass == i) %>% pull(gender)) %>% 
        filter(continent == matched_C_3 %>% dplyr::filter(subclass == i) %>% pull(continent) %>% as.character()) %>% 
        filter(BMI == DescTools::Closest(BMI, matched_C_3 %>% dplyr::filter(subclass == i) %>% pull(BMI))) %>% 
        as.data.frame() %>% 
        .[1,]
     
     
     true_matched_T_3 <- rbind(true_matched_T_3,
                               current_matched)
  }
  
  final_data_3 <- true_matched_T_3 %>% 
     left_join(matched_C_3, by = "subclass") %>% 
     dplyr::select(subclass,
                   age.x,age.y,BMI.x,BMI.y,country.x,country.y,
                   continent.x,continent.y,
                   gender.x,gender.y,
                   real_age.x,real_age.y,
                   subject_id.x,subject_id.y) 
  
  
  Correct_Continent_3 <-
     final_data_3 %>% 
     filter(continent.x == continent.y) %>% 
     filter(abs(BMI.x-BMI.y) <=3)  %>% 
     filter(abs(real_age.x-real_age.y) <=5) 
  
  not_correct <- final_data_3 %>% 
     filter(continent.x != continent.y | abs(BMI.x-BMI.y) >3 | abs(real_age.x-real_age.y) > 5) 
  
  #Forth round
  set.seed(123)
  matched_4 <- matchit(Group ~   real_age + BMI + continent + gender, 
                       data = dataset %>% 
                         filter(subject_id %!in% Correct_Continent_1$subject_id.y) %>% 
                         filter(subject_id %!in% Correct_Continent_1$subject_id.x) %>% 
                         filter(subject_id %!in% Correct_Continent_2$subject_id.y) %>% 
                         filter(subject_id %!in% Correct_Continent_2$subject_id.x), 
                       distance = "mahalanobis", method = "optimal",  within = list(BMI=3,real_age = 5),
                       ratio = 5)
  
  matched_C_4 <-  match.data(matched_4) %>% dplyr::filter(Group == "CRC") %>% as.data.frame()
  matched_T_4 <-  match.data(matched_4) %>% dplyr::filter(Group == "Control") %>% as.data.frame()
  
  
  
  true_matched_T_4 = data.frame()
  for (i in 1:nrow(matched_C_4)) {
    print(i)
    current_matched <-
      matched_T_4 %>% 
      dplyr::filter(subclass == i) %>% 
      filter(gender == matched_C_4 %>% dplyr::filter(subclass == i) %>% pull(gender)) %>% 
      filter(continent == matched_C_4 %>% dplyr::filter(subclass == i) %>% pull(continent) %>% as.character()) %>% 
      filter(BMI == DescTools::Closest(BMI, matched_C_4 %>% dplyr::filter(subclass == i) %>% pull(BMI))) %>% 
      as.data.frame() %>% 
      .[1,]
    
    
    true_matched_T_4 <- rbind(true_matched_T_4,
                              current_matched)
  }
  
  final_data_4 <- true_matched_T_4 %>% 
    left_join(matched_C_4, by = "subclass") %>% 
    dplyr::select(subclass,
                  age.x,age.y,BMI.x,BMI.y,country.x,country.y,
                  continent.x,continent.y,
                  gender.x,gender.y,
                  real_age.x,real_age.y,
                  subject_id.x,subject_id.y) 
  
  
  Correct_Continent_4 <-
    final_data_4 %>% 
    filter(continent.x == continent.y) %>% 
    filter(abs(BMI.x-BMI.y) <=3)  %>% 
    filter(abs(real_age.x-real_age.y) <=5) 
  
  not_correct <- final_data_4 %>% 
    filter(continent.x != continent.y | abs(BMI.x-BMI.y) >3 | abs(real_age.x-real_age.y) > 5)
  
  #Fifth round
  set.seed(123)
  matched_5 <- matchit(Group ~   real_age + BMI + continent + gender, 
                       data = dataset %>% 
                         filter(subject_id %!in% Correct_Continent_1$subject_id.y) %>% 
                         filter(subject_id %!in% Correct_Continent_1$subject_id.x) %>% 
                         filter(subject_id %!in% Correct_Continent_2$subject_id.y) %>% 
                         filter(subject_id %!in% Correct_Continent_2$subject_id.x), 
                       distance = "mahalanobis", method = "optimal",  within = list(BMI=3,real_age = 5),
                       ratio = 5)
  
  matched_C_5 <-  match.data(matched_5) %>% dplyr::filter(Group == "CRC") %>% as.data.frame()
  matched_T_5 <-  match.data(matched_5) %>% dplyr::filter(Group == "Control") %>% as.data.frame()
  
  
  
  true_matched_T_5 = data.frame()
  for (i in 1:nrow(matched_C_5)) {
    print(i)
    current_matched <-
      matched_T_5 %>% 
      dplyr::filter(subclass == i) %>% 
      filter(gender == matched_C_5 %>% dplyr::filter(subclass == i) %>% pull(gender)) %>% 
      filter(continent == matched_C_5 %>% dplyr::filter(subclass == i) %>% pull(continent) %>% as.character()) %>% 
      filter(BMI == DescTools::Closest(BMI, matched_C_5 %>% dplyr::filter(subclass == i) %>% pull(BMI))) %>% 
      as.data.frame() %>% 
      .[1,]
    
    
    true_matched_T_5 <- rbind(true_matched_T_5,
                              current_matched)
  }
  
  final_data_5 <- true_matched_T_5 %>% 
    left_join(matched_C_5, by = "subclass") %>% 
    dplyr::select(subclass,
                  age.x,age.y,BMI.x,BMI.y,country.x,country.y,
                  continent.x,continent.y,
                  gender.x,gender.y,
                  real_age.x,real_age.y,
                  subject_id.x,subject_id.y) 
  
  
  Correct_Continent_5 <-
    final_data_5 %>% 
    filter(continent.x == continent.y) %>% 
    filter(abs(BMI.x-BMI.y) <=3)  %>% 
    filter(abs(real_age.x-real_age.y) <=5) 
  
  not_correct <- final_data_5 %>% 
    filter(continent.x != continent.y | abs(BMI.x-BMI.y) >3 | abs(real_age.x-real_age.y) > 5)
  
   
  all_correct <-rbind(Correct_Continent_1,
                      Correct_Continent_2,
                      Correct_Continent_3,
                      Correct_Continent_4,
                      Correct_Continent_5)
  
  List = list(Matched_data = all_correct,
              Not_Correct = not_correct)
  

  return(List)
  
}


matched_CRC_Ctrl_Continents_u50 <- 
  the_best_match_CRC_continent(CRC_Ctrl_Continents_u50,Continents_under_50) 

matched_CRC_Ctrl_Continents_a50 <- 
  the_best_match_CRC_continent(CRC_Ctrl_Continents_a50,Continents_above_50) 

matched_CRC_Ctrl_Continents_u50$Matched_data %>% dim()  #Not_Correct
matched_CRC_Ctrl_Continents_a50$Matched_data %>% dim()


data.frame(c(matched_CRC_Ctrl_Countries_u50$Matched_data$subject_id.y,
             matched_CRC_Ctrl_Countries_u50$Not_Correct$subject_id.y))

data.frame(c(matched_CRC_Ctrl_Countries_a50$Matched_data$subject_id.y,
             matched_CRC_Ctrl_Countries_a50$Not_Correct$subject_id.y))

write.table(data.frame(c(matched_CRC_Ctrl_Continents_a50$Matched_data)),
            col.names = T,
            row.names = F,
            quote = FALSE,
            "~/ct_matched_CRC_CTR_a50.txt")

