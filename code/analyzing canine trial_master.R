library(readxl)
library(stringr)
library(tidyr)
library(splitstackshape)
library(plyr)
library(dplyr)
library(RISmed)

#load
main_data <- read_excel("C:/Users/tanyu/Box/SPARK/Comparative trials/Pubmed search/main data analysis copy.xlsx", sheet = 1)
parent <- read_excel("C:/Users/tanyu/Box/SPARK/Comparative trials/Pubmed search/main data analysis copy.xlsx", sheet = 2)

#analyze inclusion exclusion

  #remove rows where excluded is NA
  parent <- parent[complete.cases(parent$Excluded),]
  
  #count evaluated studies
  nrow(parent)
  
  #count accepted studies
  nrow(parent[parent$Excluded == 0,])
  
  #count removed at full text
  nrow(parent[grepl("F", parent$Excluded),])
  
  #count removed at abstract
  nrow(parent[grepl("A", parent$Excluded),])
  
  #reasons
  table(parent$Excluded)
  
#generate tumor type subdata set
  tumor_subdata <- main_data[,1:2]
  
  #separate by /
  colnames(tumor_subdata)[2] <- "tumor_location"
  tumor_subdata<- separate_rows(tumor_subdata, "tumor_location", sep = "/")
  
  #separate by : for tumor type
  tumor_subdata <- separate(tumor_subdata, "tumor_location", c("tumor_type", "tumor_location"), sep = ":")
  
  #separate by , for tumor location
  tumor_subdata <- separate_rows(tumor_subdata, "tumor_location", sep = ",")
  
  #separate by ( for sample size
  tumor_subdata <- separate(tumor_subdata, "tumor_location", into = c("tumor_location", "sample_size"), sep = "\\(|\\)")
  tumor_subdata$sample_size <- as.numeric(tumor_subdata$sample_size)
  
  #identify rows where type or location is not reported
  tumor_subdata$not_reported = 0
  tumor_subdata$not_reported[grepl("reported", tumor_subdata$tumor_type)|grepl("reported", tumor_subdata$tumor_location)] = 1
  
  #summary of how many tumor types per study
  tumor_type_count <- data.frame(table(tumor_subdata$`Pubmed UID`))
  hist(tumor_type_count$Freq)
  
  #annotate not reported in tumor count
  tumor_type_count$Freq[tumor_type_count$Var1 %in% 
                          tumor_subdata$`Pubmed UID`[tumor_subdata$not_reported ==1]] <- NA
  
  #annotate prior use 
  parent$EntrezUID <- as.character(parent$EntrezUID)
  tumor_type_count$Var1 <- as.character(tumor_type_count$Var1)
  tumor_type_count <- left_join(tumor_type_count, parent[,c(10,13)], by=c("Var1" = "EntrezUID"))
  
  tumor_type_count$priorhuman <-"0"
  tumor_type_count$priorhuman[grepl("3", tumor_type_count$`Prior Use`)] <- "1"
  
  #annotate prior use on specific type
  main_data$`Pubmed UID`<- as.character(main_data$`Pubmed UID`)
  tumor_type_count <- left_join(tumor_type_count, main_data[,c(1,16)], by=c("Var1" = "Pubmed UID"))
  
  tumor_type_count$human_unspec <- "0"
  tumor_type_count$human_unspec[tumor_type_count$priorhuman == 1 &
                                  tumor_type_count$`Rationale specific tumor type` == "N"] <- "1"
  
  #separate by tumor type
  carc <- unique(tumor_subdata[tumor_subdata$tumor_type=="Carcinoma",c(1,2)])
  sarc <- unique(tumor_subdata[tumor_subdata$tumor_type=="Sarcoma",c(1,2)])
  lym <- unique(tumor_subdata[tumor_subdata$tumor_type=="Lymphoma",c(1,2)])
  osteo <- unique(tumor_subdata[tumor_subdata$tumor_type=="Osteosarcoma",c(1,2)])
  mast <- unique(tumor_subdata[tumor_subdata$tumor_type=="Mast cell tumor",c(1,2)])
  mela <- unique(tumor_subdata[tumor_subdata$tumor_type=="Melanoma",c(1,2)])
  adeno <- unique(tumor_subdata[tumor_subdata$tumor_type=="Adenocarcinoma",c(1,2)])
  
  main_data.carc <- main_data[main_data$`Pubmed UID` %in% carc$`Pubmed UID`,]
  main_data.sarc <- main_data[main_data$`Pubmed UID` %in% sarc$`Pubmed UID`,]
  main_data.lym <- main_data[main_data$`Pubmed UID` %in% lym$`Pubmed UID`,]
  main_data.osteo <- main_data[main_data$`Pubmed UID` %in% osteo$`Pubmed UID`,]
  main_data.mast <- main_data[main_data$`Pubmed UID` %in% mast$`Pubmed UID`,]
  main_data.mela <- main_data[main_data$`Pubmed UID` %in% mela$`Pubmed UID`,]
  main_data.adeno <- main_data[main_data$`Pubmed UID` %in% adeno$`Pubmed UID`,]
  
  #prior use
  main_data.lym <- merge(x = main_data.lym, y = parent[,c(10,13)], by.x = "Pubmed UID", by.y= "EntrezUID")
  
  
  #export
  write.csv(tumor_subdata, "C:/Users/tanyu/Box/SPARK/Comparative trials/Pubmed search/analysis output/tumortypelong.csv")
  write.csv(tumor_type_count, "C:/Users/tanyu/Box/SPARK/Comparative trials/Pubmed search/analysis output/tumortypecount.csv")
  
  write.csv(main_data.carc, "C:/Users/tanyu/Box/SPARK/Comparative trials/Pubmed search/analysis output/main_datacarc.csv")
  write.csv(main_data.sarc, "C:/Users/tanyu/Box/SPARK/Comparative trials/Pubmed search/analysis output/main_datasarc.csv")
  write.csv(main_data.lym, "C:/Users/tanyu/Box/SPARK/Comparative trials/Pubmed search/analysis output/main_datalym.csv")
  write.csv(main_data.osteo, "C:/Users/tanyu/Box/SPARK/Comparative trials/Pubmed search/analysis output/main_dataosteo.csv")
  write.csv(main_data.mast, "C:/Users/tanyu/Box/SPARK/Comparative trials/Pubmed search/analysis output/main_datamast.csv")
  write.csv(main_data.mela, "C:/Users/tanyu/Box/SPARK/Comparative trials/Pubmed search/analysis output/main_datamela.csv")
  write.csv(main_data.adeno, "C:/Users/tanyu/Box/SPARK/Comparative trials/Pubmed search/analysis output/main_dataadeno.csv")
  
  #summary of how many studies study individual tumor types
  sum(tumor_type_count$Freq == 1, na.rm = TRUE)/nrow(tumor_type_count)
  
  #summary of sample size per tumor type
  plot(ecdf(tumor_subdata$sample_size))
  
#generate breed sub dataset
  breed_subdata <- main_data[,c(1,4)]
  
  #separate by , for breed
  breed_subdata <- separate_rows(breed_subdata, "Breed", sep = ",")
  
  #separate by () for sample size
  breed_subdata <- separate(breed_subdata, "Breed", into = c("Breed","sample_size"), sep = "\\(|\\)")
  breed_subdata$sample_size <- as.numeric(breed_subdata$sample_size)
  
  #identify rows where breed is not reported, listed as pure or other. First set all to lowercase
  breed_subdata$Breed <- tolower(breed_subdata$Breed)
  breed_subdata$not_reported = 0
  breed_subdata$not_reported[grepl("reported", breed_subdata$Breed)|grepl("pure", breed_subdata$Breed)|grepl("other", breed_subdata$Breed)] = 1
  
  #summary of how many breeds per study
  breed_count <- data.frame(table(breed_subdata$`Pubmed UID`))
  sum(breed_count$Freq==1)/nrow(breed_count)
  hist(breed_count$Freq)
  
  #annotate studies with at least one not reported
  breed_count$not.reported <- breed_count$Var1 %in% breed_subdata$`Pubmed UID`[breed_subdata$not_reported ==1]
  
  #adjusted breed count with not reported studies in separate bin
  breed_count$adjusted.n.breed <- breed_count$Freq
  breed_count$adjusted.n.breed[breed_count$not.reported == 1] <- NA
  
  #clean up names of columns
  colnames(breed_count)[c(1,2)] <- c("UID", "n.breed")
  
  #write
  write.csv(breed_count, "C:/Users/tanyu/Box/SPARK/Comparative trials/Pubmed search/analysis output/breedcount.csv")
  write.csv(breed_subdata, "C:/Users/tanyu/Box/SPARK/Comparative trials/Pubmed search/analysis output/breed_subdata.csv")
  
  
#create expt arm subdata
  exptarm_subdata <- main_data[,c(1,5,6)]
  
  #separate by : for intervention type
  exptarm_subdata <- separate(exptarm_subdata, "Intervention", into = c("intervention", "dose"), sep = ":")
  
  #count dose levels
  exptarm_subdata[is.na(exptarm_subdata$dose),3] <- "-"
  exptarm_subdata$doselvl <- str_count(exptarm_subdata$dose, ",") + 1
  
  #count intervention levels
  exptarm_subdata$intlvl <- str_count(exptarm_subdata$intervention, ",") +1
  
  #separate interventions by ,
  intervention_split <- exptarm_subdata[,c(1:2,5)]
  intervention_split <- separate_rows(intervention_split, "intervention", sep = ",")
  
  #copy intervention rows where dose >1
  intervention_split <- expandRows(intervention_split, "doselvl")
  
  #separate dose by , and copy dose rows
  dose_split <- exptarm_subdata[,c(1,3,6)]
  dose_split <- separate_rows(dose_split, "dose", sep = ",")
  
  dose_split <- expandRows(dose_split, "intlvl")
  
  #expand sample size by ,
  n_int_split <- exptarm_subdata[,c(1,4)]
  colnames(n_int_split)[2] <- "n_intervention"
  n_int_split <- separate_rows(n_int_split, "n_intervention", sep = ",")
  
  #cbind expanded data
  exptarm_subdata_clean <- cbind(intervention_split, dose_split, n_int_split)
  
  #verify correct order
  exptarm_subdata_clean$matched[exptarm_subdata_clean[1] == exptarm_subdata_clean[3] & exptarm_subdata_clean[5] == exptarm_subdata_clean[3]] <- 1
  
  #drop excess UID rows if there are no mismatches
  if (sum(exptarm_subdata_clean$matched != 1) == 0) {exptarm_subdata_clean <- exptarm_subdata_clean[,-c(3,5,7)] }
  
  #pull out control data
  control_subdata <- main_data[,c(1,8,9)]
  control_subdata <- control_subdata[control_subdata$Control != "-",]
  
  #add in dose column, rename, reorder
  control_subdata$dose <- rep("-", nrow(control_subdata))
  colnames(control_subdata)[c(2,3)] <- c("intervention", "n_intervention")
  control_subdata <- control_subdata[,c(1,2,4,3)]
  
  #rbind control and cleaned exptarm
  exptarm_subdata_clean <- rbind(exptarm_subdata_clean, control_subdata)
  
  #sort by pubmedid 
  exptarm_subdata_clean <- exptarm_subdata_clean[order(-exptarm_subdata_clean$`Pubmed UID`),]
  
  #numeric for sample size
  exptarm_subdata_clean$n_intervention <- as.numeric(exptarm_subdata_clean$n_intervention)
  
  #cumulative distribution plot
  plot(ecdf(exptarm_subdata_clean$n_intervention))
  
  #export
  write.csv(exptarm_subdata_clean, "C:/Users/tanyu/Box/SPARK/Comparative trials/Pubmed search/analysis output/exptarm_subdata_clean.csv")
  
#create outcome subdata
  outcome_split <- main_data[,c(1,10)]
  fig_split <- main_data[,c(1,11)]
  def_split <- main_data[,c(1,12)]
  prog_split <- main_data[,c(1,13)]
  
  #diagnose outcome-fig
  outcome_split$count <- str_count(outcome_split$Outcomes, ",")
  fig_split$count <- str_count(fig_split$Figures, ",")
  
  outcome_split$matched <- outcome_split$count == fig_split$count
  outcome_split[(outcome_split$matched == FALSE), 1]
  
  #diagnose fig-defined
  def_split$count <- str_count(def_split$`Defined progressive outcomes`, ",")
  
  fig_split$matched <- fig_split$count == def_split$count
  fig_split[(fig_split$matched == FALSE), 1]
  
  #diagnose defined-prog
  prog_split$count <- str_count(prog_split$`Prog outcome criteria`, ",")
  
  def_split$matched <- def_split$count == prog_split$count
  def_split[(def_split$matched == FALSE), 1]
  
  #separate rows for all
  colnames(def_split)[2] <- "def_prog_outcome"
  colnames(prog_split)[2] <-"prog_outcome_criteria"
  
  if (sum(outcome_split$matched != TRUE) ==0) {
    outcome_split <- separate_rows(outcome_split, "Outcomes", sep = ",")
    outcome_split <- outcome_split[,c(1,4)]
  }
  
  if (sum(fig_split$matched != TRUE) ==0){
    fig_split <- separate_rows(fig_split, "Figures", sep = ",")
    fig_split <- fig_split[,c(1,4)]
  }

  if (sum(def_split$matched != TRUE) ==0)
    {def_split <- separate_rows(def_split, "def_prog_outcome", sep = ",")
    def_split<- def_split[,c(1,4)]
  } 
  
  if (sum(prog_split$matched != TRUE) ==0){
    prog_split <- separate_rows(prog_split, "prog_outcome_criteria", sep = ",")
    prog_split <- prog_split[,c(1,3)]
  }
  
  #cbind everything
  outcome_subdata <- cbind(outcome_split, fig_split, def_split, prog_split)
  
  #verify pubmed IDs match and drop unecessary columns
  outcome_subdata$matched <- (outcome_subdata[1] == outcome_subdata[3] &
                                outcome_subdata[3] == outcome_subdata[5] &
                                outcome_subdata[5] == outcome_subdata[7])
  if(sum(outcome_subdata$matched != TRUE)==0) {outcome_subdata <- outcome_subdata[,-c(3,5,7,9)]}
  
  #annotate if solid, lymphoma, or both
  outcome_subdata <- join(outcome_subdata, main_data[,c(1,2)], type="left", match="all")
  
  outcome_subdata$tumor_type <- "S"
  outcome_subdata$tumor_type[grepl("Lymphoma", outcome_subdata$`Phenotypic disease population`)] <- "L"
  outcome_subdata$tumor_type[grepl("Lymphoma", outcome_subdata$`Phenotypic disease population`) 
                             & grepl("/", outcome_subdata$`Phenotypic disease population`)] <- "B" 
  
  #annotate publication date
  outcome_subdata$pubyear <- YearPubmed(EUtilsGet(outcome_subdata$`Pubmed UID`))
  
  #pull things that look at progression free survival by lesion size
  #NEED TO RECHECK CRITERIA
  outcome_subdata$Outcomes <- tolower(outcome_subdata$Outcomes)
  outcome_subdata$prog_outcome_criteria <- tolower(outcome_subdata$prog_outcome_criteria)
  outcome_subdata$progr.def <- "0"
  outcome_subdata$progr.def[grepl("progress", outcome_subdata$Outcomes) & 
                              outcome_subdata$def_prog_outcome =="Y"] <- '1'
  
  #pull things that used vcog or human recist criteria
  outcome_subdata$vcog <- "0"
  outcome_subdata$vcog[grepl("vcog", outcome_subdata$prog_outcome_criteria)|grepl("recist", outcome_subdata$prog_outcome_criteria) &
                         outcome_subdata$progr.def == 1] <- '1'
  
  #export
  write.csv(outcome_subdata, "C:/Users/tanyu/Box/SPARK/Comparative trials/Pubmed search/analysis output/outcome_subdata.csv")

#create cochrane subdata
  cochrane_subdata <- main_data[,c(1,17:29)]
  
  write.csv(cochrane_subdata, "C://Users/tanyu/Box/SPARK/Comparative trials/Pubmed search/analysis output/cochrane_subdata.csv")
  
#create primary outcome subdata
  primary_subdata <- main_data[,c(1,30)]
  abstract_subdata <- main_data[,c(1,31)]
  colnames(primary_subdata)[2] <- "primary.outcome"
  colnames(abstract_subdata)[2] <- "abstract.outcome"
  
  primary_subdata <- separate_rows(primary_subdata, "primary.outcome", sep = ",")
  primary_subdata <- separate(primary_subdata, "primary.outcome", into = c("primary.outcome","significant"), sep = "\\(|\\)")
  
  #not reported annotate
  primary_subdata$significant[primary_subdata$primary.outcome %in% "-" & is.na(primary_subdata$significant)] <- "Not reported"
  
  #repeat for abstract
  abstract_subdata <- separate_rows(abstract_subdata, "abstract.outcome", sep = ",")
  abstract_subdata <- separate(abstract_subdata, "abstract.outcome", into = c("abstract.outcome","significant"), sep = "\\(|\\)")
  
  #not reported annotate
  abstract_subdata$significant[abstract_subdata$abstract.outcome %in% "-" & is.na(abstract_subdata$significant)] <- "Not reported"
  
  #export
  write.csv(primary_subdata, "C:/Users/tanyu/Box/SPARK/Comparative trials/Pubmed search/analysis output/primary_subdata.csv")
  write.csv(abstract_subdata, "C:/Users/tanyu/Box/SPARK/Comparative trials/Pubmed search/analysis output/abstract_subdata.csv")
  