###### CML cDNA blanks

library(tidyverse)
library(odbc)
library(readxl)
library(stringr)


query_db<-function(query){
  file_path<-"G:/COMPUTER/Shire/read_only_user/ForReports/db1.mdb"    
  conn <- dbConnect(drv = odbc(), .connection_string = paste0("Driver={Microsoft Access Driver (*.mdb, *.accdb)};DBQ=",file_path,";"))
  res<-dbGetQuery(conn, query, stringsAsFactors = FALSE)
  dbDisconnect(conn)
  return(res)
}

######### Dates e13a2/e14a2 cDNAs for RQ prepped
get_cDNA_dates<-function(){
  
  q_cDNA_dates<-paste("SELECT DNA_Worksheet_det_result.WORKSHEET, DNA_Worksheet_det_result.LANE, DNA_Worksheet_det_result.TEST, DNA_Worksheet_det.LABNO, 
                               DNA_Worksheet_det.DATE_ACTIVATED, DNA_EXTRACTS.TUBE, DNA_EXTRACTS.EXTRACT_SHEET, DNA_EXTRACTS.EXTRACTION_METHOD, DNA_EXTRACTS.CONDITION, 
                               DNA_Extractsheet.EXTRACTION_DATE
                       FROM ((DNA_Worksheet_det_result INNER JOIN DNA_Worksheet_det ON (DNA_Worksheet_det_result.LANE = DNA_Worksheet_det.LANE) 
                        AND (DNA_Worksheet_det_result.WORKSHEET = DNA_Worksheet_det.WORKSHEET)) INNER JOIN DNA_EXTRACTS ON DNA_Worksheet_det.LABNO = DNA_EXTRACTS.LABNO) 
                        INNER JOIN DNA_Extractsheet ON DNA_EXTRACTS.EXTRACT_SHEET = DNA_Extractsheet.EXTRACTSHEET
                      WHERE (((DNA_Worksheet_det_result.TEST)='BCR-ABL1 RQ-PCR') AND ((DNA_Worksheet_det.DATE_ACTIVATED)>#1/1/2023#) AND ((DNA_EXTRACTS.EXTRACTION_METHOD)='cDNA Prep'));")
  
  #SELECT DISTINCT DNA_Extractsheet.EXTRACTION_DATE
  #FROM ((DNA_Worksheet_det_result INNER JOIN DNA_Worksheet_det ON (DNA_Worksheet_det_result.LANE = DNA_Worksheet_det.LANE) AND (DNA_Worksheet_det_result.WORKSHEET = DNA_Worksheet_det.WORKSHEET)) 
  #     INNER JOIN DNA_EXTRACTS ON DNA_Worksheet_det.LABNO = DNA_EXTRACTS.LABNO) INNER JOIN DNA_Extractsheet ON DNA_EXTRACTS.EXTRACT_SHEET = DNA_Extractsheet.EXTRACTSHEET
  #WHERE (((DNA_Worksheet_det_result.TEST)="BCR-ABL1 RQ-PCR") AND ((DNA_Worksheet_det.DATE_ACTIVATED)>#1/1/2023#) AND ((DNA_EXTRACTS.EXTRACTION_METHOD)="cDNA Prep"));
  
  
 #SELECT DNA_Worksheet_det_result.WORKSHEET, DNA_Worksheet_det_result.LANE, DNA_Worksheet_det_result.TEST, DNA_Worksheet_det.LABNO, DNA_Worksheet_det.DATE_ACTIVATED, DNA_EXTRACTS.TUBE, DNA_EXTRACTS.EXTRACT_SHEET, DNA_EXTRACTS.EXTRACTION_METHOD, DNA_EXTRACTS.CONDITION, DNA_Extractsheet.EXTRACTION_DATE
 #FROM ((DNA_Worksheet_det_result INNER JOIN DNA_Worksheet_det ON (DNA_Worksheet_det_result.WORKSHEET = DNA_Worksheet_det.WORKSHEET) AND (DNA_Worksheet_det_result.LANE = DNA_Worksheet_det.LANE)) INNER JOIN DNA_EXTRACTS ON DNA_Worksheet_det.LABNO = DNA_EXTRACTS.LABNO) INNER JOIN DNA_Extractsheet ON DNA_EXTRACTS.EXTRACT_SHEET = DNA_Extractsheet.EXTRACTSHEET
 #WHERE (((DNA_Worksheet_det_result.TEST)="BCR-ABL1 single step") AND ((DNA_Worksheet_det.DATE_ACTIVATED)>#1/1/2023#) AND ((DNA_EXTRACTS.EXTRACTION_METHOD)='cDNA Prep'));
                                                                         
  cDNA_dates<-query_db(q_cDNA_dates) %>%
                   mutate(cDNA_Date=ymd(EXTRACTION_DATE))%>%
                   arrange(cDNA_Date) %>%
                   select(cDNA_Date) %>%
                   unique()
  
  return(cDNA_dates)
}


cDNA_dates<-get_cDNA_dates() %>%
                mutate(cDNA="Yes") %>%
                `rownames<-`( NULL )
###########################################################
path<-"G:/DNA/RESULTS/Taqman/Haemato-Oncology Taqman 2013 onwards/2023/Taqman analysed spreadsheets/CML"
wks_files<-file.path(path,list.files(path, recursive = TRUE))
wks_files<-wks_files[-which(grepl("Data|[$]", wks_files))]


get_dates_rqblanks<-function(path){
  dates<-read_excel(path, sheet="7500 Paste", range = "B16:B17", col_names = FALSE) %>%
                   rename(RqBl_Date="...1") %>%
                   mutate(RqBl_Date=dmy(str_extract(RqBl_Date,"[0-9]{2}.[0-9]{2}.[0-9]{2}")))
}

rq_blank_dates<-data.frame()
for (i in 1:length(wks_files)){
   rq_blank_dates<-rbind(rq_blank_dates,get_dates_rqblanks(wks_files[i]))
}

rq_cDNA_dates<-cDNA_dates%>%
                            anti_join(rq_blank_dates, by=c("cDNA_Date"="RqBl_Date")) %>%
                            rename(Date=cDNA_Date)
                           
