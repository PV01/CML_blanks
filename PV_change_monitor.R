##### PV change request ###
##### http://www.columbia.edu/~sg3637/blog/Time_Series_Heatmaps.html


library(tidyverse)
library(odbc)




query_db<-function(query){
  file_path<-"G:/COMPUTER/Shire/read_only_user/ForReports/db1.mdb"    
  conn <- dbConnect(drv = odbc(), .connection_string = paste0("Driver={Microsoft Access Driver (*.mdb, *.accdb)};DBQ=",file_path,";"))
  res<-dbGetQuery(conn, query, stringsAsFactors = FALSE)
  dbDisconnect(conn)
  return(res)
}

get_status<-function(){

  q_status<-paste("SELECT DNALAB_INDICATION.LABNO, DNALAB_INDICATION.INDICATION, DNALAB_INDICATION.REASON, DNALAB_INDICATION.DATE_ACTIVATED, DNALAB_INDICATION.DATEREQUESTED, 
                          DNALAB_TEST.TEST, DNALAB_TEST.WORKSHEET, DNALAB_TEST.LANE, DNALAB_TEST.RESULT, DNALAB_TEST.RETEST, DNA_EXTRACTS.EXTRACT_SHEET, DNA_EXTRACTS.EXTRACTION_METHOD, 
                          DNA_EXTRACTS.TUBE, DNA_EXTRACTS.CONDITION
                   FROM (DNALAB_INDICATION INNER JOIN DNALAB_TEST ON (DNALAB_INDICATION.REASON = DNALAB_TEST.REASON) AND (DNALAB_INDICATION.INDICATION = DNALAB_TEST.INDICATION) AND 
                        (DNALAB_INDICATION.LABNO = DNALAB_TEST.LABNO)) LEFT JOIN DNA_EXTRACTS ON DNALAB_TEST.LABNO = DNA_EXTRACTS.LABNO
                  WHERE (((DNALAB_INDICATION.INDICATION) In ('R-MPD','R-BCR')) AND ((DNALAB_INDICATION.REASON) In ('Diagnosis','Follow-up')) AND 
                        ((DNALAB_INDICATION.DATEREQUESTED)>Date()-10));")
 status<-query_db(q_status) 
  
}


cDNAs_to_prep<-get_status() %>% 
                           filter(EXTRACTION_METHOD=="cDNA Prep", is.na(CONDITION)) %>%
                           count() %>%
                           pull()



RMPD_cDNAs_to_prep<-get_status() %>%
                        filter(TEST=="BCR-ABL1 single step", EXTRACTION_METHOD=="cDNA Prep", is.na(CONDITION))%>%
                        count() %>%
                        pull()
cDNAs_to_prep
RMPD_cDNAs_to_prep


  
  
  
  