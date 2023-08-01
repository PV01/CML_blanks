##### PV change request ###
##### http://www.columbia.edu/~sg3637/blog/Time_Series_Heatmaps.html


library(tidyverse)
library(odbc)
library(readxl)
library(stringr)
library(makeR)
library(ggplot2)


query_db<-function(query){
  file_path<-"G:/COMPUTER/Shire/read_only_user/ForReports/db1.mdb"    
  conn <- dbConnect(drv = odbc(), .connection_string = paste0("Driver={Microsoft Access Driver (*.mdb, *.accdb)};DBQ=",file_path,";"))
  res<-dbGetQuery(conn, query, stringsAsFactors = FALSE)
  dbDisconnect(conn)
  return(res)
}


######### RT-PCR worksheet info
get_wks_info<-function(){
  
  q_wks_info<-paste("SELECT DNA_Worksheet_det_result.WORKSHEET, DNA_Worksheet_det_result.TEST, DNA_Worksheet_det_result.LANE, DNA_Worksheet_det.LABNO, DNA_Worksheet.WORKSHEET_DATE, 
                             DNA_Worksheet_det_result.RESULT_DATE, DNA_EXTRACTS.EXTRACT_SHEET, DNA_EXTRACTS.EXTRACTION_METHOD, DNA_EXTRACTS.CONDITION, DNA_Extractsheet.EXTRACTION_DATE
                       FROM (DNA_Worksheet_det_result INNER JOIN ((DNA_Worksheet_det INNER JOIN DNA_EXTRACTS ON DNA_Worksheet_det.LABNO = DNA_EXTRACTS.LABNO) 
                             INNER JOIN DNA_Extractsheet ON DNA_EXTRACTS.EXTRACT_SHEET = DNA_Extractsheet.EXTRACTSHEET) ON (DNA_Worksheet_det_result.WORKSHEET = DNA_Worksheet_det.WORKSHEET) 
                             AND (DNA_Worksheet_det_result.LANE = DNA_Worksheet_det.LANE)) INNER JOIN DNA_Worksheet ON DNA_Worksheet_det_result.WORKSHEET = DNA_Worksheet.WORKSHEET
                      WHERE (((DNA_Worksheet_det_result.TEST)='BCR-ABL1 single step') AND ((DNA_Worksheet.WORKSHEET_DATE)>#1/1/2023#) AND ((DNA_EXTRACTS.EXTRACTION_METHOD)='cDNA Prep'))
                      ORDER BY DNA_Worksheet_det_result.WORKSHEET;")
  

  
  wks_info<-query_db(q_wks_info) %>%
                   mutate(WORKSHEET_WEEKDAY=weekdays(WORKSHEET_DATE)) %>%
                   relocate(WORKSHEET_WEEKDAY, .after=WORKSHEET_DATE)
}

wks_info<-get_wks_info()

wks_info %>%
  distinct(WORKSHEET, .keep_all=TRUE) %>%
        count(WORKSHEET_WEEKDAY)

#WORKSHEET_WEEKDAY  n
#1            Monday  9
#2          Thursday  4
#3           Tuesday 15
#4         Wednesday  3


wks_info %>%
        distinct(EXTRACTION_DATE, .keep_all=TRUE) %>%
                     count(WORKSHEET, name="cDNA_prep_dates") %>%
                     count(cDNA_prep_dates, name="Worksheets") %>%
                     select(Worksheets,cDNA_prep_dates)

#Worksheets cDNA_prep_dates
#1          2               1
#2          4               2
#3         12               3
#4          7               4
#5          5               5
#6          1               7



sample_nr<-wks_info %>%
                  distinct(WORKSHEET, LABNO, .keep_all = TRUE) %>%
                  count(WORKSHEET, name="Nr_Samples") 

hist(sample_nr$Nr_Samples)

my_hist <- hist(sample_nr$Nr_Samples)                     # Store histogram info
my_hist$counts <- cumsum(my_hist$counts)    # Change histogram counts
plot(my_hist) 




ggplot(sample_nr, aes(X=Nr_Samples))+
  geom_histogram()


#####################################################################################################
wks_info_heat<-wks_info %>%
                  mutate(WORKSHEET_DATE=as.Date(WORKSHEET_DATE)) %>%
                  select(WORKSHEET_DATE, LABNO) %>%
                  distinct() %>%
                  count(WORKSHEET_DATE)
r2g <- c("#D61818", "#FFAE63", "#FFFFBD", "#B5E384")
g2r <- c("#B5E384", "#FFFFBD", "#FFAE63", "#D61818")
  
calendarHeat(wks_info_heat$WORKSHEET_DATE, wks_info_heat$n,  ncolors = 99, color = "g2r", varname="Nr. Samples/Worksheet")
#############################################################################################


###################################
get_cDNA_wks_info<-function(){
  
  q_cDNA_wks_info<-paste("SELECT DNA_EXTRACTS.EXTRACT_SHEET, DNA_EXTRACTS.EXTRACTION_METHOD, DNA_EXTRACTS.TUBE, DNA_EXTRACTS.LABNO, DNALAB_INDICATION.INDICATION, DNALAB_INDICATION.REASON, DNA_Extractsheet.EXTRACTION_DATE
                          FROM DNALAB_INDICATION INNER JOIN (DNA_EXTRACTS INNER JOIN DNA_Extractsheet ON DNA_EXTRACTS.EXTRACT_SHEET = DNA_Extractsheet.EXTRACTSHEET) ON DNALAB_INDICATION.LABNO = DNA_EXTRACTS.LABNO
                          WHERE (((DNA_EXTRACTS.EXTRACTION_METHOD)='cDNA Prep') AND ((DNA_Extractsheet.EXTRACTION_DATE)>#3/1/2023#))
                         ORDER BY DNA_EXTRACTS.EXTRACT_SHEET;")
  
  
  
  cDNA_wks_info<-query_db(q_cDNA_wks_info)
  
}

cDNA_wks_info<-get_cDNA_wks_info()

cDNA_wks_info %>% 
             distinct(EXTRACT_SHEET, LABNO, TUBE, .keep_all = TRUE) %>%
             count(EXTRACT_SHEET, EXTRACTION_DATE) %>%
             filter(n>31)


#EXTRACT_SHEET EXTRACTION_DATE  n
#1        115387      2023-03-06 32
#2        115491      2023-03-08 32
#3        115650      2023-03-13 32
#4        117404      2023-05-02 32
#5        117642      2023-05-09 32
#6        117862      2023-05-15 32
#7        118114      2023-05-22 32
#8        119471      2023-06-26 32
#9        119732      2023-07-03 32


cDNA_wks_info %>% 
  distinct(EXTRACT_SHEET, LABNO, TUBE, .keep_all = TRUE) %>%
  count(EXTRACT_SHEET, EXTRACTION_DATE, name="Nr_cDNAs") %>%
  mutate(WEEKDAY=weekdays(EXTRACTION_DATE)) %>%
  select(WEEKDAY, Nr_cDNAs)



info<-cDNA_wks_info %>% 
           distinct(EXTRACT_SHEET, LABNO, TUBE, .keep_all = TRUE) %>%
           count(EXTRACT_SHEET, EXTRACTION_DATE, INDICATION, name="Nr_cDNAs_RMPD") %>%
           group_by(EXTRACT_SHEET, EXTRACTION_DATE) %>%
           mutate(Nr_cDNAs_TOTAL=sum(Nr_cDNAs_RMPD)) %>%
           ungroup() %>%
           filter(INDICATION=="R-MPD") %>%
           mutate(WEEKDAY=weekdays(EXTRACTION_DATE)) %>%
           select(WEEKDAY, Nr_cDNAs_RMPD, Nr_cDNAs_TOTAL) %>%
           mutate(WEEKDAY=factor(WEEKDAY, levels=c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday")))


ggplot(info, aes(x=WEEKDAY, y=Nr_cDNAs_TOTAL))+
  geom_boxplot(fill="darkgreen")

ggplot(info, aes(x=WEEKDAY, y=Nr_cDNAs_TOTAL))+
  ggeom_boxplot(aes(x=WEEKDAY, y=Nr_cDNAs_RMPD), fill="orange")


####################################################
get_cDNA_prelist<-function(){
  
  q_cDNA_prelist<-paste("SELECT DNA_EXTRACTS.EXTRACT_SHEET, DNA_EXTRACTS.LABNO, DNALAB_INDICATION.INDICATION, DNA_EXTRACTS.TUBE, DNA_EXTRACTS.EXTRACTION_METHOD, DNA_EXTRACTS.UPDATES, 
                                DNA_EXTRACTS.UPDATEDDATE, DNA_EXTRACTS.UPDATEDTIME, DNA_Extractsheet.EXTRACTION_DATE
                         FROM (DNA_Extractsheet INNER JOIN DNA_EXTRACTS ON DNA_Extractsheet.EXTRACTSHEET = DNA_EXTRACTS.EXTRACT_SHEET) INNER JOIN DNALAB_INDICATION ON DNA_EXTRACTS.LABNO = DNALAB_INDICATION.LABNO
                         WHERE (((DNA_EXTRACTS.EXTRACTION_METHOD)='cDNA prep') AND ((DNA_Extractsheet.EXTRACTION_DATE)>#2/1/2023#));")
  
  
  
  cDNA_prelist<-query_db(q_cDNA_prelist) %>%
                  mutate(EXTRACTION_DATE=ymd(EXTRACTION_DATE),
                         cDNA__Req_date=dmy(str_extract(UPDATES, "[0-9]{2}/[0-9]{2}/[0-9]{4}")),
                         cDNA__Req_Time=str_extract(UPDATES, "[0-9]{2}:[0-9]{2}")) %>%
                  select(-UPDATES, -UPDATEDDATE, -UPDATEDTIME) %>%
                  distinct(.keep_all = TRUE)
  
  
  
}


cDNA_prelist<-get_cDNA_prelist()


cDNA_prelist %>% 
           mutate(HR=as.numeric(str_extract(cDNA__Req_Time,"[0-9]{2}"))) %>%
           filter(cDNA__Req_date!=EXTRACTION_DATE, HR<14, INDICATION=="R-MPD")

#EXTRACT_SHEET     LABNO INDICATION TUBE EXTRACTION_METHOD EXTRACTION_DATE cDNA__Req_date cDNA__Req_Time HR
#2        116424 D23.11778      R-MPD    2         cDNA Prep      2023-04-03     2023-04-01          12:50 12

#"116377" "116424"


#EXTRACT_SHEET     LABNO INDICATION TUBE EXTRACTION_METHOD EXTRACTION_DATE cDNA__Req_date cDNA__Req_Time
#1         116377 D22.79021      R-MPD    2         cDNA Prep      2023-03-31     2023-03-30          17:00
#2         116377 D23.15651      R-MPD    2         cDNA Prep      2023-03-31     2023-03-30          16:38
#3         116377 D23.21852      R-MPD    2         cDNA Prep      2023-03-31     2023-03-31          12:55
#4         116377 D23.21866      R-MPD    2         cDNA Prep      2023-03-31     2023-03-31          10:48
#5         116377 D23.22452      R-ALL    2         cDNA Prep      2023-03-31     2023-03-30          18:04
#6         116377 D23.22627      R-BCR    2         cDNA Prep      2023-03-31     2023-03-31          08:55
#7         116377 D23.22721      R-BCR    2         cDNA Prep      2023-03-31     2023-03-31          08:57
#8         116377 D23.22727      R-BCR    2         cDNA Prep      2023-03-31     2023-03-31          09:00
#9         116377 D23.22733      R-BCR    2         cDNA Prep      2023-03-31     2023-03-31          09:02
#10        116377 D23.22738      R-BCR    2         cDNA Prep      2023-03-31     2023-03-31          09:05
#11        116377 D23.22740      R-BCR    2         cDNA Prep      2023-03-31     2023-03-31          08:33
#12        116377 D23.22744      R-BCR    2         cDNA Prep      2023-03-31     2023-03-31          08:43
#13        116377 D23.22803      R-BCR    2         cDNA Prep      2023-03-31     2023-03-31          08:45
#14        116377 D23.22846      R-BCR    2         cDNA Prep      2023-03-31     2023-03-31          08:47
#15        116377 D23.22853      R-BCR    2         cDNA Prep      2023-03-31     2023-03-31          08:53
#16        116377 D23.22891      R-BCR    2         cDNA Prep      2023-03-31     2023-03-31          08:55


#par(mfrow) and layout()
#https://bookdown.org/ndphillips/YaRrr/arranging-plots-with-parmfrow-and-layout.html
