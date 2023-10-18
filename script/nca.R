rm(list=ls()); gc();  gc();
if (!require("renv")  ) {install.packages("renv")}
if (!require("pacman")) {install.packages("pacman")}

# CRAN から入手可能なパッケージ
##############################
pacman::p_load(
  # 一般的なデータ管理
  ####################
  tidyverse, 
  magrittr, 
  
  # パッケージのインストールと管理
  ################################
  pacman,   # パッケージのインストール・読み込み
  renv,     # グループで作業する際のパッケージのバージョン管理  
  
  # プロジェクトとファイルの管理
  ##############################
  here,     # Rのプロジェクトフォルダを基準とするファイルパス
  rio,      # 様々なタイプのデータのインポート・エクスポート
  
  # 臨床薬理領域系の解析パッケージ
  ################################
  NonCompart, # NCA処理するバッケージ
  ncar,       # NonCompartの拡張版
  
  # CDISC ADaM関連パッケージ
  ##########################
  Tplyr,     # Rのプロジェクトフォルダを基準とするファイルパス

  # スタイルテーブル関連パッケージ
  ################################
  huxtable,  # HTML, LaTeX, RTF, 'Word', 'Excel', and 'PowerPoint'へ変換可能なスタイル
  
  # 図表関連パッケージ
  ####################
  patchwork, # 複数の図表をまとめられるパッケージ
  
  # 出力形式関連パッケージ
  ########################
  pharmaRTF  # 医薬品申請関連資料の出力用パッケージ
  )

###############################################################################################
#setwd("temp")
path   <- here::here()
adsl   <- import(paste(path , "adam", "adsl.sas7bdat" , sep = "/"))
adpc   <- import(paste(path , "adam", "adpc.sas7bdat" , sep = "/"))
adpp   <- import(paste(path , "adam", "adpp.sas7bdat" , sep = "/"))

###############################################################################################

adpc2 <- adpc %>%
  filter(  PARAMCD=="THEOPHS" 
         & TRTA %in% c("Theophylline_C01","Theophylline_C02") ) %>%
  mutate(ARELTM2 = if_else(ARELTM < 0, 0 ,ARELTM))

## NCA 320mg
nca_c01 <- adpc2 %>% 
  filter(TRTA=="Theophylline_C01") %>%
  tblNCA(  .
         , key     = "SUBJID"
         , colTime = "ARELTM2"
         , colConc = "AVAL"
         , dose    = 320
         , adm     = "Extravascular"
         , dur     = 0ß
         , doseUnit = "mg"
         , timeUnit = "h"
         , concUnit = "mg/L"
         , down     = "Linear")

## NCA 320mg
nca_c02 <- adpc2 %>% 
  filter(TRTA=="Theophylline_C02") %>%
  tblNCA(  .
         , key     = "SUBJID"
         , colTime = "ARELTM2"
         , colConc = "AVAL"
         , dose    = 640
         , adm     = "Extravascular"
         , dur     = 0
         , doseUnit = "mg"
         , timeUnit = "h"
         , concUnit = "mg/L"
         , down     = "Linear")

nca <- rbind(nca_c01, nca_c02)

nca <- nca %>% 
  mutate(TRTA = if_else(SUBJID %in% c("A01","A02","A03","A04","A05","A06"),"C01","C02"))

nca_t  <- nca %>% 
  pivot_longer(-c(SUBJID,TRTA), names_to = "PARAMCD", values_to = "AVAL") 

prec_data <- tibble::tribble(
  ~PARAMCD, ~max_int, ~max_dec,
  "CMAX"    ,   2, 1,
  "AUCLST"  ,   4, 1,
  "AUCIFO"  ,   4, 1,
  "TMAX"    ,   2, 2,
  "MRTEVIFO",   3, 1,
  "LAMZHL"  ,   2, 2,
  ) %>%
  mutate(PARAMCD = factor(PARAMCD,c("CMAX","AUCLST","AUCIFO","TMAX","LAMZHL","MRTEVIFO")))


header_data <- adsl %>%
  filter(SAFFL == "Y" & TRT01AN %in% c(1,2)) %>%
  mutate(TRTA = factor(if_else(TRT01AN == 1,"C01","C02"),c("C01","C02"))) %>%
  group_by(TRTA) %>%
  summarise(n=n()) 

nca_summary <- nca_t %>%
  filter(PARAMCD %in% c("CMAX","AUCLST","AUCIFO","TMAX","MRTEVIFO","LAMZHL")) %>%
  mutate(PARAMCD = factor(PARAMCD,c("CMAX","AUCLST","AUCIFO","TMAX","LAMZHL","MRTEVIFO"))) %>%
  tplyr_table(.,TRTA) %>%
    add_layer(
      group_desc(AVAL, by = PARAMCD) %>% 
        set_custom_summaries(
          CV = (mean(.var) / sd(.var)) * 100
         ,geometric_mean = exp(sum(log(.var[.var > 0]), na.rm=TRUE) / length(.var))
         ,CV_geo_mean = (sqrt(exp(var(log(.var[.var > 0])-1)))) * 100
        ) %>%
        set_format_strings(
           'N'            = f_str('xx'   , n)
          ,'Mean (SD)'    = f_str('a.a+1 (a.a+2)', mean, sd)
          ,'CV% mean'     = f_str('a.a+1', CV)
          ,'Geo-mean'     = f_str('a.a+1', geometric_mean)
          ,'CV% Geo-mean' = f_str('a.a+1', CV_geo_mean)
          ,'Median'       = f_str('a.a+1', median)
          ,'[Min; Max]'   = f_str('[a.a+0; a.a+0]', min, max)
        ) %>% 
        set_precision_on(AVAL) %>% 
        set_precision_by(PARAMCD) %>%
        set_precision_data(prec_data)
    ) %>%
    build()

nca_summary2 <- nca_summary %>%
    mutate(row_label1 = case_when(
      row_label1 == "CMAX" ~ "Cmax"
     ,row_label1 == "AUCLST" ~ "AUC last"
     ,row_label1 == "AUCIFO" ~ "AUC Inf"
     ,row_label1 == "TMAX"~ "tmax"
     ,row_label1 == "LAMZHL"~ "t1/2"
     ,row_label1 == "MRTEVIFO" ~ "MRT"
     ,TRUE ~ ""
    )) %>%
    apply_row_masks(row_breaks = TRUE) %>%   
    select(-starts_with("ord")) %>% 
    add_column_headers(
      paste0(" | | Theophylline(320mg)\\line(N=**C01**)| Theophylline(640mg)\\line(N=**C02**) "), 
      header_n = header_data) 

ht <- nca_summary2 %>% 
  huxtable::as_hux(., add_colnames=FALSE) %>%
  huxtable::set_bold(1, 1:ncol(.), TRUE) %>% # bold the first row
  huxtable::set_align(1, 1:ncol(.), 'center') %>% # Center align the first row 
  huxtable::set_align(2:nrow(.), 3:ncol(.), 'center') %>% # Center align the results
  huxtable::set_valign(1, 1:ncol(.), 'bottom') %>% # Bottom align the first row
  huxtable::set_bottom_border(1, 1:ncol(.), 1) %>% # Put a border under the first row
  huxtable::set_width(1.5) %>% # Set the table width
  huxtable::set_escape_contents(FALSE) %>% # Don't escape RTF syntax
  huxtable::set_col_width(c(.2, .2, .2, .2)) # Set the column widths
ht

doc <- pharmaRTF::rtf_doc(ht) %>% 
  pharmaRTF::add_titles(
    pharmaRTF::hf_line("Protocol: CDISCPILOT01", "PAGE_FORMAT: Page %s of %s", align='split', bold=TRUE, italic=TRUE),
    pharmaRTF::hf_line("Table 14-2.01", align='center', bold=TRUE, italic=TRUE),
    pharmaRTF::hf_line("Summary of Demographic and Baseline Characteristics", align='center', bold=TRUE, italic=TRUE)
  ) %>% 
  pharmaRTF::add_footnotes(
    pharmaRTF::hf_line("FILE_PATH: Source: %s", "DATE_FORMAT: %H:%M %A, %B %d, %Y", align='split', bold=FALSE, italic=TRUE)
  ) %>% 
  pharmaRTF::set_font_size(10) %>%
  pharmaRTF::set_ignore_cell_padding(TRUE) %>% 
  pharmaRTF::set_column_header_buffer(top=1)

pharmaRTF::write_rtf(doc, file='styled_example.rtf')

################################################################################################

p1 <- adpc2 %>% 
  filter(TRTA=="Theophylline_C01") %>%
  ggplot(.,aes(x=ARELTM2,y=AVAL,group=SUBJID))+
  theme_set(theme_classic()) +
  geom_line(aes(linetype = SUBJID))+
  ggtitle("Linier view") +
  xlab("Time (h)")+
  ylab("Concentration (mg/L)")+
  theme_bw()+ 
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5)) 

p2<- adpc2 %>% 
  filter(TRTA=="Theophylline_C01") %>%
  ggplot(.,aes(x=ARELTM2,y=AVAL,group=SUBJID))+ 
  theme_set(theme_classic()) +
  geom_line(aes(linetype = SUBJID))+
  geom_abline(intercept = log10(0.5), slope = 0,linetype = 2) +
  annotate("text", x=22, y=0.4, label="BLQ")+
  ggtitle("Semilogarithmic view") +
  xlab("Time (h)")+
  scale_y_continuous(trans='log10')+
  ylab("Concentration (mg/L)")+
  theme_bw()+ 
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5))

p1 + p2

p3 <- adpc2 %>% 
  filter(TRTA=="Theophylline_C02") %>%
  ggplot(.,aes(x=ARELTM2,y=AVAL,group=SUBJID))+
  theme_set(theme_classic()) +
  geom_line(aes(linetype = SUBJID))+
  ggtitle("Linier view") +
  xlab("Time (h)")+
  ylab("Concentration (mg/L)")+
  theme_bw()+ 
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5)) 

p4 <- adpc2 %>% 
  filter(TRTA=="Theophylline_C02") %>%
  ggplot(.,aes(x=ARELTM2,y=AVAL,group=SUBJID))+ 
  theme_set(theme_classic()) +
  geom_line(aes(linetype = SUBJID))+
  geom_abline(intercept = log10(0.5), slope = 0,linetype = 2) +
  annotate("text", x=22, y=0.4, label="BLQ")+
  ggtitle("Semilogarithmic view") +
  xlab("Time (h)")+
  scale_y_continuous(trans='log10')+
  ylab("Concentration (mg/L)")+
  theme_bw()+ 
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5)) 

p3 + p4