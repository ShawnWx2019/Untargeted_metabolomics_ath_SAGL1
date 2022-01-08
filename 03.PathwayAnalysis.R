###################################################
#         Prj: SAGL1
#         Assignment: Ath compound database construction
#         Author: Shawn Wang
#         Date: 2022.1.20
#         Locate: HENU                            
###################################################

## construction of ath kegg compound database


library(readxl)
library(tidyverse)
library(rvest)
sag_kegg <- read_xlsx("../03.progress/03.untarget/Step1_Summary.xlsx",sheet = 3,col_names = T)
head(sag_kegg)

clean = sag_kegg %>% 
  select(query,KEGG) %>% 
  distinct()


# get ath background by KEGG API ------------------------------------------

## 01.Get kegg ath pathway information
read_html("http://rest.kegg.jp/link/pathway/ath") %>% 
  html_text() %>% 
  read.table(text = .,sep = "\t") %>% 
  mutate(PathwayID = gsub("path:ath","map",V2)) %>% 
  select(PathwayID) %>% 
  unique() -> ath.db

## 02. map to compound id（CID）
get_cpd_map <- function(pathway){
  url = "http://rest.kegg.jp/link/cpd/"
  url_query = paste0(url,pathway)
  a = data.frame(
    PathwayID = pathway,
    KEGGID = NA
  )
  tryCatch({
    a = read_html(url_query) %>% 
      html_text() %>% 
      read.table(text = .,sep = "\t") %>% 
      mutate(
        PathwayID = gsub("path:","",V1),
        KEGGID = gsub("cpd:","",V2)
      ) %>% 
      select(PathwayID,KEGGID)
    return(a)
  }, error = function(e) {
    print(paste0(pathway,"no results"))
  }
  )
  return(a)
}

## run 
library(furrr)
future::plan(multisession)
ath.db_all = future_map_dfr(unique(ath.db$PathwayID),get_cpd_map,.progress = T)

## 03. get pathway ID to pathway name
pathway2map = read_html("http://rest.kegg.jp/list/pathway/ath") %>% 
  html_text() %>% 
  read.table(text = .,sep = "\t") %>% 
  mutate(
    PathwayID = gsub("path:ath","map",V1),
    Pathway = gsub(pattern = " - Arabidopsis thaliana \\(thale cress\\)",replacement = "",V2)
  ) %>% 
  select(PathwayID,Pathway)
## 04.construction ath kegg compound DB
ath.db_final <- ath.db_all %>% drop_na() %>% 
  inner_join(.,pathway2map,by = "PathwayID") %>% unique()

save(maize.db_final,file = "~/02.MyScript/PCK_Db_github/02.Rdata/ath.db")



# KEGG pathway analysis ---------------------------------------------------------

# for 
clean %>% 
  inner_join(.,ath.db_final,by = c("KEGG"="KEGGID")) %>% unique() -> final_result
write.table(final_result,file = "../04.result/untarget/kegg_annotation.table.xls",sep = "\t",row.names = F)



