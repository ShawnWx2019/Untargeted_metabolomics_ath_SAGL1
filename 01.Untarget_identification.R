###################################################
#         Prj: SAGL1
#         Assignment: Indentification and classification
#         Author: Shawn Wang
#         Date: 2022.1.20
#         Locate: HENU                            
###################################################


# import packages and functions -------------------------------------------


library(tidyverse)
library(PCAtools)
library(webchem)
library(classyfireR)
library(readxl)
library(furrr)
library(rvest)
options(future.rng.onMisuse = "ignore")

# my fun ------------------------------------------------------------------

## 01.get cid via compound names
shawn_get_cid = function(name) {
  Sys.sleep(0.2)
  a = tibble(
    query = name,
    cid = NA
  )
  tryCatch({
    a = get_cid(
      query = name,
      from = "name",
      domain = "compound",
      match = cid_mod 
    )
  },error=function(e) {
    cat(paste0(name," no match results\n"))
  }
  )
  return(a)
}

## 02.get detail via InChikey
classfireR2tbl = function(x){
  # make a black df
  if (is.null(x)) {
    tmp_x = data.frame(
      InchIkeys = NA,
      Level = NA,
      Classification = NA
    ) %>% as_tibble()
  } else if(is.null(x@meta$inchikey)) {  # make na values as blank
    tmp_x = data.frame(
      InchIkeys = NA,
      Level = NA,
      Classification = NA)
  } else {
    tmp_x = classification(x) %>% 
      mutate(InchIkeys = gsub("InChIKey=","",meta(x)$inchikey)) %>% 
      select(InchIkeys,Level,Classification)
  }
  return(tmp_x)
}

# get kegg via InChikey
InChikey2kegg = function(url){
  Sys.sleep(0.5)
  a = read_html(url) %>%
    html_nodes("p") %>% html_text() %>%
    gsub(pattern = '\"',replacement = "",x = .) %>%
    str_split(.,pattern = ",",n = 4,simplify = T) %>%
    as.data.frame() %>%
    mutate(
      InChIkey = gsub("searchTerm:","",V3),
      KEGG = str_extract(V4,pattern = "C\\d+")
    ) %>%
    select(InChIkey,KEGG)
  return(a)
}

classfireR2tbl = function(x){
  # make a black df
  if (is.null(x)) {
    tmp_x = data.frame(
      InchIkeys = NA,
      Level = NA,
      Classification = NA
    ) %>% as_tibble()
  } else if(is.null(x@meta$inchikey)) {  # make na values as blank
    tmp_x = data.frame(
      InchIkeys = NA,
      Level = NA,
      Classification = NA)
  } else {
    tmp_x = classification(x) %>% 
      mutate(InchIkeys = gsub("InChIKey=","",meta(x)$inchikey)) %>% 
      select(InchIkeys,Level,Classification)
  }
  return(tmp_x)
}

InChikey2kegg = function(url){
  Sys.sleep(0.5)
  a = read_html(url) %>%
    html_nodes("p") %>% html_text() %>%
    gsub(pattern = '\"',replacement = "",x = .) %>%
    str_split(.,pattern = ",",n = 4,simplify = T) %>%
    as.data.frame() %>%
    mutate(
      InChIkey = gsub("searchTerm:","",V3),
      KEGG = str_extract(V4,pattern = "C\\d+")
    ) %>%
    select(InChIkey,KEGG)
  return(a)
}
 # Identification and classification -------------------------------------------------------------------

file_input = read_xlsx("../02.data/SAGL1_metadata.xlsx",sheet = 1,col_names = T)
head(file_input)


data_info = data.frame(
  Compound_ID = file_input$Compound_ID,
  name = file_input$Name,
  mf = gsub(" ","",file_input$Formula),
  mw = file_input$`Calc  MW`,
  RT = file_input$`RT [min]`,
  CD_judge = file_input$Judge
)

Name_clean = data_info %>% 
  select(name) %>% 
  distinct() %>% 
  drop_na()

## 01. get cid and InChiKey via webchem
future::plan(multisession)

## cidmode

cid_mod = "all"
cid_tbl = future_map_dfr(Name_clean,shawn_get_cid,.progress = T) %>% 
  filter(!is.na(cid))
head(cid_tbl)
pc_prop_progress = function(cid) {
  a =  webchem::pc_prop(cid = cid,properties = c("MolecularFormula", "MolecularWeight","InChIKey", 
                                                 "IUPACName","ExactMass"))
  return(a)
}

tmp_cid_char = unique(cid_tbl$cid)
tmp_punchem_item = future_map_dfr(tmp_cid_char,pc_prop_progress,.progress = T) %>% 
  mutate(CID = as.character(CID))

# classfire  --------------------------------------------------------------

InchIKeys = tmp_punchem_item$InChIKey

shawn_get_classification = function(id){
  Sys.sleep(0.1)
  get_classification(inchi_key = id)
}

classification_list = map(InchIKeys,shawn_get_classification)



# convert as classfireR2tbl
result_tbl = purrr::map_df(classification_list,classfireR2tbl)
# convert long table to wide table
result_final = bind_rows(result_tbl) %>% 
  drop_na() %>% unique() %>% 
  pivot_wider(names_from = Level,values_from = Classification) 

result_out = cid_tbl %>% 
  left_join(.,tmp_punchem_item,by = c("cid" = "CID") )%>% 
  left_join(.,result_final,by = c("InChIKey"="InchIkeys")) %>% 
  right_join(.,data_info,by = c("query" = "name") ) %>% 
  mutate(
    MolecularWeight = as.numeric(MolecularWeight),
    mw = as.numeric(mw),
    Check_MF = case_when(
      mf == MolecularFormula ~ TRUE,
      TRUE ~ FALSE
    ),
    Check_MW = case_when(
      abs(mw - MolecularWeight) <= 5 ~  TRUE,
      TRUE ~ FALSE
    ),
    High_identical = case_when(
      Check_MF == TRUE & Check_MW == TRUE ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>% 
  distinct() %>% 
  relocate(Compound_ID,.before = query) %>% 
  as.data.frame()


# kegg annotation ---------------------------------------------------------

head(result_out)

CTS_url = "https://cts.fiehnlab.ucdavis.edu/rest/convert/InChIKey/KEGG/"
Inchikey = result_out %>%
  mutate(url = paste0(CTS_url,InChIKey)) %>%
  filter(InChIKey != "NA") %>% 
  select(url) %>% unique()
Inchikey = Inchikey$url
KEGGID = future_map_dfr(Inchikey,tryCatch({InChikey2kegg},error = function(e) { print("failed")}))
result_out_kegg <- left_join(result_out,KEGGID,by = c("InChIKey" = "InChIkey"))

head(result_out_kegg)
# fileout -----------------------------------------------------------------

file_out = "../03.progress/03.untarget/Step1_clean"

write.csv(result_out_kegg,file = paste0(file_out,".csv"),row.names = F)
