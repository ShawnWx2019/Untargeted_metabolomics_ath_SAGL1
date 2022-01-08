###################################################
#         Prj: SAGL1
#         Assignment: untarget metabolite
#         Author: Shawn
#         Date: 2022.1.20
#         Locate: HENU                            
###################################################


# packages ----------------------------------------------------------------


library(tidyverse)
library(ropls)
library(plotly)
library(ComplexHeatmap)
library(PCAtools)
library("broom")
library(xlsx)
library(readxl)
library(ropls)

# function ----------------------------------------------------------------


meta_heatmap = function(raw_data,meta_sample){
  scaled_data = raw_data %>% 
    column_to_rownames("CompoundID")
  anno_t = data.frame(
    row.names = meta_sample$SampleID,
    group = factor(meta_sample$Sample_group_data,levels = unique(meta_sample$Sample_group_data))
  )
  x = raw_data %>% 
    column_to_rownames("CompoundID")
  a = pheatmap::pheatmap(
    x,scale = "row",show_rownames = F,
    annotation_col = anno_t,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    color = colorRampPalette(c("blue", "white", "red"))(50)
  )
  return(a)
  
}

PCA_PLS_plt = function(logged_mat,comp_group,title,type = "pls"){
  ## pca with pca tools
  if (ncol(comp_group) != 2) {
    cat(paste0("Please check your sample group data \n",
               "It should only have two columns \n",
               "1st column is Sample ID, which is same as sample ids in compound peak-area matrix, second column \n",
               "2nd column is Sample group, which is the group you want to display in pca or pls-da plot"
    )
    )
  } else {
    colnames(comp_group) = c("SampleID","group")
    comp_group = comp_group %>% 
      column_to_rownames("SampleID")
    pca_result = pca(mat = logged_mat,metadata = comp_group,scale = T,removeVar = 0.1)
    biplot_result = biplot(
      pcaobj = pca_result,
      legendPosition = "top",
      colby = "group",
      title = "PCA biplot",
      subtitle = paste0("result of ", title),
      showLoadings = T,
      hline = 0,
      vline = 0
    )
    x = comp_group$group
    x.level = unique(x)
    g_factor = factor(x,levels = x.level)
    if (length(x.level < 3)){
      cross_val1 = 6
    } else {
      cross_val1 = 7
    }
    if (type == "pls") {
      plsDA_plot = opls(x = t(log10(logged_mat)),y = g_factor,crossvalI = cross_val1)
    } else if(type == "opls"){
      plsDA_plot = opls(x = t(log10(logged_mat)),y = g_factor,crossvalI = cross_val1,orthoI = 1)
    }

    result_all = list(
      pca_result = pca_result,
      biplot_result = biplot_result,
      plsDA_plot = plsDA_plot
    )
  }
}

maketrace = function(data,name,hovertext_in){
  trace <- list(
    mode = "markers",
    name = name,
    type = "scatter3d",
    x = as.numeric(filter(data,Type == name)$PC1),
    y = as.numeric(filter(data,Type == name)$PC2),
    z = as.numeric(filter(data,Type == name)$PC3),
    hovertext = paste(name,hovertext_in)
  )
}

PCA_3D_Interact = function(data,group,title,hovertext_in){
  pca.info = pca(mat = data,metadata = group,scale = T,removeVar = 0.1)
  
  pca.data = data.frame(
    sample = colnames(data),
    Type = group$group,
    pca.info$rotated
  )
  ## set layout 
  layout <- list(
    scene = list(
      xaxis = list(
        type = "linear", 
        title = "PC1"
      ), 
      yaxis = list(
        type = "linear", 
        title = "PC2"
      ), 
      zaxis = list(
        type = "linear", 
        title = "PC3"
      ), 
      camera = list(
        up = list(
          x = 0, 
          y = 0, 
          z = 1
        ), 
        eye = list(
          x = 0.7872575900178773, 
          y = 0.9534734809221759, 
          z = 1.5538414418639417
        ), 
        center = list(
          x = 0, 
          y = 0, 
          z = 0
        )
      ), 
      aspectmode = "auto", 
      aspectratio = list(
        x = 1.663724192127058, 
        y = 0.7845820579507525, 
        z = 0.7660908864950811
      )
    ), 
    title = title, 
    autosize = TRUE
  )
  ## make trace 
  gp_uniq = unique(group$group)
  trace_list = list()
  p <- plot_ly()
  for (i in 1:length(gp_uniq)) {
    trace_list[[i]] = maketrace(data = pca.data,name = gp_uniq[i],hovertext_in = hovertext_in)
    p <- add_trace(p,mode=trace_list[[i]]$mode, name=trace_list[[i]]$name, type=trace_list[[i]]$type, x=trace_list[[i]]$x, y=trace_list[[i]]$y, z=trace_list[[i]]$z,hovertext = trace_list[[i]]$hovertext)
  }
  
  layout(p,scene=layout$scene, title=layout$title, autosize=layout$autosize)
}

outliner_remove <- function(raw_mat,raw_sample,out_sample){
  clean_mat <- meta_mat %>% 
    select(-out_sample)
  clean_sample <- meta_sample %>% 
    filter(SampleID != out_sample)
  out_list = list(
    clean_mat = clean_mat,
    clean_sample = clean_sample
  )
}

DM_analysis = function(x,left_index,right_index,l_num,r_num,right,left){
  
  mat = x %>% 
    column_to_rownames("CompoundID") %>% 
    select(left_index,right_index)
  colnames(mat) = c(paste0("left","_",str_pad(c(1:l_num), 3,side = "left", "0")),
                    paste0("right","_",str_pad(c(1:r_num), 3,side = "left", "0")))
  mean_mat = mat %>% 
    rownames_to_column("CompoundID") %>% 
    pivot_longer(!CompoundID,names_to = "Sample",values_to = "Peak") %>% 
    mutate(Sample = gsub("....$","",Sample)) %>% 
    group_by(CompoundID,Sample) %>% 
    summarise(mean = mean(Peak)) %>% 
    pivot_wider(names_from = Sample,values_from = mean)
  # get.Test
  mat %>% 
    rownames_to_column("CompoundID") %>% 
    pivot_longer(!CompoundID,names_to = "Sample",values_to = "Peak") %>% 
    mutate(Sample = gsub("....$","",Sample)) %>% 
    group_by(CompoundID,Sample) %>% 
    nest() %>% 
    spread(key = Sample,value = data) %>% 
    mutate(
      t_test = map2(
        .x = left,
        .y = right,
        .f =  ~{t.test(.x$Peak, .y$Peak) %>% tidy()}
      ),
      log2fc = map2(
        .x = left,
        .y = right,
        .f = ~{log2(mean(.x$Peak)/mean(.y$Peak))}
      )
    ) %>% 
    unnest() %>% 
    select(-contains("Peak")) %>% 
    unique() -> Diff_result_all
  Diff_result_all$FDR = p.adjust(p = Diff_result_all$p.value,method = "BH",n = nrow(Diff_result_all))
  opls_da = opls(
    x = mat %>% t(),
    y = data.frame(
      row.names = colnames(mat),
      group = rep(c(left,right),times = c(l_num,r_num))
    ) %>% as.matrix(),
    orthoI = 1
  )
  VIP= getVipVn(opls_da) %>% as.data.frame()
  Diff_result_all$VIP = VIP$.
  Diff_result_clean = data.frame(
    CompoundID = Diff_result_all$CompoundID,
    left = mean_mat$left,
    right = mean_mat$right,
    log2fc = Diff_result_all$log2fc,
    pvalue = Diff_result_all$p.value,
    FDR = Diff_result_all$FDR,
    VIP = Diff_result_all$VIP
  ) 
  colnames(Diff_result_clean)[c(2,3)] = c(left,right)
  
  ## get vip
  return(Diff_result_clean)
}

vocPlotly = function(x,left,right){
  x %>% 
    mutate(group = case_when(
      pvalue <= 0.05 & VIP >= 1 & log2fc > 0 ~ "up",
      pvalue <= 0.05 & VIP >= 1 & log2fc < 0 ~ "down",
      TRUE ~ "not significant"
    ),
    log10pval = -log10(pvalue)
    )-> vocdata
  
  p<- plot_ly(
    data = vocdata,
    x = ~log2fc,
    y = ~log10pval,
    mode = "markers",
    color = ~ group,
    marker = list(size = ~VIP*5,opacity = 0.8),
    hovertext = vocdata$CompoundID
  )  
  layout(p,title = paste(left,"vs",right,sep = " "),
         yaxis = list(title = "-log10(p-value)"),
         xaxis = list(title = "log2(Fold change)")
  )
}

Vocanoplot = function(x,title){
  x %>% 
    mutate(group = case_when(
      pvalue <= 0.05 & VIP >= 1 & log2fc > 0 ~ "up",
      pvalue <= 0.05 & VIP >= 1 & log2fc < 0 ~ "down",
      TRUE ~ "not significant"
    ),
    log10pval = -log10(pvalue)
    )-> vocdata
  p = ggplot(data = vocdata,aes(x = log2fc,y = log10pval))+
    geom_point(aes(color = group,size = VIP))+
    scale_color_manual(values = c("blue","grey","red"))+
    scale_size(range = c(0.05,1))+
    theme_bw()+
    ylab(expression(paste(-log10,("p-value"))))+
    xlab(expression(paste(log2,"(fold change)")))+
    xlim(-max(vocdata$log2fc)-1,max(vocdata$log2fc)+1)+
    ylim(0,max(vocdata$log10pval)+1)+
    ggtitle(title)+
    annotate("text",x = max(vocdata$log2fc)-1,y = max(vocdata$log10pval)-1,
             label = paste0("up regular:",nrow(filter(vocdata,group == "up"))),
             color = "red",size = 3)+
    annotate("text",x = -max(vocdata$log2fc)+1,y = max(vocdata$log10pval)-1,
             label = paste0("down regular:",nrow(filter(vocdata,group == "down"))),
             color = "blue",size = 3)+
    geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)+
    theme_prism(border = T)
  
  return(p)
}


# data import -------------------------------------------------------------
data_meta <- read.csv("../03.progress/03.untarget/Step1_clean.csv")
head(data_meta)

data_mat <- read_xlsx("../02.data/SAGL1_metadata.xlsx",sheet = 3,col_names = T)
head(data_mat)
group <- read_xlsx("../02.data/SAGL1_metadata.xlsx",sheet = 2,col_names = T)
head(group)
# PCA ---------------------------------------------------------------------



## pick out high_conf compounds

## pick out high confidence compound with pubchem database selection.
high_conf = data_meta %>% 
  filter( High_identical == "TRUE"
  ) 



Comp_id_candidate = tibble(
  Compound_ID = high_conf$Compound_ID
) %>% 
  distinct()

high_mat <- data_mat %>% 
  inner_join(.,Comp_id_candidate,by = "Compound_ID") %>% 
  column_to_rownames("Compound_ID")

result1_pca_plsda = PCA_PLS_plt(logged_mat = high_mat,comp_group = group,title = "SAGL",type = "pls")

pdf(file = "../03.progress/03.untarget/PCA_with_loading.pdf",height = 8,width = 7.5)
result1_pca_plsda$biplot_result
dev.off()
pdf(file = "../03.progress/03.untarget/PCA.pdf",height = 8,width = 7.5)
biplot(
  pcaobj = result1_pca_plsda$pca_result,
  legendPosition = "top",pointSize = 8,lab = "",legendIconSize = 8,
  colby = "group",
  title = "PCA biplot",
  hline = 0,
  vline = 0
)
dev.off()



# DMA ---------------------------------------------------------------------

tmp_a = colnames(high_mat)

get_DM_index = function(sample_name,left,right) {
  l_index = sample_name[grepl(pattern = left,x = sample_name)]
  r_index = sample_name[grepl(pattern = right,x = sample_name)]
  index = list(
    l_index = l_index,
    r_index = r_index
  )
  return(index)
}

sagl_DM = function(left,right,...) {
  x_index = get_DM_index(
    sample_name = tmp_a,
    left = left,
    right = right
  )
  
  dm_mat = high_mat %>% rownames_to_column("CompoundID")
  outfile = DM_analysis(
    x = dm_mat,
    left_index = x_index$l_index,
    right_index = x_index$r_index,
    l_num = 4,
    r_num = 4,
    left = left,
    right = right
  )
  outfilename = paste(left,right,sep = "_vs_")
  write.csv(outfile,file = paste0("../03.progress/03.untarget/",outfilename,".csv"),row.names = F)
  return(outfile)
}

## kfb_vs_sagl1

DM_list = list()
DM_list2 = list()
sample = c("sagl","kfbQ","kfbT","WT")

comp = data.frame(
  left = rep(sample,times = c(3,2,1,0)),
  right = c(sample[-1],sample[c(3,4)],sample[4])
)

for (i in 1:nrow(comp)) {
  DM_list[[i]] =  sagl_DM(
    left = comp$left[i],
    right = comp$right[i]
  )
  DM_list2[[i]] = DM_list[[i]] 
  DM_list2[[i]]$group = paste0(comp$left[i]," vs ",comp$right[i])
  colnames(DM_list2[[i]])[c(2,3)] = c("left","right")
}

DM_tbl_merge = bind_rows(DM_list2)

DM_all = DM_tbl_merge %>% filter(pvalue <= 0.05 & VIP >= 1) %>% 
  left_join(high_conf,by = c("CompoundID" = "Compound_ID"))



write.csv(DM_all,file = "../04.result/untarget/DM_all_group2.csv",row.names = F)

svg(file = paste0("../03.progress/03.untarget/",outfilename,".svg"),height = 8,width = 8)
Vocanoplot(x = outfile,title = paste0(left," vs ",right))
dev.off()


