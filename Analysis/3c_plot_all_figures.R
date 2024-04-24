#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(GagnonMR)
library(forestploter)

setwd("/home/gagelo01/workspace/Projects/2024/LPL_pathway/")
source("Analysis/tosourcedrug.R")
ctp <- fread("/mnt/sda/gagelo01/Projects/2024/UKB_phewas/Data/Modified/clinical_toxicity_panel.txt")
res_multicis <- fread("Data/Modified/res_multicis_independent.txt")
df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
df_index[id=="trait-10-3", clean_variable_name := "Waist-to-hip ratio"]
dt_gene_region <- fread("Data/Modified/dt_gene_region.txt")

dt<-data.table(id=res_multicis[!grepl("^eqtl-|^prot-", id.exposure), unique(id.exposure)])
dt[, c("hgnc", "id_small"):=tstrsplit(id, "_", fixed = TRUE)]
dt_gene_region <- merge(dt_gene_region, dt[,.(id, id_small)], by = "id", all.x = TRUE)
dt_gene_region[is.na(id_small), id_small:=id]
dt_gene_region <- merge(dt_gene_region, df_index[, .(id,sample_size)], by.x = "id_small", by.y = "id", all.x = TRUE)

res_multicis <- merge(res_multicis, dt_gene_region[,.(id, hgnc, sample_size)], by.x = "id.exposure", by.y = "id", all.x = TRUE)
res_multicis <- merge(res_multicis, df_index[, .(id,clean_variable_name)], by.x = "id.outcome", by.y = "id", all.x = TRUE)
res_multicis[,c("hgnc", "idscale") := tstrsplit(id.exposure, "_")]
##Figure2. Association of LPL pathway with CAD, T2D. (forestplot)
arg1 <- expand_grid(ID = list(c("dis-13-1", "dis-24-4")), unit = "or")
arg2 <- expand_grid(ID = list(c("trait-1-1", "trait-10-3", "trait-19-5")), unit = "sd")
arg <- rbind(arg1, arg2)
arg <- expand_grid(arg, idscale=res_multicis$idscale %>% unique)
setDT(arg)
list_forest <- map(split(arg, 1:arg[,.N]), function(x) {
  dt_forest <- res_multicis[idscale == x$idscale, ]
  xscale <- df_index[id==x$idscale, clean_variable_name]
  dt_forest <- dt_forest[method == "Inverse variance weighted", ]
  dt_forest[,b := b*-1]
  dt_forest[, mechanism := ifelse(hgnc%in%"LPL", "activation", "inhibition")]
  dt_forest[, hgnc_mechanism := paste0(hgnc, " ", mechanism)]
  dt_forest[,lci := b-1.96*se]
  dt_forest[,uci := b+1.96*se]
  dt_forest[,panel :=stringr::str_pad("", 42, "right", " "),]
  dt_forest<-dt_forest[order(hgnc_mechanism, clean_variable_name)]
  dt<-data.table(dt_forest)
  
  
  p <- my_forestplotter_fancy(data = dt[id.outcome %in% x$ID[[1]], ], 
                              col.below.header = NULL, 
                              col.header = NULL, 
                              col.header.heading = NULL, 
                              col.left = c("hgnc_mechanism", "clean_variable_name"),
                              col.left.heading = c("Exposure", "Outcome"),
                              col.left.toformat = "hgnc_mechanism",
                              effect.name = if(x$unit=="sd"){"Beta (95% CI)"}else{"OR (95% CI)"}, 
                              col.right = NULL, 
                              exponentiate = if(x$unit=="sd"){FALSE}else{TRUE},
                              xlab = paste0("Effects scaled on 1-SD lower ",xscale))
  
  
  return(p)
})


panel_1 <- cowplot::plot_grid(list_forest[[1]], list_forest[[2]],
                              list_forest[[3]], list_forest[[4]],
                              labels = c('A)', 'B)', "C)", "D)"),
                              label_size = 14, ncol = 2)


ggsave(filename = paste0("Results/lplpathway_figure2.png"), plot = panel_1,
       width = 1100/72,height = 500/72,units="in",scale=1, device = "png")

#####Fig 3 Genetic mimicry analysis#######
# #####Concordance with LPL enhancement######
k1 <- res_multicis[idscale=="trait-16-39" & grepl("trait-33-", id.outcome)&method%in%c("Inverse variance weighted"), ]
scatter <- dcast(k1, id.outcome ~ hgnc, value.var = c("b", "se"))

colnom<-colnames(scatter)[grepl("^b_", colnames(scatter))]
scatter[,(colnom):=lapply(.SD, function(x) x*-1) , .SDcols = colnom]
cor(scatter[,.SD,.SDcols = colnom], method = c("pearson"))
# scatter <- merge(scatter, distinct(dt), by.x = "id.outcome", by.y = "nom")
data.scatter <- data.table(scatter)

plot_scatter <- function(scatter,
                         my_xlab = "Genetically proxied ANGPTL4 inhibition",
                         my_ylab = "Genetically proxied LPL enhancing",
                         title) {
  pearson_cor <- cor(scatter$b_x, scatter$b_y, use="complete.obs", method = "pearson")
  lm_res <- lm(formula = "b_x ~ b_y", data = scatter) %>%
    broom::tidy() %>%
    as.data.table()
  
  
  ggplot2::ggplot(data = scatter, ggplot2::aes(x = b_x,
                                               y = b_y, colour = value)) + ggplot2::geom_errorbarh(ggplot2::aes(xmin = b_x -
                                                                                                                  se_x, xmax = b_x + se_x),
                                                                                                   colour = "grey", height = 0) + ggplot2::geom_errorbar(ggplot2::aes(ymin = b_y -
                                                                                                                                                                        se_y, ymax = b_y + se_y),
                                                                                                                                                         colour = "grey", width = 0) + ggplot2::geom_point(size =2) +
    ggplot2::geom_abline(ggplot2::aes(intercept = lm_res[1,estimate], slope = lm_res[2, estimate]), linetype = "dashed") +
    
    ggplot2::scale_colour_manual(values = c(rgb(112, 54, 153, maxColorValue = 255),
                                            rgb(66, 94, 191, maxColorValue = 255),
                                            rgb(84,201,237, maxColorValue = 255),
                                            rgb(59,166,18, maxColorValue = 255),
                                            rgb(255,110,26, maxColorValue = 255),
                                            rgb(149,199,71, maxColorValue = 255),
                                            rgb(161,15,125, maxColorValue = 255),
                                            # rgb(249, 106, 27, maxColorValue = 255),
                                            # rgb(214,15,102, maxColorValue = 255),
                                            rgb(8, 161, 217, maxColorValue = 255),
                                            rgb(255,186,51, maxColorValue = 255),
                                            rgb(54, 150, 214, maxColorValue = 255))) +
    ggplot2::labs( x = my_xlab, y = my_ylab ) +
    ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2)) +
    expand_limits(x = 0, y = 0) +
    theme(
      legend.background = element_rect(fill = "white", linewidth = 4, colour = "white"),
      axis.ticks = element_line(colour = "black", size = 0.4),
      axis.title=element_text(size=14,face="bold"),
      # panel.grid.major = element_line(colour = "grey70", size = 0.2),
      panel.grid.minor = element_blank(),
      legend.key=element_blank(),
      panel.background = element_rect(fill = 'white', colour = 'white'),
      # panel.border = element_rect(colour = "grey70", fill=NA, size=1),
      axis.line.x.bottom = element_line(color = "black", size = 0.2),
      axis.line.y.left   = element_line(color = "black", size = 0.2),
      axis.line.y.right  = element_blank(),
      axis.text.y.right  = element_blank(),
      panel.border       = element_blank(),
      legend.title= element_blank(),
      legend.margin=margin(c(1,5,5,5)),
      plot.caption = element_text(hjust=0.5),
      plot.title = element_text(hjust = 0.5, size=14,face="bold"),
      legend.position=c(0.3,0.85)) +
    labs(caption = "[Each measures are scaled to 1-SD lower triglyceride levels]", size = 4) +
    annotate(geom="text", x=0.8, y=-0.8, label= paste0("~italic(r)==", round(pearson_cor, digits = 2)), parse=TRUE, size=7) +
    ggtitle(title)
}

data.scatter[,value:="metabolites"]
scatter <- data.table(data.scatter)
scatter <- scatter[,.(b_APOC3, se_APOC3, b_LPL, se_LPL, value)]
setnames(scatter, c("b_APOC3", "se_APOC3", "b_LPL", "se_LPL"), c("b_x", "se_x", "b_y", "se_y"))
fig6_a <- plot_scatter(scatter = scatter,
                       my_xlab = "Genetically proxied APOC3 inhibition",
                       my_ylab = "Genetically proxied LPL enhancement",
                       title = "APOC3 inhibition vs. LPL enhancement")

scatter <- data.table(data.scatter)
scatter <- scatter[,.(b_ANGPTL4, se_ANGPTL4, b_LPL, se_LPL, value)]
setnames(scatter, c("b_ANGPTL4", "se_ANGPTL4", "b_LPL", "se_LPL"), c("b_x", "se_x", "b_y", "se_y"))
fig6_b<- plot_scatter(scatter = scatter,
                      my_xlab = "Genetically proxied ANGPTL4 inhibition",
                      my_ylab = "Genetically proxied LPL enhancement",
                      title = "ANGPTL4 inhibition vs. LPL enhancement")

panel_6 <- cowplot::plot_grid(fig6_a, fig6_b, labels = c('A)', 'B)'),
                              label_size = 22, ncol = 2)
ggsave(plot = panel_6, filename = paste0("Results/lplpathway_figure3.png"),
       width = 1000/72,height = 500/72,units="in",scale=1, device = "png")


#####Figure 4 Association og LPL pathway targets with clinical toxicology panel########
dt_ctp <- res_multicis[idscale == "trait-16-39" & id.outcome %in% ctp$id, ]
xscale <- df_index[id=="trait-16-39", clean_variable_name]
dt_ctp <- dt_ctp[method == "Inverse variance weighted", ]
dt_ctp[,b := b*-1]
dt_ctp[, mechanism := ifelse(hgnc%in%"LPL", "activation", "inhibition")]
dt_ctp[, hgnc_mechanism := paste0(hgnc, " ", mechanism)]
dt_ctp[,lci := b-1.96*se]
dt_ctp[,uci := b+1.96*se]
dt_ctp[,panel := hgnc_mechanism]
dt_ctp <- merge(dt_ctp, ctp[,.(id, category, subcategory)], by.x = "id.outcome", by.y = "id", all.x = TRUE)
dt_ctp <- dt_ctp[order(category, subcategory)]

dt <- data.table(dt_ctp)
p <- my_forestplotter_fancy(data = dt, 
                            col.below.header = "clean_variable_name", 
                            col.header = "category", 
                            col.header.heading = "Outcome", 
                            col.left = NULL,
                            col.left.heading = NULL,
                            col.left.toformat = NULL,
                            effect.name = "Beta (95% CI)", 
                            col.right = NULL, 
                            exponentiate = FALSE,
                            annotate_pvalue=TRUE,
                            ntest_adjustment=dt[,.N],
                            xlab = "")


ggsave(filename = paste0("Results/lplpathway_figure4.png"), plot = p,
       width = 1200/72,height = 650/72,units="in",scale=1, device = "png")



####Figure5. Association of LPL pathway with diseases in Finngen. (phewas plot)
phewas <- res_multicis[method=="Inverse variance weighted" & grepl("dis-26-", id.outcome), ]
ID <- phewas[, id.exposure[which.max(sample_size)], by = "hgnc"]$V1
phewas <- phewas[id.exposure%in% ID, ]
info <- fread("/mnt/sda/gagelo01/Projects/Pipelines/Robust_efficient_drug_pipeline/Data/Modified/finngen_info.txt")
info <- merge(info, df_index[grepl("dis-26-", id), .(trait, clean_variable_name, note)], by = c("note"))
setnames(info, "trait", "phenocode")
phewas[, mechanism := ifelse(hgnc%in%"LPL", "activation", "inhibition")]
phewas[, hgnc_mechanism := paste0(hgnc, " ", mechanism)]
phewas  <- merge(phewas , info[,.(phenocode,category)], by.x="outcome",by.y="phenocode")
phewas[,b := b*-1]
ntest<- phewas[,length(unique(clean_variable_name))]
setnames(phewas, c("outcome", "pval", "b") ,c("phenotype", "p", "beta"))
phewas[,groupnum:=as.factor(category)]
phewas[,direction := ifelse(beta>=0, "25", "24")]
phewas[,description := clean_variable_name]
phewas[, value:=-log10(as.numeric(p))]

signif_adj = -log10(0.05/ntest)
phewas <- phewas[order(p),]
k<-phewas[, .SD[order(p)][1:10][p<0.05/ntest], by = "id.exposure", .SDcols = c("id.outcome", "p")]
k[,annotate:=TRUE]
phewas <- merge(phewas, k[,.(id.exposure,id.outcome,annotate)], by = c("id.exposure", "id.outcome"), all.x = TRUE)

library(MetBrewer)
source("/mnt/sda/gagelo01/Projects/small_MR_exploration/F2_finngen/Analysis/my_phenotypePlot_function.R")
list_phewas <- map(unique(phewas$hgnc), function(hgnc_toplot) {
  plot <-phenotypePlot(phewas[hgnc == hgnc_toplot,.(hgnc, groupnum, direction, p, description, phenotype, value, annotate)],
                       suggestive.line=-log10(0.05),significant.line=signif_adj,direction=T,color.palette=rep(met.brewer("Redon"),12),
                       x.group.labels=T,y.axis.interval=5, annotate.size=3, size.x.labels=8, size.y.labels=14,
                       x.axis.label="", y.axis.label="-log(p-value)",
                       title=phewas[hgnc == hgnc_toplot,hgnc_mechanism][1])
  return(plot)
  
})

names(list_phewas)<-unique(phewas$hgnc)
panel_2 <- cowplot::plot_grid(list_phewas$LPL, list_phewas$ANGPTL4,
                              list_phewas$APOC3, list_phewas$PCSK9,
                              labels = c('A)', 'B)', "C)", "D)"),
                              label_size = 14, ncol = 2)


ggsave(filename = paste0("Results/lplpathway_figure5.png"), plot = panel_2,
       width = 1100/72,height = 900/72,units="in",scale=1, device = "png")

# ##### Figure3. Score median.
# dat <- fread("Data/Modified/data_cox.txt")
# prs_lplpathway<-c("PRS_ANGPTL4_trait_16_4", "PRS_APOC3_trait_16_4", "PRS_LPL_trait_16_4")
# dat[, PRS_sum := apply(.SD, 1, function(x) sum(x, na.rm = TRUE)),.SDcols = prs_lplpathway]
# dat[,("PRS_sum"):=lapply(.SD, function(x) (x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)) ,.SDcols = "PRS_sum"]
# 
# colnom <- colnames(dat)[grepl("^PRS_", colnames(dat))]
# colnom_new <- paste0("dicho_", colnom)
# dat[, (colnom_new):=lapply(.SD, function(x) cut(x, breaks=c(-Inf, median(x, na.rm = TRUE), Inf), labels=c('low', "high"))),.SDcols = colnom]
# 
# #create the cox resluts
# colnames(dat)[grepl("^dicho",colnames(dat))]
# trans <- c(dicho_PRS_PCSK9_trait_16_2="PCSK9", "dicho_PRS_HMGCR_trait_16_2"="HMGCR", "dicho_PRS_sum"="LPL")
# prs1= "dicho_PRS_sum"
# prs2= c("dicho_PRS_PCSK9_trait_16_2", "dicho_PRS_HMGCR_trait_16_2")
# DV=c("CAD", "T2D")
# arg<-expand_grid(prs1,prs2,DV)%>% as.data.table(.)
# mylevels<-c("high-high", "high-low", "low-high", "low-low")
# res <- map(split(arg, 1:arg[,.N]), function(x) { 
#   dat[,funky := factor(paste0(get(x$prs1), "-", get(x$prs2)), levels = mylevels)]
#   res <- GagnonMR::run_coxph_wrapper(dat = dat[get(paste0(x$DV, "_toinclude"))==1,], IV = "funky", DV = x$DV, cov_inc = "")
#   res[,gene1:=trans[x$prs1]]
#   res[,gene2:=trans[x$prs2]]
#   res[1, exposure:="high-high"]
#   res[, exposure := gsub("funky", "", exposure)]
#   return(res)
# }) %>% rbindlist(., fill = TRUE)
# 
# dt_forest <- data.table(res)
# setnames(dt_forest, c("HR"), c("b"))
# colnom <- c("b","lci","uci")
# dt_forest[ , (colnom) := lapply(.SD, log) , .SDcols = colnom]
# dt_forest[,panel := outcome]
# dt_forest <- dt_forest[order(outcome, gene1, gene2, exposure)]
# dt_forest[, c("g1", "g2") := tstrsplit(exposure, "-")]
# 
# dt_forest <- map(split(dt_forest,list(dt_forest$gene2,dt_forest$panel)), function(x) {
#   k<-data.table(g2=x$gene2[1], panel=x$panel[1])
#   x<-rbind(k,x, fill = TRUE)
# }) %>% rbindlist(.,fill=TRUE)
# dictio <- c(CAD="Coronary arery disease", T2D="Type 2 diabetes")
# dt_forest[,panel:=dictio[panel]]
# p <- my_forestplotter_fancy(data = data.table(dt_forest), 
#                             col.below.header = NULL, 
#                             col.header = NULL, 
#                             col.header.heading = NULL,
#                             col.left = c("g1", "g2"),
#                             col.left.heading = c("LPL","Gene 2"),
#                             effect.name = "HR (95% CI)", 
#                             col.right = NULL, 
#                             exponentiate = TRUE,
#                             xlab = paste0(""))
# p
# 
# ggsave(filename = paste0("Results/lplpathway_figure3.png"), plot = p,
#        width = 722/72,height = 211/72,units="in",scale=1, device = "png") 

###supplementary figure 1
library(geni.plots)
arg <- dt_gene_region[!grepl("prot|eqtl", id),.(id_small, hgnc, gene_region)]
ID <- c("dis-13-1", "dis-24-4", "trait-19-5")
order_exposure <-  c("LDL cholesterol", "Triglyceride",  "Apolipoprotein b", "Coronary artery disease","Type 2 diabetes")
dt<-data.table(hgnc = c("ANGPTL4", "APOC3", "LPL", "PCSK9"), order_hgnc = 1:4)
arg<- merge(arg, dt, by = "hgnc")
arg <- arg[order(order_hgnc), ]

list_stack <- map(arg$hgnc, function(x) {
  inst <- GagnonMR:::intern_vcfpath_to_TwoSampelMR_region(vcf = paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/", c(ID, arg[hgnc==x,id_small]),"/", c(ID, arg[hgnc==x,id_small]), ".vcf.gz"), chrompos = arg[hgnc==x, gene_region], parameters = parameters)
  inst <- merge(inst, df_index[, .(id, clean_variable_name)], by.x = "id.exposure", by.y = "id")
  inst[,exposure:=clean_variable_name]
  inst_mvmr <- GagnonMR::prepare_for_mvmr(inst, inst, should_clump = FALSE, harmonise_strictness = 1)
  df_aligned<- inst_mvmr
  traits_inorder <- order_exposure[order_exposure %in% unique(inst_mvmr$exposure)] 
  stack <- stack_assoc_plot_wrapper(df_aligned = inst_mvmr, top_snp = NULL, traits_inorder = traits_inorder)
  return(stack)
})

panel_1 <- cowplot::plot_grid(list_stack[[1]], list_stack[[2]],
                              list_stack[[3]], list_stack[[4]],
                              labels = c('A)', 'B)', "C)", "D)"),
                              label_size = 14, ncol = 2)


ggsave(filename = paste0("Results/lplpathway_supfigure1.png"), plot = panel_1,
       width = 700/72,height = 1200/72,units="in",scale=1, device = "png")

message("This script finished without errors")

