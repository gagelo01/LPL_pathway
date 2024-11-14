#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(GagnonMR)
library(forestploter)
library(ggrepel)
setwd("/home/gagelo01/workspace/Projects/2024/LPL_pathway/")
source("Analysis/tosourcedrug.R")
source("/mnt/sda/gagelo01/Learning/gwasvcf/forestplotter.R")
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
  dt_forest[,group:="group1"]
  dt<-data.table(dt_forest)
  
  
  p <- plot_forest(data = dt[id.outcome %in% x$ID[[1]], ], 
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


ggsave(filename = paste0("Results/lplpathway_figure2.tiff"), plot = panel_1, dpi = 100,
       width = 1100/72,height = 550/72,units="in",scale=1, device = "tiff")

#####Fig 3 Genetic mimicry analysis#######
# #####Concordance with LPL enhancement######
k1 <- res_multicis[grepl("trait-33-", id.outcome)&method%in%c("Inverse variance weighted"), ]
scatter <- dcast(k1, id.outcome ~ hgnc, value.var = c("b", "se"))

colnom<-colnames(scatter)[grepl("^b_", colnames(scatter))]
scatter[,(colnom):=lapply(.SD, function(x) x*-1) , .SDcols = colnom]
cor(scatter[,.SD,.SDcols = colnom], method = c("pearson"))
# scatter <- merge(scatter, distinct(dt), by.x = "id.outcome", by.y = "nom")
data.scatter <- data.table(scatter)
data.scatter <- merge(data.scatter, df_index[,.(id, clean_variable_name)], by.x = "id.outcome", by.y = "id", all.x = TRUE)

plot_scatter <- function(scatter,
                         my_xlab = "Genetically proxied ANGPTL4 inhibition",
                         my_ylab = "Genetically proxied LPL enhancing",
                         title) {
  pearson_cor <- cor(scatter$b_x, scatter$b_y, use="complete.obs", method = "pearson")
  
k <- with(scatter, TwoSampleMR::mr_egger_regression(b_exp = b_x, b_out = b_y, se_exp = se_x, se_out = se_y))
lm_res <- k$mod$coefficients %>% as.data.table(.)
setnames(lm_res ,"Estimate", "estimate")
  # lm_res <- lm(formula = "b_x ~ b_y", data = scatter) %>%
  #   broom::tidy() %>%
  #   as.data.table()
  
  
  ggplot2::ggplot(data = scatter, ggplot2::aes(x = b_x,
                                               y = b_y, colour = value)) + ggplot2::geom_errorbarh(ggplot2::aes(xmin = b_x -
                                                                                                                  se_x, xmax = b_x + se_x),
                                                                                                   colour = "grey", height = 0) + ggplot2::geom_errorbar(ggplot2::aes(ymin = b_y -
                                                                                                                                                                        se_y, ymax = b_y + se_y),
                                                                                                                                                         colour = "grey", width = 0) + ggplot2::geom_point(size =2) +
    ggplot2::geom_abline(ggplot2::aes(intercept = lm_res[1,estimate], slope = lm_res[2, estimate]), linetype = "dashed") +
    geom_text_repel(
      aes(label = ifelse(toannotate, clean_variable_name, "")),  # Annotate only where to_annotate is TRUE
      vjust = -1,   # Adjust the vertical position of the text
      color = "black"  # Set text color
    ) +
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
      legend.position=c(0.2, 0.8)) +
    # labs(caption = "[Each measures are scaled to 1-SD lower triglyceride levels]", size = 4) +
    annotate(geom="text", x=0.7, y=-0.8, label= paste0("~italic(r)==", round(pearson_cor, digits = 2)), parse=TRUE, size=7) +
    ggtitle(title)
}

data.scatter[,value:="metabolites"]
data.scatter[,toannotate := FALSE]

scatter <- data.table(data.scatter)
scatter <- scatter[,.(b_ANGPTL4, se_ANGPTL4, b_LPL, se_LPL, value, clean_variable_name)]
setnames(scatter, c("b_ANGPTL4", "se_ANGPTL4", "b_LPL", "se_LPL"), c("b_x", "se_x", "b_y", "se_y"))
scatter[,toannotate := FALSE]
fig3_a<- plot_scatter(scatter = scatter,
                      my_xlab = "Genetically proxied ANGPTL4 inhibition",
                      my_ylab = "Genetically proxied LPL enhancement",
                      title = "")

scatter <- data.table(data.scatter)
scatter <- scatter[,.(b_APOC3, se_APOC3, b_LPL, se_LPL, value)]
setnames(scatter, c("b_APOC3", "se_APOC3", "b_LPL", "se_LPL"), c("b_x", "se_x", "b_y", "se_y"))
scatter[,toannotate := FALSE]
fig3_b <- plot_scatter(scatter = scatter,
                       my_xlab = "Genetically proxied APOC3 inhibition",
                       my_ylab = "Genetically proxied LPL enhancement",
                       title = "")

scatter <- data.table(data.scatter)
scatter <- scatter[,.(id.outcome,b_ANGPTL3, se_ANGPTL3, b_LPL, se_LPL, value, clean_variable_name)]
setnames(scatter, c("b_ANGPTL3", "se_ANGPTL3", "b_LPL", "se_LPL"), c("b_x", "se_x", "b_y", "se_y"))
scatter[,toannotate := FALSE]
# scatter[,toannotate := id.outcome %in% c("trait-33-5","trait-33-29", "trait-33-31")]

# k <-data.scatter[, .(id.outcome, b_LPL, b_ANGPTL3)]
# k[,diffabs := abs(b_ANGPTL3 - b_LPL)]
# k<-merge(k, df_index[,.(id, clean_variable_name)], by.x = "id.outcome", by.y = "id")
# k[order(-diffabs), ][1:20,]

fig3_c <- plot_scatter(scatter = scatter,
                       my_xlab = "Genetically proxied ANGPTL3 inhibition",
                       my_ylab = "Genetically proxied LPL enhancement",
                       title = "")

scatter <- data.table(data.scatter)
scatter <- scatter[,.(b_PCSK9, se_PCSK9, b_LPL, se_LPL, value)]
setnames(scatter, c("b_PCSK9", "se_PCSK9", "b_LPL", "se_LPL"), c("b_x", "se_x", "b_y", "se_y"))
scatter[,toannotate := FALSE]
fig3_d <- plot_scatter(scatter = scatter,
                       my_xlab = "Genetically proxied PCSK9 inhibition",
                       my_ylab = "Genetically proxied LPL enhancement",
                       title = "") +theme(legend.position="none")

scatter <- data.table(data.scatter)
scatter <- scatter[,.(b_HMGCR, se_HMGCR, b_LPL, se_LPL, value)]
setnames(scatter, c("b_HMGCR", "se_HMGCR", "b_LPL", "se_LPL"), c("b_x", "se_x", "b_y", "se_y"))
scatter[,toannotate := FALSE]
fig3_e <- plot_scatter(scatter = scatter,
                       my_xlab = "Genetically proxied HMGCR inhibition",
                       my_ylab = "Genetically proxied LPL enhancement",
                       title = "")+ theme(legend.position = "none")

scatter <- data.table(data.scatter)
scatter <- scatter[,.(b_HMGCR, se_HMGCR, b_PCSK9, se_PCSK9, value)]
setnames(scatter, c("b_HMGCR", "se_HMGCR", "b_PCSK9", "se_PCSK9"), c("b_x", "se_x", "b_y", "se_y"))
scatter[,toannotate := FALSE]
fig3_f <- plot_scatter(scatter = scatter,
                       my_xlab = "Genetically proxied HMGCR inhibition",
                       my_ylab = "Genetically proxied PCSK9 enhancement",
                       title = "")

panel_3 <- cowplot::plot_grid(fig3_a, fig3_b, fig3_c, fig3_d, fig3_e, fig3_f, labels = c('A)', 'B)', "C)", "D)", "E)", "F)"),
                              label_size = 22, ncol = 3, nrow = 2)
ggsave(plot = panel_3, filename = paste0("Results/lplpathway_figure3.tiff"),
       dpi = 100, width = 1310/72,height = 900/72,units="in",scale=1, device = "tiff")

#####Figure 4 Interaction MR analysis##############
res_logit <- fread("Data/Modified/res_logit.txt")
k <- res_logit[cov_inc == "+ age_enrollment + sex + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10" &  
                 grepl("^PRS", exposure) & str_count(string = IV,pattern = "\\*")==1, ]
k[!grepl(":", exposure), b := b*-1]

dt_forest <- k[grepl("trait_16_39", IV) & grepl("trait_16_37", IV)]
dt_forest[,c("V1", "V2"):=tstrsplit(IV, "\\*")]
trans <- c(PRS_ANGPTL3_trait_16_39 = "ANGPTL3",PRS_ANGPTL4_trait_16_39 = "ANGPTL4", PRS_APOC3_trait_16_39 = "APOC3", PRS_LPL_trait_16_39 = "LPL",
           PRS_HMGCR_trait_16_37 = "HMGCR", PRS_PCSK9_trait_16_37 = "PCSK9")
dt_forest[,V1:=trans[V1]]
dt_forest[,V2:=trans[V2]]
tobind<- data.table(dt_forest)
tobind[,todump:=V1][,V1:=V2][,V2:=todump][,todump:=NULL]
dt_forest <- rbind(dt_forest, tobind)
dt_forest <- dt_forest[V1%in%c("ANGPTL3","ANGPTL4","APOC3","LPL")&V2%in%c("HMGCR", "PCSK9")]
dictio <- c(CAD="Coronary arery disease", T2D="Type 2 diabetes")
dt_forest[,panel := outcome]
dt_forest[,panel:=dictio[panel]]
dt_forest[,lci:=b-se*1.96]
dt_forest[,uci:=b+se*1.96]
dt_forest[,group := exposure %>% fifelse(grepl(":", .), "Interaction", .) %>%
            fifelse(grepl("trait_16_39", .), "LPL pathway", .) %>%
            fifelse(grepl("trait_16_37", .), "LDLR pathway", .)]
dt_forest[,group := factor(group, levels = c("LPL pathway", "LDLR pathway", "Interaction"))]
dt_forest[,mechanism1 := fifelse(V1=="LPL", "activation", "inhibition")]
dt_forest[,mechanism2 := fifelse(V2=="LPL", "activation", "inhibition")]
dt_forest[,V1:=paste0(V1,"\n",mechanism1)]
dt_forest[,V2:=paste0(V2,"\n",mechanism2)]
dt_forest<-dt_forest[order(V1, V2, group)]


tm<-forest_theme(base_size = 10,
                 ci_Theight = 0.2, ci_col = c("#5884E5", "#9E131E", "#4CAF50"),#"#377eb8", "#4daf4a"
                 legend_position = "bottom", legend_name = "Exposures",
                 legend_value = dt_forest$group %>% unique)

plot4 <- plot_forest(data = data.table(dt_forest),
                     col.below.header = NULL,
                     col.header = NULL,
                     col.header.heading = NULL,
                     col.left = c("V1", "V2"),
                     col.left.heading = c("LPL pathway", "LDLR pathway"),
                     col.left.toformat = "V1",
                     effect.name = "OR (95% CI)",
                     col.right = NULL,
                     exponentiate = TRUE,
                     xlab = paste0(""),
                     tm = tm)

ggsave(filename = paste0("Results/lplpathway_figure4.tiff"), plot = plot4, dpi = 100,
       width = 740/72,height = 430/72,units="in",scale=1, device = "tiff")

# #####Figure 5 Association og LPL pathway targets with clinical toxicology panel########
# dt_ctp <- res_multicis[idscale == "trait-16-39" & id.outcome %in% ctp$id, ]
# xscale <- df_index[id=="trait-16-39", clean_variable_name]
# dt_ctp <- dt_ctp[method == "Inverse variance weighted", ]
# dt_ctp[,b := b*-1]
# dt_ctp[, mechanism := ifelse(hgnc%in%"LPL", "activation", "inhibition")]
# dt_ctp[, hgnc_mechanism := paste0(hgnc, " ", mechanism)]
# dt_ctp[,lci := b-1.96*se]
# dt_ctp[,uci := b+1.96*se]
# dt_ctp[,panel := hgnc_mechanism]
# dt_ctp[,group := "group1"]
# dt_ctp <- merge(dt_ctp, ctp[,.(id, category, subcategory)], by.x = "id.outcome", by.y = "id", all.x = TRUE)
# dt_ctp <- dt_ctp[order(category, subcategory)]
# 
# dt <- data.table(dt_ctp)
# p <- plot_forest(data = dt, 
#                             col.below.header = "clean_variable_name", 
#                             col.header = "category", 
#                             col.header.heading = "Outcome", 
#                             col.left = NULL,
#                             col.left.heading = NULL,
#                             col.left.toformat = NULL,
#                             effect.name = "Beta (95% CI)", 
#                             col.right = NULL, 
#                             exponentiate = FALSE,
#                             annotate_pvalue=TRUE,
#                             ntest_adjustment=dt[,.N],
#                             xlab = "")
# 
# 
# ggsave(filename = paste0("Results/lplpathway_figure5.png"), plot = p,
#        width = 1200/72,height = 680/72,units="in",scale=1, device = "png")
# 


####Figure6. Association of LPL pathway with diseases in Finngen. (phewas plot)
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
                       x.axis.label="", y.axis.label="-log(p-value)", grid_rows = NULL, grid_cols = NULL, grid_scales = "free_y",
                       title=phewas[hgnc == hgnc_toplot,hgnc_mechanism][1]) + 
    scale_y_continuous(trans = "sqrt")
  return(plot)
  
})

names(list_phewas)<-unique(phewas$hgnc)
panel_2 <- cowplot::plot_grid(list_phewas$LPL, list_phewas$ANGPTL4,
                              list_phewas$APOC3, list_phewas$ANGPTL3,
                              labels = c('A)', 'B)', "C)", "D)"),
                              label_size = 14, ncol = 2)


ggsave(filename = paste0("Results/lplpathway_figure5.tiff"), plot = panel_2, dpi = 100,
       width = 1100/72,height = 900/72,units="in",scale=1, device = "tiff")

# # ##### Figure6. Score median.
# dat <- fread("Data/Modified/data_cox.txt")
# prs_lplpathway<-c("PRS_ANGPTL4_trait_16_39", "PRS_APOC3_trait_16_39", "PRS_LPL_trait_16_39")
# dat[, PRS_sum := apply(.SD, 1, function(x) sum(x, na.rm = TRUE)),.SDcols = prs_lplpathway]
# dat[,("PRS_sum"):=lapply(.SD, function(x) (x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)) ,.SDcols = "PRS_sum"]
# 
# colnom <- colnames(dat)[grepl("^PRS_", colnames(dat))]
# colnom_new <- paste0("dicho_", colnom)
# dat[, (colnom_new):=lapply(.SD, function(x) cut(x, breaks=c(-Inf, median(x, na.rm = TRUE), Inf), labels=c('low', "high"))),.SDcols = colnom]
# 
# #create the logit resluts
# k <- colnames(dat)[grepl("^dicho",colnames(dat))]
# trans <- k %>% gsub("dicho_PRS_", "", .) %>% gsub("_.*$", "", .)
# names(trans) <- k
# prs2= c("dicho_PRS_PCSK9_trait_16_37", "dicho_PRS_HMGCR_trait_16_37")
# prs1= trad[!(trad%in%c(prs2,  "dicho_PRS_sum" ))]
# DV=c("CAD", "T2D")
# arg<-expand_grid(prs1,prs2,DV)%>% as.data.table(.)
# mylevels<-c("high-high", "high-low", "low-high", "low-low")
# res <- map(split(arg, 1:arg[,.N]), function(x) {
#   dat[,funky := factor(paste0(get(x$prs1), "-", get(x$prs2)), levels = mylevels)]
#   res<-glm(paste0(x$DV, "_censored ~ ",  "funky", x$cov_inc),data =  dat[get(paste0(x$DV, "_toinclude"))==1,], family = "binomial")
#   res <- broom::tidy(res) %>% as.data.table(.)
#   colnames(res)<-c("exposure", "b", "se", "todump", "pval")
#   res[,todump:=NULL]
#   res[1,b:=NA][1,se:=NA][1,pval:=NA]
#   # res <- GagnonMR::run_coxph_wrapper(dat = dat[get(paste0(x$DV, "_toinclude"))==1,], IV = "funky", DV = x$DV, cov_inc = "")
#   res[,gene1:=trans[x$prs1]]
#   res[,gene2:=trans[x$prs2]]
#   res[1, exposure:="high-high"]
#   res[, exposure := gsub("funky", "", exposure)]
#   res[,outcome:=x$DV]
#   return(res)
# }) %>% rbindlist(., fill = TRUE)
# 
# dt_forest <- data.table(res)
# dt_forest[,lci := b-1.96*se][,uci:=b+1.96*se]
# dt_forest[,panel := outcome]
# 
# dt_forest[, exposure_plot := exposure %>% fifelse(. == "high-high", "both score above median", .) %>% fifelse(.=="high-low", paste0(gene1, " score below median"),.)%>%
#                      fifelse(.=="low-high", paste0(gene2, " score below median"), .) %>% fifelse(.=="low-low", "both score below median", .)]
# 
# dt_forest[, group := exposure %>% fifelse(. == "high-high", "both score above median", .) %>% fifelse(.=="high-low", paste0("LPL score below median"),.)%>%
#             fifelse(.=="low-high", paste0("LDLR score below median"), .) %>% fifelse(.=="low-low", "both score below median", .)]
# # dt_forest <- dt_forest[order(outcome, gene1, gene2, exposure)]
# # dt_forest[, c("g1", "g2") := tstrsplit(exposure, "-")]
# # 
# # dt_forest <- map(split(dt_forest,list(dt_forest$gene2,dt_forest$panel)), function(x) {
# #   k<-data.table(g2=x$gene2[1], panel=x$panel[1])
# #   x<-rbind(k,x, fill = TRUE)
# # }) %>% rbindlist(.,fill=TRUE)
# dictio <- c(CAD="Coronary arery disease", T2D="Type 2 diabetes")
# dt_forest[,panel:=dictio[panel]]
# 
# ci_col <- c("#5884E5", "#9E131E", "#F4A261", "#2A9D8F", "#E9C46A","#7B2CBF","#707070")
# tm<-forest_theme(base_size = 10,
#                  ci_Theight = 0.1, ci_col = ci_col[(1:length(unique(dt_forest$group)))], #"#5884E5", "#9E131E", "#F4A261", "#2A9D8F"
#                  legend_position = "bottom", legend_name = "Adjusted for",
#                  legend_value = dt_forest$group %>% unique)
# 
# p <- plot_forest(data = data.table(dt_forest),
#                             col.below.header = NULL,
#                             col.header = NULL,
#                             col.header.heading = NULL,
#                             col.left = c("gene1", "gene2"),
#                             col.left.heading = c("Gene 1","Gene 2"),
#                             col.left.toformat = c("gene1", "gene2"),
#                             effect.name = "OR (95% CI)",
#                             col.right = NULL,
#                             exponentiate = TRUE,
#                             xlab = paste0(""),
#                  tm = tm )
# p
# 
# ggsave(filename = paste0("Results/lplpathway_figure3.png"), plot = p,
#        width = 722/72,height = 211/72,units="in",scale=1, device = "png")
# 
# ###supplementary figure 1
# library(geni.plots)
# arg <- dt_gene_region[!grepl("prot|eqtl", id),.(id_small, hgnc, gene_region)]
# ID <- c("dis-13-1", "dis-24-4", "trait-19-5")
# order_exposure <-  c("LDL cholesterol", "Triglyceride",  "Apolipoprotein b", "Coronary artery disease","Type 2 diabetes")
# dt<-data.table(hgnc = c("ANGPTL4", "APOC3", "LPL", "PCSK9"), order_hgnc = 1:4)
# arg<- merge(arg, dt, by = "hgnc")
# arg <- arg[order(order_hgnc), ]
# 
# list_stack <- map(arg$hgnc, function(x) {
#   inst <- GagnonMR:::intern_vcfpath_to_TwoSampelMR_region(vcf = paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/", c(ID, arg[hgnc==x,id_small]),"/", c(ID, arg[hgnc==x,id_small]), ".vcf.gz"), chrompos = arg[hgnc==x, gene_region], parameters = parameters)
#   inst <- merge(inst, df_index[, .(id, clean_variable_name)], by.x = "id.exposure", by.y = "id")
#   inst[,exposure:=clean_variable_name]
#   inst_mvmr <- GagnonMR::prepare_for_mvmr(inst, inst, should_clump = FALSE, harmonise_strictness = 1)
#   df_aligned<- inst_mvmr
#   traits_inorder <- order_exposure[order_exposure %in% unique(inst_mvmr$exposure)] 
#   stack <- stack_assoc_plot_wrapper(df_aligned = inst_mvmr, top_snp = NULL, traits_inorder = traits_inorder)
#   return(stack)
# })
# 
# panel_1 <- cowplot::plot_grid(list_stack[[1]], list_stack[[2]],
#                               list_stack[[3]], list_stack[[4]],
#                               labels = c('A)', 'B)', "C)", "D)"),
#                               label_size = 14, ncol = 2)
# 
# 
# ggsave(filename = paste0("Results/lplpathway_supfigure1.png"), plot = panel_1,
#        width = 700/72,height = 1200/72,units="in",scale=1, device = "png")

# #####Figure 6a######
# res_cox <- fread("Data/Modified/res_cox.txt")
# dt_forest <- res_cox[grepl("^PRS", exposure) & !grepl("\\*", IV)& cov_inc ==  "+ age_enrollment + sex + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10"]
# trans <- c(PRS_ANGPTL4_trait_16_39 = "ANGPTL4", PRS_APOC3_trait_16_39 = "APOC3", PRS_LPL_trait_16_39 = "LPL",
#            PRS_HMGCR_trait_16_37 = "HMGCR", PRS_PCSK9_trait_16_37 = "PCSK9")
# 
# setnames(dt_forest, c("HR"), c("b"))
# colnom <- c("b","lci","uci")
# dt_forest[ , (colnom) := lapply(.SD, function(x) log(1/x)) , .SDcols = colnom]
# setnames(dt_forest, c("lci","uci"), c("uci", "lci"))
# dt_forest[,panel := outcome]
# 
# # dt_forest <- map(split(dt_forest,list(dt_forest$gene2,dt_forest$panel)), function(x) {
# #   k<-data.table(g2=x$gene2[1], panel=x$panel[1])
# #   x<-rbind(k,x, fill = TRUE)
# # }) %>% rbindlist(.,fill=TRUE)
# dictio <- c(CAD="Coronary arery disease", T2D="Type 2 diabetes")
# dt_forest[,panel:=dictio[panel]]
# dt_forest[,exposure:=trans[exposure]]
# dt_forest[,pathway := fifelse(exposure %in% c("PCSK9", "HMGCR"), "LDLR pathway", "LPL pathway")]
# dt_forest[,mechanism := fifelse(exposure %in% "LPL", "activation", "inhibition")][,exposure_mechanism := paste0(exposure, " ", mechanism)]
# 
# plot6a <- my_forestplotter_fancy(data = data.table(dt_forest),
#                                  col.below.header = "exposure_mechanism",
#                                  col.header = "pathway",
#                                  col.header.heading = "GRS",
#                                  col.left = NULL,
#                                  col.left.heading = NULL,
#                                  effect.name = "HR (95% CI)",
#                                  col.right = NULL,
#                                  exponentiate = TRUE,
#                                  xlab = paste0(""))
# 
# ####Figure 6B#######
# k<-res_cox[grepl("^PRS", exposure) & outcome %in% c("CAD", "T2D") & cov_inc ==  "+ age_enrollment + sex + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10"]
# k <- k[str_count(exposure, ":")==1 & str_count(IV, "\\*")==1, ]
# k[,c("exposure1", "exposure2") := tstrsplit(exposure, ":")]
# k[,exposure1 := trans[exposure1]]
# k[,exposure2 := trans[exposure2]]
# ktobind<- data.table(k)
# ktobind[,exposure:=exposure1][,exposure1:=exposure2][,exposure2:=exposure]
# k<-rbind(k,ktobind)
# dt<-data.table(exposure1= c("ANGPTL4", "APOC3", "LPL", "HMGCR", "PCSK9"))
# dt[,exposure2 := exposure1]
# dt[,outcome:="CAD"]
# k <- rbind(k, dt, fill = TRUE)
# dt[,outcome:="T2D"]
# k <- rbind(k, dt, fill = TRUE)
# k <- k[,.(exposure1,exposure2, outcome, HR, pval)]
# bonferroni_threshold <- 0.05/20
# 
# k <- k[order(exposure1,exposure2)]
# k[, log10_pval := -log10(pval)]
# dat_plot_croix <- data.table(k)
# dat_plot_croix[,shape_point := fifelse(pval >= bonferroni_threshold, "Non-significant", NA)]
# dat_plot_croix[pval < bonferroni_threshold,c("HR", "log10_pval") := NA]
# 
# dat_plot_rond <- data.table(k)
# dat_plot_rond[,shape_point := fifelse(pval < bonferroni_threshold, "rond", NA)]
# dat_plot_rond[pval >= bonferroni_threshold,c("HR", "log10_pval") := NA]
# 
# balloon_plot = ggplot2::ggplot() +
#   ggplot2::geom_point(data = dat_plot_croix, ggplot2::aes(x = exposure1, y = exposure2, shape = factor(shape_point)), size = 2, color = "gray20") +
#   ggplot2::scale_shape_manual(name = "", values = c(4)) 
# # k[,shape_point := fifelse(pval < bonferroni_threshold, "rond", "Non-significant")]
# 
# plot6b <- balloon_plot +
#   ggplot2::geom_point(data = dat_plot_rond, ggplot2::aes(x = exposure1, y = exposure2, size = log10_pval, color = HR)) +
#   # ggplot2::geom_point(data = dat_plot_rond, ggplot2::aes(x = exposure1, y = exposure2, size = log10_pval, color = HR)) +
#   # ggplot2::ggplot() +
#   #   ggplot2::geom_point(data = k, ggplot2::aes(x = exposure1, y = exposure2, size = log10_pval, color = HR, shape = factor(shape_point))) +
#   ggplot2::facet_grid(cols = vars(outcome)) +
#   # ggplot2::scale_color_gradient2(name = "Hazard ratio",
#   #                                low = scales::muted("#5884E5"),
#   #                                mid = "white",
#   #                                high = scales::muted("#9E131E"),
#   #                                midpoint = 1
#   # ) +
#   ggplot2::scale_shape_manual(name = "", values = c(4,1)) +
#   ggplot2::scale_size(guide = "none") +
#   ggplot2::scale_color_continuous(guide = "none") +
#   # ggplot2::scale_size(name = expression(-Log[10](P)), range = c(4,7)) +
#   ggplot2::coord_fixed(clip = "off", ratio = 1) +
#   # ggplot2::guides(size = guide_legend(order = 1),
#   #                 shape = guide_legend(order = 2)) +
#   ggplot2::theme(
#     panel.grid.major.y = element_line(size = 0.25, colour = "gray60"),
#     panel.grid.major.x = element_line(size = 0.25, colour = "gray60"),
#     panel.grid.minor.y = element_blank(),
#     panel.grid.minor.x = element_blank(),
#     panel.background = element_blank(),
#     # plot.margin = margin(t = 2, r = 0.5, b = 0.5, l = 0.5, "cm"),
#     legend.position = "top",
#     legend.text = element_text(
#       color = "gray20",
#       size = 10),
#     legend.title=element_blank(),
#     legend.margin = margin(0, 0, 0, 0),
#     legend.spacing.x = unit(0, "mm"),
#     legend.spacing.y = unit(0, "mm"),
#     plot.margin=grid::unit(c(0,0,0,0), "mm"),
#     # legend.key = element_rect(fill = "transparent", colour = "transparent"),
#     # legend.key.size = unit(0.8, "cm"),
#     axis.title = element_blank(),
#     axis.line = element_line(size = 0.5, colour = "gray20"),
#     axis.ticks = element_line(size = 0.5, colour = "gray20"),
#     axis.text.y = element_text(
#       size = 10,
#       colour = "gray20"
#     ),
#     axis.text.x = element_text(
#       angle = 60,
#       size = 8,
#       hjust = 1,
#       face = "plain",
#       colour = "gray20"
#     ))
# 
# plot6b 
# #####Figure 6c#######
# dat <- fread("Data/Modified/data_cox.txt")
# colnom <- colnames(dat)[grepl("^PRS_", colnames(dat))]
# colnom_new <- paste0("dicho_", colnom)
# dat[, (colnom_new):=lapply(.SD, function(x) cut(x, breaks=c(-Inf, median(x, na.rm = TRUE), Inf), labels=c('low', "high"))),.SDcols = colnom]
# 
# #create the cox resluts
# prs1= "dicho_PRS_HMGCR_trait_16_37" 
# prs2= c("dicho_PRS_PCSK9_trait_16_37")
# DV=c("T2D")
# arg<-expand_grid(prs1,prs2,DV)%>% as.data.table(.)
# mylevels<-c("high-high", "high-low", "low-high", "low-low")
# names(trans)<-paste0("dicho_",names(trans))
# res <- map(split(arg, 1:arg[,.N]), function(x) {
#   dat[,funky := factor(paste0(get(x$prs1), "-", get(x$prs2)), levels = mylevels)]
#   res <- GagnonMR::run_coxph_wrapper(dat = dat[get(paste0(x$DV, "_toinclude"))==1,], IV = "funky", DV = x$DV, cov_inc = "+ age_enrollment + sex + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10")
#   res<-res[grepl("^funky",exposure)]
#   res[,gene1:=trans[x$prs1]]
#   res[,gene2:=trans[x$prs2]]
#   res[1, exposure:=mylevels[1]]
#   res[, exposure := gsub("funky", "", exposure)]
#   return(res)
# }) %>% rbindlist(., fill = TRUE)
# 
# 
# 
# dt_forest <- data.table(res)
# setnames(dt_forest, c("HR"), c("b"))
# colnom <- c("b","lci","uci")
# dt_forest[ , (colnom) := lapply(.SD, log) , .SDcols = colnom]
# dt_forest[,panel := outcome]
# dt_forest <- dt_forest[order(outcome, gene1, gene2, exposure)]
# dt_forest[, c("g1", "g2") := tstrsplit(exposure, "-")]
# 
# # dt_forest <- map(split(dt_forest,list(dt_forest$gene2,dt_forest$panel)), function(x) {
# #   k<-data.table(g2=x$gene2[1], panel=x$panel[1])
# #   x<-rbind(k,x, fill = TRUE)
# # }) %>% rbindlist(.,fill=TRUE)
# dictio <- c(CAD="Coronary arery disease", T2D="Type 2 diabetes")
# dt_forest[,panel:=dictio[panel]]
# 
# dt_forest[,exposure:=trans[exposure]]
# plot6c <- my_forestplotter_fancy(data = data.table(dt_forest),
#                                  col.below.header = NULL,
#                                  col.header = NULL,
#                                  col.header.heading = NULL,
#                                  col.left = c("g1", "g2"),
#                                  col.left.heading = dt_forest[1, c(gene1, gene2)],
#                                  effect.name = "HR (95% CI)",
#                                  col.right = NULL,
#                                  exponentiate = TRUE,
#                                  xlab = paste0(""))
# plot6c
# 
# 
# 
# # left_column <- cowplot::plot_grid(plot5a, plot5c, labels = c('A)', 'C)'), label_size = 14, ncol =1, nrow=2)
# # 
# # cowplot::plot_grid(left_column , plot5b, labels = c('', 'B)'), label_size = 14, ncol = 2, nrow=1)
# 
# bottom_row <- cowplot::plot_grid(plot6b, plot6c, labels = c('B)', 'C)'), label_size = 14, ncol =2, nrow=1)
# 
# panel_interaction <- cowplot::plot_grid(plot6a , bottom_row, labels = c('A)', ''), label_size = 14, ncol = 1, nrow=2,
#                                         rel_heights = c(1.5,2))
# 
# ggsave(filename = paste0("Results/lplpathway_figure6.png"), plot = panel_interaction,
#        width = 936/72,height = 531/72,units="in",scale=1, device = "png")

# res_cox <- fread("Data/Modified/res_cox.txt")
# k <- res_cox[cov_inc == "+ age_enrollment + sex + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10" &  
#                str_count(string = exposure,pattern = ":")==1 & str_count(string = IV,pattern = "\\*")==1, ]
# 
# 
# dt_forest <- k[grepl("trait_16_39", exposure) & grepl("trait_16_37", exposure)]
# dt_forest[,c("V1", "V2"):=tstrsplit(exposure, ":")]
# trans <- c(PRS_ANGPTL4_trait_16_39 = "ANGPTL4", PRS_APOC3_trait_16_39 = "APOC3", PRS_LPL_trait_16_39 = "LPL",
#            PRS_HMGCR_trait_16_37 = "HMGCR", PRS_PCSK9_trait_16_37 = "PCSK9")
# dt_forest[,V1:=trans[V1]]
# dt_forest[,V2:=trans[V2]]
# tobind<- data.table(dt_forest)
# tobind[,exposure:=V1][,V1:=V2][,V2:=exposure]
# dt_forest <- rbind(dt_forest, tobind)
# dt_forest <- dt_forest[V1%in%c("ANGPTL4","APOC3","LPL")&V2%in%c("HMGCR", "PCSK9")]
# 
# setnames(dt_forest, c("HR"), c("b"))
# colnom <- c("b","lci","uci")
# dt_forest[ , (colnom) := lapply(.SD, function(x) log(x)) , .SDcols = colnom]
# dt_forest[,panel := outcome]
# 
# # dt_forest <- map(split(dt_forest,list(dt_forest$gene2,dt_forest$panel)), function(x) {
# #   k<-data.table(g2=x$gene2[1], panel=x$panel[1])
# #   x<-rbind(k,x, fill = TRUE)
# # }) %>% rbindlist(.,fill=TRUE)
# dictio <- c(CAD="Coronary arery disease", T2D="Type 2 diabetes")
# dt_forest[,panel:=dictio[panel]]
# # dt_forest[,pathway := fifelse(exposure %in% c("PCSK9", "HMGCR"), "LDLR pathway", "LPL pathway")]
# # dt_forest[,mechanism := fifelse(exposure %in% "LPL", "activation", "inhibition")][,exposure_mechanism := paste0(exposure, " ", mechanism)]
# dt_forest<-dt_forest[order(V1, V2)]
# # dt_forest[,pval:=p.adjust(p = pval, method = "fdr")]
# 
# plot6 <- my_forestplotter_fancy(data = data.table(dt_forest),
#                                  col.below.header = NULL,
#                                  col.header = NULL,
#                                  col.header.heading = NULL,
#                                  col.left = c("V1", "V2"),
#                                  col.left.heading = c("LPL pathway", "LDLR pathway"),
#                                  col.left.toformat = "V1",
#                                  effect.name = "HR (95% CI)",
#                                  col.right = NULL,
#                                  exponentiate = TRUE,
#                                  xlab = paste0(""))
# 
# ggsave(filename = paste0("Results/lplpathway_figure6.png"), plot = plot6,
#        width = 800/72,height = 150/72,units="in",scale=1, device = "png")

#######fgure 7##############
#per APOB particule the LPL vs. LDLR pathway effect on LDL cholesterol and remnant cholesterol
forstandardise <- res_multicis[method=="Inverse variance weighted" & id.outcome %in% c("trait-19-5"),.(id.exposure, b)]
setnames(forstandardise, "b", "b_forstandardise")
res_multicis_standardise <- res_multicis[id.outcome %in% c("trait-33-133", "trait-33-51"), ][forstandardise, on = "id.exposure", nomatch = 0]
res_multicis_standardise[,c("b", "se"):= map(.SD, function(x) x/b_forstandardise),.SDcols=c("b", "se")]
res_multicis_standardise <- res_multicis_standardise[method=="Inverse variance weighted", ]
wide <- dcast(res_multicis_standardise, id.exposure ~ outcome, value.var = c("b"))
wide[,ratio:=Remnant_cholesterol_nonHDL_nonLDL_cholesterol/Total_cholesterol_in_LDL]


dt_forest <- data.table(res_multicis_standardise)
dt_forest <- merge(dt_forest, wide[,.(id.exposure, ratio)], by.x = "id.exposure", by.y= "id.exposure")
dt_forest[, ratio:=round(ratio, digits=2)]
dt_forest[,ratio:=as.character(ratio)]
dt_forest[seq(2, .N, by = 2), ratio := " "]
dt_forest[,b := b*-1]
dt_forest[, mechanism := ifelse(hgnc%in%"LPL", "activation", "inhibition")]
dt_forest[, hgnc_mechanism := paste0(hgnc, " ", mechanism)]
dt_forest[,lci := b-1.96*se]
dt_forest[,uci := b+1.96*se]
dt_forest[,panel :=stringr::str_pad(" ", 45, 'right', ' ')]
dt_forest<-dt_forest[order(hgnc_mechanism, clean_variable_name)]
dt_forest[,group:=fifelse(outcome=="Total_cholesterol_in_LDL", "LDL cholesterol", "Remnant cholesterol")]
dt_forest[,category:=fifelse(grepl("trait-16-39", id.exposure), "LPL pathway", "LDLR pathway")]
dt_forest<- dt_forest[order(category, hgnc_mechanism), ]
source("/mnt/sda/gagelo01/Learning/gwasvcf/forestplotter.R")


dt <- data.table(dt_forest)
library(forestploter)

tm<-forest_theme(base_size = 10,
                 ci_Theight = 0.4, ci_col = c("#377eb8", "#4daf4a"), #"#5884E5", "#9E131E"
                 legend_position = "bottom", legend_name = "Outcomes",
                 legend_value = dt_forest$group %>% unique)

p <- plot_forest(data = dt, 
                 col.below.header = c("hgnc_mechanism"), 
                 col.header = "category", 
                 col.header.heading = c("Exposure"), 
                 col.left = NULL,
                 col.left.heading = NULL,
                 col.left.toformat = NULL,
                 effect.name = "Effect (95% CI)", 
                 col.right = "ratio",
                 col.right.heading = "Ratio", 
                 exponentiate = FALSE,
                 xlab = paste0("Effects scaled on 1-SD lower circulating apoB levels"),
                 tm = tm)

ggsave(filename = paste0("Results/lplpathway_figure7.tiff"), plot = p, dpi = 100,
       width = 500/72,height = 260/72,units="in",scale=1, device = "tiff")

message("This script finished without errors")

