# ============================================================
# 构建 跨时期 TF→靶基因→GO通路 调控网络（完整代码，修正列匹配）
# ============================================================
library(igraph)


# 1. 手动定义边（ICM + 8-cell） ------------------------------------
edges <- data.frame(
  Source = c(
    # NR4A2 (ICM)
    "NR4A2","NR4A2","NR4A2","NR4A2","NR4A2","NR4A2","NR4A2","NR4A2","NR4A2","NR4A2",
    "NR4A2","NR4A2","NR4A2","NR4A2","NR4A2",
    # ZBTB18 (ICM)
    "ZBTB18","ZBTB18","ZBTB18","ZBTB18","ZBTB18","ZBTB18","ZBTB18","ZBTB18","ZBTB18",
    # PRDM1 (ICM)
    "PRDM1","PRDM1","PRDM1","PRDM1","PRDM1","PRDM1","PRDM1",
    # ZNF135 (8-cell)
    "ZNF135","ZNF135","ZNF135","ZNF135","ZNF135","ZNF135","ZNF135","ZNF135","ZNF135","ZNF135","ZNF135",
    # IGV 共定位 (ICM)
    "NR4A2","PRDM1"
  ),
  Target = c(
    # NR4A2 靶基因
    "Sox1","Neurog2","Foxg1","Lmx1b","Hoxa11","Hoxa13","Runx2","Tfap2a","Lep","Evx1os",
    "Kcnc3","Dgkh","S100b","Prdx4","Rhox4g",
    # ZBTB18 靶基因
    "Sox1","Neurog2","Foxg1","Lmx1b","Hoxa11","Hoxa13","S100b","Evx1os","Runx2",
    # PRDM1 靶基因
    "Igf2","Wnt7b","Tal1","Dlx4","Nrk","Uncx","Slc24a4",
    # ZNF135 靶基因
    "Atoh1","Cer1","Egr3","Fst","Lox","Pdgfa","Sycp1","Tfap2a","Kcnc4","Sox21","Thbs2",
    # 共定位
    "Otx2","Neurod1"
  ),
  Evidence = c(
    rep("HOMER motif scan", 15), rep("HOMER motif scan", 9), rep("HOMER motif scan", 7),
    rep("HOMER motif scan", 11),
    "IGV co-localization (Otx2)", "IGV co-localization (Neurod1)"
  ),
  Stage = c(rep("ICM",15),rep("ICM",9),rep("ICM",7), rep("8-cell",11), "ICM","ICM"),
  stringsAsFactors = FALSE
)

# 2. 节点功能注释 --------------------------------------------------
node_func <- c(
  NR4A2="", ZBTB18="", PRDM1="", ZNF135="",
  Sox1="neural progenitor / forebrain development",
  Neurog2="neurogenesis / neuron differentiation",
  Foxg1="forebrain development",
  Lmx1b="limb/CNS development",
  Hoxa11="embryonic pattern formation",
  Hoxa13="embryonic pattern formation",
  Runx2="osteoblast differentiation / development",
  Tfap2a="neural crest development",
  Lep="energy homeostasis / CNS development",
  Evx1os="non-coding RNA (potential regulatory)",
  Kcnc3="neuronal excitability",
  Dgkh="lipid signaling / development",
  S100b="CNS development (astrocyte)",
  Prdx4="antioxidant",
  Rhox4g="gamete / reproductive development",
  Igf2="growth factor signaling",
  Wnt7b="organ morphogenesis",
  Tal1="hematopoiesis / development",
  Dlx4="forebrain development",
  Nrk="neuronal differentiation",
  Uncx="homeobox transcription factor",
  Slc24a4="calcium homeostasis",
  Atoh1="neuron differentiation",
  Cer1="mesoderm/pattern specification",
  Egr3="early growth response / neurogenesis",
  Fst="BMP inhibitor / pattern formation",
  Lox="extracellular matrix organization",
  Pdgfa="growth factor / development",
  Sycp1="meiotic chromosome organization",
  Kcnc4="neuronal excitability",
  Sox21="neural progenitor / development",
  Thbs2="extracellular matrix / angiogenesis",
  Otx2="forebrain development (IGV validated)",
  Neurod1="neuronal differentiation (IGV validated)"
)

all_nodes <- unique(c(edges$Source, edges$Target))
nodes <- data.frame(
  Node = all_nodes,
  Type = ifelse(all_nodes %in% c("NR4A2","ZBTB18","PRDM1","ZNF135"), "TF", "Target"),
  Function = node_func[all_nodes],
  Stage = ifelse(all_nodes %in% c("ZNF135","Atoh1","Cer1","Egr3","Fst","Lox","Pdgfa","Sycp1","Kcnc4","Sox21","Thbs2"), "8-cell", "ICM"),
  row.names = NULL
)
nodes$Function[is.na(nodes$Function)] <- ""

# 3. 读入GO并构建基因→通路关系 ------------------------------------
icm_go   <- read.csv("ICM_H3K27me3_promoter_TSS_GO.csv")
icm_sig  <- subset(icm_go, p.adjust < 0.05)
eight_go <- read.csv("8_cell_H3K27me3_promoter_TSS_GO.csv")
eight_sig <- subset(eight_go, p.adjust < 0.05)

icm_genes   <- setNames(strsplit(icm_sig$geneID, "/"),   icm_sig$ID)
eight_genes <- setNames(strsplit(eight_sig$geneID, "/"), eight_sig$ID)

icm_selected   <- c("GO:0021954","GO:0007389","GO:0048562","GO:0030182","GO:0009952")
eight_selected <- c("GO:0007389","GO:0048562","GO:0035282","GO:0001822")

make_pathway_edges <- function(gene_list, go_list, selected_go, go_sig_df, stage_name) {
  res <- data.frame(Source=character(), Target=character(), Evidence=character(), Stage=character(), Type=character())
  for(goid in selected_go){
    gg <- go_list[[goid]]
    genes_in_network <- intersect(gene_list, gg)
    if(length(genes_in_network) > 0){
      go_name <- unique(go_sig_df$Description[go_sig_df$ID == goid])
      res <- rbind(res, data.frame(
        Source = genes_in_network,
        Target = go_name,
        Evidence = "functional annotation",
        Stage = stage_name,
        Type = "pathway"
      ))
    }
  }
  return(res)
}

target_icm   <- unique(edges$Target[edges$Stage == "ICM"])
target_8cell <- unique(edges$Target[edges$Stage == "8-cell"])

pathway_icm   <- make_pathway_edges(target_icm,   icm_genes,   icm_selected,   icm_sig,   "ICM")
pathway_8cell <- make_pathway_edges(target_8cell, eight_genes, eight_selected, eight_sig, "8-cell")

# 4. 合并所有边（修正列匹配） ---------------------------------------
motif_edges <- data.frame(edges,
  Type = ifelse(edges$Evidence %in% c("IGV co-localization (Otx2)","IGV co-localization (Neurod1)"), "IGV", "Motif"))
all_edges <- rbind(motif_edges, pathway_icm, pathway_8cell)

# 5. 更新节点表（加入通路节点） ------------------------------------
pathway_nodes <- data.frame(
  Node     = unique(c(pathway_icm$Target, pathway_8cell$Target)),
  Type     = "GO_BP",
  Function = "GO Biological Process",
  Stage    = NA,
  row.names = NULL
)
final_nodes <- rbind(nodes, pathway_nodes)

# 6. 保存文件 ------------------------------------------------------
write.table(all_edges,  "network_edges_all.txt",  sep="\t", quote=FALSE, row.names=FALSE)
write.table(final_nodes, "network_nodes_all.txt",  sep="\t", quote=FALSE, row.names=FALSE)
cat("网络文件已保存。\n")

# 7. 绘制预览图 ----------------------------------------------------
g <- graph_from_data_frame(all_edges, vertices = final_nodes$Node, directed = TRUE)
V(g)$color <- ifelse(final_nodes$Type == "TF", "orange", 
                     ifelse(final_nodes$Type == "GO_BP", "lightyellow", "skyblue"))
V(g)$label.cex <- 0.5
E(g)$color <- ifelse(all_edges$Type == "IGV", "red", 
                     ifelse(all_edges$Evidence == "functional annotation", "gray70", "gray50"))

pdf("Regulatory_network_ICM_8cell.pdf", width = 12, height = 10)
plot(g, layout = layout_with_fr, vertex.size = 10, edge.arrow.size = 0.3,
     main = "TF-target gene-GO regulatory network (ICM & 8-cell)")
dev.off()
cat("网络图已保存。\n")
