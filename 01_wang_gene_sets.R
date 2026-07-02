# Marker gene sets used to calculate expression-enrichment scores.
# These scores are not deconvolved cell-type proportions.

GENE_SETS <- list(
  TUMOR_EPITHELIAL = c(
    "EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "MUC1", "TACSTD2"
  ),
  TNBC_BASAL = c(
    "KRT5", "KRT14", "KRT17", "EGFR", "TP63", "SOX10", "MIA"
  ),
  T_CELL = c(
    "CD3D", "CD3E", "CD2", "TRAC", "IL7R", "LCK"
  ),
  CD8_T_CYTOTOXIC = c(
    "CD8A", "CD8B", "GZMB", "GZMH", "NKG7", "PRF1", "CXCR3", "GNLY"
  ),
  B_CELL = c(
    "MS4A1", "CD79A", "CD79B", "BANK1", "CD74", "HLA-DRA"
  ),
  PLASMA = c(
    "MZB1", "JCHAIN", "XBP1", "IGHG1", "IGKC", "SDC1"
  ),
  MACROPHAGE = c(
    "CD68", "LST1", "C1QA", "C1QB", "C1QC", "TYROBP", "FCGR3A", "AIF1"
  ),
  M2_MACROPHAGE = c(
    "CD163", "MRC1", "MSR1", "MARCO", "IL10"
  ),
  DC_ANTIGEN_PRESENTATION = c(
    "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "B2M"
  ),
  FIBROBLAST = c(
    "COL1A1", "COL1A2", "DCN", "LUM", "COL3A1", "ACTA2", "TAGLN"
  ),
  ENDOTHELIAL = c(
    "PECAM1", "VWF", "KDR", "ENG", "CDH5", "ESAM"
  ),
  IFN_RESPONSE = c(
    "IFIT1", "IFIT2", "IFIT3", "ISG15", "MX1", "OAS1", "STAT1", "IRF1"
  ),
  CHEMOKINE_CCL = c(
    "CCL2", "CCL3", "CCL4", "CCL5", "CCL18", "CCL19", "CCL20", "CCL21"
  ),
  CHEMOKINE_CXCL = c(
    "CXCL8", "CXCL9", "CXCL10", "CXCL11", "CXCL12", "CXCL13", "CXCL16"
  ),
  MIF_CD74_AXIS = c(
    "MIF", "CD74", "CXCR4", "CXCR2", "CXCR7"
  ),
  TLS_LIKE = c(
    "MS4A1", "CD79A", "CD3D", "CD3E", "CXCL13", "CCL19", "CCL21",
    "LAMP3", "CCR7", "CD40", "LTB"
  ),
  HYPOXIA = c(
    "VEGFA", "CA9", "SLC2A1", "ENO1", "LDHA", "BNIP3", "NDRG1"
  ),
  EMT = c(
    "VIM", "FN1", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "TWIST1", "MMP2"
  ),
  TLS_WANG_30 = c(
    "CD79A", "CD79B", "TNFRSF13C", "BLK", "CD22", "CD37", "MS4A1",
    "NIBAN3", "CD19", "IKZF3", "LINC00926", "CXCR5", "LTB", "SELL",
    "CCL19", "POU2AF1", "CXCL13", "FCRLA", "VPREB3", "FCMR", "AL928768.3",
    "RASGRP2", "TCF7", "RIPOR2", "RAC2", "IL16", "CCR7", "CD52", "TCL1A", "ATP2A3"
  )
)
