getEnrichment <- function(df, tag){
  ###########
  dirEnrichDirection <- paste0(dirEnrich, tag, "/")
  
  if (!dir.exists(dirEnrichDirection)){
    dir.create(dirEnrichDirection)
  }else{
    print(paste("The directory", dirEnrichDirection, "already exists"))
  }
  ##############
  
  # TODO add dbs
  
  TF <- df$TRANSFAC_and_JASPAR_PWMs
  
  CC <- df$GO_Cellular_Component_2023
  BP <- df$GO_Biological_Process_2023
  MF <- df$GO_Molecular_Function_2023
  
  KEGG <- df$KEGG_2021_Human
  Reactome <- df$Reactome_2022
  
  DisGeNET <- df$DisGeNET
  GWAS <- df$GWAS_Catalog_2023
  
  
  getEnrichmentPlot(CC, "CC", top_term, thr_pval, dirEnrichDirection)
  getEnrichmentPlot(BP, "BP", top_term, thr_pval, dirEnrichDirection)
  getEnrichmentPlot(MF, "MF", top_term, thr_pval, dirEnrichDirection)
  
  getEnrichmentPlot(KEGG, "KEGG", top_term, thr_pval, dirEnrichDirection)
  getEnrichmentPlot(Reactome, "Reactome", top_term, thr_pval, dirEnrichDirection)
  
  getEnrichmentPlot(DisGeNET, "DisGeNET", top_term, thr_pval, dirEnrichDirection)
  getEnrichmentPlot(GWAS, "GWAS", top_term, thr_pval, dirEnrichDirection)
  
  getEnrichmentPlot(TF, "TF", top_term, thr_pval, dirEnrichDirection)
  
}

