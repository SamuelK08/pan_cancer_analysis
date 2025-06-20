# run_enrichment.R

#--- load libraries ---
#use suppressPackageStartupMessages to keep the output clean
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(argparse))

#--- main function ---
run_pathway_enrichment <- function(modules_file, gmt_folder, output_folder) {
  print("--- step 2: pathway enrichment analysis ---")
  
  #--- 1. load pathways ---
  print("loading pathway gene sets (.gmt files)...")
  gmt_files <- list.files(gmt_folder, pattern = "*.gmt", full.names = TRUE)
  if (length(gmt_files) == 0) {
    stop("no .gmt files found in the specified folder.")
  }
  
  all_pathways <- do.call(c, lapply(gmt_files, gmtPathways))
  print(paste("loaded", length(all_pathways), "total pathways."))
  
  #--- 2. get universe of all genes ---
  print("defining the universe of all protein-coding genes from ensembl...")
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  genes <- getBM(attributes = "entrezgene_id",
                 filters = "biotype",
                 values = "protein_coding",
                 mart = mart)
  all_genes_universe <- as.character(genes[['entrezgene_id']])
  print(paste("gene universe defined with", length(all_genes_universe), "genes."))
  
  #--- 3. load modules ---
  print(paste("Loading modules from:", modules_file))
  modules_raw <- readLines(modules_file)
  if (length(modules_raw) == 0) {
      print("Warning: Modules file is empty. Nothing to process.")
      return() # Exit gracefully
  }
  module_list <- lapply(modules_raw, function(line) str_split(line, ", ")[[1]])
  
  #create the output directory if it doesn't exist
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  
  #--- 4. loop through each module and perform enrichment ---
  for (i in seq_along(module_list)) {
    print(paste("analyzing module #", i, "of", length(module_list)))
    
    #extract gene ids (non-mirna) from the module
    module_nodes <- module_list[[i]]
    module_genes <- module_nodes[!grepl("hsa-", module_nodes)] # Using grepl is slightly more robust in R
    
    #initialize results vectors
    pathway_names_res <- c()
    p_values <- c()
    odds_ratios <- c()
    fold_enrichments <- c()
    shared_genes_list <- c()
    jaccard_indices <- c() # Adding Jaccard as per previous advice
    
    #--- 5. loop through each pathway for the current module ---
    for (pathway_name in names(all_pathways)) {
      pathway_genes <- intersect(all_pathways[[pathway_name]], all_genes_universe)
      
      shared_genes <- intersect(module_genes, pathway_genes)
      
      if (length(shared_genes) > 0) {
        a <- length(shared_genes)
        b <- length(setdiff(pathway_genes, module_genes))
        c <- length(setdiff(module_genes, pathway_genes))
        d <- length(all_genes_universe) - (a + b + c)
        
        contingency_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
        fisher_result <- fisher.test(contingency_table, alternative = "greater")
        FE <- (a * length(all_genes_universe)) / (length(module_genes) * length(pathway_genes))
        jaccard_index <- a / (a + b + c)

        pathway_names_res <- c(pathway_names_res, pathway_name)
        p_values <- c(p_values, fisher_result$p.value)
        odds_ratios <- c(odds_ratios, fisher_result$estimate)
        fold_enrichments <- c(fold_enrichments, FE)
        shared_genes_list <- c(shared_genes_list, paste(shared_genes, collapse=", "))
        jaccard_indices <- c(jaccard_indices, jaccard_index)
      }
    }
    
    #--- 6. adjust p-values and save results ---
    if (length(p_values) > 0) {
      fdr_values <- p.adjust(p_values, method = "BH")
      
      results_df <- data.frame(
        Term = pathway_names_res,
        Jaccard = jaccard_indices,
        Odds_Ratio = odds_ratios,
        Fold_Enrichment = fold_enrichments,
        PValue = p_values,
        FDR = fdr_values,
        Genes = shared_genes_list
      )
      
      results_df <- results_df %>% arrange(desc(Fold_Enrichment))
      
      output_file_path <- file.path(output_folder, paste0("BRCA_Module_", i, "_pathwayAll_FisherResults.csv"))
      write.csv(results_df, output_file_path, row.names = FALSE)
      print(paste("saved results for module #", i, "to", output_file_path))
    } else {
      print(paste("no overlapping pathways found for module #", i))
    }
  }
  print("--- step 2 complete ---")
}

# --- execute the script with command-line arguments ---
parser <- ArgumentParser(description="Run pathway enrichment analysis on network modules.")
parser$add_argument("--modules", required=TRUE, help="Path to the input modules file.")
parser$add_argument("--gmt", required=TRUE, help="Path to the folder containing .gmt files.")
parser$add_argument("--output", required=TRUE, help="Path to the output folder for saving results.")
args <- parser$parse_args()

run_pathway_enrichment(
  modules_file = args$modules,
  gmt_folder = args$gmt,
  output_folder = args$output
)