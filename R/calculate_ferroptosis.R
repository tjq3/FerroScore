#' Calculate Ferroptosis Score for Single Samples
#'
#' This function performs differential expression analysis, constructs protein-protein interaction (PPI) networks,
#' identifies functional modules, and calculates a ferroptosis score based on iron, lipid, and glutathione metabolism.
#'
#' @param counts_data A matrix or data frame of raw counts with genes as rows and samples as columns.
#' @param sample_name Character string specifying the sample name being analyzed.
#' @param control_cols Numeric vector specifying the columns of control samples.
#' @param treated_cols Numeric vector specifying the columns of treated samples.
#' @param species Character string specifying the species, either "mouse" or "human".
#' @param ppi_file Path to the PPI network file (STRING database format).
#'
#'
#' @return A list containing:
#' \itemize{
#'   \item sample_name: The analyzed sample name
#'   \item DEG_results: Data frame of differential expression results
#'   \item iron_metabolism: Score for iron-related modules
#'   \item lipid_metabolism: Score for lipid-related modules
#'   \item gsh_metabolism: Score for glutathione-related modules
#'   \item ferroptosis_score: Composite ferroptosis score
#' }
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' result <- calculate_ferroptosis_score(
#'   counts_data = counts_matrix,
#'   sample_name = "Treatment1",
#'   control_cols = 1:3,
#'   treated_cols = 4:6,
#'   species = "human",
#'   ppi_file = "path/to/STRINGdb_network.txt"
#' )
#' }
#'
#' @export
calculate_ferroptosis_score <- function(counts_data, sample_name, control_cols, treated_cols,
                                        species = "mouse", ppi_file = NULL) {


  if (!species %in% c("mouse", "human")) {
    stop("Species must be either 'mouse' or 'human'")
  }
  if (species == "mouse") {
    org_db <- org.Mm.eg.db
    kegg_org <- "mmu"
    string_species <- 10090
  } else {
    org_db <- org.Hs.eg.db
    kegg_org <- "hsa"
    string_species <- 9606
  }


  #Differential Expression Analysis
  group <- factor(c(rep("control", length(control_cols)), rep("treated", length(treated_cols))))
  dge <- DGEList(counts = counts_data[, c(control_cols, treated_cols)], group = group)

  keep <- rowSums(cpm(dge) > 100) >= 1
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge$samples$lib.size <- colSums(dge$counts)

  dge <- calcNormFactors(dge, method = 'TMM')

  dge <- estimateDisp(dge)
  et <- exactTest(dge)

  DEG_results <- topTags(et, n = nrow(et$table))$table %>% na.omit()
  DEG_results$change <- ifelse(DEG_results$PValue < 0.05 & abs(DEG_results$logFC) > 1,
                               ifelse(DEG_results$logFC > 1, "Up", "Down"), "Not Change")

  diff_genes <- filter(DEG_results, change != "Not Change")
  if (nrow(diff_genes) < 5) {
    warning(paste("Insufficient differential genes for", sample_name, "- skipping"))
    return(list(
      sample_name = sample_name,
      DEG_results = DEG_results,
      iron_metabolism = NA,
      lipid_metabolism = NA,
      gsh_metabolism = NA,
      ferroptosis_score = NA
    ))
  }


  #Protein–Protein Interaction Network Construction
  diff_genes <- filter(DEG_results, change != "Not Change")
  diff_genes <- rownames(diff_genes)

  standardize_gene_names <- function(genes, species) {
    if (species == "mouse") {
      str_to_title(genes)
    } else {
      toupper(genes)
    }
  }
  diff_genes_std <- standardize_gene_names(diff_genes, species)

  string_db <- STRINGdb$new(version="11.5", species=string_species, score_threshold=400, input_directory="", protocol="http")
  mapped_genes <- string_db$map(data.frame(gene = diff_genes_std), "gene", removeUnmappedRows = TRUE)

  if (is.null(ppi_file)) {
    stop("PPI network file path must be provided")
  }

  ppi_data <- read.table(ppi_file, header=TRUE, stringsAsFactors=FALSE)
  ppi_network <- ppi_data[ppi_data$protein1 %in% mapped_genes$STRING_id &
                            ppi_data$protein2 %in% mapped_genes$STRING_id, ]
  ppi_graph <- graph_from_data_frame(ppi_network[, c("protein1", "protein2")], directed = FALSE)
  gene_name_map <- setNames(mapped_genes$gene, mapped_genes$STRING_id)
  V(ppi_graph)$name <- gene_name_map[V(ppi_graph)$name]


  #Module Identification
  set.seed(123)
  leiden_clusters <- cluster_leiden(ppi_graph, objective_function = "modularity", resolution = 1, n_iterations = 1000000)
  module_membership <- membership(leiden_clusters)
  V(ppi_graph)$module <- module_membership
  module_genes <- split(names(module_membership), module_membership)

  if (length(module_genes) == 0) {
    warning(paste("No modules found for", sample_name, "- skipping"))
    return(list(
      sample_name = sample_name,
      DEG_results = DEG_results,
      iron_metabolism = NA,
      lipid_metabolism = NA,
      gsh_metabolism = NA,
      ferroptosis_score = NA
    ))
  }


  #Module Function Identification
  enrichment_results <- list()
  for (i in 1:length(module_genes)) {
    genes <- module_genes[[i]]
    genes_std <- standardize_gene_names(genes, species)
    gene_entrez <- bitr(genes_std, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_db)

    go_enrich <- enrichGO(gene = gene_entrez$ENTREZID, OrgDb = org_db,
                          keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
                          pvalueCutoff = 0.05, qvalueCutoff = 0.05)
    kegg_enrich <- enrichKEGG(gene = gene_entrez$ENTREZID, organism = kegg_org,
                              pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
    enrichment_results[[paste0("Module_", i)]] <- list(GO = go_enrich, KEGG = kegg_enrich)
  }

  for (i in 1:length(enrichment_results)) {
    go_enrich <- enrichment_results[[i]]$GO
    kegg_enrich <- enrichment_results[[i]]$KEGG

    go_significant <- as.data.frame(go_enrich)[as.data.frame(go_enrich)$p.adjust < 0.05, ]
    kegg_significant <- as.data.frame(kegg_enrich)[as.data.frame(kegg_enrich)$p.adjust < 0.05, ]

    if (nrow(go_significant) == 0 && nrow(kegg_significant) == 0) {
      genes <- module_genes[[i]]
      genes_std <- standardize_gene_names(genes, species)
      gene_entrez <- bitr(genes_std, fromType = "SYMBOL",
                          toType = "ENTREZID", OrgDb = org_db)

      go_enrich_all <- enrichGO(gene = gene_entrez$ENTREZID,
                                OrgDb = org_db,
                                keyType = "ENTREZID",
                                ont = "ALL",
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05)
      enrichment_results[[paste0("Module_", i)]]$GO_ALL <- go_enrich_all
    }
  }

  wbg <- createWorkbook()
  for (i in 1:length(enrichment_results)) {
    go_enrich <- enrichment_results[[i]]$GO
    go_df <- as.data.frame(go_enrich)

    sheet_name <- paste0("GO_Module_", i)
    addWorksheet(wbg, sheetName = sheet_name)
    writeData(wbg, sheet = sheet_name, x = go_df, rowNames = FALSE)
  }
  wbk <- createWorkbook()
  for (i in 1:length(enrichment_results)) {
    kegg_enrich <- enrichment_results[[i]]$KEGG
    kegg_df <- as.data.frame(kegg_enrich)

    sheet_name <- paste0("KEGG_Module_", i)
    addWorksheet(wbk, sheetName = sheet_name)
    writeData(wbk, sheet = sheet_name, x = kegg_df, rowNames = FALSE)
  }

  module_functions <- list()

  for (i in 1:length(enrichment_results)) {
    go_enrich <- enrichment_results[[i]]$GO
    kegg_enrich <- enrichment_results[[i]]$KEGG
    go_enrich_all <- enrichment_results[[i]]$GO_ALL

    go_significant <- as.data.frame(go_enrich)[as.data.frame(go_enrich)$p.adjust < 0.05, ]
    kegg_significant <- as.data.frame(kegg_enrich)[as.data.frame(kegg_enrich)$p.adjust < 0.05, ]

    if (nrow(go_significant) > 0 || nrow(kegg_significant) > 0) {
      go_terms <- paste(go_significant$Description, collapse = "; ")
      kegg_pathways <- paste(kegg_significant$Description, collapse = "; ")
      module_summary <- paste(
        "Module", i, "_BP:\n",
        "GO Terms: ", go_terms, "\n",
        "KEGG Pathways: ", kegg_pathways, "\n"
      )
    } else if (!is.null(go_enrich_all)) {
      go_all_significant <- as.data.frame(go_enrich_all)[as.data.frame(go_enrich_all)$p.adjust < 0.05, ]
      if (nrow(go_all_significant) > 0) {
        go_all_terms <- paste(go_all_significant$Description, collapse = "; ")
        module_summary <- paste(
          "Module", i, "_ALL:\n",
          "GO Terms: ", go_all_terms, "\n"
        )
      } else {
        module_summary <- paste(
          "Module", i, "_None:\n",
          "No significant GO terms or KEGG pathways found.\n"
        )
      }
    } else {
      module_summary <- paste(
        "Module", i, "_None:\n",
        "No significant GO terms or KEGG pathways found.\n"
      )
    }

    module_functions[[paste0("Module_", i)]] <- module_summary
  }


  #Module Annotation
  module_functions <- list()
  lipid_keywords <- c("phospholipid", "lipid", "glycerol", "glycerolipid", "unsaturated fatty acid",
                      "fatty acid", "membrane lipid", "sphingolipid", "triglyceride",
                      "lipoprotein", "fat", "cholesterol", "cholesterol metabolism",
                      "glycolipid", "wax", "bile acid", "steroid", "lipogenesis", "lipolysis", "phospholipase")

  gsh_keywords <- c("reactive oxygen species", "oxidative stress", "cysteine", "methionine", "glutamate", "glycine",
                    "glutathione", "GSH", "glutathione metabolism", "glutathione peroxidase",
                    "homocysteine")

  iron_keywords <- c("iron coordination entity transport", "iron ion transport",
                     "iron", "ferroptosis", "ferritin", "transferrin", "heme",
                     "hepcidin", "ferroportin", "iron-sulfur clusters", "dcytb", "ceruloplasmin", "lactoferrin")

  module_names <- vector("character", length(enrichment_results))

  for (i in seq_along(enrichment_results)) {
    module_enrichment <- enrichment_results[[i]]
    go_terms <- tolower(module_enrichment$GO$Description)
    kegg_terms <- tolower(module_enrichment$KEGG$Description)
    all_terms <- c(go_terms, kegg_terms)
    total_terms <- length(all_terms)

    calculate_ratio <- function(terms, keywords) {
      if(total_terms == 0) return(0)
      matches <- sum(sapply(keywords, function(keyword) {
        sum(grepl(paste0("\\b", keyword, "\\b"), all_terms))
      }))
      matches / total_terms
    }

    lipid_ratio <- calculate_ratio(all_terms, lipid_keywords)
    gsh_ratio <- calculate_ratio(all_terms, gsh_keywords)
    iron_ratio <- calculate_ratio(all_terms, iron_keywords)

    if(total_terms == 0) {
      module_names[i] <- paste("Module", i)
      next
    }

    qualified_categories <- c()
    if(lipid_ratio >= 0.04) qualified_categories <- c(qualified_categories, "lipid")
    if(gsh_ratio >= 0.08) qualified_categories <- c(qualified_categories, "gsh")
    if(iron_ratio >= 0.03) qualified_categories <- c(qualified_categories, "iron")

    if(length(qualified_categories) == 0) {
      if(length(go_terms) > 0) {
        module_names[i] <- module_enrichment$GO$Description[1]
      } else if(length(kegg_terms) > 0) {
        module_names[i] <- module_enrichment$KEGG$Description[1]
      } else {
        module_names[i] <- paste("Module", i)
      }
    } else if(length(qualified_categories) == 1) {
      module_names[i] <- switch(qualified_categories[1],
                                "lipid" = "Lipid metabolism",
                                "gsh" = "GSH metabolism",
                                "iron" = "Iron metabolism")
    } else {
      ratios <- c(lipid = lipid_ratio, gsh = gsh_ratio, iron = iron_ratio)
      selected_category <- names(which.max(ratios[qualified_categories]))
      module_names[i] <- switch(selected_category,
                                "lipid" = "Lipid metabolism",
                                "gsh" = "GSH metabolism",
                                "iron" = "Iron metabolism")
    }
  }


  #Ferroptosis Network Extraction
  num_modules <- max(module_membership)
  new_adj_matrix <- matrix(0, nrow = num_modules, ncol = num_modules)

  for (edge in E(ppi_graph)) {
    from_node <- ends(ppi_graph, edge)[1]
    to_node <- ends(ppi_graph, edge)[2]

    from_module <- module_membership[from_node]
    to_module <- module_membership[to_node]

    edge_weight <- 1

    if (from_module == to_module) {
      new_adj_matrix[from_module, from_module] <- new_adj_matrix[from_module, from_module] + edge_weight
    } else {
      new_adj_matrix[from_module, to_module] <- new_adj_matrix[from_module, to_module] + edge_weight
      new_adj_matrix[to_module, from_module] <- new_adj_matrix[to_module, from_module] + edge_weight
    }
  }

  new_graph <- graph_from_adjacency_matrix(new_adj_matrix, mode = "undirected", weighted = TRUE)
  custom_module_names <- module_names
  V(new_graph)$name <- custom_module_names

  ferroptosis_modules <- which(module_names %in% c("Lipid metabolism", "GSH metabolism", "Iron metabolism"))
  module_nodes <- ferroptosis_modules
  module_subgraph <- induced_subgraph(new_graph, module_nodes)

  if (length(ferroptosis_modules) == 0) {
    return(list(
      sample_name = sample_name,
      DEG_results = DEG_results,
      iron_metabolism = NA,
      lipid_metabolism = NA,
      gsh_metabolism = NA,
      ferroptosis_score = NA
    ))
  }


  #Topological Feature Extraction
  module_nodes <- ferroptosis_modules
  module_name <- module_names[ferroptosis_modules]
  module_subgraph <- induced_subgraph(new_graph, module_nodes)

  node_degree <- degree(module_subgraph)
  node_betweenness <- betweenness(module_subgraph)
  node_clustering <- transitivity(module_subgraph, type = "local")

  node_features <- data.frame(
    Node = V(module_subgraph)$name,
    Degree = node_degree,
    Betweenness = node_betweenness,
    Clustering = node_clustering
  )

  calculate_module_metrics <- function(ppi_graph, module_membership, c) {
    module_c_nodes <- names(module_membership[module_membership == c])
    m <- ecount(ppi_graph)
    Q_C <- 0

    for (i in module_c_nodes) {
      for (j in module_c_nodes) {
        if (i < j) {
          A_ij <- ifelse(are_adjacent(ppi_graph, i, j), 1, 0)
          k_i <- degree(ppi_graph, i)
          k_j <- degree(ppi_graph, j)
          Q_C <- Q_C + (A_ij - (k_i * k_j) / (2 * m))
        }
      }
    }
    Q_C <- Q_C / (2 * m)

    module_c_subgraph <- induced_subgraph(ppi_graph, module_c_nodes)
    E_C <- ecount(module_c_subgraph)
    N_C <- vcount(module_c_subgraph)
    intra_density <- (2 * E_C) / (N_C * (N_C - 1))

    return(list(
      Module = c,
      Modularity = Q_C,
      Intra_module_Density = intra_density
    ))
  }

  module_metrics <- lapply(ferroptosis_modules, function(c) {
    calculate_module_metrics(ppi_graph, module_membership, c)
  })

  module_features <- do.call(rbind, lapply(module_metrics, function(x) {
    data.frame(
      Module = x$Module,
      Modularity = x$Modularity,
      Intra_module_Density = x$Intra_module_Density
    )
  }))

  combined_features <- lapply(seq_along(module_nodes), function(i) {
    c <- module_nodes[i]
    module_node_name <- module_name[i]

    module_metrics <- module_features[module_features$Module == c, ]
    node_metrics <- node_features[node_features$Node == module_node_name, ]

    metrics <- c(
      node_metrics$Degree,
      node_metrics$Betweenness,
      node_metrics$Clustering,
      module_metrics$Modularity,
      module_metrics$Intra_module_Density
    )

    metrics[is.na(metrics)] <- 0
    combined_score <- mean(metrics)

    return(list(
      Module = c,
      Node = module_node_name,
      Combined_Score = combined_score
    ))
  })

  combined_features_df <- do.call(rbind, lapply(combined_features, function(x) {
    data.frame(
      Module = x$Module,
      Node = x$Node,
      Combined_Score = x$Combined_Score
    )
  }))


  #Ferroptosis Score Calculation
  metabolism_scores <- list(
    lipid_metabolism = combined_features_df$Combined_Score[combined_features_df$Node == "Lipid metabolism"],
    iron_metabolism = combined_features_df$Combined_Score[combined_features_df$Node == "Iron metabolism"],
    gsh_metabolism = combined_features_df$Combined_Score[combined_features_df$Node == "GSH metabolism"]
  )

  summarize_scores <- function(scores) {
    if (is.null(scores) || length(scores) == 0) {
      return(1)
    } else {
      return(mean(scores, na.rm = TRUE))
    }
  }

  lipid_metabolism <- summarize_scores(metabolism_scores$lipid_metabolism)
  iron_metabolism <- summarize_scores(metabolism_scores$iron_metabolism)
  gsh_metabolism <- summarize_scores(metabolism_scores$gsh_metabolism)

  if (gsh_metabolism == 0) gsh_metabolism <- 1
  ferroptosis_score <- (iron_metabolism * lipid_metabolism) / gsh_metabolism

  return(list(
    sample_name = sample_name,
    DEG_results = DEG_results,
    iron_metabolism = iron_metabolism,
    lipid_metabolism = lipid_metabolism,
    gsh_metabolism = gsh_metabolism,
    ferroptosis_score = ferroptosis_score
  ))
}




