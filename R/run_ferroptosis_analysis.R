#' Run Ferroptosis Analysis for Multiple Samples
#'
#' Executes the complete ferroptosis analysis pipeline for multiple samples, including:
#' - Differential expression analysis
#' - PPI network construction
#' - Functional module detection
#' - Ferroptosis score calculation
#'
#' @param counts_data A matrix or data frame of raw counts with genes as rows and samples as columns.
#' @param meta_data Data frame containing sample metadata (must contain 'sample_id' and 'group' columns).
#' @param control_group Name of the control group in meta_data (default = "control").
#' @param species Character string specifying the species ("mouse" or "human", default = "mouse").
#' @param ppi_file Path to the PPI network file (STRING database format). If NULL, uses built-in data.
#' @param output_dir Directory to save results (default = tempdir()).
#' @param save_intermediate Logical indicating whether to save intermediate results (default = FALSE).
#' @param parallel Logical indicating whether to use parallel processing (default = FALSE).
#' @param n_cores Number of cores for parallel processing (default = 2).
#'
#' @return A list containing:
#' \itemize{
#'   \item results: List of detailed analysis results for each sample
#'   \item summary: Data frame of ferroptosis scores for all samples
#'   \item merged_DEG: Combined differential expression results
#' }
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' results <- run_ferroptosis_analysis(
#'   counts_data = counts_matrix,
#'   meta_data = sample_metadata,
#'   control_group = "control",
#'   species = "mouse"
#' )
#' }
#' @export
run_ferroptosis_analysis <- function(counts_data,
                                     meta_data,
                                     control_group = "control",
                                     species = "mouse",
                                     ppi_file = NULL,
                                     output_dir = tempdir(),
                                     save_intermediate = FALSE,
                                     parallel = FALSE,
                                     n_cores = 2) {


  #Input validation
  validate_inputs(counts_data, meta_data, control_group, species)

  #Prepare sample information
  sample_info <- prepare_sample_info(meta_data, control_group)
  control_cols <- which(sample_info$is_control)
  treated_samples <- sample_info[!sample_info$is_control, ]

  #Get PPI file
  if (is.null(ppi_file)) {
  ppi_file <- get_ppi_file(species)
  }

  #Set up parallel processing
  if (parallel) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      stop("Package 'parallel' needed for parallel processing. Please install it.")
    }
    n_cores <- min(n_cores, parallel::detectCores() - 1)
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl))
  }

  #Initialize result storage
  analysis_results <- list()
  all_deg <- list()

  #Process each sample
  for (i in seq_len(nrow(treated_samples))) {
    sample_id <- treated_samples$sample_id[i]
    sample_col <- which(colnames(counts_data) == sample_id)

    if (length(sample_col) == 0) {
      warning(paste("Sample", sample_id, "not found in counts_data - skipping"))
      next
    }

    cat("Processing sample:", sample_id, "...\n")

    #Calculate ferroptosis score
    result <- tryCatch({
      calculate_ferroptosis_score(
        counts_data = counts_data,
        sample_name = sample_id,
        control_cols = control_cols,
        treated_cols = sample_col,
        species = species,
        ppi_file = ppi_file
      )
    }, error = function(e) {
      warning(paste("Error processing sample", sample_id, ":", e$message))
      NULL
    })

    if (!is.null(result)) {
      analysis_results[[sample_id]] <- result

      #Store differential expression results
      deg_df <- result$DEG_results
      deg_df <- cbind(Gene = rownames(deg_df), deg_df)
      rownames(deg_df) <- NULL
      deg_df$Sample <- sample_id
      all_deg[[sample_id]] <- deg_df

      #Save intermediate results
      if (save_intermediate) {
        save_sample_results(result, output_dir, sample_id)
      }
    }
  }

  if (length(analysis_results) == 0) {
    stop("No samples were processed successfully")
  }

  #Create summary
  summary_table <- do.call(rbind, lapply(analysis_results, function(x) {
    data.frame(
      Sample = x$sample_name,
      Iron_metabolism = x$iron_metabolism,
      Lipid_metabolism = x$lipid_metabolism,
      GSH_metabolism = x$gsh_metabolism,
      Ferroptosis_score = x$ferroptosis_score,
      stringsAsFactors = FALSE
    )
  }))

  #Merge all differential expression results
  merged_deg <- do.call(rbind, all_deg)

  #Save summary results
  if (save_intermediate) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    utils::write.csv(merged_deg,
                     file.path(output_dir, "merged_DEG_results.csv"),
                     row.names = FALSE)
    utils::write.csv(summary_table,
                     file.path(output_dir, "ferroptosis_scores_summary.csv"),
                     row.names = FALSE)
  }

  #Return final results
  list(
    results = analysis_results,
    summary = summary_table,
    merged_DEG = merged_deg
  )
}

# Internal helper functions ------------------------------------------------------------

#' Validate input data
#' @keywords internal
validate_inputs <- function(counts_data, meta_data, control_group, species) {
  if (!is.data.frame(counts_data) && !is.matrix(counts_data)) {
    stop("counts_data must be a data frame or matrix")
  }

  if (!is.data.frame(meta_data)) {
    stop("meta_data must be a data frame")
  }

  required_cols <- c("sample_id", "group")
  if (!all(required_cols %in% colnames(meta_data))) {
    stop(paste("meta_data must contain columns:", paste(required_cols, collapse = ", ")))
  }

  if (!control_group %in% meta_data$group) {
    stop(paste("control_group", control_group, "not found in meta_data$group"))
  }

  if (!species %in% c("mouse", "human")) {
    stop("species must be either 'mouse' or 'human'")
  }

  if (ncol(counts_data) != nrow(meta_data)) {
    stop("Number of columns in counts_data must match number of rows in meta_data")
  }
}


#' Prepare sample information
#' @keywords internal
prepare_sample_info <- function(meta_data, control_group) {
  meta_data %>%
    dplyr::mutate(
      is_control = (group == control_group),
      sample_id = as.character(sample_id)
    ) %>%
    dplyr::arrange(desc(is_control)) %>%  # Controls first
    dplyr::select(sample_id, is_control)
}


#' Get PPI file path
#' @keywords internal
get_ppi_file <- function(species) {
  if (species == "mouse") {
    ppi_file <- system.file("extdata", "10090.protein.links.v11.5.txt.gz", package = "FerroScore")
  } else {
    ppi_file <- system.file("extdata", "9606.protein.links.v11.5.txt.gz", package = "FerroScore")
  }

  if (!file.exists(ppi_file)) {
    stop(paste("Built-in PPI file not found for species", species,
               "\nPlease provide ppi_file parameter or install the package with the required data files."))
  }

  return(ppi_file)
}


#' Save sample-level results
#' @keywords internal
save_sample_results <- function(result, output_dir, sample_id) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Save DEG results
  deg_file <- file.path(output_dir, paste0(sample_id, "_DEG.csv"))
  utils::write.csv(result$DEG_results, deg_file, row.names = FALSE)

  # Save scores
  scores <- data.frame(
    sample = sample_id,
    iron_score = result$iron_metabolism,
    lipid_score = result$lipid_metabolism,
    gsh_score = result$gsh_metabolism,
    ferroptosis_score = result$ferroptosis_score
  )
  score_file <- file.path(output_dir, paste0(sample_id, "_scores.csv"))
  utils::write.csv(scores, score_file, row.names = FALSE)
}


