#### Extract meta data for RNA-seq files using ARCHS4 ####
### [Optional]: Can also extract meta data for RNA-seq studies using this
#### For connecting python with R ####
### library(reticulate)
### use_python("/env/bin/python3", required = TRUE)
### py_config()
### py_module_available("archs4py")
### a4 <- import("archs4py")
### file = "/app/data/ARCHS4/human_gene_v2.5.h5"
### a4$ls(file)

# pData_RNA <- lapply(rna_studies, function( geo ) {
#   file = "/app/data/ARCHS4/human_gene_v2.5.h5"
#   a4$meta$series(file, geo)
# })

#### Generating expression matrices ####
generate_exprs_mtx <- function(DNA = NULL, RNA = NULL, dna_studies = list(), rna_mtx = list(), pData_RNA = list()) {
  
  #### Handling duplicated genes ####
  handling_duplicates <- function(data_table) {
    if (is.list(data_table$Genes)) {
      data_table[, Genes := vapply(Genes, function(x) {
        if (is.null(x) || length(x) == 0 || all(is.na(x))) {
          return(NA_character_)
        } else if (is.list(x)) {
          return(as.character(unlist(x)))
        } else {
          return(as.character(x))
        }
      }, FUN.VALUE = character(1))]
    }
    
    data_table <- na.omit(data_table[(Genes != "" & !is.na(Genes)), ])
    
    if (any(duplicated(data_table$Genes))) {
      exprs_data.table <- data_table[, lapply(.SD, mean), by = Genes, .SDcols = !'Genes']
    } else {
      exprs_data.table <- data_table
    }
    
    exprs_mtx <- as.matrix(exprs_data.table[, -1])
    rownames(exprs_mtx) <- exprs_data.table$Genes
    
    return(exprs_mtx)
  }
  
  #### For DNA microarray data sets ####
  if (!is.null(DNA) && DNA == TRUE) {
    if (!is.list(dna_studies))
      stop("dna_studies needs to be a list")
    
    if (!all(sapply(dna_studies, is.character)))
      stop("Elements in dna_studies must be character GEO IDs, e.g. list(study1 = 'GSEXXXX')")
    
    list_of_geo_data <- lapply(dna_studies, function(file) {
      tryCatch({
        gse_obj <- GEOquery::getGEO(GEO = file, AnnotGPL = TRUE)
        if (inherits(gse_obj, "ExpressionSet")) {
          list(gse_obj)
        } else {
          gse_obj
        }
      }, error = function(e) {
        warning(paste("Error retrieving", file, ":", e$message))
        NULL
      })
    })
    
    names(list_of_geo_data) <- names(dna_studies)
    
    out_dna <- vector("list", length(dna_studies))
    names(out_dna) <- names(dna_studies)
    
    for (s in seq_along(list_of_geo_data)) {
      es_list <- list_of_geo_data[[s]]
      if (is.null(es_list)) next
      
      per_platform_mtx <- vector("list", length(es_list))
      
      for (k in seq_along(es_list)) {
        geo <- es_list[[k]]
        if (is.null(geo)) next
        
        mtx <- Biobase::exprs(geo)
        fData <- Biobase::fData(geo)
        
        gene_cols <- c("Gene symbol", "Gene Symbol", "Symbol", "ILMN_Gene", "symbol", "GB_ACC")
        present <- intersect(gene_cols, colnames(fData))
        
        if (length(present) == 0L) {
          warning(sprintf("No gene annotation columns for study %s platform %d; using rownames.",
                          names(dna_studies)[s], k))
          fData$Genes <- rownames(fData)
        } else {
          if ("Gene symbol" %in% present) {
            selected_col <- "Gene symbol"
          } else if ("Gene Symbol" %in% present) {
            selected_col <- "Gene Symbol"
          } else if ("Symbol" %in% present) {
            selected_col <- "Symbol"
          } else if ("ILMN_Gene" %in% present) {
            selected_col <- "ILMN_Gene"
          } else if ("symbol" %in% present) {
            selected_col <- "symbol"
          } else {
            selected_col <- "GB_ACC"
          }
          
          if (selected_col == "GB_ACC") {
            probe_ids <- as.character(fData[[selected_col]])
            
            probe_ids_base <- sub("\\.\\d+$", "", probe_ids)
            
            
            symbols <- tryCatch(
              AnnotationDbi::mapIds(
                x = org.Hs.eg.db,
                keys = probe_ids_base,
                column = "SYMBOL",
                keytype = "ACCNUM",
                multiVals = "first"
              ),
              error = function(e) {
                warning(sprintf(
                  "mapIds ACCNUM failed for study %s platform %d: %s. Using raw GB_ACC as Genes.",
                  names(dna_studies)[s], k, e$message
                ))
                rep(NA_character_, length(probe_ids))
              }
            )
            
            if (all(is.na(symbols))) {
              warning(sprintf(
                "GB_ACC in study %s platform %d does not map as ACCNUM; using raw GB_ACC as Genes.",
                names(dna_studies)[s], k
              ))
              fData$Genes <- probe_ids
            } else {
              fData$Genes <- symbols
            }
          } else {
            fData$Genes <- as.character(fData[[selected_col]])
          }
        }
        
        genes_vector <- as.character(fData$Genes)
        
        df <- data.frame(Genes = genes_vector, mtx, check.names = FALSE)
        dt <- data.table::as.data.table(df)
        
        mtx2 <- handling_duplicates(dt)
        
        gene_names <- rownames(mtx2)
        exprs_mtx <- apply(mtx2, 2, as.numeric)
        log_mtx <- log(exprs_mtx + 1)
        log_mtx[is.nan(log_mtx)] <- 0
        
        norm_mtx <- limma::normalizeBetweenArrays(log_mtx, method = "quantile")
        rownames(norm_mtx) <- gene_names
        per_platform_mtx[[k]] <- norm_mtx
      }
      
      per_platform_mtx <- Filter(Negate(is.null), per_platform_mtx)
      if (length(per_platform_mtx) == 0L) next
      
      if (length(per_platform_mtx) == 1L) {
        out_dna[[s]] <- per_platform_mtx[[1]]
      } else {
        gene_lists <- lapply(per_platform_mtx, rownames)
        common_genes <- Reduce(intersect, gene_lists)
        
        if (length(common_genes) == 0L) {
          warning(sprintf(
            "No overlapping genes found across platforms for study %s; returning first platform matrix.",
            names(dna_studies)[s]
          ))
          out_dna[[s]] <- per_platform_mtx[[1]]
        } else {
          per_platform_common <- lapply(per_platform_mtx, function(m) {
            m_common <- m[common_genes, , drop = FALSE]
            m_common[common_genes, , drop = FALSE]
          })
          combined_mtx <- do.call(cbind, per_platform_common)
          rownames(combined_mtx) <- common_genes
          out_dna[[s]] <- combined_mtx
        }
      }
      
      pdf(paste0("boxplot_", names(dna_studies)[s], ".pdf"), width = 8, height = 8)
      boxplot(out_dna[[s]], outline = FALSE,
              col = "skyblue",
              main = paste0("Log transformation/Quantile\n", names(dna_studies)[s]))
      dev.off()
    }
    return(out_dna)
  }
  
  
  #### For RNA-Seq data sets ####
  if (!is.null(RNA) && RNA == TRUE) {
    if (!is.list(rna_mtx))
      stop("rna_mtx needs to be a list of GEO accession IDs or GEO ExpressionSets")
    
    list_of_norm_RNA_mtx <- vector("list", length(rna_mtx))
    list_of_dds <- vector("list", length(rna_mtx))
    diagnostics <- vector("list", length(rna_mtx))
    names(list_of_norm_RNA_mtx) <- names(rna_mtx)
    names(list_of_dds) <- names(rna_mtx)
    names(diagnostics) <- names(rna_mtx)
    
    for (i in seq_along(rna_mtx)) {
      diag_i <- list(status = "ok", messages = character(0))
      gse_in <- rna_mtx[[i]]
      gse_obj <- NULL
      
      # Accept char GEO id or ExpressionSet or list of ExpressionSets
      if (is.character(gse_in)) {
        gse_try <- tryCatch(GEOquery::getGEO(GEO = gse_in, AnnotGPL = TRUE),
                            error = function(e) e)
        if (inherits(gse_try, "error")) {
          diag_i$status <- "error"
          diag_i$messages <- c(diag_i$messages, paste0("getGEO failed: ", gse_try$message))
          diagnostics[[i]] <- diag_i
          next
        }
        # getGEO may return a list (platforms) or an ExpressionSet
        if (is.list(gse_try) && length(gse_try) >= 1) {
          gse_obj <- gse_try[[1]]
        } else if (inherits(gse_try, "ExpressionSet")) {
          gse_obj <- gse_try
        } else {
          diag_i$status <- "error"
          diag_i$messages <- c(diag_i$messages, "getGEO returned unexpected object")
          diagnostics[[i]] <- diag_i
          next
        }
      } else if (inherits(gse_in, "ExpressionSet")) {
        gse_obj <- gse_in
      } else if (is.list(gse_in) && length(gse_in) >= 1 && inherits(gse_in[[1]], "ExpressionSet")) {
        gse_obj <- gse_in[[1]]
      } else {
        diag_i$status <- "error"
        diag_i$messages <- c(diag_i$messages, "Unsupported rna_mtx element (expect GEO id or ExpressionSet)")
        diagnostics[[i]] <- diag_i
        next
      }
      
      # extract matrix and pData safely
      mtx_raw <- tryCatch(Biobase::exprs(gse_obj), error = function(e) e)
      if (inherits(mtx_raw, "error")) {
        diag_i$status <- "error"
        diag_i$messages <- c(diag_i$messages, paste0("exprs() failed: ", mtx_raw$message))
        diagnostics[[i]] <- diag_i
        next
      }
      
      # Ensure matrix shape and coercion to numeric WITHOUT producing logical mode
      mtx <- tryCatch({
        mtx_mat <- as.matrix(mtx_raw)
        # Coerce column-wise to numeric, preserve dimnames
        rn <- rownames(mtx_mat)
        cn <- colnames(mtx_mat)
        mtx_num <- apply(mtx_mat, 2, function(col) as.numeric(as.character(col)))
        # If single column apply returns vector, coerce back to matrix
        if (is.null(dim(mtx_num))) mtx_num <- matrix(mtx_num, ncol = 1)
        dimnames(mtx_num) <- list(rn, cn)
        mtx_num
      }, warning = function(w) {
        # still try to recover
        mtx_mat <- as.matrix(mtx_raw)
        rn <- rownames(mtx_mat); cn <- colnames(mtx_mat)
        mtx_num <- matrix(as.numeric(as.character(mtx_mat)), nrow = nrow(mtx_mat), ncol = ncol(mtx_mat))
        dimnames(mtx_num) <- list(rn, cn)
        mtx_num
      }, error = function(e) e)
      
      if (inherits(mtx, "error")) {
        diag_i$status <- "error"
        diag_i$messages <- c(diag_i$messages, paste0("Failed to coerce expression to numeric: ", mtx$message))
        diagnostics[[i]] <- diag_i
        next
      }
      
      pData <- tryCatch(Biobase::pData(gse_obj), error = function(e) e)
      if (inherits(pData, "error")) {
        diag_i$status <- "error"
        diag_i$messages <- c(diag_i$messages, paste0("pData() failed: ", pData$message))
        diagnostics[[i]] <- diag_i
        next
      }
      
      # ensure condition column
      if (!("condition" %in% colnames(pData))) {
        if ("characteristics_ch1" %in% colnames(pData)) {
          pData$condition <- pData$characteristics_ch1
        } else {
          pData$condition <- factor(rep("unknown", nrow(pData)))
        }
      }
      
      # match samples by names present in data
      common_samples <- intersect(colnames(mtx), rownames(pData))
      if (length(common_samples) == 0) {
        diag_i$status <- "error"
        diag_i$messages <- c(diag_i$messages, "No common sample names between expr matrix and pData (check colnames/rownames).")
        diagnostics[[i]] <- diag_i
        next
      }
      
      mtx_fil <- mtx[, common_samples, drop = FALSE]
      pData_fil <- pData[common_samples, , drop = FALSE]
      
      # Final numeric enforcement (matrix of numeric values)
      mtx_fil <- matrix(as.numeric(mtx_fil), nrow = nrow(mtx_fil), ncol = ncol(mtx_fil),
                        dimnames = list(rownames(mtx_fil), colnames(mtx_fil)))
      
      # If values are not raw integer counts, round for DESeq2 (safe)
      is_integer_counts <- all(abs(mtx_fil - round(mtx_fil)) < 1e-6)
      if (!is_integer_counts) {
        count_matrix <- round(mtx_fil)
        diag_i$messages <- c(diag_i$messages, "Input non-integer; rounded to integers for DESeq2.")
      } else {
        count_matrix <- mtx_fil
      }
      
      # filter genes with zero counts
      keep_genes <- rowSums(count_matrix, na.rm = TRUE) > 0
      if (sum(keep_genes) == 0) {
        diag_i$status <- "error"
        diag_i$messages <- c(diag_i$messages, "No genes with non-zero counts after filtering.")
        diagnostics[[i]] <- diag_i
        next
      }
      count_matrix <- count_matrix[keep_genes, , drop = FALSE]
      
      # Build DESeq dataset safely
      dds_1 <- tryCatch(DESeqDataSetFromMatrix(count_matrix, colData = pData_fil, design = ~ condition),
                        error = function(e) e)
      if (inherits(dds_1, "error")) {
        diag_i$status <- "error"
        diag_i$messages <- c(diag_i$messages, paste0("DESeqDataSetFromMatrix failed: ", dds_1$message))
        diagnostics[[i]] <- diag_i
        next
      }
      
      dds_2 <- tryCatch(DESeq(dds_1), error = function(e) e)
      if (inherits(dds_2, "error")) {
        diag_i$status <- "error"
        diag_i$messages <- c(diag_i$messages, paste0("DESeq() failed: ", dds_2$message))
        diagnostics[[i]] <- diag_i
        next
      }
      
      norm_counts <- tryCatch(counts(dds_2, normalized = TRUE), error = function(e) e)
      if (inherits(norm_counts, "error")) {
        diag_i$status <- "error"
        diag_i$messages <- c(diag_i$messages, paste0("counts(..., normalized=TRUE) failed: ", norm_counts$message))
        diagnostics[[i]] <- diag_i
        next
      }
      
      # store results
      list_of_dds[[i]] <- dds_2
      list_of_norm_RNA_mtx[[i]] <- norm_counts
      diagnostics[[i]] <- diag_i
    } # end for
    
    # write boxplots for those that succeeded (named)
    for (i in seq_along(list_of_norm_RNA_mtx)) {
      if (is.null(list_of_norm_RNA_mtx[[i]]) || length(list_of_norm_RNA_mtx[[i]]) == 0) next
      pdf(paste0("boxplot_", ifelse(is.null(names(rna_mtx)[i]), paste0("rna_", i), names(rna_mtx)[i]), ".pdf"),
          width = 8, height = 8)
      boxplot(list_of_norm_RNA_mtx[[i]], outline = FALSE, col = "skyblue",
              main = paste0("Median of Ratios\n", ifelse(is.null(names(rna_mtx)[i]), paste0("rna_", i), names(rna_mtx)[i])))
      dev.off()
    }
    
    return(list(dds_objs = list_of_dds, norm_counts = list_of_norm_RNA_mtx, diagnostics = diagnostics))
  }
  
}

####
####

#### Finding common genes across studies #####
find_common_genes <- function(DNA = FALSE, RNA = FALSE, DNA_mtx = list(), RNA_mtx = list()) {
  if(DNA && !RNA) {
    return(Reduce(intersect, lapply(DNA_mtx, rownames)))
  }
  if(RNA && !DNA) {
    if(length(RNA_mtx) == 2) stop("Use RNA_mtx$norm_counts")
    return(Reduce(intersect, lapply(RNA_mtx, rownames)))
  }
  if(DNA && RNA) {
    list_of_matrices <- c(DNA_mtx, RNA_mtx)
    return(Reduce(intersect, lapply(list_of_matrices, rownames)))
  }
  return(NULL)
}

####
####

#### Meta-analysis function ####
meta_results <- function(list_of_studies) {

list_of_es <- lapply(list_of_studies, function( study ) {
  effect.sizes(study)
  })
  
  summary <- combine.effect.sizes(list_of_es)
  
  summary$g <- summary$g[!apply(is.na(summary$g) | is.infinite(summary$g), 1, any), ]
  summary$se.g <- summary$se.g[!apply(is.na(summary$se.g) | is.infinite(summary$se.g), 1, any), ]
  g_genes <- rownames(summary$g)
  
  summary$pooled.estimates$genes <- rownames(summary$pooled.estimates)
  summary$pooled.estimates <- summary$pooled.estimates %>%
    dplyr::filter(genes %in% g_genes)
  
  summary$pooled.estimates <- summary$pooled.estimates %>%
    dplyr::filter(!apply(is.na(.) | is.infinite(as.matrix(.)), 1, any))
  
  summary$pooled.estimates <- as.data.frame(lapply(summary$pooled.estimates, function(x) ifelse(x == 0, 1e-200, x)))
  rownames(summary$pooled.estimates) <- g_genes
  
  g <- summary$g
  se.g <- summary$se.g
  
  pool    <- summary$pooled.estimates[, "summary" ]
  names(pool) <- rownames(g)
  
  se.pool <- summary$pooled.estimates[, "se.summary" ]
  names(se.pool) <- rownames(g)
  
  x.label <- "Standardized Mean Difference (log2 scale)"
  
  consistent_genes <- c()
  
  for (gene in rownames(g)) {
    g_gene <- g[gene, ]
    if ((all(g_gene > 0) & pool[gene] > 0.5) | (all(g_gene < 0) & pool[gene] < -0.5)) {
      consistent_genes <- c(consistent_genes, gene)
    }
  }
  
  for (gene in consistent_genes) {
    g_gene <- g[gene, ]
    se.g_gene <- se.g[gene, ]
    study_names <- gsub("_g", "", names(g_gene))
    
    pdf(file = paste0(gene, ".pdf"), height = 3.5, width = 3.5)
    par(cex = 0.65)
    metaplot(g_gene, se.g_gene,
             labels = study_names,
             summn = pool[gene],
             sumse = se.pool[gene],
             sumnn = 1/se.pool[gene]^2,
             summlabel = "Summary effect",
             xlab = x.label,
             ylab = "",
             main = bquote(italic(.(gene))),
             colors = meta.colors(box = "violetred", lines = "plum", summary = "mediumpurple", 
                                  text = "black", axes = "black", zero = "black"),  
             boxsize = 1,
             lty.random = 1,  
             lwd.random = 2,   
             zero = 0,  
             col.zero = "black", 
             lty.zero = 3)
    dev.off()
  }
  
  return(list(summary = summary,
              genes = consistent_genes))
}

####
####

#### Function to store meta-analysis results ####
generate_list_for_meta_analysis <- function(DNA = FALSE, RNA = FALSE, dna_mtx = list(), rna_mtx = list(), pData = list(), study = list()) {
  
  #### Creating empty list to store meta-analysis results ####
  list_of_studies <- list()
  
  #### Preparing data for meta-analysis ####
  if(DNA && !RNA) {
    common_genes <- find_common_genes(DNA = TRUE, DNA_mtx = dna_mtx)
    for(i in 1:length(study)) {
      
      name <- names(study)[i]
      mtx <- dna_mtx[[i]]
      p <- pData[[i]]
      
      common_samples <- intersect(colnames(mtx), rownames(p))
      mtx_fil <- mtx[, common_samples, drop = FALSE]
      p_fil <- p[common_samples, , drop = FALSE]
      
      #### Dynamically creating list ####
      list_of_studies[[name]] <- list(
        expr = mtx_fil[common_genes, ],
        pheno = p_fil$condition,
        keys = common_genes,
        class = as.numeric(ifelse(p_fil$condition == levels(p_fil$condition)[1], 0, 1))
      )
    }
  }
  
  #### Preparing data for meta-analysis ####
  if(RNA && !DNA) {
    common_genes <- find_common_genes(RNA = TRUE, RNA_mtx = RNA_matrices$norm_counts)
    for(i in 1:length(study)) {
      name <- names(study)[i]
      mtx <- rna_mtx[[i]]
      p <- pData[[i]]
      
      common_samples <- intersect(colnames(mtx), rownames(p))
      mtx_fil <- mtx[, common_samples, drop = FALSE]
      p_fil <- p[common_samples, , drop = FALSE]
     
      #### Dynamically creating list ####
      list_of_studies[[name]] <- list(
        expr = mtx_fil[common_genes, ],
        pheno = p_fil$condition,
        keys = common_genes,
        class = as.numeric(ifelse(p_fil$condition == levels(p_fil$condition)[1], 0, 1))
      )
    }
  }
  
  
  if(DNA && RNA) {
    common_genes <- find_common_genes(DNA = TRUE, DNA_mtx = DNA_matrices, RNA = TRUE, RNA_mtx = RNA_matrices$norm_counts)
    for(i in 1:length(study)) {
      
      name <- names(study)[i]
      list_of_mtx <- c(DNA_mtx, RNA_mtx)
      mtx <- list_of_mtx[[i]]
      p <- pData[[i]]
      
      common_samples <- intersect(colnames(mtx), rownames(p))
      mtx_fil <- mtx[, common_samples, drop = FALSE]
      p_fil <- p[common_samples, , drop = FALSE]
      
      #### Dynamically creating list ####
      list_of_studies[[name]] <- list(
        expr = mtx_fil[common_genes, ],
        pheno = p_fil$condition,
        keys = common_genes,
        class = as.numeric(ifelse(p_fil$condition == levels(p_fil$condition)[1], 0, 1))
      )
    }
  }
  
  return(list(studies = list_of_studies,
              results = meta_results(list_of_studies)))
}

####
####

#### Return pvalues and qvalues ####
get_p_q_values <- function(meta_results) {
  
  cleaned_studies <- list(
    studies = lapply(meta_results$studies, function(study) {
      
      valid_rows <- complete.cases(study$class, study$pheno)
      expr_clean <- study$expr[, valid_rows]
      
      valid_genes <- apply(expr_clean, 1, function(gene_values) {
        all(!is.na(gene_values) & !is.nan(gene_values) & !is.infinite(gene_values))
      })
      expr_clean <- expr_clean[valid_genes, ]
      
      study_cleaned <- list(
        expr = expr_clean,
        pheno = study$pheno[valid_rows],
        keys = study$keys[valid_genes],
        class = study$class[valid_rows]
      )
      return(study_cleaned)
    }),
    
    results = meta_results$results
  )
  
  list.of.sigs <- lapply(cleaned_studies$studies, function( study ) {
    fisher_pval <- ttest.Pvalues(study)
    genes <- cleaned_studies$results$genes
    fisher_pval <- fisher_pval %>% dplyr::filter(keys%in%genes)
  })
  
  output.Fisher <- as.data.frame(sum.of.logs(list.of.sigs)) # For fisher p-value
  adj_pval <- adjust.fisher(output.Fisher=output.Fisher)[, c("F.Qval.up", "F.Qval.down")]
  
  genes <- cleaned_studies$results$genes
  pool <- cleaned_studies$results$summary$pooled.estimates[genes, "summary"]
  se.pool <- cleaned_studies$results$summary$pooled.estimates[genes, "se.summary" ]
  
  es_pval <- numeric()
  es_qval <- numeric()
  
  for(i in 1:nrow(output.Fisher)) {
    if(output.Fisher$F.pval.up[i] <= 0.05) {
      es_pval <- c(es_pval, output.Fisher$F.pval.up[i])
    }
    if(output.Fisher$F.pval.down[i] <= 0.05) {
      es_pval <- c(es_pval, output.Fisher$F.pval.down[i])
    }
  }
  
  for(i in 1:nrow(adj_pval)) {
    if(adj_pval$F.Qval.up[i] <= 0.05) {
      es_qval <- c(es_qval, adj_pval$F.Qval.up[i])
    }
    if(adj_pval$F.Qval.down[i] <= 0.05) {
      es_qval <- c(es_qval, adj_pval$F.Qval.down[i])
    }
  }
  
  sig <- data.frame(Genes = genes,
                    Pooled_es = pool,
                    Pooled_es_se = se.pool,
                    es_pval = es_pval,
                    es_qval = es_qval)
  return(sig)
}
