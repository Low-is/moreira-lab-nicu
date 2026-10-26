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
#   file = "C:/Users/RANDOLPHL/Documents/ARCHS4/human_gene_v2.5.h5"
#   a4$meta$series(file, geo)
# })


#### Creating ENSEMBL annotation database ####
DB_PATH <- "convertIDs.db"
db <- dbConnect(SQLite(), DB_PATH)
dbExecute(db, "CREATE TABLE IF NOT EXISTS gene_annotations (
          ensembl_id TEXT,
          gene_symbol TEXT,
          entrez_id TEXT,
          species TEXT,
        PRIMARY KEY (ensembl_id, species)
);")

fetch_ensembl <- function(dataset) {
  mart <- useMart("ensembl", dataset=dataset)
  getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id"),
    mart=mart,
    useCache = FALSE
  )
}

#### Download *every* ENSEMBL dataset ####
ensembl_root <- useEnsemblGenomes(biomart = "ensembl")
datasets <- listDatasets(ensembl_root)


generate_exprs_mtx <- function(DNA = NULL, RNA = NULL, dna_files = list(), rna_mtx = list(), pData_RNA = list()) {
  #### Changes 5-29-2025: Changing 'dna_studies' to 'dna_files' for naming consistency ####
  #### Additions 5-29-2025: Going to add argument for selecting species type prior to gene annotations ####
  
  #### Handling duplicated genes ####
  handling_duplicates <- function(data_table) {
    
    data_table <- na.omit(data_table[Genes != "" & !is.na(Genes)])
    
    
    if (any(duplicated(data_table$Genes))) {
      exprs_data.table <- data_table[, lapply(.SD, mean), by = Genes, .SDcols =! 'Genes']
    }
    
    exprs_mtx <- as.matrix(exprs_data.table[, -1])
    rownames(exprs_mtx) <- exprs_data.table$Genes
    
    return(exprs_mtx)
  }
  
  #### For DNA microarray data sets ####
  if(!is.null(DNA) && DNA == TRUE) {
    
    #### Checking if dna_files is a list ####
    if(!is.list(dna_files)) stop("dna_files needs to be a list")
    
    #### Checking if all elements in the list are character strings ####
    if(!all(sapply(dna_files, is.character))) stop("Elements in dna_files needs to be a list of character strings\nExample: list(dataset1 = 'GSEXXXX', dataset2 = 'GSEXXXX') ")
    
    ### Extracting list of DNA microarray expression matrices ####
    list_of_geo_data <- lapply(dna_files, function( file ) {
      tryCatch({
        geo_data <- getGEO(GEO = file, AnnotGPL = TRUE)[[1]]
        if(is.null(geo_data) || length(geo_data) == 0) stop(paste("Failed to retrieve GEO dataset:", file))
        return(geo_data)
      }, error = function(e) {
        warning(paste("Error retrieving", file, ":", e$message))
        return(NULL) ### Return NULL in case of failure
      })
    })
    
    #### Extracting expression matrices ####
    list_of_dna_mtx <- lapply(list_of_geo_data, function( geo ) {
      if(!is.null(geo)) {
        return(exprs(geo))
      } else {
        return(NULL)
      }
    })
    
    #### Extracting feature data ####
    list_of_fData <- lapply(list_of_geo_data, function( geo ){
      if(!is.null(geo)) {
        return(fData(geo))
      } else {
        return(NULL)
      }
    })
    
    #### Annotating probe IDs ####
    for(i in 1:length(list_of_fData)) {
      fData <- list_of_fData[[i]]
      gene_cols <- c("Gene symbol", "ILMN_Gene", "GB_ACC")
      matched_col <- intersect(gene_cols, colnames(fData))
      
      if(length(matched_col) > 0) {
        selected_col <- matched_col[1]  ### Use the first matched column
        
        #### Ensure probe IDs are character strings ####
        probe_ids <- as.character(fData[[selected_col]])
        
        #### Check if conversion is needed ####
        if(!(selected_col %in% c("Gene symbol", "ILMN_Gene"))) {
          keytype <- ifelse(selected_col == "GB_ACC", "ACCNUM")
          
          #### Convert probe IDs to gene symbols ####
          fData$Genes <- mapIds(
            x = org.Hs.eg.db,
            keys = probe_ids,
            column = "SYMBOL",
            keytype = keytype,
            multiVals = "first"
          )
        } else {
          #### If already gene symbols, just assign them ####
          fData$Genes <- fData[[selected_col]]
        }
      }
      
      list_of_fData[[i]] <- fData  ### Save updated feature data back to list
    }
    
    #### Converting expression matrices to data.table ####
    list_of_dna_data_table <- list()
    
    for (i in 1:length(list_of_dna_mtx)) {
      mtx <- list_of_dna_mtx[[i]]
      fData <- list_of_fData[[i]]
      
      genes_vector <- as.character(fData$Genes)
      
      # Create data.table for each entry
      df <- data.frame(Genes = genes_vector, mtx)
      
      data_table <- as.data.table(df)
      
      # Append the data_table to the list
      list_of_dna_data_table[[i]] <- data_table
    }
    
    list_of_dna_mtx <- lapply(list_of_dna_data_table, function ( tbl ) {
      mtx <- handling_duplicates( tbl )
      return(mtx)
    })
    
    ##### Quantile normalization #####
    list_of_dna_mtx <- lapply(list_of_dna_mtx, function( mtx ) {
      gene_names <- rownames( mtx )
      
      exprs_mtx <- apply(mtx, 2, as.numeric)
      
      log_mtx <- log(exprs_mtx + 1)
      
      if(any(is.nan(log_mtx)) == TRUE) {
        log_mtx[is.nan(log_mtx)] <- 0
      }
      norm_mtx <- normalizeBetweenArrays(log_mtx, method = "quantile")
      
      #### Re-assigning gene names to the normalized matrix ####
      rownames(norm_mtx) <- gene_names
      return(norm_mtx)
    })
    
    #### Saving boxplots of normalized data ####
    for(i in 1:length(list_of_dna_mtx)) {
      pdf(paste0("boxplot_", names(dna_files)[i], ".pdf"), width = 8, height = 8)
      boxplot <- boxplot(list_of_dna_mtx[[i]], outline = FALSE,
                         col = "skyblue", main = paste0("Log transformation/Quantile", "\n", names(dna_files)[i]))
      dev.off()
    }
    
    return( list_of_dna_mtx )
  }
  
  #### For RNA-Seq data sets ####
  #### Changes 5-29-2025: Changing 'rna_mtx' to 'rna_files' for naming consistency ####
  if(!is.null(RNA) && RNA == TRUE) {
    #### Checking if rna_mtx is a list ####
    if(!is.list(rna_files)) stop("rna_files needs to be a list")
    
    #### Checking if pData_RNA has "condition" column ####
    if(!is.list(pData_RNA)) stop("pData_RNA needs to be a list")

    #### Checking there is a column named "condition"
    if(!all(sapply(pData_RNA, function(x) "condition" %in% colnames(x)))) {
      stop("Each element in pData_RNA must have a 'condition' column")
    }
    
    #### Creating empty list ####
    list_of_norm_RNA_mtx <- list()
    list_of_dds <- list()
    
    for(i in 1:length(list_of_RNA_mtx)) {
      mtx <- list_of_RNA_mtx[[i]]
      pData <- pData_RNA[[i]]
      
      dds_1 <- DESeqDataSetFromMatrix(round(mtx), colData = pData, design = ~ condition)
      dds_2 <- DESeq( dds_1 )
      
      norm_counts <- counts(dds_2, normalized = TRUE)
      
      list_of_dds[[i]] <- dds_2
      list_of_norm_RNA_mtx[[i]] <- norm_counts
    }
    
    #### Saving boxplots of normalized data ####
    for(i in 1:length(list_of_norm_RNA_mtx)) {
      pdf(paste0("boxplot_", names(rna_files)[i], ".pdf"), width = 8, height = 8)
      boxplot <- boxplot(list_of_norm_RNA_mtx[[i]], outline = FALSE,
                         col = "skyblue", main = paste0("Median of Ratios", "\n", names(rna_files)[i]))
      dev.off()
    }
    
    return( list(dds_objs = list_of_dds,
                 norm_counts = list_of_norm_RNA_mtx) )
  }
}

#### Testing function ####
#### 11-19-2025 ####
dna_files <- list(GSE3140 = "GSE3140",
                  GSE4607 = "GSE4607",
                  GSE8121 = "GSE8121",
                  GSE9692 = "GSE9692",
                  GSE10361 = "GSE10361",
                  GSE11755 = "GSE11755",
                  GSE13904 = "GSE13904",
                  GSE15495 = "GSE15495",
                  GSE26378 = "GSE26378",
                  GSE30119 = "GSE30119",
                  GSE25504 = "GSE25504",
                  GSE53471 = "GSE53471",
                  GSE61477 = "GSE61477",
                  GSE66099 = "GSE66099",
                  GSE69686 = "GSE69686",
                  GSE70461 = "GSE70461",
                  GSE101880 = "GSE101880",
                  GSE119217 = "GSE119217",
                  GSE145227 = "GSE145227",
                  GSE138712 = "GSE138712",
                  GSE100159 = "GSE100159",
                  GSE150707 = "GSE150707",
                  GSE236877 = "GSE236877",
                  GSE279448 = "GSE279448")
DNA_matrices <- generate_exprs_mtx(DNA = TRUE,  RNA = NULL, dna_studies = dna_files)
# RNA_matrics <- generate_exprs_mtx(DNA = NULL, RNA = TRUE, rna_mtx = rna_mtx, pData_RNA = pData_RNA)


list_of_pData <- list()

for(i in seq_along(dna_files)) {
  
  gse <- getGEO(dna_files[[i]], AnnotGPL = TRUE)[[1]]
  
  pData <- pData(gse)
  
  list_of_pData[[i]] <- pData
}

all_pData <- do.call(rbind, list_of_pData)