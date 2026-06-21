#' @title DiseaseCIO Client
#' @description An R6 client to interact with the DiseaseCIO database. It allows users to query metadata, filter datasets, batch-download files, and load data directly into R memory for downstream visualization and analysis.
#' @importFrom R6 R6Class
#' @importFrom data.table fread
#' @importFrom fst read.fst
#' @importFrom base64enc base64decode
#' @importFrom httr GET http_error http_status progress
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual theme_bw labs geom_vline geom_hline geom_boxplot geom_jitter theme_classic
#' @importFrom ggpubr stat_compare_means
#' @importFrom pROC roc ggroc
#' @importFrom rlang .data
#' @export
DiseaseCIO <- R6Class("DiseaseCIO",
                     private = list(
                       api_token = "aHR0cDovLzIxOC44LjI0MS4yNDg6MzgzOC9EaXNUeFJFU1AvZG93bmxvYWQv",
                       feature_prefix_map = c(
                         "Gene Expression" = "Gene",
                         "Pathway Activity" = "Pathway",
                         "Immune Infiltration(RNA)" = "Immune_Infiltration_RNA",
                         "TF Regulon Activity" = "Transcription_Factor",
                         "Secreted Protein Activity"= "Secreted_Protein_Activity",
                         "Cell Death Mode" = "Cell_Death",
                         "Metabolic Flux Score" = "Metabolic_Flux",
                         "Ligand-Receptor Pair" = "LR_Pair",

                         "Taxonomy Abundance"= "Taxonomy_Abundance",
                         "Alpha Diversity" = "Alpha_Diversity",
                         "Microbial Function" = "Micro_Pathway",
                         "Enzymatic Activity" = "Enzymatic_Activity",
                         "Antibiotic Resistance Gene" = "Antibiotic_Resistance_Gene",
                         "Virulence Factor" = "Virulence_Factor",

                         "CpG-level Methylation" = "CpGmethy",
                         "Gene-level Methylation" = "Genemethy",
                         "Epigenomic Regulation" = "Epigenomic_Regulation",
                         "Immune Infiltration(DNAm)" = "Immune_Infiltration_DNAm",
                         "Epigenetic Clock" = "Epigenetic_clock",
                         "Episcore" = "Episcore"
                       ),
                       make_request = function(endpoint) {
                         dl_url = paste0(rawToChar(base64decode(private$api_token)), endpoint)
                         return(dl_url)
                       },
                       download_file = function(url, destfile, progress=TRUE) {
                         if (progress) {
                           response <- httr::GET(url,
                                                 httr::write_disk(destfile, overwrite = TRUE),
                                                 httr::progress(),
                                                 httr::config(max_recv_speed_large = 1000000))
                         }else{
                           response <- httr::GET(url,
                                                 httr::write_disk(destfile, overwrite = TRUE),
                                                 httr::config(max_recv_speed_large = 1000000))
                         }
                         if(httr::http_error(response)) {
                           stop("Download failed: ", httr::http_status(response)$message)
                         }
                         return(invisible(response))
                       },
                       query_file = function(url,progress=FALSE,file_type="csv") {
                         temp_file <- tempfile()
                         sig = tryCatch({
                           private$download_file(url = url,destfile = temp_file,progress=progress)
                         }, error = function(e) {
                           return(NULL)
                         })
                         if (is.null(sig)) {
                           stop("Connection error!")
                         }
                         if (file_type=="csv") {
                           temp <- fread(temp_file,data.table = FALSE,header = TRUE)
                         }else if (file_type=="fst") {
                           temp <- read.fst(temp_file)
                         }
                         unlink(temp_file)
                         return(temp)
                       }
                     ),
                     public = list(
                       #' @field metadata A data.frame containing the full metadata manifest downloaded from the server.
                       metadata = NULL,
                       #' @field sub_table A data.frame containing the filtered subset of the metadata.
                       sub_table = NULL,
                       #' @field local_data A nested list storing the loaded clinical data, feature matrices, and differential results in R memory.
                       local_data = NULL,

                       #' @description Create a new database client and the fetch full metadata.
                       #' @param manifest String. The local file path for the metadata manifest. If the file already exists locally, it will be read directly; otherwise, it will be downloaded from the server. Defaults to "Download_Manifest.csv".
                       #' @return A new \code{DiseaseCIO} object.
                       initialize = function(manifest="Download_Manifest.csv") {
                         if (!file.exists(manifest)) {
                           message("Connecting to server and fetching metadata table...")
                           meta_url <- private$make_request("Download_Manifest.csv")

                           private$download_file(url = meta_url, destfile = manifest)

                         }else{
                           message("Local file found, reading...")
                         }
                         self$metadata <- fread(manifest,data.table = FALSE,header = TRUE,showProgress=FALSE)
                         self$sub_table <- self$metadata
                       },

                       #' @description Filter the metadata table based on specific biological or clinical requirements.
                       #' @param dataset Character vector. The dataset identifiers from the DiseaseCIO database (e.g., c("D100001")).
                       #' @param omic Character vector. The omics data types to retain (e.g., c("Transcriptomics")"").
                       #' @param disease Character vector. The primary disease types to filter by (e.g., c("Leukemia")).
                       #' @param disease_sub Character vector. The specific disease subtypes (e.g., c("Acute myeloid leukemia")).
                       #' @param therapy Character vector. The therapy categorys applied to the cohorts (e.g., c("Targeted Therapy")).
                       #' @param treatment Character vector. The treatment regimens (e.g., c("Ruxolitinib")).
                       #' @param sampling_location Character vector. The tissue origins or sampling locations(e.g., c("Tissue","Bone Marrow")).
                       #' @param feature Character vector. The feature types to retain (e.g., c("Gene Expression","Clinical Data", "Corpus", "QA")).
                       #' @param file_type Character vector. The file types to retain (e.g., c("Meta Info", "Feature Matrix", "AI-ready Data")).
                       #' @param min_size Integer. The minimum sample size required for a dataset to be retained. Defaults to 0.
                       #' @return Returns the modified R6 object invisibly, allowing for method chaining.
                       filter_metadata = function(dataset=NULL,omic = NULL,disease = NULL,
                                                  disease_sub = NULL,therapy = NULL,
                                                  treatment = NULL,sampling_location = NULL,
                                                  feature = NULL, file_type = NULL,min_size=0) {
                         if (is.null(self$metadata)) stop("Metadata is empty. Please re-initialize.")
                         temp_df <- self$metadata

                         if (!is.null(dataset)) {
                           temp_df <- temp_df[temp_df$Dataset %in% dataset, ]
                         }
                         if (!is.null(omic)) {
                           temp_df <- temp_df[grepl(paste0(omic,collapse = "|"),temp_df$Omics), ]
                         }
                         if (!is.null(disease)) {
                           temp_df <- temp_df[grepl(paste0(disease,collapse = "|"),temp_df$Disease_Type), ]
                         }
                         if (!is.null(disease_sub)) {
                           temp_df <- temp_df[grepl(paste0(disease_sub,collapse = "|"),temp_df$Disease_Subtype), ]
                         }
                         if (!is.null(therapy)) {
                           temp_df <- temp_df[grepl(paste0(therapy,collapse = "|"),temp_df$Therapy_Category), ]
                         }
                         if (!is.null(treatment)) {
                           temp_df <- temp_df[grepl(paste0(treatment,collapse = "|"),temp_df$Treatment_Regimen), ]
                         }
                         if (!is.null(sampling_location)) {
                           temp_df <- temp_df[grepl(paste0(sampling_location,collapse = "|"),temp_df$Sampling_Location), ]
                         }
                         if (!is.null(feature)) {
                           temp_df <- temp_df[grepl(paste0(feature,collapse = "|"),temp_df$Feature_Type), ]
                         }
                         if (!is.null(file_type)) {
                           temp_df <- temp_df[grepl(paste0(file_type,collapse = "|"),temp_df$File_Type), ]
                         }
                         temp_df = temp_df[temp_df$Sample_Size>=min_size,]

                         row.names(temp_df)=NULL
                         self$sub_table <- temp_df
                         message(sprintf("Filtered down to %d records.", nrow(self$sub_table)))

                         return(invisible(self))
                       },
                       #' @description Search module function of DiseaseCIO.
                       #' @param feature_type String. The type of feature (e.g., "Gene Expression").You can view all feature types by using `unique(client$metadata$Feature_Type)`.
                       #' @param feature_name String. The name of feature (e.g., "CTLA4").
                       #' @param disease Character vector. The primary disease types to filter by (e.g., c("Leukemia")).
                       #' @param therapy Character vector. The therapy categorys applied to the cohorts (e.g., c("Targeted Therapy")).
                       #' @param treatment Character vector. The treatment regimens (e.g., c("Ruxolitinib")).
                       #' @return A \code{data.frame} object.
                       search_DCIO = function(feature_type,feature_name,disease=NULL,therapy=NULL,treatment=NULL,
                                              threshold = c("FDR<0.05", "FDR<0.2")) {
                         threshold <- match.arg(threshold)
                         if (length(feature_type)!=1) stop("Only one feature_type id can be input!")
                         if (length(feature_name)!=1) stop("Only one feature_name id can be input!")

                         message("Connecting to DiseaseCIO...")
                         feature = private$feature_prefix_map[feature_type]
                         if (length(feature)!=1) stop("Please enter the correct feature_type!")
                         search_all <- switch(threshold,
                                             "FDR<0.05" = paste0("searchall/",feature,"_SearchAll_out_fdr0.05.fst"),
                                             "FDR<0.2"  = paste0("searchall/",feature,"_SearchAll_out_fdr0.2.fst"))
                         metadata <- private$query_file(private$make_request(search_all),file_type = "fst")

                         message("Searching data...")
                         metadata = metadata[metadata$ID == feature_name,]
                         if (!is.null(disease)) {
                           metadata <- metadata[grepl(paste0(disease,collapse = "|"),metadata$Disease_Type), ]
                         }
                         if (!is.null(therapy)) {
                           metadata <- metadata[grepl(paste0(therapy,collapse = "|"),metadata$Therapy_Category), ]
                         }
                         if (!is.null(treatment)) {
                           metadata <- metadata[grepl(paste0(treatment,collapse = "|"),metadata$Treatment_Regimen), ]
                         }
                         if (nrow(metadata)==0) stop("No dataset found!")
                         rownames(metadata) = NULL
                         return(metadata)
                       },
                       #' @description Browse module function of DiseaseCIO.
                       #' @param dataset String. The dataset identifier (e.g., "D100001").
                       #' @param feature_type String. The type of feature (e.g., "Gene Expression").You can view all feature types by using `unique(client$metadata$Feature_Type)`.
                       #' @param feature_name String. The name of feature (e.g., "CTLA4").
                       #' @return A \code{list} object.
                       browse_DCIO = function(dataset,feature_type,feature_name) {
                         if (length(dataset)!=1) stop("Only one dataset id can be input!")
                         if (length(feature_type)!=1) stop("Only one feature_type id can be input!")
                         if (length(feature_name)!=1) stop("Only one feature_name id can be input!")

                         message("Connecting to DiseaseCIO...")
                         metadata <- private$query_file(private$make_request("Download_Manifest.csv"))

                         message("Querying data...")
                         metadata = metadata[metadata$File_Type=="Differential Results",]

                         metadata = metadata[metadata$Dataset == dataset,]
                         metadata = metadata[metadata$Feature_Type == feature_type,]
                         file_diff = metadata$File_URL
                         if (length(file_diff)!=1) stop("No file found!")

                         diff_tb = private$query_file(gsub("csv$","fst",file_diff),file_type = "fst")
                         diff_tb = diff_tb[diff_tb$ID==feature_name,]

                         if (nrow(diff_tb)==0) stop("No feature found!")
                         if (nrow(diff_tb)!=1) stop("Matching error, please check your input!")
                         res = as.list(diff_tb)
                         return(res)
                       },

                       #' @description Batch download files in the filtered \code{sub_table} to a local directory.
                       #' @param output_dir String. The target local directory path to save the downloaded files. Defaults to "downloaded_data".
                       #' @param skip_existing Logical. Whether to skip downloading files that already exist in the target directory. Defaults to TRUE.
                       #' @return Returns the modified R6 object invisibly.
                       download_files = function(output_dir = "downloaded_data",skip_existing=TRUE) {
                         if (nrow(self$sub_table) == 0) stop("The filtered table is empty. Nothing to download.")

                         for (i in seq_len(nrow(self$sub_table))) {
                           row_data <- self$sub_table[i, ]
                           data_id <- row_data$Dataset
                           f_type <- gsub(" ","_",row_data$File_Type)
                           f_name <- row_data$File_Name

                           temp_dir <- ifelse(f_type=="Meta_Info",paste0(output_dir,"/",data_id),paste0(output_dir,"/",data_id,"/",f_type))
                           file_name <- paste0(temp_dir,"/",f_name)
                           if (skip_existing) {
                             if (file.exists(file_name)) {
                               message(sprintf("File already exists: %s...", f_name))
                               next
                             }
                           }

                           if (!dir.exists(temp_dir)) dir.create(temp_dir, recursive = TRUE)

                           dl_url <- private$make_request(sprintf("%s/%s",data_id, f_name))

                           message(sprintf("Downloading %s...", f_name))
                           private$download_file(url = dl_url,destfile = file_name)
                         }
                         message(" \n Batch download complete.")
                         return(invisible(self))
                       },
                       #' @description Fetch the filtered files and load them into an R list structure.
                       #' @param local_dir String. The directory path containing locally downloaded files.
                       #' @return Returns the modified R6 object invisibly. The loaded data is accessible via the \code{local_data} field.
                       load_to_memory = function(local_dir="downloaded_data") {
                         if (nrow(self$sub_table) == 0) stop("The filtered table is empty. Nothing to load.")

                         result_list <- list()

                         for (i in seq_len(nrow(self$sub_table))) {
                           row_data <- self$sub_table[i, ]
                           data_id <- row_data$Dataset
                           f_type <- gsub(" ","_",row_data$File_Type)
                           f_name <- row_data$File_Name
                           feature <- row_data$Feature_Type

                           if (!is.null(local_dir)) {
                             temp_dir <- ifelse(f_type=="Meta_Info",paste0(local_dir,"/",data_id),paste0(local_dir,"/",data_id,"/",f_type))
                             file_name <- paste0(temp_dir,"/",f_name)
                             if (file.exists(file_name)) {
                               message(sprintf("File found, loading %s into memory...", f_name))
                             }else{
                               stop(sprintf("The file was not found: %s", f_name))
                             }
                           }else{
                             message("Loading path cannot be empty!")
                           }

                           tryCatch({
                             if (f_type=="Meta_Info") {
                               result_list[[data_id]][[f_type]] <- fread(file_name, data.table = FALSE, header = TRUE,showProgress=FALSE)
                             }else{
                               result_list[[data_id]][[f_type]][[feature]] <- fread(file_name, data.table = FALSE, header = TRUE,showProgress=FALSE)
                             }
                           }, error = function(e) {
                             warning(sprintf("Failed to load %s: %s", data_id, e$message))
                           })
                         }
                         message(" \n All selected files successfully loaded into memory.")
                         self$local_data = result_list
                         return(invisible(self))
                       },
                       #' @description Generate a Volcano Plot for differential analysis results.
                       #' @param dataset String. The dataset identifier (e.g., "D100001").
                       #' @param feature_type String. The type of feature (e.g., "Gene Expression").
                       #' @param logX_col String. The column name representing the log2 Fold Change (or LogOR) in the differential matrix. Defaults to "LogOR".
                       #' @param pval_col String. The column name representing the P-value or FDR in the differential matrix. Defaults to "P.value".
                       #' @param fc_cutoff Numeric. The absolute threshold for fold change significance. Defaults to 1.0.
                       #' @param p_cutoff Numeric. The threshold for P-value significance. Defaults to 0.05.
                       #' @return A \code{ggplot} object.
                       plot_volcano = function(dataset, feature_type, logX_col = "LogOR", pval_col = "P.value",
                                               fc_cutoff = 1.0, p_cutoff = 0.05) {

                         if (is.null(self$local_data[[dataset]])) stop("The dataset is not found!")
                         diff_df <- self$local_data[[dataset]]$Differential_Results[[feature_type]]
                         if (is.null(diff_df)) stop("No difference matrix found for this feature type!")

                         diff_df$Significance <- "Not Significant"
                         diff_df$Significance[diff_df[[logX_col]] > fc_cutoff & diff_df[[pval_col]] < p_cutoff] <- "Up"
                         diff_df$Significance[diff_df[[logX_col]] < -fc_cutoff & diff_df[[pval_col]] < p_cutoff] <- "Down"

                         p <- ggplot2::ggplot(diff_df, ggplot2::aes(x = .data[[logX_col]],
                                                                    y = -log10(.data[[pval_col]]),
                                                                    color = .data[["Significance"]])) +
                           ggplot2::geom_point(alpha = 0.8, size = 1.5) +
                           ggplot2::scale_color_manual(values = c("Up" = "#d73027", "Down" = "#4575b4", "Not Significant" = "grey80")) +
                           ggplot2::geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "black") +
                           ggplot2::geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "black") +
                           ggplot2::theme_bw() +
                           ggplot2::labs(title = sprintf("%s (%s)", dataset, feature_type),
                                         x = logX_col, y = "-log10(P-value)")

                         return(p)
                       },

                       #' @description Generate a Boxplot to visualize feature expression across different clinical groups.
                       #' @param dataset String. The dataset identifier (e.g., "D100001").
                       #' @param feature_type String. The type of feature matrix (e.g., "Gene Expression").
                       #' @param feature_name String. The specific feature name to plot (e.g., "CD274").
                       #' @param group_col String. The column name in the clinical metadata used for grouping samples. Defaults to "Response".
                       #' @param stat_method String. The Method for comparing means.
                       #' @return A \code{ggplot} object.
                       plot_boxplot = function(dataset, feature_type, feature_name, group_col="Response",stat_method="wilcox.test") {

                         clin_df <- self$local_data[[dataset]]$Meta_Info
                         if (is.null(clin_df)) stop("No clinical information found!")
                         feat_mat <- self$local_data[[dataset]]$Feature_Matrix[[feature_type]]
                         if (is.null(feat_mat)) stop("No matrix found for this feature type!")

                         if (!feature_name %in% feat_mat[,1]) stop(sprintf("Feature not found in matrix: %s", feature_name))
                         if (!group_col %in% colnames(clin_df)) stop(sprintf("Grouping column not found in clinical information: %s", group_col))

                         common_samples <- intersect(clin_df$Sample, colnames(feat_mat))
                         if (length(common_samples) == 0) stop("There is no matching sample ID between clinical information and feature matrix!")

                         plot_data <- data.frame(
                           Sample = common_samples,
                           Value = as.numeric(feat_mat[feat_mat[,1]==feature_name, common_samples]),
                           Group = as.character(clin_df[match(common_samples,clin_df$Sample), group_col])
                         )

                         plot_data <- plot_data[!is.na(plot_data$Group), ]

                         p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Group, y = Value, fill = Group)) +
                           ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
                           ggplot2::geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
                           ggplot2::theme_classic() +
                           ggplot2::labs(title = sprintf("Feature: %s", feature_name),
                                         subtitle = sprintf("Dataset: %s", dataset),
                                         y = "Feature Value", x = group_col) +
                           ggplot2::theme(legend.position = "none")+
                           ggpubr::stat_compare_means(
                             method = stat_method,
                             label = "p.format",
                             label.x = 1.5
                           )

                         return(p)
                       },

                       #' @description Generate an ROC Curve to evaluate the predictive performance of a specific feature.
                       #' @param dataset String. The dataset identifier (e.g., "D100001").
                       #' @param feature_type String. The type of feature matrix (e.g., "Gene Expression").
                       #' @param feature_name String. The specific feature name to evaluate (e.g., "CD274").
                       #' @param group_col String. The column name in the clinical metadata representing the true binary response. Defaults to "Response".
                       #' @param positive_class String. The label in the grouping column that represents the positive class (e.g., "R" for Responder). Defaults to "R".
                       #' @return A \code{ggplot} object.
                       plot_roc = function(dataset, feature_type, feature_name, group_col="Response", positive_class="R") {
                         if (is.null(self$local_data[[dataset]])) stop("The dataset is not found!")
                         clin_df <- self$local_data[[dataset]]$Meta_Info
                         if (is.null(clin_df)) stop("No clinical information found!")
                         feat_mat <- self$local_data[[dataset]]$Feature_Matrix[[feature_type]]
                         if (is.null(feat_mat)) stop("No matrix found for this feature type!")

                         if (!feature_name %in% feat_mat[,1]) stop(sprintf("Feature not found in matrix: %s", feature_name))
                         if (!group_col %in% colnames(clin_df)) stop(sprintf("Grouping column not found in clinical information: %s", group_col))

                         common_samples <- intersect(clin_df$Sample, colnames(feat_mat))
                         if (length(common_samples) == 0) stop("There is no matching sample ID between clinical information and feature matrix!")

                         plot_data <- data.frame(
                           Sample = common_samples,
                           Value = as.numeric(feat_mat[feat_mat[,1]==feature_name, common_samples]),
                           Group = as.character(clin_df[match(common_samples,clin_df$Sample), group_col])
                         )
                         plot_data <- plot_data[!is.na(plot_data$Group), ]

                         if (!positive_class %in% plot_data$Group) stop("The specified 'positive_class' cannot be found!")
                         plot_data$Label <- ifelse(plot_data$Group == positive_class, 1, 0)

                         roc_obj <- pROC::roc(response = plot_data$Label, predictor = plot_data$Value, quiet = TRUE)
                         auc_val <- round(roc_obj$auc, 3)

                         p <- pROC::ggroc(roc_obj, color = "#d73027", size = 1) +
                           ggplot2::theme_minimal() +
                           ggplot2::geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "grey50") +
                           ggplot2::labs(title = sprintf("Feature: %s", feature_name),
                                         subtitle = sprintf("AUC = %s", auc_val),
                                         x = "Specificity", y = "Sensitivity")

                         return(p)
                       },

                       #' @description Generate a Pie Chart for categorical clinical information.
                       #' @param dataset String. The dataset identifier (e.g., "D100001").
                       #' @param clinical_col String. The column name in the clinical metadata representing the categorical variable (e.g., "Response", "Sex").
                       #' @return A \code{ggplot} object.
                       plot_pie = function(dataset, clinical_col) {
                         if (is.null(self$local_data[[dataset]])) stop("The dataset is not found!")
                         clin_df <- self$local_data[[dataset]]$Meta_Info
                         if (is.null(clin_df)) stop("No clinical information found!")
                         if (!clinical_col %in% colnames(clin_df)) stop(sprintf("Column not found in clinical information: %s", clinical_col))

                         clin_vector <- na.omit(clin_df[[clinical_col]])
                         if (length(clin_vector) == 0) stop("The specified clinical column only contains NA values!")

                         plot_data <- as.data.frame(table(clin_vector))
                         colnames(plot_data) <- c("Category", "Count")
                         plot_data$Percentage <- prop.table(plot_data$Count) * 100

                         plot_data$Label <- sprintf("%s\n(%.1f%%)", plot_data$Category, plot_data$Percentage)

                         p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = "", y = Count, fill = Category)) +
                           ggplot2::geom_bar(stat = "identity", width = 1, color = "white", size = 0.5) +
                           ggplot2::coord_polar(theta = "y", start = 0) +
                           ggplot2::geom_text(ggplot2::aes(label = Label),
                                              position = ggplot2::position_stack(vjust = 0.5),
                                              size = 4, color = "white", fontface = "bold") +
                           ggplot2::theme_void() +
                           ggplot2::labs(title = sprintf("Distribution of %s", clinical_col),
                                         subtitle = sprintf("Dataset: %s (N = %d)", dataset, sum(plot_data$Count)),
                                         fill = clinical_col) +
                           ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                                          plot.subtitle = ggplot2::element_text(hjust = 0.5))
                         return(p)
                       }
                     )
)
