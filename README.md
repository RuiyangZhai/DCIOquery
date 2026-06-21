# DCIOquery

**DCIOquery** is an intuitive, R6-based client package designed to interact with the DiseaseCIO database. It streamlines the entire bioinformatics workflow from data retrieval to visualization. 

With `DiseaseCIO`, you can easily:
* **Query & Filter**: Search the global metadata manifest using specific biological and clinical criteria (e.g., disease type, omics data, treatment regimens).
* **Search & Browse**: Use the DiseaseCIO database search and browsing functions.
* **Batch Download**: Automatically download selected clinical data, feature matrices, and differential analysis results to your local machine.
* **In-Memory Analysis**: Load targeted datasets directly into the R environment as structured lists.
* **Built-in Visualization**: Generate high-quality, publication-ready plots (Volcano plots, Boxplots, and ROC curves) with simple, one-line commands.

## Installation

You can install the development version of DTXRquery from GitHub using the `devtools` or `remotes` package:

```R
# install.packages("devtools")
devtools::install_github("RuiyangZhai/DTXRquery")
```

## Quick Start / Examples
Here is a typical workflow demonstrating how to initialize the client, filter data, load it into memory, and generate plots.

#### 1. Initialize the Client
Creating a new `DisTxRESP` object will automatically fetch the latest data manifest from the server.
```R
library(DTXRquery)

# Initialize the client
client <- DisTxRESP$new()
# Connecting to server and fetching metadata table...
# Downloading: 310 kB

# View the full metadata
head(client$metadata)

```
#### 2. Filter and Load Data
You can chain methods together to filter the required data sets, download them and load them into memory for analysis.
```R
# Filter and screen specific omics, diseases, drugs, etc
client$filter_metadata(omic = c("Transcriptomics"),
                       disease = c("Inflammatory bowel disease"),
                       intervention = c("Golimumab"),
                       feature = c("Clinical Data","Gene Expression"))
# Filtered down to 6 records.

# Check the filtered subset table
head(client$sub_table)
# Dataset           Omics               Disease_Type    Disease_Subtype        Treatment_Regimen Intervention Sampling_Location Sample_Size    Feature_Type                        File_Name            File_Type
# 1 DTXR100447 Transcriptomics Inflammatory bowel disease Ulcerative Colitis Immunomodulatory Therapy    Golimumab            Tissue          84   Clinical Data             DTXR100447_pdata.csv            Meta Info
# 2 DTXR100447 Transcriptomics Inflammatory bowel disease Ulcerative Colitis Immunomodulatory Therapy    Golimumab            Tissue          84 Gene Expression DTXR100447_Differential_Gene.csv Differential Results
# 3 DTXR100447 Transcriptomics Inflammatory bowel disease Ulcerative Colitis Immunomodulatory Therapy    Golimumab            Tissue          84 Gene Expression              DTXR100447_Gene.csv       Feature Matrix
# 4 DTXR100448 Transcriptomics Inflammatory bowel disease Ulcerative Colitis Immunomodulatory Therapy    Golimumab            Tissue          59   Clinical Data             DTXR100448_pdata.csv            Meta Info
# 5 DTXR100448 Transcriptomics Inflammatory bowel disease Ulcerative Colitis Immunomodulatory Therapy    Golimumab            Tissue          59 Gene Expression DTXR100448_Differential_Gene.csv Differential Results
# 6 DTXR100448 Transcriptomics Inflammatory bowel disease Ulcerative Colitis Immunomodulatory Therapy    Golimumab            Tissue          59 Gene Expression              DTXR100448_Gene.csv       Feature Matrix

# Batch download files
your_save_dir = "my_data"
client$download_files(your_save_dir)
# Downloading DTXR100447_pdata.csv...
# Downloading: 390 B     
# Downloading DTXR100447_Differential_Gene.csv...
# Downloading: 500 kB     
# Downloading DTXR100447_Gene.csv...
# Downloading: 13 MB     
# Downloading DTXR100448_pdata.csv...
# Downloading: 440 B     
# Downloading DTXR100448_Differential_Gene.csv...
# Downloading: 410 kB     
# Downloading DTXR100448_Gene.csv...
# Downloading: 9.2 MB     
# 
# Batch download complete.

# Filter specific data sets and feature types
client$filter_metadata(
  dtxr = c("DTXR600031","DTXR500068"),
  feature = c("Clinical Data","Microbial Abundance")
)
# Filtered down to 6 records.

# Batch download files and load into memory
client$download_files(your_save_dir)$load_to_memory(your_save_dir)
# File found, loading DTXR500068_pdata.csv into memory...
# File found, loading DTXR500068_Differential_Microbial_Composition.csv into memory...
# File found, loading DTXR500068_Microbial_Composition.csv into memory...
# File found, loading DTXR600031_pdata.csv into memory...
# File found, loading DTXR600031_Differential_Microbial_Composition.csv into memory...
# File found, loading DTXR600031_Microbial_Composition.csv into memory...
# 
# All selected files successfully loaded into memory.
```

#### 3. Data Visualization
Once the data is loaded into memory, you can use the built-in visualization methods.

__Pie Plot__
Pie chart for visualizing clinical baseline characteristics：
```R
# View the proportion of different response groups in the dataset
pie_plt <- client$plot_pie(dataset = "DTXR600031",clinical_col = "Response")

# View gender ratio of patients
# pie_plt <- client$plot_pie(dataset = "DTXR600031",clinical_col = "Sex")
print(pie_plt)
```
<div align=center>
  <img src="https://github.com/RuiyangZhai/img/blob/main/DTXRquery/pie_plt.png?raw=true" width="600">
</div>

__Volcano Plot__

Visualize differential expression results quickly:
```R
# Plot a volcano plot for a specific dataset
volcano_plt <- client$plot_volcano(
  dataset = "DTXR600031", 
  feature_type = "Microbial Abundance",
  logX_col = "LogOR",
  pval_col = "P.value",
  fc_cutoff = 1.0, 
  p_cutoff = 0.05
)
print(volcano_plt)
```
<div align=center>
  <img src="https://github.com/RuiyangZhai/img/blob/main/DTXRquery/volcano_plt.png?raw=true" width="600">
</div>

__Boxplot__

Compare the expression or values of a specific feature across clinical groups:
```R
# Plot the abundance of Akkermansia across treatment response groups
box_plt <- client$plot_boxplot(
  dataset = "DTXR600031", 
  feature_type = "Microbial Abundance", 
  feature_name = "g__Akkermansia", 
  group_col = "Response"
)
print(box_plt)
```
<div align=center>
  <img src="https://github.com/RuiyangZhai/img/blob/main/DTXRquery/box_plt.png?raw=true" width="600">
</div>

__ROC Curve__

Evaluate the predictive performance of a biomarker:
```R
# Generate an ROC curve to assess Akkermansia as a predictor for Responder ("R")
roc_plt <- client$plot_roc(
  dataset = "DTXR600031", 
  feature_type = "Microbial Abundance", 
  feature_name = "g__Akkermansia", 
  group_col = "Response", 
  positive_class = "R"
)
print(roc_plt)
```
<div align=center>
  <img src="https://github.com/RuiyangZhai/img/blob/main/DTXRquery/roc_plt.png?raw=true" width="600">
</div>

## Contact
Any technical question please contact Ruiyang Zhai (zhairuiyang@foxmail.com) or Te Ma (mate.compbio@foxmail.com).
