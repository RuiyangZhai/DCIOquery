# DCIOquery

**DCIOquery** is an intuitive, R6-based client package designed to interact with the DiseaseCIO database. It streamlines the entire bioinformatics workflow from data retrieval to visualization. 

With `DiseaseCIO`, you can easily:
* **Query & Filter**: Search the global metadata manifest using specific biological and clinical criteria (e.g., disease type, omics data, treatment regimens).
* **Search & Browse**: Use the DiseaseCIO database search and browsing functions.
* **Batch Download**: Automatically download selected clinical data, feature matrices, and differential analysis results to your local machine.
* **In-Memory Analysis**: Load targeted datasets directly into the R environment as structured lists.
* **Built-in Visualization**: Generate high-quality, publication-ready plots (Volcano plots, Boxplots, and ROC curves) with simple, one-line commands.

## Installation

You can install the development version of DCIOquery from GitHub using the `devtools` or `remotes` package:

```R
# install.packages("devtools")
devtools::install_github("RuiyangZhai/DCIOquery")
```

## Quick Start / Examples
Here is a typical workflow demonstrating how to initialize the client, filter data, load it into memory, and generate plots.

#### 1. Initialize the Client
Creating a new `DCIOquery` object will automatically fetch the latest data manifest from the server.
```R
library(DCIOquery)

# Initialize the client
client <- DCIOquery$new()
# Connecting to server and fetching metadata table...
# Downloading: 310 kB

# View the full metadata
head(client$metadata)

```
#### 2. Search & Browse
You can utilize the search and browsing functionalities of the DiseaseCIO database locally.
```R
##Search
search_res <- client$search_DCIO("Taxonomy Abundance","g__Barnesiella",threshold = "FDR<0.2")
# Connecting to DiseaseCIO...
# Searching data...
print(search_res)

##Browse
browse_res <- client$browse_DCIO("D310042","Gene-level Methylation","ART3")
# Connecting to DiseaseCIO...
# Querying data...
print(search_res)
```
#### 3. Filter and Load Data
You can chain methods together to filter the required data sets, download them and load them into memory for analysis.
```R
# Filter and screen specific omics, diseases, drugs, etc
client$filter_metadata(omic = c("Transcriptomics"),
                       disease = c("Inflammatory bowel disease"),
                       treatment = c("Golimumab"),
                       feature = c("Clinical Data","Gene Expression"))
# Filtered down to 6 records.

# Check the filtered subset table
head(client$sub_table)

# Batch download files
your_save_dir = "my_data"
client$download_files(your_save_dir)
# Downloading D110204_pdata.csv...
# Downloading: 360 B     
# Downloading D110204_Differential_Gene.csv...
# Downloading: 450 kB     
# Downloading D110204_Gene.csv...
# Downloading: 13 MB     
# Downloading D110205_pdata.csv...
# Downloading: 420 B     
# Downloading D110205_Differential_Gene.csv...
# Downloading: 390 kB     
# Downloading D110205_Gene.csv...
# Downloading: 9.2 MB     
# 
# Batch download complete.

# Filter specific data sets and feature types
client$filter_metadata(
  dataset = c("D610044","D530001"),
  feature = c("Clinical Data","Taxonomy Abundance")
)
# Filtered down to 6 records.

# Batch download files and load into memory
client$download_files(your_save_dir)$load_to_memory(your_save_dir)
# File found, loading D530001_pdata.csv into memory...
# File found, loading D530001_Differential_Taxonomy_Abundance.csv into memory...
# File found, loading D530001_Taxonomy_Abundance.csv into memory...
# File found, loading D610044_pdata.csv into memory...
# File found, loading D610044_Differential_Taxonomy_Abundance.csv into memory...
# File found, loading D610044_Taxonomy_Abundance.csv into memory...
# 
# All selected files successfully loaded into memory.
```

#### 4. Data Visualization
Once the data is loaded into memory, you can use the built-in visualization methods.

__Pie Plot__
Pie chart for visualizing clinical baseline characteristics：
```R
# View the proportion of different response groups in the dataset
pie_plt <- client$plot_pie(dataset = "D610044",clinical_col = "Response")

# View gender ratio of patients
# pie_plt <- client$plot_pie(dataset = "D610044",clinical_col = "Sex")
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
  dataset = "D610044", 
  feature_type = "Taxonomy Abundance",
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
  dataset = "D610044", 
  feature_type = "Taxonomy Abundance", 
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
  dataset = "D610044", 
  feature_type = "Taxonomy Abundance", 
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
