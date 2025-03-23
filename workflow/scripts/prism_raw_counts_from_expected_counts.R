####Change to inlcude only 6175 genes from CosMx panel 


suppressMessages(library('tidyverse'))
suppressMessages(library('dplyr'))

### SNAKEMAKE I/O ###
raw_expected_counts <- snakemake@input[["raw_expected_counts"]]
ccle_default_line <- snakemake@input[["ccle_default_line"]]
protein_coding_genes <- snakemake@input[["protein_coding_genes"]]
cosmx_genes <- snakemake@input[["cosmx_genes"]]
cell_lines_annotation  <- snakemake@input[["cell_lines_annotation"]]
raw_gene_counts <- snakemake@output[["raw_gene_counts"]]

### SNAKEMAKE PARAMS ###
coding_genes <- as.logical(snakemake@params["coding_genes_only"])

## Load CCLE raw counts
ccle_counts <- data.table::fread(raw_expected_counts)

## Load cell line equivalences
cell_lines <- data.table::fread(ccle_default_line)

## Keep RNA cell line equivalences
cell_lines <- cell_lines %>%
  filter(ProfileType == "rna") %>%
  rename(V1 = ProfileID) %>%
  select(V1, ModelID) %>%
  unique()

## Convert CCLE raw counts to a matrix where each gene is a row and each sample, 
## a column.
## Input is RSEM expected counts (thus I have to round up to units)
genes <- colnames(ccle_counts)[-1]
ccle_counts <- ccle_counts %>%
  merge(cell_lines, by = "V1") %>%
  column_to_rownames("ModelID") %>%
  select(-V1) %>%
  t() %>%
  as.data.frame()
rownames(ccle_counts) <- genes

## Most of the genes with NA values are ERCC-0000X spike-ins
ccle_counts <- na.omit(ccle_counts)

## BROAD format for gene names is HUGO (ENSEMBL)
ccle_counts <- ccle_counts %>%
  mutate(gene_variance = unname(apply(ccle_counts, MARGIN = 1, FUN = var)),
         Gene = str_remove(rownames(ccle_counts), pattern = " \\(.*"),
         Ensembl = str_extract(rownames(ccle_counts), 
                               pattern = "ENSG[0-9]+")) %>%
  ## Get rid of remaining spikes
  filter(!grepl(Gene, pattern = "ERCC-", fixed = TRUE))

## If coding_genes = TRUE, just keep protein coding genes
if (coding_genes) {
  hgnc_coding_genes <- data.table::fread(protein_coding_genes) %>%
    rename(Gene = symbol, Ensembl = ensembl_gene_id) %>%
    select(Gene, Ensembl) %>%
    unique()
  ccle_counts <- ccle_counts %>%
    select(-Gene) %>%
    merge(hgnc_coding_genes)
}

## Keep the most variable gene when two ENSEMBL ids point to the same HGNC
ccle_counts <- ccle_counts %>%
  group_by(Gene) %>%
  filter(gene_variance == max(gene_variance)) %>%
  column_to_rownames("Gene") %>%
  select(-gene_variance, -Ensembl)

## Round counts
ccle_counts <- round(as.matrix(ccle_counts), digits = 0)
cosmx_genes_list <- read.csv(cosmx_genes)
cosmx_genes_list <- cosmx_genes_list$Display.Name

gene_mapping <- c("HLA-DRB1" = "HLA-DRB",
                  "MT-CO1" = "COX1",
                  "MT-CO2" = "COX2",
                  "TMT1A" = "METTL7A",
                  "NHERF1" = "SLC9A3R1",
                  "IFT70A" = "TTC30A",
                  "FCGR1A" = "FCGR1A/BP",
                  "TNXB" = "TNXA/B")

# Replace row names in ccle_counts using the gene mapping
new_rownames <- ifelse(rownames(ccle_counts) %in% names(gene_mapping),
                       gene_mapping[rownames(ccle_counts)], 
                       rownames(ccle_counts))

# Apply the new row names
rownames(ccle_counts) <- new_rownames

# Define the gene mapping list
slashed_genes <- list(
  "CTAG1A/B" = c("CTAG1A", "CTAG1B"),
  "FKBP1A/C" = c("FKBP1A", "FKBP1C"),
  "HBA1/2" = c("HBA1", "HBA2"),
  "IFNA4/10/17" = c("IFNA4", "IFNA10", "IFNA17"),
  "IFNL2/3" = c("IFNL2", "IFNL3"),
  "KIR2DL1/3" = c("KIR2DL1", "KIR2DL3"),
  "KLRC1/2" = c("KLRC1", "KLRC2"),
  "LILRB3/A6" = c("LILRB3", "LILRA6"),
  "MAGEA3/6" = c("MAGEA3", "MAGEA6"),
  "MZT2A/B" = c("MZT2A", "MZT2B"),
  "SMN1/2" = c("SMN1", "SMN2"),
  "TPSAB1/2" = c("TPSAB1", "TPSB2"),
  "XCL1/2" = c("XCL1", "XCL2"),
  "CT45A" = c("CT45A1", "CT45A2", "CT45A3", "CT45A5", "CT45A6", "CT45A7", "CT45A8", "CT45A9", "CT45A10"),
  "MHC I" = c("HLA-A", "HLA-B", "HLA-C")
)

# Iterate over the slashed_genes list
for (new_gene in names(slashed_genes)) {
  genes_to_combine <- slashed_genes[[new_gene]]  # Genes to combine
  
  # Step 1: Create a new row with NAs if it doesn't exist
  if (!(new_gene %in% rownames(ccle_counts))) {
    # Create a row of NAs and add it to the data frame
    ccle_counts <- rbind(ccle_counts, rep(NA, ncol(ccle_counts)))
    rownames(ccle_counts)[nrow(ccle_counts)] <- new_gene
  }
  
  # Step 2: Subset the rows of ccle_counts to be combined
  rows_to_combine <- ccle_counts[genes_to_combine, , drop = FALSE]
  
  # Step 3: Combine the rows (e.g., sum them)
  combined_row <- colSums(rows_to_combine, na.rm = TRUE)
  
  # Step 4: Replace the NAs in the new row with the combined row
  ccle_counts[new_gene, ] <- combined_row
  
  # Step 5: Remove the old rows that were combined
  ccle_counts <- ccle_counts[!rownames(ccle_counts) %in% genes_to_combine, ]
}

#remove rownames thaat are not present in cosmx
common_genes <- intersect(rownames(ccle_counts), cosmx_genes_list)

# Subset ccle_counts to keep only the common rows
ccle_counts <- ccle_counts[common_genes, , drop = FALSE]

cell_lines_annotation <- read.csv(cell_lines_annotation)

nb_ccle <- cell_lines_annotation$DepMap_ID

# Find the common columns between the column names of ccle_counts and nb_ccle
common_columns <- intersect(colnames(ccle_counts), nb_ccle)

# Subset ccle_counts to keep only the common columns
ccle_counts <- ccle_counts[, common_columns, drop = FALSE]

## Save object
saveRDS(ccle_counts, file = raw_gene_counts)
