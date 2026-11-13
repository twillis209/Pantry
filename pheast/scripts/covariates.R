library(impute)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
source(file.path(script_dir, "fusion_twas", "utils", "plink_utils.R"))

load_geno <- function(filename) {
    geno <- read_plink(filename)
    geno_mat <- t(geno$bed)
    colnames(geno_mat) <- geno$fam[, 2]
    geno_mat
}

# PCs of the rows (features) of df, or NULL if fewer than 2 rows have non-zero
# variance: n_pcs is capped at one less than the number of PCs, so 1 feature gives
# none (upstream selected 1:0 and returned the feature itself as a PC).
get_PCs <- function(df, n_pcs) {
    if (nrow(df) < 2) {
        return(NULL)
    }
    if (sum(is.na(df)) > 0) {
        df <- impute.knn(df)$data # Expects samples as columns
    }
    df <- df[apply(df, 1, var) != 0, , drop = FALSE]
    if (nrow(df) < 2) {
        return(NULL)
    }
    pca <- prcomp(t(df), center = TRUE, scale = TRUE)
    n_pcs <- min(n_pcs, ncol(pca$x) - 1)
    pcs <- round(pca$x[, 1:n_pcs, drop = FALSE], 6)
    pcs_df <- data.frame(ID = colnames(pcs), t(pcs))
    colnames(pcs_df) <- c("ID", rownames(pcs))  # Column names starting with digits get 'fixed' and must be changed back
    pcs_df
}

args <- commandArgs(trailingOnly = TRUE)
GENO_PREFIX <- args[1]
BED_FILE <- args[2]
N_GENO_PCS <- as.integer(args[3])
N_PHENO_PCS <- as.integer(args[4])
OUT_FILE <- args[5]

pheno <- read.delim(BED_FILE, check.names = FALSE, row.names = 4)[, -(1:3)]
if (ncol(pheno) < 2) stop("Computing covariate PCs requires more than 1 sample.")
pheno_pcs <- get_PCs(pheno, N_PHENO_PCS)

geno <- load_geno(GENO_PREFIX)
geno <- geno[, colnames(pheno)]
geno_pcs <- get_PCs(geno, N_GENO_PCS)
if (is.null(geno_pcs)) stop("Computing genotype PCs requires at least 2 variants with non-zero variance.")
geno_pcs$ID <- paste("geno", geno_pcs$ID, sep = "_")

if (is.null(pheno_pcs)) {
    # Expected for small modalities: no phenotype PCs exist, so none are used.
    message(
        "covariates.R: no phenotype PCs (", nrow(pheno), " phenotype(s), ",
        "fewer than 2 with non-zero variance); writing genotype PCs only"
    )
    covars <- geno_pcs
} else {
    pheno_pcs$ID <- paste("pheno", pheno_pcs$ID, sep = "_")
    stopifnot(identical(colnames(geno_pcs), colnames(pheno_pcs)))
    covars <- rbind(geno_pcs, pheno_pcs)
}
write.table(covars, OUT_FILE, sep = "\t", quote = FALSE, row.names = FALSE)
