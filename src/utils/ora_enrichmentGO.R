library(clusterProfiler)
library(org.Hs.eg.db)
library(arrow)
library(GOSemSim)



#' Get the gene universe from the edgelist Parquet file
#'
#' Reads only the column names of the Parquet file (one column per gene), so
#' the data is never loaded in memory. Index columns are dropped.
#'
#' @param path Character. Path to the Parquet file whose columns are genes.
#' @return Character vector of unique gene symbols.
#' @details Defensive fix: the `MT.` prefix is restored to `MT-` in case an
#'   upstream tool replaced the dash with a dot.
gene_universe_from_total_measured <- function(path) {

    genes <- setdiff(open_dataset(path)$schema$names, c("", "Unnamed: 0"))
    genes <- sub("^MT\\.", "MT-", genes)

    return(unique(genes))
}

#' Compare GO enrichment across communities
#'
#' Runs GO Biological Process over-representation analysis for every
#' community with `clusterProfiler::compareCluster()`. Communities with fewer
#' than `min_genes` genes are dropped.
#'
#' @param universe Character vector. Background gene universe (gene symbols).
#' @param community_path Character. Path to the CSV with columns `gene` and
#'   `community`.
#' @param min_genes Integer. Minimum number of genes for a community to be
#'   tested.
#' @param rds_path Character or `NULL`. If given, the result is saved there
#'   with `saveRDS()` so it can be reused for plotting. The parent directory
#'   is created if missing. Default `NULL` (nothing is saved).
#' @return A `compareClusterResult` object.
#' @details Uses `keyType = "SYMBOL"`, BH adjustment, `pvalueCutoff = 0.05`
#'   and `qvalueCutoff = 1`.
compare_clusters_GO <- function(universe, community_path, min_genes, rds_path = NULL) {

    df_communities <- read.csv(community_path)

    df_communities$gene <- sub("^MT\\.", "MT-", df_communities$gene)

    communities_list <- split(df_communities$gene, df_communities$community)

    communities_list <- communities_list[lengths(communities_list) >= min_genes]

    result <- compareCluster(
        geneClusters  = communities_list,
        fun           = "enrichGO",
        OrgDb         = org.Hs.eg.db,
        keyType       = "SYMBOL",
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 1,
        universe      = universe
    )

    if (!is.null(rds_path)) {
        dir.create(dirname(rds_path), recursive = TRUE, showWarnings = FALSE)
        saveRDS(result, rds_path)
    }

    return(result)
}

#' Simplify GO enrichment results by semantic similarity
#'
#' Removes redundant GO terms with `clusterProfiler::simplify()`, keeping the
#' term with the lowest adjusted p-value in each group of similar terms.
#'
#' @param go_results A `compareClusterResult` object, for example from
#'   [compare_clusters_GO()].
#' @param cutoff Numeric. Similarity threshold above which terms are
#'   considered redundant. Default `0.7`.
#' @param sem_data A `GOSemSimDATA` object (from `GOSemSim::godata()`).
#'   Default `NULL`.
#' @param rds_path Character or `NULL`. If given, the simplified result is
#'   saved there with `saveRDS()`. The parent directory is created if missing.
#'   Default `NULL` (nothing is saved).
#' @return The simplified `compareClusterResult` object, invisibly.
#' @details Similarity is measured with the Wang method.
simplify_results <- function(go_results, cutoff = 0.7, sem_data = NULL, rds_path = NULL) {
    simplified_results <- clusterProfiler::simplify(
        go_results,
        cutoff     = cutoff,
        by         = "p.adjust",
        select_fun = min,
        measure    = "Wang",
        semData    = sem_data
    )

    if (!is.null(rds_path)) {
        dir.create(dirname(rds_path), recursive = TRUE, showWarnings = FALSE)
        saveRDS(simplified_results, rds_path)
    }

    invisible(simplified_results)
}
