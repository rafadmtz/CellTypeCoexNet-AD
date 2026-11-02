library(igraph)
library(aricode)

#' Compare community detection algorithms by pairwise AMI
#'
#' Runs seven community detection algorithms on a graph and computes the
#' Adjusted Mutual Information (AMI) between every pair of partitions. Each
#' algorithm is summarized by its mean AMI against the other six.
#'
#' Modularity (Q) and the number of communities per algorithm are printed
#' as a quick sanity check.
#'
#' @param G An undirected igraph graph.
#' @param seed Integer. Random seed set before running the algorithms.
#'   Default `2`.
#' @return Named numeric vector with the mean AMI of each algorithm, sorted
#'   from highest to lowest.
#' @details Uses `aricode::AMI()`. 1 means identical partitions and values
#'   near 0 mean independent partitions.
compare_community_ami <- function(G, seed = 2) {

  set.seed(seed)

  communities <- list()
  communities[["infomap"]]     <- cluster_infomap(G, nb.trials = 1000)
  communities[["fast_greedy"]] <- cluster_fast_greedy(G)
  communities[["louvain"]]     <- cluster_louvain(G)
  communities[["walktrap"]]    <- cluster_walktrap(G)
  communities[["eigen"]]       <- cluster_leading_eigen(G)
  communities[["prop"]]        <- cluster_label_prop(G)
  communities[["betweenness"]] <- cluster_edge_betweenness(G)

  cat("Q values for each algorithm\n")
  print(round(sapply(communities, modularity), 3))

  cat("Length of communities:\n")
  print(sapply(communities, length))

  community_vectors <- lapply(communities, function(x) as.integer(membership(x)))
  stopifnot(all(sapply(community_vectors, length) == vcount(G)))

  ami_matrix <- matrix(
    NA_real_,
    nrow = length(communities),
    ncol = length(communities),
    dimnames = list(names(communities), names(communities))
  )

  for (p in combn(names(community_vectors), 2, simplify = FALSE)) {
    ami_value <- AMI(community_vectors[[p[1]]], community_vectors[[p[2]]])
    ami_matrix[p[1], p[2]] <- ami_value
    ami_matrix[p[2], p[1]] <- ami_value
  }

  ami_by_algorithm <- sort(rowMeans(ami_matrix, na.rm = TRUE), decreasing = TRUE)

  return(ami_by_algorithm)
}


#' Measure AMI stability of community detection across random seeds
#'
#' Runs each algorithm (louvain, prop, infomap) with several random seeds and
#' computes the mean pairwise AMI between its partitions. Higher values mean
#' the algorithm gives more similar results across seeds.
#'
#' @param G An undirected igraph graph.
#' @param n_seeds Integer. Number of random seeds. Default `10`.
#' @param base_seed Integer. Run `i` uses seed `base_seed + i`. Default `2`.
#' @return Named numeric vector with the mean AMI of each algorithm.
#' @details Infomap uses 100 trials per run.
compare_stability <- function(G, n_seeds = 10, base_seed = 2) {

  algorithms <- list(
    louvain = cluster_louvain,
    prop    = cluster_label_prop,
    infomap = function(g) cluster_infomap(g, nb.trials = 100)
  )

  runs <- list(louvain = list(), prop = list(), infomap = list())

  for (i in seq_len(n_seeds)) {
    set.seed(base_seed + i)
    for (algorithm in names(algorithms)) {
      runs[[algorithm]][[i]] <- as.integer(membership(algorithms[[algorithm]](G)))
    }
  }

  ami_by_algorithm <- sapply(names(runs), function(algorithm) {
    values <- c()
    for (p in combn(seq_len(n_seeds), 2, simplify = FALSE)) {
      values <- c(values, AMI(runs[[algorithm]][[p[1]]], runs[[algorithm]][[p[2]]]))
    }
    mean(values)
  })

  return(ami_by_algorithm)
}


#' Run Infomap and export its communities
#'
#' Runs Infomap with 1000 trials, prints the modularity (Q) and the number of
#' communities, and writes the gene-to-community assignment to a CSV file.
#'
#' @param G An undirected igraph graph. Vertex names are gene identifiers.
#' @param outdir Character. Directory where the CSV is written.
#' @param tag Character. Label added to the output file name, for example a
#'   group name.
#' @return The `communities` object from `cluster_infomap()`.
#' @details Output file: `infomap_communities_<tag>.csv`, with columns `gene`
#'   and `community`.
run_infomap <- function(G, outdir, tag) {

  communities <- cluster_infomap(G, nb.trials = 1000)

  cat("Q value for:", tag, "\n")
  print(round(modularity(communities), 3))

  cat("Length of communities:", length(communities), "\n")

  df <- data.frame(
    gene      = names(membership(communities)),
    community = as.integer(membership(communities))
  )

  out_file <- file.path(outdir, paste0("infomap_communities_", tag, ".csv"))
  write.csv(df, out_file, row.names = FALSE)

  return(communities)
}