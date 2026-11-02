library (dplyr)

jaccard_enrichments <- function(enrich1, enrich2){

    v_1 <- c(enrich1$ID)
    v_2 <- c(enrich2$ID)

    cat("enrichment1: ", length(unique(v_1)), "terms\n") 
    cat("enrichment2: ", length(unique(v_2)), "terms\n") 

    return(length(intersect(v_1, v_2)) / length(union(v_1, v_2)))
}

jac_com <- function(com1, com2){
    return(length(intersect(com1, com2)) / length(union(com1, com2)))
}

jaccard_communities <- function(enrich1, enrich2){

    communities1 <- split(enrich1$ID, enrich1$Cluster)
    communities2 <- split(enrich2$ID, enrich2$Cluster)

    cat("communities1: ", length(communities1), "\n")
    cat("communities2: ", length(communities2), "\n")

    m <- matrix(
        NA_real_, 
        nrow = length(communities1), 
        ncol = length(communities2),
        dimnames = list(
            names(communities1), 
            names(communities2)
        )
    )

    for (i in names(communities1)){
        for (j in names(communities2)){
            m[i,j] <- jacc_com(comunidades1[[i]], comunidades2[[j]])
        }   
    }
    return(m)
}

resume_functional_matrix_jac_com <- function(m, threshold = 0.5){

    j_max <- apply(m, 1, max)

    n_communities <- nrow(m)
    n_exact <- sum(apply(m, 1, function(row) any(row == 1)))
    n_threshold <- sum(j_max >= threshold)

    return(data.frame(
        communities = n_communities,
        J_1 = n_exactas,
        J_threshold = n_umbral,
        mean_J_max = mean(j_max)
    ))
}