# Copyright (c) 2020, ETH Zurich

#' Allows the user to define the rate at which geographic clusters accumulate differentiation
#' with each other.
#'
#' @details This function determines the increase in divergence between separated clusters of a species. This function
#' should return either (i) a single value if there is an homogeneous divergence, or (ii) a matrix indicating the divergence that
#' should be accumulated between specific pairwise geographic clusters.
#' 
#' The function can either return a single value or a full cluster by cluster matrix. If only one value is returned it will be used 
#' to increment divergence between any given distinct cluster pairs. If a matrix is returned it has to be in the dimension of
#' cluster x cluster, in which case the divergence values will be increased according to the cluster membership of any cell pairs.
#'
#' For every time step, the divergence between geographic clusters can increase by a defined number. The divergence values can be 
#' scaled optionally using the species or landscape information. For instance, the divergence between clusters could be higher under
#' warmer temperature, or difference in ecological traits could promote faster divergence between clusters.
#' 
#' Oppositely, for every time-step, if cluster are merged their divergence is reduced by one (1). 
#'
#' @param species the species of the current time step
#' @param cluster_indices an index vector indicating the cluster every occupied site is part of
#' @param landscape the landscape of the current time step
#' @param config the config of the simulation
#'
#' @return a single value or a matrix of divergences between all clusters occurring in clusters_indices
#' @export
get_divergence_factor <- function(species, cluster_indices, landscape, config){
  stop("this function documents the user function interface only, do not use it!")
}

#' Orchestrates the speciation of any species alive in the simulation
#'
#' @param config the current config object
#' @param data the current data object
#' @param vars the current vars object
#' @param cluster a PSOCK cluster as returned by makeCluster()
#'
#' @return an expanded species list including all newly created species
#' @importFrom parallel parLapply
#' @noRd
loop_speciation <- function(config, data, vars, cluster = NULL) {
  if (config$gen3sis$general$verbose >= 3) {
    cat(paste("entering speciation module \n"))
  }

  if (!is.null(cluster)) {
    speciation_list <- parLapply(
      cluster,
      data$all_species,
      speciate,
      config, data$landscape, data$distance_matrix
    )
  } else {
    speciation_list <- lapply(
      data$all_species,
      speciate,
      config, data$landscape, data$distance_matrix
    )
  }

  n_new_sp_ti_total <- sum(sapply(speciation_list, function(x) x$n_new_sp))
  ## NOTE (Adrián Castro Insua): I am not completely sure about the
  ## specific function of these two variables, but I set them to their
  ## corresponding values since they were already there
  vars$n_new_sp_ti <- vars$n_new_sp_ti + n_new_sp_ti_total
  vars$n_sp_added_ti <- vars$n_sp_added_ti + n_new_sp_ti_total

  # Append new species to phylogeny table
  new_sp_phy <- list()
  new_sp_phy$Ancestor <- unlist(
    lapply(
      speciation_list,
      function(x) rep(as.numeric(x$parent$id), x$n_new_sp)
    )
  )
  new_sp_phy$Descendent <- seq(
    length(speciation_list) + 1,
    length.out = n_new_sp_ti_total
  )
  new_sp_phy$Speciation.Time <- rep(vars$ti, n_new_sp_ti_total)
  new_sp_phy$Extinction.Time <- rep(vars$ti, n_new_sp_ti_total)
  new_sp_phy$Speciation.Type <- rep("Genetic", n_new_sp_ti_total)

  data$phy <- rbind(data$phy, as.data.frame(new_sp_phy))

  # Create new species
  has_new_sp <- unlist(lapply(speciation_list, function(x) x$n_new_sp > 0))
  speciated <- speciation_list[has_new_sp]
  speciated_n_sp <- sapply(speciated, function(x) x$n_new_sp)

  starting_ids <- length(speciation_list) +
    cumsum(c(1, head(speciated_n_sp, -1)))

  new_species <- vector("list", length = n_new_sp_ti_total)

  new_sp_idx <- 0
  for (speciation_idx in seq_along(speciated)) {
    sp <- speciated[[speciation_idx]]
    ## Create all descendant species of this parent species
    new_sp_ids <- seq(starting_ids[speciation_idx], length.out = sp$n_new_sp)
    for (i in seq_along(new_sp_ids)) {
      new_sp_idx <- new_sp_idx + 1
      tep_clu_gen_index <- i + 1
      new_sp <- create_species_from_existing(
        sp$parent,
        new_sp_ids[i],
        names(sp$parent$abundance[
          sp$gen_clusters == tep_clu_gen_index
        ]),
        config
      )
      new_species[[new_sp_idx]] <- new_sp
    }

    if (config$gen3sis$general$verbose >= 3) {
      cat(paste("[!]   Welcome, Strange Thing   [!] \n"))
      cat(paste(sp$n_new_sp, "speciation event(s) happened \n"))
    }
  }

  ## Restrict parent species distribution if speciation has occurred
  speciation_list[has_new_sp] <- lapply(speciated, function(sp) {
    sp$parent <- limit_species_to_cells(
      species = sp$parent,
      cells = names(sp$parent$abundance[sp$gen_clusters == 1])
    )

    sp
  })

  # Update list of species
  data$all_species <- lapply(speciation_list, function(x) x$parent)
  # Append new species to species list
  data$all_species <- append(data$all_species, new_species)

  if (config$gen3sis$general$verbose >= 3) {
    cat(paste("exiting speciation module \n"))
  }
  if (config$gen3sis$general$verbose >= 3 && vars$n_sp_added_ti > 0) {
    cat(paste(vars$n_sp_added_ti, "new species created \n"))
  }
  return(list(config = config, data = data, vars = vars))
}

##' Simulates the speciation process for a given living species
##'
##' @param species a `gen3sis_species` object
##' @param config the current config object
##' @param landscape a `gen3sis_landscape` object
##' @param distance_matrix the distance matrix that was created with
##'   the landscape
##'
##' @return a list with the following elements:
##'
##' - parent: the parent species (possibly with a modified divergence
##' matrix)
##'
##' - gen_clusters: vector of genetic clusters, indicating to which
##' cluster each site/population belongs (or NULL if there was no
##' speciation)
##'
##' - n_new_sp: number of new species generated
##' @noRd
speciate <- function(species, config, landscape, distance_matrix) {
  ## Check if the species is alive
  if (!length(species[["abundance"]])) {
    return(
      list(
        parent = species,
        gen_clusters = NULL,
        n_new_sp = 0
      )
    )
  }

  # define occupied cells by species
  species_presence <- names(species[["abundance"]])

  ## calling RCPP function to define physical clusters
  if (length(species_presence) == 1) { # check if only one cell is occupied
    clu_geo_spi_ti <- 1
  } else {
    distances <- config$gen3sis$dispersal$get_dispersal_values(
      length(species_presence), species, landscape, config
    )
    permutation <- sample(seq_along(species_presence))
    clu_geo_spi_ti <- Tdbscan_variable(
      distance_matrix[
        species_presence[permutation],
        species_presence[permutation],
        drop = FALSE
      ],
      distances,
      1
    )

    clu_geo_spi_ti <- clu_geo_spi_ti[order(permutation)]
  }

  gen_dist_spi <- decompress_divergence(species[["divergence"]])
  # update genetic distances
  ifactor <- config$gen3sis$speciation$get_divergence_factor(
    species, clu_geo_spi_ti, landscape, config
  )
  gen_dist_spi <- update_divergence(
    gen_dist_spi,
    clu_geo_spi_ti,
    ifactor = ifactor
  )

  gen_dist_spi <- compress_divergence(gen_dist_spi)

  species[["divergence"]] <- gen_dist_spi

  clu_gen_spi_ti_c <- Tdbscan(
    gen_dist_spi$compressed_matrix,
    config$gen3sis$speciation$divergence_threshold,
    1
  )
  clu_gen_spi_ti <- clu_gen_spi_ti_c[gen_dist_spi$index]
  n_new_sp <- max(clu_gen_spi_ti) - 1

  gen_clusters <- if (n_new_sp > 0) {
    clu_gen_spi_ti
  } else {
    NULL
  }

  list(
    parent = species,
    gen_clusters = gen_clusters,
    n_new_sp = n_new_sp
  )
}

#' Updates a given divergence matrix
#'
#' @param gen_dist_spi a divergence matrix
#' @param clu_geo_spi_ti a cluster index
#' @param ifactor the divergence factor by which the clusters distances are to be increased
#'
#' @return an updated divergence matrix
#' @noRd
update_divergence <- function(divergence, cluster_indices, ifactor) {
  #udpate genetic distance
  clusters <- unique(cluster_indices)
  if( length(ifactor) == 1 ) {
    # scalar ifactor
    divergence <- divergence + ifactor
    dfactor <- 1+ifactor
  } else {
    # matrix ifactor
    divergence <- divergence + ifactor[cluster_indices, cluster_indices]
    dfactor <- 1
  }
  for ( i in clusters ){
    #in case they belong to same clusters, subtract -2 (for the default case), to that final diference is -1 given previous addition!
    divergence[cluster_indices == i, cluster_indices == i] <-
      divergence[cluster_indices == i, cluster_indices== i] - dfactor
  }
  #setting -1 to zero. Genetic differences can not be negative
  divergence[divergence < 0] <- 0
  ##end updating genetic distance##
  return(divergence)
}


#' Updates the total number of species
#'
#' @param config the current config object
#' @param data the current data list
#' @param vars the current vars list
#'
#' @return the updated vals list
#' @noRd
update1.n_sp.all_geo_sp_ti <- function(config, data, vars) {
  # update number of species
  vars$n_sp <- vars$n_sp+vars$n_sp_added_ti
  return(list(config = config, data = data, vars = vars))
}


#' Updates the total number of species alive
#'
#' @param config the current config object
#' @param data the current data list
#' @param vars the current vars list
#'
#' @return the updated vals list
#' @noRd
update2.n_sp_alive.geo_sp_ti <- function(config, data, vars) {
  # update number of species alive
  vars$n_sp_alive <- sum( sapply(data$all_species, function(sp){ifelse(length(sp[["abundance"]]), 1, 0) }))
  return(list(config = config, data = data, vars = vars))
}
