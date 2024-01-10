# Copyright (c) 2020, ETH Zurich


#' Allows the user to define the ecological consequences for species
#' within each site, defining thus species survival and abundance
#'
#' @details The arguments of the function allows to apply abiotic and
#'   biotic ecological rules to species in each site. Based on those
#'   rules, the function updates the abundance of each species in each
#'   site. If the abundance is null, the species is absent or
#'   extinct. Ecology can account for local environmental conditions,
#'   the abundance of species, and/or their traits.
#'
#' @param abundance a named vector of abundances with one abundance
#'   value per species
#' @param traits a named matrix containing the species traits, one row
#'   per species
#' @param local_environment the environmental values for the given
#'   site
#' @param config the config of the simulation
#'
#' @return an abundance vector with the new abundance values for every
#'   species. An abundance value of 0 indicates species death, any
#'   other value indicates survival.
#' @export
apply_ecology <- function(abundance, traits, local_environment, config) {
  stop("this function documents the user function interface only, do not use it.")
}



#' Orchestrates for applying the ecology function to all sites
#'
#' @details The ecology is applied on a per site basis over all
#'   species occurring in each site.  Therefore, this function
#'   iterates over all sites and collects the abundance and traits of
#'   any species occurring there.  It then calls the user supplied
#'   apply_ecology function to this collection and apply ecology to
#'   each site.
#'
#' @param config the general config of the simulation
#' @param data the general data list
#' @param vars the general variables list
#' @param cluster a PSOCK cluster as returned by parallel::makeCluster()
#'
#' @return returns the standard val(config, data, vars) list
#' @importFrom parallel parLapply
#' @noRd
loop_ecology <- function(config, data, vars, cluster = NULL) {
  # skip ecology function if config$exp$enable_eco_mec is FALSE
  if (config$gen3sis$general$verbose >= 3) {
    cat(paste("entering ecology module @ time", vars$ti, "\n"))
  }

  all_cells <- rownames(data$landscape$environment)
  abund_matrix <- do.call(cbind, lapply(data$all_species, function(sp) {
    sp_abund <- numeric(length(all_cells))
    sp_abund[all_cells %in% names(sp$abundance)] <- sp$abundance
    sp_abund
  }))
  rownames(abund_matrix) <- all_cells

  # Take ids that have at least one species
  occupied_cells <- rownames(abund_matrix)[rowSums(abund_matrix) > 0]

  ## Check if the maximum number of species per cell has been exceeded
  max_n_sp_idi <- config$gen3sis$general$max_number_of_coexisting_species
  richness <- rowSums(abund_matrix > 0)
  sites_with_excess <- occupied_cells[richness > max_n_sp_idi]
  if (length(sites_with_excess) > 0) {
    vars$flag <- "max_number_of_coexisting_species"
    ## Print only the first site that exceeds the maximum number of
    ## co-occurring species
    paste0(
      "Maximum number of species per cell reached. Specifically ",
      richness[richness > max_n_sp_idi][1], "(>", max_n_sp_idi,
      ") species @ t", vars$ti,
      " site id ", sites_with_excess[1]
    )
    return(list(config = config, data = data, vars = vars))
  }

  ## Get environment and species (abundance, traits) information for
  ## all occupied sites
  sites_with_info <- lapply(
    occupied_cells,
    get_site_info,
    abund_matrix, data$all_species, data$landscape$environment, config
  )

  if (!is.null(cluster)) {
    new_abund_matrix <- parLapply(
      cluster,
      sites_with_info,
      update_site_abundances,
      length(data$all_species), config
    )
  } else {
    new_abund_matrix <- lapply(
      sites_with_info,
      update_site_abundances,
      length(data$all_species), config
    )
  }
  ## The columns of this matrix are sites, the rows are species. It is
  ## an abundance matrix of the occupied sites, not a global matrix.
  new_abund_matrix <- do.call(cbind, new_abund_matrix)
  colnames(new_abund_matrix) <- occupied_cells

  ## Update abundances in the list of all species
  for (spi in seq_len(nrow(new_abund_matrix))) {
    sp_new_abund <- new_abund_matrix[spi, ]
    sp_new_abund <- sp_new_abund[sp_new_abund > 0]
    data$all_species[[spi]]$abundance <- sp_new_abund
  }

  ## Update species distributions
  data$all_species <- lapply(data$all_species, function(sp) {
    limit_species_to_cells(sp, names(sp$abundance)[sp$abundance != 0])
  })

  if (config$gen3sis$general$verbose >= 3) {
    cat(paste("exiting ecology module @ time", vars$ti, "\n"))
  }
  return(list(config = config, data = data, vars = vars))
}

##' Get data that is required to run the ecology function for a site
##'
##' For the given site, retrieve some of the data that is required to
##' run the ecology function: a vector of abundances and a matrix of
##' traits of all species present, and the environmental variables for
##' that site.
##'
##' @param site_id the site identifier
##' @param abund_matrix the global species abundance matrix
##' @param all_species a list with all the species in the simulation
##' @param environment a data.frame containing the landscape
##'   environmental variables, with sites as rows and environmental
##'   variables as columns
##' @param config the general config of the simulation
##'
##' @return a list with the following elements:
##'
##' - site_id: the provided site ID
##'
##' - abundance: abundance vector of the species present in the site
##'
##' - traits: trait matrix of the species present in the site, with
##' species as rows and traits as columns
##'
##' - local_env: environmental variables of the site
##' @noRd
get_site_info <- function(site_id, abund_matrix, all_species,
                          environment, config) {
  coo_sp <- which(abund_matrix[site_id, ] > 0)
  ## Abundance vector for this site
  abundance <- abund_matrix[site_id, coo_sp]
  names(abundance) <- coo_sp

  ## Create trait matrix of the species present in this cell
  traits <- matrix(
    nrow = length(coo_sp),
    ncol = length(config$gen3sis$general$trait_names)
  )
  for (i in seq_along(coo_sp)) {
    traits[i, ] <- all_species[[coo_sp[i]]]$traits[
      site_id, config$gen3sis$general$trait_names
    ]
  }
  colnames(traits) <- config$gen3sis$general$trait_names
  rownames(traits) <- coo_sp

  local_env <- environment[site_id, , drop = FALSE]

  list(
    site_id = site_id, abundance = abundance,
    traits = traits, local_env = local_env
  )
}

##' Update species abundances at one site by applying the ecology
##' function specified in the configuration file
##'
##' @param site_with_info a list containing all information that is
##'   necessary from a site to apply the ecology function,
##'   i.e. abundances from the species that are present, species
##'   traits, and local environmental variables
##' @param n_sp_total total number of species, including extinct ones
##' @param config the general config of the simulation
##'
##' @return A vector of updated species abundances for a site
##' @noRd
update_site_abundances <- function(site_with_info, n_sp_total, config) {
  new_abd <- config$gen3sis$ecology$apply_ecology(
    site_with_info$abundance,
    site_with_info$traits,
    site_with_info$local_env,
    config
  )

  new_abd_site <- numeric(n_sp_total)
  new_abd_site[as.integer(names(site_with_info$abundance))] <- new_abd
  new_abd_site
}
