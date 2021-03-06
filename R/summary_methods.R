# Summary methods for mizer package

# Copyright 2012 Finlay Scott and Julia Blanchard.
# Copyright 2018 Gustav Delius and Richard Southwell.
# Development has received funding from the European Commission's Horizon 2020 
# Research and Innovation Programme under Grant Agreement No. 634495 
# for the project MINOUW (http://minouw-project.eu/).
# Distributed under the GPL 3 or later 
# Maintainer: Gustav Delius, University of York, <gustav.delius@york.ac.uk>

# Soundtrack: The Definitive Lead Belly

#' Get diet of predator resolved by prey species
#'
#' Calculates the rate at which a predator of a particular species and size
#' consumes biomass of a prey species, where for this purpose we treat the
#' plankton like another prey species. 
#' 
#' This function performs the same integration as
#' \code{getAvailEnergy()} but does not aggregate over prey species, and
#' multiplies by (1-feeding_level) to get the consumed biomass rather than the
#' available biomass. Outside the range of sizes for a predator species the
#' returned rate is zero.
#'
#' @param params A MizerParams object
#' @param proportion If TRUE (default) the function returns the diet as a
#'   proportion of the total consumption rate. If FALSE it returns the 
#'   consumption rate in grams.
#' @param n An array (species x size) with the abundance density of fish
#' @param n_pp A vector with the abundance of plankton
#' @param n_bb A vector with the abundance of benthos
#' @param n_aa A vector with the abundance of algae
#' 
#' @return An array (predator species  x predator size x (prey species + 3 background spectra) )
#' @export
#' 
getDiet <- function(params, 
                    n, 
                    n_pp, 
                    n_bb, 
                    n_aa,
                    proportion = TRUE) {
 
  # The code is based on that for getAvailEnergy(), but on Nov5 2019 has been modified based on the latest mizer master getDiet() by Gustav
  species <- params@species_params$species
  no_sp <- length(species)
#@  no_sp <- dim(n)[1]
#@  no_w <- dim(n)[2]
#@  no_w_full <- length(n_pp)
  no_w <- length(params@w)
  no_w_full <- length(params@w_full)
  
  diet <- array(0, dim = c(no_sp, no_w, no_sp + 3),
                dimnames = list("predator" = species,
                                "w" = dimnames(n)$w,
                                "prey" = c(as.character(species), "plankton", "benthos", "algae")))
  # idx_sp are the index values of object@w_full such that
  # object@w_full[idx_sp] = object@w
  idx_sp <- (no_w_full - no_w + 1):no_w_full
  
  # If the feeding kernel does not have a fixed predator/prey mass ratio
  # then the integral is not a convolution integral and we can not use fft.
  
  if (length(params@ft_pred_kernel_e) == 1) {
    # pred_kernel is predator species x predator size x prey size
    # We want to multiply this by the prey abundance, which is
    # prey species by prey size, sum over prey size. We use matrix
    # multiplication for this. Then we multiply 1st and 3rd 
    ae <- matrix(params@pred_kernel[, , idx_sp, drop = FALSE],
                 ncol = no_w) %*%
      t(sweep(n, 2, params@w * params@dw, "*"))
      diet[, , 1:no_sp] <- ae #new#@
#@    dim(ae) <- c(no_sp, no_w, no_sp)
#@    # We multiply by interaction matrix, choosing the correct dimensions
#@    diet[, , 1:no_sp] <- sweep(ae, c(1, 3), params@interaction, "*")
    
    #### CHECK THIS PART ####
    # Eating the plankton: On May 2 this is corrected NOT to multiple by availability at this stage 
    diet[, , no_sp + 1] <- rowSums(sweep(
      params@pred_kernel, 3, params@dw_full * params@w_full * n_pp, "*"), dims = 2)
    #"*", check.margin = FALSE), dims = 2)
    # Eating the benthos
    diet[, , no_sp + 2] <- rowSums(sweep(
      params@pred_kernel, 3, params@dw_full * params@w_full * n_bb, "*"), dims = 2)
    #     "*", check.margin = FALSE), dims = 2)
    # Eating the algae
    diet[, , no_sp + 3] <- rowSums(sweep(
      params@pred_kernel, 3, params@dw_full * params@w_full * n_aa, "*"), dims = 2)
    # "*", check.margin = FALSE), dims = 2)
    #### 
    
  } 
    else 
  {
    prey <- matrix(0, nrow = no_sp + 3, ncol = no_w_full) #AA - replace +1 with +3
    prey[1:no_sp, idx_sp] <- sweep(n, 2, params@w * params@dw, "*")
    prey[no_sp + 1, ] <- n_pp * params@w_full * params@dw_full
    prey[no_sp + 2, ] <- n_bb * params@w_full * params@dw_full
    prey[no_sp + 3, ] <- n_aa * params@w_full * params@dw_full
    
    ft <- array(rep(params@ft_pred_kernel_e, times = no_sp + 3) *  #AA - replace +1 with +3
                  rep(mvfft(t(prey)), each = no_sp),
                dim = c(no_sp, no_w_full, no_sp + 3))              #AA - replace +1 with +3
    # We now have an array predator x wave number x prey
    # To Fourier transform back we need a matrix of wave number x everything
    ft <- matrix(aperm(ft, c(2, 1, 3)), nrow = no_w_full)
    ae <- array(Re(mvfft(ft, inverse = TRUE) / no_w_full), 
                dim = c(no_w_full, no_sp, no_sp + 3))              #AA - replace +1 with +3
    ae <- ae[idx_sp, , , drop = FALSE]
    ae <- aperm(ae, c(2, 1, 3))
    # Due to numerical errors we might get negative or very small entries that
    # should be 0
    ae[ae < 1e-18] <- 0
    diet[, , 1:(no_sp + 3)] <- ae
    
  }
  
  #print("diet dims")
  #print(dim(diet))
  #print(diet[1,c(100:120),19])
  ## On May2: only now multiply by the interaction and availability 
  inter <- cbind(params@interaction, params@species_params$avail_PP, params@species_params$avail_BB, params@species_params$avail_AA)
  
  #  diet.temp  <- sweep(diet[, , 1:(no_sp + 3), drop = FALSE],
  #                                         c(1, 3), inter, "*")
  #  print(diet.temp[1,c(100:120),19])  
  
  diet[, , 1:(no_sp + 3)] <- sweep(sweep(diet[, , 1:(no_sp + 3), drop = FALSE],
                                         c(1, 3), inter, "*"), 
                                   c(1, 2), params@search_vol, "*")
  #  print("diet dims")
  #  print(dim(diet))
  #  print(diet[1,c(100:120),19]) 
  
  
  # Correct for satiation and keep only entries corresponding to fish sizes
  f <- getFeedingLevel(params, n, n_pp, n_bb, n_aa) ##AA
  fish_mask <- n > 0
  diet <- sweep(diet, c(1, 2), (1 - f) * fish_mask, "*")
  if (proportion) {
    total <- rowSums(diet, dims = 2)
    diet <- sweep(diet, c(1, 2), total, "/")
    diet[is.nan(diet)] <- 0
  }
  return(diet)
}



#' Calculate the SSB of species
#' 
#' Calculates the spawning stock biomass (SSB) through time of the species in
#' the \code{MizerSim} class. SSB is calculated as the total mass of all mature
#' individuals.
#' 
#' @param sim An object of class \code{MizerSim}.
#'   
#' @return An array containing the SSB (time x species)
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' getSSB(sim)
#' }
getSSB <- function(sim) {
    ssb <- apply(sweep(sweep(sim@n, c(2,3), sim@params@psi,"*"), 3, 
                       sim@params@w * sim@params@dw, "*"), c(1, 2), sum) 
    return(ssb)
}


#' Calculate the total biomass of each species within a size range at each time 
#' step.
#' 
#' Calculates the total biomass through time of the species in the
#' \code{MizerSim} class within user defined size limits. The default option is
#' to use the whole size range. You can specify minimum and maximum weight or
#' length range for the species. Lengths take precedence over weights (i.e. if
#' both min_l and min_w are supplied, only min_l will be used).
#' 
#' @param sim An object of class \code{MizerSim}.
#' @param ... Other arguments to select the size range of fish to be used
#'   in the calculation (min_w, max_w, min_l, max_l).
#'
#' @return An array containing the biomass (time x species)
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' getBiomass(sim)
#' getBiomass(sim, min_w = 10, max_w = 1000)
#' }
getBiomass <- function(sim, ...) {
    size_range <- get_size_range_array(sim@params, ...)
    biomass <- apply(sweep(sweep(sim@n, c(2, 3), size_range, "*"), 3,
                           sim@params@w * sim@params@dw, "*"), c(1, 2), sum)
    return(biomass)
}


#' Calculate the total abundance in terms of numbers of species within a size range
#'
#' Calculates the total numbers through time of the species in the
#' \code{MizerSim} class within user defined size limits. The default option is
#' to use the whole size range You can specify minimum and maximum weight or
#' lengths for the species. Lengths take precedence over weights (i.e. if both
#' min_l and min_w are supplied, only min_l will be used)
#' 
#' @param sim An object of class \code{MizerSim}.
#' @param ... Other arguments to select the size range of the species to be used
#'   in the calculation (min_w, max_w, min_l, max_l).
#'
#' @return An array containing the total numbers (time x species)
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' getN(sim)
#' getN(sim, min_w = 10, max_w = 1000)
#' }
getN <- function(sim, ...) {
    size_range <- get_size_range_array(sim@params, ...)
    n <- apply(sweep(sweep(sim@n, c(2, 3), size_range, "*"), 3,
                     sim@params@dw, "*"), c(1, 2), sum)
    return(n)
}


#' Calculate the total yield per gear and species
#'
#' Calculates the total yield per gear and species at each simulation
#' time step.
#'
#' @param sim An object of class \code{MizerSim}.
#'
#' @return An array containing the total yield (time x gear x species)
#' @export
#' @seealso \code{\link{getYield}}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' getYieldGear(sim)
#' }
getYieldGear <- function(sim) {
    biomass <- sweep(sim@n, 3, sim@params@w * sim@params@dw, "*")
    f_gear <- getFMortGear(sim)
    yield_species_gear <- apply(sweep(f_gear, c(1, 3, 4), biomass, "*"),
                                c(1, 2, 3), sum)
    return(yield_species_gear)
}


#' Calculate the total yield of each species
#'
#' Calculates the total yield of each species across all gears at each
#' simulation time step.
#'
#' @param sim An object of class \code{MizerSim}.
#'
#' @return An array containing the total yield (time x species)
#' @export
#' @seealso \code{\link{getYieldGear}}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=10)
#' y <- getYield(sim)
#' }
getYield <- function(sim) {
    # biomass less the first time step
    yield_gear_species <- getYieldGear(sim)
    return(apply(yield_gear_species, c(1, 3), sum))
}


# Helper function that returns an array (no_sp x no_w) of boolean values indicating whether that size bin is within
# the size limits specified by the arguments
# If min_l or max_l are supplied they take precendence over the min_w and max_w
# But you can mix min_l and max_w etc
# Not exported
get_size_range_array <- function(params, min_w = min(params@w), 
                                 max_w = max(params@w), 
                                 min_l = NULL, max_l = NULL, ...) {
    no_sp <- nrow(params@species_params)
    if (!is.null(min_l) | !is.null(max_l))
        if (any(!c("a","b") %in% names(params@species_params)))
            stop("species_params slot must have columns 'a' and 'b' for length-weight conversion")
    if (!is.null(min_l))
        min_w <- params@species_params$a * min_l ^ params@species_params$b
    else min_w <- rep(min_w,no_sp)
    if (!is.null(max_l))
        max_w <- params@species_params$a * max_l ^ params@species_params$b
    else max_w <- rep(max_w,no_sp)
    if (!all(min_w < max_w))
        stop("min_w must be less than max_w")
    min_n <- aaply(min_w, 1, function(x) params@w >= x, .drop = FALSE)
    max_n <- aaply(max_w, 1, function(x) params@w <= x, .drop = FALSE)
    size_n <- min_n & max_n
    # Add dimnames?
    dimnames(size_n) <- list(sp = params@species_params$species, w = signif(params@w,3)) 
    return(size_n)
}

# TODO: Check documentation for summary
#### summary for MizerParams ####
#' Summarize MizerParams object 
#'
#' Outputs a general summary of the structure and content of the object
#' @param object A \code{MizerParams} object.
#' @param ... Other arguments (currently not used).
#'
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears,inter)
#' summary(params)
#' }
setMethod("summary", signature(object = "MizerParams"), function(object, ...) {
    cat("An object of class \"", as.character(class(object)), "\" \n", sep = "")
    cat("Community size spectrum:\n")
    cat("\tminimum size:\t", signif(min(object@w)), "\n", sep = "")
    cat("\tmaximum size:\t", signif(max(object@w)), "\n", sep = "")
    cat("\tno. size bins:\t", length(object@w), "\n", sep = "")
    # Length of background? 
    cat("Background size spectrum:\n")
    cat("\tminimum size:\t", signif(min(object@w_full)), "\n", sep = "")
    cat("\tmaximum size:\t", signif(max(object@w_full)), "\n", sep = "")
    cat("\tno. size bins:\t", length(object@w_full), "\n", sep = "")
    # w range - min, max, number of w
    # w background min max
    # no species and names and wInf,  - not all these wMat, beta, sigma
    # no gears, gear names catching what
    cat("Species details:\n")
    #cat("\tSpecies\t\tw_inf\n")
    #	for (i in 1:nrow(object@species_params))
    #	    cat("\t",as.character(object@species_params$species)[i], "\t\t ",signif(object@species_params$w_inf[i],3), "\n", sep = "")
    print(object@species_params[,c("species","w_inf","w_mat","beta","sigma")])
    cat("Fishing gear details:\n")
    cat("\tGear\t\t\tTarget species\n")
    for (i in 1:dim(object@catchability)[1]){
        cat("\t",dimnames(object@catchability)$gear[i], "\t\t",dimnames(object@catchability)$sp[object@catchability[i,]>0], "\n", sep=" ") 
    }
})


#### summary for MizerSim ####
#' Summarize MizerSim object 
#'
#' Outputs a general summary of the structure and content of the object
#' @param object A \code{MizerSim} object.
#' @param ... Other arguments (currently not used).
#'
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears,inter)
#' sim <- project(params, effort=1, t_max=5)
#' summary(sim)
#' }

setMethod("summary", signature(object = "MizerSim"), function(object, ...){
    cat("An object of class \"", as.character(class(object)), "\" \n", sep = "")
    cat("Parameters:\n")
    summary(object@params)
    cat("Simulation parameters:\n")
    # Need to store t_max and dt in a description slot? Or just in simulation time parameters? Like a list?
    cat("\tFinal time step: ", max(as.numeric(dimnames(object@n)$time)), "\n", sep = "")
    cat("\tOutput stored every ", as.numeric(dimnames(object@n)$time)[2] - as.numeric(dimnames(object@n)$time)[1], " time units\n", sep = "")
})


#' Calculate the proportion of large fish
#' 
#' Calculates the proportion of large fish through time in the \code{MizerSim}
#' class within user defined size limits. The default option is to use the whole
#' size range. You can specify minimum and maximum size ranges for the species
#' and also the threshold size for large fish. Sizes can be expressed as weight
#' or size. Lengths take precedence over weights (i.e. if both min_l and min_w
#' are supplied, only min_l will be used). You can also specify the species to
#' be used in the calculation. This method can be used to calculate the Large
#' Fish Index. The proportion is based on either abundance or biomass.
#' 
#' @param sim An object of class \code{MizerSim}.
#' @param species numeric or character vector of species to include in the
#'   calculation.
#' @param threshold_w the size used as the cutoff between large and small fish.
#'   Default value is 100.
#' @param threshold_l the size used as the cutoff between large and small fish.
#' @param biomass_proportion a boolean value. If TRUE the proportion calculated
#'   is based on biomass, if FALSE it is based on numbers of individuals.
#'   Default is TRUE.
#' @param ... Other arguments to select the size range of the species to be used
#'   in the calculation (min_w, max_w, min_l, max_l).
#'   
#' @return An array containing the proportion of large fish through time
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=10)
#' getProportionOfLargeFish(sim)
#' getProportionOfLargeFish(sim, species=c("Herring","Sprat","N.pout"))
#' getProportionOfLargeFish(sim, min_w = 10, max_w = 5000)
#' getProportionOfLargeFish(sim, min_w = 10, max_w = 5000, threshold_w = 500)
#' getProportionOfLargeFish(sim, min_w = 10, max_w = 5000,
#'     threshold_w = 500, biomass_proportion=FALSE)
#' }
getProportionOfLargeFish <- function(sim, 
                                     species = 1:nrow(sim@params@species_params), 
                                     threshold_w = 100, threshold_l = NULL, 
                                     biomass_proportion=TRUE, ...) {
    check_species(sim,species)
    # This args stuff is pretty ugly - couldn't work out another way of using ...
    args <- list(...)
    args[["params"]] <- sim@params
    total_size_range <- do.call("get_size_range_array", args = args)
    args[["max_w"]] <- threshold_w
    args[["max_l"]] <- threshold_l
    large_size_range <- do.call("get_size_range_array", args = args)
    w <- sim@params@w
    if (!biomass_proportion) # based on abundance numbers
        w[] <- 1
    total_measure <- apply(sweep(sweep(sim@n[,species,,drop=FALSE],c(2,3),total_size_range[species,,drop=FALSE],"*"),3,w * sim@params@dw, "*"),1,sum)
    upto_threshold_measure <- apply(sweep(sweep(sim@n[,species,,drop=FALSE],c(2,3),large_size_range[species,,drop=FALSE],"*"),3,w * sim@params@dw, "*"),1,sum)
    #lfi = data.frame(time = as.numeric(dimnames(sim@n)$time), proportion = 1-(upto_threshold_measure / total_measure))
    #return(lfi)
    return(1 - (upto_threshold_measure / total_measure))
}


#' Calculate the mean weight of the community
#'
#' Calculates the mean weight of the community through time.
#' This is simply the total biomass of the community divided by the abundance in numbers.
#' You can specify minimum and maximum weight or length range for the species. Lengths take precedence over weights (i.e. if both min_l and min_w are supplied, only min_l will be used).
#' You can also specify the species to be used in the calculation.
#'
#' @param sim An object of class \code{MizerSim}
#' @param species numeric or character vector of species to include in the calculation
#' @param ... Other arguments for the \code{getN} and \code{getBiomass} methods such as \code{min_w}, \code{max_w} \code{min_l} and \code{max_l}.
#'
#' @return A vector containing the mean weight of the community through time
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=10)
#' getMeanWeight(sim)
#' getMeanWeight(sim, species=c("Herring","Sprat","N.pout"))
#' getMeanWeight(sim, min_w = 10, max_w = 5000)
#' }
getMeanWeight <- function(sim, species = 1:nrow(sim@params@species_params), ...){
    check_species(sim, species)
    n_species <- getN(sim, ...)
    biomass_species <- getBiomass(sim, ...)
    n_total <- apply(n_species[, species, drop = FALSE], 1, sum)
    biomass_total <- apply(biomass_species[, species, drop = FALSE], 1, sum)
    return(biomass_total / n_total)
}


#' Calculate the mean maximum weight of the community
#'
#' Calculates the mean maximum weight of the community through time. This can be
#' calculated by numbers or biomass. The calculation is the sum of the w_inf *
#' abundance of each species, divided by the total abundance community, where
#' abundance is either in biomass or numbers. You can specify minimum and
#' maximum weight or length range for the species. Lengths take precedence over
#' weights (i.e. if both min_l and min_w are supplied, only min_l will be used).
#' You can also specify the species to be used in the calculation.
#'
#' @param sim An object of class \code{MizerSim}.
#' @param species numeric or character vector of species to include in the calculation.
#' @param measure The measure to return. Can be 'numbers', 'biomass' or 'both'
#' @param ... Other arguments for the \code{getN} and \code{getBiomass} methods such as \code{min_w}, \code{max_w} \code{min_l} and \code{max_l}.
#'
#' @return A matrix or vector containing the mean maximum weight of the community through time
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=10)
#' getMeanMaxWeight(sim)
#' getMeanMaxWeight(sim, species=c("Herring","Sprat","N.pout"))
#' getMeanMaxWeight(sim, min_w = 10, max_w = 5000)
#' }
getMeanMaxWeight <- function(sim, species = 1:nrow(sim@params@species_params), 
                             measure = "both", ...) {
    if (!(measure %in% c("both","numbers","biomass"))) {
        stop("measure must be one of 'both', 'numbers' or 'biomass'")
    }
    check_species(sim, species)
    n_species <- getN(sim, ...)
    biomass_species <- getBiomass(sim, ...)
    n_winf <- apply(sweep(n_species, 2, sim@params@species_params$w_inf,"*")[,species,drop=FALSE], 1, sum)
    biomass_winf <- apply(sweep(biomass_species, 2, sim@params@species_params$w_inf,"*")[,species,drop=FALSE], 1, sum)
    mmw_numbers <- n_winf / apply(n_species, 1, sum)
    mmw_biomass <- biomass_winf / apply(biomass_species, 1, sum)
    if (measure == "numbers")
        return(mmw_numbers)
    if (measure == "biomass")
        return(mmw_biomass)
    if (measure == "both")
        return(cbind(mmw_numbers, mmw_biomass)) 
}


#' Calculate the slope of the community abundance
#'
#' Calculates the slope of the community abundance through time by performing a linear regression on the logged total numerical abundance at weight and logged weights (natural logs, not log to base 10, are used).
#' You can specify minimum and maximum weight or length range for the species. Lengths take precedence over weights (i.e. if both min_l and min_w are supplied, only min_l will be used).
#' You can also specify the species to be used in the calculation.
#'
#' @param sim An object of class \code{MizerSim}.
#' @param species Numeric or character vector of species to include in the calculation.
#' @param biomass Boolean. If TRUE (default), the abundance is based on biomass, if FALSE the abundance is based on numbers. 
#' @param ... Optional parameters include
#'   \itemize{
#'     \item min_w Minimum weight of species to be used in the calculation.
#'     \item max_w Maximum weight of species to be used in the calculation.
#'     \item min_l Minimum length of species to be used in the calculation.
#'     \item max_l Maximum length of species to be used in the calculation.
#'   }
#'
#' @return A data frame with slope, intercept and R2 values.
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=40, dt = 1, t_save = 1)
#' # Slope based on biomass, using all species and sizes
#' slope_biomass <- getCommunitySlope(sim)
#' # Slope based on numbers, using all species and sizes
#' slope_numbers <- getCommunitySlope(sim, biomass=FALSE)
#' # Slope based on biomass, using all species and sizes between 10g and 1000g
#' slope_biomass <- getCommunitySlope(sim, min_w = 10, max_w = 1000)
#' # Slope based on biomass, using only demersal species and sizes between 10g and 1000g
#' dem_species <- c("Dab","Whiting","Sole","Gurnard","Plaice","Haddock", "Cod","Saithe")
#' slope_biomass <- getCommunitySlope(sim, species = dem_species, min_w = 10, max_w = 1000)
#' }
getCommunitySlope <- function(sim, species = 1:nrow(sim@params@species_params),
                              biomass = TRUE, ...) {
    check_species(sim, species)
    size_range <- get_size_range_array(sim@params, ...)
    # set entries for unwanted sizes to zero and sum over wanted species, giving
    # array (time x size)
    total_n <-
        apply(sweep(sim@n, c(2, 3), size_range, "*")[, species, , drop = FALSE],
              c(1, 3), sum)
    # numbers or biomass?
    if (biomass)
        total_n <- sweep(total_n, 2, sim@params@w, "*")
    # previously unwanted entries were set to zero, now set them to NA
    # so that they will be ignored when fitting the linear model
    total_n[total_n <= 0] <- NA
    # fit linear model at every time and put result in data frame
    slope <- adply(total_n, 1, function(x, w) {
        summary_fit <- summary(lm(log(x) ~ log(w)))
        out_df <- data.frame(
            slope = summary_fit$coefficients[2, 1],
            intercept = summary_fit$coefficients[1, 1],
            r2 = summary_fit$r.squared
        )
    }, w = sim@params@w)
    dimnames(slope)[[1]] <- slope[, 1]
    slope <- slope[, -1]
    return(slope)
}


# internal
check_species <- function(object, species){
    if (!(is(species,"character") | is(species,"numeric")))
        stop("species argument must be either a numeric or character vector")
    if (is(species,"character"))
        check <- all(species %in% dimnames(object@n)$sp)  
    if (is(species,"numeric"))
        check <- all(species %in% 1:dim(object@n)[2])
    if (!check)
        stop("species argument not in the model species. species must be a character vector of names in the model, or a numeric vector referencing the species")
    return(check)
}

