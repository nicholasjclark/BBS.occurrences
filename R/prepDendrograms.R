#'Create a multivariate functional similarity dendrogram for species
#'using lists of trait datasets
#'
#'
#'@param func.datasets A \code{list} containing dataframes of trait variables, where
#'rownames of each dataframe represent the species' names and columns represent trait variables.
#'@param prep.types A \code{list} describing which prepping method to apply to each dataset in
#'\code{func.datasets}. Each prep method must match one of the following:
#'\code{'prep.fuzzy'} (for proportional variables), \code{'prep.binary'} (for binary variables),
#'or \code{'Q'} (for continuous variables). For more information about prep methods, see
#'\code{\link[ade4]{dist.ktab}}
#'
#'@seealso \code{\link[ade4]{dist.ktab}}
#'
#'@details Multiple functional trait datasets are incorporated with equal weights to
#'create a mixed-variable coefficient of distance, which is a generalized Gower's coefficient
#'of distance. This matrix is then scaled to range \code{[0, 1]} and returned.
#'
#'@export
#'
prepDendrograms = function(func.datasets, prep.types){

#### Prep the trait datasets for the Gower's distance calculations ####
    prep.dat <- lapply(seq_along(func.datasets), function(j){

      if(prep.types[j] == 'prep.fuzzy') {
        ade4::prep.fuzzy(func.datasets[[j]],
                   col.blocks = ncol(func.datasets[[j]])) }

      else if(prep.types[j] == 'prep.binary') {
        ade4::prep.binary(func.datasets[[j]],
                    col.blocks = ncol(func.datasets[[j]])) }

      else if(prep.types[j] == 'Q') {
        data.frame(func.datasets[[j]]) }
    })
    ktab1 <- ade4::ktab.list.df(prep.dat)

  #### Create a vector of prep types to match ade4 syntax ####
  type <- vector()
  for(i in 1:length(prep.types)){
    type[i] <- if(prep.types[i] == "prep.binary"){'B'}
    else if(prep.types[i] == "prep.fuzzy"){'F'}
    else {'Q'}
  }

  # Use the vector of types to create a vector of scaling options
  option <- vector()
  for(i in 1:length(prep.types)){
    option[i] <- if(prep.types[i] == "prep.binary"){'noscale'}
    else if(prep.types[i] == "prep.fuzzy"){'scaledBYrange'}
    else {'scaledBYsd'}
  }

  #### Calculate the Gower's distance matrix ####
  distrait <- ade4::dist.ktab(ktab1, type = type, option = option)

  # Scale the matrix to range [0, 1] and return
  ecol.dist <- as.matrix(distrait)
  output <- ecol.dist / max(ecol.dist, na.rm = T)
  colnames(output) <- rownames(func.datasets[[1]])
  rownames(output) <- rownames(func.datasets[[1]])

  return(output)
}

