#'Determines if a species is considered "valid" by looking at its latin and
#'common names. We also add a BirdTree.org column to facilitate extraction of
#'phylogenetic and functional trait information
#'
#'@importFrom magrittr %>%
#'
#'@param species.df A \code{data.frame} containing results from a
#'\code{\link[pullBBS]{pullBBS}} object. Note that an internet connection is needed for
#'downloading the latest version of the BirdTree / BirdLife taxonomy file
#'
#'@seealso \code{\link[pullBBS]{pullBBS}}
#'
#'@details 'Species' names from BBS observations may not represent valid species. This
#'function uses pattern matching in attempts to classify these vague species into valid ones.
#'It will also match names to BirdTree.org / BirdLife taxonomy to facilitate
#'collection of species' phylogenetic and functional traits.
#'
#'@export
#'
validateSpecies = function(species.df){

  #### Matching names to USGS format to remove non-valid species ####
  # Add leading zero for three-digit AOUs to match USGS format
  species.df$AOU <- withr::with_options(c(scipen = 999),
                                      stringr::str_pad(species.df$AOU, 4, pad = "0"))
  species.df$AOU.char <- as.character(species.df$AOU)

  # Read in species names and AOU codes from the USGS dataset
  data("USGS.bird.names")
  bird.spec.names <- USGS.bird.names

  # Replace missing species names in USGS dataset
  for(i in 1:length(bird.spec.names$Scientific.Name)){
    bird.spec.names$Scientific.Name[i] <- ifelse(bird.spec.names$Scientific.Name[i] == "",
                                                bird.spec.names$Comments[i],
                                                bird.spec.names$Scientific.Name[i])
  }

  # Remove subspecies names from three-name trinomials
  bird.spec.names$Scientific.Name <- sapply(
    strsplit(bird.spec.names$Scientific.Name, " "),
    function(x){
      paste(x[1:2], collapse = " ")
    }
  )

  # Join scientific names to BBS data based on character AOU column
  species.df <- species.df %>%
    dplyr::left_join(bird.spec.names, c("AOU.char" = "Species.Number"))

  #### Matching with BirdTree.org taxonomy for collecting phylogenetic / functional traits ####
  # Load the datafile containing known synonyms and their verified BirdTree names
  data("BirdTree.syns")
  BirdTree.syns <- BirdTree.syns[, c("Synonym", "BirdTree.name")]

  # Download the latest version of the BirdTree / BirdLife taxonomy
  temp <- tempfile()
  download.file('https://data.vertlife.org/birdtree/BLIOCPhyloMasterTax.csv',
                temp)
  True.birdtree.tax <- read.csv(temp)

  # Add column for BirdTree species name and remove species that don't match to BirdTree
  species.df$Scientific.Name <- stringr::str_replace(species.df$Scientific.Name,
                                                     " ", "_")

  # Check all bird names against known BirdTree tip labels and against synonyms;
  # replace with BirdTree names or leave as 'unknown'
  species.df = species.df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(BirdTree.species = ifelse(Scientific.Name %in% True.birdtree.tax$TipLabel,
                                  True.birdtree.tax$TipLabel[True.birdtree.tax$TipLabel %in% Scientific.Name],
                                  "unknown")) %>%
    dplyr::mutate(BirdTree.species = ifelse(Scientific.Name %in% BirdTree.syns$Synonym,
                                  BirdTree.syns$BirdTree.name[BirdTree.syns$Synonym %in% Scientific.Name],
                                  BirdTree.species))

  #### Filter out the unknowns and return ####
  species.df <- species.df[!species.df$BirdTree.species == "unknown",]
  return(species.df)
}

