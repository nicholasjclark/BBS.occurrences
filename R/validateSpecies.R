#'Determines if a species is "valid" by looking at its latin and
# common names. Also adds a column for BirdTree.org ID
#'
#'
#'@export
#'
validateSpecies = function(species.df){
  require(devtools)
  require(dplyr)

  #Add leading zero for three-digit AOUs to match USGS format
  species.df$AOU <- withr::with_options(c(scipen = 999),
                                      stringr::str_pad(species.df$AOU, 4, pad = "0"))
  species.df$AOU.char <- as.character(species.df$AOU)

  #Read in species names and AOU codes from the USGS dataset
  data("USGS.bird.names")
  bird.spec.names <- USGS.bird.names

  #Replace missing species names in USGS dataset
  for(i in 1:length(bird.spec.names$Scientific.Name)){
    bird.spec.names$Scientific.Name[i] <- ifelse(bird.spec.names$Scientific.Name[i] == "",
                                                bird.spec.names$Comments[i],
                                                bird.spec.names$Scientific.Name[i])
  }

  #Remove subspecies names from three-name trinomials
  bird.spec.names$Scientific.Name <- sapply(
    strsplit(bird.spec.names$Scientific.Name, " "),
    function(x){
      paste(x[1:2], collapse = " ")
    }
  )

  #Join scientific names to BBS data based on character AOU column
  species.df<-species.df %>%
    dplyr::left_join(bird.spec.names, c("AOU.char" = "Species.Number"))

  #Add column for BirdTree species name and remove species that don't match to BirdTree
  data("BirdTree.syns")
  BirdTree.syns <- BirdTree.syns[, c("Synonym", "BirdTree.name")]

  temp <- tempfile()
  download.file('https://data.vertlife.org/birdtree/BLIOCPhyloMasterTax.csv',
                temp)
  True.birdtree.tax <- read.csv(temp)
  species.df$Scientific.Name <- stringr::str_replace(species.df$Scientific.Name,
                                                     " ", "_")

  #Check all bird names against known Birdtree tip labels and avian synonyms; replace with BirdTree names
  species.df = species.df %>% rowwise() %>%
    dplyr::mutate(BirdTree.species = ifelse(Scientific.Name%in%True.birdtree.tax$TipLabel,
                                  True.birdtree.tax$TipLabel[True.birdtree.tax$TipLabel %in% Scientific.Name],
                                  "unknown")) %>%
    dplyr::mutate(BirdTree.species = ifelse(Scientific.Name%in% BirdTree.syns$Synonym,
                                  BirdTree.syns$BirdTree.name[BirdTree.syns$Synonym %in% Scientific.Name],
                                  BirdTree.species))

  species.df <- species.df[!species.df$BirdTree.species == "unknown",]
  return(species.df)

}

