#'Clean downloaded NOAA storm data and return the top proportion of most damaging
#'storms
#'
#'@importFrom magrittr %>%
#'
#'@param Stormdata.raw A \code{data.table} object of raw NOAA storm data, downloaded from
#'climage.gov
#'
#'@param damage.threshold A positive numeric value stating the percentage of most damaging storms
#'to return. For example, if \code{damage.threshold == 0.25}, then only storms registering in
#'the top 25 percent of total property and crop damages will be returned. Default is to return
#'the top third of most damaging storms \code{damage.threshold == 0.33}
#'
#'@return A \code{data.frame} object
#'
#'@export
#'
cleanStormData = function(Stormdata.raw, damage.threshold){
  if(missing(damage.threshold)){
    damage.threshold <- 0.33
  }

#### Change lat and long data to the correct scale and direction ####
  cat('Scaling coordinates data... \n')
  Storm.data.raw$LATITUDE <- Storm.data.raw$LATITUDE / 100
  Storm.data.raw$LONGITUDE <- (Storm.data.raw$LONGITUDE / 100) * -1
  Storm.data.raw$LONGITUDE_ <- (Storm.data.raw$LONGITUDE_ / 100) * -1
  Storm.data.raw$LATITUDE_E <- Storm.data.raw$LATITUDE_E / 100
  Storm.data.raw$LONGITUDE_E <- Storm.data.raw$LONGITUDE_
  Storm.data.raw$LONGITUDE_ <- NULL

#### Clean damage columns and separate state / county names ####
  cat('Creating adjusted damage columns, cleaning state names... \n')
  utils::data(county.fips, package = "maps")
  Storm.data.trim = Storm.data.raw %>%
    dplyr::mutate(Prop.damage = ifelse(grepl("k", PROPDMGEXP, ignore.case = TRUE, perl = TRUE),
                                       PROPDMG * 1000,
                                       ifelse(grepl("m", PROPDMGEXP, ignore.case = TRUE,
                                                    perl = TRUE),PROPDMG * 1000000,
                                              ifelse(grepl("", PROPDMGEXP, ignore.case = TRUE,
                                                           perl = TRUE), 0, NA))),
                  Crop.damage = ifelse(grepl("k", CROPDMGEXP, ignore.case = TRUE, perl = TRUE),
                                       CROPDMG * 1000,
                                       ifelse(grepl("m", CROPDMGEXP, ignore.case = TRUE,
                                                    perl = TRUE), CROPDMG * 1000000,
                                              ifelse(grepl("", PROPDMGEXP, ignore.case = TRUE,
                                                           perl = TRUE), 0, NA))),
                  Total.damage = Prop.damage + Crop.damage) %>%
    dplyr::mutate(Year = lubridate::year(lubridate::mdy_hms(BGN_DATE)),
                  Month = lubridate::month(lubridate::mdy_hms(BGN_DATE)),
                  County = tolower(COUNTYNAME),
                  State = tolower(STATE),
                  Storm.type = EVTYPE,
                  Storm.length = LENGTH,
                  Storm.width = WIDTH,
                  Latitude = LATITUDE,
                  Longitude = LONGITUDE,
                  Fatalities = as.numeric(FATALITIES),
                  Injuries = as.numeric(INJURIES)) %>%
    dplyr::select(Prop.damage:Injuries) %>%
    dplyr::left_join(county.fips %>%
                       tidyr::separate_("polyname", c("State", "County"), sep = ",") %>%
                       dplyr::mutate_(County = ~ stringr::str_replace(County, ":.+", "")) %>%
                       dplyr::distinct_() %>%
                       purrr::set_names(c('County.code', 'State', 'County'))) %>%
    dplyr::mutate_at(dplyr::vars(County, State), stringr::str_to_title)

#### Use pattern matching to categorise storm types ####
  # Adapted from code:https://rpubs.com/roelicaal/NOAA_Storm_Database

  cat('Categorising storm types using pattern matching... \n')
Storm.data.trim = Storm.data.trim %>%
  dplyr::mutate(Storm.type = dplyr::case_when(
    grepl("TSTM|THUNDERSTORM|LIGHTNING", Storm.type, perl=TRUE, ignore.case=TRUE) ~ 'thunderstorm',
    grepl("WIND|MICROBURST|(?=.*MICRO)(?=.*BURST)", Storm.type, perl=TRUE, ignore.case=TRUE) ~ 'wind',
    grepl("HAIL", Storm.type, ignore.case=TRUE) ~ 'hail',
    grepl("FLOOD|FLD", Storm.type, perl=TRUE, ignore.case=TRUE) ~ 'flood',
    grepl("TORNADO|FUNNEL|WATERSPOUT", Storm.type, perl=TRUE, ignore.case=TRUE) ~ "tornado",
    grepl("SLEET|SNOW", Storm.type, perl=TRUE, ignore.case=TRUE) ~ "sleet and snow",
    grepl("RAIN", Storm.type, ignore.case=TRUE) ~ "rain",
    grepl("SURF|TIDE|SURGE|RIP|CURRENT", Storm.type, perl=TRUE, ignore.case=TRUE) ~ "surftide",
    grepl("ICE|FREEZ|FROST|FROZEN|COLD|CHILL", Storm.type, perl=TRUE, ignore.case=TRUE) ~ "cold",
    grepl("BLIZZARD|(?=.*ICE)(?=.*STORM)|(?=.*SNOW)(?=.*STORM)|(?=.*WINTER)(?=.*STORM)|(?=.*LAKE)(?=.*EFFECT)",
          Storm.type, perl=TRUE, ignore.case=TRUE) ~ "blizzard",
    grepl("DUST", Storm.type, ignore.case=TRUE) ~ "dust",
    grepl("WILDFIRE|(?=.*WILD)(?=.*FIRE)|(?=.*FOREST)(?=.*FIRE)", Storm.type, perl=TRUE, ignore.case=TRUE) ~ "fire",
    grepl("HEAT|WARM", Storm.type, perl=TRUE, ignore.case=TRUE) ~ "heat",
    grepl("(?=.*DRY)(?=.*WARM)|DROUGHT", Storm.type, perl=TRUE, ignore.case=TRUE) ~ "drought",
    grepl("FOG", Storm.type, ignore.case=TRUE) ~ "fog",
    grepl("HURRICANE|TYPHOON|(?=.*TROPICAL)(?=.*STORM)(?=.*depression)",
          Storm.type, perl=TRUE, ignore.case=TRUE) ~ "tropical storm",
    grepl("LANDSLIDE", Storm.type, ignore.case=TRUE) ~ "landslide",
    grepl("AVALANCHE", Storm.type, ignore.case=TRUE) ~ "avalanche",
    TRUE ~ 'other'))

#### Keep storms with above average damage ratings ####
  Storm.data <- Storm.data.trim %>%
    dplyr::filter(grepl("^thunderstorm$|^wind$|^hail$|^tornado$|^rain$|^fire$|^drought$|^tropical storm$",
                        Storm.type, perl = TRUE)) %>%
    dplyr::filter(Total.damage > quantile(Total.damage, damage.threshold)) %>%
    as.data.frame()

Storm.data$Storm.type <- as.factor(as.character(Storm.data$Storm.type))
colnames(Storm.data) <- stringr::str_to_title(colnames(Storm.data))

rm(Storm.data.trim)

#### Final cleaning of state and county names ####
cat('Final processing of State and County names... \n')
storm_data_z = Storm.data
rm(Storm.data)
storm_data_z$State <- state.name[match(toupper(storm_data_z$State), state.abb)]
small_data <- storm_data_z %>%
  dplyr::tbl_df() %>%
  dplyr::select(-County.code) %>%
  dplyr::mutate_(State = ~ stringr::str_to_lower(State),
                 County = ~ stringr::str_to_lower(County),
                 County = ~stringr::str_replace_all(County, "[.'``]", ""))

# Match `County` to county name in `county.fips`
utils::data(county.fips, package = "maps")
county.fips %>%
  tidyr::separate_("polyname", c("State", "County"), sep = ",") %>%
  dplyr::mutate_(County = ~ stringr::str_replace(County, ":.+", "")) %>%
  dplyr::distinct_() %>%
  purrr::set_names(c('County.code', 'State', 'County')) -> county.fips

a <- small_data %>%
  dplyr::left_join(county.fips) %>%
  dplyr::mutate_(County.code = ~ ifelse(State == "district of columbia", 11001, County.code),
                 County.code = ~ ifelse(State == "virginia" & County == "chesapeake", 51550, County.code),
                 County.code = ~ ifelse(State == "alabama" & County == "dekalb", 1049, County.code)) %>%
  dplyr::filter_(~ !is.na(County.code))

matched_data <- a %>%
  dplyr::mutate_(County.code = ~ ifelse(County.code == 49049, NA, County.code)) %>%
  dplyr::distinct()

#### Match by the cleaned county and state names and return ####
storm_data_z <- storm_data_z %>%
  dplyr::select(-County.code, -State, -County) %>%
  dplyr::left_join(matched_data) %>%
  dplyr::distinct() %>%
  dplyr::filter(!is.na(State) & !is.na(County))

return(storm_data_z)
}
