# Setup -------------------------------------------------------------------

library("JMbayes2") # remove this when real data is available
library("dplyr")
library("DescTools")

# Trapezoid function ------------------------------------------------------

import_data <- function(input_dataframe) {
  cols_to_use <- c("id", "years", "year", "serBilir") # alternative include albumin instead of serBilir, both are NA free
  data <- input_dataframe[cols_to_use] %>%
    dplyr::mutate(
      whole_year = ceiling(year)
    )

  # add a column with the whole year
  return(data)
}

data <- import_data(JMbayes2::pbc2) # For the PBC2 dataset the years = 0 value is the same as years[0] + 1 in most cases.

# Subset data to make it work for one person
data_sub <- dplyr::filter(data, id %in% c(3))

calculate_trapezoid_auc <- function(dataframe, max_timespan = 10) {
  for (id in dataframe$id) {

  vector_size <- nrow(dataframe)
  vector_hold <- vector(mode = "list", length = vector_size)

  if(is.na(dataframe$serBili[1])) {
    from_point <- 2
    vector_hold[1] <- NaN
  } else {
    from_point <- 1
  }

  for (i in from_point:vector_size) {
    vector_hold[[i]] <- DescTools::AUC(
      x = dataframe$year,
      y = dataframe$serBilir,
      from = dataframe$year[from_point],
      to = dataframe$serBilir[i],
      method = "trapezoid"
    )
  }

  dataframe$auc <- vector_hold

  # Collapse datframe to whole years only
  dataframe <- dataframe %>%
  dplyr::select(id, years, year, serBilir, whole_year, auc) %>%
    dplyr::group_by(id, years, whole_year) %>%
    dplyr::summarise(
      last_meassure_year = last(year),
      serBilir = last(serBilir),
      auc = last(auc)
    ) %>%
    dplyr::rename(year = whole_year)

  return(dataframe)
  }
  }

# debugonce(calculate_trapezoid_auc)
data_test <- calculate_trapezoid_auc(data_sub)
head(data_test)

dplyr::group_by(data_test, id, years) %>%
  summarise(
    year = first(year),
    y = last(serBilir),
    whole_year = max(whole_year),
    auc = last(auc)
)

data_test %>%
  # filter(!(year == 0 & is.na(serBilir))) %>%
  select(id, years, year, serBilir, whole_year) %>%
  group_by(id, years, whole_year) %>%
  summarise(
    year = last(year),
    serBilir = last(serBilir)
  )

calculate_trapezoid_auc <- function(dataframe, max_timespan = 10) {
  dataframe <- data[3:111, ] %>% dplyr::mutate(whole_year = ceiling(year))



  x <- dataframe %>%
    dplyr::mutate(whole_year = ceiling(year)) %>%
    dplyr::group_by(id, years)

  simple <- x %>%
    group_by(id) %>%
    summarise(max_time = max(ceiling(years)), min_year = min(whole_year), max_year = (max(whole_year)))

  for (i in 2:max_timespan) {
    auc_out <- DescTools::AUC(
      x = x$year,
      y = x$serBilir,
      from = x$year[2],
      to = x$year[i],
      method = "trapezoid"
    )
    print(auc_out)
  }


  result <- dataframe %>%
    dplyr::group_by(id) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      AUC = ifelse(
        whole_year > 0,
        DescTools::AUC(
        x = .$year,
        y = .$serBilir,
        from = year[1]
      ),
      NA
      )
    )
result

  print(x)
  DescTools::AUC(
    x$year[1:4],
    x$serBilir[1:4],
    method = "trapezoid"
  )
}

# Assumptions

# Simple; uniform grid

# Complex: non-uniform grid
