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
data_sub <- dplyr::filter(data, id %in% c(2,3))

calculate_trapezoid_auc <- function(dataframe, max_timespan = 10) {
  output_dataframe <- data.frame()
  for (unique_id in unique(dataframe$id)) {
    # create a temp dataframe
    temp_dataframe <- dataframe %>% dplyr::filter(id == unique_id)
    vector_size <- nrow(temp_dataframe)
    vector_hold <- vector(mode = "list", length = vector_size)

    # Check if the first point is NA
    # TODO: Might need to be updated to overwrite the first point if NA to the second point?!
    if (is.na(dataframe$serBilir[1])) {
      from_point <- 2
      vector_hold[1] <- NaN
    } else {
      from_point <- 1
    }

    # Iteratively loop through each set of x-y coordinates from: from_point to i (max length of data for patient)
    for (i in from_point:vector_size) {
      vector_hold[[i]] <- DescTools::AUC(
        x = temp_dataframe$year[from_point:i],
        y = temp_dataframe$serBilir[from_point:i],
        method = "trapezoid"
      )
    }

    temp_dataframe$auc <- vector_hold

    # Collapse datframe to whole years only
    temp_dataframe <- temp_dataframe %>%
      dplyr::select(id, years, year, serBilir, whole_year, auc) %>%
      dplyr::group_by(id, years, whole_year) %>%
      dplyr::summarise(
        last_meassure_year = last(year),
        serBilir = last(serBilir),
        auc = last(auc)
      ) %>%
      dplyr::rename(year = whole_year)

    output_dataframe <- rbind(output_dataframe, temp_dataframe)
    # return(dataframe)
  }
  return(output_dataframe)
}

data_test <- calculate_trapezoid_auc(data)
