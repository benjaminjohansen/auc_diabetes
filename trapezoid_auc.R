# Setup -------------------------------------------------------------------

library("JMbayes2") # remove this when real data is available
library("dplyr")
library("DescTools")
library("tidyverse")

# Trapezoid function ------------------------------------------------------

#' Dataloader function. Update cols_to_use to use your own user defined columns
#'
#' @param input_dataframe A dataframe containing at least a time variable (fx years) and an outcome variable (fx serBilir)
#'
#' @return A dataframe with your loaded data
#'
#' @examples import_data(JMbayes2::pbc2)
import_data <- function(input_dataframe) {
  cols_to_use <- c("id", "years", "year", "serBilir") # alternative include albumin instead of serBilir, both are NA free
  data <- input_dataframe[cols_to_use] %>%
    dplyr::mutate(
      whole_year = ceiling(year)
    )
  # add a column with the whole year
  return(data)
}


#' Calculate AUC cumulated over time stamps for a given x-point (time) and a y-point (outcome variable).
#' Collapses the output based on the last calculated time in the time span
#'
#' @param dataframe A dataframe in long-format containing as minimum a column with x-points and y-points.
#' @param max_timespan (int) The max time span of the time column.
#'
#' @return Returns a dataframe in long format with calculated AUC scores.
#'
#' @examples
calculate_trapezoid_auc <- function(dataframe, max_timespan = 10) {
  # Remove unnecessary timespans
  dataframe <- dataframe %>% dplyr::filter(whole_year <= 10)

  output_list <- vector("list", length(unique(dataframe$id)))
  names(output_list) <- unique(dataframe$id)

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
    temp_dataframe$auc[from_point:vector_size] <- sapply(from_point:vector_size, function(i) {
      DescTools::AUC(
        x = temp_dataframe$year[from_point:i],
        y = temp_dataframe$serBilir[from_point:i],
        method = "trapezoid"
      )
    })

    # Collapse datframe to whole years only
    temp_dataframe <- temp_dataframe %>%
      dplyr::select(id, years, year, serBilir, whole_year, auc) %>%
      dplyr::group_by(id, years, whole_year) %>%
      dplyr::summarise(
        last_meassure_year = last(year),
        serBilir = last(serBilir),
        auc = last(auc),
        .groups = "drop"
      ) %>%
      dplyr::rename(year = whole_year)

    output_list[[as.character(unique_id)]] <- temp_dataframe
  }
  output_dataframe <- do.call(rbind, output_list)

  return(output_dataframe)
}

# TODO make a melt function to turn into wide format

#' Converts a dataframe from long to wide format. Assumes to have columns named 0 to n.
#'
#' @param long_dataframe input dataframe in long-format
#' @param max_timespan How far the ordering should go
#'
#' @return A wide-format dataframe
#'
#' @examples wide_dataframe <- create_wide_dataframe(long_dataframe = data_test, max_timespan = max_time)
create_wide_dataframe <- function(long_dataframe, max_timespan = 10) {
  ordered_cols <- c("id", 0:max_timespan)

  wide_dataframe <- long_dataframe %>%
    dplyr::select(id, year, auc) %>%
    tidyr::pivot_wider(names_from = year, values_from = auc)

  reordered_dataframe <- wide_dataframe %>% dplyr::select(all_of(ordered_cols))

  return(reordered_dataframe)
}


# Execute functions -------------------------------------------------------
# variables
max_time <- 10

data <- import_data(JMbayes2::pbc2) # For the PBC2 dataset the years = 0 value is the same as years[0] + 1 in most cases.

data_test <- calculate_trapezoid_auc(dataframe = data, max_timespan = max_time)

wide_dataframe <- create_wide_dataframe(long_dataframe = data_test, max_timespan = max_time)
