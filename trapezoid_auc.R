# Setup -------------------------------------------------------------------

library("JMbayes2") # remove this when real data is available
library("dplyr")
library("DescTools")
library("tidyverse")
library("ggplot2")

project_path <- file.path("C:/SDCA/git/trap_auc") #update to your project path
results <- file.path(project_path, "results")
plots <- file.path(project_path, "plots")

setwd(project_path)

# Functions ------------------------------------------------------

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

#' Converts a dataframe from long to wide format. Assumes to have columns named 0 to n.
#'
#' @param long_dataframe input dataframe in long-format
#' @param max_timespan How far the ordering should go
#'
#' @return A wide-format dataframe
#'
#' @examples wide_dataframe <- create_wide_dataframe(long_dataframe = data_test, max_timespan = max_time)
create_wide_dataframe <- function(long_dataframe, time_col, id_col, outcome_col, max_timespan = 10) {
  ordered_cols <- c("id", 0:max_timespan)

  wide_dataframe <- long_dataframe %>%
    dplyr::select({{ id_col }}, {{ time_col }}, {{ outcome_col }}) %>%
    tidyr::pivot_wider(names_from = {{ time_col }}, values_from = {{ outcome_col }})

  reordered_dataframe <- wide_dataframe %>% dplyr::select(all_of(ordered_cols))

  return(reordered_dataframe)
}

#' Calculate AUC for a dataframe given the x (time) points and the outcome (y) points.
#'
#' @param input_dataframe Long dataframe containing the x, and y points
#' @param max_timespan The maximum time to calculate AUC for
#' @param y_outcome The column containing the y points
#' @param recorded_time The column containing the x points
#' @param end_of_follow_up The column containing the end of follow up time
#' @param wide If the output should be in wide format, then TRUE, if set to FALSE then long format
#'
#' @return A dataframe in either long or wide format with calculated AUC
#' @examples
create_auc <- function(input_dataframe, max_timespan = 10, y_outcome, recorded_time, end_of_follow_up, wide = TRUE) {
  # Setup variables, empty tibbles and vectors for holding data
  # Get unique IDs
  unique_ids <- unique(input_dataframe$id)

  # Tibble to hold the output
  output_tibble <- dplyr::tibble()


  # Add a column to the input data with whole years. These will be used to split intervals.
  input_dataframe <- input_dataframe %>%
    dplyr::mutate(
      whole_year = ceiling(year)
    )


  # Main script
  # 1) Split by ID
  # 2) If no recording at time 0, add a first row for the ID at time 0
  # 3) For each interval,
  #   calculate the AUC.
  #   Copy points nearest to the interval.
  #   Stop when reaching stop criterion

  for (id_unique in unique_ids) {
    # Subset the data by ID
    dataframe_subset <- input_dataframe %>% dplyr::filter(id == id_unique)

    ## Create variables to use
    # Create a copy of the last row with data
    last_row <- dataframe_subset %>% dplyr::slice(n())

    # Set the last recorded time
    latest_recording <- ceiling(last_row[[recorded_time]])

    # Set the stop criterion:
    # The minimum of either the max_timespan (function input) or the floored value of end_of_follow_up column.
    stop_criterion <- min(max_timespan, floor(last_row[[end_of_follow_up]]))

    if (stop_criterion > 0) {
      # Use the interval as a counter
      interval <- 1

      # Empty vectors to hold the results
      id_vec <- vector()
      x_time_vec <- vector()
      y_outcome_vec <- vector()


      # If no data recorded at the start, create a row to hold the first recorded value.
      if (dataframe_subset$whole_year[1] > 0) {
        new_row <- dataframe_subset[1, ]
        new_row$whole_year <- 0
        new_row[[recorded_time]] <- 0
        dataframe_subset <- rbind(new_row, dataframe_subset)
      }

      ## Start the loop
      # As long as the interval is less than or equal to the stop criterion, keep calculating AUC
      while (interval <= stop_criterion) {
        temp_auc <- dplyr::tibble()

        # Subset data for the given interval
        auc_data <- dataframe_subset %>%
          dplyr::filter(whole_year <= interval)

        # Create a copy of the last row - will be needed
        last_auc_row <- auc_data %>% dplyr::slice(n())

        # If the time of the last row is less than the current interval time and is NOT 0, then add a row
        # to the subsetted data
        if (last_auc_row[[recorded_time]] < interval & last_auc_row[[recorded_time]] != 0) {
          last_auc_row[[recorded_time]] <- interval # ceiling(last_auc_row[[recorded_time]])
          temp_auc <- rbind(auc_data, last_auc_row)
        } else {
          temp_auc <- auc_data
        }

        # Calculate the AUC
        auc <- DescTools::AUC(temp_auc[[recorded_time]], temp_auc[[y_outcome]])

        # Add the ID, the interval and the calculated AUC to vectors
        id_vec <- append(id_vec, id_unique)
        x_time_vec <- append(x_time_vec, sprintf("Interval 0:%s", interval))
        y_outcome_vec <- append(y_outcome_vec, auc)

        interval <- interval + 1
      }

      # Bind the outcome to the tibble when the over a subject is done
      new_rows <- tibble(id = id_vec, x_time = x_time_vec, y_outcome = y_outcome_vec)
      output_tibble <- bind_rows(output_tibble, new_rows)
    } else {
      # Bind a row with ID and NA values if the stop criterion is less than or equal to 0
      new_rows <- tibble(id = id_unique, x_time = NA, y_outcome = NA)
      output_tibble <- bind_rows(output_tibble, new_rows)
    }
  }

  if (wide == TRUE) {
    output_tibble <- tidyr::pivot_wider(output_tibble, names_from = x_time, values_from = y_outcome)
    col_names <- names(output_tibble)
    if ("NA" %in% col_names) {
      output_tibble <- dplyr::select(output_tibble, -"NA")
    }
    print(col_names)
  }

  return(output_tibble)
}


# Test cases --------------------------------------------------------------
max_time <- 5

test_simple <- function() {
  input_test <- tidyr::tibble(
    id = c(1, 1, 1, 1, 2, 2, 2, 2),
    years = c(6, 6, 6, 6, 4.6, 4.6, 4.6, 4.6),
    year = c(0.5, 1.5, 2, 3.5, 1.5, 3, 3.5, 4.5),
    serBilir = c(1, 2, 2, 1, 1, 2, 1, 0.5)
  )

  output_test <- data.frame(
    id = c(1, 2),
    AUC_T1 = c(1, NA),
    AUC_T2 = c(3, 2),
    AUC_T3 = c(5, 3.75),
    AUC_T4 = c(5.75, 5),
    AUC_T5 = c(6.75, NA)
  )

  test <- create_auc(input_test,
    y_outcome = "serBilir",
    recorded_time = "year",
    end_of_follow_up = "years",
    max_timespan = max_time
  )

  return(test)
}

head(test_simple())


# Test on JMbayes
test_bayes <- function() {
  test_bayes <- import_data(JMbayes2::pbc2) # %>% dplyr::filter(id == 10)
  # debugonce(create_auc)
  bayes_output <- create_auc(test_bayes,
    y_outcome = "serBilir",
    recorded_time = "year",
    end_of_follow_up = "years",
    max_timespan = max_time,
    wide = TRUE
  )

  return(bayes_output)
}

x <- test_bayes()
head(x)
