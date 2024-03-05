# Setup -------------------------------------------------------------------

library("JMbayes2") # remove this when real data is available
library("dplyr")
library("DescTools")
library("tidyverse")
library("ggplot2")

project_path <- file.path("C:/SDCA/git/trap_auc")
results <- file.path(project_path, "results")
plots <- file.path(project_path, "plots")

setwd(project_path)

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

#' Extend a dataframe in time. Meaning, if no records at baseline is made, the first entry is copied to baseline, with time a 0
#' Extend dataframe up to end-of-follow up if missing data. Will copy the last observed outcome.
#'
#' @param input_dataframe Dataframe that should be extended
#' @param outcome_variable The column containing the outcome variable (y)
#' @param time_variable The column containing the time vairable (x)
#' @param end_of_study_column The end of study column
#'
#' @return
#' @export
#'
#' @examples
extend_dataframe <- function(input_dataframe, outcome_variable, time_variable, end_of_study_column) {
  # Get last year and last recorded y-outcome:
  max_year <- floor(input_dataframe[[end_of_study_column]][1])
  last_y <- tail(input_dataframe[[outcome_variable]], 1)

  # Create a new x_time variable for calculating AUC
  output_dataframe <- input_dataframe %>%
    group_by(whole_year) %>%
    mutate(x_time := ifelse(
      !!sym(time_variable) == max(!!sym(time_variable)),
      whole_year,
      !!sym(time_variable)
    ))

  # If no data recorded at the start, create a row to hold the first recorded value.
  if (output_dataframe$whole_year[1] != 0) {
    new_row <- output_dataframe[1, ]
    new_row$whole_year <- 0
    new_row[[time_variable]] <- 0
    new_row$x_time <- 0
    output_dataframe <- rbind(new_row, output_dataframe)
  }

  # update the dataframe with the new time_variable
  # output_dataframe$x_time = !!sym(time_variable)
  output_dataframe <- output_dataframe %>%
    group_by(whole_year) %>%
    # Find the last row in each whole year
    dplyr::slice(n()) %>%
    # add x_time
    dplyr::mutate(!!sym(time_variable) := ceiling(!!sym(time_variable))) %>%
    rbind(output_dataframe, .) %>%
    dplyr::arrange(!!sym(time_variable)) %>%
    dplyr::distinct()

  # Keep adding new rows until the end of follow up
  new_rows <- list()
  current_year <- max(output_dataframe$whole_year)
  while (current_year < max_year) {
    current_year <- current_year + 1
    new_row <- output_dataframe[nrow(output_dataframe), ]
    new_row$whole_year <- current_year
    new_row[[time_variable]] <- current_year
    new_row[[outcome_variable]] <- last_y
    new_rows[[length(new_rows) + 1]] <- new_row
  }
  output_dataframe <- rbind(output_dataframe, do.call(rbind, new_rows))

  return(output_dataframe)
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
# TODO: Update to take columns. See extend dataframe for scafolding/inspiration
calculate_trapezoid_auc <- function(dataframe, max_timespan = NULL) {
  # Remove unnecessary timespans
  if (max_timespan) {
    dataframe <- dataframe %>% dplyr::filter(whole_year <= max_timespan)
  }

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


# Setting up a test case --------------------------------------------------
# TODO: Agree on how we calculate AUC in different time intervals.
# debugonce(extend_dataframe)
id_2 <- data %>% dplyr::filter(id == 2)
id_4 <- data %>% dplyr::filter(id == 4 & whole_year > 0)

test_id <- import_data(pbc2) %>% dplyr::filter(id == 2)
test_out <- extend_dataframe(test_id, outcome_variable = "serBilir", time_variable = "year", end_of_study_column = "years")


x <- extend_dataframe(input_dataframe = id_2, outcome_variable = "serBilir", time_variable = "year", end_of_study_column = "years")
y <- extend_dataframe(input_dataframe = id_4, outcome_variable = "serBilir", time_variable = "year", end_of_study_column = "years")

plt <- ggplot() +
  geom_point(data = test_id, aes(x = year, y = serBilir), color = "blue", alpha = .5, size = 4) +
  # geom_point(data = test_out, aes(x = year, y = serBilir), size = 4, alpha = .5) +
  # geom_line(data = test_out, aes(x = year, y = serBilir)) +
  scale_x_continuous(breaks = seq(0, max(x$whole_year), by = 1)) +
  theme_minimal()

plt

ggplot2::ggsave(plot = plt, filename = file.path(project_path, sprintf("plots/%s_new_algo.png", today())))

id_4 <- data %>% dplyr::filter(id == 4)
id_4$serBilir[1] <- NaN
id_4$whole_year

id_4 %>% dplyr::mutate(lag_year = lag(year), lead_year = lead(year), lag_y = lag(serBilir), lead_y = (serBilir))

id_4 %>%
  dplyr::group_by(id, years, whole_year, serBilir, year) %>%
  dplyr::mutate(lag = lag(serBilir), lead(serBilir)) %>%
  dplyr::summarise(first_in_year = min(year), last_in_year = max(year))

# function to get points to calculate AUC
get_time_interval <- function(input_dataframe, x, y, x_start, x_end) {
  output_data <- dplyr::select(input_dataframe[[x]], inputdataframe[[y]]) %>%
    dplyr::filter(between(input_dataframe[[x]]), x_start, x_end)

  return(output_data)
}

dplyr::select(id_4, year, serBilir, whole_year) %>%
  dplyr::group_by(whole_year) %>%
  dplyr::summarise(min_y = min(serBilir), max_y = max(serBilir), min_year = min(year), max_year = max(year)) %>%
  dplyr::mutate(
    floor_year = floor(min_year),
    ceil_year = ceiling(max_year)
  )

max(id_4$whole_year)
min(id_4$whole_year)


dplyr::filter(between(input_dataframe[[x]]), x_start, x_end)

get_time_interval(id_4, "year", "serBilir", 0, 2)
DescTools::AUC(c(0, 0.55, 1, 1.001, 1.02, 1.99), c(1.6, 1.6, 1.6, 1.7, 1.7, 3.2))

plot(c(0, 0.55, 1, 1.001, 1.02, 1.99), c(1.6, 1.6, 1.6, 1.7, 1.7, 3.2), ylim = c(0, 4))

# Testing for itterations
# Update year to whole_year if year is the maximum within each whole_year group
id_2 <- dplyr::filter(data, id == 2 & whole_year > 0)
id_2

# Get last year and last recorded y-outcome:
max_year <- floor(df$years[1])
last_y <- tail(df$serBilir, 1)

# update the last row of recorded data. Set the last year to whole_year
df <- id_2 %>%
  # group_by(whole_year) %>%
  mutate(test = ifelse(
    year == max(year),
    whole_year,
    year
  ))
# update all rows
df <- id_4 %>%
  group_by(whole_year) %>%
  mutate(test = ifelse(
    year == max(year),
    whole_year,
    year
  ))
df

my_function <- function(input_data, time_variable) {
  df <- input_data %>%
    group_by(whole_year) %>%
    mutate(!!sym(time_variable) := ifelse(
      !!sym(time_variable) == max(!!sym(time_variable)),
      whole_year,
      !!sym(time_variable)
    ))
  return(df)
}

x <- my_function(id_4, "year")

# If no data recorded at the start, create a row to hold the first recorded value.
if (df$whole_year[1] != 0) {
  new_row <- df[1, ]
  new_row$whole_year <- 0
  new_row$serBilir <- df$serBilir[1]
  new_row$year <- 0
  df <- rbind(new_row, df)
}

# Keep adding new rows until the end of follow up
while (max(df$whole_year) < max_year) {
  new_row <- df[nrow(df), ]
  new_row$whole_year <- max(df$whole_year) + 1
  new_row$year <- max(df$whole_year) + 1
  new_row$serBilir <- last_y
  df <- rbind(df, new_row)
}

df$auc[1:nrow(df)] <- sapply(1:nrow(df), function(i) {
  DescTools::AUC(
    x = df$year[1:i],
    y = df$serBilir[1:i],
    method = "trapezoid"
  )
})

df %>%
  group_by(id, whole_year) %>%
  summarise(auc = max(auc))


# Test cases --------------------------------------------------------------

input_test <- tidyr::tibble(
  id = c(1, 1, 1, 1, 2, 2, 2, 2),
  years = c(6, 6, 6, 6, 4.6, 4.6, 4.6, 4.6),
  year = c(0.5, 1.5, 2, 3.5, 1.5, 3, 3.5, 4.5),
  serBilir = c(1, 2, 2, 1, 1, 2, 1, 0.5)
)

create_auc <- function(input_dataframe, max_timespan = 10, y_outcome, recorded_time, end_of_follow_up) {
  input_dataframe <- input_dataframe %>%
    dplyr::mutate(
      whole_year = ceiling(year)
    )

  # Get unique IDs
  unique_ids <- unique(input_dataframe$id)

  output_data <- dplyr::tibble()

  output_tibble <- dplyr::tibble()

  for (id_unique in unique_ids) {
    dataframe_subset <- input_dataframe %>% dplyr::filter(id == id_unique)

    # If no data recorded at the start, create a row to hold the first recorded value.
    if (dataframe_subset$whole_year[1] > 0) {
      new_row <- dataframe_subset[1, ]
      new_row$whole_year <- 0
      new_row[[recorded_time]] <- 0
      dataframe_subset <- rbind(new_row, dataframe_subset)
    }

    # Create a copy of the last row with data
    last_row <- dataframe_subset %>% dplyr::slice(n())
    # Set the stop criterion:
    # The minimum of either the max_timespan (function input) or the floored value of end_of_follow_up column.
    stop_criterion <- min(max_timespan, floor(last_row[[end_of_follow_up]]))

    # Set the last recorded time
    latest_recording <- ceiling(last_row[[recorded_time]])

    # Calculate AUC on the fly
    interval <- 1

    id_vec <- vector()
    x_time_vec <- vector()
    y_outcome_vec <- vector()
    # for (year_interval in dataframe_subset$whole_year) {
    # while (interval <= max(dataframe_subset$whole_year)) {
    while (interval <= stop_criterion) {
      temp_auc <- dplyr::tibble()

      auc_data <- dataframe_subset %>%
        dplyr::filter(whole_year <= interval)
      last_auc_row <- auc_data %>% dplyr::slice(n())

      if (last_auc_row[[recorded_time]] < interval & last_auc_row[[recorded_time]] != 0) {
        last_auc_row[[recorded_time]] <- interval #ceiling(last_auc_row[[recorded_time]])
        temp_auc <- rbind(auc_data, last_auc_row)
      } else {
        temp_auc <- auc_data
      }

      auc <- DescTools::AUC(temp_auc[[recorded_time]], temp_auc[[y_outcome]])

      id_vec <- append(id_vec, id_unique)
      x_time_vec <- append(x_time_vec, sprintf("Interval 0:%s",interval))
      y_outcome_vec <- append(y_outcome_vec, auc)
      # id_vec <- append(id_vec, id_unique)
      # auc_vector <- rbind(auc_vector, c(id_unique, auc, max(temp_auc$whole_year)))
      # output_vector <- temp_auc %>%
      #   dplyr::select(all_of("id", "whole_year", "auc")) %>%
      #   dplyr::group_by("id") %>% summarise(whole_year = max(whole_year),
      #                                       auc = max(auc)
      #                                       ) %>% rbind(auc_vector, rbind)
      interval <- interval + 1
    }

    new_rows <- tibble(id = id_vec, x_time = x_time_vec, y_outcome = y_outcome_vec)
    # output_tibble <- dplyr::tibble(id = id_vec, interval_upper_bound = x_time_vec, y_outcome = y_outcome_vec)
    output_tibble <- bind_rows(output_tibble, new_rows)
    # add rows to the end, meaning: to either max_timespan or two the last of years!


    # Loop over the dataframe subset as long as the latest_recording <= stop criterion
    # Then add the new row to the dataframe subset
    while (latest_recording <= stop_criterion) {
      new_row <- last_row
      new_row$year <- latest_recording
      new_row$whole_year <- latest_recording
      # Update the subset
      dataframe_subset <- rbind(new_row, dataframe_subset)
      # Update latest_recording
      latest_recording <- latest_recording + 1
    }

    output_data <- rbind(dataframe_subset)
  }

  print(output_tibble)

  return(output_tibble)
}

# debugonce(create_auc)
test <- create_auc(input_test,
  y_outcome = "serBilir",
  recorded_time = "year",
  end_of_follow_up = "years",
  max_timespan = 5
)
test %>% slice(n())

# calculate AUC for T=1,2,3,4,5

output_test <- data.frame(
  id = c(1, 2),
  AUC_T1 = c(1, NA),
  AUC_T2 = c(3, 2),
  AUC_T3 = c(5, 3.75),
  AUC_T4 = c(5.75, 5),
  AUC_T5 = c(6.75, NA)
)

plot(input_test$year, input_test$serBilir, col = rep(1:2))

plt <- ggplot(data = input_test, aes(x = year, y = serBilir, colour = as.factor(id))) +
  geom_point(size = 5, alpha = 0.9) +
  geom_line() +
  ylim(0, NA) +
  xlim(0, NA) +
  scale_color_grey() +
  theme_minimal()

plt

ggsave(plt, filename = file.path(plots, sprintf("%s_test_cases.pdf", today())), device = "pdf")
