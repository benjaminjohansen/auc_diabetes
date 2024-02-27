# Setup -------------------------------------------------------------------

library("JMbayes2") # remove this when real data is available
library("dplyr")
library("DescTools")
library("tidyverse")
library("ggplot2")

project_path <- file.path("C:/SDCA/git/trap_auc")

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

  # Keep adding new rows until the end of follow up
  new_rows <- list()
  current_year <- max(output_dataframe$whole_year)
  while (current_year < max_year) {
    current_year <- current_year + 1
    new_row <- output_dataframe[nrow(output_dataframe), ]
    new_row$whole_year <- current_year
    new_row$x_time <- current_year
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
  if(max_timespan){
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

# debugonce(extend_dataframe)
id_2 <- data %>% dplyr::filter(id == 2)
id_4 <- data %>% dplyr::filter(id == 4)

x <- extend_dataframe(input_dataframe = id_2, outcome_variable = "serBilir", time_variable = "year", end_of_study_column = "years")

plt <- ggplot() +
  geom_point(data = id_2, aes(x = year, y = serBilir), color = "blue", alpha = .5, size = 4) +
  geom_point(data = x, aes(x = x_time, y = serBilir), size = 4, alpha = .5) +
  geom_line(data = x, aes(x = x_time, y = serBilir)) +
  scale_x_continuous(breaks = seq(0, max(x$whole_year), by = 1)) +
  theme_classic()

ggplot2::ggsave(plot = plt, filename = file.path(project_path, sprintf("plots/%s_point_moving.png", today())))

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
