library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(data.table)
# library(purrr)


# ----------------------- Helper Functions to Implement ------------------------

#' Read the expression data "csv" file.
#'
#' Function to read microarray expression data stored in a csv file. The
#' function should return a sample x gene tibble, with an extra column named
#' "subject_id" that contains the geo accession ids for each subject.
#'
#' @param filename (str): the file to read.
#'
#' @return
#' @export
#'
#' @examples expr_mat <- read_expression_table('example_intensity_data.csv')
read_expression_table <- function(filename) {
  
  data <- readr::read_delim(filename, delim = ' ') %>% 
    t() %>%
    as_tibble(rownames="subject_id")
  data[1,1] <- 'subject_id'
  colnames(data) <- data[1,]
  data <- data[-1,]

  return (data)
}

#' Replaces all '.' in a string with '_'
#'
#' @param str String to operate upon.
#'
#' @return reformatted string.
#' @export
#'
#' @examples
#' period_to_underscore("foo.bar")
#' "foo_bar"
period_to_underscore <- function(str) {
  str <- str_replace_all(str, "[.]", "_")
  return (str)
}


# rename variables:
# Age_at_diagnosis to Age
# SixSubtypesClassification to Subtype
# normalizationcombatbatch to Batch

#' Rename and select specified columns.
#'
#' Function to rename Age_at_diagnosis, SixSubtypesClassification, and
#' normalizationcombatbatch columns to Age, Subtype, and Batch, respectively. A
#' subset of the data should be returned, containing only the Sex, Age, TNM_Stage,
#' Tumor_Location, geo_accession, KRAS_Mutation, Subtype, and Batch columns.
#'
#' @param data (tibble) metadata information for each sample
#'
#' @return (tibble) renamed and subsetted tibble
#' @export
#'
#' @examples rename_and_select(metadata)
#' 
#' 
rename_and_select <- function(data) {
  
  data
  data <- rename(data, Subtype = "SixSubtypesClassification", Age = "Age_at_diagnosis", Batch = "normalizationcombatbatch")
  data <- select(data, 'Sex', 'Age', 'TNM_Stage', 'Tumor_Location', 'geo_accession', 'KRAS_Mutation', 'Subtype', 'Batch')
  
  return (data)
}


#' Create new "Stage" column containing "stage " prefix.
#'
#' Creates a new column "Stage" with elements following a "stage x" format, where
#' `x` is the cancer stage data held in the existing TNM_Stage column. Stage
#' should have a factor data type.
#'
#' @param data  (tibble) metadata information for each sample
#'
#' @return (tibble) updated metadata with "Stage" column
#' @export
#'
#' @examples metadata <- stage_as_factor(metadata)
stage_as_factor <- function(data) {
  
  data <- mutate(data, Stage = paste0('stage ', TNM_Stage))
  data$Stage <- as.factor(data$Stage)
  
  return (data)
}


#' Calculate age of samples from a specified sex.
#'
#' @param data (tibble) metadata information for each sample
#' @param sex (str) which sex to calculate mean age. Possible values are "M"
#' and "F"
#'
#' @return (float) mean age of specified samples
#' @export
#'
#' @examples mean_age_by_sex(metadata, "F")
mean_age_by_sex <- function(data, sex) {
  
  data <- group_by(data, Sex) %>%
    summarise(Avg = mean(Age))
  
  val <- data$Avg[data$Sex==sex]
  
  return (val)
}


#' Calculate average age of samples within each cancer stage. Stages should be
#' from the newly created "Stage" column.
#'
#' @param data (tibble) metadata information for each sample
#'
#' @return (tibble) summarized tibble containing average age for all samples from
#' each stage.
#' @export
#'
#' @examples age_by_stage(data)
age_by_stage <- function(data) {
  
  data <- group_by(data, Stage) %>%
    summarise(mean(Age))
  
  return (data)
}

#' Create a cross tabulated table for Subtype and Stage using dplyr methods.
#'
#' @param data (tibble) metadata information for each sample
#'
#' @return (tibble) table where rows are the cancer stage of each sample, and the
#' columns are each cancer subtype. Elements represent the number of samples from
#' the corresponding stage and subtype. If no instances of a specific pair are
#' observed, a zero entry is expected.
#' @export
#'
#' @examples cross_tab <- dplyr_cross_tab(metadata)
subtype_stage_cross_tab <- function(data) {
  
  crosstab_data <- data %>%
    group_by(Stage, Subtype) %>%
    tally() %>%
    spread(Subtype, n) %>%
    replace_na(list(C3=0, C4=0))
  
  return (crosstab_data)
}

#' Summarize average expression and probe variability over expression matrix.
#'
#' @param exprs An (n x p) expression matrix, where n is the number of samples,
#' and p is the number of probes.
#'
#' @return A summarized tibble containing `main_exp`, `variance`, and `probe`
#' columns documenting average expression, probe variability, and probe ids,
#' respectively.
summarize_expression <- function(exprs) {
  
  exprs <- exprs %>% 
    mutate_at(c(2:54676), as.numeric)
  
  mean_probes <- select(exprs, where(is.numeric)) %>% colMeans()
  xvar_probes <- cbind(lapply(exprs[2:54676], FUN = var))
  mean_var_probe = tibble(mean_exp = mean_probes, variance = xvar_probes, probe = names(mean_probes))
  
  return (mean_var_probe)
}