# Clinical-Statistical-Data-Analysis
This team project aimed to determine the minimally effective dose (MED) of Drug X
# Ensure required libraries are installed and loaded
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("tidyr")) install.packages("tidyr")
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("nlme")) install.packages("nlme")
if (!require("MASS")) install.packages("MASS")
if (!require("car")) install.packages("car")
if (!require("PopED")) install.packages("PopED")
if (!require("openxlsx")) install.packages("openxlsx")
if (!require("broom")) install.packages("broom")
if (!require("nlstools")) install.packages("nlstools")
if (!require("gridExtra")) install.packages("gridExtra")
if (!require("FSA")) install.packages("FSA")

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(nlme)
library(MASS)
library(car)
library(PopED)
library(openxlsx)
library(broom)
library(nlstools)
library(gridExtra)
library(FSA)

# Upload the dataset and check data structures 
file_path <- "C:/Users/halaf/Downloads/Clinical data analysis_data_5.csv"
data <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)
# Data Structure
str(data)
# Check the column names in the data
colnames(data)

# Data preparation 
# Calculate change in effect (DV_change)
# Filter for baseline and end-of-study and rename columns
baseline <- data %>%
  filter(TIME == 0) %>%
  dplyr::select(ID, DV, STUDY_ARM) %>%
  rename(DV_baseline = DV)

end_of_study <- data %>%
  filter(TIME == 240) %>%
  dplyr::select(ID, DV) %>%
  rename(DV_end = DV)

# Merge baseline and end-of-study data
change_data <- merge(baseline, end_of_study, by = "ID")
change_data <- change_data %>%
  mutate(DV_change = DV_end - DV_baseline)

# Summarize changes by study arm
summary_stats <- change_data %>%
  group_by(STUDY_ARM) %>%
  summarise(
    mean_change = mean(DV_change, na.rm = TRUE),
    sd_change = sd(DV_change, na.rm = TRUE),
    n = n()
  )

print(summary_stats)

# Boxplot of DV Change by Study Arm
ggplot(change_data, aes(x = as.factor(STUDY_ARM), y = DV_change)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  labs(title = "Change in DV by Study Arm",
       x = "Study Arm",
       y = "Change in DV (ΔDV)") +
  theme_minimal()

# Assess Normality of Data
# Perform the Shapiro-Wilk test for normality
shapiro.test(change_data$DV_change)
# Visual Inspection
# Histogram
hist(change_data$DV_change, main = "Histogram of DV_change", xlab = "DV_change", col = "lightblue", border = "black")

# QQ Plot
qqnorm(change_data$DV_change)
qqline(change_data$DV_change, col = "red")

# Check Homogeneity of Variance
# Levene's Test
levene_result <- leveneTest(DV_change ~ as.factor(STUDY_ARM), data = change_data)
print(levene_result)

# Kruskal-Wallis Test to compare medians
kruskal_result <- kruskal.test(DV_change ~ as.factor(STUDY_ARM), data = change_data)
print(kruskal_result)

# Perform Dunn's Test with Holm's correction for all pairwise comparisons
dunn_result <- dunnTest(DV_change ~ as.factor(STUDY_ARM), data = change_data, method = "holm")

# Extract the results and filter for comparisons against Study Arm 1
dunn_result_holm <- dunn_result$res
dunn_result_holm_filtered <- dunn_result_holm[dunn_result_holm$Comparison %in% c("1 - 2", "1 - 3", "1 - 4", "1 - 5", "1 - 6"), ]

# Print the filtered results
print(dunn_result_holm_filtered)

# Function to determine MED with strict and relaxed criteria
determine_MED <- function(dunn_results, summary_stats, threshold = 200, p_value_threshold = 0.05) {
  # Merge Dunn's test results with summary statistics to get mean differences
  dunn_results <- dunn_results %>%
    mutate(
      Study_Arm = as.numeric(str_extract(Comparison, "\\d+$")),
      Mean_Difference = summary_stats$mean_change[match(Study_Arm, summary_stats$STUDY_ARM)] - summary_stats$mean_change[summary_stats$STUDY_ARM == 1]
    )
  
  # Step 1: Apply strict criteria
  med_candidates <- dunn_results %>%
    mutate(
      Significant = P.adj < p_value_threshold,
      Exceeds_Threshold = Mean_Difference > threshold
    ) %>%
    filter(Significant == TRUE & Exceeds_Threshold == TRUE) %>%
    arrange(Study_Arm)  # Sort by dose
  
  if (nrow(med_candidates) > 0) {
    # Return the smallest Study_Arm that meets criteria
    MED <- med_candidates$Study_Arm[1]
    cat("✅ The Minimally Effective Dose (MED) is STUDY_ARM:", MED, "\n")
    return(MED)
  } else {
    cat("⚠ No dose meets the strict MED criteria. Relaxing thresholds...\n")
    
    # Step 2: Relax criteria - Select highest Mean_Difference & lowest p-value
    relaxed_candidates <- dunn_results %>%
      arrange(desc(Mean_Difference), P.adj) %>%
      slice(1)  # Select the top candidate
    
    if (nrow(relaxed_candidates) > 0) {
      MED <- relaxed_candidates$Study_Arm
      cat("🔍 Using relaxed criteria, the closest MED is STUDY_ARM:", MED, 
          "with Mean_Difference =", relaxed_candidates$Mean_Difference, 
          "and p-value =", relaxed_candidates$P.adj, "\n")
      return(MED)
    } else {
      cat("❌ No suitable dose found under relaxed criteria. Returning NA.\n")
      return(NA)
    }
  }
}

# Call the function on Dunn's test results
MED <- determine_MED(dunn_result_holm_filtered, summary_stats)

# Output the final MED value
cat("Final MED Value:", MED, "\n")
