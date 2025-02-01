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

#AUC claculation
  #PATIENT FACTORS
    AGE <- data$AGE
    WT <- data$WT
    HT <- data$HT
    FEMALE <- data$FEMALE

  #PK and PD
    DOSE <- data$DOSE
    AUC <- data$AUC
    CMAX <- data$CMAX
    DV <- data$DV

  #Calculating BMI AND TO THE DATAFRAME
    data$BMI <- data$WT / ((data$HT/100)^2)

  #DV_delta calculation
    data <- data %>% 
      group_by(ID) %>% 
      mutate(DV_delta = DV[TIME == 240] - DV[TIME == 0])
    boxplot(DV_delta ~ DOSE, data = data, main = "Effect (DV change from start to end of study) by Dose", xlab = "Dose", ylab = "DV change")

    DV_delta <- data %>% 
      group_by(ID) %>% 
      mutate(DV_delta = DV[TIME == 240] - DV[TIME == 0])
    
    # Filter out rows with missing data in DV_change
    data <- data %>% drop_na(DV_delta)
    
    # DV_0 and DV_240
        Data_0 <- subset(data, TIME == 0, select = c(ID,AUC,DOSE,CMAX, AGE,WT,FEMALE, HT,DOSE,TAU, STUDY_ARM, side_effect, BMI, DV_delta))
        Data_240 <- subset(data, TIME == 240, select = c(ID,AUC,DOSE,CMAX, AGE,WT,FEMALE, HT,DOSE,TAU, STUDY_ARM, side_effect, BMI, DV_delta))

  # Drop rows with missing data in DV_change
     data <- data %>% drop_na(DV_delta)

  #inspect
    names(data)

  # Check the class of the converted variables
    sapply(data[c("AGE", "WT","DOSE", "HT", "AUC", "CMAX", "DV","BMI", "DV_delta")], class)



# Step 2: Exploratory Data Analysis
  # Summary statistics and scatter plots
    summary(data)
    pairs(~ DV_delta + CMAX + AUC + AGE + WT + DOSE, data = data)

  # Summary statistics of the key variables
    summary(data)



# Step 3: Linear Model
  #LINEAR REGRESSION
    #Patient Factors
    model_1 <- lm(HT ~ FEMALE, data = data)
    model_2 <- lm(WT ~ FEMALE, data = data)
    model_3 <- lm(WT ~ HT, data = data)
    model_4 <- lm(WT ~ AGE, data = data)
    model_5 <- lm(WT ~ HT+AGE, data = data)
    model_6 <- lm(WT ~ HT*AGE, data = data)
    model_7 <- lm(BMI ~ AGE, data = data)
    
    #PK parameters
    # Filter the data for specific DOSE values
    filtered_data <- data[data$DOSE %in% c(5, 10, 20, 150, 200), ]
    
    model_8 <- lm(CMAX ~ DOSE, data = data) #plot for extracting the MED
    
    model_9 <- lm(AUC ~ CMAX, data = data)
    model_10 <- lm(AUC ~ DOSE, data = data) #plot for extracting the MED
    
    #Patient factors with PK parameters (CMAX, AUC, DOSE)
    model_11 <- lm(AUC ~ BMI, data = data)
    model_12 <- lm(AUC ~ AGE, data = data)
    model_13 <- lm(DOSE ~ AGE, data = data)
    model_14 <- lm(CMAX ~ AGE, data = data)
    model_15 <- lm(AUC ~ FEMALE, data = data)
    
    #Patient factors and PD
    model_16 <- lm(DV_delta ~ BMI, data = data)
    model_17 <- lm(DV_delta ~ AGE, data = data)
    model_18 <- lm(DV_delta ~ BMI+AGE, data = data)
    model_19 <- lm(DV_delta ~ BMI*AGE, data = data)
    model_20 <- lm(DV_delta ~ FEMALE, data = data)
    model_21 <- lm(DV_delta ~ BMI+AGE+FEMALE, data = data)
    
    #PK and PD
    model_22 <- lm(DV_delta ~ AUC, data = filtered_data) #plot for extracting auc THAT ACHIEVES MED
    model_23 <- lm(DV_delta ~ CMAX, data = data) #plot
    model_24 <- lm(DV_delta ~ DOSE, data = filtered_data)
    model_25 <- lm(DV_delta ~ AUC+BMI+FEMALE, data = filtered_data)
    
    #filter out an outlier in PD
    data_filtered <- data[data$ID != 117, ]
    
    #check
    data_filtered
    
    model_26 <- lm(DV_delta ~ AUC, data = data_filtered)
    model_27 <- lm(DV_delta ~ CMAX, data = data_filtered)
    model_28 <- lm(DV_delta ~ DOSE, data = data_filtered)
    model_29 <- lm(DV_delta ~ AUC+BMI+FEMALE, data = data_filtered)
    
    # Perform the regression and get summary statistics
      summary(model_22) 
      summary(model_10) 
    
    # Example plot, replace with plot_X
    plot_22 <- ggplot(data_filtered, aes(x = AUC, y = DV_delta)) +
      geom_point() +  # Scatter plot of the original data points
      geom_smooth(method = "lm", color = "blue") 
    plot_22
    
    plot_10 <- ggplot(data, aes(x = DOSE, y = AUC)) +
      geom_point() +  # Scatter plot of the original data points
      geom_smooth(method = "lm", color = "blue") 
    
    # Plotting the Linear Regression
    plot_22 <- ggplot(data_filtered, aes(x = AUC, y = DV_delta)) +
      geom_point() +  # Scatter plot of the original data points
      geom_smooth(method = "lm", color = "blue", se = TRUE) +  # Show uncertainty with confidence intervals
      labs(title = "Linear Model: DV_delta ~ AUC",
           x = "AUC", y = "DV_delta") +
      theme_minimal()
    
    print(plot_22)
    
    #AUC estimation
      # model_22 
      # DV <- (-1.96 + 3.163 * (data_filtered$AUC))
    
    # AUC estimation using model coefficients 
      # coef(model_22)[1] = intercept, coef(model_22)[2] = AUC
    DV <- coef(model_22)[1] + coef(model_22)[2] * data_filtered$AUC
  
    #MED ESTIMATION
    model_10

    #model_22 <- lm(DV_delta ~ AUC, data = data)
    #model_10 <- lm(AUC ~ DOSE, data = data) #plot for extracting the MED
    
    
    # Arrange all plots in a grid layout (3 rows and 5 columns)
    combined_plots <- grid.arrange(plot_22, plot_10, 
                                   nrow=2)
    
    # Save the combined plots as a PNG file with high resolution
    ggsave("combined_regression_plots_22+10.png", combined_plots,
           width = 12, height = 20, dpi = 600)
     
    # NON-linear models
    nls_fit_2 <- nls(DV_delta ~ E0 + ((EMAX * (AUC^HILL))/(EC50^HILL + AUC^HILL)),  data=data_filtered, 
                     start = list( E0= -200,EMAX = 1100, 
                                   EC50 = 100,HILL=1), trace = TRUE)
    summary(nls_fit_2)
 
    #VERIFICATION
    AIC(model_22, nls_fit_2)
   
    
    # Create a function to generate plots
    create_plot <- function(model, x_var, y_var) {
      ggplot(f_data_f[!is.na(data[[x_var]]) & !is.na(f_data[[y_var]]), ], aes_string(x = x_var, y = y_var)) +
        geom_point(size = 3, shape = 21, fill = "lightblue", color = "black", stroke = 0.5) +
        geom_smooth(method = "lm", se = TRUE, fill = "lightblue", color = "blue") +
        labs(title = paste(y_var, "~", x_var),
             x = x_var,
             y = y_var) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.grid.major = element_line(size = 0.5),
          panel.grid.minor = element_blank(),
          axis.title.x = element_text(face = "bold"),
          axis.title.y = element_text(face = "bold")
        )
    }
    
    
    #PLOTTING ALL GRAPHS TOGETHER
    
    plot_22 + geom_smooth(method="lm",se = FALSE) +
      geom_line(aes(y=fitted(nls_fit_2)),col="red",lty="dotted",lwd=1)
    
    
    #Calculation of AUC
      # AUC <- ((200*(30.15^3.98))/(264.37-200))^(1/3.98)
      # AUC
    # AUC at MED = 40.08598
    
    #Extracting from AUC values
      # Dose <- 2.334*40.09
      # Dose
    #Dose = 93.57
    params <- coef(model_22)
    
    
    # Assign parameter values from the model output
    E0 <- params["E0"]
    EMAX <- params["EMAX"]
    EC50 <- params["EC50"]
    HILL <- params["HILL"]
    
    # Define DV_delta 
    DV_delta <- 200  
    
    # Compute AUC dynamically using extracted values
    AUC <- ((DV_delta * EC50^HILL) / (EMAX - DV_delta))^(1 / HILL)

    # Print AUC
    cat("Predicted AUC for MED:", AUC, "\n")
    
    
    #Calcule MED: predicted dose at the calculated AUC for delta_DV = 200
    model_10 <- lm(AUC ~ DOSE, data = data) #plot for extracting the MED
    model_10_summary <- summary(model_10)
    
    plot_10 <- ggplot(data, aes(x = DOSE, y = AUC)) +
      geom_point() +  # Scatter plot of the original data points
      geom_smooth(method = "lm", color = "blue")
    
    
    predicted_MED <- (AUC/(coef(model_10)[2]))  
    cat("Predicted dose for MED:", predicted_MED, "\n")
    
    #Dose = 154.3665 
   
    # Extract the necessary components
    estimate_dose <- coef(model_10)["DOSE"]
    std_error_dose <- summary(model_10)$coefficients["DOSE", "Std. Error"]
    df <- df.residual(model_10)  # Residual degrees of freedom
    t_critical <- qt(0.975, df)  # t-critical value for 95% CI
    
    # Calculate the confidence interval
    ci_lower <- estimate_dose - t_critical * std_error_dose
    ci_upper <- estimate_dose + t_critical * std_error_dose
    
    # Print the confidence interval
    cat("95% Confidence Interval for DOSE coefficient: [", ci_lower, ",", ci_upper, "]\n")
    

    upperdose <- AUC/ci_lower
    
    lowerdose <- AUC/ci_upper
    cat("Predicted dose for MED with 95% confidence interval: [", lowerdose, "-", upperdose, "]\n")

    
    lowMED <- predicted_MED - t_critical * std_error_dose
    
    highMED <- predicted_MED + t_critical * std_error_dose
    
    cat("Predicted dose for MED with 95% confidence interval: [", lowMED, "-", highMED, "]\n")
    
# Step 6: Investigate Side Effects
    # Analyze side effects rate for each dose group
    file_path <- "C:/Users/halaf/Downloads/Clinical data analysis_data_5.csv"
    data <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)
    side_effects <- data %>%
      group_by(STUDY_ARM) %>%
      summarize(SideEffectRate = mean(side_effect, na.rm = TRUE))
    
    # Compare predicted MED side effects to placebo
    placebo_rate <- side_effects$SideEffectRate[side_effects$STUDY_ARM == 1]
    med_side_effect_rate <- side_effects$SideEffectRate[which.min(abs(side_effects$STUDY_ARM - upperdose))]
    cat("Placebo side effect rate:", placebo_rate, "\n")
    cat("Predicted MED side effect rate:", med_side_effect_rate, "\n")
  
    
    # Step 7: Uncertainty Analysis Using Optimal Design Software
    
    # Install and load the PopED package
    if (!requireNamespace("PopED", quietly = TRUE)) install.packages("PopED")
    library(PopED)
    
    # Define the Emax model function
    ff_emax <- function(model_switch, xt, parameters, poped.db) {
      with(as.list(parameters), {
        DOSE <- xt
        y <- BASE + (EMAX * DOSE^(GAMMA)) / (ED50^(GAMMA) + DOSE^(GAMMA))
        return(list(y = y, poped.db = poped.db))
      })
    }
    
    # Define parameter generation function
    sfg_emax <- function(x, a, bpop, b, bocc) {
      parameters <- c(
        EMAX = bpop[1],  # Maximum effect
        ED50 = bpop[2],  # Dose at 50% of Emax
        GAMMA = bpop[3], # Hill coefficient
        BASE = bpop[4]   # Baseline effect
      )
      return(parameters)
    }
    
    # Create the PopED database with adjusted parameters
    poped_db <- create.poped.database(
      
      # The model
      ff_fun = ff_emax,
      fg_fun = sfg_emax,
      fError_fun = feps.add.prop,  # Additive and proportional error structure
      
      # Model parameters (fixed effects)
      bpop = c(
        EMAX = 545.0635  , # Hypothetical maximum effect
        ED50 = 30.14669  , # Dose for half-maximal effect
        GAMMA = 3.983012 , # Hill coefficient
        BASE = 1.875232    # Baseline effect
        
      ),
      
      # Residual unexplained variability
      sigma = c(
        PROP = 0.2^2,  # Proportional error variance
        ADD = 3^2      # Additive error variance
      ),
      
      # Design: Sample size and doses
      groupsize = 20,  # Individuals per dose group
      xt = rbind(
        group_1 = c(0),    # Placebo
        group_2 = c(5),   # Low dose
        group_3 = c(10),   # Low dose
        group_4 = c(20),  # Low dose
        group_4 = c(150),  # Medium dose
        group_6 = c(200)   # High dose
      ),
      
      # Timepoints for effect measurements
      ourzero = 0         # Baseline
    )
    
    # Debugging wrapper function
    evaluate_design_with_debugging <- function(poped_db) {
      tryCatch({
        # Evaluate the design
        design_eval <- evaluate_design(poped_db)
        cat("Design Evaluation Completed Successfully:\n")
        print(design_eval)
        
        # Extract and display FIM statistics
        fim_stats <- design_eval$poped.db$FIM_stats
        if (is.null(fim_stats)) stop("FIM statistics could not be computed.")
        
        cat("\nExpected Standard Errors (SE) of Model Parameters:\n")
        print(fim_stats)
        
        # Validate FIM
        if (any(is.na(design_eval$poped.db$FIM))) {
          cat("\nWarning: FIM contains NA values. Consider revising model or design.\n")
        }
      }, error = function(e) {
        cat("\nError in evaluating design or computing FIM:\n", e$message, "\n")
        cat("This could be due to poorly scaled parameters or invalid design inputs.\n")
      })
    }
    
    # Run the evaluation with debugging
    evaluate_design_with_debugging(poped_db)
    
    # Visualize model prediction and uncertainty
    library(ggplot2)
    plot_model_prediction(poped_db, PI = TRUE) + 
      xlab("Dose (mg)") + 
      ylab("Effect (DV)") + 
      ggtitle("Predicted Dose-Response Curve with Uncertainty")
    
    
    # DESIGN OPTIMIZATION
    # First we need to define the design space in our database:
    poped_db_2 <- create.poped.database(poped_db,minxt=0,maxxt=200)
    
    # Now we can optimize the doses of the design and plot the design:
    output <- poped_optim(poped_db_2, opt_xt=TRUE)
    
    summary(output)
    plot_model_prediction(output$poped.db, PI=T) + xlab("Dose") 
    
    # How many individuals would you need in the original design to match the 
    # information gain in the optimized design?
    optimize_n_eff(poped_db,output$ofv)
    
    
    
    
    
    
    
    
    
    


