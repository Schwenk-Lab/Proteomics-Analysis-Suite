# ------------------------------------------------------------------------------
# Script:        StudySpecFuncs.R
# Author:        Erik Kling Åhlén, eahle@kth.se
# Date:          2025-10-08
# Purpose:       Functions specific for the thesis and the associated dataset
# ------------------------------------------------------------------------------

################################################################################
# Create Ordinal Categorical variables from Continous variables
################################################################################
createCategoricalVars <- function(sinfo, sample_id = "sample_id") {
  
  # Create categorical columns for age, BMI, alcohol consumption and dense area
  sinfo <- sinfo %>%
    dplyr::mutate(
      alcoholCat = dplyr::case_when(
        alcohol_gram_week == 0 ~ "0",
        alcohol_gram_week > 0 & alcohol_gram_week <= 25 ~ ">0-≤25",
        alcohol_gram_week > 25 & alcohol_gram_week <= 50 ~ ">25-≤50",
        alcohol_gram_week > 50 & alcohol_gram_week <= 75 ~ ">50-≤75",
        alcohol_gram_week > 75 & alcohol_gram_week <= 100 ~ ">75-≤100",
        alcohol_gram_week > 100 & alcohol_gram_week <= 150 ~ ">100-≤150",
        alcohol_gram_week > 150 & alcohol_gram_week <= 250 ~ ">150-≤250",
        alcohol_gram_week > 250 ~ ">250"
      ),
      ageCat = dplyr::case_when(
        age < 50 ~ "≤40s",
        age >= 50 & age < 60 ~ "50s",
        age >= 60 & age < 70 ~ "60s",
        age >= 70 ~ "≥70s",
      ),
      denseCat = dplyr::case_when(
        stratus_densearea_cm2 < 14.12 ~ "<33rd",
        stratus_densearea_cm2 >= 14.12 & stratus_densearea_cm2 < 34.88 ~ "Mid",
        stratus_densearea_cm2 >= 34.88 ~ ">67th"
      ),
      bmiCat = dplyr::case_when(
        bmi >= 0 & bmi < 18.5 ~ "Underweight",
        bmi >= 18.5 & bmi < 25 ~ "Normal",
        bmi >= 25 & bmi < 30 ~ "Overweight",
        bmi >= 30 ~ "Obese"
      )
    ) %>%
    # Convert variables to ordered factors:
    dplyr::mutate(
      alcoholCat = factor(alcoholCat,
                          levels = c("0", ">0-≤25", ">25-≤50", ">50-≤75",
                                     ">75-≤100", ">100-≤150", ">150-≤250",
                                     ">250"), ordered = TRUE),
      ageCat = factor(ageCat,
                      levels = c("≤40s", "50s", "60s", "≥70s"), ordered = TRUE),
      bmiCat = factor(bmiCat,
                      levels = c("Underweight", "Normal", "Overweight",
                                 "Obese"), ordered = TRUE),
      denseCat = factor(denseCat,
                        levels = c("<33rd", "Mid", ">67th"), ordered = TRUE)
    ) %>%
    dplyr::mutate(
      alcoholCat = droplevels(alcoholCat),
      ageCat = droplevels(ageCat),
      bmiCat = droplevels(bmiCat),
      denseCat = droplevels(denseCat),
      plasmaAge = droplevels(plasmaAge),
    ) %>%
    as.data.frame()
  
  rownames(sinfo) <- sinfo[[sample_id]]
  assign("select.sinfo", select.sinfo, envir = .GlobalEnv)
  
  return(sinfo)
}
process_data <- function(sinfo, binfo, npx) {
  # Make sure sample_id and individual_id are characters
  sinfo$sample_id <- as.character(sinfo$sample_id)
  sinfo$individual_id <- as.character(sinfo$individual_id)
  npx$sample_id <- as.character(npx$sample_id)
  
  # Make x_BC numeric
  sinfo$x_BC <- as.numeric(as.character(sinfo$x_BC))
  
  # Remove failed samples
  sinfo <- sinfo[sinfo$sample_type != "fail", ]
  npx <- npx[npx$sample_id %in% sinfo$sample_id, ]
  
  # Create CancerTime column based on blood_draw_date and bc_invasive_1stdiagdate
  sinfo <- sinfo %>%
    mutate(
      # Convert the date columns from character to Date
      blood_date = as.Date(blood_draw_date),
      diag_date  = as.Date(bc_invasive_1stdiagdate),
      
      # Compute the difference in years (approximate)
      diff_years = as.numeric(diag_date - blood_date) / 365
    ) %>%
    mutate(
      CancerTime = case_when(
        x_BC == 0 ~ 0,  # No cancer samples get 0 regardless of dates
        
        # When blood_draw_date comes before the diagnosis date:
        blood_date < diag_date & diff_years > 2  ~ 2,  # more than 2 years difference
        blood_date < diag_date & diff_years <= 2 ~ 1,   # 2 years or less difference
        
        # When the diagnosis date comes before the blood draw date:
        diag_date < blood_date ~ 3,
        
        # Any other cases become NA (e.g., equal dates)
        TRUE ~ NA_real_
      )
    )
  
  # Create categorical variables for year, month and season
  sinfo$year <- as.numeric(format(sinfo$blood_date, "%Y"))
  sinfo$plasmaAge <- as.numeric(2020 - sinfo$year)
  sinfo$month <- as.numeric(format(sinfo$blood_date, "%m"))
  sinfo$season <- ifelse(sinfo$month %in% 3:5, "Spring",
                         ifelse(sinfo$month %in% 6:8, "Summer",
                                ifelse(sinfo$month %in% 9:11, "Fall", "Winter")))
  
  
  # Create categorical variables for alcohol, age, bmi     
  sinfo <- sinfo %>%
    mutate(
      alcoholCat = case_when(
        alcohol_gram_week == 0 ~ "0",
        alcohol_gram_week > 0  & alcohol_gram_week <= 25  ~ ">0-≤25",
        alcohol_gram_week > 25 & alcohol_gram_week <= 50  ~ ">25-≤50",
        alcohol_gram_week > 50 & alcohol_gram_week <= 75  ~ ">50-≤75",
        alcohol_gram_week > 75 & alcohol_gram_week <= 100 ~ ">75-≤100",
        alcohol_gram_week > 100 & alcohol_gram_week <= 150 ~ ">100-≤150",
        alcohol_gram_week > 150 & alcohol_gram_week <= 250 ~ ">150-≤250",
        alcohol_gram_week > 250 ~ ">250"
      ),
      ageCat = case_when(
        age < 50  ~ "<=40s",
        age >= 50 & age < 60  ~ "50s",
        age >= 60 & age < 70  ~ "60s",
        age >= 70 ~ "70s=<"
      ),
      bmiCat = case_when(
        bmi >= 0  & bmi < 18.5  ~ "Underweight",
        bmi >= 18.5 & bmi < 25   ~ "Normal",
        bmi >= 25 & bmi < 30  ~ "Overweight",
        bmi >= 30 ~ "Obese"
      ),
      denseCat = dplyr::case_when(
        stratus_densearea_cm2 < 14.12  ~ "<33rd",
        stratus_densearea_cm2 >= 14.12 & stratus_densearea_cm2 < 34.88  ~ "Mid",
        stratus_densearea_cm2 >= 34.88 ~ ">67th"
      )
      
    ) %>%
    # Convert to factors (with ordering where appropriate)
    mutate(
      year = factor(year, levels = c(2010, 2011, 2012, 2013)),
      plasmaAge = factor(plasmaAge, levels = c(7, 8, 9, 10), ordered = TRUE),
      month = factor(month,levels = c(1,2,3,4,5,6,7,8,9,10,11,12), 
                     labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", 
                                "Aug", "Sep", "Oct", "Nov", "Dec") ),
      season = factor(season, levels = c("Spring", "Summer","Fall", "Winter")),
      
      alcoholCat = factor(alcoholCat, levels = c("0", ">0-≤25", ">25-≤50", ">50-≤75", 
                                                 ">75-≤100", ">100-≤150", ">150-≤250", 
                                                 ">250"), ordered = TRUE),
      ageCat = factor(ageCat, levels = c("<=40s", "50s", "60s", "70s=<"), ordered = TRUE),
      bmiCat = factor(bmiCat, levels = c("Underweight", "Normal", "Overweight", 
                                         "Obese"), ordered = TRUE),
      denseCat = factor(denseCat, 
                        levels = c("<33rd","Mid",">67th"), ordered = TRUE),
      
      menopause_status = factor(menopause_status, levels = c(1, 2, 3), 
                                labels = c("pre-menopause", 
                                           "peri-menopause", 
                                           "post-menopause")),
      
      mht_status = factor(mht_status, levels = c(0, 1, 2),
                          labels = c("never", "previous", "current")),
      
      smoking_status = factor(smoking_status, levels = c(0, 1, 2),
                              labels = c("never","previous", "current")),
      
      x_BC = factor(x_BC, levels = c(0, 1), labels = c("control", "cancer")),
      
      CancerTime = factor(CancerTime, levels = c(0, 1, 2),
                          labels = c("control","<2years",">2years"))
    ) %>%
    as.data.frame()  # In case the result is a tibble and you prefer a classic dataframe
  
  
  # Set sample_id as rownames and remove sample_id column from npx
  rownames(sinfo) <- sinfo$sample_id
  rownames(npx) <- npx$sample_id
  npx <- npx[-1]  # Remove sample_id column from npx
  
  
  # Region selection
  regions <- c("stockholm", "skane", "all")
  cat("\nPlease select a region to process (1 = stockholm, 2 = skane, 3 = all, 0 = exit):\n")
  region_choice <- as.numeric(readline())
  
  if (region_choice == 0) {
    cat("\nNo region selected. Exiting script.\n")
    return(NULL) # Exit if no region is selected
  } else if (region_choice %in% 1:3) {
    selected_region <- regions[region_choice]
    region_filter <- if (selected_region == "stockholm") "stockholm" else if (selected_region == "skane") "skane" else NULL  # NULL means "all"
    
    # Notify region selection
    cat(paste("\nRegion selected:", selected_region, "\n"))
    
    # Filter sinfo and binfo by region if applicable
    filtered_sinfo <- if (!is.null(region_filter)) sinfo %>% filter(region == region_filter) else sinfo
    filtered_binfo <- if (!is.null(region_filter)) binfo %>% filter(region == region_filter) else binfo
    filtered_npx <- npx[rownames(npx) %in% rownames(filtered_sinfo), ]
    
    # Cancer selection
    cancer_options <- c("cancer", "no cancer", "all")
    cat("\nPlease select cancer status to process (1 = cancer, 2 = no cancer, 3 = all, 0 = exit):\n")
    cancer_choice <- as.numeric(readline())
    
    if (cancer_choice == 0) {
      cat("\nNo cancer status selected. Exiting script.\n")
      return(NULL)
    } else if (cancer_choice %in% 1:3) {
      selected_cancer <- cancer_options[cancer_choice]
      
      # Notify cancer selection
      cat(paste("\nCancer status selected:", selected_cancer, "\n"))
      
      # Filter sinfo and npx by cancer status
      filtered_sinfo <- if (selected_cancer == "cancer") {
        filtered_sinfo %>% filter(x_BC == 1)
      } else if (selected_cancer == "no cancer") {
        filtered_sinfo %>% filter(x_BC == 0)
      } else {
        filtered_sinfo
      }
      filtered_npx <- filtered_npx[rownames(filtered_npx) %in% rownames(filtered_sinfo), ]
      
      # Panel selection
      panels <- c("CAM", "IMONC", "all")
      cat("\nPlease select a panel to process (1 = CAM, 2 = IMONC, 3 = all, 0 = exit):\n")
      panel_choice <- as.numeric(readline())
      
      if (panel_choice == 0) {
        cat("\nNo panel selected. Exiting script.\n")
        return(NULL)
      } else if (panel_choice %in% 1:3) {
        selected_panel <- panels[panel_choice]
        
        # Notify panel selection
        cat(paste("\nPanel selected:", selected_panel, "\n"))
        
        # Filter binfo by panel (if not "all") and by region if applicable
        filtered_binfo <- if (selected_panel != "all") filtered_binfo %>% filter(panel == selected_panel) else filtered_binfo
        
        # Replace invalid characters in protein names and column names
        filtered_binfo$protein_name <- gsub("[ /]", "-", filtered_binfo$protein_name)
        colnames(filtered_npx) <- gsub("\\.", "-", colnames(filtered_npx))
        
        # For regions "stockholm" or "skane", set rownames of filtered_binfo to its protein_name column.
        if (selected_region %in% c("stockholm", "skane")) {
          rownames(filtered_binfo) <- filtered_binfo$protein_name
        }
        
        # Matching protein names between binfo and npx
        matching_protein_names <- intersect(filtered_binfo$protein_name, colnames(filtered_npx))
        
        # Identify any mismatches (optional)
        mismatched_protein_names <- setdiff(filtered_binfo$protein_name, colnames(filtered_npx))
        
        # Subset npx based on matching proteins
        if (length(matching_protein_names) > 0) {
          filtered_npx <- filtered_npx[, matching_protein_names, drop = FALSE]
        } else {
          cat("\nNo matching protein names found for npx. Subsetting skipped.\n")
        }
        
        # New prompt: remove samples with CancerTime == 2 if desired
        cat("\nDo you want to remove samples where the individual received the cancer diagnosis less than 2 years after the blood sample draw?\n")
        cat("(Enter 1 for Yes, 0 for No):\n")
        remove_choice <- as.numeric(readline())
        
        if (remove_choice == 1) {
          # Remove samples with CancerTime == 1 from filtered_sinfo
          filtered_sinfo <- filtered_sinfo[filtered_sinfo$CancerTime != 2, ]
          # Remove corresponding samples from filtered_npx based on rownames
          filtered_npx <- filtered_npx[rownames(filtered_npx) %in% rownames(filtered_sinfo), ]
          
          # Update the variables in the Global Environment
          assign(paste0(selected_region, if(selected_cancer != "all") paste0(".", selected_cancer) else "", 
                        if(selected_panel != "all") paste0(".", selected_panel) else "", ".sinfo"), 
                 filtered_sinfo, envir = .GlobalEnv)
          assign(paste0(selected_region, if(selected_cancer != "all") paste0(".", selected_cancer) else "", 
                        if(selected_panel != "all") paste0(".", selected_panel) else "", ".npx"), 
                 filtered_npx, envir = .GlobalEnv)
          
          cat("\nSamples with CancerTime == 1 have been removed.\n")
        } else {
          cat("\nNo samples were removed based on CancerTime.\n")
        }
        
        # Ask for a frequency cutoff for protein removal
        cat("\nEnter a cutoff value between 0 and 1 for frequency below LOD \n (proteins with freq_below_lod above this cutoff will be removed):\n")
        cutoff <- as.numeric(readline())
        if (is.na(cutoff) || cutoff < 0 || cutoff > 1) {
          cat("Invalid cutoff provided; using default value of 0.9.\n")
          cutoff <- 0.9
        }
        # Here selected_proteins are those whose freq_below_lod is below the cutoff
        selected_proteins <- filtered_binfo$protein_name[filtered_binfo$freq_below_lod > cutoff]
        cat("\n", length(selected_proteins), "proteins found with freq_below_lod above the cutoff.\n")
        if (length(selected_proteins) > 0) {
          filtered_binfo <- filtered_binfo[!(filtered_binfo$protein_name %in% selected_proteins), ]
          filtered_npx <- filtered_npx[, !(colnames(filtered_npx) %in% selected_proteins), drop = FALSE]
          cat("Proteins removed from binfo and npx based on the cutoff.\n")
        } else {
          cat("No proteins were removed based on the cutoff.\n")
        }
        
        
        # Variable naming: create a prefix using region, cancer, and panel
        variable_name_prefix <- paste0(
          selected_region, 
          if (selected_cancer != "all") paste0(".", selected_cancer) else "",
          if (selected_panel != "all") paste0(".", selected_panel) else ""
        )
        
        # Assign the filtered data to the global environment
        assign(paste0(variable_name_prefix, ".sinfo"), filtered_sinfo, envir = .GlobalEnv)
        assign(paste0(variable_name_prefix, ".binfo"), filtered_binfo, envir = .GlobalEnv)
        assign(paste0(variable_name_prefix, ".npx"), filtered_npx, envir = .GlobalEnv)
        
        cat("\nProcessing complete.\n")
        
      } else {
        cat("\nInvalid panel selection. Exiting script.\n")
      }
    } else {
      cat("\nInvalid cancer status selected. Exiting script.\n")
    }
  } else {
    cat("\nInvalid region selected. Exiting script.\n")
  }
}
