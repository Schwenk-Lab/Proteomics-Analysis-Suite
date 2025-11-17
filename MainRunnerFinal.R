# ------------------------------------------------------------------------------
# Script:        MainRunnerFinal.R
# Author:        Erik Kling Åhlén, eahle@kth.se
# Date:          2025-10-08
# Purpose:       Main Menu runner function for the Interactive Analysis Suite
# ------------------------------------------------------------------------------


################################################################################
# Main Runner Function
################################################################################

run_main_menu <- function(sinfo, binfo, ptx) {
  # Log start time
  start_time <- Sys.time()
  
  cat("\n=====================================================\n")
  cat("Welcome to the Plasma Proteomics Analysis Suite\n")
  cat("=====================================================\n")
  
  # Main menu loop:
  repeat { 
    cat("\n=== Main Menu: ===\n")
    cat("  1 = Extract Data (specific for the study)\n")
    cat("  2 = Set Global Seed\n")
    cat("  3 = QC Workflow\n")
    cat("  4 = Stratification\n")
    cat("  5 = Preprocessing\n")
    cat("  6 = Clustering Parameter Optimization\n")
    cat("  7 = Clustering and Analysis\n")
    cat("  0 = Quit and show run summary\n")
    
    choice <- my_readline("Enter your choice [0-7]: ")
    
    if (choice == "1"){
      process_data(sinfo, binfo, ptx)
    } else if (choice == "2"){
      set_global_seed()
    } else if (choice == "3"){
      QC_main_menu()
    } else if (choice == "4"){
      stratification_menu()
    } else if (choice == "5"){
      preprocess_main_menu()
    } else if (choice == "6"){
      clustering_parameter_optimization_menu()
    } else if (choice == "7"){
      analysis_main_menu()
    } else if (choice == "0"){
      cat("Exiting main menu...\n")
      break
    } else {
      cat("Invalid option. Please choose a number from 0 to 7")
    }
  }
  # Record end time and calculate total run duration
  end_time <- Sys.time()
  total_time <- difftime(end_time, start_time, units = "secs")
  extracted_log <- extract_user_inputs()
  cat("\n=============================\n")
  cat("Run Summary:\n")
  cat("Total run time (seconds):", total_time, "\n")
  cat("Captured user inputs stored in extracted_log:\n")
  
  # Return a list with results if needed for later use
  invisible(list(total_runtime = total_time, inputs = user_inputs))
}

