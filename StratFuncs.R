# ------------------------------------------------------------------------------
# Script:        StratFuncs.R
# Author:        Erik Kling Åhlén, eahle@kth.se
# Date:          2025-10-08
# Purpose:       Functions for examining dominant traits and Stratification
# ------------------------------------------------------------------------------

################################################################################
# Stratification Menu Function:
################################################################################
stratification_menu <- function() {
  repeat {
    cat("\n=== Stratification Menu ===\n")
    cat("  1 = Create combined trait variables\n")
    cat("  2 = Mosaic plots\n")
    cat("  3 = Create tree of dominant trait combinations\n")
    cat("  4 = Stratification based on selected traits\n")
    cat("  5 = Drop unused factor levels\n")
    cat("  0 = Exit Stratification Menu\n")
    strat_choice <- read_choice("Enter your choice [0-5]: " )
    
    if (is.na(strat_choice) || strat_choice == 0) {
      cat("Exiting Stratification Menu.\n")
      break
    }
    else if (strat_choice == 1) {
      combine_traits_with_coding(sinfo = select.sinfo, reproduce = NULL)
    } 
    else if (strat_choice == 2) {
      interactive_mosaic_plot(data = select.sinfo)
    }
    else if (strat_choice == 3) {
      run_tree_builder(df = select.sinfo, param_tree_builder = NULL)
    }
    else if (strat_choice == 4) {
      create_stratus(data = select.sinfo, reproduce = NULL)
    }
    else if (strat_choice == 5) {
      clean_factor_columns(df = select.sinfo)
    }
  }
  }

################################################################################
# Create combined traits with optional coding of names:
################################################################################
combine_traits_with_coding <- function(sinfo, reproduce = NULL) {
  if (is.null(reproduce)) {
    
    reproduce <- data.frame(
      traits = I(list(NULL)),
      coded = NA,
      stringsAsFactors = FALSE
    )
    
    sel_input <- choose_traits(sinfo)
    sel_split <- unlist(strsplit(sel_input, split = ","))
    sel_split <- trimws(sel_split)
    reproduce$traits[[1]] <- sel_split
    
    # Determine if inputs are numeric or names
    if (all(suppressWarnings(!is.na(as.numeric(sel_split))))) {
      sel_indices <- as.numeric(sel_split)
      if (any(sel_indices < 1 | sel_indices > ncol(sinfo))) {
        stop("One or more columns are out of bounds.\n")
      }
      selected_df <- sinfo[, sel_indices, drop = FALSE]
    } else {
      if (!all(sel_split %in% colnames(sinfo))) {
        stop("Some provided columns do not exist in dataframe")
      }
      selected_df <- sinfo[, sel_split, drop = FALSE]
    }
  }
  # Ensure selected columns are factors
  selected_df[] <- lapply(selected_df, function(x) {
    if (!is.factor(x)) as.factor(x) else x
  })
  # Create combined factor using interaction
  combined_factor <- interaction(selected_df, sep = "+", drop = TRUE)
  
  # One-hot encoding of combined factor
  dummy_matrix <- model.matrix(~ combined_factor - 1,
                               data = data.frame(combined_factor = combined_factor))
  
  # Prompt if to use coded variable names
  code_choice <- readline(prompt = "Do you want to use coded variable names? (Y/N): ")
  code_choice <- tolower(trimws(code_choice))
  reproduce$coded <- code_choice
  
  if (code_choice %in% c("y", "yes")) {
    # Generate codes in order
    levs <- levels(combined_factor)
    num_codes <- length(levs)
    codes <- paste0("C", formatC(1:num_codes, width = 2, flag = "0"))
    
    # Create a dataframe containing the code mapping
    mapping <- data.frame(Code = codes, Combination = levs, stringsAsFactors = FALSE)
    
    # Replace dummy matrix colnames
    colnames(dummy_matrix) <- codes
    cat("Coded variable names selected\n")
  } else {
    cat("Using full combined names for dummy columns. \n")
    # Remove prefix:
    colnames(dummy_matrix) <- sub("^combined_factor", "", colnames(dummy_matrix))
  }
  # Remove Empty columns
  valid_cols <- which(colSums(dummy_matrix) > 0)
  if (length(valid_cols) < ncol(dummy_matrix)) {
    removed <- setdiff(seq_len(ncol(dummy_matrix)), valid_cols)
    cat("Removing dummy columns with no group members:",
        paste(colnames(dummy_matrix)[removed], collapse = ", "), "\n")
    dummy_matrix <- dummy_matrix[, valid_cols, drop = FALSE]
  }
  if (code_choice %in% c("y", "yes")) {
    mapping <- mapping[mapping$Code %in% colnames(dummy_matrix), , drop = FALSE]
    assign("code_mapping", mapping, envir = .GlobalEnv)
    cat("coded mapping assigned to global variable 'code_mapping'\n")
  }
  else {
    sel_split <- reproduce$traits[[1]]
    selected_df <- sinfo[, sel_split, drop = FALSE]
    selected_df[] <- lapply(selected_df, function(x) {
      if (!is.factor(x)) as.factor(x) else x
    })
    # Create combined factor using interaction
    combined_factor <- interaction(selected_df, sep = "+", drop = TRUE)
    
    # One-hot encoding of combined factor
    dummy_matrix <- model.matrix(~ combined_factor - 1,
                                 data = data.frame(combined_factor = combined_factor))
    
    # Prompt if to use coded variable names
    code_choice <- reproduce$coded
    code_choice <- tolower(trimws(code_choice))
    if (code_choice %in% c("y", "yes")) {
      # Generate codes in order:
      levs <- levels(combined_factor)
      num_codes <- length(levs)
      codes <- paste0("C", formatC(1:num_codes, width = 2, flag = "0"))
      
      # Create mapping df
      mapping <- data.frame(Code = codes, Combination = levs, stringsAsFactors = FALSE)
      
      # Replace dummy matrix colnames
      colnames(dummy_matrix) <- codes
    } else {
      # Optionally remove prefix:
      colnames(dummy_matrix) <- sub("^combined_factor", "", colnames(dummy_matrix))
    }
    # Remove Empty columns
    valid_cols <- which(colSums(dummy_matrix) > 0)
    if (length(valid_cols) < ncol(dummy_matrix)) {
      removed <- setdiff(seq_len(ncol(dummy_matrix)), valid_cols)
      dummy_matrix <- dummy_matrix[, valid_cols, drop = FALSE]
    }
    if (code_choice %in% c("y", "yes")) {
      mapping <- mapping[mapping$Code %in% colnames(dummy_matrix), , drop = FALSE]
      assign("code_mapping", mapping, envir = .GlobalEnv)
    }
  }
  result_df <- cbind(sinfo, dummy_matrix)
  assign("select.sinfo", result_df, envir = .GlobalEnv)
  assign("CombinedTraitsReproduce", reproduce, envir = .GlobalEnv)
  return(result_df)
}


################################################################################
# Stratification Function:
################################################################################

create_stratus <- function(data, reproduce = NULL) {
  # If data frame does not have row names, assign sequential numbers
  if (is.null(rownames(data))) {
    rownames(data) <- as.character(seq_len(nrow(data)))
  }
  if (is.null(reproduce)) {
    reproduce <- data.frame(
      Columns = NA,
      Values = NA,
      stringsAsFactors = FALSE
    )
    # Select column to use for stratification
    all_columns <- names(data)
    cat("Please select a column from the following list:\n")
    for (i in seq_along(all_columns)) {
      cat(i, ": ", all_columns[i], "\n")
    }
    column_choice <- as.numeric(readline(prompt = "Enter the number of your chosen column: "))
    if (is.na(column_choice) || column_choice < 1 || column_choice > length(all_columns)) {
      stop("Invalid selection for column")
    }
    selected_column <- all_columns[column_choice]
    reproduce$Columns <- selected_column
    
    # Extract unique values from selected column
    unique_vals <- as.character(unique(data[[selected_column]]))
    
    # Helper to extract numeric part following underscore
    extract_num <- function(x) {
      num <- sub(".*_([^_]+)$", "\\1", x)
      as.numeric(num)
    }
    
    numeric_values <- suppressWarnings(sapply(unique_vals, extract_num))
    if (all(is.na(numeric_values))) {
      ordered_values <- sort(unique_vals)
    } else {
      numeric_values[is.na(numeric_values)] <- Inf
      order_indices <- order(numeric_values)
      ordered_values <- unique_vals[order_indices]
    }
    
    # List ordered values and prompt to select one or more.
    cat("Values available in column", selected_column, ":\n")
    for (i in seq_along(ordered_values)) {
      cat(i, ": ", ordered_values[i], "\n")
    }
    cat("Enter the number(s) corresponding to your desired value(s) separated by comma:\n")
    value_input <- readline(prompt = "Your selection: ")
    
    # Split the input on commas remove spaces and convert to numeric
    value_choices <- as.numeric(unlist(strsplit(value_input, ",")))
    if (any(is.na(value_choices)) || length(value_choices) < 1 ||
        any(value_choices < 1) || any(value_choices > length(ordered_values))) {
      stop("Invalid selection for value(s).")
    }
    selected_values <- ordered_values[value_choices]
    reproduce$Values <- selected_values
    filtered_data <- data[data[[selected_column]] %in% selected_values, ]
    cat("Stratification based on: ", selected_values, "is complete\n")
  } else { 
    selected_column <- reproduce$Columns
    selected_values <- reproduce$Values
    filtered_data <- data[data[[selected_column]] %in% selected_values, ]
    }
  stratusClust <- rownames(filtered_data)
  bup.ptx <<- select.ptx
  bup.sinfo <<- select.sinfo
  select.ptx <<- select.ptx[stratusClust,]
  select.sinfo <<- select.sinfo[stratusClust,]
  
  if(exists("unscaled.ptx")){
    bup.unscaled <<- unscaled.ptx
    unscaled.ptx <<- unscaled.ptx[stratusClust,]
  }
  if(exists("unadjusted.ptx")) {
    bup.unadjusted <<- unadjusted.ptx
    unadjusted.ptx <<- unadjusted.ptx[stratusClust,]
  }
  
  if (selected_column == "region") {
    bup.binfo <<- select.binfo
    select.binfo <<- select.binfo[select.binfo$region == selected_values,]
  }
  
}

################################################################################
# Interactive mosaic plot function:
################################################################################

interactive_mosaic_plot <- function(data) {
  # List columns with numbered indices
  cat("Columns in the data:\n")
  colNames <- names(data)
  for (i in seq_along(colNames)) {
    cat(i,": ", colNames[i], "\n", sep = "")
  }
  cat("Enter the numbers (comma separated) of the categorical columns")
  cat_input <- readline()
  # Split input at commas, trim whitespace, convert to numeric
  cat_indices <- as.numeric(unlist(strsplit(cat_input, split = ",")))
  cat_vars <- colNames[cat_indices]
  
  # Prompt for number corresponding to the fixed variable of interest
  cat("Enter the number corresponding to the fixed variable: ")
  fixed_input <- readline()
  fixed_index <- as.numeric(fixed_input)
  fixed_var <- colNames[fixed_index]
  
  # Show choices
  cat("\nSelected categorical columns: ", paste(cat_vars, collapse = ", "), "\n")
  cat("Selected fixed variable: ", fixed_var, "\n\n")
  
  plot_cat_mosaic(data, cat_vars = cat_vars, fixed_var = fixed_var)
}

################################################################################
# Function for removing unused factor levels remaining after stratification:
################################################################################

clean_factor_columns <- function(df) {
  # Identify all factor columns
  factor_cols <- names(df)[sapply(df, is.factor)]
  
  # Drop unused levels from each factor column
  for (col in factor_cols) {
    df[[col]] <- droplevels(df[[col]])
  }
  # Assign cleaned data to global environment
  assign("select.sinfo", df, envir = .GlobalEnv)
}

################################################################################
# Bulid a hierarchical tree with branches based on dominant trait representation
################################################################################

build_exclusive_trait_tree <- function(data, trait_cols, filter_trait, filter_level,
                                       cluster_col = NULL, cluster_value = NULL,
                                       filter_cluster = FALSE, top_n = 2,
                                       min_prop = 0.05) {
  # Optional cluster filtering
  if (filter_cluster) {
    if (is.null(cluster_col) || is.null(cluster_value)) {
      stop("If filter_cluster = TRUE, you must provide cluster_col and cluster_value.")
    }
    data <- data %>% filter(.data[[cluster_col]] == cluster_value)
  }
  # Ensure traits are character
  data <- data %>% 
    mutate(across(all_of(trait_cols), as.character))
  
  # Root subset:
  root_subset <- data %>% filter(.data[[filter_trait]] == filter_level)
  
  # Recursive function for exclusive split:
  recurse_branch <- function(subset_df, used_traits, path_traits, box_root_df) {
    remaining_traits <- setdiff(trait_cols, used_traits)
    if (length(remaining_traits) == 0 || nrow(subset_df) == 0) return(NULL)
    
    # Count co-traits in subset
    counts <- subset_df %>%
      pivot_longer(cols = all_of(remaining_traits), names_to = "trait",
                   values_to = "level") %>%
      group_by(trait, level) %>%
      summarise(n = n(), .groups = "drop") %>%
      arrange(desc(n))
    
    # Keep only top_n with proportion ≥ min_prop
    box_root_size <- nrow(box_root_df)
    counts <- counts %>%
      mutate(prop_box = n / box_root_size)
    top_traits <- counts %>% 
      filter(prop_box >= min_prop) %>%
      slice_max(order_by = n, n = top_n)
    
    if (nrow(top_traits) == 0) return(NULL)
    
    # Exclusive branching
    remaining_df <- subset_df
    branches <- list()
    
    for (i in seq_len(nrow(top_traits))) {
      tr <- top_traits$trait[i]
      lv <- top_traits$level[i]
      branch_subset <- remaining_df[remaining_df[[tr]] == lv, , drop = FALSE]
      remaining_df <- remaining_df[remaining_df[[tr]] != lv, , drop = FALSE]
      
      # Percentage calculation with fixed denominator and cumulative numerator
      if (length(path_traits) == 0) {
        # First trait in box
        pct <- 100
        new_path_traits <- list(c(tr, lv))
      } else {
        # Cumulative numerator: Trait2 & Trait3 etc...
        cumulative_traits <- append(path_traits, list(c(tr, lv)))
        filtered <- box_root_df
        for (cond in cumulative_traits[-1]) {
          filtered <- filtered[filtered[[cond[1]]] == cond[2], , drop = FALSE]
        }
        pct <- (nrow(filtered) / nrow(box_root_df)) * 100
        new_path_traits <- cumulative_traits
      }
      label_lines <- c()
      for (j in seq_along(new_path_traits)) {
        cond <- new_path_traits[[j]]
        if (j == 1) {
          pct_line <- 100
        } else {
          filtered <- box_root_df
          for (k in 2:j) {
            cond2 <- new_path_traits[[k]]
            filtered <- filtered[filtered[[cond2[1]]] == cond2[2], , drop = FALSE]
          }
          pct_line <- (nrow(filtered) / nrow(box_root_df)) * 100
        }
        label_lines <- c(label_lines,
                         paste0(cond[1], "=", cond[2]),
                         paste0("(", round(pct_line, 1), "%)"))
      }
      branch_label <- paste(label_lines, collapse = "\n")
      branches[[length(branches) + 1]] <- list(
        label = branch_label,
        children = recurse_branch(branch_subset, c(used_traits, tr),
                                  new_path_traits, box_root_df)
      )
    }
    return(branches)
  }
  # Build tree starting from root
  tree <- list(
    label = paste0(filter_trait, "=", filter_level, "\n(100%)"),
    children = recurse_branch(root_subset, filter_trait, 
                              list(c(filter_trait, filter_level)),
                              root_subset)
  )
  
  return(tree)
}


################################################################################
# Helper function to create a tree using data.tree
################################################################################
list_to_tree <- function(x, root_name = "Root") {
  build_node <- function(node_list) {
    node <- Node$new(node_list$label)
    if (!is.null(node_list$children)) {
      for (child in node_list$children) {
        node$AddChildNode(build_node(child))
      }
    }
    return(node)
  }
  build_node(x)
}
################################################################################
# Helper functions for plotting tree results
################################################################################

tree_to_nodes <- function(node, depth = 0, parent = NA) {
  this <- tibble(label = node$label, depth = depth, parent = parent)
  if (is.null(node$children)) return(this)
  children <- map_dfr(node$children, ~ tree_to_nodes(.x, depth + 1, node$label))
  bind_rows(this, children)
}

tree_to_edges <- function(node) {
  if (is.null(node$children)) return(tibble(from = character(), to = character()))
  edges <- map_dfr(node$children, ~ tibble(from = node$label, to = .x$label))
  bind_rows(edges, map_dfr(node$children, tree_to_edges))
}

################################################################################
#  Assign positions for tree function with label formatting
################################################################################
assign_positions <- function(nodes_df, edges_df, wrap_width = 20, h_gap = 3, v_gap = 4) {
  # Helper: format a single label into "trait\n(percent)" pairs, joined by "\n"
  format_label <- function(lbl) {
    # Extract all chunks  without nested parentheses
    chunks <- stringr::str_extract_all(lbl, "[^()]+\\([^()]+\\)")[[1]]
    if (length(chunks) == 0) return(lbl)
    
    # For each chunk, make two lines: left part, then (percent)
    lines <- purrr::map_chr(chunks, function(ch) {
      left <- stringr::str_trim(sub("\\([^()]*\\)$", "", ch))            # text before (...)
      perc <- stringr::str_trim(sub(".*\\(([^()]*)\\)$", "(\\1)", ch))   # (...) as-is
      paste0(left, "\n", perc)
    })
    paste(lines, collapse = "\n")
  }
  
  # Helper: wrap each line separately to preserve hard breaks
  wrap_hardlines <- function(text, width) {
    lines <- unlist(strsplit(text, "\n", fixed = TRUE))
    wrapped_lines <- purrr::map_chr(lines, ~ stringr::str_wrap(.x, width = width))
    paste(wrapped_lines, collapse = "\n")
  }
  
  # Build label_wrapped robustly for each node
  nodes_df <- nodes_df %>%
    mutate(
      label_clean   = vapply(label, format_label, character(1)),
      label_wrapped = vapply(label_clean, wrap_hardlines, character(1), width = wrap_width),
      # Ensure it's plain character and not factor
      label_wrapped = as.character(label_wrapped)
    )
  
  # Tree topology helpers
  get_children <- function(parent_label) {
    edges_df %>% filter(from == parent_label) %>% pull(to)
  }
  
  positions <- list()
  x_cursor <- 0
  
  # Recursive placement with depth-dependent vertical spacing
  place_node <- function(label, depth) {
    children <- get_children(label)
    
    if (length(children) == 0) {
      wrapped <- nodes_df$label_wrapped[nodes_df$label == label]
      lines <- strsplit(wrapped, "\n", fixed = TRUE)[[1]]
      width <- max(nchar(lines)) * 0.12
      x_center <- x_cursor + width / 2
      positions[[label]] <<- c(x_center, -depth * (v_gap + depth * 0.5))
      x_cursor <<- x_cursor + width + h_gap
      return(width)
    } else {
      child_widths <- purrr::map_dbl(children, ~ place_node(.x, depth + 1))
      total_width <- sum(child_widths) + h_gap * (length(children) - 1)
      x_center <- mean(sapply(children, function(ch) positions[[ch]][1]))
      positions[[label]] <<- c(x_center, -depth * (v_gap + depth * 0.5))
      return(total_width)
    }
  }
  
  # Place roots
  roots <- setdiff(nodes_df$label, edges_df$to)
  purrr::map_dbl(roots, ~ place_node(.x, 0))
  
  # Build plotting data
  pos_df <- tibble::tibble(
    label = names(positions),
    x = sapply(positions, `[[`, 1),
    y = sapply(positions, `[[`, 2)
  )
  
  # Upward shift to use empty space above root
  y_min <- min(pos_df$y)
  pos_df$y <- pos_df$y - (y_min / 2)
  
  nodes_df <- nodes_df %>% left_join(pos_df, by = "label")
  edges_plot <- edges_df %>%
    left_join(nodes_df %>% select(from = label, x.from = x, y.from = y), by = "from") %>%
    left_join(nodes_df %>% select(to = label, x.to = x, y.to = y), by = "to")
  
  list(nodes = nodes_df, edges = edges_plot)
}


################################################################################
# Main plotting function for trait tree
################################################################################
plot_tree_overlap_free <- function(tree, wrap_width = 20, h_gap = 3, v_gap = 4,
                                   text_size = 2.8, padding = 0.15) {
  nodes_df <- tree_to_nodes(tree)
  edges_df <- tree_to_edges(tree)
  
  layout <- assign_positions(nodes_df, edges_df, wrap_width, h_gap, v_gap)
  
  # Ensure label_wrapped is plain character with \n intact
  layout$nodes$label_wrapped <- as.character(layout$nodes$label_wrapped)
  
  # Create Plot
  ggplot() +
    geom_segment(data = layout$edges,
                 aes(x = x.from, y = y.from, xend = x.to, yend = y.to),
                 arrow = arrow(length = unit(3, "mm")),
                 lineend = "round") +
    geom_label(data = layout$nodes,
               aes(x = x, y = y, label = label_wrapped),
               size = text_size,
               label.padding = unit(padding, "lines"),
               label.size = 0.15,
               lineheight = 1.0) +
    theme_void() +
    coord_cartesian(clip = "off") +
    theme(plot.margin = ggplot2::margin(80, 80, 80, 80))
}

################################################################################
# Runner function for tree builder
################################################################################

run_tree_builder <- function(df, param_tree_builder = NULL) {
  
  if (!is.null(param_tree_builder)) {
    message("Using pre-specified values — skipping interactive prompts.")
    
    traits <- param_tree_builder$traits
    filter_trait <- param_tree_builder$filter_trait
    filter_level <- param_tree_builder$filter_level
    filter_cluster <- param_tree_builder$filter_cluster
    cluster_col <- param_tree_builder$cluster_col
    cluster_value <- param_tree_builder$cluster_value
    top_n <- param_tree_builder$top_n
    min_prop <- param_tree_builder$min_prop
    plot_choice <- param_tree_builder$plot_choice
    
  } else {
    # Interactive mode
    traits <- choose_traits(df)
    filter_info <- choose_filter(df)
    cluster_info <- choose_cluster_filter(df)
    
    filter_trait <- filter_info$filter_trait
    filter_level <- filter_info$filter_level
    filter_cluster <- cluster_info$filter_cluster
    cluster_col <- cluster_info$cluster_col
    cluster_value <- cluster_info$cluster_value
    
    top_n <- as.integer(readline(prompt = "Enter top_n (number of branches from every branch): "))
    if (is.na(top_n) || top_n < 1) stop("Invalid top_n value.")
    
    min_prop <- as.numeric(readline(prompt = "Enter min_prop (minimum proportion, e.g., 0.05 for 5%): "))
    if (is.na(min_prop) || min_prop <= 0 || min_prop > 1) stop("Invalid min_prop value.")
    
    cat("\nChoose plotting method:\n")
    cat("[1] ggplot\n")
    cat("[2] DiagrammeR\n")
    plot_choice <- as.integer(readline(prompt = "Enter choice: "))
    
    # Prompt if to save selections
    save_choice <- readline(prompt = "Save these selections for reuse? (y/n): ")
    if (tolower(save_choice) == "y") {
      param_tree_builder <- list(
        traits = traits,
        filter_trait = filter_trait,
        filter_level = filter_level,
        filter_cluster = filter_cluster,
        cluster_col = cluster_col,
        cluster_value = cluster_value,
        top_n = top_n,
        min_prop = min_prop,
        plot_choice = plot_choice
      )
      assign("param_tree_builder", param_tree_builder, envir = .GlobalEnv)
      message("Selections saved in variable 'param_tree_builder'.")
    }
  }
  
  # Build the tree plot
  tree <- build_exclusive_trait_tree(
    data = df,
    trait_cols = traits,
    filter_trait = filter_trait,
    filter_level = filter_level,
    cluster_col = cluster_col,
    cluster_value = cluster_value,
    filter_cluster = filter_cluster,
    top_n = top_n,
    min_prop = min_prop
  )
  
  # Plot Results
  if (plot_choice == 1) {
    if (!exists("plot_tree_overlap_free")) {
      stop("Function plot_tree_overlap_free() not found in environment.")
    }
    p_tree <<- plot_tree_overlap_free(tree)
    print(p_tree)
    
  } else if (plot_choice == 2) {
    root_node <- list_to_tree(tree)
    SetGraphStyle(root_node, rankdir = "TB")
    g <- DiagrammeR::render_graph(ToDiagrammeRGraph(root_node))
    print(g)
    
  } else {
    stop("Invalid plotting choice.")
  }
  
  invisible(tree)
}
