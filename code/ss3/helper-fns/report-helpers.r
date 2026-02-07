extract_mfcl_quantities = function(rep_file, par_file, model_id, output_dir, fy = 1979, quarterly = FALSE) {
        # Read MFCL report and par files
        rep = read.MFCLRep(rep_file)
        par = read.MFCLPar(par_file,first.yr=fy)

        # Extract years as numeric from the model
        years = as.numeric(dimnames(popN(rep))$year)
        start_yr = min(years)
        end_yr = max(years)
        
        # Extract SSB with proper handling of MFCL output structure
        ssb_values = as.data.table(ssb(rep))
        
        # Aggregate SSB by year (summing or averaging across seasons based on quarterly flag)
        # For annual SSB, select "all" ages and aggregate seasons
        if(quarterly) {
            # If quarterly model, take average of quarters
            annual_ssb = ssb_values[age == "all", .(
                Value = mean(value)
            ), by = .(yr = as.numeric(year), unit, area, iter)]
        } else {
            # If not quarterly (or annual), sum across seasons
            annual_ssb = ssb_values[age == "all", .(
                Value = sum(value)
            ), by = .(yr = as.numeric(year), unit, area, iter)]
        }

        annual_ssb = annual_ssb[, .(yr, Value)]

        # Add placeholder standard deviations (adjust as needed)
        annual_ssb[, StdDev := Value * 0.001]  # Using 0.001 CV as placeholder
        
        # Format SSB data.table
        ssb_dt = annual_ssb[, .(
            id = rep(model_id, .N),
            type = "ssb",
            yr = yr,
            Value = Value,
            StdDev = StdDev
        )]
        
        # Extract dynamic depletion components
        # Get SSB with fishing (already extracted above)
        ssb_fishing = annual_ssb[, .(Yr = yr, SSB = Value)]
        
        # Get unfished SSB (dynamic B0)
        # Using adultBiomass_nofish or equivalent function
        ssb_nofishing_values = as.data.table(adultBiomass_nofish(rep))
        
        # Handle similar multi-dimensional structure for sbf0 output
        if(quarterly) {
            # If quarterly model, take average of quarters
            annual_sbf0 = ssb_nofishing_values[age == "all", .(
                SSB_nofishing = mean(value)
            ), by = .(Yr = as.numeric(year), unit, area, iter)]
        } else {
            # If not quarterly, sum across seasons
            annual_sbf0 = ssb_nofishing_values[age == "all", .(
                SSB_nofishing = sum(value)
            ), by = .(Yr = as.numeric(year), unit, area, iter)]
        }

        annual_sbf0 = annual_sbf0[, .(Yr, SSB_nofishing)]
        
        # Merge and format depletion data.table
        dep_dt = merge(ssb_fishing, annual_sbf0, by = "Yr")
        dep_dt[, `:=`(
            id = rep(model_id, .N),
            Era = ifelse(Yr <= end_yr, "TIME", "FORE")
        )]
        
        # Select and order columns to match SS output format
        dep_dt = dep_dt[, .(id, Yr, Era, SSB, SSB_nofishing)]

        # Get population numbers
        pop_values = as.data.table(popN(rep))
        
        if(quarterly) {
            # If quarterly model, take average of quarters
            annual_numbers = pop_values[, .(
                Value = mean(value)/1000
            ), by = .(yr = as.numeric(year), unit, area, iter)]
        } else {
            # If not quarterly, sum across seasons
            annual_numbers = pop_values[, .(
                Value = sum(value)/1000
            ), by = .(yr = as.numeric(year), unit, area, iter)]
        }
        
        annual_numbers = annual_numbers[, .(yr, Value)]

        # Add placeholder standard deviations
        annual_numbers[, StdDev := Value * 0.001]
        
        # Format numbers data.table
        numbers_dt = annual_numbers[, .(
            id = rep(model_id, .N),
            type = "numbers",
            yr = yr,
            Value = Value,
            StdDev = StdDev
        )]
        
        # Extract recruitment from youngest age class in population numbers
        pop_values = as.data.table(popN(rep))
        
        # Get the youngest age class
        # MFCL popN typically has ages as numeric values in the 'age' column
        youngest_age = min(as.numeric(pop_values$age[pop_values$age != "all"]))
        rec_values = pop_values[age == as.character(youngest_age)]
        
        # Aggregate by year based on quarterly setting
        if(quarterly) {
            annual_rec = rec_values[, .(
                Value = mean(value)/1000  # Convert to thousands of fish
            ), by = .(yr = as.numeric(year), unit, area, iter)]
        } else {
            annual_rec = rec_values[, .(
                Value = sum(value)/1000  # Convert to thousands of fish
            ), by = .(yr = as.numeric(year), unit, area, iter)]
        }
        
        annual_rec = annual_rec[, .(yr, Value)]
        
        # Use slightly higher CV for recruitment as it's typically more variable
        annual_rec[, StdDev := Value * 0.001]  # Using 0.1 CV as placeholder
        
        # Format recruitment data.table
        rec_dt = annual_rec[, .(
            id = rep(model_id, .N),
            type = "rec",
            yr = yr,
            Value = Value,
            StdDev = StdDev
        )]
        
        # Write output files
        fwrite(ssb_dt, file = file.path(output_dir, "ssb.csv"))
        fwrite(numbers_dt, file = file.path(output_dir, "N.csv"))
        fwrite(dep_dt, file = file.path(output_dir, "dyn_dep.csv"))
        fwrite(rec_dt, file = file.path(output_dir, "rec.csv"))
        
        return(list(ssb = ssb_dt, dep = dep_dt, N = numbers_dt, rec = rec_dt))
    }

    get_population_numbers = function(model_dir, model_name, beginning_of_year = FALSE, exclude_age0 = FALSE) {
        # Extract structured population data from Stock Synthesis output
        model = r4ss::SS_output(dir = model_dir, verbose = FALSE)
        
        # Retrieve age-structured abundance matrix
        natage = model$natage
        
        # Specify temporal reference point for population estimates
        time_point = ifelse(beginning_of_year, "B", "M")
        natage_filtered = natage[natage$"Beg/Mid" == time_point & !(natage$Era %in% c("VIRG", "FORE")),  ]
        
        # Identify age columns based on column names rather than positions
        age_cols = grep("^[0-9]+$", names(natage))
        
        # Implement age-0 exclusion if specified
        if(exclude_age0) {
            age_cols = age_cols[-1]  # Remove first age column (age-0)
        }
        
        # Extract temporal sequence for population trajectory
        years = unique(natage_filtered$Yr)
        total_numbers = numeric(length(years))
        
        # Compute total abundance for each temporal stratum
        for(i in seq_along(years)) {
            year_data = natage_filtered[natage_filtered$Yr == years[i], ]
            # Sum across all ages and any morphs/sexes present
            age_values = as.matrix(year_data[, age_cols])[4,] # hard code for last season
            total_numbers[i] = sum(age_values, na.rm = TRUE)
        }
        
        # Formulate standardized output data structure
        results = data.table::data.table(
            id = rep(model_name, length(years)),  # User-specified model identifier
            type = rep("numbers", length(years)),  # Type specification as "numbers"
            yr = years,
            Value = total_numbers,
            StdDev = rep(0.001, length(years))  # Placeholder uncertainty quantification
        )

        # Write output files
            data.table::fwrite(results, file = file.path(model_dir, "N.csv"))
        
        return(list(N=results))
    }

# Function to create comparative time series plots
plot_model_comparison_ts = function(model_ids, model_stem, categories = c("SSB", "Depletion (D)", "N", "Recruitment"), 
                                  show_se = TRUE, n_col = 2, model_labels = NULL, custom_colors = NULL,
                                  save_plot = FALSE, save_dir = NULL, plot_name = "model_comparison", 
                                  width = 10, height = 8, dpi = 300, file_type = "png") {
    
    # Validate inputs
    if(length(model_ids) < 1 || length(categories) < 1) {
      stop("Please provide at least one model ID and category")
    }

    # Handle model_stem as vector or single path
    if(length(model_stem) == 1) {
      model_stem = rep(model_stem, length(model_ids))
    } else if(length(model_stem) != length(model_ids)) {
      stop("If model_stem is a vector, it must have the same length as model_ids")
    }
    
    # If model_labels is not provided, use model_ids as labels
    if(is.null(model_labels)) {
      model_labels = model_ids
    } else if(length(model_labels) != length(model_ids)) {
      stop("If model_labels is provided, it must have the same length as model_ids")
    }
    
    # Create vector of model specifications with labels
    selected_models = data.frame(
      id = model_ids,
      path = model_stem,
      label = model_labels,
      stringsAsFactors = FALSE
    )
    
    # Initialize plot data list
    plot_dt.list = list("SSB" = NA, "Depletion (D)" = NA, "N" = NA, "Recruitment" = NA)
    
    # Process SSB if selected
    if("SSB" %in% categories) {
      plot_dt.list[["SSB"]] = rbindlist(lapply(1:nrow(selected_models), function(i) {
        dt = fread(file.path(selected_models$path[i], selected_models$id[i], "ssb.csv"))
        dt[, id := selected_models$id[i]]
        dt[, label := selected_models$label[i]]
        return(dt)
      })) %>%
        setnames(., c("Value", "StdDev"), c("value", "sd")) %>%
        .[, cv := sd/value] %>%
        .[, l95 := exp(log(value) - 2*sqrt(log(cv^2 + 1)))] %>%
        .[, u95 := exp(log(value) + 2*sqrt(log(cv^2 + 1)))] %>%
        .[, .(id, label, type, yr, value, l95, u95)] %>%
        .[, type := "SSB"] %>%
        .[, type := factor(type, levels = c("SSB", "Depletion (D)", "N", "Recruitment"))]
    }

    # Process N if selected
    if("N" %in% categories) {
      plot_dt.list[["N"]] = rbindlist(lapply(1:nrow(selected_models), function(i) {
        dt = fread(file.path(selected_models$path[i], selected_models$id[i], "N.csv"))
        dt[, id := selected_models$id[i]]
        dt[, label := selected_models$label[i]]
        return(dt)
      })) %>%
        setnames(., c("Value", "StdDev"), c("value", "sd")) %>%
        .[, cv := sd/value] %>%
        .[, l95 := exp(log(value) - 2*sqrt(log(cv^2 + 1)))] %>%
        .[, u95 := exp(log(value) + 2*sqrt(log(cv^2 + 1)))] %>%
        .[, .(id, label, type, yr, value, l95, u95)] %>%
        .[, type := "N"] %>%
        .[, type := factor(type, levels = c("SSB", "Depletion (D)", "N", "Recruitment"))]
    }
    
    # Process Depletion if selected
    if("Depletion (D)" %in% categories) {
      plot_dt.list[["Depletion (D)"]] = rbindlist(lapply(1:nrow(selected_models), function(i) {
        dt = fread(file.path(selected_models$path[i], selected_models$id[i], "dyn_dep.csv"))
        dt[, id := selected_models$id[i]]
        dt[, label := selected_models$label[i]]
        return(dt)
      })) %>%
        .[Era %in% c("TIME")] %>%
        setnames(., c("Yr"), c("yr")) %>%
        .[, value := SSB/SSB_nofishing] %>%
        .[, l95 := value] %>% # Note: No uncertainty quantified in depletion from dyn_dep.csv
        .[, u95 := value] %>%
        .[, type := "Depletion (D)"] %>%
        .[, .(id, label, type, yr, value, l95, u95)] %>%
        .[, type := factor(type, levels = c("SSB", "Depletion (D)", "N", "Recruitment"))]
    }
    
    # Process Recruitment if selected - MODIFIED TO SHOW RELATIVE RECRUITMENT
    if("Recruitment" %in% categories) {
      plot_dt.list[["Recruitment"]] = rbindlist(lapply(1:nrow(selected_models), function(i) {
        dt = fread(file.path(selected_models$path[i], selected_models$id[i], "rec.csv"))
        dt[, id := selected_models$id[i]]
        dt[, label := selected_models$label[i]]
        
        # Calculate mean recruitment for this model to rescale
        mean_rec = mean(dt$Value, na.rm = TRUE)
        
        # Rescale recruitment to relative values (mean = 1)
        dt[, Value := Value / mean_rec]
        dt[, StdDev := StdDev / mean_rec]  # Scale standard deviation proportionally
        
        return(dt)
      })) %>%
        setnames(., c("Value", "StdDev"), c("value", "sd")) %>%
        .[, cv := sd/value] %>%
        .[, l95 := exp(log(value) - 2*sqrt(log(cv^2 + 1)))] %>%
        .[, u95 := exp(log(value) + 2*sqrt(log(cv^2 + 1)))] %>%
        .[, .(id, label, type, yr, value, l95, u95)] %>%
        .[, type := "Recruitment"] %>%
        .[, type := factor(type, levels = c("SSB", "Depletion (D)", "N", "Recruitment"))]
    }
    
    # Combine all available data
    plot_dt = rbindlist(plot_dt.list[which(sapply(plot_dt.list, length) > 1)]) %>% 
      .[order(id, type, yr)] %>%
      .[,label:=factor(label,levels=model_labels)]
    
    # Validate that we have data to plot
    if(nrow(plot_dt) == 0) {
      stop("No data available for the selected models and categories")
    }
    
    # Create plot
    p = plot_dt %>%
      ggplot() +
      ylim(0, NA) +
      xlab("Year") +
      facet_wrap(~type, scales = "free_y", ncol = n_col)
    
    # Modify y-axis label for Recruitment to indicate relative values
    p = p + ylab("Metric") +
      geom_hline(yintercept = 0)
    
    # Add horizontal line at y=1 for Recruitment panel to indicate mean level
    if("Recruitment" %in% categories){
      p = p + geom_hline(data = data.frame(type = "Recruitment", y = 1), 
                        aes(yintercept = y), 
                        linetype = "dashed", 
                        color = "darkgray")
    }

    # Add uncertainty ribbons if requested
    if(show_se) {
      p = p + geom_ribbon(aes(x = yr, ymin = l95, ymax = u95, group = id, fill = label), alpha = 0.25)
    }
    
    # Add lines for mean values
    p = p + geom_path(aes(x = yr, y = value, group = id, color = label), linewidth = 1.5)
    p = p + geom_path(data=plot_dt[label==tail(model_labels,n=1)], aes(x = yr, y = value, group = id, color = label), linewidth = 1.5)
    
    # Apply color scales - either custom colors or viridis palette
    if(!is.null(custom_colors)) {
      # Validate that enough colors are provided
      if(length(custom_colors) < length(unique(plot_dt$label))) {
        warning("Not enough custom colors provided. Falling back to viridis palette.")
        p = p + viridis::scale_color_viridis("Model", begin = 0.1, end = 0.8, 
                                            direction = 1, option = "H", discrete = TRUE) +
          viridis::scale_fill_viridis("Model", begin = 0.1, end = 0.8, 
                                    direction = 1, option = "H", discrete = TRUE)
      } else {
        # Use custom colors
        p = p + scale_color_manual("Model", values = custom_colors) +
          scale_fill_manual("Model", values = custom_colors)
      }
    } else {
      # Use default viridis palette
      p = p + viridis::scale_color_viridis("Model", begin = 0.1, end = 0.8, 
                                          direction = 1, option = "H", discrete = TRUE) +
        viridis::scale_fill_viridis("Model", begin = 0.1, end = 0.8, 
                                  direction = 1, option = "H", discrete = TRUE)
    }
    
    # Apply consistent theme elements
    p = p + theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
            panel.grid.major = element_line(color = 'gray70', linetype = "dotted"),
            panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
            strip.background = element_rect(fill = "white"),
            legend.key = element_rect(fill = "white"))
    
    # Save the plot if requested
    if(save_plot) {
      # Validate save directory
      if(is.null(save_dir)) {
        save_dir = getwd()
        warning(paste("No save directory specified. Saving to current working directory:", save_dir))
      } else {
        # Create directory if it doesn't exist
        if(!dir.exists(save_dir)) {
          dir.create(save_dir, recursive = TRUE)
          message(paste("Created directory:", save_dir))
        }
      }
      
      # Construct full file path with extension
      file_path = file.path(save_dir, paste0(plot_name, ".", file_type))
      
      # Save using appropriate function based on file type
      if(file_type == "pdf") {
        pdf(file_path, width = width, height = height)
        print(p)
        dev.off()
      } else {
        ggsave(filename = file_path, plot = p, width = width, height = height, 
              dpi = dpi, units = "in")
      }
      
      message(paste("Plot saved to:", file_path))
      return(invisible(NULL))
    } else {
      return(p)
    }
  }

# Function to create comparative length-based selectivity plots
plot_model_comparison_selex = function(model_ids, model_stem, 
                                      n_col = 2, model_labels = NULL, custom_colors = NULL,
                                      save_plot = FALSE, save_dir = NULL, plot_name = "selex_comparison", 
                                      width = 10, height = 8, dpi = 300, file_type = "png") {
  
  # Validate inputs
  if(length(model_ids) < 1) {
    stop("Please provide at least one model ID")
  }

  # Handle model_stem as vector or single path
  if(length(model_stem) == 1) {
    model_stem = rep(model_stem, length(model_ids))
  } else if(length(model_stem) != length(model_ids)) {
    stop("If model_stem is a vector, it must have the same length as model_ids")
  }
  
  # If model_labels is not provided, use model_ids as labels
  if(is.null(model_labels)) {
    model_labels = model_ids
  } else if(length(model_labels) != length(model_ids)) {
    stop("If model_labels is provided, it must have the same length as model_ids")
  }
  
  # Create vector of model specifications with labels
  selected_models = data.frame(
    id = model_ids,
    path = model_stem,
    label = model_labels,
    stringsAsFactors = FALSE
  )
  
  # Read and process selectivity data from each model
  selex_list = lapply(1:nrow(selected_models), function(i) {
    # Read selectivity data
    selex_file = file.path(selected_models$path[i], selected_models$id[i], "selex_l.csv")
    if (!file.exists(selex_file)) {
      warning(paste("File not found:", selex_file))
      return(NULL)
    }
    
    dt = fread(selex_file)
    dt[, model_id := selected_models$id[i]]
    dt[, model_label := selected_models$label[i]]
    return(dt)
  })
  
  # Combine all selectivity data
  selex_data = rbindlist(selex_list, use.names = TRUE, fill = TRUE)
  
  # Check if we have data
  if(nrow(selex_data) == 0) {
    stop("No selectivity data found for the specified models")
  }
  
  # Read fleet summary data for each model to get fleet names
  fleet_list = lapply(1:nrow(selected_models), function(i) {
    # Read fleet summary data
    fleet_file = file.path(selected_models$path[i], selected_models$id[i], "fleet_summary.csv")
    if (!file.exists(fleet_file)) {
      warning(paste("Fleet summary file not found:", fleet_file))
      return(NULL)
    }
    
    dt = fread(fleet_file)
    dt[, model_id := selected_models$id[i]]
    return(dt)
  })
  
  # Combine all fleet data
  fleet_data = rbindlist(fleet_list, use.names = TRUE, fill = TRUE)
  
  # Check if we have fleet data
  if(nrow(fleet_data) == 0) {
    warning("No fleet summary data found, using original Fleet_name")
  } else {
    # Map fleet names from fleet_summary to selex_data
    selex_data = merge(
      selex_data,
      fleet_data[, .(fleet, fleetname, model_id)],
      by.x = c("Fleet", "model_id"),
      by.y = c("fleet", "model_id"),
      all.x = TRUE
    )
    
    # Use fleetname from fleet_summary if available, otherwise keep original Fleet_name
    selex_data[!is.na(fleetname), Fleet_name := fleetname]
  }
  
  # Check if fleet names are consistent across models
  fleet_consistency = selex_data[, .(
    unique_names = length(unique(Fleet_name))
  ), by = .(Fleet)]
  
  inconsistent_fleets = fleet_consistency[unique_names > 1, Fleet]
  if(length(inconsistent_fleets) > 0) {
    warning(paste("Inconsistent fleet names found for fleet(s):", 
                 paste(inconsistent_fleets, collapse = ", ")))
  }
  
  # Add sex label for better readability
  selex_data[Sex == 1, sex_label := "Female"]
  selex_data[Sex == 2, sex_label := "Male"]
  selex_data[is.na(sex_label), sex_label := "NA"]
  
  # For each Fleet, if there's only one sex, change the label to "Aggregated"
  fleet_sex_counts <- selex_data[, .(sex_count = uniqueN(Sex)), by = .(Fleet)]
  single_sex_fleets <- fleet_sex_counts[sex_count == 1, Fleet]
  
  if(length(single_sex_fleets) > 0) {
    selex_data[Fleet %in% single_sex_fleets, sex_label := "Aggregated"]
  }
  
  # Convert labels to factors to ensure consistent ordering
  selex_data[, model_label := factor(model_label, levels = model_labels)]
  # Sort unique sex labels alphabetically for factor levels
  sex_levels <- sort(unique(selex_data$sex_label))
  selex_data[, sex_label := factor(sex_label, levels = sex_levels)]
  
  # Create plot
  p = ggplot(selex_data, aes(x = variable, y = value, color = model_label, linetype = sex_label)) +
    geom_line(linewidth = 1.2) +
    facet_wrap(~ Fleet_name, scales = "free_x", ncol = n_col) +
    ylim(0, 1) +
    xlab("Length") +
    ylab("Selectivity") +
    labs(color = "Model", linetype = "Sex") +
    theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
          panel.grid.major = element_line(color = 'gray70', linetype = "dotted"),
          panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
          strip.background = element_rect(fill = "white"),
          legend.key = element_rect(fill = "white"))
  
  # Apply color scales - either custom colors or viridis palette
  if(!is.null(custom_colors)) {
    # Validate that enough colors are provided
    if(length(custom_colors) < length(unique(selex_data$model_label))) {
      warning("Not enough custom colors provided. Falling back to viridis palette.")
      p = p + viridis::scale_color_viridis("Model", begin = 0.1, end = 0.8, 
                                         direction = 1, option = "H", discrete = TRUE)
    } else {
      # Use custom colors
      p = p + scale_color_manual("Model", values = custom_colors)
    }
  } else {
    # Use default viridis palette
    p = p + viridis::scale_color_viridis("Model", begin = 0.1, end = 0.8, 
                                       direction = 1, option = "H", discrete = TRUE)
  }
  
  # Save the plot if requested
  if(save_plot) {
    # Validate save directory
    if(is.null(save_dir)) {
      save_dir = getwd()
      warning(paste("No save directory specified. Saving to current working directory:", save_dir))
    } else {
      # Create directory if it doesn't exist
      if(!dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE)
        message(paste("Created directory:", save_dir))
      }
    }
    
    # Construct full file path with extension
    file_path = file.path(save_dir, paste0(plot_name, ".", file_type))
    
    # Save using appropriate function based on file type
    if(file_type == "pdf") {
      pdf(file_path, width = width, height = height)
      print(p)
      dev.off()
    } else {
      ggsave(filename = file_path, plot = p, width = width, height = height, 
            dpi = dpi, units = "in")
    }
    
    message(paste("Plot saved to:", file_path))
    return(invisible(p))
  } else {
    return(p)
  }
}


# Function to create comparative CPUE index plots across multiple models
plot_index_comparison = function(model_ids, model_stem, 
                               show_se = TRUE, 
                               show_fit = TRUE, 
                               use_log_scale = FALSE,
                               apply_varadj = TRUE,
                               filter_lambda = TRUE,
                               n_col = 2, 
                               model_labels = NULL, 
                               custom_colors = NULL,
                               save_plot = FALSE, 
                               save_dir = NULL, 
                               plot_name = "cpue_comparison", 
                               width = 10, 
                               height = 8, 
                               dpi = 300, 
                               file_type = "png") {
  
  # Required packages
  required_packages <- c("data.table", "ggplot2", "viridis")
  for(pkg in required_packages) {
    if(!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is needed for this function to work. Please install it."))
    }
  }
  
  # Validate inputs
  if(length(model_ids) < 1) {
    stop("Please provide at least one model ID")
  }
  
  # Handle model_stem as vector or single path
  if(length(model_stem) == 1) {
    model_stem = rep(model_stem, length(model_ids))
  } else if(length(model_stem) != length(model_ids)) {
    stop("If model_stem is a vector, it must have the same length as model_ids")
  }
  
  # If model_labels is not provided, use model_ids as labels
  if(is.null(model_labels)) {
    model_labels = model_ids
  } else if(length(model_labels) != length(model_ids)) {
    stop("If model_labels is provided, it must have the same length as model_ids")
  }
  
  # Create vector of model specifications with labels
  selected_models = data.frame(
    id = model_ids,
    path = model_stem,
    label = model_labels,
    stringsAsFactors = FALSE
  )
  
  # Check if cpue.csv exists for all models
  file_exists = sapply(1:nrow(selected_models), function(i) {
    file.exists(file.path(selected_models$path[i], selected_models$id[i], "cpue.csv"))
  })
  
  if(sum(file_exists) == 0) {
    stop("The file cpue.csv does not exist for any of the models you selected.")
  } else if(sum(file_exists) < length(file_exists)) {
    warning("The file cpue.csv does not exist for some models. These models will be excluded.")
    selected_models = selected_models[file_exists, ]
  }
  
  # Read and compile CPUE data
  cpue_data_list = lapply(1:nrow(selected_models), function(i) {
    # Read CPUE data
    cpue_file = file.path(selected_models$path[i], selected_models$id[i], "cpue.csv")
    cpue_dt = fread(cpue_file)
    
    # Rename columns to be consistent with shiny app code
    setnames(cpue_dt, 
             c("Fleet", "Fleet_name", "Time", "Obs", "Exp", "SE"), 
             c("fleet", "fleetname", "time", "obs", "exp", "se"),
             skip_absent = TRUE)
    
    # Add model identifiers
    cpue_dt[, id := selected_models$id[i]]
    cpue_dt[, label := selected_models$label[i]]
    
    # Round observed values for consistent ID generation
    cpue_dt[, obs := round(obs, digits=3)]
    
    # Generate ID for grouping by model and fleet combination
    cpue_dt[, id3 := as.numeric(as.factor(paste0(id, "_", fleetname)))]
    
    return(cpue_dt)
  })
  
  # Combine all CPUE data
  plot_dt = rbindlist(cpue_data_list)
  
  # Read variance adjustment and lambda values
  varlambda_dt_list = lapply(1:nrow(selected_models), function(i) {
    model_id = selected_models$id[i]
    model_path = selected_models$path[i]
    
    # Read variance adjustments
    varadj_file = file.path(model_path, selected_models$id[i], "varadj.csv")
    if(file.exists(varadj_file)) {
      tmp_var = fread(varadj_file)
      if(nrow(tmp_var) == 0) {
        # Create empty table if file exists but is empty
        tmp_var = data.table(
          fleet = unique(plot_dt[id == model_id]$fleet),
          id = model_id,
          factor = 1,
          value = 0
        )[, .(id, factor, fleet, value)]
      } else {
        # Filter for factor 1
        tmp_var = tmp_var[factor == 1]
        
        # Add missing fleets with zero value
        tmp_var2 = data.table(
          fleet = unique(plot_dt[id == model_id]$fleet),
          id = model_id,
          factor = 1,
          value = 0
        )[, .(id, factor, fleet, value)]
        
        tmp_var2 = tmp_var2[!(fleet %in% tmp_var$fleet)]
        tmp_var = rbind(tmp_var, tmp_var2)[order(fleet)]
      }
    } else {
      # Create default table if file doesn't exist
      tmp_var = data.table(
        fleet = unique(plot_dt[id == model_id]$fleet),
        id = model_id,
        factor = 1,
        value = 0
      )[, .(id, factor, fleet, value)]
    }
    
    # Rename columns
    tmp_var = tmp_var[, .(id, fleet, value)]
    setnames(tmp_var, "value", "varadj")
    
    # Read lambdas
    lambda_file = file.path(model_path, selected_models$id[i], "lambdas.csv")
    if(file.exists(lambda_file)) {
      tmp_lambda = fread(lambda_file)
      if("like_comp" %in% names(tmp_lambda)) {
        tmp_lambda = tmp_lambda[like_comp == 1][, .(id, fleet, value)]
      } else {
        # Create default table if file doesn't have expected structure
        tmp_lambda = data.table(
          fleet = unique(plot_dt[id == model_id]$fleet),
          id = model_id,
          value = 0
        )[, .(id, fleet, value)]
      }
      
      # Add missing fleets with zero value
      tmp_lambda2 = data.table(
        fleet = unique(plot_dt[id == model_id]$fleet),
        id = model_id,
        value = 0
      )[, .(id, fleet, value)]
      
      tmp_lambda2 = tmp_lambda2[!(fleet %in% tmp_lambda$fleet)]
      tmp_lambda = rbind(tmp_lambda, tmp_lambda2)[order(fleet)]
    } else {
      # Create default table if file doesn't exist
      tmp_lambda = data.table(
        fleet = unique(plot_dt[id == model_id]$fleet),
        id = model_id,
        value = 0
      )[, .(id, fleet, value)]
    }
    
    # Rename columns
    setnames(tmp_lambda, "value", "lambda")
    
    # Merge variance adjustments and lambdas
    return(merge(tmp_var, tmp_lambda, by = c("id", "fleet")))
  })
  
  # Combine all variance adjustment and lambda data
  varlambda_dt = rbindlist(varlambda_dt_list)
  
  # Merge with CPUE data
  plot_dt = merge(plot_dt, varlambda_dt, by = c("id", "fleet"))
  
  # Apply variance adjustment if requested
  if(!apply_varadj) {
    plot_dt[, se := se - varadj]
  }
  
  # Check consistency of observed values across models for the same fleet/index
  consistency_check = plot_dt[, .(
    obs_values = list(unique(obs)),
    se_values = list(unique(se)),
    models = list(unique(label))
  ), by = .(fleetname, time)]
  
  # Check if observed values are consistent
  obs_inconsistent = consistency_check[sapply(obs_values, length) > 1]
  if(nrow(obs_inconsistent) > 0) {
    error_msg = paste0(
      "Observed CPUE values differ across models for the same fleet and time point:\n",
      paste(
        apply(obs_inconsistent[, .(fleetname, time)], 1, function(row) {
          return(paste0("  - Fleet: ", row[1], ", Time: ", row[2]))
        }),
        collapse = "\n"
      )
    )
    stop(error_msg)
  }
  
  # Check if standard errors are consistent, but only enforce if show_se is TRUE
  se_inconsistent = consistency_check[sapply(se_values, length) > 1]
  if(nrow(se_inconsistent) > 0) {
    inconsistency_msg = paste0(
      "Standard errors (post variance adjustment) differ across models for the same fleet and time point:\n",
      paste(
        apply(se_inconsistent[, .(fleetname, time)], 1, function(row) {
          return(paste0("  - Fleet: ", row[1], ", Time: ", row[2]))
        }),
        collapse = "\n"
      )
    )
    
    if(show_se) {
      # If we're showing SE bars, this is an error
      stop(inconsistency_msg)
    } else {
      # If we're not showing SE bars, just warn
      warning(inconsistency_msg)
    }
  }
  
  # Create a version of the data with just the first model's observations
  # This ensures we only plot the observations once
  obs_dt = plot_dt[, .SD[1], by = .(fleetname, time)]

  # Make sure obs_dt has the same calculated columns as plot_dt
  if(use_log_scale) {
    obs_dt[, lse := obs - se]
    obs_dt[, use := obs + se]
  } else {
    obs_dt[, lse := exp(log(obs) - se)]
    obs_dt[, use := exp(log(obs) + se)]
  }
  
  # Apply log transformation if requested
  if(use_log_scale) {
    plot_dt[, obs := log(obs)]
    plot_dt[, exp := log(exp)]
    plot_dt[, lse := obs - se]
    plot_dt[, use := obs + se]
    ylab_txt = "Index (log-scale)"
    yint = 0
  } else {
    plot_dt[, lse := exp(log(obs) - se)]
    plot_dt[, use := exp(log(obs) + se)]
    ylab_txt = "Index"
    yint = 1
  }
  
  # Filter by lambda if requested
  if(filter_lambda) {
    plot_dt = plot_dt[lambda == 1]
  }
  
  # Check if we have data to plot
  if(nrow(plot_dt) == 0) {
    stop("No data available after applying filters")
  }
  
  # Set label as factor with correct order
  plot_dt[, label := factor(label, levels = model_labels)]
  
  # Create base plot
  p = ggplot(plot_dt[order(id3, time)])  

  # Set y-axis limits for non-log scale
  if(!use_log_scale) {
    p = p + ylim(0, NA)
  }
  
  # Add common elements
  p = p + 
    xlab("Time") +
    ylab(ylab_txt) +
    geom_hline(yintercept = yint, linetype = "dashed", linewidth = 1) +
    facet_wrap(~fleetname, ncol = n_col)
  
  # Add error bars if requested
  if(show_se) {
    p = p + geom_errorbar(data = obs_dt, aes(x = time, ymin = lse, ymax = use), 
                       width = 0.2, linewidth = 0.7, color = "gray60")
  }
  
  # Add observed points
  p = p + geom_point(data = obs_dt, aes(x = time, y = obs), 
                   shape = 21, size = 3, stroke = 0.5, fill = "gray90")
  
  # Add fitted lines if requested
  if(show_fit) {
    p = p + geom_line(aes(x = time, y = exp, color = label), linewidth = 1)
  }
  
  # Apply color scales - either custom colors or viridis palette
  color_count = length(unique(plot_dt$label))
  color_legend = "Model"
  
  if(!is.null(custom_colors)) {
    # Validate that enough colors are provided
    if(length(custom_colors) < color_count) {
      warning("Not enough custom colors provided. Falling back to viridis palette.")
      p = p + 
        viridis::scale_color_viridis(color_legend, begin = 0.1, end = 0.8, 
                                   direction = 1, option = "H", discrete = TRUE) +
        viridis::scale_fill_viridis(color_legend, begin = 0.1, end = 0.8, 
                                  direction = 1, option = "H", discrete = TRUE)
    } else {
      # Use custom colors
      p = p + 
        scale_color_manual(color_legend, values = custom_colors) +
        scale_fill_manual(color_legend, values = custom_colors)
    }
  } else {
    # Use default viridis palette
    p = p + 
      viridis::scale_color_viridis(color_legend, begin = 0.1, end = 0.8, 
                                 direction = 1, option = "H", discrete = TRUE) +
      viridis::scale_fill_viridis(color_legend, begin = 0.1, end = 0.8, 
                                direction = 1, option = "H", discrete = TRUE)
  }
  
  # Apply consistent theme elements
  p = p + theme(
    panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
    panel.grid.major = element_line(color = 'gray70', linetype = "dotted"),
    panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
    strip.background = element_rect(fill = "white"),
    legend.key = element_rect(fill = "white")
  )
  
  # Save the plot if requested
  if(save_plot) {
    # Validate save directory
    if(is.null(save_dir)) {
      save_dir = getwd()
      warning(paste("No save directory specified. Saving to current working directory:", save_dir))
    } else {
      # Create directory if it doesn't exist
      if(!dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE)
        message(paste("Created directory:", save_dir))
      }
    }
    
    # Construct full file path with extension
    file_path = file.path(save_dir, paste0(plot_name, ".", file_type))
    
    # Save using appropriate function based on file type
    if(file_type == "pdf") {
      pdf(file_path, width = width, height = height)
      print(p)
      dev.off()
    } else {
      ggsave(filename = file_path, plot = p, width = width, height = height, 
             dpi = dpi, units = "in")
    }
    
    message(paste("Plot saved to:", file_path))
    invisible(p)
  } else {
    return(p)
  }
}



# Function to create comparative composition plots (length or weight) across multiple models
plot_composition_comparison = function(model_ids, model_stem, 
                                    comp_type = "length", # Can be "length" or "weight"
                                    show_fit = TRUE, 
                                    n_col = 2, 
                                    free_y_scale = TRUE,  # Option to have free y scales
                                    model_labels = NULL, 
                                    custom_colors = NULL,
                                    save_plot = FALSE, 
                                    save_dir = NULL, 
                                    plot_name = NULL, 
                                    width = 10, 
                                    height = 8, 
                                    dpi = 300, 
                                    file_type = "png") {
  
  # Required packages
  required_packages <- c("data.table", "ggplot2", "viridis")
  for(pkg in required_packages) {
    if(!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is needed for this function to work. Please install it."))
    }
  }
  
  # Set comp_type to lowercase for consistent handling
  comp_type = tolower(comp_type)
  
  # Validate comp_type
  if(!comp_type %in% c("length", "weight")) {
    stop("comp_type must be either 'length' or 'weight'")
  }
  
  # Set file names and plot titles based on comp_type
  if(comp_type == "length") {
    file_name = "comp_len.csv"
    kind_filter = "LEN"
    if(is.null(plot_name)) plot_name = "length_composition_comparison"
    x_label = "Length bin"
    y_label = "Proportion"
    plot_title = "Length Composition Comparison"
  } else { # weight
    file_name = "comp_size.csv"
    kind_filter = "SIZE"
    if(is.null(plot_name)) plot_name = "weight_composition_comparison"
    x_label = "Weight bin"
    y_label = "Proportion"
    plot_title = "Weight Composition Comparison"
  }
  
  # Validate inputs
  if(length(model_ids) < 1) {
    stop("Please provide at least one model ID")
  }
  
  # Handle model_stem as vector or single path
  if(length(model_stem) == 1) {
    model_stem = rep(model_stem, length(model_ids))
  } else if(length(model_stem) != length(model_ids)) {
    stop("If model_stem is a vector, it must have the same length as model_ids")
  }
  
  # If model_labels is not provided, use model_ids as labels
  if(is.null(model_labels)) {
    model_labels = model_ids
  } else if(length(model_labels) != length(model_ids)) {
    stop("If model_labels is provided, it must have the same length as model_ids")
  }
  
  # Create vector of model specifications with labels
  selected_models = data.frame(
    id = model_ids,
    path = model_stem,
    label = model_labels,
    stringsAsFactors = FALSE
  )
  
  # Check if composition file exists for all models
  file_exists = sapply(1:nrow(selected_models), function(i) {
    file.exists(file.path(selected_models$path[i], selected_models$id[i], file_name))
  })
  
  if(sum(file_exists) == 0) {
    stop(paste0("The file ", file_name, " does not exist for any of the models you selected."))
  } else if(sum(file_exists) < length(file_exists)) {
    warning(paste0("The file ", file_name, " does not exist for some models. These models will be excluded."))
    selected_models = selected_models[file_exists, ]
  }
  
  # Read and compile composition data
  comp_data_list = lapply(1:nrow(selected_models), function(i) {
    # Read composition data
    comp_file = file.path(selected_models$path[i], selected_models$id[i], file_name)
    comp_dt = fread(comp_file)
    
    # Filter by Kind
    comp_dt = comp_dt[Kind == kind_filter]
    
    # Add model identifiers
    comp_dt[, model_id := selected_models$id[i]]
    comp_dt[, model_label := selected_models$label[i]]
    
    # Round observed values for consistent ID generation
    comp_dt[, Obs := round(Obs, digits=6)]
    
    return(comp_dt)
  })
  
  # Combine all composition data
  plot_dt = rbindlist(comp_data_list, fill = TRUE)
  
  # Check if we have any data after filtering
  if(nrow(plot_dt) == 0) {
    stop(paste0("No ", comp_type, " composition data found with Kind = ", kind_filter))
  }
  
  # Calculate bin widths for each fleet
  bin_info = plot_dt[, .(
    Bin = sort(unique(Bin))
  ), by = .(Fleet_name)]
  
  # Calculate bin widths as differences between consecutive bins
  # For the last bin, use the width of the previous bin
  bin_info[, bin_width := c(diff(Bin), tail(diff(Bin), 1)), by = .(Fleet_name)]
  
  # Calculate bin centers and edges for plotting
  bin_info[, bin_center := Bin + bin_width/2]
  
  # Merge bin widths back to the main data table
  plot_dt = merge(plot_dt, bin_info, by = c("Fleet_name", "Bin"))
  
  # Check consistency of observed values across models for the same fleet/bin
  consistency_check = plot_dt[, .(
    obs_values = list(unique(Obs)),
    models = list(unique(model_label))
  ), by = .(Fleet_name, Bin)]
  
  # Check if observed values are consistent
  obs_inconsistent = consistency_check[sapply(obs_values, length) > 1]
  if(nrow(obs_inconsistent) > 0) {
    error_msg = paste0(
      "Observed ", comp_type, " composition values differ across models for the same fleet and bin:\n",
      paste(
        apply(obs_inconsistent[, .(Fleet_name, Bin)], 1, function(row) {
          return(paste0("  - Fleet: ", row[1], ", Bin: ", row[2]))
        }),
        collapse = "\n"
      )
    )
    stop(error_msg)
  }
  
  # Create a version of the data with just the first model's observations
  # This ensures we only plot the observations once
  obs_dt = plot_dt[, .SD[1], by = .(Fleet_name, Bin)]
  
  # Check if we have data to plot
  if(nrow(plot_dt) == 0) {
    stop("No data available after applying filters")
  }
  
  # Set model_label as factor with correct order
  plot_dt[, model_label := factor(model_label, levels = model_labels)]
  
  # Use Fleet_name as facet label
  plot_dt[, facet_label := Fleet_name]
  obs_dt[, facet_label := Fleet_name]
  
  # Combine data for plotting
  combined_obs = obs_dt[, .(Fleet_name, Bin, bin_width, bin_center, Obs, facet_label)]
  combined_exp = plot_dt[, .(Fleet_name, Bin, bin_width, bin_center, Exp, model_label, facet_label)]
  
  # Create a combined plot with facets
  combined_plot = ggplot() +
    # Add bars for observed data using geom_rect
    geom_rect(data = combined_obs, 
              aes(xmin = Bin, xmax = Bin + bin_width, 
                  ymin = 0, ymax = Obs),
              fill = "gray80", color = "gray40", alpha = 0.7) +
    # Add fitted lines if requested
    {if(show_fit) geom_line(data = combined_exp, 
                           aes(x = bin_center, y = Exp, color = model_label, group = model_label),
                           linewidth = 1)} +
    # Facet by fleet, with option for free y scales
    facet_wrap(~ facet_label, ncol = n_col, scales = if(free_y_scale) "free_y" else "fixed") +
    # Labels
    labs(x = x_label, y = y_label, title = plot_title)
  
  # Apply color scales - either custom colors or viridis palette
  color_count = length(unique(plot_dt$model_label))
  color_legend = "Model"
  
  if(!is.null(custom_colors)) {
    # Validate that enough colors are provided
    if(length(custom_colors) < color_count) {
      warning("Not enough custom colors provided. Falling back to viridis palette.")
      combined_plot = combined_plot + 
        viridis::scale_color_viridis(color_legend, begin = 0.1, end = 0.8, 
                                   direction = 1, option = "H", discrete = TRUE)
    } else {
      # Use custom colors
      combined_plot = combined_plot + scale_color_manual(color_legend, values = custom_colors)
    }
  } else {
    # Use default viridis palette
    combined_plot = combined_plot + 
      viridis::scale_color_viridis(color_legend, begin = 0.1, end = 0.8, 
                                 direction = 1, option = "H", discrete = TRUE)
  }
  
  # Apply consistent theme elements
  combined_plot = combined_plot + theme(
    panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
    panel.grid.major = element_line(color = 'gray70', linetype = "dotted"),
    panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
    strip.background = element_rect(fill = "white"),
    legend.key = element_rect(fill = "white")
  )
  
  # Save the plot if requested
  if(save_plot) {
    # Validate save directory
    if(is.null(save_dir)) {
      save_dir = getwd()
      warning(paste("No save directory specified. Saving to current working directory:", save_dir))
    } else {
      # Create directory if it doesn't exist
      if(!dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE)
        message(paste("Created directory:", save_dir))
      }
    }
    
    # Construct full file path with extension
    file_path = file.path(save_dir, paste0(plot_name, ".", file_type))
    
    # Save using appropriate function based on file type
    if(file_type == "pdf") {
      pdf(file_path, width = width, height = height)
      print(combined_plot)
      dev.off()
    } else {
      ggsave(filename = file_path, plot = combined_plot, width = width, height = height, 
             dpi = dpi, units = "in")
    }
    
    message(paste("Plot saved to:", file_path))
    invisible(combined_plot)
  } else {
    return(combined_plot)
  }
}

# Function to create comparative catch plots across multiple models
plot_catch_comparison = function(model_ids, model_stem,
                               fleets = NULL,        # New parameter for fleet selection
                               show_fit = TRUE,
                               n_col = 2,
                               free_y_scale = TRUE,  # Option to have free y scales
                               model_labels = NULL,
                               custom_colors = NULL,
                               save_plot = FALSE,
                               save_dir = NULL,
                               plot_name = NULL,
                               width = 10,
                               height = 8,
                               dpi = 300,
                               file_type = "png") {
  
  # Required packages
  required_packages <- c("data.table", "ggplot2", "viridis")
  for(pkg in required_packages) {
    if(!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is needed for this function to work. Please install it."))
    }
  }
  
  # Set file names and plot titles
  file_name = "catch.csv"
  if(is.null(plot_name)) plot_name = "catch_comparison"
  x_label = "Year"
  y_label = "Catch"
  plot_title = "Catch Comparison"
  
  # Validate inputs
  if(length(model_ids) < 1) {
    stop("Please provide at least one model ID")
  }
  
  # Handle model_stem as vector or single path
  if(length(model_stem) == 1) {
    model_stem = rep(model_stem, length(model_ids))
  } else if(length(model_stem) != length(model_ids)) {
    stop("If model_stem is a vector, it must have the same length as model_ids")
  }
  
  # If model_labels is not provided, use model_ids as labels
  if(is.null(model_labels)) {
    model_labels = model_ids
  } else if(length(model_labels) != length(model_ids)) {
    stop("If model_labels is provided, it must have the same length as model_ids")
  }
  
  # Create vector of model specifications with labels
  selected_models = data.frame(
    id = model_ids,
    path = model_stem,
    label = model_labels,
    stringsAsFactors = FALSE
  )
  
  # Check if catch file exists for all models
  file_exists = sapply(1:nrow(selected_models), function(i) {
    file.exists(file.path(selected_models$path[i], selected_models$id[i], file_name))
  })
  
  if(sum(file_exists) == 0) {
    stop(paste0("The file ", file_name, " does not exist for any of the models you selected."))
  } else if(sum(file_exists) < length(file_exists)) {
    warning(paste0("The file ", file_name, " does not exist for some models. These models will be excluded."))
    selected_models = selected_models[file_exists, ]
  }
  
  # Read and compile catch data
  catch_data_list = lapply(1:nrow(selected_models), function(i) {
    # Read catch data
    catch_file = file.path(selected_models$path[i], selected_models$id[i], file_name)
    catch_dt = fread(catch_file)
    
    # Add model identifiers
    catch_dt[, model_id := selected_models$id[i]]
    catch_dt[, model_label := selected_models$label[i]]
    
    # Round values for consistent ID generation
    catch_dt[, Obs := round(Obs, digits=6)]
    catch_dt[, Exp := round(Exp, digits=6)]
    
    # Round Time to avoid floating point issues
    catch_dt[, Time := round(Time, digits=2)]
    
    return(catch_dt)
  })
  
  # Combine all catch data
  plot_dt = rbindlist(catch_data_list, fill = TRUE)
  
  # Check if we have any data after filtering
  if(nrow(plot_dt) == 0) {
    stop("No catch data found")
  }
  
  # Get fleet names (if not in data, use Fleet number)
  if(!"Fleet_name" %in% names(plot_dt)) {
    # Try to create fleet names if not available
    plot_dt[, Fleet_name := paste("Fleet", Fleet)]
  }
  
  # Filter by fleet if the fleets parameter is provided
  if(!is.null(fleets)) {
    # Determine if fleets are specified by name or ID
    if(is.character(fleets)) {
      # Check if all specified fleet names exist
      missing_fleets <- fleets[!fleets %in% unique(plot_dt$Fleet_name)]
      if(length(missing_fleets) > 0) {
        warning(paste("The following specified fleet names were not found:", 
                     paste(missing_fleets, collapse=", ")))
      }
      
      # Filter by fleet name
      plot_dt <- plot_dt[Fleet_name %in% fleets]
    } else {
      # Check if all specified fleet IDs exist
      missing_fleets <- fleets[!fleets %in% unique(plot_dt$Fleet)]
      if(length(missing_fleets) > 0) {
        warning(paste("The following specified fleet IDs were not found:", 
                     paste(missing_fleets, collapse=", ")))
      }
      
      # Filter by fleet ID
      plot_dt <- plot_dt[Fleet %in% fleets]
    }
    
    # Check if we have any data after fleet filtering
    if(nrow(plot_dt) == 0) {
      stop("No data available after filtering for the specified fleets")
    }
  }
  
  # Set the order of fleet names based on the order in the fleets parameter
  if(!is.null(fleets)) {
    if(is.character(fleets)) {
      # For character fleet names, use the order as provided
      valid_fleets <- fleets[fleets %in% unique(plot_dt$Fleet_name)]
      plot_dt[, Fleet_name := factor(Fleet_name, levels = valid_fleets)]
    } else {
      # For numeric fleet IDs, we need to match to fleet names
      fleet_mapping <- unique(plot_dt[, .(Fleet, Fleet_name)])
      fleet_mapping <- fleet_mapping[Fleet %in% fleets]
      setkey(fleet_mapping, Fleet)
      # Order by the provided fleet IDs
      ordered_names <- fleet_mapping[J(fleets[fleets %in% unique(plot_dt$Fleet)]), Fleet_name]
      # Set as factor with levels in the right order
      plot_dt[, Fleet_name := factor(Fleet_name, levels = ordered_names)]
    }
  }
  
  # Calculate the appropriate bin width based on the smallest interval in the data
  # First, get the unique Time values and sort them
  unique_times = sort(unique(plot_dt$Time))
  
  # Calculate the differences between consecutive times
  time_diffs = diff(unique_times)
  
  # Use the smallest non-zero difference as our bin width
  # If no non-zero differences, default to 1
  non_zero_diffs = time_diffs[time_diffs > 0]
  if(length(non_zero_diffs) > 0) {
    bin_width = min(non_zero_diffs)
  } else {
    bin_width = 1
  }
  
  # Round to 6 decimal places to avoid floating point issues
  bin_width = round(bin_width, 6)
  
  # Add bin information to the data
  plot_dt[, bin_start := Time]
  plot_dt[, bin_end := Time + bin_width]
  plot_dt[, bin_center := Time + bin_width/2]
  
  # Check consistency of observed values across models for the same fleet/time
  consistency_check = plot_dt[, .(
    obs_values = list(unique(Obs)),
    models = list(unique(model_label))
  ), by = .(Fleet, Time)]
  
  # Check if observed values are consistent
  obs_inconsistent = consistency_check[sapply(obs_values, length) > 1]
  if(nrow(obs_inconsistent) > 0) {
    warning_msg = paste0(
      "Observed catch values differ across models for the same fleet and time:\n",
      paste(
        apply(obs_inconsistent[, .(Fleet, Time)], 1, function(row) {
          return(paste0("  - Fleet: ", row[1], ", Time: ", row[2]))
        }),
        collapse = "\n"
      ),
      "\nUsing the first model's observations for plotting."
    )
    warning(warning_msg)
  }
  
  # Create a version of the data with just the first model's observations
  # This ensures we only plot the observations once
  obs_dt = plot_dt[, .SD[1], by = .(Fleet, Fleet_name, Time)]
  
  # Make sure obs_dt has the same fleet name factor levels as plot_dt
  if("Fleet_name" %in% names(plot_dt) && is.factor(plot_dt$Fleet_name)) {
    obs_dt[, Fleet_name := factor(Fleet_name, levels = levels(plot_dt$Fleet_name))]
  }
  
  # Check if we have data to plot
  if(nrow(plot_dt) == 0) {
    stop("No data available after applying filters")
  }
  
  # Set model_label as factor with correct order
  plot_dt[, model_label := factor(model_label, levels = model_labels)]
  
  # Use Fleet_name as facet label
  plot_dt[, facet_label := Fleet_name]
  obs_dt[, facet_label := Fleet_name]
  
  # Combine data for plotting
  combined_obs = obs_dt[, .(Fleet, Fleet_name, Time, bin_start, bin_end, bin_center, 
                           Obs, facet_label)]
  combined_exp = plot_dt[, .(Fleet, Fleet_name, Time, bin_start, bin_end, bin_center, 
                            Exp, model_label, facet_label)]
  
  # Create a combined plot with facets
  combined_plot = ggplot() +
    # Add bars for observed data using geom_rect
    geom_rect(data = combined_obs, 
              aes(xmin = bin_start, xmax = bin_end, 
                  ymin = 0, ymax = Obs),
              fill = "gray80", color = "gray40", alpha = 0.7) +
    # Add fitted lines if requested
    {if(show_fit) geom_line(data = combined_exp, 
                           aes(x = bin_center, y = Exp, color = model_label, group = model_label),
                           linewidth = 1)} +
    # Facet by fleet, with option for free y scales
    facet_wrap(~ facet_label, ncol = n_col, scales = if(free_y_scale) "free_y" else "fixed") +
    # Labels
    labs(x = x_label, y = y_label, title = plot_title)
  
  # Apply color scales - either custom colors or viridis palette
  color_count = length(unique(plot_dt$model_label))
  color_legend = "Model"
  
  if(!is.null(custom_colors)) {
    # Validate that enough colors are provided
    if(length(custom_colors) < color_count) {
      warning("Not enough custom colors provided. Falling back to viridis palette.")
      combined_plot = combined_plot + 
        viridis::scale_color_viridis(color_legend, begin = 0.1, end = 0.8, 
                                   direction = 1, option = "H", discrete = TRUE)
    } else {
      # Use custom colors
      combined_plot = combined_plot + scale_color_manual(color_legend, values = custom_colors)
    }
  } else {
    # Use default viridis palette
    combined_plot = combined_plot + 
      viridis::scale_color_viridis(color_legend, begin = 0.1, end = 0.8, 
                                 direction = 1, option = "H", discrete = TRUE)
  }
  
  # Apply consistent theme elements
  combined_plot = combined_plot + theme(
    panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
    panel.grid.major = element_line(color = 'gray70', linetype = "dotted"),
    panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
    strip.background = element_rect(fill = "white"),
    legend.key = element_rect(fill = "white")
  )
  
  # Save the plot if requested
  if(save_plot) {
    # Validate save directory
    if(is.null(save_dir)) {
      save_dir = getwd()
      warning(paste("No save directory specified. Saving to current working directory:", save_dir))
    } else {
      # Create directory if it doesn't exist
      if(!dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE)
        message(paste("Created directory:", save_dir))
      }
    }
    
    # Construct full file path with extension
    file_path = file.path(save_dir, paste0(plot_name, ".", file_type))
    
    # Save using appropriate function based on file type
    if(file_type == "pdf") {
      pdf(file_path, width = width, height = height)
      print(combined_plot)
      dev.off()
    } else {
      ggsave(filename = file_path, plot = combined_plot, width = width, height = height, 
             dpi = dpi, units = "in")
    }
    
    message(paste("Plot saved to:", file_path))
    invisible(combined_plot)
  } else {
    return(combined_plot)
  }
}

# Function to create comparative biological plots (growth, weight, maturity)
plot_model_comparison_bio = function(model_ids, model_stem, 
                                   n_col = 2, model_labels = NULL, custom_colors = NULL,
                                   show_variability = TRUE, 
                                   save_plot = FALSE, save_dir = NULL, plot_name = "bio_comparison", 
                                   width = 12, height = 10, dpi = 300, file_type = "png") {
  
  # Validate inputs
  if(length(model_ids) < 1) {
    stop("Please provide at least one model ID")
  }

  # Handle model_stem as vector or single path
  if(length(model_stem) == 1) {
    model_stem = rep(model_stem, length(model_ids))
  } else if(length(model_stem) != length(model_ids)) {
    stop("If model_stem is a vector, it must have the same length as model_ids")
  }
  
  # If model_labels is not provided, use model_ids as labels
  if(is.null(model_labels)) {
    model_labels = model_ids
  } else if(length(model_labels) != length(model_ids)) {
    stop("If model_labels is provided, it must have the same length as model_ids")
  }
  
  # Create vector of model specifications with labels
  selected_models = data.frame(
    id = model_ids,
    path = model_stem,
    label = model_labels,
    stringsAsFactors = FALSE
  )
  
  # Read and process biological data from each model
  bio_list = lapply(1:nrow(selected_models), function(i) {
    # Read bio_age.csv data
    bio_file = file.path(selected_models$path[i], selected_models$id[i], "bio_age.csv")
    if (!file.exists(bio_file)) {
      warning(paste("File not found:", bio_file))
      return(NULL)
    }
    
    dt = fread(bio_file)
    dt[, model_id := selected_models$id[i]]
    dt[, model_label := selected_models$label[i]]
    return(dt)
  })
  
  # Combine all biological data
  bio_data = rbindlist(bio_list, use.names = TRUE, fill = TRUE)
  
  # Check if we have data
  if(nrow(bio_data) == 0) {
    stop("No biological data found for the specified models")
  }
  
  # Add sex label for better readability
  bio_data[sex == 1, sex_label := "Female"]
  bio_data[sex == 2, sex_label := "Male"]
  bio_data[is.na(sex_label), sex_label := "NA"]

  if(uniqueN(bio_data$sex)==1){
    bio_data[sex == 1, sex_label := "Aggregated"]
  }
  
  # Convert labels to factors to ensure consistent ordering
  bio_data[, model_label := factor(model_label, levels = model_labels)]
  # Sort unique sex labels alphabetically for factor levels
  sex_levels <- sort(unique(bio_data$sex_label))
  bio_data[, sex_label := factor(sex_label, levels = sex_levels)]
  
  # Calculate confidence intervals for length using len_sd
  bio_data[, len_upper := len + 1.96 * len_sd]
  bio_data[, len_lower := len - 1.96 * len_sd]
  # Ensure lower bound is not negative
  bio_data[len_lower < 0, len_lower := 0]
  
  # Create long format data for faceting on bio variable
  bio_long <- melt(bio_data, 
                   id.vars = c("model_id", "model_label", "id", "sex", "sex_label", "age"),
                   measure.vars = c("len", "wt", "mat"),
                   variable.name = "bio_var",
                   value.name = "value")
  
  # Add confidence intervals for length to long format
  bio_long_with_ci <- bio_long[, .(model_id, model_label, id, sex, sex_label, age, bio_var, value)]
  
  # Create separate data for the confidence intervals (only for length)
  ci_data <- bio_data[, .(model_id, model_label, id, sex, sex_label, age, 
                          bio_var = "len", 
                          upper = len_upper, 
                          lower = len_lower)]
  
  # Add readable labels for facets
  bio_long_with_ci[bio_var == "len", bio_label := "Length at Age (cm)"]
  ci_data[bio_var == "len", bio_label := "Length at Age (cm)"]
  bio_long_with_ci[bio_var == "wt", bio_label := "Weight at Age (kg)"]
  bio_long_with_ci[bio_var == "mat", bio_label := "Maturity at Age"]
  
  # Convert bio_label to factor with desired order
  bio_var_levels <- c("Length at Age (cm)", "Weight at Age (kg)", "Maturity at Age")
  bio_long_with_ci[, bio_label := factor(bio_label, levels = bio_var_levels)]
  ci_data[, bio_label := factor(bio_label, levels = bio_var_levels)]
  
  # Create the plot
  p = ggplot() +
    # Add lines for each model and sex combination
    geom_line(data = bio_long_with_ci, 
              aes(x = age, y = value, color = model_label, linetype = sex_label),
              linewidth = 1.2) +
    # Add facets for different biological variables with free y scales
    facet_wrap(~ bio_label, scales = "free_y", ncol = n_col) +
    # Set y-axis to start at 0
    ylim(0, NA) +
    xlab("Age (years)") +
    ylab("Value") +
    labs(color = "Model", linetype = "Sex") +
    theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
          panel.grid.major = element_line(color = 'gray70', linetype = "dotted"),
          panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
          strip.background = element_rect(fill = "white"),
          legend.key = element_rect(fill = "white"))
  
  # Add variability ribbons for length if requested, but only in the length at age facet
  if(show_variability) {
    # Filter the ci_data to include only bio_var == "len" for the ribbon
    length_ci_data <- ci_data[bio_var == "len"]
    
    p = p + geom_ribbon(data = length_ci_data,
                        aes(x = age, ymin = lower, ymax = upper, 
                            fill = model_label, group = interaction(model_label, sex_label)),
                        alpha = 0.2, show.legend = FALSE)
  }
  
  # Apply color scales - either custom colors or viridis palette
  if(!is.null(custom_colors)) {
    # Validate that enough colors are provided
    if(length(custom_colors) < length(unique(bio_long_with_ci$model_label))) {
      warning("Not enough custom colors provided. Falling back to viridis palette.")
      p = p + viridis::scale_color_viridis("Model", begin = 0.1, end = 0.8, 
                                         direction = 1, option = "H", discrete = TRUE) +
             viridis::scale_fill_viridis("Model", begin = 0.1, end = 0.8, 
                                       direction = 1, option = "H", discrete = TRUE)
    } else {
      # Use custom colors
      p = p + scale_color_manual("Model", values = custom_colors) +
             scale_fill_manual("Model", values = custom_colors)
    }
  } else {
    # Use default viridis palette
    p = p + viridis::scale_color_viridis("Model", begin = 0.1, end = 0.8, 
                                       direction = 1, option = "H", discrete = TRUE) +
           viridis::scale_fill_viridis("Model", begin = 0.1, end = 0.8, 
                                     direction = 1, option = "H", discrete = TRUE)
  }
  
  # Save the plot if requested
  if(save_plot) {
    # Validate save directory
    if(is.null(save_dir)) {
      save_dir = getwd()
      warning(paste("No save directory specified. Saving to current working directory:", save_dir))
    } else {
      # Create directory if it doesn't exist
      if(!dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE)
        message(paste("Created directory:", save_dir))
      }
    }
    
    # Construct full file path with extension
    file_path = file.path(save_dir, paste0(plot_name, ".", file_type))
    
    # Save using appropriate function based on file type
    if(file_type == "pdf") {
      pdf(file_path, width = width, height = height)
      print(p)
      dev.off()
    } else {
      ggsave(filename = file_path, plot = p, width = width, height = height, 
            dpi = dpi, units = "in")
    }
    
    message(paste("Plot saved to:", file_path))
    return(invisible(p))
  } else {
    return(p)
  }
}
