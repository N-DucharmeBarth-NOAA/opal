# Nicholas Ducharme-Barth
# 2026/01/22
# R code to fix selectivity parameters at bounds from 01-bet-base model
# Reads Report.sso, identifies parameters flagged HI/LO, and expands bounds

# Copyright (c) 2026 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# load packages
	library(data.table)
	library(magrittr)
	library(r4ss)

#_____________________________________________________________________________________________________________________________
# define paths
	proj_dir = this.path::this.proj()
	dir_model = file.path(proj_dir,"model-files")
    dir_base_stock_synthesis = file.path(dir_model,"ss3","01-bet-base")
    dir_helper_fns = file.path(proj_dir,"code","ss3","helper-fns")

#________________________________________________________________________________________________________________________________________________________________________________________________________
# source helper functions
    sapply(file.path(dir_helper_fns,(list.files(dir_helper_fns))),source)

#_____________________________________________________________________________________________________________________________
# read in baseline stock synthesis files from 01-bet-base
    tmp_starter = SS_readstarter(file=file.path(dir_base_stock_synthesis,"starter.ss_new"),verbose=FALSE)
    tmp_ctl = SS_readctl(file=file.path(dir_base_stock_synthesis,"control.ss_new"),datlist = file.path(dir_base_stock_synthesis,"data_echo.ss_new"))
    tmp_data = SS_readdat(file=file.path(dir_base_stock_synthesis,"data_echo.ss_new"))
    tmp_forecast = SS_readforecast(file=file.path(dir_base_stock_synthesis,"forecast.ss_new"),verbose=FALSE)

#_____________________________________________________________________________________________________________________________
# read Report.sso to identify parameters at bounds
    report_file = file.path(dir_base_stock_synthesis,"Report.sso")
    report_lines = readLines(report_file)
    
    # Find parameters flagged as HI or LO in Report (parameter report section)
    # Look for Status column containing "HI" or "LO" as separate word
    flag_idx = grep("\\s(HI|LO)\\s", report_lines)
    
    if(length(flag_idx)==0){
      stop("No parameters flagged HI/LO found in Report.sso. Check bounds manually.")
    }
    
    # Extract parameter names and bounds from Report.sso lines
    flagged_params = data.table()
    for(idx in flag_idx){
      line = report_lines[idx]
      parts = strsplit(trimws(line), "\\s+")[[1]]
      if(length(parts) >= 10){
        param_name = parts[2]
        param_value = as.numeric(parts[3])
        param_phase = as.numeric(parts[5])
        min_bound = as.numeric(parts[6])
        max_bound = as.numeric(parts[7])
        init_value = as.numeric(parts[8])
        status = parts[10]
        flagged_params = rbind(flagged_params, 
          data.table(name=param_name, value=param_value, phase=param_phase, lo=min_bound, hi=max_bound, init=init_value, status=status))
      }
    }
    
    message(sprintf("Found %d parameters at bounds:", nrow(flagged_params)))
    print(flagged_params)

#_____________________________________________________________________________________________________________________________
# expand bounds for problematic parameters
    # Strategy: 
    # - For ascend_se/descend_se with bounds -7 to 7: fix if init value outside [-6, 6]
    # - For end_logit with bounds -999 to 9: change lower bound to -9 if being estimated
    # - For Size_95%width with bounds 0 to 500: change to 0 to 50
    
    for(i in seq_len(nrow(flagged_params))){
      param_name = flagged_params$name[i]
      old_lo = flagged_params$lo[i]
      old_hi = flagged_params$hi[i]
      init_val = flagged_params$init[i]
      param_phase = flagged_params$phase[i]
      
      new_lo = old_lo
      new_hi = old_hi
      new_phase = param_phase
      
      if(old_lo == -7 && old_hi == 7){
        # ascend_se or descend_se parameter: fix if init value outside [-6, 6]
        if(init_val < -6 || init_val > 6){
          new_phase = -abs(param_phase)  # fix at current value
          message(sprintf("  Fixing %s at init value %.3f (outside [-6, 6])", param_name, init_val))
        } else {
          # Parameter is within reasonable bounds, skip it
          next
        }
      } else if(old_lo == 0 && old_hi == 500){
        # Size_95%width parameter: change bounds to 0 to 50
        new_lo = 0
        new_hi = 50
      } else if(grepl("end_logit", param_name) && (old_lo == -999 || old_lo == -999.0) && param_phase > 0){
        # end_logit parameter being estimated: change lower bound from -999 to -9
        new_lo = -9
        new_hi = old_hi
      } else if(old_lo == -999 || old_lo == -999.0){
        # Skip start_logit and other logit parameters with bounds -999 to 9
        next
      } else {
        next
      }
      
      # Map Report.sso parameter name to selex_parms rowname
      # Report format: Size_DblN_PARTYPE_FLEETNAME(id) or Size_inflection_FLEETNAME or Size_95%width_FLEETNAME
      # selex_parms format: SizeSel_P_N_FLEETNAME(id)
      
      # Extract fleet name and number from param_name
      # e.g., "Size_DblN_ascend_se_F04_LL.EQUAT(4)" or "Size_95%width_F11_DOM.HL(11)"
      fleet_match = regmatches(param_name, regexpr("F\\d+_[^(]+", param_name))
      fleet_num_match = regmatches(param_name, regexpr("\\(\\d+\\)", param_name))
      
      if(length(fleet_match) == 0 || length(fleet_num_match) == 0){
        warning(sprintf("Could not parse fleet from parameter %s", param_name))
        next
      }
      
      fleet_name = fleet_match
      fleet_num = as.numeric(gsub("[()]", "", fleet_num_match))
      
      # Determine parameter number (P_N) based on parameter type
      param_type = NA
      if(grepl("peak", param_name)) param_type = 1
      else if(grepl("top_logit", param_name)) param_type = 2
      else if(grepl("ascend_se", param_name)) param_type = 3
      else if(grepl("descend_se", param_name)) param_type = 4
      else if(grepl("start_logit", param_name)) param_type = 5
      else if(grepl("end_logit", param_name)) param_type = 6
      else if(grepl("inflection", param_name)) param_type = 1  # for logistic pattern
      else if(grepl("95%width", param_name)) param_type = 2   # for logistic pattern
      
      if(is.na(param_type)){
        warning(sprintf("Could not determine parameter type for %s", param_name))
        next
      }
      
      # Build the rowname as it appears in size_selex_parms
      rowname_target = sprintf("SizeSel_P_%d_%s(%d)", param_type, fleet_name, fleet_num)
      
      # Find and update in size_selex_parms
      param_idx = which(rownames(tmp_ctl$size_selex_parms) == rowname_target)
      if(length(param_idx) > 0){
        # Update bounds and/or phase
        tmp_ctl$size_selex_parms[param_idx, 1] = new_lo
        tmp_ctl$size_selex_parms[param_idx, 2] = new_hi
        tmp_ctl$size_selex_parms[param_idx, 7] = new_phase
        
        if(new_phase != param_phase){
          message(sprintf("  %s [%s]: PHASE %d -> %d", 
            param_name, rowname_target, param_phase, new_phase))
        } else {
          message(sprintf("  %s [%s]: LO %f -> %f; HI %f -> %f", 
            param_name, rowname_target, old_lo, new_lo, old_hi, new_hi))
        }
      } else {
        warning(sprintf("Parameter %s (looking for %s) not found in size_selex_parms", param_name, rowname_target))
      }
    }

#_____________________________________________________________________________________________________________________________
# Also check ALL selectivity parameters with bounds -7 to 7 (not just flagged ones)
    message("\nChecking all selectivity parameters with bounds -7 to 7 for extreme initial values...")
    
    for(i in 1:nrow(tmp_ctl$size_selex_parms)){
      row_name = rownames(tmp_ctl$size_selex_parms)[i]
      old_lo = tmp_ctl$size_selex_parms[i, 1]
      old_hi = tmp_ctl$size_selex_parms[i, 2]
      init_val = tmp_ctl$size_selex_parms[i, 3]
      param_phase = tmp_ctl$size_selex_parms[i, 7]
      
      # Check if this parameter has bounds -7 to 7 and is being estimated
      if(old_lo == -7 && old_hi == 7 && param_phase > 0){
        # If init value is outside [-6, 6], fix the parameter
        if(init_val < -6 || init_val > 6){
          new_phase = -abs(param_phase)
          tmp_ctl$size_selex_parms[i, 7] = new_phase
          message(sprintf("  Fixing %s: init=%.3f, phase %d -> %d", 
                         row_name, init_val, param_phase, new_phase))
        }
      }
    }

#_____________________________________________________________________________________________________________________________
# update virgin recruitment
    tmp_ctl$SR_parms["SR_LN(R0)","INIT"] = 15

#_____________________________________________________________________________________________________________________________
# setup lambdas (likelihood weights) for all data types
    tmp_lambdas_surv = data.table(like_comp=rep(1,length(tmp_ctl$fleetnames)),fleet=1:length(tmp_ctl$fleetnames)) %>%
                       .[,phase:=1] %>%
                       .[,value:=0] %>%
                       .[,sizefreq_method:=0] %>%
                       .[fleet %in% c(15),value:=1]
    
    # Get unique fleets from length composition data
    lencomp_fleets = unique(tmp_data$lencomp$fleet)
    tmp_lambdas_lf = data.table(like_comp=rep(4,length(tmp_ctl$fleetnames)),fleet=1:length(tmp_ctl$fleetnames)) %>%
                     .[,phase:=1] %>%
                     .[,value:=0] %>%
                     .[,sizefreq_method:=0] %>%
                     .[fleet %in% lencomp_fleets,value:=1]
    
    # Get unique fleets from generalized size comp (weight) data
    sizefreq_fleets = unique(tmp_data$sizefreq_data_list[[1]]$fleet)
    tmp_lambdas_gs = data.table(like_comp=rep(6,length(sizefreq_fleets)),fleet=sizefreq_fleets) %>%
                     .[,phase:=1] %>%
                     .[,value:=1] %>%
                     .[,sizefreq_method:=1]
    
    # Combine all lambda components
    tmp_ctl$lambdas = as.data.frame(rbind(tmp_lambdas_surv,tmp_lambdas_lf,tmp_lambdas_gs))
    tmp_ctl$N_lambdas = nrow(tmp_ctl$lambdas)

#_____________________________________________________________________________________________________________________________
# create new directory for fixed stock synthesis model
    dir_fixed_stock_synthesis = file.path(dir_model,"ss3","02-fix-sel")
    dir.create(dir_fixed_stock_synthesis,recursive=TRUE)

#_____________________________________________________________________________________________________________________________
# write out all files to new directory
    SS_writestarter(tmp_starter, dir=dir_fixed_stock_synthesis, overwrite=TRUE)
    SS_writedat(tmp_data, outfile=file.path(dir_fixed_stock_synthesis,"data.ss"), overwrite=TRUE)
    SS_writectl(tmp_ctl, file.path(dir_fixed_stock_synthesis,"control.ss"), overwrite=TRUE)
    SS_writeforecast(tmp_forecast, dir=dir_fixed_stock_synthesis, overwrite=TRUE)

#_____________________________________________________________________________________________________________________________
# run the model
    file.copy(from=paste0(proj_dir,"/executables/stock-synthesis/3.30.24.1/ss3_win.exe"), to=dir_fixed_stock_synthesis)
    run(dir=dir_fixed_stock_synthesis, exe="ss3_win.exe", show_in_console=TRUE, skipfinished=FALSE)
