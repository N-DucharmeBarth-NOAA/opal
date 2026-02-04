# Extract spawning biomass and depletion from MFCL Report (.rep) file
# Returns data.table with columns: model, year, ts, season, ssb, ssb_se, depletion, depletion_se
# When quarterly=TRUE: ts is sequential (1,2,3,4,5,...), season is 1-4 per year
# When quarterly=FALSE: ts is sequential by year (1,2,3,...), season is 1 for all

# library(FLR4MFCL)
# library(data.table)
# library(magrittr)
# library(this.path)
# rep_file = file.path(this.path::this.proj(), "model-files", "mfcl", "v11", "plot-10.par.rep")
# result = extract_mfcl_biomass(rep_file, model_name = "MFCL-v11", quarterly = TRUE)
# head(result)

extract_mfcl_biomass = function(rep_file, model_name = NULL, quarterly = TRUE) {
	if(!file.exists(rep_file)) {
		stop(sprintf("Report file not found: %s", rep_file))
	}
	
	if(is.null(model_name)) {
		model_name = basename(dirname(rep_file))
	}
	
	rep = read.MFCLRep(rep_file)
	
	# Extract fished SSB using FLR4MFCL ssb() function
	ssb_values = as.data.table(ssb(rep))
	
	if(quarterly) {
		ssb_fished = ssb_values[age == "all", .(
			Value = mean(value, na.rm = TRUE)
		), by = .(yr = as.numeric(year), season = as.numeric(season), unit, area, iter)]
		
		# Calculate sequential timestep: ts = (years_from_start) * 4 + season
		ssb_fished[, ts := (yr - min(yr)) * 4L + season]
		ssb_fished = ssb_fished[, .(yr, ts, season, Value)]
	} else {
		ssb_fished = ssb_values[age == "all", .(
			Value = mean(value, na.rm = TRUE)
		), by = .(yr = as.numeric(year), unit, area, iter)]
		
		# Annual timestep indexes by year
		ssb_fished[, ts := yr - min(yr) + 1L]
		ssb_fished[, season := 1L]
		ssb_fished = ssb_fished[, .(yr, ts, season, Value)]
	}
	
	# Extract unfished SSB (dynamic B0) using adultBiomass_nofish()
	ssb_nofishing = as.data.table(adultBiomass_nofish(rep))
	
	if(quarterly) {
		ssb_unfished = ssb_nofishing[age == "all", .(
			Value = mean(value, na.rm = TRUE)
		), by = .(yr = as.numeric(year), season = as.numeric(season), unit, area, iter)]
		
		# Calculate sequential timestep: ts = (years_from_start) * 4 + season
		ssb_unfished[, ts := (yr - min(yr)) * 4L + season]
		ssb_unfished = ssb_unfished[, .(yr, ts, season, Value)]
	} else {
		ssb_unfished = ssb_nofishing[age == "all", .(
			Value = mean(value, na.rm = TRUE)
		), by = .(yr = as.numeric(year), unit, area, iter)]
		
		# Annual timestep indexes by year
		ssb_unfished[, ts := yr - min(yr) + 1L]
		ssb_unfished[, season := 1L]
		ssb_unfished = ssb_unfished[, .(yr, ts, season, Value)]
	}
	
	# Merge fished and unfished, calculate depletion
	biomass_dt = merge(ssb_fished, ssb_unfished, by = c("yr", "ts", "season"), suffixes = c("_fished", "_unfished")) %>%
		.[, .(
			model = model_name,
			year = yr,
			ts,
			season,
			ssb = Value_fished,
			ssb_se = NA_real_,
			depletion = Value_fished / Value_unfished,
			depletion_se = NA_real_
		)]
	
	setorderv(biomass_dt, c("year"))
	return(biomass_dt)
}
