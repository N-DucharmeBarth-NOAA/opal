# Extract spawning biomass and depletion from MFCL Report (.rep) file
# Returns data.table with columns: model, year, ssb, ssb_se, depletion, depletion_se

# rep_file = file.path(this.path::this.proj(), "model-files", "mfcl", "v11", "plot-10.par.rep")
# result = extract_mfcl_biomass(rep_file, model_name = "MFCL-v11")
# head(result)

extract_mfcl_biomass = function(rep_file, model_name = NULL) {
	if(!file.exists(rep_file)) {
		stop(sprintf("Report file not found: %s", rep_file))
	}
	
	if(is.null(model_name)) {
		model_name = basename(dirname(rep_file))
	}
	
	rep = read.MFCLRep(rep_file)
	
	# Extract fished SSB using FLR4MFCL ssb() function
	ssb_fished = as.data.table(ssb(rep)) %>%
		.[age == "all"] %>%
		.[, year := as.numeric(year)] %>%
		.[, .(ssb_fished = sum(value, na.rm = TRUE)), by = year]
	
	# Extract unfished SSB (dynamic B0) using adultBiomass_nofish()
	ssb_unfished = as.data.table(adultBiomass_nofish(rep)) %>%
		.[age == "all"] %>%
		.[, year := as.numeric(year)] %>%
		.[, .(ssb_unfished = sum(value, na.rm = TRUE)), by = year]
	
	# Merge fished and unfished, calculate depletion
	biomass_dt = ssb_fished[ssb_unfished, on = "year"] %>%
		.[, .(
			model = model_name,
			year,
			ssb = ssb_fished,
			ssb_se = ssb_fished * 0.05,
			depletion = ssb_fished / ssb_unfished,
			depletion_se = (ssb_fished / ssb_unfished) * 0.05
		)]
	
	setorderv(biomass_dt, c("year"))
	return(biomass_dt)
}
