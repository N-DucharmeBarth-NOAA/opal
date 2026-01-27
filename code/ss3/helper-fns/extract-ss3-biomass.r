# Extract spawning biomass and depletion from SS3 Report.sso
# Returns data.table with columns: model, year, ssb, ssb_se, depletion, depletion_se

# model_dir = file.path(this.path::this.proj(), "model-files", "ss3", "01-bet-base")
# result = extract_ss3_biomass(model_dir)
# head(result)

extract_ss3_biomass = function(model_dir) {
	if(!file.exists(file.path(model_dir, "Report.sso"))) {
		stop(sprintf("Report.sso not found in %s", model_dir))
	}
	
	rep = SS_output(dir = model_dir, verbose = FALSE, printstats = FALSE)
	model_name = basename(model_dir)
	
	if(is.null(rep$derived_quants)) {
		stop("No derived_quants found in Report.sso")
	}
	
	quants = as.data.table(rep$derived_quants) %>%
		.[grep("^SSB_|^Bratio_", Label)]
	
	# Extract year only for labels matching pattern SSB_<number> or Bratio_<number>
	ssb_rows = quants[grep("^SSB_\\d+$", Label)]
	ssb_rows[, c("model", "year", "ssb", "ssb_se") := .(
		model_name,
		as.numeric(sub("SSB_", "", Label)),
		Value,
		StdDev
	)]
	
	depletion_rows = quants[grep("^Bratio_\\d+$", Label)]
	depletion_rows[, c("model", "year", "depletion", "depletion_se") := .(
		model_name,
		as.numeric(sub("Bratio_", "", Label)),
		Value,
		StdDev
	)]
	
	out = ssb_rows[, .(model, year, ssb, ssb_se)][
		depletion_rows[, .(model, year, depletion, depletion_se)],
		on = c("model", "year"), nomatch = 0
	]
	
	setorderv(out, c("year"))
	return(out)
}
