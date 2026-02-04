# Extract spawning biomass and depletion from SS3 Report.sso
# Returns data.table with columns: model, year, ts, season, ssb, ssb_se, depletion, depletion_se
# When quarterly=TRUE: ts is sequential (1,2,3,4,5,...), season is 1-4 per year, year increments per 4 timesteps
# When quarterly=FALSE: ts is sequential by year (1,2,3,...), season is 1 for all, ts and year correspond

# library(r4ss)
# library(data.table)
# library(magrittr)
# library(this.path)
# model_dir = file.path(this.path::this.proj(), "model-files", "ss3", "01-bet-base")
# result = extract_ss3_biomass(model_dir, start_year = 1952, quarterly = TRUE)
# head(result)

extract_ss3_biomass = function(model_dir, start_year = 1952, quarterly = TRUE) {
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
	
	# Extract only numeric year labels (SSB_1, SSB_2, etc., not SSB_Virgin)
	ssb_rows = quants[grep("^SSB_\\d+$", Label)]
	ts_vals = as.numeric(sub("SSB_", "", ssb_rows$Label))
	
	# Always compute quarterly structure first
	ssb_rows[, c("model", "ts", "season", "year", "ssb", "ssb_se") := .(
		model_name,
		ts_vals,
		((ts_vals - 1L) %% 4L) + 1L,
		start_year + ((ts_vals - 1L) %/% 4L),
		Value,
		StdDev
	)]
	
	# If not quarterly, aggregate within years by averaging
	if(!quarterly) {
		ssb_rows = ssb_rows[, .(
			ssb = mean(ssb, na.rm = TRUE),
			ssb_se = mean(ssb_se, na.rm = TRUE)
		), by = .(model, year)][, `:=`(
			ts = year - start_year + 1L,
			season = 1L
		)]
		ssb_rows = ssb_rows[, .(model, year, ts, season, ssb, ssb_se)]
	}
	
	depletion_rows = quants[grep("^Bratio_\\d+$", Label)]
	ts_vals_dep = as.numeric(sub("Bratio_", "", depletion_rows$Label))
	
	# Always compute quarterly structure first
	depletion_rows[, c("model", "ts", "season", "year", "depletion", "depletion_se") := .(
		model_name,
		ts_vals_dep,
		((ts_vals_dep - 1L) %% 4L) + 1L,
		start_year + ((ts_vals_dep - 1L) %/% 4L),
		Value,
		StdDev
	)]
	
	# If not quarterly, aggregate within years by averaging
	if(!quarterly) {
		depletion_rows = depletion_rows[, .(
			depletion = mean(depletion, na.rm = TRUE),
			depletion_se = mean(depletion_se, na.rm = TRUE)
		), by = .(model, year)][, `:=`(
			ts = year - start_year + 1L,
			season = 1L
		)]
		depletion_rows = depletion_rows[, .(model, year, ts, season, depletion, depletion_se)]
	}
	
	out = merge(
		ssb_rows[, .(model, year, ts, season, ssb, ssb_se)],
		depletion_rows[, .(model, year, ts, season, depletion, depletion_se)],
		by = c("model", "year", "ts", "season"),
		all.x = TRUE
	)
	
	setorderv(out, c("ts"))
	return(out)
}
