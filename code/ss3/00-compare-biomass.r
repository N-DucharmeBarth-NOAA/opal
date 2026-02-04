# Nicholas Ducharme-Barth
# 2026/01/27
# R code to compare spawning biomass and depletion across SS3 and MFCL models

# Copyright (c) 2026 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# load packages
	library(data.table)
	library(magrittr)
	library(FLR4MFCL)
	library(r4ss)
	library(ggplot2)

#_____________________________________________________________________________________________________________________________
# define paths
	proj_dir = this.path::this.proj()
	dir_model = file.path(proj_dir,"model-files")
	dir_ss3 = file.path(dir_model,"ss3")
	dir_mfcl = file.path(dir_model,"mfcl")
	dir_helper_fns_ss3 = file.path(proj_dir,"code","ss3","helper-fns")
	dir_helper_fns_mfcl = file.path(proj_dir,"code","mfcl","helper-fns")

#________________________________________________________________________________________________________________________________________________________________________________________________________
# source helper functions
	sapply(file.path(dir_helper_fns_ss3,(list.files(dir_helper_fns_ss3))),source)
	sapply(file.path(dir_helper_fns_mfcl,(list.files(dir_helper_fns_mfcl))),source)

#_____________________________________________________________________________________________________________________________
# extract SS3 biomass from specified models
	ss3_models = c("01-bet-base", "02-fix-sel")
	ss3_biomass = rbindlist(lapply(
		file.path(dir_ss3, ss3_models),
		function(x) extract_ss3_biomass(x, quarterly = FALSE)
	))

#_____________________________________________________________________________________________________________________________
# extract MFCL biomass
	mfcl_rep_file = file.path(dir_mfcl, "v11", "plot-10.par.rep")
	mfcl_biomass = extract_mfcl_biomass(mfcl_rep_file, model_name = "MFCL-v11", quarterly = FALSE)

#_____________________________________________________________________________________________________________________________
# combine SS3 and MFCL data
	all_biomass = rbind(ss3_biomass, mfcl_biomass)

#_____________________________________________________________________________________________________________________________
# plot spawning biomass comparison
	p_ssb = ggplot(all_biomass, aes(x = year, y = ssb, color = model, fill = model)) +
		geom_line(size = 0.8) +
		geom_ribbon(aes(ymin = ssb - 1.96*ssb_se, ymax = ssb + 1.96*ssb_se, color = NULL), 
		            alpha = 0.2) +
		labs(title = "Spawning Biomass Comparison", x = "Year", y = "Spawning Biomass (t)", 
		     color = "Model", fill = "Model") +
		theme_minimal() +
		theme(legend.position = "right", panel.grid.major = element_line(color = "gray90"))
	
	print(p_ssb)

#_____________________________________________________________________________________________________________________________
# plot depletion comparison
	p_depl = ggplot(all_biomass, aes(x = year, y = depletion, color = model, fill = model)) +
		geom_line(size = 0.8) +
		geom_ribbon(aes(ymin = depletion - 1.96*depletion_se, ymax = depletion + 1.96*depletion_se, color = NULL), 
		            alpha = 0.2) +
		geom_hline(yintercept = 0.25, linetype = "dashed", color = "red", size = 0.6) +
		geom_hline(yintercept = 0.35, linetype = "dashed", color = "orange", size = 0.6) +
		labs(title = "Depletion (B/B0) Comparison", x = "Year", y = "Depletion (B/B0)", 
		     color = "Model", fill = "Model") +
		theme_minimal() +
		theme(legend.position = "right", panel.grid.major = element_line(color = "gray90"))
	
	print(p_depl)

