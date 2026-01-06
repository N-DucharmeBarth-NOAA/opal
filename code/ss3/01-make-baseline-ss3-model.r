

# Nicholas Ducharme-Barth
# 2026/01/06
# R code to make a baseline stock synthesis model version of the mfcl/v11 model

# Copyright (c) 2026 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# load packages
	library(data.table)
	library(magrittr)
	library(FLR4MFCL)
    library(frqit)
    library(r4ss)

#_____________________________________________________________________________________________________________________________
# define paths
	proj_dir = this.path::this.proj()
	dir_model = paste0(proj_dir,"/model-files/")
    dir_base_mfcl = paste0(dir_model,"mfcl/v11/")
    dir_base_stock_synthesis = paste0(dir_model,"ss3/00-swpo-mls-base-file/")
    dir_helper_fns = paste0(proj_dir,"/code/helper-fns/")

#________________________________________________________________________________________________________________________________________________________________________________________________________
# source helper functions
    sapply(paste0(dir_helper_fns,(list.files(dir_helper_fns))),source)

#_____________________________________________________________________________________________________________________________
# read in baseline mfcl files
    base_frq = parse_frq(paste0(dir_base_mfcl,"bet.frq"))
    base_ini = read.MFCLIni(paste0(dir_base_mfcl,"bet.ini"), nseasons=4)
    base_par = read.MFCLPar(paste0(dir_base_mfcl,"10.par"), first.yr=1952)
    base_rep = read.MFCLRep(paste0(dir_base_mfcl,"plot-10.par.rep"))

#_____________________________________________________________________________________________________________________________
# run new version of stock synthesis
    file.copy(from=paste0(proj_dir,"/executables/stock-synthesis/3.30.24.1/ss3_win.exe"),to=dir_base_stock_synthesis)
    run(dir=dir_base_stock_synthesis,exe="ss3_win.exe")

#_____________________________________________________________________________________________________________________________
# read in baseline stock synthesis files
    tmp_starter = SS_readstarter(file=paste0(dir_base_stock_synthesis,"starter.ss_new"),verbose=FALSE)
    tmp_ctl = SS_readctl(file=paste0(dir_base_stock_synthesis,"control.ss_new"),datlist = paste0(dir_base_stock_synthesis,"data_echo.ss_new"))
    tmp_data = SS_readdat(file=paste0(dir_base_stock_synthesis,"data_echo.ss_new"))
    tmp_forecast = SS_readforecast(file=paste0(dir_base_stock_synthesis,"forecast.ss_new"),verbose=FALSE)

#_____________________________________________________________________________________________________________________________
# create new directory for stock synthesis mls files
    dir_bet_stock_synthesis_base = paste0(dir_model,"ss3/01-bet-base/")
    dir.create(dir_bet_stock_synthesis_base,recursive=TRUE)
