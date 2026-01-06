

# Nicholas Ducharme-Barth
# 2024/02/03
# Wrapper function for SS_plots to make html dashboard from reduced outputs

# Copyright (c) 2024 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

# dir = to_dir

ss_model_html_viewer = function(dir,html_dir=file.path(this.path::this.proj(),"html-dashboard"),linux_ws=is_linux_os(),verbose=TRUE,catchasnumbers = TRUE)
{   
    A = proc.time()
    # copy all files from dir to html dir
    dir.create(html_dir,recursive=TRUE)
    file.copy(from=file.path(dir,list.files(dir)), to=file.path(html_dir,list.files(dir)), overwrite = TRUE)

    # remake output
    if(file.exists(file.path(html_dir,"data.ss_new"))){
        file.rename(from=file.path(html_dir,"data.ss_new"),to=file.path(html_dir,"data_echo.ss_new"))
    }

    exe_path = readLines(file.path(html_dir,"executable.txt"))
    exec_version = strsplit(exe_path,"/")[[1]][1]
    
    if(exec_version=="3.24U"|exec_version=="3.30.22.1"){
        make_ss_output(html_dir,par="ss3.par",linux_ws=linux_ws,verbose=TRUE) 
    } else {
        make_ss_output(html_dir,linux_ws=linux_ws,verbose=TRUE) 
    }

    tmp_report = SS_output(html_dir,verbose = FALSE,printstats = FALSE) 

    tmp_report$parameters = read.table(file.path(html_dir,"html_parameters.txt"))
    row.names(tmp_report$parameters) = readLines(file.path(html_dir,"html_parameters_rownames.txt"))

    tmp_report$recruitpars = read.table(file.path(html_dir,"html_recruitpars.txt"))
    row.names(tmp_report$recruitpars) = readLines(file.path(html_dir,"html_recruitpars_rownames.txt"))

    tmp_report$estimated_non_dev_parameters = read.table(file.path(html_dir,"html_estimated_non_dev_parameters.txt"))
    row.names(tmp_report$estimated_non_dev_parameters) = readLines(file.path(html_dir,"html_estimated_non_dev_parameters_rownames.txt"))

    tmp_report$derived_quants = read.table(file.path(html_dir,"html_derived_quants.txt"))
    row.names(tmp_report$derived_quants) = readLines(file.path(html_dir,"html_derived_quants_rownames.txt"))

    summary = fread(file.path(html_dir,"summary.csv"))
    tmp_report$log_det_hessian = summary$log_det_hessian[1]
    tmp_report$maximum_gradient_component = summary$mgc[1]
    tmp_report$inputs$covar = TRUE

    SS_plots(replist = tmp_report,dir=html_dir,html=FALSE,catchasnumbers = catchasnumbers)
    SS_html(tmp_report, filenotes = NULL, plotdir = file.path(html_dir,"plots/"), verbose = TRUE,openfile = FALSE)
    B = proc.time()
    time = round((B-A)[3],digits=2)
    if(verbose){
        print(paste0("Model html/plots created in ",time," seconds."))
    }
}
