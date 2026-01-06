

# Nicholas Ducharme-Barth
# 2024/01/31
# Make SS output

# Copyright (c) 2024 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

# dir = to_dir
# starter="starter.ss_new"
# control="control.ss_new"
# data="data.ss_new"
# forecast="forecast.ss_new"
# par="ss3.par"
# exe="ss3"
# calc_hessian=FALSE
# files_to_transfer="all"
# linux_ws=TRUE

make_ss_output = function(dir,
                          starter="starter.ss_new",
                          control="control.ss_new",
                          data="data_echo.ss_new",
                          forecast="forecast.ss_new",
                          par="ss3.par",
                          files_to_transfer="all",
                          executable_dir_stem = file.path(this.path::this.proj(),"executables","stock-synthesis"),
                          calc_hessian=FALSE,
                          linux_ws=is_linux_os(),
                          verbose=TRUE)
{
    A = proc.time()
    # create tmp dir
        tmp_dir = file.path(dir,"tmp/")
        dir.create(tmp_dir,recursive=TRUE)
    # copy files
        file.copy(from=file.path(dir,c(starter,control,data,forecast,par)),to=file.path(tmp_dir,c("starter.ss","control.ss","data.ss","forecast.ss",par)))
    # read the executable.txt file & copy
    # add check to make sure copied executable is compatible with linux_ws
        exe_path = readLines(file.path(dir,"executable.txt"))
        exec_name = tail(strsplit(exe_path,"/")[[1]],n=1)
        if(linux_ws & length(grep("win.exe",exec_name))==1){
            exe_path = gsub("win.exe","linux",exe_path)
            exec_name = gsub("win.exe","linux",exec_name)
        }
        if(!linux_ws & length(grep("linux",exec_name))==1){
            exe_path = gsub("linux","win.exe",exe_path)
            exec_name = gsub("linux","win.exe",exec_name)
        }
        file.copy(from=file.path(executable_dir_stem,exe_path), to=tmp_dir, overwrite = TRUE)

    # give permissions
        if(linux_ws){
            file.rename(from=file.path(tmp_dir,exec_name),to=file.path(tmp_dir,"ss3"))
            system(paste0("cd ",tmp_dir,"; chmod 777 ss3"))
            exe = "ss3"
        } else{
            file.rename(from=file.path(tmp_dir,exec_name),to=file.path(tmp_dir,"ss3.exe"))
            exe = "ss3.exe"
        }

    # get initial file names
        original_files = unique(c(starter,control,data,forecast,par,exe,list.files(tmp_dir)))
    
    # modify starter
        tmp_starter = SS_readstarter(file=file.path(tmp_dir,"starter.ss"),verbose=FALSE)
        tmp_starter$init_values_src = 1
        SS_writestarter(tmp_starter,dir=tmp_dir,file="starter.ss",overwrite=TRUE,verbose=FALSE)
        
    # run ss
    if(calc_hessian){
        run(dir=tmp_dir,exe=exe,extras="-maxfn 0 -phase 9999",show_in_console=FALSE,skipfinished =FALSE,verbose=TRUE)
    } else {
        run(dir=tmp_dir,exe=exe,extras="-maxfn 0 -phase 9999 -nohess",show_in_console=FALSE,skipfinished =FALSE,verbose = TRUE)
    }

    # transfer files to original dir
        new_files = list.files(tmp_dir)
        dir_files = list.files(dir)
        if(length(files_to_transfer)==1&files_to_transfer=="all")
        {
            files_out = NA
            for(i in 1:length(dir_files)){
                if(!(new_files[i] %in% dir_files)){
                    files_out = c(files_out,new_files[i])
                }
            }
            files_out = as.vector(na.omit(files_out))
        } else {
            files_out = files_to_transfer
        }
        file.copy(from=file.path(tmp_dir,files_out),to=dir)

    # clean-up
        unlink(tmp_dir, recursive = TRUE)
    B=proc.time()
    time = round((B-A)[3],digits=2)
    if(verbose){
        print(paste0("Model output created in ",time," seconds."))
    }
}
