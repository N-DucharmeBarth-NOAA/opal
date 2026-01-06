

ss_transfer = function(from_dir,to_dir,executable_dir_stem,linux_ws=is_linux_os())
{
    #________________________________________________________________________________________________________________________________________________________________________________________________________
    # transfer files
        ss_files_to_transfer = c("control.ss_new","data_echo.ss_new","forecast.ss_new","starter.ss_new","ss3.par")
        if(linux_ws)
        {
            exec_version = "3.30.22.1" # only because linux version of 3.30.23.1 does not work on codespaces
            exec_name = "ss3_linux"
        } else{
            exec_version = "3.30.23.1"
            exec_name = "ss3_win.exe"
        }

        dir.create(to_dir,recursive=TRUE)
        exec = file.path(executable_dir_stem,exec_version,exec_name)
        exec_save = file.path(exec_version,exec_name)

        # transfer files
        file.copy(from=file.path(from_dir,ss_files_to_transfer), to=file.path(to_dir,ss_files_to_transfer), overwrite = TRUE)
        file.copy(from=exec, to=to_dir, overwrite = TRUE)
        writeLines(exec_save,con=file.path(to_dir,"executable.txt")) 

    # rename
        file.rename(from=file.path(to_dir,ss_files_to_transfer),to=file.path(to_dir,gsub("_echo","",gsub("_new","",ss_files_to_transfer))))
        if(linux_ws){
            file.rename(from=file.path(to_dir,exec_name),to=file.path(to_dir,"ss3"))
            system(paste0("cd ",to_dir,"; chmod 777 ss3"))
        } else {
            file.rename(from=file.path(to_dir,exec_name),to=file.path(to_dir,"ss3.exe"))
        }
}
