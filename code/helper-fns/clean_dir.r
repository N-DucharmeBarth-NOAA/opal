
clean_dir = function(to_dir,keep_files=c("ss.par","control.ss","starter.ss_new","forecast.ss_new","control.ss_new","data_echo.ss_new","data.ss_new","ss3.par","warning.sso","executable.txt")){
    all_files = list.files(to_dir)
    clean_files = all_files[which(!(all_files %in% c(all_files[grep("html_",all_files,fixed=TRUE)],all_files[grep(".csv",all_files,fixed=TRUE)],keep_files)))]
    print(paste0(length(clean_files)," files removed taking up ",round(sum(unname(sapply(file.path(to_dir,clean_files),file.size)))/1024^2,digits=2)," MB."))
    unname(sapply(file.path(to_dir,clean_files),unlink))
    print(paste0("Directory size is now ",round(sum(unname(sapply(file.path(to_dir,list.files(to_dir)),file.size)))/1024^2,digits=2)," MB."))
}

