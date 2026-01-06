

# Nicholas Ducharme-Barth
# 2024/02/15
# Summarize model fleet structure
# Number
# Name
# Type
# Units
# Selex
# Mirror
# Catch start
# Catch end
# Average catch n
# Average catch mt
# index start
# index end
# Len comp start
# Len comp end
# Len comp total N
# Weight comp start
# Weight comp end
# Weight comp total N

# Copyright (c) 2024 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

# tmp_report = SS_output(dir,verbose = FALSE,printstats = FALSE)
summarize_fleets = function(dir,tmp_report)
{
    # read ctl
    tmp_ctl = SS_readctl(file=file.path(dir,"control.ss"),datlist = file.path(dir,"data.ss"))
    # read data
    tmp_data = SS_readdat(file=file.path(dir,"data.ss"))

    tmp_catch = as.data.table(tmp_data$catch) %>%
                .[catch>0,.(catch_start=min(year),catch_end=max(year)),by=fleet]
    
    tmp_avg_catch = as.data.table(tmp_report$catch) %>%
                    .[,.(catch_kn=mean(sel_num),catch_mt=mean(sel_bio)),by=Fleet] %>%
                    setnames(.,"Fleet","fleet")
    
    tmp_cpue = as.data.table(tmp_data$CPUE) %>%
                .[,.(index_start=min(year),index_end=max(year)),by=index] %>%
                setnames(.,"index","fleet")
    
    tmp_lencomp = as.data.table(tmp_data$lencomp) %>%
                  .[,.(lencomp_start=min(year),lencomp_end=max(year),lencomp_nsamp=sum(Nsamp)),by=fleet]

    # check if any size comp in terms of length
    if(length(which(tmp_data$scale_per_method >= 3))>0)
    {
        tmp_method = which(tmp_data$scale_per_method >= 3)
        tmp_method_dt.list = as.list(rep(NA,length(tmp_method)))
        for(i in 1:length(tmp_method))
        {
            tmp_method_dt.list[[i]] = as.data.table(tmp_data$sizefreq_data_list[[tmp_method[i]]]) %>%
                                    .[,.(lencomp_start=min(year),lencomp_end=max(year),lencomp_nsamp=sum(Nsamp)),by=fleet] 
        }   
        tmp_method_dt = rbindlist(tmp_method_dt.list)
        tmp_lencomp = rbind(tmp_lencomp,tmp_method_dt) %>%
                      .[,.(lencomp_start=min(lencomp_start),lencomp_end=max(lencomp_end),lencomp_nsamp=sum(lencomp_nsamp)),by=fleet]
    }

    if(length(which(tmp_data$scale_per_method < 3))>0)
    {
        tmp_method = which(tmp_data$scale_per_method < 3)
        tmp_method_dt.list = as.list(rep(NA,length(tmp_method)))
        for(i in 1:length(tmp_method))
        {
            tmp_method_dt.list[[i]] = as.data.table(tmp_data$sizefreq_data_list[[tmp_method[i]]]) %>%
                                    .[,.(wcomp_start=min(year),wcomp_end=max(year),wcomp_nsamp=sum(Nsamp)),by=fleet] 
        }   
        tmp_wcomp = rbindlist(tmp_method_dt.list)
    } else {
        tmp_wcomp = data.table(fleet=1:nrow(tmp_data$fleetinfo)) %>%
                    .[,wcomp_start:=NA] %>%
                    .[,wcomp_end:=NA] %>%
                    .[,wcomp_nsamp:=NA]
    }

    # make table
    dt = as.data.table(tmp_data$fleetinfo) %>%
        .[,fleet:=1:.N] %>%
        .[,.(fleetname,fleet,type,area,units)] %>%
        .[,selex:=tmp_ctl$size_selex_types$Pattern] %>%
        .[,mirror:=tmp_ctl$size_selex_types$Special] %>%
        merge(.,tmp_catch,by="fleet",all.x=TRUE) %>%
        merge(.,tmp_avg_catch,by="fleet",all.x=TRUE) %>%
        merge(.,tmp_cpue,by="fleet",all.x=TRUE) %>%
        merge(.,tmp_lencomp,by="fleet",all.x=TRUE) %>%
        merge(.,tmp_wcomp,by="fleet",all.x=TRUE)

    # write table
    fwrite(dt,file=file.path(dir,"fleet_summary.csv"))
}
