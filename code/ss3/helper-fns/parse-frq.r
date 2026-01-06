    # define new S4 class so it doesn't break downstream code
    setClass("simpleFrq",
            representation(
                n_fisheries = "numeric",
                lf_range = "numeric",
                cateffpen = "data.frame",
                lnfrq = "data.frame",
                wtfrq = "data.frame"
            ),
            prototype=prototype(
                n_fisheries = numeric(),
                lf_range = numeric(),
                cateffpen = data.frame(),
                lnfrq = data.frame(),
                wtfrq = data.frame()
            ))
    
    simpleFrq = function() {return(new("simpleFrq"))}

    # define slot generics
    setGeneric("n_fisheries", function(x) standardGeneric("n_fisheries"))
    setGeneric("n_fisheries<-", function(x, value) standardGeneric("n_fisheries<-"))

    setGeneric("lf_range", function(x) standardGeneric("lf_range"))
    setGeneric("lf_range<-", function(x, value) standardGeneric("lf_range<-"))

    setGeneric("cateffpen", function(x) standardGeneric("cateffpen"))
    setGeneric("cateffpen<-", function(x, value) standardGeneric("cateffpen<-"))

    setGeneric("lnfrq", function(x) standardGeneric("lnfrq"))
    setGeneric("lnfrq<-", function(x, value) standardGeneric("lnfrq<-"))

    setGeneric("wtfrq", function(x) standardGeneric("wtfrq"))
    setGeneric("wtfrq<-", function(x, value) standardGeneric("wtfrq<-"))

    # define methods
    setMethod("n_fisheries", "simpleFrq", function(x) x@n_fisheries)
    setMethod("n_fisheries<-", "simpleFrq", function(x, value) {
    x@n_fisheries <- value
    x
    })

    setMethod("lf_range", "simpleFrq", function(x) x@lf_range)
    setMethod("lf_range<-", "simpleFrq", function(x, value) {
    x@lf_range <- value
    x
    })

    setMethod("cateffpen", "simpleFrq", function(x) x@cateffpen)
    setMethod("cateffpen<-", "simpleFrq", function(x, value) {
    x@cateffpen <- value
    x
    })

    setMethod("lnfrq", "simpleFrq", function(x) x@lnfrq)
    setMethod("lnfrq<-", "simpleFrq", function(x, value) {
    x@lnfrq <- value
    x
    })

    setMethod("wtfrq", "simpleFrq", function(x) x@wtfrq)
    setMethod("wtfrq<-", "simpleFrq", function(x, value) {
    x@wtfrq <- value
    x
    })

parse_frq = function(frq_path){
    require(data.table)
    require(magrittr)

    # read lines
        frq_lines = trimws(unname(sapply(readLines(frq_path),function(x)gsub("\t"," ",x))))
        frq_lines = trimws(unname(sapply(frq_lines,function(x)gsub("\\s+"," ",x))))

    # define new object
        outFrq = simpleFrq()

    # n_fisheries (numeric)
        pointer = grep("# Region Fisheries diffusion tag groups",frq_lines,fixed=TRUE)
        n_fisheries(outFrq) = as.numeric(strsplit(frq_lines[pointer+1],"\\s+")[[1]][2])
        
    # lf_range (named vector)
    # "Datasets"    "LFIntervals" "LFFirst"     "LFWidth"     "LFFactor" "WFIntervals" "WFFirst"     "WFWidth"     "WFFactor"
        pointer = grep("LFIntervals",frq_lines,fixed=TRUE)
        lf_range(outFrq) = as.numeric(strsplit(frq_lines[pointer+1],"\\s+")[[1]])
        names(lf_range(outFrq)) = c("Datasets","LFIntervals","LFFirst","LFWidth","LFFactor","WFIntervals","WFFirst","WFWidth","WFFactor")

    # extract fisheries instances and put in giant matrix
        pointer = grep("#YR MON WEEK FISH CATCH EFFORT CV",frq_lines,fixed=TRUE)
        if(length(pointer)==0){
            pointer = grep("# age_nage age_age1",frq_lines,fixed=TRUE) + 1
        }
        fil = lapply(frq_lines[(pointer+1):(pointer+lf_range(outFrq)["Datasets"])],function(x)as.numeric(strsplit(x,"\\s+")[[1]]))
        fil_length = sapply(fil,length)

        fi_mat = matrix(NA,nrow=length(fil_length),ncol=max(fil_length))
        for(i in seq_along(fil_length)){
            fi_mat[i,1:fil_length[i]] = fil[[i]]
        }

    # cateffpen (data.frame)
    # "year"    "month"   "week"    "fishery" "catch"   "effort"  "penalty"
        cateffpen(outFrq) = as.data.frame(fi_mat[,1:7])
        colnames(cateffpen(outFrq)) = c("year","month","week","fishery","catch","effort","penalty")

    # lnfrq (data.frame)
    # "year" "month" "week" "fishery" lfbins ...
    # 0 or greater in column 8
        lf_index = which(fi_mat[,8]>=0)
        mfcl_bin_lower = seq(from=lf_range(outFrq)[3],by=lf_range(outFrq)[4],length.out=lf_range(outFrq)[2])

        lnfrq(outFrq) = as.data.frame(fi_mat[lf_index,c(1:7,7+(1:lf_range(outFrq)[2]))])
        colnames(lnfrq(outFrq)) = c("year","month","week","fishery","catch","effort","penalty",mfcl_bin_lower)
        lnfrq(outFrq) = lnfrq(outFrq)[,-c(5:7)]

    # wtfrq (data.frame)
    # "year" "month" "week" "fishery" wfbins ...
    # if has values beyond 7 + LFIntervals and has lf then wf 
    # OR
    # column 8 is -1 and column 9 is 0 or greater
    wf_index_a = which((fi_mat[,9]>=0 & fi_mat[,8]== -1))
    wf_index_b = which(fi_mat[,8]>=0 & fi_mat[,7+lf_range(outFrq)[2]+1]>0)

    mfcl_wtbin_lower = seq(from=lf_range(outFrq)[7],by=lf_range(outFrq)[8],length.out=lf_range(outFrq)[6])

    wtfrq_a = as.data.frame(fi_mat[wf_index_a,c(1:7,8+(1:lf_range(outFrq)[6]))])
    wtfrq_b = as.data.frame(fi_mat[wf_index_b,c(1:7,7+lf_range(outFrq)[2]+(1:lf_range(outFrq)[6]))])

    wtfrq_df = rbind(wtfrq_a,wtfrq_b)
    colnames(wtfrq_df) = c("year","month","week","fishery","catch","effort","penalty",mfcl_wtbin_lower)
    wtfrq_df = wtfrq_df[,-c(5:7)]
    wtfrq_dt = as.data.table(wtfrq_df) %>%
                .[order(fishery,year,month)]
    wtfrq(outFrq) = as.data.frame(wtfrq_dt)

    return(outFrq)
}
