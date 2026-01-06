

# Nicholas Ducharme-Barth
# 2025/01/09
# Helper functions to define selex pars

# Copyright (c) 2024 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

make_size_selex_par_24 = function(peak,
top_logit,
ascend_se,
descend_se,
start_logit,
end_logit)
{
    return(data.table("LO"=c(peak[1],top_logit[1],ascend_se[1],descend_se[1],start_logit[1],end_logit[1]),
    "HI"=c(peak[2],top_logit[2],ascend_se[2],descend_se[2],start_logit[2],end_logit[2]),
    "INIT"=c(peak[3],top_logit[3],ascend_se[3],descend_se[3],start_logit[3],end_logit[3]),
    "PRIOR"=c(peak[3],top_logit[3],ascend_se[3],descend_se[3],start_logit[3],end_logit[3]),
    "PR_SD"=rep(99,6),
    "PR_type"=rep(0,6),
    "PHASE"=c(peak[4],top_logit[4],ascend_se[4],descend_se[4],start_logit[4],end_logit[4]),
    "env_var&link"=rep(0,6),
    "dev_link"=rep(0,6),
    "dev_minyr"=rep(0,6),
    "dev_maxyr"=rep(0,6),
    "dev_PH"=rep(0,6),
    "Block"=rep(0,6),
    "Block_Fxn"=rep(0,6)))
}


make_size_selex_par_5 = function()
{
    return(data.table("LO"=rep(-99,2),
    "HI"=rep(10,2),
    "INIT"=rep(-1,2),
    "PRIOR"=rep(1,2),
    "PR_SD"=rep(99,2),
    "PR_type"=rep(0,2),
    "PHASE"=rep(-99,2),
    "env_var&link"=rep(0,2),
    "dev_link"=rep(0,2),
    "dev_minyr"=rep(0,2),
    "dev_maxyr"=rep(0,2),
    "dev_PH"=rep(0,2),
    "Block"=rep(0,2),
    "Block_Fxn"=rep(0,2)))
}

make_age_selex_par_11 = function(max_age=10)
{
    return(data.table("LO"=c(0,10),
    "HI"=c(10,100),
    "INIT"=c(0,max_age),
    "PRIOR"=rep(0,2),
    "PR_SD"=rep(99,2),
    "PR_type"=rep(0,2),
    "PHASE"=rep(-99,2),
    "env_var&link"=rep(0,2),
    "dev_link"=rep(0,2),
    "dev_minyr"=rep(0,2),
    "dev_maxyr"=rep(0,2),
    "dev_PH"=rep(0,2),
    "Block"=rep(0,2),
    "Block_Fxn"=rep(0,2)))
}

make_size_selex_par_1 = function(inflection,
width)
{
    return(data.table("LO"=c(inflection[1],width[1]),
    "HI"=c(inflection[2],width[2]),
    "INIT"=c(inflection[3],width[3]),
    "PRIOR"=c(inflection[3],width[3]),
    "PR_SD"=rep(99,2),
    "PR_type"=rep(0,2),
    "PHASE"=c(inflection[4],width[4]),
    "env_var&link"=rep(0,2),
    "dev_link"=rep(0,2),
    "dev_minyr"=rep(0,2),
    "dev_maxyr"=rep(0,2),
    "dev_PH"=rep(0,2),
    "Block"=rep(0,2),
    "Block_Fxn"=rep(0,2)))
}
