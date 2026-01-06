

# Nicholas Ducharme-Barth
# 2024/02/15
# Set operating system

# Copyright (c) 2024 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

is_linux_os = function()
{
    system = Sys.info()["sysname"]
    if(system=="Linux"){
        return(TRUE)
    } else if(system=="Windows"){
        return(FALSE)
    } else {
        return(NA)
    }
}