pkgname <- "SADmoment"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('SADmoment')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("SADmoment-package")
### * SADmoment-package

flush(stderr()); flush(stdout())

### Name: SADmoment-package
### Title: Extrapolate SAD and calculate parameters of best fitting
###   distribution for subareas with several others
### Aliases: SADmoment-package SADmoment

### ** Examples

# Data from Barro Colorado Island (BCI) is used here for the example, 
# within total 50 ha area all the individuals 
# with at least 1cm (Diameter at Breast Height) (dbh) are considered.

SAD_Moment(dbci[,3],dbci[,4],dbci[,1],ext.rate=1.1,dplot=TRUE)

# In the following we used all the individuals 
# with at least 10cm (Diameter at Breast Height) (dbh) from the BCI data.
#R> name<-SAD_Moment(X,Y,Spnames,ext.rate=1.2,dplot=FALSE)
#R> summary(name)
#            Length Class     Mode   
#Data        6      -none-    list   
#Initial     4      -none-    list   
#Extrapolate 7      -none-    list   
#Time        5      proc_time numeric
#R> summary(name$Initial)
#               Length Class  Mode   
#Area            1     -none- numeric
#Moments        20     -none- numeric
#XY.Ratio        1     -none- numeric
#Number.of.bins  1     -none- numeric
#R> summary(name$Extrapolate)
#                                       Length Class  Mode   
#Area                                    1     -none- numeric
#Moments                                20     -none- numeric
#Number.of.bins                          1     -none- numeric
#Coefficients.of.Tchebychev.Polynomials 49     -none- numeric
#Tchebychev.Polynomials.Values          98     -none- numeric
#Tchebychev.Moments.Values               7     -none- numeric
#Extrapolate.SAD                        12     -none- numeric
#R> summary(name$Time)
#   user  system elapsed 
#  4.632   2.009   6.690 
#R>print(name$Data$SAD)
# [1] 27 11 20 30 29 38 30 21 20  5  4  1
#R> print(name$Extrapolate$Extrapolate.SAD)
# [1] 32 13 15 24 31 33 30 24 17 10  5  2
# 
#R> LLfun(X,Splist)
#Enter sub-area size: 9800
#     [,1]                                                                        
#[1,] "For the sub-area ( 9800 ) the best fitting distribution is : Poisson Gamma"
#[2,] "alpha= 1.2089269105249 and beta= 1.2089269105249"
#R> LLfun(X,Y,Splist)
#Enter sub-area size: 200000
#     [,1]                                                                      
#[1,] "For the sub-area ( 2e+05 ) the best fitting distribution is : Log-Normal"
#[2,] "mu= 3.51478121939837 and sigma= 2.2582385148143"
#



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
