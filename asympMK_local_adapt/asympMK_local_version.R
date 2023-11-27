#
#
#		asymptoticMK: Asymptotic McDonald-Kreitman Test web service
#
#		By Benjamin C. Haller and Philipp W. Messer
#		Copyright (C) 2017 Philipp Messer.
#
#		This web service should be available at http://benhaller.com/messerlab/asymptoticMK.html
#		The Github repository for asymptoticMK is at https://github.com/MesserLab/asymptoticMK
#
#

#  This file is part of asymptoticMK.
#
#  asymptoticMK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
#  asymptoticMK is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License along with asymptoticMK.  If not, see <http://www.gnu.org/licenses/>.


#######################################################################################
#
#	To run asymptoticMK locally, first source all the functions between this comment
#	and the next comment, below, that is in the same box style.  Then, read on from
#	that second comment.  Everything between the two box comments may be considered a
#	black box; you should not need to understand it or modify it to use asymptoticMK.
#
#######################################################################################

# Get a CI using Monte Carlo simulation based upon a fitted model.  This is necessary because
# getting confidence intervals for non-linear models is a complicated business, apparently.
# Thanks to Andrej-Nikolai Spiess (http://www.dr-spiess.de) for this code.
# See: https://www.r-bloggers.com/predictnls-part-1-monte-carlo-simulation-confidence-intervals-for-nls-models/
# Or, if that link goes stale: http://stats.stackexchange.com/a/251501/141766

######################################################################################

#### asympMK customization - MJ27112023

#### this code is a customized version of the asympMK script that fits asymptote locally
#### the original code does not fit an asymptote in the case of a exponential fit 
#### as exponential fit is still treated as a linear fit thereby giving elevated alpha results

#### !! run the following two functions as they are (predictNLS, fitMKmodel) !! ####
predictNLS <- function(
    object, 
    newdata,
    level = 0.95, 
    nsim = 10000,
    ...
){
  require(MASS, quietly = TRUE)
  
  ## get right-hand side of formula
  RHS <- as.list(object$call$formula)[[3]]
  EXPR <- as.expression(RHS)
  
  ## all variables in model
  VARS <- all.vars(EXPR)
  
  ## coefficients
  COEF <- coef(object)
  
  ## extract predictor variable    
  predNAME <- setdiff(VARS, names(COEF))  
  
  ## take fitted values, if 'newdata' is missing
  if (missing(newdata)){
    newdata <- eval(object$data)[predNAME]
    colnames(newdata) <- predNAME
  }
  
  ## check that 'newdata' has same name as predVAR
  if (names(newdata)[1] != predNAME) stop("newdata should have name '", predNAME, "'!")
  
  ## get parameter coefficients
  COEF <- coef(object)
  
  ## get variance-covariance matrix
  VCOV <- vcov(object)
  
  ## augment variance-covariance matrix for 'mvrnorm' 
  ## by adding a column/row for 'error in x'
  NCOL <- ncol(VCOV)
  ADD1 <- c(rep(0, NCOL))
  ADD1 <- matrix(ADD1, ncol = 1)
  colnames(ADD1) <- predNAME
  VCOV <- cbind(VCOV, ADD1)
  ADD2 <- c(rep(0, NCOL + 1))
  ADD2 <- matrix(ADD2, nrow = 1)
  rownames(ADD2) <- predNAME
  VCOV <- rbind(VCOV, ADD2) 
  
  ## iterate over all entries in 'newdata' as in usual 'predict.' functions
  NR <- nrow(newdata)
  respVEC <- numeric(NR)
  seVEC <- numeric(NR)
  varPLACE <- ncol(VCOV)   
  
  ## define counter function
  counter <- function (i){
    if (i%%10 == 0) 
      cat(i)
    else cat(".")
    if (i%%50 == 0) 
      cat("\n")
    flush.console()
  }
  
  outMAT <- NULL 
  
  for (i in 1:NR){
    #counter(i)		# show a counter for lengthy fits; commented out to reduce noise here...
    
    ## get predictor values and optional errors
    predVAL <- newdata[i, 1]
    if (ncol(newdata) == 2) predERROR <- newdata[i, 2] else predERROR <- 0
    names(predVAL) <- predNAME  
    names(predERROR) <- predNAME  
    
    ## create mean vector for 'mvrnorm'
    MU <- c(COEF, predVAL)
    
    ## create variance-covariance matrix for 'mvrnorm'
    ## by putting error^2 in lower-right position of VCOV
    newVCOV <- VCOV
    newVCOV[varPLACE, varPLACE] <- predERROR^2
    
    ## create MC simulation matrix
    simMAT <- mvrnorm(n = nsim, mu = MU, Sigma = newVCOV, empirical = TRUE)
    
    ## evaluate expression on rows of simMAT
    EVAL <- try(eval(EXPR, envir = as.data.frame(simMAT)), silent = TRUE)
    if (inherits(EVAL, "try-error")) stop("There was an error evaluating the simulations!")
    
    ## collect statistics
    PRED <- data.frame(predVAL)
    colnames(PRED) <- predNAME   
    FITTED <- predict(object, newdata = data.frame(PRED))
    MEAN.sim <- mean(EVAL, na.rm = TRUE)
    SD.sim <- sd(EVAL, na.rm = TRUE)
    MEDIAN.sim <- median(EVAL, na.rm = TRUE)
    MAD.sim <- mad(EVAL, na.rm = TRUE)
    QUANT <- quantile(EVAL, c((1 - level)/2, level + (1 - level)/2))
    RES <- c(FITTED, MEAN.sim, SD.sim, MEDIAN.sim, MAD.sim, QUANT[1], QUANT[2])
    outMAT <- rbind(outMAT, RES)
  }
  
  colnames(outMAT) <- c("fit", "mean", "sd", "median", "mad", names(QUANT[1]), names(QUANT[2]))
  rownames(outMAT) <- NULL
  
  #cat("\n")	# commented out along with the call to counter() above
  
  return(outMAT)  
}

# core code for two-step nls2() model fit at a given level of precision (res)
fitMKmodel <- function(alpha_trimmed, f_trimmed, res)
{
  require(nls2)
  
  mod <- tryCatch({
    st <- expand.grid(const_a=seq(-1,1,length.out=res + 1), const_b=seq(-1,1,length.out=res), const_c=seq(1,10,length.out=res + 1))
    nls2(alpha_trimmed ~ const_a + const_b * exp(-const_c* f_trimmed), start=st, algorithm="brute-force", control=nls.control(maxiter=NROW(st)))
  },
  error=function(cond) {})
  
  if (length(mod) == 0)
    return(NULL)
  
  mod2 <- tryCatch({
    nls2(alpha_trimmed ~ const_a + const_b * exp(-const_c* f_trimmed), start = mod, control=nls.control(maxiter=200))
  },
  error=function(cond) {})
  
  if (length(mod2) == 0)
    return(NULL)
  
  return(mod2)
}

###################################################################

#### ! declaring the constants here ! ####

#inputs for thaliana

#sfs containing table
df = read.delim("input_files/zam_poly_wga.txt", sep = "\t")

#declaring the divergence constants
d = 223926
d0 = 518897

#delcaring the trimming points for the analysis
xlow = 0.02
xhigh = 0.98

#### !run everything in between this and the next comment as is! ####

require(nls2)
require(MASS)

if (is.na(d0) || is.null(d0)){stop("Malformed d0 (must be numeric).")}
if (is.na(d) || is.null(d)){stop("Malformed d (must be numeric).")}
if (is.na(xlow) || is.null(xlow)){stop("Malformed xlow (must be numeric).")}
if (is.na(xhigh) || is.null(xhigh)){stop("Malformed xhigh (must be numeric).")}

#	Bounds-check response variables
#
if (d0 <= 0){stop("d0 must greater than zero.")}
if (d <= 0){stop("d must greater than zero.")}
if ((xlow < 0.0) || (xlow > 1.0)){stop("xlow must be in the interval [0,1].")}
if ((xhigh < 0.0) || (xhigh > 1.0)){stop("xhigh must be in the interval [0,1].")}
if (xlow >= xhigh){stop("xlow must be less than xhigh.")}

#	Read in the file and check its format
#
if (NCOL(df) != 3){stop("Dataframe df does not contain exactly three tab-separated columns.")}
if (NROW(df) <= 0){stop("Dataframe df contains no data rows.")}

cols <- names(df)

suppressWarnings(	# the goal is to generate NAs here, so we don't want to see the warnings...
  if (!is.na(as.numeric(cols[1])) || !is.na(as.numeric(cols[2])) || !is.na(as.numeric(cols[3]))){
    stop("Dataframe df has a numeric column name; probably the required header row is missing.")
})

f <- df[[1]]
p <- df[[2]]
p0 <- df[[3]]

if (!is.numeric(f))
  stop("The first column of the dataframe df, frequency, is not numeric.")
if (!is.numeric(p))
  stop("The second column of the dataframe df, p, is not numeric.")
if (!is.numeric(p0))
  stop("The third column of the dataframe df, p0, is not numeric.")
if (any(is.na(f)))
  stop("The first column of the dataframe df, frequency, contains NA values (not allowed).")
if (any(is.na(p)))
  stop("The second column of the dataframe df, p, contains NA values (not allowed).")
if (any(is.na(p0)))
  stop("The third column of the dataframe df, p0, contains NA values (not allowed).")
if (any(f < 0.0) || any(f > 1.0))
  stop("The first column of the dataframe df, frequency, contains values out of the required range [0,1].")
if (any(p < 0))		# note that zero is allowed, although not recommended
  stop("The second column of the dataframe df, p, contains values < 0 (not allowed).")
if (all(p == 0))		# not all can be zero, however
  stop("The second column of the dataframe df, p, contains all values == 0 (not allowed).")
if (any(p0 <= 0))
  stop("The third column of the dataframe df, p0, contains values <= 0 (not allowed).")

if (NROW(df) < 3)
  stop("At least three data rows are required, to constrain the fit.")

#	Compute alpha values and trim
#
alpha <- 1 - (d0/d) * (p/p0)
cutoff_f1 <- xlow
cutoff_f2 <- xhigh

trim <- ((f >= cutoff_f1) & (f <= cutoff_f2))

if (sum(trim) < 3)
  stop("At least three data points are required after trimming the frequency range, to constrain the fit.")

f_trimmed <- f[trim]
alpha_trimmed <- alpha[trim]

#	Compute the original McDonald-Kreitman alpha; we decided to use the trimmed data for this.
#
alpha_nonasymp <- 1 - (d0/d) * (sum(p[trim])/sum(p0[trim]))			# from trimmed data
#alpha_nonasymp <- 1 - (d0/d) * (sum(p)/sum(p0))						# from untrimmed data

#	Fit models
#
mod1 <- fitMKmodel(alpha_trimmed, f_trimmed, 10)

if (length(mod1) == 0){
  # try a deeper scan for a decent fit
  mod1 <- fitMKmodel(alpha_trimmed, f_trimmed, 20)
}

tryCatch({
  mod2 <- lm(alpha_trimmed ~ f_trimmed)
},
error=function(cond) {})

linear_better <- FALSE

if ((length(mod1) == 0) || (AIC(mod2) < AIC(mod1))){
  
  #edit MJ-27112023 - changing linear_better to FALSE to always fit the exponential model
  linear_better <- FALSE
}

if (!linear_better){
  # if we're leaning toward the exponential model, check for ridiculously wide confidence intervals; sometimes
  # we should reject the exponential model for that reason, because it is basically just a linear model with
  # a "cheat" of a swing up or down to fit one additional data point perfectly, which is lame :->
  ci_pred <- predictNLS(mod1, newdata=data.frame(f_trimmed=1.0))
  alpha_1_low <- ci_pred[6]
  alpha_1_high <- ci_pred[7]
  
  if ((alpha_1_low < -100) || (alpha_1_high > 100)){linear_better <- TRUE}
}

# Prepare for output and plotting
full_seq <- seq(from=min(f), to=max(f), by=0.001)
trimmed_seq <- seq(from=min(f_trimmed), to=max(f_trimmed), by=0.001)

if (linear_better){
  alpha_1_est <- predict(mod2, newdata=data.frame(f_trimmed=1.0))
  ci_pred <- predict(mod2, newdata=data.frame(f_trimmed=1.0), interval="confidence")	# we want confidence, not prediction
  alpha_1_low <- ci_pred[2]
  alpha_1_high <- ci_pred[3]
  const_a <- coef(mod2)["(Intercept)"]
  const_b <- coef(mod2)["f_trimmed"]
  const_c <- NA
  
  full_predicts <- predict(mod2, newdata=data.frame(f_trimmed=full_seq))
  trimmed_predicts <- predict(mod2, newdata=data.frame(f_trimmed=trimmed_seq))
  fit_color <- "red"
}else{
  alpha_1_est <- predict(mod1, newdata=data.frame(f_trimmed=1.0))
  const_a <- coef(mod1)["const_a"]
  const_b <- coef(mod1)["const_b"]
  const_c <- coef(mod1)["const_c"]
  
  full_predicts <- predict(mod1, newdata=data.frame(f_trimmed=full_seq))
  trimmed_predicts <- predict(mod1, newdata=data.frame(f_trimmed=trimmed_seq))
  fit_color <- "red"
}


#	BEGIN OUTPUT
#
result_df <- data.frame(model=(if ((length(mod1) == 0) || linear_better) "linear" else "exponential"), a=const_a, b=const_b, c=const_c, alpha_asymptotic=alpha_1_est, CI_low=alpha_1_low, CI_high=alpha_1_high, alpha_original=alpha_nonasymp, row.names=NULL)


#### ! the result_df contains the relevant information for making the asymptote plot !####

#################################################

#clearing the environment and keeping only the result_df and alpha_trimmed - contains the class-specific alpha
rm(list = setdiff(ls(), c("result_df", "alpha_trimmed", "dmel_fit", "aratha_fit")))

#declaring the SFS classes here
sfs_class = seq(0.02, 0.98, by = 0.01)

#calculating the asymptote_fit
asympt_fit = lapply(sfs_class, function(x) {result_df$a + result_df$b*exp(-result_df$c*x)})
asympt_fit = unlist(asympt_fit)

#storing the information in a dataframe
alpha_fit = data.frame(alpha = alpha_trimmed, asympt_fit = asympt_fit, sfs_class = sfs_class)

# #removing outputs here for per species if required
# aratha_fit = alpha_fit
# hsap_fit = alpha_fit
# dmel_fit = alpha_fit

#importing ggplot2
library(ggplot2)

#### plotting a single asymptote ####

#plotting the asymptote fit
asymp_fit = ggplot(alpha_fit)+geom_line(aes(x = sfs_class, y = asympt_fit))+geom_point(aes(x = sfs_class, y = alpha), size = 2)+
  geom_hline(yintercept = result_df$alpha_asymptotic) +
  scale_x_continuous(breaks = seq(0,1,by=0.1))+
  geom_hline(yintercept = result_df$CI_low, colour = "red", size = 0.5, linetype = "dashed")+
  geom_hline(yintercept = result_df$CI_high, colour = "red", size = 0.5, linetype = "dashed")+
  geom_vline(xintercept = 0.1, colour = "black", size = 0.5)+ geom_vline(xintercept = 0.9, colour = "black", size = 0.5)+
  theme_linedraw()+ xlab("SFS classes")+ylab("alpha")+theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
                                          axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18))

#saving to a file
ggsave(x = asymp_fit, filename = "plots/asympMK_fit_dmel.tiff", device = "tiff", width = 10, height = 10)

####plotting multiple asymptotes ####
# 
# #merging multiple output files
# allspecies_asymp = data.frame(asymp_fit = c(aratha_fit$asympt_fit, dmel_fit$asympt_fit, hsap_fit$asympt_fit),
#                               species = c(rep("aratha", 97), rep("dmel", 97), rep("hsap", 97)),
#                               sfs_class = c(aratha_fit$sfs_class, dmel_fit$sfs_class, hsap_fit$sfs_class))
# 
# #plotting
# ggplot(data = allspecies_asymp)+geom_line(aes(x = sfs_class, y = asymp_fit, colour = species), size = 2)+
#   theme_linedraw()

