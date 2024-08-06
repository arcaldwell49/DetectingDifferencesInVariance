library(simstudy)
library(skedastic)
library(car)
library(metafor)
gen_id <- defData(varname = "pre", dist = "normal", formula = 100, variance = 20^2,
                      id = "id")
gen_id <- defData(gen_id, varname = "t0", dist = "normal", formula = 0,
                  variance = 1)
gen_id <- defData(gen_id, varname = "rx", dist = "trtAssign",
               formula = "1;1")
gen_id <- defData(gen_id, varname = "post", dist = "normal",
                  formula = "pre - (5+t0) * rx",
                  variance = 5)
set.seed(282721)

dt_id <- genData(10000, gen_id)
dt_id$delta = dt_id$post - dt_id$pre




# Functions -----
p_from_z = function(x,
                    alternative = "two.sided"){
  if(alternative == "two.sided"){
    2*pnorm(-abs(unlist(x)))
  } else  if (alternative  == "greater"){
    pnorm(x, lower.tail = FALSE)
  } else if (alternative  == "less"){
    pnorm(x, lower.tail = TRUE)
  } else{
    stop("alternative must be two.sided, greater, or less")
  }
  
}

## tests ------
# F-test for variances # var.test()
# barlett # bartlett.test()
# levene #  car::levene
# fligner-killeen # fligner.test()
# Breusch-Pagan-Koenker # skedastic::breusch_pagan()
# glejser (model based) # skedastic::glejser()
glejser.test <- function(model, Z_var, dataset){
  # model is the lmer object
  # Z_var is the explanatory variable to regress the abs values on
  # dataset is the full data
  
  ares = auxresponse = abs(residuals(model))
  n = length(ares)
  Z = as.matrix(dataset[Z_var])
  if(length(Z_var == 1)){
    Z <- cbind(1, Z)
  }
  q = ncol(Z)-1
  auxres <- stats::lm.fit(Z, ares)$residuals
  #sum_lm <- summary(lm(as.formula(paste0("ares ~ ", Z_var)), data = dataset))
  sigma_hatsq <- sum(residuals(model)^2)/n
  teststat <- (sum(auxresponse^2) - n * mean(auxresponse)^2 - 
                 sum(auxres^2))/(sigma_hatsq * (1 - 2/pi))
  pval <- stats::pchisq(teststat, df = q, lower.tail = FALSE)
  
  statistic = teststat
  names(statistic) = "Chi-squared"
  null1 = 0
  names(null1) = "heteroskedasticity"
  
  rval <- list(statistic = statistic,
               parameter = q,
               p.value = as.numeric(pval),
               null.value = null1,
               alternative = "greater",
               method = "Glejser Test for Heteroskedasticity",
               data.name = as.character(model$call)[2])
  class(rval) <- "htest"
  return(rval)
  
}
# difference in variance (normal approximation test)
diff_var_test = function(sd1, n1,
                         sd2, n2,
                         alternative = c("two.sided",
                                         "greater",
                                         "less")){
  alternative = match.arg(alternative)
  var1 = sd1^2
  var1_se <- var1*(sqrt(2/(n1 - 1)))
  var2 = sd2^2
  var2_se <- var2*(sqrt(2/(n2 - 1)))
  
  est_diff <- var1-var2
  
  est_diff_SE <- sqrt(var1_se^2 + var2_se^2) 
  teststat <- est_diff/est_diff_SE
  pval = p_from_z(teststat, alternative = alternative)
  
  statistic = teststat
  names(statistic) = "z"
  null1 = 0
  names(null1) = "difference in variances"
  estimate = est_diff
  names(estimate)= "differnce in variance"
  std_err = sqrt(est_diff_SE)
  sum_stats = paste0("SD1 = ", sd1, ", SD2 = ", sd2)
  
  rval <- list(statistic = statistic,
               p.value = as.numeric(pval),
               stderr = std_err,
               null.value = null1,
               alternative = alternative,
               method = "Difference in Variances",
               data.name = sum_stats)
  class(rval) <- "htest"
  return(rval)
}
# SDir (normal based test from Hopkins variances estimate)
diff_sdir_test = function(sd1, n1,
                         sd2, n2,
                         alternative = c("two.sided",
                                         "greater",
                                         "less")){

  alternative = match.arg(alternative)
  var1 = sd1^2
  var1_se <- var1*(sqrt(2/(n1 - 1)))
  var2 = sd2^2
  var2_se <- var2*(sqrt(2/(n2 - 1)))
  df1 = n1 -1
  df2 = n2 - 1
  if((var1-var2)< 0){
    sd_ir = sqrt(var1-var2)
    
  } else{
    sd_ir = -1*sqrt(var2-var1)
    
  }
  est_diff <- var1-var2
  
  est_diff_SE <- sqrt(2 * (sd1^4/df1 + sd2^4/df2)) 
  teststat <- est_diff/est_diff_SE
  pval = p_from_z(teststat, alternative = alternative)
  
  statistic = teststat
  names(statistic) = "z"
  null1 = 0
  names(null1) = "Standard Deviation of the Individual Response"
  estimate = sd_ir
  names(estimate)= "Standard Deviation of the Individual Response"
  std_err = sqrt(est_diff_SE)
  sum_stats = paste0("SD1 = ", sd1, ", SD2 = ", sd2)
  
  rval <- list(statistic = statistic,
               p.value = as.numeric(pval),
               stderr = std_err,
               null.value = null1,
               alternative = alternative,
               method = "Difference in Variances",
               data.name = sum_stats)
  class(rval) <- "htest"
  return(rval)
}
# ratio of variances (normal based test?) logVR

log_vr_test = function(sd1, n1,
                       sd2, n2,
                       alternative = c("two.sided",
                                       "greater",
                                       "less")){
  
  # old code
  # var1 = sd1^2
  # var1_se <- var1*(sqrt(2/(n1 - 1)))
  # var2 = sd2^2
  # var2_se <- var2*(sqrt(2/(n2 - 1)))
  # yi <- log(sqrt(var1)/sqrt(var2)) + 1/(2 * (n1 - 1)) - 1/(2 * 
  #                                                                     (n2 - 1))
  # 
  # vi <- 1/(2 * (n1 - 1)) + 1/(2 * (n2 - 1))
  # 
  # teststat <- yi/sqrt(vi)
  alternative = match.arg(alternative)
  # metafor code
  es_est = metafor::escalc(
  measure = "VR",
  n1i = n1,
  n2i = n2,
  sd1i = sd1,
  sd2i = sd2
  )

teststat <- es_est$yi/sqrt(es_est$vi)
pval = p_from_z(teststat, 
                alternative = alternative)

statistic = teststat
names(statistic) = "z"
null1 = 0
names(null1) = "log ratio of variances"
estimate = es_est$yi
names(estimate)= "log VR"
std_err = sqrt(es_est$vi)
sum_stats = paste0("SD1 = ", sd1, ", SD2 = ", sd2)

rval <- list(statistic = statistic,
             p.value = as.numeric(pval),
             stderr = std_err,
             null.value = null1,
             alternative = alternative,
             method = "log Ratio of Variances",
             data.name = sum_stats)
class(rval) <- "htest"
return(rval)
}

# coefficient of variation approach
log_cvr_test = function(sd1, n1, mean1,
                       sd2, n2, mean2,
                       
                       alternative = c("two.sided",
                                       "greater",
                                       "less")){
  alternative = match.arg(alternative)
  es_est = metafor::escalc(
     measure = "CVR",
     n1i = n1,
     n2i = n2,
     sd1i = sd1,
     sd2i = sd2,
     m1i = mean1,
     m2i = mean2,
  )
  
  teststat <- es_est$yi/sqrt(es_est$vi)
  pval = p_from_z(teststat, 
                  alternative = alternative)
  
  statistic = teststat
  names(statistic) = "z"
  null1 = 0
  names(null1) = "log ratio of the coefficients of variation"
  estimate = es_est$yi
  names(estimate)= "log VR"
  std_err = sqrt(es_est$vi)
  sum_stats = paste0("CV1 = ", sd1,"/",mean1, ", CV2 = ", CV2,"/",mean2)
  
  rval <- list(statistic = statistic,
               p.value = as.numeric(pval),
               stderr = std_err,
               null.value = null1,
               alternative = alternative,
               method = "log Ratio of CVs",
               data.name = sum_stats)
  class(rval) <- "htest"
  return(rval)
}


breusch_pagan(mtcars_lm,
              auxdesign =  as.matrix(mtcars$wt),
              koenker = FALSE)

test1=glejser(mtcars_lm,
              auxdesign =  as.matrix(mtcars$wt),
              statonly=FALSE)
glejser.test(mtcars_lm,
             "wt",
             mtcars)
