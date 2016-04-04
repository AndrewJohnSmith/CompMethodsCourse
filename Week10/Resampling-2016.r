# Most code borrowed from:

# Tom Carsey (carsey@unc.edu) and Jeff Harden (jjharden@unc.edu)

#---------------------------
# Sampling variability 
#----------------------------

# What exactly are standard errors? Answer: an estimate of the sampling variability of some statistic of interest.

# Let's say we draw a sample of 500 observations from a population with a mean of 5 and a standard deviation of 10. The true standard error of that mean is the population SD divided by the square root of N.
N <- 1000
true.mean <- 5
true.sd <- 50
true.se <- true.sd/sqrt(N)
a.sample <- rnorm(N, true.mean, true.sd)

# We can calculate the mean of our sample and the standard deviation of our sample. This tells us our ESTIMATE of the population mean and our ESTIMATE of the population standard deviation. To get an ESTIMATE of the sampling variability of that mean we just calculated---the standard error of the mean---we use the formula s/sqrt(N), where s is our ESTIMATE of the population standard deviation and N is the sample size.
sample.mean <- mean(a.sample)
sample.sd <- sd(a.sample)
sample.se <- sample.sd/sqrt(N)

# All of these are close, but not exactly equal to the true values. Now let's check to see that in repeated samples, the sample means and SDs converge to the true values.

rep.mean <- numeric(10000)
rep.sd <- numeric(10000)
rep.se <- numeric(10000)

for (i in 1:10000){
new.sample <- rnorm(N, true.mean, true.sd)
rep.mean[i] <- mean(new.sample)
rep.sd[i] <- sd(new.sample)
rep.se[i] <- sd(new.sample)/sqrt(N)
}

par(mfrow = c(3, 1))
hist(rep.mean, breaks = 50)
abline(v = true.mean, col = "red", lwd = 4)

hist(rep.sd, breaks = 50)
abline(v = true.sd, col = "red", lwd = 4)

hist(rep.se, breaks = 50)
abline(v = true.se, col = "red", lwd = 4)

# This looks pretty good. The point of all this is to show that a standard error is an ESTIMATE of the sampling variability of some statistic. In this case, that statistic was a simple mean. It could also be an OLS coefficient (which, of course, is also a mean).

### Confidence intervals: Drawing a picture of the coefficient's variability with the standard error

# Definition: If the procedure for computing a 95% confidence interval is used over and over, 95% of the time the interval will contain the true parameter value. In other words, you have a 95% chance of creating a 95% CI for a parameter that contains the true value. However, once you've done it, your CI either covers the parameter or it doesn't.  

# Let's calculate a 95 % CI from the original sample we drew. To do this, we must assume that the population is distributed normally (in this case we know it is). Then we add/subtract 1.96*the sample SE to the sample mean. Why do we choose 1.96? Because in a normal distribution, 95% of the density lies between -1.96 and 1.96.

ci.lo <- sample.mean - 1.96*sample.se
ci.up <- sample.mean + 1.96*sample.se
ci.lo; ci.up

# This includes the true parameter (5), so this time it worked. But how do we know if this method for calculating the CI is correct? We need to know that most of the time (i.e., 95%), the confidence interval does include the true parameter. We can examine that from the simulation we did earlier. As this graph shows, sometimes the CI does not include 5:

par(mar = c(5, 5, 3, 3))
plot(seq(1, 100, length = 100), seq(-4, 16, length = 100), type = "n", axes = FALSE, xlab = "Samples", ylab = "True Values")
box()
axis(1, at = seq(0, 100, 10))
axis(2, at = c(round(true.mean - 1.96*true.se, digits = 2), true.mean, round(true.mean + 1.96*true.se, digits = 2)), las = 2)
abline(h = true.mean, lwd = 2, col = "red")
abline(h = true.mean + 1.96*true.se, lwd = 2, lty = 2, col = "red")
abline(h = true.mean - 1.96*true.se, lwd = 2, lty = 2, col = "red")
for (i in 1:100){
points(i, rep.mean[i], lwd = 2, col = ifelse(rep.mean[i] + 1.96*rep.se[i] > true.mean & rep.mean[i] - 1.96*rep.se < true.mean, "steelblue", "forestgreen"))
segments(i, rep.mean[i] + 1.96*rep.se[i], i, rep.mean[i] - 1.96*rep.se, lwd = 2, col = ifelse(rep.mean[i] + 1.96*rep.se[i] > true.mean & rep.mean[i] - 1.96*rep.se < true.mean, "steelblue", "forestgreen"))
}

# More formally, we can calculate a "coverage probability." To do this, we construct a confidence interval for each run in the simulation as shown in the graph above, note whether that confidence interval includes the true parameter, and then calculate what proportion of the simulated CIs includes the true parameter. This proportion should be 0.95 if the method is working correctly.

cp <- ifelse(rep.mean - 1.96*rep.se < true.mean & rep.mean + 1.96*rep.se > true.mean, 1, 0)
coverage.probability <- mean(cp)
coverage.probability

#------------------------------------------------------------------------------
# Resampling methods: bootstrap, randomization tests
#------------------------------------------------------------------------------

library(Design)
library(boot) 
library(bootstrap)

# Global variables
set.seed(111007)
n <- 1000
a <- 0
b1 <- .5
b2 <- .25

# runif generates random deviates
# n = number of observations. If length(n) > 1, the length is taken to be the number required
# How many observations are we using??
x1 <- runif(n, min = -1, max = 1)  # Creates a vector of random values between -1 and 1
x2 <- runif(n, min = -1, max = 1)  # Creates another vector of random values between -1 and 1
y <- a + b1*x1 + b2*x2 + rnorm(n, 0, 2) # Creates a vector of random values based on the terms above


# First example: Create a vector of names of people
class.names <- c("Santiago", "Derek", "Allan", "Jun", "Karina", "Hanna", "John", "Elizabeth", "Louis", "Chelsea", "Florian", "Tom", "Jeff")

sample(class.names, replace = TRUE)  # Sample the class names with replacement

#***********************************
# Bootstrapping function
#***********************************
bootstrap <- function(model, boot.reps){ # Inputs: the name of the model and the number of bootstrap replications
boot.est <- matrix(NA, nrow = boot.reps, ncol = model$rank) # This step creates an empty matrix for the bootstrap estimates
for(i in 1:boot.reps){ # Start bootstrap loop
	sample <- sample(length(model$residuals), replace = TRUE) # Select observations to go in bootstrap sample
  datai <- NULL # Initiate bootstrap sample data
	for(j in 1:length(sample)){ # Bind selected observations...
		datai <- rbind(datai, model$model[sample[j], ]) #...by row to make a matrix of data from the bootstrap sample
	} # End bootstrap sample data
	boot.est[i , ] <- coef(lm(y ~ x1 + x2, data = datai)) # Run the model on the bootstrap sample, then collect coefficients in boot.est matrix
cat("Completed", i, "of", boot.reps, "bootstrap replications", "\n")
} # End bootstrap loop
return(list(boot.est = as.matrix(boot.est), boot.vcv = cov(boot.est))) # Return the matrix of coefficient estimates and the bootstrap covariance matrix
} # End function
#***********************************

model.1 <- lm(y ~ x1 + x2)   # Run a regression model with two predictors

system.time( # Check how long it takes
model.1.boot <- bootstrap(model.1, 100)  # Runs the bootstrap function of the above regression model
)

ses <- cbind(sqrt(diag(vcov(model.1))), sqrt(diag(model.1.boot[[2]])))
colnames(ses) <- c("Conventional SE", "SP Bootstrapped SE")
ses   # Outputs the standard errors of the regular regression vs. the bootstrapped regression

# The next line creates the standard regression coefficients and confidence intervals as well as those generated by the bootstrapped model
cis <- cbind(model.1$coef, model.1$coef - 1.96*ses[ , 1], model.1$coef + 1.96*ses[ , 1], model.1$coef - 1.96*ses[ , 2], model.1$coef + 1.96*ses[ , 2], quantile(model.1.boot[[1]][ , 2], .025), quantile(model.1.boot[[1]][ , 2], .975))
colnames(cis) <- c("Coef", "Conventional CI - Lower", "Conventional CI - Upper", "SP Bootstrap CI - Lower", "SP Bootstrap CI - Upper", "NP Bootstrap CI - Lower", "NP Bootstrap CI - Upper")
cis

par(mar = c(4, 7, 2, 2))
plot(seq(.2, .9, length = 3), seq(-.5, .5, length = 3), type = "n", axes = FALSE, xlab = "", ylab = "")
box()
points(model.1$coef[2], -.25)
points(model.1$coef[2], 0)
points(model.1$coef[2], .25)
segments(cis[2, 2], -.25, cis[2, 3], -.25)
segments(cis[2, 4], 0, cis[2, 5], 0)
segments(cis[2, 6], .25, cis[2, 7], .25)
axis(2, at = c(-.25, 0, .25), labels = c("NP Boot", "SP Boot", "Conv. SE"), las = 2)

#------------------------------------------------------------
# Bootstrapping with the boot package
#------------------------------------------------------------
# First you need to tell it what you want to bootstrap with a function
betas <- function(formula, data, indices) {
  d <- data[indices, ] # allows boot to select sample
  fit <- lm(formula, data = d)
  return(fit$coef)
} 

# Then you use the boot() command
system.time(
boot.results2 <- boot(data = as.data.frame(cbind(y, x1, x2)), statistic = betas, R = 1000, formula = y ~ x1 + x2)
)
boot.results2

# The results produce values for each term of the model
# In this case the intercept, x1, and x2

# Plot the bootstrap results for each term of the model
plot(boot.results2, index = 1)
plot(boot.results2, index = 2)
plot(boot.results2, index = 3)

# Confidence intervals
boot.ci(boot.results2, type = c("norm", "perc"), index = 1)
boot.ci(boot.results2, type = c("norm", "perc"), index = 2)
boot.ci(boot.results2, type = c("norm", "perc"), index = 3)

# Compare all these results to the bootstrap function used previously (above)

#**************************************
# Randomization test function
#**************************************

rand.test <- function(data, variable, reps){ # Inputs: the data, the variable to be randomized, and the number of replications
t.statistics <- numeric(reps) # Empty vector to store randomized t-statistics 
for(i in 1:reps){ # Start randomization loop
new.var <- sample(variable, replace = FALSE) # Reshuffle the values of the randomized variable by sampling without replacement
ifelse(variable == x1, new.data <- data.frame(y, new.var, x2), new.data <- data.frame(y, x1, new.var)) # Create a new data matrix with the reshuffled variable. The ifelse() command adjusts to whether x1 or x2 was selected.
colnames(new.data) <- c("y", "x1", "x2") # Name the variables of the new data
new.model <- lm(y ~ x1 + x2, data = new.data) # Estimate the model on the new data
ifelse(variable == x1, t.statistics[i] <- new.model$coef[2]/sqrt(vcov(new.model)[2, 2]), t.statistics[i] <- new.model$coef[3]/sqrt(vcov(new.model)[3, 3])) # Store the t-statistic from the new model. The ifelse() command is again used to adjust to the selection of x1 or x2.
cat("Completed", i, "of", reps, "replications", "\n")
} # End randomization loop
return(t.statistics) # Return the t-statistics
}
#***************************************

#----------------------------------------
# Conduct a randomization test
#----------------------------------------

system.time(
x1.rand.test <- rand.test(data.frame(y, x1, x2), x1, 500)
)
system.time(
x2.rand.test <- rand.test(data.frame(y, x1, x2), x2, 500)
)

# Plot the distribution of t-statistics with 95% rejection regions
d <- density(x1.rand.test)
d.hi <- d$x > quantile(x1.rand.test, .975) 
d.lo <- d$x < quantile(x1.rand.test, .025) 

d2 <- density(x2.rand.test)
d.hi2 <- d2$x > quantile(x2.rand.test, .975) 
d.lo2 <- d2$x < quantile(x2.rand.test, .025) 

par(mfrow = c(2, 1))
plot(d, xlim = c(-8, 8), axes = FALSE)
polygon(x = c(d$x[d.hi == TRUE], rev(d$x[d.hi == TRUE])), y = c(d$y[d.hi == TRUE], rep(0, times = length(d$y[d.hi == TRUE]))), col = "gray70")
polygon(x = c(d$x[d.lo == TRUE], rev(d$x[d.lo == TRUE])), y = c(d$y[d.lo == TRUE], rep(0, times = length(d$y[d.lo == TRUE]))), col = "gray70")
abline(v = model.1$coef[2]/sqrt(vcov(model.1)[2, 2]), lty = 2, col = "red")
axis(1, seq(-8, 8, 1))
box()

plot(d2, xlim = c(-8, 8))
polygon(x = c(d2$x[d.hi2 == TRUE], rev(d2$x[d.hi2 == TRUE])), y = c(d2$y[d.hi2 == TRUE], rep(0, times = length(d2$y[d.hi2 == TRUE]))), col = "gray70")
polygon(x = c(d2$x[d.lo2 == TRUE], rev(d2$x[d.lo2 == TRUE])), y = c(d2$y[d.lo2 == TRUE], rep(0, times = length(d2$y[d.lo2 == TRUE]))), col = "gray70")
abline(v = model.1$coef[3]/sqrt(vcov(model.1)[3, 3]), lty = 2, col = "red")
axis(1, seq(-8, 8, 1))
box()

