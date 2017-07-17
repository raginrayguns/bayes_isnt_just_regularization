library(dplyr)

# Formula for moments of the posterior distribution when sampling
# from a Cauchy distribution with scale 1 and unknown location,
# from Hanson, Wolf 1996 "Estimators for the Cauchy Distribution", 
# doi:10.1007/978-94-015-8729-7_20

I <- function(k,x)
{
  stopifnot(0 <= k && k< 2*length(x)-1)
  stopifnot(isTRUE(all.equal(k%%1,0)))
  Re(sum(sapply(1:length(x), function(i) (x[i]+1i)^k *
                  prod(sapply((1:length(x))[-i], function(j)
                    (1/((x[i]-x[j])^2+4))*(1-2i/(x[i]-x[j])))))))
}
posterior.moment <- function(k,x)
{
  stopifnot(I(0,x) > 1e-20)
  I(k,x)/I(0,x)
}

posterior.mean <- function(x) posterior.moment(1,x)

# Function to find a local maximum of the likelihood function
local.mle <- function(x, init, tol=10^-4)
{
  curr <- init
  repeat
  {
    # Calculate weights that downweight values far from the
    # current estimate
    unnormalized.weights <- 1/((x - curr)^2 + 1)
    # Normalize the weights
    w <- unnormalized.weights / sum(unnormalized.weights)
    # Next value is a weighted average of the data
    nextval <- sum(w*x)
    
    # Stop if we're within tolerance
    if (abs(nextval-curr) < tol)
    {
      break
    }
    
    curr <- nextval
  }
  
  return(curr)
}

# Find MLE by doing a local search around each data point--that should
# work, right? I hope so
iterative.mle <- function(x)
{
  candidates <- sapply(x, function(xi) local.mle(x, xi))
  scores <- sapply(candidates, function(m) sum(log(dcauchy(x, location=m))))
  return(candidates[which.max(scores)])
}

set.seed(1020)
samples <- replicate(10^5, rcauchy(5), simplify=FALSE)
postmeans <- sapply(samples, posterior.mean)
mles <- sapply(samples, iterative.mle)
# Save the results
save(samples, postmeans, mles, file="comparison.Rdata")
#load("comparison.Rdata")

# Categorize into examples that were especially good for the posterior mean,
# and especially bad
results <- tibble(run=1:10^5, postmeans=postmeans, mles=mles) %>% mutate(advantage=mles^2 - postmeans^2)
good <- results %>% filter(advantage > quantile(advantage, .999))
bad <- results %>% filter(-advantage > quantile(-advantage, .999))

# Plot good and bad ones for browsing
for (i in good$run)
{
  x <- samples[[i]]
  plot(Vectorize(function(m) sum(dcauchy(x, location=m, log=TRUE))), xlim=2*(range(x) - median(x)) + median(x), main=i)
  abline(v=postmeans[i], col='blue')
  abline(v=mles[i], col='red')
}
for (i in bad$run)
{
  x <- samples[[i]]
  plot(Vectorize(function(m) sum(dcauchy(x, location=m, log=TRUE))), xlim=2*(range(x) - median(x)) + median(x))
  abline(v=postmeans[i], col='blue')
  abline(v=mles[i], col='red')
}

# The one I'm using for my illustration
png("comparison.png", width=540, height=540)
i <- 97729
x <- samples[[i]]
plot(Vectorize(function(m) sum(dcauchy(x, location=m, log=TRUE))),
     xlim=2*(range(x) - median(x)) + median(x),
     main="Comparison of ML/MAP and posterior mean",
     xlab="mu", ylab="likelihood")
abline(v=mles[i], col='orange')
abline(v=postmeans[i], col='cyan3')
abline(v=0, col="black")
legend(12, -21.66604,
       c("ML/MAP", "posterior mean", "truth"),
       c("orange", "cyan3", "black"))
dev.off()

# Comparison of rmse
c(sqrt(mean(postmeans^2)), sqrt(mean(mles^2)))
  