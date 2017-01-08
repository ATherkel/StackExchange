
## ---- initialize ----

## Seed for replication
set.seed(3)


## Gamma parameters
a <- 0.02
b <- 0.0025

## Length of (X_1, X_2, ..., X_n)
n <- c(10^(1:5),5e5,1e6,2e6)
## Number of simulations for each vector
nsim <- 1e2

## Estimate size of vector
size <- sum(n*nsim*8);
class(size) <- "object_size";
format(size,"Mb")

## 
sim <- lapply(n, function(num){
    matrix(rgamma(num  * nsim,a,b),nrow = num)
})


## ---- find_maxima ----

max_func <- function(num){
    numsim <- nrow(num)
    
    ## Max for each simulation run
    mymax <- sapply(1:ncol(num), function(i){
        max(num[,i])
    })
    
    ## Embrechts et al. (1997) norming constants
    c <- 1/b
    d <- c*(log(numsim)+(a-1)*log(log(numsim))-lgamma(a))
    
    mean1 <- mean(mymax)
    mean2 <- mean((mymax-d)/c)
    
    ## Norming constants (Von Mises functions)
    d3 <- qgamma(1-1/numsim,a,b)
    c3 <- b^(-1)*(1+(a-1)/(b*d3))
    
    mean3 <- mean((mymax-d3)/c3)
    
    ## Output
    c(meanreg = mean1, meannorm_em = mean2, meannorm_vm = mean3)
}

simmax <- lapply(sim,max_func)

## Extract the three vectors
simmax_reg <- sapply(simmax,function(x) x["meanreg"])
simmax_norm_em <- sapply(simmax,function(x) x["meannorm_em"])
simmax_norm_vm <- sapply(simmax,function(x) x["meannorm_vm"])


## Norming constants (Von Mises functions)
d3 <- function(x) qgamma(1-1/x,a,b)
c3 <- function(x)  b^(-1)*(1+(a-1)/(b*d3(x)))


## ---- plotting ----

par(mfrow = c(2,1))
par(mar = c(3,2,2,2)+0.1)


## Plot 1
plot(n,simmax_reg,main = bquote(paste(
    "Max of Gamma,  ",M[n])))
curve(c3(x)*(-digamma(1)) + d3(x),
      add = TRUE,lty = 2,col = 2)
legend("bottomright",legend = c(
    as.expression(bquote(paste(c[n]*gamma+d[n], "  ")))),
       lty = 2,col = c(2))


## Plot 2
plot(n,simmax_norm_vm,main = bquote(paste("Normed max of Gamma,  ",
                                    c[n]^-1*(M[n] - d[n]))),
     ylab = 'simmax_norm["mean",]')
abline(h = -digamma(1),lty = 2,col = 2)
abline(h = 1,lty = 2,col = 2)
legend("bottomright",legend = bquote(
    paste(gamma %~~%.(round(-digamma(1),3)),"  ","a")),
    lty = 2,col = 2)


par(mfrow = c(1,1))

 


