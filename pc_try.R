rm(list = ls())
#Load modified functions
source("try_pc.R")
source("try_skeleton.R")
source("try_udag2pdag.R")
source("try_udag2pdagRelaxed.R")
library(pcalg)
library(Rgraphviz)
library(COUNT)
library(plyr)


# Data 
#a = rnorm(n)
#b =    rnorm(n)
#c =   0.3*a+ 0.45*b +  rnorm(n)
#d =  0.45*c + 0.65*b+   rnorm(n)
# e = 0.6*c  +  rnorm(n)
# f = .82*b + 1.23*e +  rnorm(n)
# g = 0.45* e + rnorm(n)
# h = 0.23 * f + .29 * d + rnorm(n)
# i = 0.19* f + 0.22* d +  rnorm(n)
# j = 0.78* i + 0.47*h + rnorm(n)
#unscaleddat = matrix(c(a,b,c,d), nrow = n)
#dat = scale(unscaleddat)

#SS = list(C = cor(dat), n = nrow(dat))
#  model = try_pc(SS, indepTest =gaussCItest,  labels = c("A","B","C","D"), alpha = 0.01, verbose = TRUE)

#Discretize continuous data 
# dd = matrix(0L, nrow(dat), ncol(dat))
# for (i in 1:ncol(dat))
# {  dd[,i] <- discretize(dat[,i], method="frequency", categories = 4, onlycuts = FALSE)
# }
# ddat = dd - 1

#Binary Data
n = 10000
a <- rbinom(n, 1, 0.4)
e <- rbinom(n, 1, 0.4)
b <- as.integer(e|rbinom(n, 1, 0.7))
c <- as.integer((a|b)&rbinom(n,1,0.25))
d <- as.integer(as.integer((c&e))| (ed <- rbinom(n,1,0.01)))
dat = matrix(c(a,b,c,d,e), nrow = n)
SS = list(dm = dat, adaptDF = FALSE)
model = fci(SS, indepTest = binCItest, labels = c("A","B","C","D","E"), alpha = 0.01)


# #Discrete Data
# a <- sample(0:2, n, TRUE) 
# b <- sample(0:3, n, TRUE) 
# c <- a + b + 2*sample(0:1,n, TRUE) 
# d <-  c + b + sample(0:1,n,TRUE)
# dat = matrix(c(a,b,c,d), nrow = n)


#SS = list(dm = dat, nlev = c(7,4,6,3),adaptDF = FALSE)
#model = try_pc(SS, indepTest = disCItest, alpha = 0.05,labels = c("Day","Time","Delay","Lanes") , verbose = FALSE)

# Getting model from the data
#dev.off()
quartz()
plot(model)


#Adjacency Matrix of the graph
# adjm <- matrix(0L,ncol(dat), ncol(dat))
# q <- model@graph@edgeL
# sepset <- model@sepset
# for (i in c(1:ncol(dat)))
#    {
#   adjm[i, as.integer(q[[i]]$edges)] <- 1
#    }
 
# #2D Density Estimation
# x1 = runif(100)
# x2 = runif(100)
# x = cbind(x1,x2)
# dens <- kde2d(x1, x2)
# a = 0.4544; b = 0.766767
# index1 <- which(abs(dens$x-a)==min(abs(dens$x - a)))
# if(dens$x[[index1]]>a){index1<- index1 -1}
# index2 <- which(abs(dens$y-b)==min(abs(dens$y - b)))
# if(dens$y[[index2]]>b){index2<- index2 -1}
# pdf <- dens$z[index1,index2]

#Binary probability estimation
#Conditional Probability
 q <- count(dat)
 ind_num = which(q$x.3 ==1 & q$x.4 ==1)
 ind_den = which(q$x.3 ==1 )
 num = sum(q$freq[ind_num])
 den = sum(q$freq[ind_den])
 cond_prob = num/den
 print(cond_prob)

 #Do Probability
  x2_0 <- sum(q$freq[which(q$x.2 ==0 )])/n
  l1 = (sum(q$freq[which(q$x.4 ==1 & q$x.3 ==1 & q$x.2 ==0 )]))/(sum(q$freq[which( q$x.3 ==1 & q$x.2 ==0 )]))
 l2 = (sum(q$freq[which(q$x.4 ==1 & q$x.3 ==1 & q$x.2 ==1 )]))/(sum(q$freq[which( q$x.3 ==1 & q$x.2 ==1 )]))
 do_prob = l1*x2_0 + l2*(1-x2_0)
 print(do_prob)
 
#do-probability for multivariate gaussian (discretization approximation of integral)
# #MEAN 
# num = 0
# denom = 0
# for (i in 1:1000)
# {
#   k[i] = rnorm(1)
#   num = num + dnorm(k[i])*condMVN(mean = U, sigma = E, dependent.ind = c(4), given.ind = c(3,2), X.given = c(1,k[i]), check.sigma=TRUE)$condMean
#   denom = denom + dnorm(k[i])
# }
# num/(denom)
# #VARIANCE ( First sample and then find the sample variance)
# j = NULL
# for (i in 1:10000)
# {
#   k[i] = rnorm(1) #sample X2
#   condmean = condMVN(mean = U, sigma = E, dependent.ind = c(4), given.ind = c(3,2), X.given = c(1,k[i]), check.sigma=TRUE)$condMean
#   condvar = condMVN(mean = U, sigma = E, dependent.ind = c(4), given.ind = c(3,2), X.given = c(1,k[i]), check.sigma=TRUE)$condVar
#   
#   j[i] <- rnorm(1, mean = condmean, sd = condvar^(0.5))
# }
# var(j)