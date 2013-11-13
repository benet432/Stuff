
mini <- FALSE

#============================== Setup for running on Gauss... ==============================#

args <- commandArgs(TRUE)

cat("Command-line arguments:\n")
print(args)

####
# sim_start ==> Lowest possible dataset number
###

###################
sim_start <- 1000
###################

if (length(args)==0){
  sim_num <- sim_start + 1
  set.seed(121231)
} else {
  # SLURM can use either 0- or 1-indexing...
  # Lets use 1-indexing here...
  sim_num <- sim_start + as.numeric(args[1])
  sim_seed <- (762*(sim_num-1) + 121231)
}

cat(paste("\nAnalyzing dataset number ",sim_num,"...\n\n",sep=""))

# Find r and s indices:

jobnum=sim_num-1000
if (jobnum %% 50 ==0)
{
  s_index=jobnum/50
  r_index=50
}else
{
  s_index=as.integer(jobnum/50)+1
  r_index=jobnum %% 50
}

#============================== Run the simulation study ==============================#

# Load packages:
library(BH)
library(bigmemory.sri)
library(bigmemory)
library(biganalytics)


# mini or full?
if (mini){
	rootfilename <- "blb_lin_reg_mini"
} else {
	rootfilename <- "blb_lin_reg_data"
}

# Attach big.matrix :
dat= attach.big.matrix(dget("/home/pdbaines/data/blb_lin_reg_data.desc"), backingpath="/home/pdbaines/data/")

# Remaining BLB specs:

n=nrow(dat)
s = 5
r = 50
gamma = 0.7
b=round(n^gamma)

#sample b rows for the five subsets
set.seed(s_index)
index1=sample(1:n, b, replace = FALSE)

#sample of size n
set.seed(jobnum)
index2=rmultinom(1, n, prob=rep(1/b, b))

# Extract the subset:
data=data.frame(dat[index1, ])

# Fit lm:
coefficient=lm(X1001~.-1, data, weights=index2)$coefficient

# Output file & Save estimates to file:
write.table(coefficient, file = paste0("output/","coef_",sprintf("%02d",s_index),"_",sprintf("%02d",r_index),".txt"))