
# code base of the algorithm proposed in the paper "Bi-objective parameter setting problem of a Genetic Algorithm: An empirical study on Traveling Salesperson Problem"
# the authors are Yavuzhan Akdurana, Erdi Dasdemir, Murat Caner Testik from Department of Industrial Engineering, Hacettepe University, 06800, Ankara, Turkey
# you may contact edasdemir@hacettepe.edu.tr for your questions.
########################################################################################################################################################


# set seed
# setting a seed will ensure that we can obtain the same results every time we run the algorithm
# if the seed is not set, we may get different results every time as the algorithm inherits randomness
set.seed(70646)
#############################################################


# set the path
path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)
#############################################################

# import libraries
library(GA)
library(mco)
library(nsga2R)
library(compiler)
library(combinat)
library(tictoc)
library(FrF2)
library(colorspace)
library(gtools)


# read tsp instance
# correct path of the tsp instance should be set here
# or this file should be  directly located in the corresponding instance folder
Data <- read.csv("sample_data.csv", header = FALSE) # SET THE PATH !

# set if the tyep of input data is EUC_2D or MATRIX (see the paper and its supplementary document for the details)
input_type = "EUC_2D"

if (input_type == "EUC_2D"){
  # function that creates a distance matrix for TSP instances having EUC_2D type input (coordinate information of the nodes in (x, y) format is given)
  # matrix data defines distances between locations
  DistanceMatrix <- function(x){
    DM <- matrix(nrow = length(x[,1]),ncol = length(x[,1]))
    for (i in 1:length(x[,1])) {
      for (j in 1:length(x[,1])) {
        DM[i,j] <- sqrt((x[i,1]-x[j,1])^2+(x[i,2]-x[j,2])^2)
      }
    }
    return(DM)
  }
  
  D <- DistanceMatrix(Data) # distance matrix for EUC_2D problems
}else{
  D <- Data # MATRIX problems have data in the form of distance matrix. The data read can be used directly.
}
#############################################################


# define functions
# genetic algorithm function to find the fitness value of the corresponding setting (population member)
functionGA<- function (x){#GA algorithm function
 
   tic()  # Starts the timer and stores the start time
  
  GA<- ga(type = "permutation", fitness = fitness, distMatrix = D,min = 1, max = length(Data[,1]),
          
          #Cromosome structure
          selection = Fact1[x[1]], #Selection factor
          crossover = Fact2[x[2]], #Crossover factor
          mutation = Fact3[x[3]],  #Mutation factor 
          popSize=x[4],            #Population size factor
          pcrossover=x[5],         #Crossover probability factor
          pmutation=x[6],          #Mutation probability factor
          maxiter=x[7])            #Number of generation factor
  
  exactime <- toc() #Notes the current timer and computes elapsed time since the matching call to tic(). 
  
  exactime <- as.numeric(exactime$toc - exactime$tic) #Run time of the function
  
  return(c(GA@fitnessValue,exactime)) #F1 and F2 values
}


# crossover function 
# cross_no and no_cross are number of genes to be crossover
# cross_prob and prob_cross are crossover probability
# parent_chromosome is chromosomes of the parents
# it creates new offspring by determining how many and which genes will cross over, based on the probability of crossover.
mixedCrossover<- function(parent_chromosome,cprob,nocros){
  
  popSize=nrow(parent_chromosome)  # Population size
  varNo=ncol(parent_chromosome)   # Number of genes
  child <- parent_chromosome      # Keeps the chromosomes of the children's parents
  c=1
  childAll=c()
  
  for (c in 1:(popSize/2)){ 
    
    childPairs<- sample(nrow(child),2) #Randomly selected two parents
    
    parent1<-child[childPairs[1],] # First parent
    parent2<-child[childPairs[2],] # Second parent
    
    if(runif(1)<cprob){
      rndpos<-sort(sample(1:varNo,nocros)) #Sample random gene positions between 1-7 to exchange.
      
      child1<-parent1[rndpos] #child1 = parent1 what values are available in rndpos
      child2<-parent2[rndpos] #child2 = parent2 what values are available in rndpos
      
      parent1[rndpos]<-child2 #We are replacing the values in rndpos of parent1 with values from parent2 in the same position (child2).
      parent2[rndpos]<-child1 #We are replacing the values in rndpos of parent2 with values from parent1 in the same position (child1).
    }
    
    child<-child[-childPairs,] #Removes chromosomes from parents that have crossed over from the list
    childAll<-rbind(childAll,parent1,parent2) #Chromosomes of all parents and new children
    
  }
  rownames(childAll)=NULL
  return(childAll)
}


# mutation function
# there are 7 "if statements" in the function her as 7 parameters are considered in the paper.
# please adjust the number of "if statements" according to your case.
# each if decides whether a mutation will be implemented to the corresponding gene (parameter) or not.
# m_prob = p_mut -> mutation probability of each gene
# parent_chromosome -> chromosomes of the parents
mixedMutation<- function(parent_chromosome,mprob){
  popSize=nrow(parent_chromosome)  # population size
  varNo=ncol(parent_chromosome)    # number of genes
  child = parent_chromosome        # keeps the chromosomes of the children's parents

  childMut=t(apply(child,1, function (m) { 

    
    probGenes=runif(varNo) # generate 7 random probability values from a uniform distribution by number of genes
    
    applyMut=which(probGenes<mprob) # if the probability value of the gene is less than the mutation probability value, the mutation occurs.
    
    if(1 %in% applyMut){
      m[1]<-sample((1:4)[1:4!= m[1]],1) # it takes the value 2 if it is currently 1, and 1 if it is 2.
    }
    
    if(2 %in% applyMut){
      m[2]<-sample((1:4)[1:4!= m[2]],1) # it takes the value 2 if it is currently 1, and 1 if it is 2.
    }
    
    if(3 %in% applyMut){
      m[3]<-sample((1:5)[1:5!= m[3]],1) # it takes the value 2 if it is currently 1, and 1 if it is 2.
    }
    
    if(4 %in% applyMut){
      m[4]<-sample(Fact4[Fact4!= m[4]],1) # selects a new value from the range of values that will not be the same as the current value.
    } 
    
    if(5 %in% applyMut){
      m[5]<-sample(Fact5[Fact5!= m[5]],1) # selects a new value from the range of values that will not be the same as the current value.
    } 
    
    if(6 %in% applyMut){
      m[6]<-sample(Fact6[Fact6!= m[6]],1) # selects a new value from the range of values that will not be the same as the current value.
    } 
    
    if(7 %in% applyMut){
      m[7]<-sample(Fact7[Fact7!= m[7]],1) # selects a new value from the range of values that will not be the same as the current value.
    } 
    
    m
    
  }))
  return(childMut)
}

# function to calculate tour length
tourLength <- function(tour, distMatrix) {
  tour <- c(tour, tour[1])
  route <- embed(tour, 2)[,2:1]
  sum(distMatrix[route])
}


# fitness function to find tour length
# the GA function works as a maximization. To convert it into a minimization problem, its exponent is written as -1
fitness <- function(tour, ...) 1/tourLength(tour, ...) #The GA function works as a maximization. To convert it into a minimization problem, its exponent is written as -1
################################################################################################################################


# set parameters of the multi-objective algorithm (NSGAII)
popSize <- 4 # Initial population Size 
tourSize <- 2 # Size of tournament
generations <- 4 # Number of generations
pcross <- 0.9 # Crossover probability
crosNo <- 3 # Number of genes to be crossovered
objDim <- 2 # Number of objectives considered, f1: f value, f2:time
nfact <- 7 # Number of factors considered
pmut <- 1-((1-0.2)^(1/nfact)) # Mutation probability of each gene--> The overall prob of a chromosome to be crossovered is 0.2 here.
################################################################################################################################


#-------------Genetic algorithm factors----------------
Fact1 <- c("gaperm_lrSelection","gaperm_nlrSelection","gaperm_rwSelection","gaperm_tourSelection") #Selection Operator
Fact2 <- c("gaperm_cxCrossover","gaperm_pmxCrossover","gaperm_oxCrossover","gaperm_pbxCrossover") #Crossover Operato
Fact3 <- c("gaperm_simMutation", "gaperm_ismMutation","gaperm_swMutation","gaperm_dmMutation","gaperm_scrMutation") #Mutation Operator
Fact4 <- c(50,100,150,200,300,400,500,600,700,800,900,1000) #Population Size
Fact5 <- c(0.1,0.25,0.5,0.75,0.9,1) #Crossover Probability
Fact6 <- c(0.1,0.25,0.5,0.75,0.9,1) #Mutation Probability
Fact7 <- c(50,100,150,200,300,400,500,600,700,800,900,1000) #Number of generations


# start timer for multi-objective algorithm
ptm <- proc.time() #proc.time determines how much real and CPU time (in seconds) the currently running R process has already taken


# generate a random parent population----------------
parent <- t(sapply(1:popSize,function (i) {sample(1:4,2,replace=TRUE)}))
parent <- cbind(parent,(sapply(1:popSize,function (i) {sample(1:5,1,replace=TRUE)})))
parent <- cbind(parent, sample(Fact4,popSize,replace=TRUE), sample(Fact5,popSize,replace = TRUE), sample(Fact6,popSize,replace = TRUE), sample(Fact7,popSize,replace=TRUE))


# add objective function columns to the parent matrix----------------
parentf1f2 <- apply(parent,1,functionGA) # The exponent was written as -1 to turn it into a minimization problem. The process was repeated to return it to its old value.
parentf1f2[1,] <- (1/parentf1f2[1,])
parentf1f2 <- round(parentf1f2,6) # To round values like 1.0000000000000000000012, 1.0000000000000000000000011
parent <- cbind(parent,t(parentf1f2))


# plot of initial population----------------
mypath <- file.path(path,paste("0myplot_", "initial", ".png", sep = ""))
png(file=mypath)
mytitle = paste("Initial population")
plot(parent[,(nfact+1):(nfact+2)], col=1, pch=1, xlab="f1 value (function)", ylab="f2 value (time)",main = mytitle, yaxs="i",xaxs="i")
dev.off()


# sorting non-dominated values----------------
ranking <- fastNonDominatedSorting(parent[,(nfact+1):(nfact+2)]) #A fast approach to sort non-dominated solutions into different nondomination levels
rnkIndex <- integer(popSize) 
i <- 1
while (i <= length(ranking)) {#creating index numbers of ranks
  rnkIndex[ranking[[i]]] <- i
  i <- i + 1
} 
parent <- cbind(parent,rnkIndex)#assigns index numbers of rows


# crowding ranking ----------------
maxGlobalf <- apply(parent[,(nfact+1):(nfact+2)], 2, max) #Global maximum f1 and f2 values
minGlobalf <- apply(parent[,(nfact+1):(nfact+2)], 2, min) #Global minimum f1 and f2 values
objRange <- maxGlobalf-minGlobalf #objective function range
cd <- crowdingDist4frnt(parent,ranking,objRange) #crowding distance function
parent <- cbind(parent,apply(cd,1,sum)) #adding result of crowding distance function to population record


InitialPopulation <- parent # we record the initial population


##########################ITERATIONS##################################
iter<-1
for (iter in 1:generations){
  
  cat("---------------generation---------------",iter,"starts")

  #-------------Tournament Selection----------------  
  matingPool <- tournamentSelection(parent,popSize,tourSize)
  
  #-------------Apply Crossover---------------- 
  childAfterX<- mixedCrossover(matingPool[,1:nfact],pcross,crosNo)
  
  #-------------Apply Mutation----------------   
  childAfterM<- mixedMutation(childAfterX, pmut)
  
  #-------------Apply GA to find objective function values----------------
  cf1f2 <- apply(childAfterM,1,functionGA)
  cf1f2[1,] <- (1/cf1f2[1,]) # The exponent was written as -1 to turn it into a minimization problem. The process was repeated to return it to its old value.
  cf1f2 <- round(cf1f2,6) # To round values like 1.0000000000000000000012, 1.0000000000000000000000011
  
  childAfterM<-cbind(childAfterM,t(cf1f2))# Add f1 and f2 values to the offspring matrix (childAfterM)
  
  parentNext <- rbind(parent[,1:(nfact+objDim)],childAfterM) # Combine the parent with the childAfterM
  
  #-------------Plot of results----------------
  mypath <- file.path(path,paste("0myplot_", iter, ".png", sep = ""))
  png(file=mypath)
  mytitle = paste("Initial Population of generation ",iter)
  plot(parentNext[,(nfact+1):(nfact+objDim)], col=iter+1, pch=1, xlab="f1 (objective)", ylab="f2 (time)",main = mytitle, yaxs="i",xaxs="i")
  dev.off()
  
  #-------------sorting non-dominated values----------------
  ranking <- fastNonDominatedSorting(parentNext[,(nfact+1):(nfact+2)])
  rnkIndex <- integer(popSize)
  i <- 1
  while (i <= length(ranking)) {
    rnkIndex[ranking[[i]]] <- i
    i <- i + 1
  } 
  parentNext <- cbind(parentNext,rnkIndex)
  

  #-------------Crowding ranking ----------------
  maxf<-apply(parentNext[,(nfact+1):(nfact+objDim)], 2, max)
  minf<-apply(parentNext[,(nfact+1):(nfact+objDim)], 2, min)
  objRange<-maxf-minf
  maxGlobalf<-apply(rbind(maxGlobalf,maxf),2,max) ###!!!!HOCAM BURAYI KULLANMAMIÞIZ maxf ve minf KULLANMIÞIZ
  minGlobalf<-apply(rbind(maxGlobalf,minf),2,min) ###!!!!HOCAM BURAYI KULLANMAMIÞIZ maxf ve minf KULLANMIÞIZ
  cd <- crowdingDist4frnt(parentNext,ranking,objRange)
  parentNext <- cbind(parentNext,apply(cd,1,sum))

  
  parentNext.sort <- parentNext[order(parentNext[,nfact+objDim+1],-parentNext[,nfact+objDim+2]),]#Population records are ordered by objective functions
  
  
  # choose the first 'popSize' rows for next generation 
  parent <- parentNext.sort[1:popSize,]
  
  #-------------Plot of results----------------
  mypath <- file.path(path,paste("0myplot_", iter,"_2", ".png", sep = ""))
  png(file=mypath)
  mytitle = paste("Final Population of generation ",iter)
  plot(parent[,(nfact+1):(nfact+2)], col=iter+1, pch=1, xlab="f1 (objective)", ylab="f2 (time)",main = mytitle, yaxs="i",xaxs="i")
  dev.off() 
  
  
}

#-------------Time for Multiobjective Algorithm ends----------------
SimTime<-proc.time() - ptm
SimTimeF<-SimTime[1] + SimTime[2]

#-------------Recording results----------------
write.table(parent, file.path(path,paste("parent.txt")),sep = "\t")
write.table(parent[which(parent[,(nfact+objDim+1)]==1),], file.path(path,paste("parent_nondominated.txt")),sep = "\t")
plot(parent[,(nfact+1):(nfact+2)], col=iter+1, pch=1, xlab="f1 (objective)", ylab="f2 (time)",main = mytitle, yaxs="i",xaxs="i")

