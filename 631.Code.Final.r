#STAT 631 Project
# ALberta's data on COVID-19 gathered from https://covid19stats.alberta.ca/ (04-13-2020)
#STAT 631 Project - Abraham Hussein
# Data importation
library(deSolve)
Data = read.csv("covid19dataexport.csv")

# sort by date
newdata <- Data[order(Data$Date.reported),] 
# Given I(t) and R(t)
It = c()
Rt = c()
for (i in 1:1732) {
  It = c(It, ifelse(newdata$Case.type[i] == "Confirmed",1,0))
if (newdata$Case.status[i] == "Recovered") {
  k = 1
} else if (newdata$Case.status[i] == "Died") {
  k= -1
} else {
  k = 0
  }
  Rt = c(Rt,k)
}
# Updating Dataframe
newdata$Case.type = It

newdata$Case.status = Rt

#Adding missing dates to the levels of data set
levels(newdata$Date.reported) = c(levels(newdata$Date.reported),"2020-03-07","2020-03-08")
#Dates 03-06 -> 04-13
Date = c(levels(newdata$Date.reported))
Date = sort(Date)
levels(Date) = c(Date)


### Generate Groups S(t), I(t), & R(t)
### Preamble
count = 0
pcount = 0
Daily.growth = c()
Probable.growth = c()
St = c()
popsize = 4371316 # Total population as Sept 2019 is 4,371,316 
                  # reported from plo


for (i in 1:39) {
  
# COunting the new cases
  new.cases = sum(ifelse(newdata$Date.reported == Date[i] & newdata$Case.type == 1,1,0))
  count = count  + new.cases
  Daily.growth = c(Daily.growth, count)
  
# Counting the probable casses
  prb.cases = sum(ifelse(newdata$Date.reported == Date[i] & newdata$Case.type == 0,1,0))
  pcount = pcount + prb.cases
  Probable.growth = c(Probable.growth, pcount)
# S(t) based on CONFIRMED cases
  St = c(St,popsize-count)
}

# Recovered vectors
Dt = c(rep(0,8),rep(1,3),rep(2,2),rep(3,2),rep(5,2),7,10,16,16,
      21,26,29,30,rep(35,2),rep(39,2),41,44,45,rep(46,2),rep(48,5))
Rt = c(rep(0,8),rep(1,4),3,4,rep(5,2),12,17,21,26,34,53,74,94,119,
       143,183,209,242,288,344,409,476,549,632,710,771,823,875)

# Figure 1 - Cases in Alberta for COVID-19
plot(1:39,Daily.growth,  ylim=c(0,1230),
     ylab = "Cases", col=2, xlab = "Days (since March 05, 2020)")
lines(1:39,Daily.growth,type = 'b', main = "Cases in Alberta for COVID-19", 
      ylab = "Cases", col=2, xlab = "Days (since March 05, 2020)")
lines(1:39,Probable.growth,type='l' ,lty=2, col=3)
lines(1:39,Probable.growth,type='b' ,lty=2, col=3)
legend(2,1200,legend = c(" I(t) - Confirmed Cases","P(t) - Probable Cases"), col = c(2:3), lty =c(1:2))

# Table 1 - Cases
Workingdata = data.frame(cbind(Date,c(1:39),Daily.growth,Probable.growth,Rt,Dt))
colnames(Workingdata) = c("Date", "Days", "Case (C)", "Case (P)", "Recovered", "Died")
write.csv(Workingdata, "table.csv")
colnames(Workingdata) = c("Date", "Days", "Case.C", "Case.P", "Recovered", "Died")



# Building our model
# a(t) + i(t) + p(t) + r(t) + d(t) = 1
attach(Workingdata)
#Preamble 
# Susceptible s(t)
s <- function(t, P = Probable.growth, C = Daily.growth) {
  s = popsize - P[t]  - C[t]
  s = s / popsize
  s
}
# Infected i(t)
i <- function(t,P = Probable.growth, C = Daily.growth) {
  I = C[t] + P[t]
  I = I / popsize
  I
}
# Recovered r(t)
r <- function(t, H = Rt, D = Dt) {
  R = H[t] + D[t]
  R = R / popsize
  R
}

parameter <- function(popsize = popsize, C, P, H, D){

  beta = c(0)
  Gamm = c(0)
  pred = c(0)
for (t in 2:39) {
  beta[t] = ( s(t-1) - s(t) ) / ( s(t-1) * i(t-1) )
  Gamm[t] = ( r(t) - r(t-1) ) / i(t-1)
  pred[t] = i(t-1) + beta[t] * ( s(t-1) * i(t-1) ) - Gamm[t] * i(t-1)
}
  param = data.frame(cbind(beta,Gamm,pred))
  
  return(param)
}

firstrun = parameter(C = Daily.growth, P = Probable.growth, H = Rt, D = Dt)
params = c(beta = mean(firstrun$beta),
           gamma =mean(firstrun$Gamm))
plot(1:39,firstrun$beta, type = 'l')
plot(1:39,firstrun$Gamm, type = 'l')


# SIR Model
library(deSolve)
#Function originally implemented from [3] 
SIR.model <-function(t, x, param){
  S = x[1]
  I = x[2]
  R = x[3]
  
  beta = params["beta"]
  gamma = params["gamma"]
  dS = -beta*S*I
  dI = beta*S*I-gamma*I
  dR = gamma*I
  dX = c(dS, dI, dR)
  list(dX)
}

params = c(beta = mean(firstrun$beta),
          gamma =mean(firstrun$Gamm))

times = 1:150

xstart = c(S=s(1), I=i(1), R=r(1))


out = as.data.frame(
  ode(func=SIR.model,
      y=xstart,
      times=times,
      parms=params)
  )

# Output
plot(S~time,data=out,col=3,type='l',xlab = "Days" )    
lines(I~time,data=out,col=2,type='l')
lines(R~time,data=out,col=4,type='l')
points(1:39,s(1:39), col = 3)
points(1:39,i(1:39), col = 2)
points(1:39,r(1:39), col = 4)
legend(120,0.8,legend = c("Susceptible","Infected","Removed"), col = c(2:4),lty=c(1,1,1))


# Simulations for the Spread of the virus
Network.Spatial.Simulation <- function(Contact, Beta, Gamma,n , times, points = T){
# Initialization
  
  pop = n^2
  I_0 =  sample(pop,1)
  Recovered = rep(0,pop) #Dummy Variable
  MatriX  = matrix(rep(0,n),ncol = n, nrow = n) # Network
  address = expand.grid(x = 1:n, y = 1:n)

# Updating Matrix to Infected status
Cough <- function(I_0, M = MatriX){
  M[address[I_0,1], address[I_0,2]] = 1
  return(M)
}

# Updating Matrix for neighbor infections
Cough.Neighbor <-function(r,c, M = MatriX, N, chance){
  #Neighbor 1
  if (is.na(N[1]) == F & N[1] != 1) {
    if (chance > runif(1)){  #Beta being infectivity rate and test for transmission
      M[r+1, c+1] = 1 }} 
  #Neighbor 2
  if (is.na(N[2]) == F & N[2] != 1) {
    if ( chance > runif(1)){  
      M[r+1, c] = 1 }}
  #Neighbor 3
  if (is.na(N[3]) == F & N[3] != 1) {
    if (chance > runif(1)){  
      M[r+1, c-1] = 1 }}
  #Neighbor 4
  if (is.na(N[4]) == F & N[4] != 1) {
    if (chance > runif(1)){  
      M[r, c+1] = 1 }}
  #Neighbor 5
  if (is.na(N[5]) == F & N[5] != 1) {
    if (chance > runif(1)){  
      M[r, c-1] = 1 }}
  #Neighbor 6
  if (is.na(N[6]) == F & N[6] != 1) {
    if (chance > runif(1)){  #Beta being infectivity rate and test for transmission
      M[r-1, c+1] = 1 }} 
  #Neighbor 7
  if (is.na(N[7]) == F & N[7] != 1) {
    if (chance > runif(1)){  
      M[r-1, c] = 1 }}
  #Neighbor 8
  if (is.na(N[8]) == F & N[8] != 1) {
    if (chance > runif(1)){  
      M[r-1, c-1] = 1 }}
  return(M)
}

# Immediate neighbors - Function from mrip [2] 
# Returns column of neighbors for position i
Neighbor = function(mat,col) {
    m2 = cbind(NA,rbind(NA,mat,NA),NA)
    addresses = expand.grid(x = 1:n, y = 1:n)
    ret = c()
    for(i in 1:-1)
    for(j in 1:-1)
    if(i!=0 || j !=0)
    ret = rbind(ret,m2[addresses$x+i+1+nrow(m2)*(addresses$y+j)]) 
    return(ret[,col])
}

# ALgorithim
  MatriX = Cough(I_0) #Ground zero
  sim.S = c() #Susceptible
  sim.I = c() #Infected
  sim.R = c() #Removed
# Initial plot
  if (points == T) {
  plot(address$y, address$x,  col=ifelse(Cough(I_0),2,8),
  pch =ifelse(Cough(I_0),19,20), main = "Spread of COVID-19", xlab = paste("Day ",0), ylab="")
  }
  
for (t in 1:times) {
  Infected = which(MatriX[,] == 1) #Assess infected individuals
  for (j in 1:length(Infected)) {
    r = address[Infected[j],1]
    c = address[Infected[j],2]
    
# Amount of neighbors being infected visits
    Hello.neighbor = floor(Contact[1] + rnorm(1,0,Contact[2]))
    if (Hello.neighbor > 8) {
      Hello.neighbor = 8             #Ensuring amount is within 0 and 8.
    } else if (Hello.neighbor < 0) {
      Hello.neighbor = 0
    }
    chance = Beta * (Hello.neighbor/8) #probability of infecting immediate neighbor
    Nextdoors = Neighbor(MatriX,Infected[j]) #Vector of Neighbors (0 - S, 1 - I, 2 - R)
    MatriX = Cough.Neighbor(r,c,M= MatriX,N = Nextdoors, chance)
  }
# Curves
    Result = as.numeric(MatriX)
  
# Recovery - **Possible Memory Leak**
  for (c in 1:pop) {
  if (isTRUE(Gamma > runif(1)) == TRUE & Result[c] == 1) { #Ensuring cured person has disease 
    Recovered[c] = 1
  }
  }
  sim.S = c(sim.S, pop - sum(Result) )
  sim.I = c(sim.I,       sum(Result) - sum(Recovered))
  sim.R = c(sim.R,                     sum(Recovered))
  
##  Graphics  ##
  
if (points == T) {
# Network Display of Disease Spread
  inf.color = ifelse(Result == 1,2,8)
  inf.shape = ifelse(Result == 1,15,15)
  cre.color = ifelse(Recovered == 1,4,inf.color)
  cre.shape = ifelse(Recovered == 1,15,inf.shape)
  plot(address$y,address$x, col=inf.color, pch = inf.shape,
       main = "Spread of COVID-19", xlab = paste("Day ",t), ylab = "")
  points(address$y,address$x, col=cre.color, pch = cre.shape)
}
}
  if (points == F) {
# S I R Density Curves
    plot (1:times, sim.S/pop, col="green", type = "l",
          ylim = c(0,1))
    lines(1:times, sim.I/pop, col="red")
    lines(1:times, sim.R/pop, col="blue")
  }
  
  
  Final = data.frame(cbind(sim.S/pop, sim.I/pop, sim.R/pop))
  colnames(Final) = c("S","I","R")
  
  return(Final)
}

# Based on Netherlands findings as The cases for two Countries are close [1]
Contact = c(mean = 13.85,
            sd = 10.54)

# Reduction of Contacts
Contact.half = Contact / 2
# Social Distancing measaures
Distance = c(0,1)



# simulation time and size have been reduced for your convenience:
# Simulation Without Social Distancing
est = Network.Spatial.Simulation (Contact,Beta = params[1],Gamma = params[2],n = 20,times = 20, points = T)

plot (est$S, col="green", type = "l", main= "Implementation of Self-Quarantine ",
      xlab= "Days", ylab = "Density",ylim = c(0,1))
lines(est$I, col="red")
lines(est$R, col="blue")
legend(0.82*length(est$S),0.9,legend = c("Susceptible","Infected","Removed"), col = c(3,2,4),lty=c(1,1,1))

# Simulation With Contact reduction
est = Network.Spatial.Simulation (Contact.half ,Beta = params[1],Gamma = params[2],n = 20,times = 20, points = T)

plot (est$S, col="green", type = "l", main= "Implementation of Self-Quarantine ",
      xlab= "Days", ylab = "Density",ylim = c(0,1))
lines(est$I, col="red")
lines(est$R, col="blue")
legend(0.82*length(est$S),0.9,legend = c("Susceptible","Infected","Removed"), col = c(3,2,4),lty=c(1,1,1))

# Simulation With Self Isolation
est = Network.Spatial.Simulation (Distance,Beta = params[1],Gamma = params[2],n = 30,times = 250, points = T)

plot (est$S, col="green", type = "l", main= "Implementation of Self-Quarantine ",
      xlab= "Days", ylab = "Density",ylim = c(0,1))
lines(est$I, col="red")
lines(est$R, col="blue")
legend(0.82*length(est$S),0.9,legend = c("Susceptible","Infected","Removed"), col = c(3,2,4),lty=c(1,1,1))






#[1] : https://doi.org/10.1186/s12889-018-5709-x
#[2] : https://stackoverflow.com/questions/29105175/find-neighbouring-elements-of-a-matrix-in-r
#[3] : https://kingaa.github.io/thid/odes/ODEs_in_R.pdf
#[4] : https://kingaa.github.io/clim-dis/parest/parest.html