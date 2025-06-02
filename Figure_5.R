### Diamond Model Community Disassembly Graph
### Author: Joe Brennan

### Figure 5

require(deSolve)
require(igraph)

# A. Create disassembly graph.
# Code adapted from Shoemaker et al, 2020

# ----------------------------------------------------------------------------------------------------
# Function to run model with both species for overall dynamics
VassFox_Cvar <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    
    dP = - (M_P * P) + ( ( (J_P * P) * ( (O_P_C1 * C1) + ( (1 - O_P_C1) * C2) ) ) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) ) 
    dC1 = - (M_C1 * C1) + ( (O_C1_R * J_C1 * C1 * R) / (R + R_0_1) ) - ( (O_P_C1 * J_P * P * C1) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
    dC2 = - (M_C2 * C2) + ( (O_C2_R * J_C2 * C2 * R) / (R + R_0_2) ) - ( ((1-O_P_C1) * J_P * P * C2) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
    dR = r * R * (1 - (R / K)) - ( (O_C1_R * J_C1 * C1 * R) / (R + R_0_1) ) - ( (O_C2_R * J_C2 * C2 * R) / (R + R_0_2) )
    
    return(list(c(dP, dC1, dC2, dR)))
    
  })
}

# ----------------------------------------------------------------------------------------------------
# parameters
# resource intrinsic rate of growth
r = 1.0
# resource carrying capacity
K = 1.0
# consumer 1 ingestion rate
J_C1 = 0.8036
# consumer 2 ingestion rate
J_C2 = 0.7
# predator ingestion rate
J_P = 0.4
# predator mortality rate
M_P = 0.08
# half saturation constant
R_0_1 = 0.16129
R_0_2 = 0.9
C_0 = 0.5
# preference coefficient
O_P_C1 = 0.92
O_C1_R = 1.0
O_C2_R = 0.98
# mortality rate of competitors
M_C1 = 0.4
M_C2 = 0.2 

# create array containing all parameters
pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R)

# Start finding all coexisting communities and evaluating invasion growth rate of each species at each coexisting state

# Empty set
# species abundance at equilibria
P = 0
C1 = 0
C2 = 0
R = 0

# get IGR at respective species abundance
empty.P = - (M_P) + ( ( (J_P) * ( (O_P_C1 * C1) + ( (1 - O_P_C1) * C2) ) ) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) ) 
empty.C1 = - (M_C1) + ( (O_C1_R * J_C1 * R) / (R + R_0_1) ) - ( (O_P_C1 * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
empty.C2 = - (M_C2) + ( (O_C2_R * J_C2 * R) / (R + R_0_2) ) - ( ((1-O_P_C1) * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
empty.R = r * (1 - (R / K)) - ( (O_C1_R * J_C1 * C1) / (R + R_0_1) ) - ( (O_C2_R * J_C2 * C2) / (R + R_0_2) )

empty.IGRs = c(empty.P, empty.C1, empty.C2, empty.R)

####################################################################################

# at R equilibrium. Analytically solvable. R* = K.
# species abundance at equilibria
P = 0
C1 = 0
C2 = 0
R = K

# get IGR at respective species abundance
Req.P = - (M_P) + ( ( (J_P) * ( (O_P_C1 * C1) + ( (1 - O_P_C1) * C2) ) ) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) ) 
Req.C1 = - (M_C1) + ( (O_C1_R * J_C1 * R) / (R + R_0_1) ) - ( (O_P_C1 * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
Req.C2 = - (M_C2) + ( (O_C2_R * J_C2 * R) / (R + R_0_2) ) - ( ((1-O_P_C1) * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
Req.R = r * (1 - (R / K)) - ( (O_C1_R * J_C1 * C1) / (R + R_0_1) ) - ( (O_C2_R * J_C2 * C2) / (R + R_0_2) ) # should equal 0

Req.IGRs = c(Req.P, Req.C1, Req.C2, Req.R)

####################################################################################

# At R & C1 eq. Use simulations to approximate IGRs

# initial conditions
State = c(P=0, C1=1, C2=0, R=1)

# timesteps
time <- seq(0,100000,by=0.1) 

# simulate dynamics
out <- ode(y = State, times = time, func = VassFox_Cvar, parms = pars)

# end of burn in period
invade_start_time <- 10000 #index to start at

# store timeseries output of IGRs
C1Req.P.ts = rep(NA, length(time)-invade_start_time)
C1Req.C1.ts = rep(NA, length(time)-invade_start_time)
C1Req.C2.ts = rep(NA, length(time)-invade_start_time)
C1Req.R.ts = rep(NA, length(time)-invade_start_time)

# keep track of number of calculations
counter = 1

# loop through each time step after burn in period and calculate instantaneous growth rate for each species
for (t in invade_start_time:length(time)){
  
  P = out[t,2]
  C1 = out[t,3]
  C2 = out[t,4]
  R = out[t,5]
  
  C1Req.P.ts[counter] = - (M_P) + ( ( (J_P) * ( (O_P_C1 * C1) + ( (1 - O_P_C1) * C2) ) ) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) ) 
  C1Req.C1.ts[counter] = - (M_C1) + ( (O_C1_R * J_C1 * R) / (R + R_0_1) ) - ( (O_P_C1 * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
  C1Req.C2.ts[counter] = - (M_C2) + ( (O_C2_R * J_C2 * R) / (R + R_0_2) ) - ( ((1-O_P_C1) * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
  C1Req.R.ts[counter] = r * (1 - (R / K)) - ( (O_C1_R * J_C1 * C1) / (R + R_0_1) ) - ( (O_C2_R * J_C2 * C2) / (R + R_0_2) )
  
  counter <- counter + 1
  
}

# average across all instantaneous growth rates to calculate invasion growth rate
C1Req.P = mean(C1Req.P.ts)
C1Req.C1 = mean(C1Req.C1.ts) # should roughly equal 0 since at asymptotic dynamics
C1Req.C2 = mean(C1Req.C2.ts)
C1Req.R = mean(C1Req.R.ts) # should roughly equal 0 since at asymptotic dynamics

C1Req.IGRs = c(C1Req.P, 0, C1Req.C2, 0) # set the ones rougly equal to 0 to zero

####################################################################################

# At R & C2 eq. Use simulations to approximate IGRs

# initial conditions
State = c(P=0, C1=0, C2=1, R=1)

# timesteps
time <- seq(0,100000,by=0.1) 

# simulate dynamics
out <- ode(y = State, times = time, func = VassFox_Cvar, parms = pars)

# burn in period length
invade_start_time <- 10000 #index to start at

# store timeseries output of IGRs
C2Req.P.ts = rep(NA, length(time)-invade_start_time)
C2Req.C1.ts = rep(NA, length(time)-invade_start_time)
C2Req.C2.ts = rep(NA, length(time)-invade_start_time)
C2Req.R.ts = rep(NA, length(time)-invade_start_time)

# track number of iterations through IGR computations
counter = 1

# loop through each time step after burn in period and calculate instantaneous growth rate of each species
for (t in invade_start_time:length(time)){
  
  P = out[t,2]
  C1 = out[t,3]
  C2 = out[t,4]
  R = out[t,5]
  
  C2Req.P.ts[counter] = - (M_P) + ( ( (J_P) * ( (O_P_C1 * C1) + ( (1 - O_P_C1) * C2) ) ) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) ) 
  C2Req.C1.ts[counter] = - (M_C1) + ( (O_C1_R * J_C1 * R) / (R + R_0_1) ) - ( (O_P_C1 * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
  C2Req.C2.ts[counter] = - (M_C2) + ( (O_C2_R * J_C2 * R) / (R + R_0_2) ) - ( ((1-O_P_C1) * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
  C2Req.R.ts[counter] = r * (1 - (R / K)) - ( (O_C1_R * J_C1 * C1) / (R + R_0_1) ) - ( (O_C2_R * J_C2 * C2) / (R + R_0_2) )
  
  counter <- counter + 1
  
}

# average over all instantaneous growth rates for invasion growth rate
C2Req.P = mean(C2Req.P.ts)
C2Req.C1 = mean(C2Req.C1.ts) 
C2Req.C2 = mean(C2Req.C2.ts) # should roughly equal 0 since at asymptotic dynamics
C2Req.R = mean(C2Req.R.ts) # should roughly equal 0 since at asymptotic dynamics

C2Req.IGRs = c(C2Req.P, C2Req.C1, 0, 0) # set the ones rougly equal to 0 to zero

####################################################################################

# At P, C1, R eq. Use simulations to approximate IGRs

# initial conditions
State = c(P=1, C1=1, C2=0, R=1)

# timesteps
time <- seq(0,100000,by=0.1) 

# simulate dynamics
out <- ode(y = State, times = time, func = VassFox_Cvar, parms = pars)

# length of burn in period
invade_start_time <- 10000 #index to start at

# store timeseries output of IGRs
PC1Req.P.ts = rep(NA, length(time)-invade_start_time)
PC1Req.C1.ts = rep(NA, length(time)-invade_start_time)
PC1Req.C2.ts = rep(NA, length(time)-invade_start_time)
PC1Req.R.ts = rep(NA, length(time)-invade_start_time)

# intiialize counter
counter = 1

# loop through each time step after burn in period and calculate instantaneous growth rate for each species 
for (t in invade_start_time:length(time)){
  
  P = out[t,2]
  C1 = out[t,3]
  C2 = out[t,4]
  R = out[t,5]
  
  PC1Req.P.ts[counter] = - (M_P) + ( ( (J_P) * ( (O_P_C1 * C1) + ( (1 - O_P_C1) * C2) ) ) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) ) 
  PC1Req.C1.ts[counter] = - (M_C1) + ( (O_C1_R * J_C1 * R) / (R + R_0_1) ) - ( (O_P_C1 * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
  PC1Req.C2.ts[counter] = - (M_C2) + ( (O_C2_R * J_C2 * R) / (R + R_0_2) ) - ( ((1-O_P_C1) * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
  PC1Req.R.ts[counter] = r * (1 - (R / K)) - ( (O_C1_R * J_C1 * C1) / (R + R_0_1) ) - ( (O_C2_R * J_C2 * C2) / (R + R_0_2) )
  
  counter <- counter + 1
  
}

# Average over the instaneous growth rates to get invasion growth rate. 
PC1Req.P = mean(PC1Req.P.ts) # should roughly equal 0 since at asymptotic dynamics
PC1Req.C1 = mean(PC1Req.C1.ts) # should roughly equal 0 since at asymptotic dynamics
PC1Req.C2 = mean(PC1Req.C2.ts) 
PC1Req.R = mean(PC1Req.R.ts) # should roughly equal 0 since at asymptotic dynamics

PC1Req.IGRs = c(0, 0, PC1Req.C2, 0) # set the ones rougly equal to 0 to zero

####################################################################################

# At P, C1, C2, R eq. Use simulations to approximate IGRs

# initial conditions
State = c(P=1, C1=1, C2=1, R=1)

# timesteps
time <- seq(0,100000,by=0.1) 

# simulate dynamics
out <- ode(y = State, times = time, func = VassFox_Cvar, parms = pars)

# length of burn in period
invade_start_time <- 10000 #index to start at

# store timeseries output of IGRs
PC1C2Req.P.ts = rep(NA, length(time)-invade_start_time)
PC1C2Req.C1.ts = rep(NA, length(time)-invade_start_time)
PC1C2Req.C2.ts = rep(NA, length(time)-invade_start_time)
PC1C2Req.R.ts = rep(NA, length(time)-invade_start_time)

# initialize counter
counter = 1

# loop through each time step after burn in period and calculate instantaneous growth rate for each species
for (t in invade_start_time:length(time)){
  
  P = out[t,2]
  C1 = out[t,3]
  C2 = out[t,4]
  R = out[t,5]
  
  PC1C2Req.P.ts[counter] = - (M_P) + ( ( (J_P) * ( (O_P_C1 * C1) + ( (1 - O_P_C1) * C2) ) ) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) ) 
  PC1C2Req.C1.ts[counter] = - (M_C1) + ( (O_C1_R * J_C1 * R) / (R + R_0_1) ) - ( (O_P_C1 * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
  PC1C2Req.C2.ts[counter] = - (M_C2) + ( (O_C2_R * J_C2 * R) / (R + R_0_2) ) - ( ((1-O_P_C1) * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
  PC1C2Req.R.ts[counter] = r * (1 - (R / K)) - ( (O_C1_R * J_C1 * C1) / (R + R_0_1) ) - ( (O_C2_R * J_C2 * C2) / (R + R_0_2) )
  
  counter <- counter + 1
  
}

# average over instantaneous growth rates to calculate invasion growth rate
PC1C2Req.P = mean(PC1C2Req.P.ts) # should roughly equal 0 since at asymptotic dynamics
PC1C2Req.C1 = mean(PC1C2Req.C1.ts) # should roughly equal 0 since at asymptotic dynamics
PC1C2Req.C2 = mean(PC1C2Req.C2.ts) # should roughly equal 0 since at asymptotic dynamics
PC1C2Req.R = mean(PC1C2Req.R.ts) # should roughly equal 0 since at asymptotic dynamics

PC1C2Req.IGRs = c(0, 0, 0, 0) # set the ones roughly equal to 0 to zero

############################################################################

# Put all these into a matrix - invasion scheme
Invasion.Scheme = rbind(empty.IGRs, Req.IGRs, C1Req.IGRs, C2Req.IGRs, PC1Req.IGRs, PC1C2Req.IGRs)

# function to calculate the disassembly graph
Disassembly.graph = function(IS,spp.names=c()){ 
  
  n=dim(IS)[1] # number of communities
  temp.composition=list() # a list to hold the composition of the communities
  composition=list()
  edge.label = c() # list of edge labels detailing which invader leads to transition
  for(i in 1:n)temp.composition[[i]]=which(IS[i,]==0)
  for(i in 1:n)composition[[i]]=spp.names[temp.composition[[i]]]
  
  
  # Create matrix of zeroes, will be adjacency matrix for disassembly graph
  mat = matrix(data = 0, nrow = nrow(IS), ncol = nrow(IS))
  label.edge.mat = mat
  
  # Loop through each community
  for (i in 1:nrow(IS)){
    
    # extract community composition for community i
    temp.community = which(IS[i,]==0)
    
    # Identify which community are subsets of S
    #index = which(rowSums(IS[-i,setdiff(1:ncol(IS), temp.community)]!=0)==length(setdiff(1:ncol(IS), temp.community)))
    if (length(setdiff(1:ncol(IS), temp.community))>1) candidates = which(rowSums(IS[-i,setdiff(1:ncol(IS), temp.community)]!=0)==length(setdiff(1:ncol(IS), temp.community)))
    if (length(setdiff(1:ncol(IS), temp.community))==1) candidates = which((IS[-i,setdiff(1:ncol(IS), temp.community)]!=0)==length(setdiff(1:ncol(IS), temp.community)))
    if ((i == nrow(IS)) & (all(IS[nrow(IS),]==rep(0,ncol(IS))))) candidates = 1:nrow(IS)
    
    # loop through each candidate and ensure it is a possible disassembled community
    for (j in 1:length(candidates)){
      
      # Extract community composition of the j-th candidate community
      T.temp = which(IS[candidates[j],]==0)
      
      # consider case with only one species difference
      if (length(setdiff(temp.community,T.temp))==1){
        mat[i,candidates[j]] = 1
        label.edge.mat[i,candidates[j]] = setdiff(temp.community,T.temp)
      }
      
      #consider case with multiple species difference
      if (length(setdiff(temp.community, T.temp))>1 & sum(IS[candidates[j], setdiff(temp.community, T.temp)]>0)==1){
        
        mat[i,candidates[j]] = 1
        label.edge.mat[i,candidates[j]] = setdiff(temp.community, T.temp)[which(IS[candidates[j], setdiff(temp.community, T.temp)]>0)]
        
      }
      
      # consider case where multiple species difference but all neg IGR
      if (length(setdiff(temp.community, T.temp))>1 & sum(IS[candidates[j],  setdiff(temp.community, T.temp)]>0)==0){
        
        mat[i, candidates[j]] = 1
        label.edge.mat[i,candidates[j]] = paste(setdiff(temp.community, T.temp), collapse = ',')
        
      }
      
    }
    
  }
  
  # turn edge label matrix into array for graphing
  edge.label = c()
  
  # add each row in matrix into single array
  for (i in 1:nrow(label.edge.mat)){
    edge.label = c(edge.label, label.edge.mat[i,])
  }
  
  # only include entries that are actual labels, i.e. all non-zero numbers
  edge.label = spp.names[edge.label]
  
  return(list(adj.mat=mat,edge.label = edge.label, IS = IS, composition=composition))
  
}

# Create disassembly graph for this specific model
out = Disassembly.graph(Invasion.Scheme, spp.names = c('P', 'C1', 'C2', 'R'))

# plotting function for disassembly graph
plot.disassembly.graph=function(out,clear.out=c(),bend.factor=0.75,cols,impermanent.weight=0.172549,multiple.invasion.edge.weight=0.5,vlx=1,elx=1){
  
  IG=out$adj.mat
  IS=out$IS
  composition = out$composition
  number.species = 1:ncol(IS)
  out$edge.label <- paste0("-", out$edge.label)
  
  # create graph
  g=graph_from_adjacency_matrix(IG)
  # get edge vertex info: column 1 = tail, 2 = tip
  edges_data_frame <- get.data.frame(g, what = "edges")
  
  # Compute the number of species in each community
  community_sizes <- sapply(out$composition, length)
  
  # Find edges where species richness decreases by more than 1
  edges_large_drop <- edges_data_frame[community_sizes[edges_data_frame$from] - community_sizes[edges_data_frame$to] > 1, ]
  edges_large_drop$label = out$edge.label[as.numeric(rownames(edges_large_drop))]
  edges_large_drop$fromcomp = out$composition[edges_large_drop$from]
  edges_large_drop$tocomp = out$composition[edges_large_drop$to]
  
  index = setdiff(1:length(out$edge.label), as.numeric(rownames(edges_large_drop)))
  
  out$edge.label = '' # no edge labels for this
  
  # create edge widths based on one versus multiple invasion
  for(i in 1:nrow(edges_data_frame)){
    E(g)$weight[i] = multiple.invasion.edge.weight
    E(g)$color[i] = "gray"  # Default to gray (single species loss)
    
    SS = composition[[edges_data_frame[i, 1]]]
    TT = composition[[edges_data_frame[i, 2]]]
    
    # If edge represents a multiple species loss, make it yellow
    if (any(edges_large_drop$from == edges_data_frame[i, 1] & edges_large_drop$to == edges_data_frame[i, 2])) {
      E(g)$color[i] = cols[2]  # Yellow for multiple loss
    }
    
    if (is.element(edges_data_frame[i, 2], clear.out)) {
      E(g)$color[i] = NA
    }
  }
  # create the reordered version of the graph (by edge weight that is used ultimately for the plotting
  g2 <- make_graph(as.vector(t(get.edgelist(g)[order(E(g)$weight),])))
  # Compute the number of species in each community
  community_sizes <- sapply(composition, length)
  
  # Order nodes by number of species
  sorted_indices <- order(community_sizes)
  
  # Recompute xvals and yvals
  k <- length(composition)
  xvals <- numeric(k)
  yvals <- numeric(k)
  
  lengths <- table(community_sizes)  # Counts of each community size
  max_length <- max(lengths)
  
  counter <- 1
  for (i in 0:max(community_sizes)) {
    num_in_level <- lengths[as.character(i)]
    if (is.na(num_in_level)) next
    
    for (j in 1:num_in_level) {
      xvals[counter] <- j + (max_length - num_in_level) / 2  # Center alignment
      x.temp <- 2 * xvals[counter] / max_length - 1
      yvals[counter] <- i * 1.2  
      # Height based on species count
      counter <- counter + 1
    }
  }
  
  # move 1,2 to left a little
  xvals[5] = xvals[5]-0.17
  yvals[5] = yvals[5]-0.25
  
  # community names
  v.name=c(expression(symbol("\306")))
  for(i in 2:k){
    temp=composition[[i]]
    v.name=c(v.name,paste(unlist(temp),collapse=", "))
  }
  for(i in clear.out)v.name[i]=""
  
  V(g)$name=v.name
  
  vertex.frame.cols=rep(cols[6],k)
  vertex.label.cols=rep(cols[6],k)
  vertex.shapes=rep("circle",k)
  vcols = rep('white', k)
  for(i in clear.out){
    vcols[i]=NA
    vertex.frame.cols[i]=NA
    vertex.label.cols[i]=cols[5]
    for(j in which(edges_data_frame$to==i)){E(g)$color[j]=NA}
  }
  
  # assign properties to the reordered graph g2
  E(g2)$color <- E(g)$color[order(E(g)$weight)]
  E(g2)$weight<-E(g)$weight[order(E(g)$weight)]
  E(g2)$weight[as.numeric(rownames(edges_large_drop))] = 5
  
  V(g2)$name <-v.name
  par(mar=c(0,0,0,0))
  plot(g2,vertex.color=vcols,edge.width=E(g2)$weight,
       vertex.size = 25,
       vertex.frame.color=vertex.frame.cols,vertex.shape=vertex.shapes,
       vertex.label.color=vertex.label.cols,layout=cbind(xvals,yvals),
       vertex.label.cex=vlx, edge.label.font=2, edge.arrow.size=0.9, edge.label = out$edge.label,edge.label.cex=elx,edge.label.color=cols[9], edge.arrow.size = 1.5) 
}

cols = c('white', '#FFDB00', 'white', 'white', 'white', 'black', '#FFDB00', 'black')
plot.disassembly.graph(out, cols=cols)

# ----------------------------------------------------------------------------------------------------

# B. Simulate loss of P from full community

# Model
VassFox_Cvar <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    
    dP = - (M_P * P) + ( ( (J_P * P) * ( (O_P_C1 * C1) + ( (1 - O_P_C1) * C2) ) ) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) ) 
    dC1 = - (M_C1 * C1) + ( (O_C1_R * J_C1 * C1 * R) / (R + R_0_1) ) - ( (O_P_C1 * J_P * P * C1) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
    dC2 = - (M_C2 * C2) + ( (O_C2_R * J_C2 * C2 * R) / (R + R_0_2) ) - ( ((1-O_P_C1) * J_P * P * C2) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
    dR = r * R * (1 - (R / K)) - ( (O_C1_R * J_C1 * C1 * R) / (R + R_0_1) ) - ( (O_C2_R * J_C2 * C2 * R) / (R + R_0_2) )
    
    return(list(c(dP, dC1, dC2, dR)))
    
  })
}

# ----------------------------------------------------------------------------------------------------
# parameters
# resource intrinsic rate of growth
r = 1.0
# resource carrying capacity
K = 1.0
# consumer 1 ingestion rate
J_C1 = 0.8036
# consumer 2 ingestion rate
J_C2 = 0.7
# predator ingestion rate
J_P = 0.4
# predator mortality rate
M_P = 0.08
# half saturation constant
R_0_1 = 0.16129
R_0_2 = 0.9
C_0 = 0.5
# preference coefficient
O_P_C1 = 0.92
O_C1_R = 1.0
O_C2_R = 0.98
# mortality rate of competitors
M_C1 = 0.4
M_C2 = 0.2 

# put parameters into a list
pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R)

# initial conditions
State = c(P=1, C1=1, C2=1, R=1)

# timesteps
time <- seq(0,500,by=0.1) 

# simulate dynamics
out <- ode(y = State, times = time, func = VassFox_Cvar, parms = pars)

# Make P extinct
out[nrow(out),2]=0 

# continue simulating under new state
newout = ode(y = out[nrow(out),2:5], times = seq(500,1000, by=0.1), func = VassFox_Cvar, parms = pars)

# merge time series
new.out = rbind(out, newout)
colnames(new.out) = c('Time', 'P', 'C_1', 'C_2', 'R')

# read in libraries needed
require(tidyverse)
require(scales)

# create dataframe for plotting
df <- new.out[3000:7000,] %>%
  as.data.frame() %>%
  pivot_longer(cols = -Time, names_to = "Species", values_to = "Abundance")%>%
  mutate(Time = rescale(Time, to = c(0, 400), from = c(300, 700)))

# plot
ggplot(df, aes(x = Time, y = Abundance, color = Species)) +
  scale_color_manual(
    values = c("P" = "orange", "C_1" = "blue", "C_2" = "black", "R" = "purple"),
    labels = c(expression("C"[1]), expression("C"[2]), 'P', 'R')
  ) +
  geom_line(size = 1.2) +
  theme_classic() +
  ylim(c(0,1.5)) +
  labs(x = "Time",
       y = "Density") +
  theme(legend.position = c(0.85, 0.85),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.background = element_rect(color = "black", fill = NA),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 14),  
        axis.text.y = element_text(size = 14))

# ----------------------------------------------------------------------------------------------------

# C. Simulate loss of C2 from full community

# initial conditions
State = c(P=1, C1=1, C2=1, R=1)

# timesteps
time <- seq(0,500,by=0.1) 

# simulate dynamics
out <- ode(y = State, times = time, func = VassFox_Cvar, parms = pars)

# Make P extinct
out[nrow(out),4]=0 

# simulate under new state
newout = ode(y = out[nrow(out),2:5], times = seq(500,1000, by=0.1), func = VassFox_Cvar, parms = pars)

# merge time series
new.out = rbind(out, newout)
colnames(new.out) = c('Time', 'P', 'C_1', 'C_2', 'R')

# create dataframe for plotting
df <- new.out[3000:7000,] %>%
  as.data.frame() %>%
  pivot_longer(cols = -Time, names_to = "Species", values_to = "Abundance")%>%
  mutate(Time = rescale(Time, to = c(0, 400), from = c(300, 700)))

# plot
ggplot(df, aes(x = Time, y = Abundance, color = Species)) +
  scale_color_manual(
    values = c("P" = "orange", "C_1" = "blue", "C_2" = "black", "R" = "purple"),
    labels = c(expression("C"[1]), expression("C"[2]), 'P', 'R')
  ) +
  geom_line(size = 1.2) +
  theme_classic() +
  ylim(c(0,1.5)) +
  labs(x = "Time",
       y = "Density") +
  theme(legend.position = c(0.85, 0.85),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.background = element_rect(color = "black", fill = NA),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 14),  
        axis.text.y = element_text(size = 14))
