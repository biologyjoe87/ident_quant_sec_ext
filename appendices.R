### Code to create supplementary figures
### Author: Joe Brennan

library(igraph)
library(gRbase)
library(viridis)
library(deSolve)
library(tidyverse)

#################################################################
### GENERAL FUNCTIONS

# Create Invasion Graph
Invasion.Graph.Function = function(IS) {
  n = dim(IS)[1]
  composition = list()
  richness = numeric(n)
  
  # Determine community compositions
  for(i in 1:n) {
    temp = which(abs(IS[i,]) < tolerance)
    composition[[i]] = temp
    richness[i] = length(temp)
  }
  
  # Build adjacency matrix
  IG = matrix(0, n, n)
  for(i in 1:n) {
    for(j in 1:n) {
      if(i != j) {
        b = composition[[j]]
        a = composition[[i]]
        c = setdiff(b, a)
        c2 = setdiff(a, b)
        
        condition.1 = 1
        if(length(c) > 0) condition.1 = min(IS[i,c])
        
        condition.2 = -1
        if(length(c2) > 0) condition.2 = max(IS[j,c2])
        
        if((condition.1 > tolerance) & (condition.2 < -tolerance)) {
          IG[i,j] = 1
        }
      }
    }
  }
  
  return(list(
    IG = IG,
    richness = richness,
    composition = composition,
    IS = IS
  ))
}

# analyze invasion graph for annual plant
Invasion.Graph.Analysis.Function = function(IG.list) {
  n = dim(IG.list$IS)[1]
  composition = IG.list$composition
  richness = IG.list$richness
  IG = IG.list$IG
  IS = IG.list$IS
  
  # Check necessary permanence conditions and acyclic conditions for all communities. Create the partial.ordering.matrix for the communities along the way
  permanence.condition = matrix(TRUE, n, 1)
  acyclic = matrix(FALSE, n, 1)
  partial.ordering.matrix=IG*0
  for(i in 1:n) {
    temp.commmunity.list=c() # keep track of the communities that are subset of community i
    for(j in 1:n) {
      if((j != i) && is_subsetof(composition[[j]], composition[[i]])) {
        temp.commmunity.list=c(temp.commmunity.list,j)
        partial.ordering.matrix[i,j]=1 # j is a subset i 
        if((max(IS[j,composition[[i]]]) <= tolerance)) {
          permanence.condition[i] = FALSE
        }
      }
    }
    IG.temp=IG[temp.commmunity.list,temp.commmunity.list]
    if(length(IG.temp)<2){acyclic[i]=TRUE}else{
      g.temp = graph_from_adjacency_matrix(IG[temp.commmunity.list,temp.commmunity.list]) 
      acyclic[i] = is_dag(g.temp)}
  }
  
  # Find minus-i communities
  minus.i = list()
  k = dim(IS)[2]
  for(i in 1:k) {
    temp = c()
    for(j in 1:n) {
      if((!is.element(i, composition[[j]])) && (max(IS[j,-i]) <= 0)) {
        temp = c(temp, j)
      }
    }
    minus.i[[i]] = temp
  }
  
  return(list(
    IG = IG,
    acyclic = acyclic,
    permanence.condition = permanence.condition,
    minus.i = minus.i,
    richness = richness,
    composition = composition,
    IS = IS,
    partial.ordering.matrix=partial.ordering.matrix
  ))
}

# graph plotter
Graph_Plotter <- function(input, highlight_edges = NULL, highlight_color = "red", highlight_width = 4,vertex_size=30,vertex_label_cex=1,default_edge_width=2,default_edge_color="darkgray",edge_arrow_size=0.5,color_palette = "viridis",arc_height=0.75,vertex_label_col="black",vertex_frame_color="black") {
  # Create graph from adjacency matrix
  g <- graph_from_adjacency_matrix(input$IG, mode = "directed")
  
  # Get unique richness values and sort them
  unique_richness <- sort(unique(input$richness))
  
  # Create base y-coordinates based on richness levels
  y_coords <- match(input$richness, unique_richness)
  
  # Create x-coordinates and adjust y-coordinates with arc
  x_coords <- numeric(length(input$richness))
  
  # Arc height parameter - adjust this to control the curve
  arc_height <- arc_height  # Maximum deviation from the baseline
  
  for (level in unique_richness) {
    # Find vertices at this richness level
    vertices_at_level <- which(input$richness == level)
    n_vertices <- length(vertices_at_level)
    
    if (n_vertices >1) {
      # Space vertices equally between 0 and 1
      level_x_coords <- seq(0, 1, length.out = n_vertices)
      
      # Create upward parabolic arc adjustment for y-coordinates
      arc_adjustment <- -arc_height * (4 * level_x_coords * (1 - level_x_coords))
      
      # Adjust y-coordinates for vertices at this level
      y_coords[vertices_at_level] <- y_coords[vertices_at_level] + arc_adjustment
      
    } else {
      # If only one vertex at this level, put it in the center with no arc
      level_x_coords <- 0.5
      y_coords[vertices_at_level] <- y_coords[vertices_at_level]- arc_height
      
    }
    
    # Assign x coordinates to the vertices at this level
    x_coords[vertices_at_level] <- level_x_coords
  }
  
  # Combine into layout matrix
  layout_matrix <- cbind(x_coords, y_coords)
  
  # Normalize heights to [0,1] for viridis color mapping
  normalized_heights <- (input$height - min(input$height)) / 
    (max(input$height) - min(input$height))
  
  # extend highlight_colors to a vector if needed
  if(length(highlight_color)==1)highlight_color=rep(highlight_color,length(highlight_edges))
  
  # Set up edge attributes for highlighting
  edge_colors <- rep(default_edge_color, ecount(g))
  edge_widths <- rep(default_edge_width, ecount(g))
  
  # If highlight_edges is provided, update edge attributes
  if (!is.null(highlight_edges)) {
    # highlight_edges should be a list of vertex index pairs
    temp.counter=1
    for(edge in highlight_edges) {
      # Find the edge index
      edge_id <- get.edge.ids(g, edge)
      if(edge_id > 0) {  # if edge exists
        edge_colors[edge_id] <- highlight_color[temp.counter]
        edge_widths[edge_id] <- highlight_width
        temp.counter=temp.counter+1
      }
    }
  }
  # Create the plot
  par(mar=c(0,0,0,0))
  plot(g,
       layout = layout_matrix,
       vertex.color = rep('white', nrow(input$IG)),
       vertex.size = vertex_size,
       vertex.label = input$vertex.labels,
       vertex.label.color = vertex_label_col,
       vertex.frame.color=vertex_frame_color,
       vertex.label.cex = vertex_label_cex,
       vertex.label.dist = 0,
       vertex.label.degree = pi/2,
       edge.arrow.size = edge_arrow_size,
       edge.color = edge_colors,
       edge.width = edge_widths,
       main = "")
}

condensed_directed_graph <- function(adj_matrix) {
  # Find strongly connected components
  components <- find_equivalence_classes(adj_matrix)
  
  # Add singleton vertices (those not in any cycle)
  n <- nrow(adj_matrix)
  vertices_in_components <- unique(unlist(components))
  singletons <- setdiff(1:n, vertices_in_components)
  components <- c(components, as.list(singletons))
  
  # Create new adjacency matrix
  n_new <- length(components)
  condensed_adj <- matrix(0, nrow=n_new, ncol=n_new)
  
  # For each pair of components
  for(i in 1:n_new) {
    for(j in 1:n_new) {
      if(i != j) {
        # Check if there's any edge from component i to component j
        for(v1 in components[[i]]) {
          for(v2 in components[[j]]) {
            if(adj_matrix[v1,v2] > 0) {
              condensed_adj[i,j] <- 1
              break
            }
          }
          if(condensed_adj[i,j] > 0) break
        }
      }
    }
  }
  
  return(list(
    condensed_graph = condensed_adj,
    component_map = components
  ))
}

find_equivalence_classes <- function(adj_matrix) {
  n <- nrow(adj_matrix)
  index <- 0
  stack <- integer(0)
  on_stack <- logical(n)
  indices <- integer(n)
  lowlink <- integer(n)
  components <- list()
  
  tarjan <- function(v) {
    index <<- index + 1
    indices[v] <<- index
    lowlink[v] <<- index
    stack <<- c(stack, v)
    on_stack[v] <<- TRUE
    
    # Consider successors of v
    for(w in which(adj_matrix[v,] > 0)) {
      if(indices[w] == 0) {  # Successor w has not yet been visited
        tarjan(w)
        lowlink[v] <<- min(lowlink[v], lowlink[w])
      } else if(on_stack[w]) {  # Successor w is in stack and hence in the current SCC
        lowlink[v] <<- min(lowlink[v], indices[w])
      }
    }
    
    # If v is a root node, pop the stack and generate an SCC
    if(lowlink[v] == indices[v]) {
      component <- integer(0)
      repeat {
        w <- stack[length(stack)]
        stack <<- stack[-length(stack)]
        on_stack[w] <<- FALSE
        component <- c(component, w)
        if(w == v) break
      }
      if(length(component) > 1) {  # Only add components with more than one vertex
        components <<- c(components, list(sort(component)))
      }
    }
  }
  
  # Initialize data structures
  indices[] <- 0
  
  # Find SCCs from each unvisited vertex
  for(v in 1:n) {
    if(indices[v] == 0) {
      tarjan(v)
    }
  }
  
  return(components)
}

my_topo_sort = function(Adj) {
  # Get number of vertices
  n = dim(Adj)[1]
  
  # Initialize ordered list for vertices
  L = c()
  
  # Initialize counter for level sets
  counter = 0
  
  # Initialize vector for level set values
  V = c()
  
  # Initialize set of remaining vertices
  remaining = 1:n
  
  # Process vertices until none remain
  while(length(remaining > 0)) {
    # Find vertices with no incoming edges in remaining subgraph
    if(length(remaining) > 1) {
      S = which(colSums(Adj[remaining, remaining]) == 0)
    } else {
      S = 1
    }
    
    # Add found vertices to ordered list
    L = c(L, remaining[S])
    
    # Assign current level to these vertices
    V = c(V, rep(counter, length(S)))
    
    # Increment level counter
    counter = counter + 1
    
    # Remove processed vertices from remaining set
    remaining = setdiff(remaining, remaining[S])
    length(remaining)  # Note: This value isn't used
  }
  
  return(list(V=V, L=L))
}

map_condensed_to_original <- function(values, component_map) {
  # Check input validity
  if(length(values) != length(component_map)) {
    stop("Length of values must match number of components")
  }
  
  # Get number of vertices in original graph
  n_original <- max(unlist(component_map))
  
  # Create output vector
  original_values <- numeric(n_original)
  
  # Map each value to original vertices
  for(i in seq_along(component_map)) {
    original_values[component_map[[i]]] <- values[i]
  }
  
  return(original_values)
}

###################################################################
### INVASION GRAPH FOR ANNUAL PLANT COMMUNITY

# Create matrix of invasion growth rates for the annual plant model
plant.IS = function(A,g,s,lambda, tolerance=1e-14){
  
  k=dim(A)[1] # number of species
  C=list() # list to hold all the feasible (i.e. non-negative entries) equilibria of the model
  C[[1]]=rep(0,k) # the origin as the first equilibrium Note: for replicator equation will need the first k equilibria to be the vertices (1,0,0,..0), (0,1,0,...,0) etc
  no.C=1 # counter to keep track of the number of feasible equilibria
  # find all the feasible equilibria
  for(i in 1:k){ #loop for the number of species supported by the equilibrium. Note: for replicator this will start with 2 instead of 1
    temp=combn(1:k,i) # a matrix of all the possible configurations of communities with i species
    k2=dim(temp)[2] # the number of configurations
    for(j in 1:k2){ # loop through the configurations
      I=temp[,j] # pull out the j-th configuration
      xtemp=solve(a = A[I,I],b= lambda[I] * g[I] / (1 - s[I]*(1-g[I]))-1, tol = tolerance)/g[I] # solve for the equilibrium restricted to the species in I
      if(min(xtemp)>0){ # make sure all entries are positive i.e. a new feasible equilibrium
        no.C=no.C+1 # update counter
        xtemp2=C[[1]] # grab the zero vector to fill out the feasible equilibrium for all k species Note: this should be changed to the zero vector for the replicator
        xtemp2[I]=xtemp # fill in the non-zero entries
        C[[no.C]]=xtemp2 # set as the next element in the list
      }
    }
  }
  # create the invasion scheme matrix
  IS=matrix(NA,length(C),k) # matrix to hold all the per-capita growth rates
  for(i in 1:length(C)){
    IS[i,]=log(s*(1-g) + g*lambda /(1 + ((A %*% (g*C[[i]]))))) # i-th row corresponds to all the per-capita growth rates at the i-th feasible equilibrium
  }
  IS[which(abs(IS)<tolerance)]=0 # set to zero based on the tolerance
  return(IS)
  
}

# Empirically-parametrized values from Van Dyke et al
A = rbind(c(0.27290462, 0.9415389, 0.5590646),
          c(0.13755487, 0.9483990,0.1599514),
          c(0.08483079, 1.0221670, 0.4161495))
g = c(0.2240310, 0.6666667, 0.6029647)
s = c(0.2625463, 0.0450000, 0.6215016)
lambda = c(2351.5366, 745.9973, 654.6758)

# Species names
abbr.spp = c('A. wrangelianus', 'H. murinum', 'P. erecta')

# Get the invasion scheme - matrix of invasion growth rates
tolerance = 1e-8

IS = plant.IS(A,g,s,lambda)
IG.out = Invasion.Graph.Function(IS)
finalout = Invasion.Graph.Analysis.Function(IG.out)

input=c()
input$IG=finalout$IG
input$richness=finalout$richness
input$vertex.labels=sapply(finalout$composition, paste, collapse = ",")
composition=finalout$composition

result <- condensed_directed_graph(input$IG)

topo_temp<-my_topo_sort(result$condensed_graph)
condensed_values<-(topo_temp$V)
condensed_values[topo_temp$L] <-(topo_temp$V)
original_values <- map_condensed_to_original(condensed_values, result$component_map)
input$height<-original_values

Graph_Plotter(input,
              highlight_width = 3,vertex_size=17,vertex_label_cex = 0.75, color_palette = 'terrain')

##########################################################################################
### INVASION GRAPH FOR GRASSLAND SYSTEM

# empirically parameterized values from Geijzendorffer et al
spring.b = c(100.77, 45.19, 103.43, 67.09, 112.71, 76.52)/1000
spring.A = -rbind(c(1.28, -0.15, -0.10, -0.22, 0.16, -0.12),
                  c(-0.13, 0.41, 0.12, 0.09, -0.03, 0),
                  c(0.38, 0.30, 0.41, 0.36, 0.15, 0.11),
                  c(-0.95, 0.07, 0.36, 0.31, -3.6, 0.34),
                  c(0.39, 0.15, 0.82, 0.72, 0.8, -0.78),
                  c(0.02, 0.18, 0.05, 0.13, 0.57, 0.56))

# only focus on subset of the species
interest.spp = c(1,2,3,5,6)

A = spring.A[interest.spp, interest.spp]
b = spring.b[interest.spp]

LV.Invasion.Scheme = function(A, b, tolerance = 1e-8) {
  k = dim(A)[1] # number of columns is number of species
  C = list() # initialize list of coexisting commmunities
  C[[1]] = rep(0, k) # include the empty set
  no.C = 1
  
  # Find all non-negative equilibria
  for(i in 1:k) {
    temp = combn(1:k, i)
    k2 = dim(temp)[2]
    for(j in 1:k2) {
      I = temp[,j]
      xtemp = solve(a = A[I,I], b = -b[I])
      if(min(xtemp) > tolerance) {
        no.C = no.C + 1
        xtemp2 = C[[1]]
        xtemp2[I] = xtemp
        C[[no.C]] = xtemp2
      }
    }
  }
  
  # Calculate per-capita growth rates
  IS = matrix(NA, length(C), k)
  for(i in 1:length(C)) {
    IS[i,] = A %*% C[[i]] + b
  }
  IS[which(abs(IS) < tolerance)] = 0
  
  return(IS)
}

# Create invasion scheme for particular empirically derived model
IS = LV.Invasion.Scheme(A,b)

IG.out = Invasion.Graph.Function(IS)
finalout = Invasion.Graph.Analysis.Function(IG.out)

input=c()
input$IG=finalout$IG
input$richness=finalout$richness
input$vertex.labels=sapply(finalout$composition, paste, collapse = ",")
composition=finalout$composition

result <- condensed_directed_graph(input$IG)

topo_temp<-my_topo_sort(result$condensed_graph)
condensed_values<-(topo_temp$V)
condensed_values[topo_temp$L] <-(topo_temp$V)
original_values <- map_condensed_to_original(condensed_values, result$component_map)
input$height<-original_values

Graph_Plotter(input,
              highlight_width = 3,vertex_size=17,vertex_label_cex = 0.75, color_palette = 'terrain')

####################################################################################
### INVASION GRAPH FOR DIAMOND MODEL

# ----------------------------------------------------------------------------------------------------
# Function to run model with both species for overall dynamics, code from Shoemaker et al, 2020
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

# Put all these into a matrix - invasion scheme
Invasion.Scheme = rbind(empty.IGRs, Req.IGRs, C1Req.IGRs, C2Req.IGRs, PC1Req.IGRs, PC1C2Req.IGRs)

IG.out = Invasion.Graph.Function(Invasion.Scheme)
finalout = Invasion.Graph.Analysis.Function(IG.out)

input=c()
input$IG=finalout$IG
input$richness=finalout$richness
input$vertex.labels=sapply(finalout$composition, paste, collapse = ",")
composition=finalout$composition

result <- condensed_directed_graph(input$IG)

topo_temp<-my_topo_sort(result$condensed_graph)
condensed_values<-(topo_temp$V)
condensed_values[topo_temp$L] <-(topo_temp$V)
original_values <- map_condensed_to_original(condensed_values, result$component_map)
input$height<-original_values

Graph_Plotter(input,
              highlight_width = 3,vertex_size=17,vertex_label_cex = 0.75, color_palette = 'terrain')

############################################################################################
### SIMULATE ANNUAL PLANT DYNAMICS

# Empirically-parametrized values from Van Dyke et al
A = rbind(c(0.27290462, 0.9415389, 0.5590646),
          c(0.13755487, 0.9483990,0.1599514),
          c(0.08483079, 1.0221670, 0.4161495))
g = c(0.2240310, 0.6666667, 0.6029647)
s = c(0.2625463, 0.0450000, 0.6215016)
lambda = c(2351.5366, 745.9973, 654.6758)

# Set initial conditions
init = c(3000, 500, 190)

# length of time
iterations = 150

# initialize matrix to contain the time series
mat = matrix(NA, nrow=iterations, ncol=length(init))
mat[1,]=init

# Simulate
for (i in 2:nrow(mat)){
  mat[i,] = s*(1-g)*mat[i-1,] + lambda*g*mat[i-1,] / (1 + A %*% (g*mat[i-1,]))
  if (i==75) mat[i,2]=0 # extinction event
}

# prepare simulated values for plotting
mat = cbind(mat, 1:iterations)
colnames(mat) = c("A. wrangelianus", "H. murinum", "P. erecta", 'Iterations')

# create dataframe for plotting
df <- mat %>%
  as.data.frame() %>%
  pivot_longer(cols = -Iterations, names_to = "Species", values_to = "Abundance")

# plot
ggplot(df, aes(x = Iterations, y = Abundance, color = Species)) +
  scale_color_manual(
    values = c("A. wrangelianus" = "orange", "H. murinum" = "blue", "P. erecta" = "black")) +
  geom_line(size = 1.2) +
  theme_classic() +
  labs(x = "Time",
       y = "Density") +
  theme(legend.position = c(0.8, 0.85),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.background = element_rect(color = "black", fill = NA),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 14),  
        axis.text.y = element_text(size = 14))

###########

# Set initial conditions
init = c(0, 0, 2000)

# length of time
iterations = 600

# initialize matrix to contain the time series
mat = matrix(NA, nrow=iterations, ncol=length(init))
mat[1,]=init

# Simulate
for (i in 2:nrow(mat)){
  mat[i,] = s*(1-g)*mat[i-1,] + lambda*g*mat[i-1,] / (1 + A %*% (g*mat[i-1,]))
  if (i==50) mat[i,2]=1 # extinction event
  if (i==250) mat[i,1]=1
}

# prepare simulated values for plotting
mat = cbind(mat, 1:iterations)
colnames(mat) = c("A. wrangelianus", "H. murinum", "P. erecta", 'Iterations')

# create dataframe for plotting
df <- mat %>%
  as.data.frame() %>%
  pivot_longer(cols = -Iterations, names_to = "Species", values_to = "Abundance")

# plot
ggplot(df, aes(x = Iterations, y = Abundance, color = Species)) +
  scale_color_manual(
    values = c("A. wrangelianus" = "orange", "H. murinum" = "blue", "P. erecta" = "black")) +
  geom_line(size = 1.2) +
  theme_classic() +
  labs(x = "Time",
       y = "Density") +
  theme(legend.position = c(0.8, 0.85),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.background = element_rect(color = "black", fill = NA),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 14),  
        axis.text.y = element_text(size = 14))

##############################################################################
### MCT DECOMPOSITION OF OTHER SECONDARY EXTINCTION IN ANNUAL PLANT MODEL

# empirically parameterized values from Van Dyke et al
A = rbind(c(0.27290462, 0.9415389, 0.5590646),
          c(0.13755487, 0.9483990,0.1599514),
          c(0.08483079, 1.0221670, 0.4161495))
g = c(0.2240310, 0.6666667, 0.6029647)
s = c(0.2625463, 0.0450000, 0.6215016)
lambda = c(2351.5366, 745.9973, 654.6758)

# invasion growth rate for annual plant model
IGR.plant <- function(A,s,g,lambda,res.eq){ # function to get IGR for gLV
  
  IGR.row = log(s*(1-g) + g*lambda /(1 + ((A %*% (g*res.eq)))))
  
  return(IGR.row)
  
}

# specify resident community and invading species
res = 3
inv = 1

# get equilibrium density for original function given resident community
res.eq = numeric(nrow(A))
res.eq[res] = solve(a = A[res,res],b= lambda[res] * g[res] / (1 - s[res]*(1-g[res]))-1)/g[res]

# get Invasion Growth Rate of specific model
orig.IGR = IGR.plant(A,s,g,lambda,res.eq)[inv]

# set interactions of interest to 0. In our case, competition between species.
temp.A=A
temp.A[1,3]<-0

# find baseline value i.e., invasion growth rate when interactions are not present
eps.0 = IGR.plant(temp.A,s,g,lambda,res.eq)

#initialize dataframe to contain change in invasion growth rate
eps.df = as.data.frame(eps.0); names(eps.df)='eps.0'

# impact of competition from H. murinum on P. erecta
temp.A[1,3] <- A[1,3]
eps.13 = IGR.plant(temp.A, s, g, lambda, res.eq) - eps.0

# add to dataframe
eps.df = cbind(eps.df, eps.13)

# check it worked. Should read TRUE.
all(abs(IGR.plant(A,s,g,lambda,res.eq)-rowSums(eps.df))<1e-15)

# find delta values
delta.interest.values = eps.df[inv,]-eps.df[res,]
colnames(delta.interest.values) = gsub('^eps','delta',names(delta.interest.values))

# Create dataframe for plotting
df <- data.frame(
  Mechanism = c("Baseline and other weak contributions", 
                'Competition from P. erecta on A. wrangelianus'),
  IGR = as.numeric(delta.interest.values)
)

# fix order
df$Mechanism <- factor(df$Mechanism, levels = rev(c("Baseline and other weak contributions", 
                                                    'Competition from P. erecta on A. wrangelianus')))

# Create cumulative sum 
df$cumulative_IGR_backwards <- cumsum(df$IGR)  # Direct cumulative sum calculation

# specify resident community and invading species
res = 2
inv = 1

# get equilibrium density for original function given resident community
res.eq = numeric(nrow(A))
res.eq[res] = solve(a = A[res,res],b= lambda[res] * g[res] / (1 - s[res]*(1-g[res]))-1)/g[res]

# IGR for annual plant model
IGR.plant <- function(A,s,g,lambda,res.eq){ # function to get IGR for gLV
  
  IGR.row = log(s*(1-g) + g*lambda /(1 + ((A %*% (g*res.eq)))))
  
  return(IGR.row)
  
}

# get Invasion Growth Rate of specific model
orig.IGR = IGR.plant(A,s,g,lambda,res.eq)[inv]

# set interactions of interest to 0. In our case, competition between species.
temp.A=A
temp.A[1,2]<-0

# find baseline value i.e., invasion growth rate when interactions are not present
eps.0 = IGR.plant(temp.A,s,g,lambda,res.eq)

#initialize dataframe to contain change in invasion growth rate
eps.df = as.data.frame(eps.0); names(eps.df)='eps.0'

# impact of competition from H. murinum on P. erecta
temp.A[1,2] <- A[1,2]
eps.12 = IGR.plant(temp.A, s, g, lambda, res.eq) - eps.0

# add to dataframe
eps.df = cbind(eps.df, eps.12)

# check it worked. Should read TRUE.
all(abs(IGR.plant(A,s,g,lambda,res.eq)-rowSums(eps.df))<1e-15)

# find delta values
delta.values = eps.df[inv,]-colSums(eps.df[res,])/(length(res))
colnames(delta.values) = gsub('^eps','delta',names(delta.values))

# Only consider contributions whose effect is at least 0.01 in magnitude
indices = which(abs(delta.values)>0.01)
delta.interest.values = delta.values[indices]

# add all weak contributions to baseline
delta.interest.values[1] = delta.interest.values[1] + sum(delta.values[abs(delta.values)<0.01])

# Create dataframe for plotting
inv.df <- data.frame(
  Mechanism = c("Baseline and other weak contributions", 
                'Competition from H. murinum on A. wrangelianus'), 
  IGR = as.numeric(delta.interest.values)
)

# fix order
inv.df$Mechanism <- factor(inv.df$Mechanism, levels = rev(c("Baseline and other weak contributions", 
                                                            'Competition from H. murinum on A. wrangelianus')))

# Create cumulative sum 
inv.df$cumulative_IGR_backwards <- cumsum(inv.df$IGR)  # Direct cumulative sum calculation

# prepare dataframe for plotting
total.df = as.data.frame(inv.df$Mechanism)
colnames(total.df) = 'Mechanism'
total.df$minusi.comm.IGR = inv.df$IGR

trueindices = which(matchindices)
igr_indices = seq_along(df$IGR)

alter[trueindices] = df$IGR[seq_along(trueindices)]
total.df$secext.comm.IGR = alter

secextmech <- data.frame(
  Mechanism = "Competition from P. erecta on A. wrangelianus",
  minusi.comm.IGR = 0,
  secext.comm.IGR = -6.31604
)
total.df = rbind(total.df, secextmech)

require(ggplot2)
require(dplyr)
require(tidyr)
require(stringr)

# Reshape the dataframe to a long format
total_long <- total.df %>%
  pivot_longer(cols = c(minusi.comm.IGR, secext.comm.IGR), 
               names_to = "IGR_Type", values_to = "IGR")

# Compute Cumulative IGR (sum for each IGR_Type)
cumulative_igr <- total_long %>%
  group_by(IGR_Type) %>%
  summarize(IGR = sum(IGR), .groups = "drop") %>%
  mutate(Mechanism = "Cumulative IGR")  # Add "Cumulative IGR" as a new row

# Append Cumulative IGR row to the original dataframe
total_long <- bind_rows(total_long, cumulative_igr)

# Ensure "Cumulative IGR" is at the top and "Baseline" is at the bottom
mechanism_levels <- c("Cumulative IGR", 'Interaction between competition from A. wrangelianus and H. murinum on P. erecta', setdiff(unique(total.df$Mechanism), c("Cumulative IGR", "Baseline and other weak contributions",'Interaction between competition from A. wrangelianus and H. murinum on P. erecta')), "Baseline and other weak contributions")
total_long$Mechanism <- factor(total_long$Mechanism, levels = mechanism_levels, ordered = TRUE)

# Get the position (index) of "Cumulative IGR" for the vertical line
cumulative_igr_position <- which(levels(total_long$Mechanism) == "Cumulative IGR")

# Plot
ggplot(total_long, aes(x = Mechanism, y = IGR, fill = IGR_Type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black") +
  coord_flip() +  
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.5) +
  geom_vline(xintercept = cumulative_igr_position + 0.5, linetype = "dotted", color = "black", size = 0.7) +
  theme_classic() +  
  labs(x = "Mechanism", y = "Invasion Growth Rate", fill = "Community") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 27)) + 
  scale_fill_manual(values = c("minusi.comm.IGR" = "#fec488", "secext.comm.IGR" = "#e44f64"),  
                    labels = c("minusi.comm.IGR" = "-P. erecta", "secext.comm.IGR" = "-A. wrangelianus")) +
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18, face = "bold"),   
    plot.title = element_text(size = 18, face = "bold"), 
    legend.position = "top",
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18)
  )
