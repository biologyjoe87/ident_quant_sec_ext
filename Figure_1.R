### Figure 1
### Author: Joe Brennan

library('tidyverse')
library(igraph)

# First bit of code from Van Dyke et al, 2022 code

# # read in data
# s_g_data <- read.csv("./data/s_g_data.csv")
# 
# # Extract parameters
# s = s_g_data$s
# g = s_g_data$g
# 
# final_output <- read.csv("./data/final_output_nls_boot_1000.csv")
# final_output$treatment <- factor(final_output$treatment, levels = c(1, 2))
# 
# # Calculate the mean lambda values for each species at different treatments
# lambda_summary <- final_output %>%
#   group_by(focal, treatment) %>%
#   summarize(mean_lambda = mean(lambda, na.rm = TRUE))
# 
# lambda.ambient = lambda_summary$mean_lambda[lambda_summary$treatment==1]
# 
# nls_boot_pairs <- read.csv("./output/nls_boot_pairs_1000_full_model.csv")
# 
# A.ambient = matrix(nls_boot_pairs$alpha[nls_boot_pairs$treatment==1], nrow=6,ncol=6)
# A.ambient.names = matrix(paste0(nls_boot_pairs$focal[nls_boot_pairs$treatment==1],'.',nls_boot_pairs$competitor[nls_boot_pairs$treatment==1]), nrow=6, ncol=6)

# Define colors for disassembly graph
cols = c('white', '#FFDB00', 'white', 'white', 'white', 'black', '#FFDB00', 'black')

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

# Focus on community of interest
# A = A.ambient[c(1,3,4), c(1,3,4)]
# g = g[c(1,3,4)]
# s = s[c(1,3,4)]
# lambda = lambda.ambient[c(1,3,4)]

A = rbind(c(0.27290462, 0.9415389, 0.5590646),
          c(0.13755487, 0.9483990,0.1599514),
          c(0.08483079, 1.0221670, 0.4161495))
g = c(0.2240310, 0.6666667, 0.6029647)
s = c(0.2625463, 0.0450000, 0.6215016)
lambda = c(2351.5366, 745.9973, 654.6758)

# Species names
abbr.spp = c('A. wrangelianus', 'H. murinum', 'P. erecta')

# Get the invasion scheme - matrix of invasion growth rates
IS = plant.IS(A,g,s,lambda)

# Create community disassembly graph
alternate.disassembly.graph = function(IS){ 
  
  n=dim(IS)[1] # number of communities
  composition=list()
  edge.label = c() # list of edge labels detailing which invader leads to transition
  for(i in 1:n)composition[[i]]=which(IS[i,]==0)

  
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
  index = which(edge.label!=0)
  edge.label = edge.label[index]
  
  return(list(adj.mat=mat,edge.label = edge.label, IS = IS, composition=composition))
  
}

# Create disassembly graph for the specific model
out=alternate.disassembly.graph(IS)

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
  xvals[5] = xvals[5]-0.5
  
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

# plot specific disassemly graph
plot.disassembly.graph(out,cols=cols)

# Create legend
species_numbers <- seq_along(abbr.spp)
legend_labels <- paste(species_numbers, abbr.spp, sep=": ")
legend(x=0.45, y=1.2, legend = legend_labels, col="black", text.col="black", bty="o", cex=1.2, title = 'Species')

############################################################################################################################################

# Figure 3B. Simulation showing secondary extinction

# Secondary Extinction of P. erecta when A. wrangelanius is lost

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
  if (i==75) mat[i,1]=0 # extinction event
}

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

###############################################################################################################

# Fig 3c. Show A. wrangelanius and H. murinum is -P. erecta community

# set intiial conditions
init = c(3000, 500, 190)

# set number of iterations
iterations = 100

# initialize matrix to hold time series
mat = matrix(NA, nrow=iterations, ncol=length(init))
mat[1,]=init

# simulate the dynamics
for (i in 2:nrow(mat)){
  mat[i,] = s*(1-g)*mat[i-1,] + lambda*g*mat[i-1,] / (1 + A %*% (g*mat[i-1,]))
  if (i==50) mat[i,3]=0
}

mat = cbind(mat, 1:iterations)
colnames(mat) = c("A. wrangelianus", "H. murinum", "P. erecta", 'Iterations')

# create dataframe for plotting
df <- as.data.frame(mat) %>%
  pivot_longer(cols = -Iterations, names_to = "Species", values_to = "Abundance") 

# plot
ggplot(df, aes(x = Iterations, y = Abundance, color = Species)) +
  scale_color_manual(
    values = c("A. wrangelianus" = "orange", "H. murinum" = "blue", "P. erecta" = "black")) +
  geom_line(size = 1.2) +
  theme_classic() +
  labs(x = "Time",
       y = "Density") +
  theme(legend.position = c(0.8, 0.75),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.background = element_rect(color = "black", fill = NA),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 14),  
        axis.text.y = element_text(size = 14))
