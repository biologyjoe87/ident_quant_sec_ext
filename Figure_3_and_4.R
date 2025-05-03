### secondary extinction for perennial plant system lotka volterra
### Author: Joe Brennan

library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)

# Figure 3a. Disassembly graph for the grassland community.

spring.b = c(100.77, 45.19, 103.43, 67.09, 112.71, 76.52)/1000
spring.A = -rbind(c(1.28, -0.15, -0.10, -0.22, 0.16, -0.12),
                  c(-0.13, 0.41, 0.12, 0.09, -0.03, 0),
                  c(0.38, 0.30, 0.41, 0.36, 0.15, 0.11),
                  c(-0.95, 0.07, 0.36, 0.31, -3.6, 0.34),
                  c(0.39, 0.15, 0.82, 0.72, 0.8, -0.78),
                  c(0.02, 0.18, 0.05, 0.13, 0.57, 0.56))

normalized.spring.A = spring.A

for (i in 1:nrow(spring.A)){
  
  for (j in 1:ncol(spring.A)){
    
    normalized.spring.A[i,j] = (-spring.A[i,j] * spring.b[j] ) / ( spring.A[j,j] * spring.b[i] )
    
  }
  
}

interest.spp = c(1,2,3,5,6)

# A = normalized.spring.A[interest.spp, interest.spp]
# b = rep(1,5)

A = spring.A[interest.spp, interest.spp]
b = spring.b[interest.spp]

# Define species names
full.spp = c('Agrostis stolonifera', 'Lolium perenne',  'Phleum pratense', 'Trifolium pratense', 'Trifolium repens')
abbr.spp = c('A. stolonifera', 'L. perenne',  'P. pratense', 'T. pratense', 'T. repens')

# Create a matrix containing all invasion growth rates at each coexisting state. Rows represent coexisting communties while columns are species in the model. Thus, the i-th row and j-th column contains species j's invasion growth rate at the i-th community.
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

# Create the disassembly graph for this model. 
Disassembly.graph = function(IS){ 
  
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

# Get disassembly graph for this specific model
out = Disassembly.graph(IS)

# Plotting of the disassembly graph. Code modified from Hofbauer & Schreiber, 2022.
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
plot.disassembly.graph(out, cols=cols, vlx=0.82)
species_numbers <- seq_along(abbr.spp)
legend_labels <- paste(species_numbers, abbr.spp, sep=": ")
legend(x=0.6, y=1.2, legend = legend_labels, col="black", text.col="black", bty="o", cex=1.2, title = 'Species')

####################################################################################################################################

# Figure 2b. Simulate the dynamics depicting the secondary extinction

# Plotter for LV systems that allow for extinction events
sequential_extinction_GLV_plotter=function(A,b,tf=100,dt=0.1,loss.seq=c(),threshold.density=0.01,theta=1,graph=FALSE, comm.start = c())
{
  lv=function(t,x,parms){
    with(parms,{
      dx<-x*(b+A%*%x)
      list(dx)
    })
  }
  k=length(loss.seq)
  times=seq(0,tf,by=dt)
  parms=list(b=b,A=A)
  sol = solve(A[comm.start, comm.start], -b[comm.start]) # Start at equilibria
  sol[sol<0]=0
  initialx=rep(0,length(b))
  initialx[comm.start] = sol
  out=ode(y=initialx,times=times,func=lv,parms=parms)
  x=out[,2:(1+length(b))]
  for(i in 1:length(loss.seq)){
    initialx=out[length(times),2:(1+length(b))]
    initialx[which(initialx<threshold.density)]=0
    initialx[loss.seq[i]]=0
    out=ode(y=initialx,times=times,func=lv,parms=parms)
    x=rbind(x,out[,2:(1+length(b))])
  }
  #par(mar=c(4.5,4.5,0,1))
  ts=dt*(1:length(x[,1]))
  if(graph)matplot(ts,x[,comm.start]^theta,type="l",lty=1,lwd=4,bty="n",xlab="time t",ylab=expression(paste("densities ",x[,comm.start])));legend("topright", colnames(x[,comm.start]),col=seq_len(ncol(x[,comm.start])),cex=0.8,fill=seq_len(ncol(x[,comm.start])))
  return(list(x=x,ts=ts))
}

# simulate removal
ext.out = sequential_extinction_GLV_plotter(A=A, b=b, tf=1000, dt=0.1, loss.seq = c(5), comm.start=1:5)

# Create dataframe for plotting dynamics
df <- as.data.frame(ext.out$x) %>%
  mutate(time = ext.out$ts) %>%
  pivot_longer(cols = -time, names_to = "Species", values_to = "density")

species_colors <- c('darkorange', 'darkgreen', 'purple', 'darkblue', 'cornflowerblue')

ggplot(df, aes(x = time, y = density, color = Species)) +
  scale_color_manual(
    values = species_colors,
    labels = abbr.spp
  ) +
  geom_line(size = 1.2) +
  ylim(c(0,0.2)) +
  theme_classic() +
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

###################################################################################################################

# Figure 3c. Show -T. pratense community

# Simulate dynamcis where P. erecta goes extinct
ext.out = sequential_extinction_GLV_plotter(A=A, b=b, tf=1000, dt=0.1, loss.seq = c(4), comm.start=1:5)

# Plotting
df <- as.data.frame(ext.out$x) %>%
  mutate(time = ext.out$ts) %>%
  pivot_longer(cols = -time, names_to = "Species", values_to = "density")

species_colors <- c('darkorange', 'darkgreen', 'purple', 'darkblue', 'cornflowerblue')

ggplot(df, aes(x = time, y = density, color = Species)) +
  scale_color_manual(
    values = species_colors,
    labels = abbr.spp
  ) +
  geom_line(size = 1.2) +
  ylim(c(0,0.2)) + 
  theme_classic() +
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

####################################################################################################################################################

### Figure 4
### Author: Joe Brennan

# Figure 2a. Decomposition of T. pratense

library(ggplot2)
library(stringr)

# Invasion growth rate calculations
IGR <- function(A,b,res.eq){
  
  IGR.row = A%*%res.eq + b
  
  return(IGR.row)
  
}

# Decomposition function for LV systems. 
decomposition.particular.spp <- function(A, b, res, inv, avg = F, spp.compare){
  
  # get equilibrium density for original function given resident community
  res.eq = numeric(nrow(A))
  res.eq[res] = solve(A[res,res], -b[res])
  
  # get Invasion Growth Rate
  orig.IGR = IGR(A,b,res.eq)[inv]
  
  # Choose indices. Do not choose ones on diagonal (conspecific interactions) and terms associated with non-resident and non-invader species
  opt = sort(union(inv,res))
  indices = which(A!=diag(A) & !(col(A)==inv) & A %in% A[opt,opt] & A!=0)
  index.row.col = arrayInd(indices, dim(A))
  
  # Calculate average values across all values of interest
  if (avg==T) avg.val = mean(A[indices])
  if (avg==F) avg.val = 0
  
  # create temporary matrix, temp.A, to hold the changed version of A
  temp.A = A
  temp.A[indices] = avg.val
  
  # calculate epsilon naught and respective delta
  eps.0 = IGR(temp.A,b,res.eq)
  delta.e.star = eps.0[inv] - sum(eps.0[spp.compare])/(length(spp.compare))
  
  # Create lists we will store values in 
  eps_names = c('eps.0')
  eps_values = c(eps.0[inv])
  delta_names_eq = c('delta.0')
  delta_values_eq = c(delta.e.star)
  
  # loop through all other possibilities
  for (i in 1:length(indices)){
    
    # identify current term of interest
    current.index = indices[i]
    current.row.col = index.row.col[i,]
    
    # restore original values of respective term
    temp.A[current.index] = A[current.index]
    
    # store temp epsilon and temp delta
    temp.epsilon = IGR(temp.A, b, res.eq) - eps.0
    temp.delta.eq = temp.epsilon[inv] - sum(temp.epsilon[spp.compare])/(length(spp.compare))
    
    # create variable name
    eps.var.name = paste0('eps.',current.row.col[1],current.row.col[2])
    delta.var.name = paste0('delta.',current.row.col[1],current.row.col[2])
    
    # assign outputs
    assign(eps.var.name, temp.epsilon[inv])
    assign(delta.var.name, temp.delta.eq)
    
    # store variable names and values
    eps_names = c(eps_names, eps.var.name)
    eps_values = c(eps_values, temp.epsilon[inv])
    delta_names_eq = c(delta_names_eq, delta.var.name)
    delta_values_eq = c(delta_values_eq, temp.delta.eq)
    
    # turn value back
    if (avg==T) temp.A[current.index] = avg.val
    if (avg==F) temp.A[current.index] = 0
    
  }
  
  if (abs(orig.IGR - sum(delta_values_eq))<1e-8) print('eq woohoo') else print('eq oh no')
  
  # Create final arrays :D
  final.delta.array=setNames(delta_values_eq, delta_names_eq)
  final.eps.array = setNames(eps_values, eps_names)
  
  # make dataframes 
  delta_df_eq = data.frame(Delta = names(final.delta.array), Delta.Value = as.numeric(final.delta.array))
  eps_df = data.frame(Epsilon = names(final.eps.array), Epsilon.Value = as.numeric(final.eps.array))
  
  df = cbind(eps_df, delta_df_eq)
  
  return(list(IGR = orig.IGR,decomposition=df))
}

# Do decomposition of P. lanceolata at -A. stolonifera community
decomp = decomposition.particular.spp(A,b,res=c(1,2,3), inv=4, avg=F, spp.compare = c(1,2,3)) # only compare IGR to resident species who have competitive effect on species 4

# Extract non-zero contributions
delta.interest = (decomp$decomposition$Delta.Value[decomp$decomposition$Delta.Value!=0])

# Create dataframe for plotting
df <- data.frame(
  Mechanism = (c('Baseline',
                 'Facilitation from A. stolonifera on L. perenne',
                 'Competition from A. stolonifera on P. pratense',
                 'Competition from A. stolonifera on T. pratense',
                 'Facilitation from L. perenne on A. stolonifera',
                 'Competition from L. perenne on P. pratense',
                 'Competition from L. perenne on T. pratense',
                 'Facilitation from P. pratense on A. stolonifera',
                 'Competition from P. pratense on L. perenne',
                 'Competition from P. pratense on T. pratense')),
  IGR = delta.interest
)

# Fix the ordering
df$Mechanism <- factor(df$Mechanism, levels = c('Baseline',
                                                'Facilitation from A. stolonifera on L. perenne',
                                                'Competition from A. stolonifera on P. pratense',
                                                'Competition from A. stolonifera on T. pratense',
                                                'Facilitation from L. perenne on A. stolonifera',
                                                'Competition from L. perenne on P. pratense',
                                                'Competition from L. perenne on T. pratense',
                                                'Facilitation from P. pratense on A. stolonifera',
                                                'Competition from P. pratense on L. perenne',
                                                'Competition from P. pratense on T. pratense'))

#########################################################################################################

# Perform decomposition
inv.decomp = decomposition.particular.spp(A,b,res=c(1,2,3,5),inv=4, avg=F, spp.compare= c(1,2,3))

# Extract non-zero values
delta.interest = (inv.decomp$decomposition$Delta.Value[inv.decomp$decomposition$Delta.Value!=0])

# Create dataframe for plotting
inv.df <- data.frame(
  Mechanism = (c('Baseline',
                 'Facilitation from A. stolonifera on L. perenne',
                 'Competition from A. stolonifera on P. pratense',
                 'Competition from A. stolonifera on T. pratense',
                 'Facilitation from L. perenne on A. stolonifera',
                 'Competition from L. perenne on P. pratense',
                 'Competition from L. perenne on T. pratense',
                 'Facilitation from P. pratense on A. stolonifera',
                 'Competition from P. pratense on L. perenne',
                 'Competition from P. pratense on T. pratense',
                 'Facilitation from T. repens on A. stolonifera',
                 'Competition from T. repens on P. pratense',
                 'Facilitation from T. repens on T. pratense')),
  IGR = delta.interest
)

# Fix order of contributions
inv.df$Mechanism <- factor(inv.df$Mechanism, levels = c('Baseline',
                                                        'Facilitation from A. stolonifera on L. perenne',
                                                        'Competition from A. stolonifera on P. pratense',
                                                        'Competition from A. stolonifera on T. pratense',
                                                        'Facilitation from L. perenne on A. stolonifera',
                                                        'Competition from L. perenne on P. pratense',
                                                        'Competition from L. perenne on T. pratense',
                                                        'Facilitation from P. pratense on A. stolonifera',
                                                        'Competition from P. pratense on L. perenne',
                                                        'Competition from P. pratense on T. pratense',
                                                        'Facilitation from T. repens on A. stolonifera',
                                                        'Competition from T. repens on P. pratense',
                                                        'Facilitation from T. repens on T. pratense'))


##################

total.df = as.data.frame(inv.df$Mechanism)
colnames(total.df) = 'Mechanism'
total.df$minusi.comm.IGR = inv.df$IGR

matchindices = inv.df$Mechanism %in% df$Mechanism
alter = numeric(length(matchindices))

trueindices = which(matchindices)
igr_indices = seq_along(df$IGR)

alter[trueindices] = df$IGR[seq_along(trueindices)]
total.df$secext.comm.IGR = alter

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

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
mechanism_levels <- c("Cumulative IGR", setdiff(unique(total.df$Mechanism), c("Cumulative IGR", "Baseline")), "Baseline")
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
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30)) + 
  scale_fill_manual(values = c("minusi.comm.IGR" = "#fec488", "secext.comm.IGR" = "#e44f64"),  
                    labels = c("minusi.comm.IGR" = "-T. pratense", "secext.comm.IGR" = "-T. repens")) +
  
  # Add y-axis break from 3 to 10
  #scale_y_break(c(3, 10), scales = "fixed") +  
  
  # Manually adjust tick marks before and after the break
  # scale_y_continuous(
  #   breaks = c(seq(-2, 3, by = 1), seq(10, max(total_long$IGR), by = 1))
  # ) +
  
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18, face = "bold"),   
    plot.title = element_text(size = 18, face = "bold"), 
    legend.position = "top",
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    
    # Remove the right y-axis (duplicate y-axis labels)
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank()
  )