### Figure 2. Decomposition of invasion growth rates. 
### Author: Joe Brennan

# First chunk of code is modified from Van Dyke et al, 2022

# Read in parameter data
# s_g_data <- read.csv("./data/s_g_data.csv") #seed survival and germination data
# 
# # extract parameters
# s = s_g_data$s
# g = s_g_data$g
# 
# final_output <- read.csv("./output/final_output_nls_boot_1000.csv")
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
# 
# # Get parameters for species of interest
# A = A.ambient[c(1,3,4), c(1,3,4)]
# g = g[c(1,3,4)]
# s = s[c(1,3,4)]
# lambda = lambda.ambient[c(1,3,4)]

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
res = 2
inv = 3

# get equilibrium density for original function given resident community
res.eq = numeric(nrow(A))
res.eq[res] = solve(a = A[res,res],b= lambda[res] * g[res] / (1 - s[res]*(1-g[res]))-1)/g[res]

# get Invasion Growth Rate of specific model
orig.IGR = IGR.plant(A,s,g,lambda,res.eq)[inv]

# set interactions of interest to 0. In our case, competition between species.
temp.A=A
temp.A[3,2]<-0

# find baseline value i.e., invasion growth rate when interactions are not present
eps.0 = IGR.plant(temp.A,s,g,lambda,res.eq)

#initialize dataframe to contain change in invasion growth rate
eps.df = as.data.frame(eps.0); names(eps.df)='eps.0'

# impact of competition from H. murinum on P. erecta
temp.A[3,2] <- A[3,2]
eps.32 = IGR.plant(temp.A, s, g, lambda, res.eq) - eps.0

# add to dataframe
eps.df = cbind(eps.df, eps.32)

# check it worked. Should read TRUE.
all(abs(IGR.plant(A,s,g,lambda,res.eq)-rowSums(eps.df))<1e-15)

# find delta values
delta.interest.values = eps.df[inv,]-eps.df[res,]
colnames(delta.interest.values) = gsub('^eps','delta',names(delta.interest.values))

# Create dataframe for plotting
df <- data.frame(
  Mechanism = c("Baseline and other weak contributions", 
                'Competition from H. murinum on P. erecta'),
  IGR = as.numeric(delta.interest.values)
)

# fix order
df$Mechanism <- factor(df$Mechanism, levels = rev(c("Baseline and other weak contributions", 
                                                    'Competition from H. murinum on P. erecta')))

# Create cumulative sum 
df$cumulative_IGR_backwards <- cumsum(df$IGR)  # Direct cumulative sum calculation

# specify resident community and invading species
res = c(1,2)
inv = 3

# get equilibrium density for original function given resident community
res.eq = numeric(nrow(A))
res.eq[res] = solve(a = A[res,res],b= lambda[res] * g[res] / (1 - s[res]*(1-g[res]))-1)/g[res]

# IGR for annual plant model
IGR.plant <- function(A,s,g,lambda,res.eq){ # function to get IGR for gLV
  
  IGR.row = log(s*(1-g) + g*lambda /(1 + ((A %*% (g*res.eq)))))
  
  return(IGR.row)
  
}

# function to find all strict subsets of an array
all_subsets <- function(x) {
  
  
  # Create a list to store the subsets
  subsets <- list()
  
  # Iterate over all possible subset sizes
  for (i in 1:(length(x)-1)) {
    subsets <- c(subsets, combn(x, i, simplify = FALSE))
  }
  
  if (length(x)==1) subsets = list(x)
  
  return(subsets)
}

# get Invasion Growth Rate of specific model
orig.IGR = IGR.plant(A,s,g,lambda,res.eq)[inv]

# set interactions of interest to 0. In our case, competition between species.
temp.A=A
temp.A[1,2]<-temp.A[2,1]<-temp.A[3,1]<-temp.A[3,2]<-0

# find baseline value i.e., invasion growth rate when interactions are not present
eps.0 = IGR.plant(temp.A,s,g,lambda,res.eq)

# find all species interactions combinations
subsets = all_subsets(1:4)

#initialize dataframe to contain change in invasion growth rate
eps.df = as.data.frame(eps.0); names(eps.df)='eps.0'

# loop through all subsets
for (i in 1:length(subsets)){
  
  # identify traits we are restoring
  current.subset = subsets[[i]]
  
  # find eps to subtract
  int.effect.subtraction = c(1)
  if(i-1>0)for (j in 1:(i-1)){
    if (all(subsets[[j]] %in% current.subset)) int.effect.subtraction = c(int.effect.subtraction, (j+1))
  }
  
  #restore respective traits
  if (1 %in% current.subset) temp.A[1,2]=A[1,2]
  if (2 %in% current.subset) temp.A[2,1]=A[2,1]
  if (3 %in% current.subset) temp.A[3,2]=A[3,2]
  if (4 %in% current.subset) temp.A[3,1]=A[3,1]
  
  # find invasion growth rate of current system
  r = IGR.plant(temp.A,s,g,lambda,res.eq)
  
  # subtract effects to extract main effect
  temp.eps = r - rowSums(cbind(eps.df[,int.effect.subtraction], rep(0,3)))
  
  # add to df
  eps.df[,(i+1)] = temp.eps; names(eps.df)[i+1]=paste0('eps.', paste0(current.subset, collapse='.'))
  
  # set values back to 0
  temp.A[1,2]<-temp.A[2,1]<-temp.A[3,1]<-temp.A[3,2]<-0
  
}

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
                'Competition from H. murinum on A. wrangelianus', 
                "Competition from A. wrangelianus on H. murinum", 
                "Competition from H. murinum on P. erecta",
                'Competition from A. wrangelianus on P. erecta',
                'Interaction between competition from A. wrangelianus and H. murinum on P. erecta'),
  IGR = as.numeric(delta.interest.values)
)

# fix order
inv.df$Mechanism <- factor(inv.df$Mechanism, levels = rev(c("Baseline and other weak contributions", 
                                                    'Competition from H. murinum on A. wrangelianus', 
                                                    "Competition from A. wrangelianus on H. murinum", 
                                                    "Competition from H. murinum on P. erecta",
                                                    'Competition from A. wrangelianus on P. erecta',
                                                    'Interaction between competition from A. wrangelianus and H. murinum on P. erecta')))

# Create cumulative sum 
inv.df$cumulative_IGR_backwards <- cumsum(inv.df$IGR)  # Direct cumulative sum calculation

# prepare dataframe for plotting
total.df = as.data.frame(inv.df$Mechanism)
colnames(total.df) = 'Mechanism'
total.df$minusi.comm.IGR = inv.df$IGR

# find shared mechanisms between communities
matchindices = inv.df$Mechanism %in% df$Mechanism
alter = numeric(length(matchindices))

trueindices = which(matchindices)
igr_indices = seq_along(df$IGR)

alter[trueindices] = df$IGR[seq_along(trueindices)]
total.df$secext.comm.IGR = alter

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
