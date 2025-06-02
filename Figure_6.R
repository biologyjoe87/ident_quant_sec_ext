### Figure 6. Decomposition of invasion growth rates
### Author: Joe Brennan

# Code adapted from Shoemaker et al, 2020

require(deSolve)

# A. decomposition of C2 in -C2 comm

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

# timesteps
time <- seq(0,100000,by=0.1) 

# Make P extinct
State <- c(P = 1, C1 = 1, C2 = 0, R = 1)

# create array of parameters
pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R)

# simulate dynamics
out <- ode(y = State, times = time, func = VassFox_Cvar, parms = pars)

# length of burn in period
invade_start_time <- 10000

# initialize matrix to hold invasion growth rates
C2_ldgr <- matrix(data=NA, nrow=(length(time)-invade_start_time), ncol=1)
C1_resident <- C2_ldgr
counter <- 1

# loop through each time step after burn in period and calculate invasion growth rate
for (t in invade_start_time:length(time)){
  
  P = out[t,2]
  C1 = out[t,3]
  C2 = out[t,4]
  
  R = out[t,5]
  
  C2_ldgr[counter] <- - (M_C2) + ( (O_C2_R * J_C2 * R) / (R + R_0_2) ) - ( ((1-O_P_C1) * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
  C1_resident[counter] <- - (M_C1) + ( (O_C1_R * J_C1 * R) / (R + R_0_1) ) - ( (O_P_C1 * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
  
  counter <- counter + 1
  
}

# calculate invasion growth rate of competitors. Resident competitor should be roughly equal to 0.
C2_r_bar <- mean(C2_ldgr)-mean(C1_resident)

# ----------------------------------------------------------------------------------------------------
# calculate delta_0. IGR under no resource fluctuation, no predator preference, and no C1 intake of resource

# no preference between prey items
O_P_C1 <- 0.5

# average resource abundance
avg_R = mean(out[,5])
R = avg_R

# no resource intake from resident competitor
J_C1 = 0

# intialize matrix to store baseline values
C2_epsilon_0 <- matrix(data=NA, nrow=(length(time)-invade_start_time), ncol=1)
C1_resident_epsilon_0 <- C2_epsilon_0
counter <- 1

# calculate instantaneous growth rate of both competitors after the burn in period.
for (t in invade_start_time:length(time)){
  
  P = out[t,2]
  C1 = out[t,3]
  C2 = out[t,4]
  
  C2_epsilon_0[counter] <- - (M_C2) + ( (O_C2_R * J_C2 * R) / (R + R_0_2) ) - ( ((1-O_P_C1) * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
  C1_resident_epsilon_0[counter] <- - (M_C1) + ( (O_C1_R * J_C1 * R) / (R + R_0_1) ) - ( (O_P_C1 * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
  
  counter <- counter + 1
  
}

# calculate mean impact of contribution to invading competitor relative the resident
C2_delta_0 <- mean(C2_epsilon_0)-mean(C1_resident_epsilon_0)

# ----------------------------------------------------------------------------------------------------
# calculate delta_R, with R fluctuations but no predator pref and no intake of resource from C1

# initialize matrices to hold values
C2_epsilon_R <- matrix(data=NA, nrow=(length(time)-invade_start_time), ncol=1)
C1_resident_epsilon_R <- C2_epsilon_R
counter <- 1

# no preference between prey items
O_P_C1 <- 0.5

# No resource uptake from resident competitor
J_C1 = 0

# loop through each time step after the burn in period and calculate the instantaneous growth rate of the invading and resident competitors. 
for (t in invade_start_time:length(time)){
  
  # allow R values to change - do not keep average
  R = out[t,5]
  
  P = out[t,2]
  C1 = out[t,3]
  C2 = out[t,4]
  
  C2_epsilon_R[counter] <- - (M_C2) + ( (O_C2_R * J_C2 * R) / (R + R_0_2) ) - ( ((1-O_P_C1) * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
  C1_resident_epsilon_R[counter] <- - (M_C1) + ( (O_C1_R * J_C1 * R) / (R + R_0_1) ) - ( (O_P_C1 * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
  
  counter <- counter + 1
  
}

# calculate mean impact of resource fluctuations on invading competitor relative to the resident competitor
C2_delta_R <- mean(C2_epsilon_R)-mean(C1_resident_epsilon_R) - C2_delta_0

# ----------------------------------------------------------------------------------------------------
# calculate delta_Omega
# constant resource and no intake of resource from C1 but predator preference

# preference between prey items
O_P_C1 <- 0.92

# no resource uptake from resident competitor
J_C1 = 0

# Keep R values constant
R = avg_R

# initialize matrices to store values
C2_epsilon_omega <- matrix(data=NA, nrow=(length(time)-invade_start_time), ncol=1)
C1_resident_epsilon_omega <- C2_epsilon_omega
counter <- 1

# loop through each time step after the burn in period and calculate the instantaneous growth rate of the invading and resident competitors. 
for (t in invade_start_time:length(time)){
  
  P = out[t,2]
  C1 = out[t,3]
  C2 = out[t,4]
  
  C2_epsilon_omega[counter] <- - (M_C2) + ( (O_C2_R * J_C2 * R) / (R + R_0_2) ) - ( ((1-O_P_C1) * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
  C1_resident_epsilon_omega[counter] <- - (M_C1) + ( (O_C1_R * J_C1 * R) / (R + R_0_1) ) - ( (O_P_C1 * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
  
  counter <- counter + 1
  
}

# calculate mean impact of predator preference on invading competitor relative to the resident competitor
C2_delta_omega <- mean(C2_epsilon_omega)-mean(C1_resident_epsilon_omega) - C2_delta_0

# ----------------------------------------------------------------------------------------------------
# calculate delta_JC1
# constant resource and no pred pref but resource intake of C1

# no preference between prey items
O_P_C1 <- 0.5 

# Resident competitor uptake rate of R
J_C1 = 0.8036

# Keep R constant
R = avg_R

# initialize matrices to hold values
C2_epsilon_JC1 <- matrix(data=NA, nrow=(length(time)-invade_start_time), ncol=1)
C1_resident_epsilon_JC1 <- C2_epsilon_JC1
counter <- 1

# evaluate instantaneous growth rate at each time step after burn in period for both resident and invading competitor
for (t in invade_start_time:length(time)){
  
  P = out[t,2]
  C1 = out[t,3]
  C2 = out[t,4]
  
  C2_epsilon_JC1[counter] <- - (M_C2) + ( (O_C2_R * J_C2 * R) / (R + R_0_2) ) - ( ((1-O_P_C1) * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
  C1_resident_epsilon_JC1[counter] <- - (M_C1) + ( (O_C1_R * J_C1 * R) / (R + R_0_1) ) - ( (O_P_C1 * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
  
  counter <- counter + 1
  
}

# calculate the mean impact of resource intake from resident competitor on invading competitor relative to the resident competitor
C2_delta_JC1 <- mean(C2_epsilon_JC1)-mean(C1_resident_epsilon_JC1) - C2_delta_0

# ----------------------------------------------------------------------------------------------------
# calculate delta_R_omega
# fluctuating resource and pred pref but no resource intake of C1

# preference between prey items
O_P_C1 <- 0.92

# no resource uptake from resident competitor
J_C1 = 0

# initialize matrices to hold values
C2_epsilon_R_omega <- matrix(data=NA, nrow=(length(time)-invade_start_time), ncol=1)
C1_resident_epsilon_R_omega <- C2_epsilon_R_omega
counter <- 1

# calculate instantaneous growth rate of both competitors for each time step after the burn in period
for (t in invade_start_time:length(time)){
  
  P = out[t,2]
  C1 = out[t,3]
  C2 = out[t,4]
  
  R = out[t,5]
  
  C2_epsilon_R_omega[counter] <- - (M_C2) + ( (O_C2_R * J_C2 * R) / (R + R_0_2) ) - ( ((1-O_P_C1) * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
  C1_resident_epsilon_R_omega[counter] <- - (M_C1) + ( (O_C1_R * J_C1 * R) / (R + R_0_1) ) - ( (O_P_C1 * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
  
  counter <- counter + 1
  
}

# Calculate mean interactive effect between resource fluctuations and predator preference for invading species relative to the resident competitor
C2_delta_R_omega <- mean(C2_epsilon_R_omega)-mean(C1_resident_epsilon_R_omega) - (C2_delta_0 + C2_delta_R + C2_delta_omega) # should equal zero

# ----------------------------------------------------------------------------------------------------
# calculate delta_R_JC1
# fluctuating resource and resource intake by C1 but no pred pref

# No preference between prey items
O_P_C1 <- 0.5

# Resource uptake from resident competitor
J_C1 = 0.8036

# initialize matrices to store values
C2_epsilon_R_JC1 <- matrix(data=NA, nrow=(length(time)-invade_start_time), ncol=1)
C1_resident_epsilon_R_JC1 <- C2_epsilon_R_JC1
counter <- 1

# calculate instantaneous growth rate of both competitor species for each time step after the burn in period
for (t in invade_start_time:length(time)){
  
  P = out[t,2]
  C1 = out[t,3]
  C2 = out[t,4]
  
  R = out[t,5]
  
  C2_epsilon_R_JC1[counter] <- - (M_C2) + ( (O_C2_R * J_C2 * R) / (R + R_0_2) ) - ( ((1-O_P_C1) * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
  C1_resident_epsilon_R_JC1[counter] <- - (M_C1) + ( (O_C1_R * J_C1 * R) / (R + R_0_1) ) - ( (O_P_C1 * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
  
  counter <- counter + 1
  
}

# calculate mean impact of interaction between resource fluctuations and resource uptake by resident competitor on the invading competitor relative to the resident competitor
C2_delta_R_JC1 <- mean(C2_epsilon_R_JC1)-mean(C1_resident_epsilon_R_JC1) - (C2_delta_0 + C2_delta_R + C2_delta_JC1) 

# ----------------------------------------------------------------------------------------------------
# calculate delta_omega_JC1
# and resource intake by C1 and pred pref but no fluc in R

# preference between prey items
O_P_C1 <- 0.92

# Resource uptake by resident competitor
J_C1 = 0.8036

# Constant resource density
R = avg_R

# initialize matrices to hold values
C2_epsilon_omega_JC1 <- matrix(data=NA, nrow=(length(time)-invade_start_time), ncol=1)
C1_resident_epsilon_omega_JC1 <- C2_epsilon_omega_JC1
counter <- 1

# loop through each time step after burn in period and calculate instantaneous growth rate for both competitor species
for (t in invade_start_time:length(time)){
  
  P = out[t,2]
  C1 = out[t,3]
  C2 = out[t,4]
  
  C2_epsilon_omega_JC1[counter] <- - (M_C2) + ( (O_C2_R * J_C2 * R) / (R + R_0_2) ) - ( ((1-O_P_C1) * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
  C1_resident_epsilon_omega_JC1[counter] <- - (M_C1) + ( (O_C1_R * J_C1 * R) / (R + R_0_1) ) - ( (O_P_C1 * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
  
  counter <- counter + 1
  
}

# calculate mean impact of interactions between predator preference and resource uptake by  resident competitor to the invading competitor relative to the resident competitor.
C2_delta_omega_JC1 <- mean(C2_epsilon_omega_JC1)-mean(C1_resident_epsilon_omega_JC1) - (C2_delta_0 + C2_delta_omega + C2_delta_JC1)  # should be zero

# Calculate three way interaction
C2_delta_R_omega_JC1 = C2_r_bar-(C2_delta_0+C2_delta_R+C2_delta_omega+C2_delta_JC1 + C2_delta_R_omega + C2_delta_R_JC1 + C2_delta_omega_JC1) # should be zero

# ----------------------------------------------------------------------------------------------------

# Load ggplot2 library
require(ggplot2)

# only consider non-zero
C2_final_mechanisms <- c(C2_delta_R_JC1, C2_delta_JC1, C2_delta_omega, C2_delta_R, C2_delta_0)

# Create dataframe for plotting
inv.df <- data.frame(
  Mechanism = c("Interaction between resource fluctuations and competitor 1 resource intake", 
                'Competitor 1 resource intake', 
                "Predator preference", 
                "Resource fluctuations", 
                "Baseline"),
  IGR = C2_final_mechanisms
)

# fix order of interactions
inv.df$Mechanism <- factor(inv.df$Mechanism, levels = c("Interaction between resource fluctuations and competitor 1 resource intake", 
                                                'Competitor 1 resource intake', 
                                                "Predator preference", 
                                                "Resource fluctuations", 
                                                "Baseline"))

# calculate cumulative sum
inv.df$cumulative_IGR_backwards <- rev(cumsum(rev(inv.df$IGR)))

# # plot
# ggplot(inv.df, aes(x = Mechanism, y = IGR)) +
#   geom_bar(stat = "identity", width = 0.7, fill = "lightgrey", color = "black") +
#   coord_flip() +  
#   geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.5) +
#   geom_line(aes(x = as.numeric(Mechanism), y = cumulative_IGR_backwards, group = 1), color = "black", size = 1) +  # Add cumulative sum line starting from baseline
#   # geom_text(aes(x = as.numeric(Mechanism), y = cumulative_IGR_backwards, label = round(cumulative_IGR_backwards, 2)), 
#   #           color = "black", hjust = -0.8, size = 4.5) +  # Add labels next to each point on the cumulative sum line
#   geom_point(aes(x = as.numeric(Mechanism)[1], y = cumulative_IGR_backwards[1]), color = "red", size = 3) +  # Add red dot for cumulative sum of contributions representing the value of the invasion growth rate
#   theme_classic() +  
#   labs(x = "Mechanism", y = "Invasion Growth Rate") +
#   ggtitle(expression("Decomposition of C"[2]~"at -C"[2]~'Community'))+
#   scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) + 
#   scale_y_continuous(breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6)) +
#   theme(
#     axis.text = element_text(size = 12),
#     axis.title = element_text(size = 14, face = "bold"),   
#     plot.title = element_text(size = 14, face = "bold"), 
#     legend.position = "none" 
#   )

#----------------------------------------------------------------------------------------------------

# B. decomposition of C2 in -P comm

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

# timesteps
time <- seq(0,100000,by=0.1) 

# initial conditions
State <- c(P = 0, C1 = 1, C2 = 0, R = 1)

# create array for parameters
pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R)

# simulate dynamics
out <- ode(y = State, times = time, func = VassFox_Cvar, parms = pars)

# length of burn in period
invade_start_time <- 10000 #index to start invasions at

# initialize matrices to store values
C2_ldgr <- matrix(data=NA, nrow=(length(time)-invade_start_time), ncol=1)
C1_resident <- C2_ldgr
counter <- 1

# loop through each time step after burn in period and calculate instantaneous growth rate of both competitor species
for (t in invade_start_time:length(time)){
  
  P = out[t,2]
  C1 = out[t,3]
  C2 = out[t,4]
  
  R = out[t,5]
  
  C2_ldgr[counter] <- - (M_C2) + ( (O_C2_R * J_C2 * R) / (R + R_0_2) ) - ( ((1-O_P_C1) * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
  C1_resident[counter] <- - (M_C1) + ( (O_C1_R * J_C1 * R) / (R + R_0_1) ) - ( (O_P_C1 * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
  
  counter <- counter + 1
  
}

# average over all instantaneous growth rates for each species to calculate invaison growth rate. The resident's should be numerical zero. Calculate their difference. 
C2_r_bar <- mean(C2_ldgr)-mean(C1_resident)

# ----------------------------------------------------------------------------------------------------
# calculate delta_0. IGR under no resource fluctuation and no C1 intake of resource

# Constant resource abundance
avg_R = mean(out[,5])
R = avg_R

# No resource uptake by resident competitor
J_C1 = 0

# initialize matrices to hold values
C2_epsilon_0 <- matrix(data=NA, nrow=(length(time)-invade_start_time), ncol=1)
C1_resident_epsilon_0 <- C2_epsilon_0
counter <- 1

# loop through each time step after burn in period and calculate instantaneous growth rate of both competitors. 
for (t in invade_start_time:length(time)){
  
  P = out[t,2]
  C1 = out[t,3]
  C2 = out[t,4]
  
  C2_epsilon_0[counter] <- - (M_C2) + ( (O_C2_R * J_C2 * R) / (R + R_0_2) ) - ( ((1-O_P_C1) * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
  C1_resident_epsilon_0[counter] <- - (M_C1) + ( (O_C1_R * J_C1 * R) / (R + R_0_1) ) - ( (O_P_C1 * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
  
  counter <- counter + 1
  
}

# Average across instanteous growth rates to find invasion growth rates and take their difference. 
C2_delta_0 <- mean(C2_epsilon_0)-mean(C1_resident_epsilon_0)

# ----------------------------------------------------------------------------------------------------

# calculate delta_R, with R fluctuations but no C1 resource intake

# initialize matrices to hold values
C2_epsilon_R <- matrix(data=NA, nrow=(length(time)-invade_start_time), ncol=1)
C1_resident_epsilon_R <- C2_epsilon_R
counter <- 1

# no resource uptake by resident competitor
J_C1 = 0

# loop through each time step after burn in period and calculate instantaneous growth rate of both competitors 
for (t in invade_start_time:length(time)){
  
  # allow dynamics of resource 
  R = out[t,5]
  
  P = out[t,2]
  C1 = out[t,3]
  C2 = out[t,4]
  
  C2_epsilon_R[counter] <- - (M_C2) + ( (O_C2_R * J_C2 * R) / (R + R_0_2) ) - ( ((1-O_P_C1) * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
  C1_resident_epsilon_R[counter] <- - (M_C1) + ( (O_C1_R * J_C1 * R) / (R + R_0_1) ) - ( (O_P_C1 * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
  
  counter <- counter + 1
  
}

# calculate mean impact of resource fluctuations on invading species relative to the resident competitor.
C2_delta_R <- mean(C2_epsilon_R)-mean(C1_resident_epsilon_R) - C2_delta_0

# ----------------------------------------------------------------------------------------------------

# calculate delta_JC1, with C1 resource intake but no fluctuations in R

# initialize matrices to store values
C2_epsilon_JC1 <- matrix(data=NA, nrow=(length(time)-invade_start_time), ncol=1)
C1_resident_epsilon_JC1 <- C2_epsilon_JC1
counter <- 1

# resource uptake rate of resident competitor
J_C1 = 0.8036

# constant resource density
R = avg_R

# loop through each time step after burn in period and calculate instantaneous growth rate of both competitors
for (t in invade_start_time:length(time)){
  
  P = out[t,2]
  C1 = out[t,3]
  C2 = out[t,4]
  
  C2_epsilon_JC1[counter] <- - (M_C2) + ( (O_C2_R * J_C2 * R) / (R + R_0_2) ) - ( ((1-O_P_C1) * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
  C1_resident_epsilon_JC1[counter] <- - (M_C1) + ( (O_C1_R * J_C1 * R) / (R + R_0_1) ) - ( (O_P_C1 * J_P * P) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
  
  counter <- counter + 1
  
}

# calculate mean impact of resource uptake of resident species on invading competitor relative to the resident competitor. 
C2_delta_JC1 <- mean(C2_epsilon_JC1)-mean(C1_resident_epsilon_JC1) - C2_delta_0

# calculate two way intercations between resource fluctuations and resource uptake by resident competitor
C2_delta_R_JC1 = C2_r_bar - (C2_delta_0 + C2_delta_R + C2_delta_JC1)

# store mean values
C2_final_mechanisms <- c(C2_delta_R_JC1, C2_delta_JC1, C2_delta_R, C2_delta_0)

# Create dataframe for plotting
df <- data.frame(
  Mechanism = c("Interaction between resource fluctuations and competitor 1 resource intake", 
                'Competitor 1 resource intake', 
                "Resource fluctuations", 
                "Baseline"),
  IGR = C2_final_mechanisms
)

# fix order of interactions
df$Mechanism <- factor(df$Mechanism, levels = c("Interaction between resource fluctuations and competitor 1 resource intake", 
                                                'Competitor 1 resource intake', 
                                                "Resource fluctuations", 
                                                "Baseline"))

# create cumulative sum
df$cumulative_IGR_backwards <- rev(cumsum(rev(df$IGR)))

# # plot
# ggplot(df, aes(x = Mechanism, y = IGR)) +
#   geom_bar(stat = "identity", width = 0.7, fill = "lightgrey", color = "black") +
#   coord_flip() +  
#   geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.5) +
#   geom_line(aes(x = as.numeric(Mechanism), y = cumulative_IGR_backwards, group = 1), color = "black", size = 1) +  # Add cumulative sum line starting from baseline
#   # geom_text(aes(x = as.numeric(Mechanism), y = cumulative_IGR_backwards, label = round(cumulative_IGR_backwards, 2)), 
#   #           color = "black", hjust = -0.7, size = 4.5) +  # Add labels next to each point on the cumulative sum line
#   geom_point(aes(x = as.numeric(Mechanism)[1], y = cumulative_IGR_backwards[1]), color = "red", size = 3) +  # Add red dot for cumulative sum of contributions representing the value of the invasion growth rate
#   theme_classic() +  
#   labs(x = "Mechanism", y = "Invasion Growth Rate") +
#   ggtitle(expression("Decomposition of C"[2]~"at -P community")) +
#   scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) + 
#   scale_y_continuous(breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6)) +
#   theme(
#     axis.text = element_text(size = 12),
#     axis.title = element_text(size = 14, face = "bold"),   
#     plot.title = element_text(size = 14, face = "bold"), 
#     legend.position = "none" 
#   )

###################################################################################################################

# create dataframe for final plotting
total.df = as.data.frame(inv.df$Mechanism)
colnames(total.df) = 'Mechanism'
total.df$minusi.comm.IGR = inv.df$IGR

# find mechanisms shared between two communities
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
mechanism_levels <- c("Cumulative IGR", setdiff(unique(total.df$Mechanism), c("Cumulative IGR", "Baseline and other weak contributions")), "Baseline and other weak contributions")
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
  scale_x_discrete(labels = function(x) str_wrap(x, width = 23)) + 
  scale_fill_manual(values = c("minusi.comm.IGR" = "#fec488", "secext.comm.IGR" = "#e44f64"),  
                    labels = c(expression(-C[2] ~ "   "), expression(-P))) +
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18, face = "bold"),   
    plot.title = element_text(size = 18, face = "bold"), 
    legend.position = "top",
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18),
  )
