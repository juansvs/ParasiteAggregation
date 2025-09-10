#--- set parameters

## arena setup
Na = 5 # N animals
N_patches =  78*78 # N patches
total_time = 60*24*1 # total simulation time (in minutes) 
h_max =  400 # maximum sward height (logistic growth)
h0 = 50 # ungrazeable portion of sward
s0 = 2000 # units of feces deposited per defecation event
B = 1 # bite size

## patch parameters
gamma = 0.00004 # sward growth rate
epsilon = 5e-5 # development L -> L3
omega = 1e-4 # death rate of pre-infective larvae
rho = 1.5e-5 # death rate of L3 larvae
phi = 1.776e-5 # feces decay rate

## animal parameters
beta = 0.1  # per-capita bite rate
mu_f = 5 # livestock feces avoidance parameter
zeta = 2e-5 # death of immature larvae in host (this not specified in the paper)
chi = 3e-5 # maturation rate of larvae in host
sig = 1.9e-8 # rate of resistance loss
psi = 0.25 # resistance gain coefficient 1
eta = 0.025 # resistance gain coefficient 2
tau = 2e-5 # death rate of larvae in host
lambda = 2 # rate of egg production
Lambda = 0 # anorexia coefficient
nu = 0.015 # intrinsic search/movement rate
theta = 0.4 # probability of L3 establishment as adults
f_dep = 1 # defecation rate
alpha = 10 #  power law search coefficient


# initial state for infective larvae. In original paper 24000 larvae distributed
# across 20 patches in a 78x78 arena
larvae_patches <- sample(1:N_patches, 20)
L_counts <- sample(larvae_patches, 24000, replace = T) |> table()
L <- rep(0, N_patches)
L[as.numeric(names(L_counts))] <- L_counts

# Patch variables
h = rep(h_max / 2, N_patches) # Sward height in each patch
f = rep(0, N_patches)       # Livestock faeces
l = rep(0, N_patches)       # uninfective larvae in pasture
L = L      # infective larvae (L3) in pasture
# Host variables
animal_locations = sample(1:N_patches, Na, replace = TRUE)  # Place animals randomly on the grid
r = rep(0, Na)              # Immune response of each host
a = rep(0, Na)              # Number of immature parasites in host
A = rep(0, Na)              # Number of adult parasites in host  
eg = rep(0, Na)             # Number of eggs in host
s = rep(0, Na)              # Stomach content per animal
