## To do: modify rates maybe create a setup function

#--- set parameters
pars <- list(
  ## arena setup
  Na = 5 ,            # N animals
  N_patches =  78*78, # N patches
  total_time = 60*12*1,# total simulation time (in minutes) 
  h_max =  400 ,       # maximum sward height (logistic growth)
  h0 = 50 ,           # ungrazeable portion of sward
  s0 = 2000 ,         # units of feces deposited per defecation event
  B = 1 ,             # bite size
  
  
  ## patch parameters
  gamma = 0.00004 ,# sward growth rate
  epsilon = 5e-5 , # development L -> L3
  omega = 1e-4 ,   # death rate of pre-infective larvae
  rho = 1.5e-5 ,   # death rate of L3 larvae
  phi = 1.776e-5,  # feces decay rate
  
  ## animal parameters
  beta = 0.1,      # per-capita bite rate
  mu_f = 5  ,      # livestock feces avoidance parameter
  zeta = 2e-5 ,    # death of immature larvae in host (this not specified in the paper)
  chi = 3e-5  ,    # maturation rate of larvae in host
  sig = 1.9e-8,    # rate of resistance loss
  psi = 0.25  ,    # resistance gain coefficient 1
  eta = 0.025 ,# resistance gain coefficient 2
  tau = 2e-5 ,# death rate of larvae in host
  lambda = 2 ,# rate of egg production
  Lambda = 0 ,# anorexia coefficient
  theta = 0.4 ,# probability of L3 establishment as adults
  f_dep = 1 ,# defecation rate
  
  ## movement process
  nu = 0.015 ,# intrinsic search/movement rate
  alpha = 10 #  power law search coefficient
)

# initial state for infective larvae. In original paper 24000 larvae distributed
# across 20 patches in a 78x78 arena
larvae_patches <- sample(1:N_patches, 20)
L_counts <- sample(larvae_patches, 24000, replace = T) |> table()
L <- rep(0, N_patches)
L[as.numeric(names(L_counts))] <- L_counts

S <- list(
  # Patch variables
  h = rep(h_max / 2, N_patches), # Sward height in each patch
  f = rep(0, N_patches) ,      # Livestock faeces
  l = rep(0, N_patches) ,      # uninfective larvae in pasture
  L = L      ,# infective larvae (L3) in pasture
  # Host variables
  animal_locations = sample(1:N_patches, Na, replace = TRUE) , # Place animals randomly on the grid
  r = rep(0, Na)  ,            # Immune response of each host
  a = rep(0, Na) ,             # Number of immature parasites in host
  A = rep(0, Na),              # Number of adult parasites in host  
  eg = rep(0, Na),             # Number of eggs in host
  s = rep(0, Na)              # Stomach content per animal
)
## Event rates
# sward growth
growth_rates <- gamma * h * (1 - h / h_max)
# development of larvae in patches
dev_l <- epsilon * l
# death of pre-infective larvae
death_l <- omega * l
# death of infective larvae
death_L <- rho * L
# Fecal decay rates for each patch
f_decay <- phi * f

# Grazing rates for each animal
grazing_rates <- beta * (h[animal_locations] - h0) * exp(-mu_f * f[animal_locations] *(a+A)^Lambda)
# death of immature adults in host
death_a <- zeta * a
# development into adult parasites
dev_a <- chi * a
# death of adult parasites
death_A <- tau * A
# gain of immunity due to parasite burden
immun_gain <- (a + A) * eta
# loss of immunity
immun_loss <- sig * r
# egg production
egg_prod <- lambda * A/2
# Defecation rates for each animal
defecation_rates <- f_dep*(s-s0)*as.numeric(s>s0)
# Movement rates for each animal
movement_rates <- sapply(animal_locations, mov_rate, hj = h, nu = nu, alpha = alpha, rw = sqrt(N_patches), cl = sqrt(N_patches))
