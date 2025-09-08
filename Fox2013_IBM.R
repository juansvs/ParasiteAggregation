#### TO DO ####

# Incorporate parasite migration (e.g. Gloworm)

#
#
#### --- Main Simulation Function --- ####
#
run_model <- function(
    pars,
    # Na, N_patches, total_time, beta, gamma, h_max, h0, s0, mu_f, 
    # lambda_f, nu, alpha, theta, 
    seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Initialize state variables
  # Patch variables
  h <- rep(h_max / 2, N_patches) # Sward height in each patch
  f <- rep(0, N_patches)       # Livestock faeces
  l <- rep(0, N_patches)       # uninfective larvae in pasture
  L <- rep(0, N_patches)       # infective larvae (L3) in pasture
  # Host variables
  animal_locations <- sample(1:N_patches, Na, replace = TRUE)  # Place animals randomly on the grid
  r <- rep(0, Na)              # Immune response of each host
  a <- rep(0, Na)              # Number of immature parasites in host
  A <- rep(0, Na)              # Number of adult parasites in host  
  eg <- rep(0, Na)             # Number of eggs in host
  s <- rep(0, Na)              # Stomach content per animal


  # Set up output dataframe
  current_time <- 0
  n_feces_investig <- 0
  time_series <- data.frame(time = numeric(total_time%/%30+1), avg_height = 0, avg_contamination = 0, feces_invest = 0)
  time_series[1,] <- data.frame(time = current_time, avg_height = mean(h), avg_contamination = mean(f), feces_invest = n_feces_investig)
  ix <- 2
  
  # Main simulation loop
  while (current_time < total_time) {
    # Create a bunch of random numbers at once to save computation time
    rnums <- matrix(runif(2e7), ncol = 2)
    for (y in 1:1e7) {
      # 1. Calculate event rates for all possible events
      rates <- get_event_rates(states, pars)
      total_rate <- sum(rates$all_rates)
      
      # 2. Determine time to next event (exponential distribution)
      # if (total_rate == 0) {
      #   break # No events left to occur
      # }
      delta_t <- -log(rnums[y,1])/total_rate
      
      # 3. Choose which event occurs
      event_probs <- cumsum(rates$all_rates/total_rate)
      event_index <- 1+sum(rnums[y,2]>event_probs)
      
      # 4. Update state variables based on the chosen event
      event_type <- rates$event_types[event_index]
      event_patch_or_animal_idx <- rates$event_indices[event_index]
      
      # Apply updates based on event type
      states <- update_state(event_type)
      
      # 5. Advance time and record state every 30 minutes
      new_time <- current_time+delta_t
      record_state <- all(floor(new_time)>floor(current_time),floor(new_time)%%30==0)
      if(record_state) {
        time_series[ix,] <- data.frame(time = current_time, avg_height = mean(h), avg_contamination = mean(f), feces_invest = n_feces_investig)
        ix <- ix+1
      }
      current_time <- new_time
    }
  }
  return(time_series)
}

#
#### --- Helper Functions --- ####
#
get_event_rates <- function(states, pars) {
  h <- states$h;  f <- states$f;  l <- states$l
  L <- states$L;  s <- states$s;  a <- states$a
  A <- states$A;  eg <- states$eg;  r <- states$r
  animal_locations <- states$animal_locations
  Na <- pars$Na;  N_patches <- pars$N_patches
  h_max <- pars$h_max;  beta <- pars$beta
  h0 <- pars$h0;  mu_f <- pars$mu_f
  lambda_f <- pars$lambda_f
  s0 <- pars$s0;  nu <- pars$nu
  alpha <- pars$alpha;  gamma <- pars$gamma
  theta <- pars$theta;  omega <- pars$omega
  eta <- pars$eta;  rho <- pars$rho
  phi <- pars$phi;  chi <- pars$chi
  tau <- pars$tau;  sig <- pars$sig
  Lambda <- pars$Lambda;  zeta <- pars$zeta

  all_rates <- numeric(N_patches*5+Na*8+N_patches*Na)
  event_types <- c(rep(c("growth","development_l","death_l",
                       "death_L","f_decay"),each = N_patches),
                   rep(c("bite","death_a",
                       "development_a","death_A","immunity_gain",
                       "immunity_loss","egg_prod","defecation"), each = Na),
                   rep("movement",N_patches*Na))
  event_indices <- c(rep(1:N_patches,5), rep(1:Na,8),c(1:Na, each = N_patches*Na))
  destination <- c(rep(NA,N_patches*5+Na*8), rep(1:N_patches, Na))

  # 1. Grazing rates for each animal
  grazing_rates <- beta * (h[animal_locations] - h0) * exp(-mu_f * f[animal_locations] *(a+A)*Lambda)

  # 2. Defecation rates for each animal
  defecation_rates <- f_dep*(s-s0)*as.numeric(s>s0)

  # 3. Movement rates for each animal
  movement_rates <- lapply(animal_locations, mov_rate, hj = h, nu = nu, alpha = alpha, rw = sqrt(N_patches), cl = sqrt(N_patches))

  # 4. Sward growth rates for each patch
  growth_rates <- gamma * h * (1 - h / h_max)

  # 5. Faecal decay rates for each patch
  f_decay_rates <- lambda_f * f

  #
  
  return(list(all_rates = all_rates, event_types = event_types, event_indices = event_indices, destinations = destination))
}


# function to calculate the rates at which an individual in patch i moves to
# patch j
mov_rate <- function(cell_curr, hj, nu, alpha, rw, cl) {
  index_matrix <- matrix(nrow = rw, ncol = cl)
  rowmat <- row(index_matrix)
  colmat <- col(index_matrix)
  # get x,y coordinates
  curr_row <- row(index_matrix)[cell_curr]
  curr_col <- col(index_matrix)[cell_curr]
  # get distances in x and y from the rest of the grid
  xdist <- curr_col-colmat
  ydist <- curr_row-rowmat
  # power law using the Euclidean distance, and the alpha value
  F_i <- sqrt(xdist^2+ydist^2)^-alpha
  # substitute the current cell value by 0
  F_i[cell_curr] <- 0
  z_i <- sum(F_i)
  rate <- nu/z_i*F_i*hj
  return(rate)
}

# function to update state space based on event type. 
update_state <- function(event_type) {
  if (event_type == "growth") {
    h[event_patch_or_animal_idx] <- h[event_patch_or_animal_idx] + 1
  } else if (event_type == "development_l") {
    l[event_patch_or_animal_idx] <- l[event_patch_or_animal_idx]-1
    L[event_patch_or_animal_idx] <- L[event_patch_or_animal_idx]+1
  } else if(event_type == "death_l") {
    l[event_patch_or_animal_idx] <- l[event_patch_or_animal_idx]-1
  } else if(event_type == "death_L") {
    L[event_patch_or_animal_idx] <- L[event_patch_or_animal_idx]-1
  } else if (event_type == "f_decay") {
    f[event_patch_or_animal_idx] <- f[event_patch_or_animal_idx] - 1
  } else if (event_type == "bite") {
    animal_idx    <- event_patch_or_animal_idx
    patch_idx     <- animal_locations[animal_idx]
    h[patch_idx]  <- h[patch_idx]  - 1   # reduce patch sward height
    s[animal_idx] <- s[animal_idx] + 1  # increase stomach content
    l[patch_idx]  <- l[patch_idx]  - B/h[patch_idx]*l[patch_idx] # reduce number of larve in patch
    L[patch_idx]  <- L[patch_idx]  - B/h[patch_idx]*L[patch_idx] # reduce number of larve in patch
    a[animal_idx] <- a[animal_idx] + theta*r[animal_idx]*(B/h[patch_idx])*L[patch_idx] # increase number of larvae in host
    r[animal_idx] <- r[animal_idx] + phi*B*L[patch_idx]/h[patch_idx] # update host resistance
  } else if(event_type == "death_a") {
    a[event_patch_or_animal_idx] <- a[event_patch_or_animal_idx] - 1
  } else if(event_type == "development_a") {
    a[event_patch_or_animal_idx] <- a[event_patch_or_animal_idx] - 1
    A[event_patch_or_animal_idx] <- A[event_patch_or_animal_idx] + 1
  } else if(event_type == "death_A") {
    A[event_patch_or_animal_idx] <- A[event_patch_or_animal_idx] - 1
  } else if(event_type == "immunity_gain") {
    r[event_patch_or_animal_idx] <- r[event_patch_or_animal_idx] + 1
  } else if(event_type == "immunity_loss") {
    r[event_patch_or_animal_idx] <- r[event_patch_or_animal_idx] - 1
  } else if(event_type == "egg_prod") {
    eg[event_patch_or_animal_idx] <- eg[event_patch_or_animal_idx] + 1
  } else if (event_type == "defecation") {
    animal_idx <- event_patch_or_animal_idx
    patch_idx <- animal_locations[animal_idx]
    # Heaviside function (Theta(s_k - s0))
    if (s[animal_idx] >= s0) {
      f[patch_idx] <- f[patch_idx] + s0
      s[animal_idx] <- s[animal_idx] - s0
    }
  } else if (event_type == "movement") {
    dest_patch <- rates$destinations[event_index]
    animal_idx <- event_patch_or_animal_idx
    # current_patch <- animal_locations[animal_idx]
    # c[current_patch] <- c[current_patch]-1
    # c[dest_patch] <- c[dest_patch]+1
    animal_locations[animal_idx] <- dest_patch
  } 
}

