#
# --- Main Simulation Function ---
#
run_grazing_model <- function(
    Na, N_patches, total_time, beta, gamma, h_max, h0, s0, mu_f, mu_w,
    lambda_w, lambda_f, nu, alpha, seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Initialize state variables
  h <- rep(h_max / 2, N_patches) # Sward height in each patch
  f <- rep(0, N_patches)        # Livestock faeces
  w <- rep(0, N_patches)        # Wildlife faeces (assumed external/static for this example)
  c <- rep(0, N_patches)        # Animal count per patch
  s <- rep(0, Na)              # Stomach content per animal

  # Place animals randomly on the grid
  animal_locations <- sample(1:N_patches, Na, replace = TRUE)
  c[animal_locations] <- c[animal_locations]+1

  current_time <- 0
  n_feces_investig <- 0
  time_series <- data.frame(time = numeric(total_time%/%5+1), avg_height = 0, avg_contamination = 0, feces_invest = 0)
  time_series[1,] <- data.frame(time = current_time, avg_height = mean(h), avg_contamination = mean(f), feces_invest = n_feces_investig)
  ix <- 1
  
  # Main simulation loop
  while (current_time < total_time) {

    # 1. Calculate event rates for all possible events
    rates <- get_event_rates(
      h, c, f, w, s, animal_locations, Na, N_patches,
      h_max, beta, h0, mu_f, mu_w, lambda_w, lambda_f, s0, nu, alpha
    )
    total_rate <- sum(rates$all_rates)

    # 2. Determine time to next event (exponential distribution)
    # if (total_rate == 0) {
    #   break # No events left to occur
    # }
    rnums <- runif(2)
    delta_t <- -log(rnums[1])/total_rate

    # 3. Choose which event occurs
    event_probs <- cumsum(rates$all_rates/total_rate)
    event_index <- 1+sum(rnums[2]>event_probs)

    # 4. Update state variables based on the chosen event
    event_type <- rates$event_types[event_index]
    event_patch_or_animal_idx <- rates$event_indices[event_index]
    
    # Apply updates based on event type
    if (event_type == "growth") {
      h[event_patch_or_animal_idx] <- h[event_patch_or_animal_idx] + 1
    } else if (event_type == "grazing") {
      animal_idx <- event_patch_or_animal_idx
      patch_idx <- animal_locations[animal_idx]
      if (h[patch_idx] > h0) {
        h[patch_idx] <- h[patch_idx] - 1
        s[animal_idx] <- s[animal_idx] + 1
      }
      if(f[patch_idx]>0) n_feces_investig <- n_feces_investig + 1
    } else if (event_type == "f_decay") {
      if (f[event_patch_or_animal_idx] > 0) {
        f[event_patch_or_animal_idx] <- f[event_patch_or_animal_idx] - 1
      }
    } else if (event_type == "w_decay") {
      if (w[event_patch_or_animal_idx] > 0) {
        w[event_patch_or_animal_idx] <- w[event_patch_or_animal_idx] - 1
      }
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
      current_patch <- animal_locations[animal_idx]
      c[current_patch] <- c[current_patch]-1
      c[dest_patch] <- c[dest_patch]+1
      animal_locations[animal_idx] <- dest_patch
    }

    # 5. Advance time and record state every 5 minutes
    new_time <- current_time+delta_t
    if(floor(new_time)%%5==0) {
      time_series[ix,] <- data.frame(time = current_time, avg_height = mean(h), avg_contamination = mean(f), feces_invest = n_feces_investig)
      ix <- ix+1
    }
    current_time <- new_time
    
  }
  return(time_series)
}

#
# --- Helper Function for Rate Calculation ---
#
get_event_rates <- function(h, c, f, w, s, animal_locations, Na, N_patches,
                           h_max, beta, h0, mu_f, mu_w, lambda_w, lambda_f, s0, nu, alpha) {

  all_rates <- c()
  event_types <- c()
  event_indices <- c()
  destination <- c()

  # 1. Grazing rates for each animal
  grazing_rates <- beta * (h[animal_locations] - h0) * exp(-mu_f * f[animal_locations] - mu_w * w[animal_locations])
  all_rates <- c(all_rates, grazing_rates)
  event_types <- c(event_types, rep("grazing", Na))
  event_indices <- c(event_indices, 1:Na)
  destination <- c(destination, rep(NA, Na))
  
  # 2. Defecation rates for each animal
  defecation_rates <- f_dep*(s-s0)*as.numeric(s>s0)
  all_rates <- c(all_rates, defecation_rates)
  event_types <- c(event_types, rep("defecation", Na))
  event_indices <- c(event_indices, 1:Na)
  destination <- c(destination, rep(NA, Na))
  
  # 3. Movement rates for each animal
  movement_rates <- lapply(animal_locations, mov_rate, hj = h, nu = nu, alpha = alpha, rw = sqrt(N_patches), cl = sqrt(N_patches))
  all_rates <- c(all_rates, unlist(movement_rates))
  event_types <- c(event_types, rep("movement", Na*N_patches))
  event_indices <- c(event_indices, rep(1:Na, each = N_patches))
  destination <- c(destination, rep(1:N_patches, Na))

  # 4. Sward growth rates for each patch
  growth_rates <- gamma * h * (1 - h / h_max)
  all_rates <- c(all_rates, growth_rates)
  event_types <- c(event_types, rep("growth", N_patches))
  event_indices <- c(event_indices, 1:N_patches)
  destination <- c(destination, rep(NA, N_patches))
  
  # 5. Faecal decay rates for each patch
  f_decay_rates <- lambda_f * f
  w_decay_rates <- lambda_w * w
  all_rates <- c(all_rates, f_decay_rates, w_decay_rates)
  event_types <- c(event_types, rep("f_decay", N_patches), rep("w_decay", N_patches))
  event_indices <- c(event_indices, 1:N_patches, 1:N_patches)
  destination <- c(destination, rep(NA, 2*N_patches))
  
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
