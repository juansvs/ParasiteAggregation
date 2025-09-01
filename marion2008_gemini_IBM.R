#
# --- Main Simulation Function ---
#
run_grazing_model <- function(
    Na, N_patches, total_time, h_max, beta, h0, mu_f, mu_w,
    lambda_w, lambda_f, s0, nu, alpha, seed = NULL
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
  for (loc in animal_locations) {
    c[loc] <- c[loc] + 1
  }

  current_time <- 0
  time_series <- data.frame(time = 0, avg_height = mean(h), avg_stomach = mean(s))

  # Main simulation loop
  while (current_time < total_time) {

    # 1. Calculate event rates for all possible events
    rates <- get_event_rates(
      h, c, f, w, s, animal_locations, Na, N_patches,
      h_max, beta, h0, mu_f, mu_w, lambda_w, lambda_f, s0, nu, alpha
    )
    total_rate <- sum(rates$all_rates)

    # 2. Determine time to next event (exponential distribution)
    if (total_rate == 0) {
      break # No events left to occur
    }
    delta_t <- rexp(1, rate = total_rate)

    # 3. Choose which event occurs
    event_index <- sample(
      1:length(rates$all_rates),
      size = 1,
      prob = rates$all_rates / total_rate
    )

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
      animal_idx <- event_patch_or_animal_idx
      current_patch <- animal_locations[animal_idx]
      # Select new patch based on search kernel
      
      # For simplicity, we'll implement the global search model from the text
      # A better implementation would pre-calculate probabilities
      possible_moves <- setdiff(1:N_patches, current_patch)
      
      # Calculate power law probabilities for all other patches
      distances <- sqrt((floor((current_patch - 1) / sqrt(N_patches)) - floor((possible_moves - 1) / sqrt(N_patches)))^2 + 
                        ((current_patch - 1) %% sqrt(N_patches) - (possible_moves - 1) %% sqrt(N_patches))^2)
      probs <- distances^(-alpha)
      
      if (sum(probs) > 0) {
          new_patch <- sample(possible_moves, size = 1, prob = probs)
          c[current_patch] <- c[current_patch] - 1
          c[new_patch] <- c[new_patch] + 1
          animal_locations[animal_idx] <- new_patch
      }
    }

    # 5. Advance time and record state
    current_time <- current_time + delta_t
    time_series <- rbind(time_series, data.frame(time = current_time, avg_height = mean(h), avg_stomach = mean(s)))
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

  # 1. Grazing rates for each animal
  grazing_rates <- beta * (h[animal_locations] - h0) * exp(-mu_f * f[animal_locations] - mu_w * w[animal_locations])
  grazing_rates[grazing_rates < 0] <- 0 # Ensure no negative rates
  all_rates <- c(all_rates, grazing_rates)
  event_types <- c(event_types, rep("grazing", Na))
  event_indices <- c(event_indices, 1:Na)

  # 2. Defecation rates for each animal
  defecation_rates <- (nu * (s >= s0)) / s0 # Simplified since intake and faeces are in same units
  all_rates <- c(all_rates, defecation_rates)
  event_types <- c(event_types, rep("defecation", Na))
  event_indices <- c(event_indices, 1:Na)

  # 3. Movement rates for each animal
  movement_rates <- rep(nu, Na)
  all_rates <- c(all_rates, movement_rates)
  event_types <- c(event_types, rep("movement", Na))
  event_indices <- c(event_indices, 1:Na)

  # 4. Sward growth rates for each patch
  growth_rates <- h * (1 - h / h_max)
  all_rates <- c(all_rates, growth_rates)
  event_types <- c(event_types, rep("growth", N_patches))
  event_indices <- c(event_indices, 1:N_patches)

  # 5. Faecal decay rates for each patch
  f_decay_rates <- lambda_f * f
  w_decay_rates <- lambda_w * w
  all_rates <- c(all_rates, f_decay_rates, w_decay_rates)
  event_types <- c(event_types, rep("f_decay", N_patches), rep("w_decay", N_patches))
  event_indices <- c(event_indices, 1:N_patches, 1:N_patches)

  return(list(all_rates = all_rates, event_types = event_types, event_indices = event_indices))
}