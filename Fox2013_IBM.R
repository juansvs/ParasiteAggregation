#### TO DO ####

# Incorporate parasite migration (e.g. Gloworm)
# check functions of resistance: rates lambda, tau, and change theta

#
#
#### --- Main Simulation Function --- ####
#
run_model <- function(seed = NULL, outf) {
  starttime <- Sys.time()
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Set up output dataframe
  current_time <- 0
  time_series <- data.frame(time = numeric(total_time%/%30+1), avg_height = 0, avg_a = 0, avg_A = 0, sd_A = 0, avg_l = 0, avg_L = 0)
  time_series[1,] <- c(current_time, mean(h), mean(a), mean(A), sd(A), mean(l), mean(L))
  ix <- 2
  cat(c("time","sim_time", "avg_h", "avg_a", "avg_A", "sd_A", "avg_l", "avg_L", "\n"), sep = "\t", file = outf)
  cat(c(0,current_time, mean(h), mean(a), mean(A), mean(l), mean(L)),"\n", sep = "\t", file = outf, append = T)
  
  #rates list object
  event_db <- list(event_types = c(rep(c("growth","development_l","death_l",
                                           "death_L","f_decay"),each = N_patches),
                                     rep(c("bite","death_a",
                                           "development_a","death_A","immunity_gain",
                                           "immunity_loss","egg_prod","defecation"), each = Na),
                                     rep("movement",N_patches*Na)),
                     event_indices = c(rep(1:N_patches,5), rep(1:Na,8),rep(1:Na, each = N_patches)),
                     destination = c(rep(NA,N_patches*5+Na*8), rep(1:N_patches, Na)))
  
  
  # Main simulation loop
  while (current_time < total_time) {
    # Create a bunch of random numbers at once to save computation time
    rnums <- matrix(runif(2e4), ncol = 2)
    for (y in 1:nrow(rnums)) {
      # 1. Calculate event rates for all possible events
      all_rates <- if(current_time > 0) get_event_rates(type = 'sort',event_index,event_db) else get_event_rates(type='exact')
      total_rate <- sum(all_rates)
      
      # 2. Determine time to next event (exponential distribution)
      # if (total_rate == 0) {
      #   break # No events left to occur
      # }
      delta_t <- -log(rnums[y,1])/total_rate
      
      # 3. Choose which event occurs
      event_probs <- cumsum(all_rates/total_rate)
      event_index <- 1+sum(rnums[y,2]>event_probs)
      
      # 4. Update state variables based on the chosen event
      update_state(event_index, event_db)
      
      # 5. Advance time and record state every 30 minutes
      new_time <- current_time+delta_t
      record_state <- all(floor(new_time)>floor(current_time),floor(new_time)%%30==0)
      if(record_state) {
        time_series[ix,] <- c(new_time, mean(h), mean(a), mean(A), sd(A), mean(l), mean(L))
        elapsedtime <- as.numeric(difftime(Sys.time(), starttime), units = "hours")
        cat(c(elapsedtime, new_time, mean(h), mean(a), mean(A), sd(A), mean(l), mean(L)), "\n", sep = "\t", file = outf, append = T)
        ix <- ix+1
      }
      current_time <- new_time
      if(current_time>=total_time) break
    }
  }
  return(time_series)
}

#
#### --- Helper Functions --- ####
#
get_event_rates0 <- function(pars = NULL, S = NULL) {
  with(c(pars,S),{ 
    # sward growth
    growth <- gamma * h * (1 - h / h_max)
    # development of larvae in patches
    dev_l <- epsilon * l
    # death of pre-infective larvae
    death_l <- omega * l
    # death of infective larvae
    death_L <- rho * L
    # Fecal decay rates for each patch
    f_decay <- phi * f
    
    # Grazing rates for each animal
    grazing <- beta * (h[animal_locations] - h0) * exp(-mu_f * f[animal_locations] *(a+A)^Lambda)
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
    defecation <- f_dep*(s-s0)*as.numeric(s>s0)
    # Movement rates for each animal
    movement <- sapply(animal_locations, mov_rate, hj = h, nu = nu, alpha = alpha, rw = sqrt(N_patches), cl = sqrt(N_patches))
    
    
    #create output lists
    rates <- list(growth = growth, dev_l = dev_l, death_l = death_l, death_L = death_L, f_decay = f_decay, 
                grazing = grazing, 
                death_a = death_a, dev_a = dev_a, death_A = death_A, immun_gain = immun_gain, immun_loss = immun_loss, egg_prod = egg_prod, 
                defecation = defecation, 
                movement = movement)
    times <- lapply(rates, rexp)
    
    return(list(rates = rates, 
                times = times))
  }
  )
}


# function to get the event rates using the Optimized method (sorting method)
# that does not recalculate every rate, but rather only those affected by the
# latest event
get_event_rates_opt <- function(event_type, event_index, rates, pars, S) {
  with(c(pars, S, rates), {
    # only update based on previous event
    if(event_type == "growth") {
      rates$growth[event_index] <- gamma * h[event_index] * (1 - h[event_index] / h_max)
      times$growth[event_index] <- rexp(1,rates$growth[event_index])
      # check which animals are in the patch
      animals_in_patch <- which(animal_locations==event_index)
      rates$grazing[animals_in_patch] <- beta * (h[event_index] - h0) * exp(-mu_f * f[event_index] *(a[animals_in_patch]+A[animals_in_patch])^Lambda)
      times$grazing[animals_in_patch] <- rexp(length(animals_in_patch),rates$grazing[animals_in_patch])
    } else if (event_type=="dev_l") {
      rates$dev_l[event_index]   <- epsilon * l[event_index]
      times$dev_l[event_index]   <- rexp(1,rates$dev_l[event_index])
      rates$death_l[event_index] <- omega * l[event_index]
      times$death_l[event_index] <- rexp(1,rates$death_l[event_index])
      rates$death_L[event_index] <- rho * L[event_index]
      times$death_L[event_index] <- rexp(1,rates$death_L[event_index])
    } else if (event_type=="death_l") {
      rates$death_l[event_index] <- omega * l[event_index]
      times$death_l[event_index] <- rexp(1,rates$death_l[event_index])
      rates$dev_l[event_index]   <- epsilon * l[event_index]
      times$dev_l[event_index]   <- rexp(1,rates$dev_l[event_index])
    } else if (event_type=="death_L") {
      rates$death_L[event_index] <- rho * L[event_index]
      times$death_L[event_index] <- rexp(1,rates$death_L[event_index])
    } else if (event_type=="f_decay") {
      rates$f_decay[event_index] <- phi * f_decay[event_index]
      times$f_decay[event_index] <- rexp(1, rates$f_decay[event_index])
      animals_in_patch <- which(animal_locations==event_index)
      rates$grazing[animals_in_patch] <- beta * (h[event_index] - h0) * exp(-mu_f * f[event_index] *(a[animals_in_patch]+A[animals_in_patch])^Lambda)
      times$grazing[animals_in_patch] <- rexp(length(animals_in_patch),rates$grazing[animals_in_patch])
    } else if (event_type=="bite") {
      patch <- animal_locations[event_index]
      rates$grazing[event_index] <- beta * (h[patch] - h0) * exp(-mu_f * f[patch] *(a[event_index]+A[event_index])^Lambda)
      times$grazing[event_index] <- rexp(1,rates$grazing[event_index])
      rates$growth[patch]  <- gamma * h[patch] * (1 - h[patch] / h_max)
      times$growth[patch] <- rexp(1, rates$growth[patch])
      rates$death_a[event_index] <- zeta * a[event_index]
      times$death_a[event_index] <- rexp(1,rates$death_a[event_index])
      rates$dev_a[event_index]   <- chi * a[event_index]
      times$dev_a[event_index] <- rexp(1,rates$dev_a[event_index])
      rates$dev_l[event_index]   <- epsilon * l[patch]
      times$dev_l[event_index] <- rexp(1,rates$dev_l[event_index])
      rates$death_l[event_index] <- omega * l[patch]
      times$death_l[event_index] <- rexp(1,rates$death_l[event_index])
      rates$death_L[event_index] <- rho * L[patch]
      times$death_L[event_index] <- rexp(1,rates$death_L[event_index])
    } else if (event_type=="death_a") {
      rates$death_a[event_index] <- zeta * a[event_index]
      times$death_a[event_index] <- rexp(1,rates$death_a[event_index])
      rates$dev_a[event_index]   <- chi * a[event_index]
      times$dev_a[event_index] <- rexp(1,rates$dev_a[event_index])
      rates$immun_gain[event_index] <- (a[event_index] + A[event_index]) * eta
      times$immun_gain[event_index] <- rexp(1, rates$immun_gain[event_index])
    } else if (event_type=="development_a") {
      rates$death_a[event_index]  <- zeta * a[event_index]
      times$death_a[event_index] <- rexp(1,rates$death_a[event_index])
      rates$dev_a[event_index]    <- chi * a[event_index]
      times$dev_a[event_index] <- rexp(1,rates$dev_a[event_index])
      rates$death_A[event_index]  <- tau * A[event_index]
      times$death_A[event_index] <- rexp(1, rates$death_A[event_index])
      rates$egg_prod[event_index] <- lambda * A[event_index]/2
      times$egg_prod[event_index] <- rexp(1, rates$egg_prod[event_index])
      rates$immun_gain[event_index] <- (a[event_index] + A[event_index]) * eta
      times$immun_gain[event_index] <- rexp(1, rates$immun_gain[event_index])
    } else if (event_type=="death_A") {
      rates$death_A[event_index]   <- tau * A[event_index]
      times$death_A[event_index] <- rexp(1, rates$death_A[event_index])
      rates$egg_prod[event_index]  <- lambda * A[event_index]/2
      times$egg_prod[event_index] <- rexp(1, rates$egg_prod[event_index])
      rates$immun_gain[event_index] <- (a[event_index] + A[event_index]) * eta
      rates$immun_gain[event_index] <- (a[event_index] + A[event_index]) * eta
    } else if (event_type=="immunity_loss") {
      rates$immun_loss[event_index] <- sig * r[event_index]
      times$immun_loss[event_index] <- rexp(1, rates$immun_loss[event_index])
    } else if (event_type=="defecation") {
      rates$defecation[event_index] <- f_dep*(s[event_index]-s0)*as.numeric(s[event_index]>s0)
      times$defecation[event_index] <- rexp(1, rates$defecation[event_index])
      patch <- animal_locations[event_index]
      rates$grazing[event_index] <- beta * (h[patch] - h0) * exp(-mu_f * f[patch] *(a[event_index]+A[event_index])^Lambda)
      times$grazing[event_index] <- rexp(1,rates$grazing[event_index])
      rates$dev_l[event_index]   <- epsilon * l[patch]
      times$dev_l[event_index] <- rexp(1,rates$dev_l[event_index])
      rates$death_l[event_index] <- omega * l[patch]
      times$death_l[event_index] <- rexp(1,rates$death_l[event_index])
    } else if (event_type=="movement") {
      rates$movement[,event_index] <- mov_rate(animal_locations[event_index], hj = h, nu = nu, alpha = alpha, rw = sqrt(N_patches), cl = sqrt(N_patches))
      times$movement[,event_index] <- rexp(length(rates$movement[,event_index]), rates$movement[,event_index])
    }
    return(list(rates = rates, times = times))
  }
  )
}

tauleap <- function(rates, tau) {
  # calculate the number of times each event happens in the interval tau
  events_N <- lapply(rates, \(x) rpois(length(x), x*tau))
  # update states accordingly
  
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
update_state <- function(event_index, event_db) {
  event_type <- event_db$event_types[event_index]
  patch_or_animal_idx <- event_db$event_indices[event_index]
  if (event_type == "growth") {
    h[patch_or_animal_idx] <<- h[patch_or_animal_idx] + 1
  } else if (event_type == "development_l") {
    l[patch_or_animal_idx] <<- l[patch_or_animal_idx]-1
    L[patch_or_animal_idx] <<- L[patch_or_animal_idx]+1
  } else if(event_type == "death_l") {
    l[patch_or_animal_idx] <<- l[patch_or_animal_idx]-1
  } else if(event_type == "death_L") {
    L[patch_or_animal_idx] <<- L[patch_or_animal_idx]-1
  } else if (event_type == "f_decay") {
    f[patch_or_animal_idx] <<- f[patch_or_animal_idx] - 1
  } else if (event_type == "bite") {
    animal_idx    <- patch_or_animal_idx
    patch_idx     <- animal_locations[animal_idx]
    h[patch_idx]  <<- h[patch_idx]  - 1   # reduce patch sward height
    s[animal_idx] <<- s[animal_idx] + 1  # increase stomach content
    l[patch_idx]  <<- l[patch_idx]  - B/h[patch_idx]*l[patch_idx] # reduce number of larve in patch
    L[patch_idx]  <<- L[patch_idx]  - B/h[patch_idx]*L[patch_idx] # reduce number of larve in patch
    a[animal_idx] <<- a[animal_idx] + theta*(B/h[patch_idx])*L[patch_idx] # increase number of larvae in host
    r[animal_idx] <<- r[animal_idx] + psi*B*L[patch_idx]/h[patch_idx] # update host resistance
  } else if(event_type == "death_a") {
    a[patch_or_animal_idx] <<- a[patch_or_animal_idx] - 1
  } else if(event_type == "development_a") {
    a[patch_or_animal_idx] <<- a[patch_or_animal_idx] - 1
    A[patch_or_animal_idx] <<- A[patch_or_animal_idx] + 1
  } else if(event_type == "death_A") {
    A[patch_or_animal_idx] <<- A[patch_or_animal_idx] - 1
  } else if(event_type == "immunity_gain") {
    r[patch_or_animal_idx] <<- r[patch_or_animal_idx] + 1
  } else if(event_type == "immunity_loss") {
    r[patch_or_animal_idx] <<- r[patch_or_animal_idx] - 1
  } else if(event_type == "egg_prod") {
    eg[patch_or_animal_idx] <<- eg[patch_or_animal_idx] + 1
  } else if (event_type == "defecation") {
    animal_idx <- patch_or_animal_idx
    patch_idx <- animal_locations[animal_idx]
    # Heaviside function (Theta(s_k - s0))
    if (s[animal_idx] >= s0) {
      eg[animal_idx] <<- eg[animal_idx] - s0/s[animal_idx]*eg[animal_idx]
      l[patch_idx] <<- l[patch_idx] + s0/s[animal_idx]*eg[animal_idx]
      s[animal_idx] <<- s[animal_idx] - s0
      f[patch_idx] <<- f[patch_idx] + s0
    }
  } else if (event_type == "movement") {
    dest_patch <- event_db$destination[event_index]
    animal_idx <- patch_or_animal_idx
    # current_patch <- animal_locations[animal_idx]
    # c[current_patch] <- c[current_patch]-1
    # c[dest_patch] <- c[dest_patch]+1
    animal_locations[animal_idx] <<- dest_patch
  } 
}

