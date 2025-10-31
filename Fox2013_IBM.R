#### TO DO ####

# Incorporate parasite migration (e.g. Gloworm)
# check functions of resistance: rates lambda, tau, and change theta

#
#
#### --- Main Simulation Function --- ####
#
run_model <- function(seed = NULL, tstep = 30, outf, pars, S, method) {
  # start timer to keep track of how long the sim is taking
  starttime <- Sys.time()
  # set seed if specified
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Set up output dataframes
  current_time <- 0
  time_series <- data.frame(time = numeric(pars$total_time%/%tstep+1), avg_height = 0, avg_a = 0, avg_A = 0, sd_A = 0, avg_l = 0, avg_L = 0, avg_s = 0, avg_f = 0)
  time_series[1,] <- c(current_time, mean(S$h), mean(S$a), mean(S$A), sd(S$A), mean(S$l), mean(S$L), mean(S$s), mean(S$f))
  ix <- 2
  cat(c("time","sim_time", "avg_h", "avg_a", "avg_A", "sd_A", "avg_l", "avg_L","avg_s", "avg_f", "\n"), sep = "\t", file = outf)
  cat(c(0,current_time, mean(S$h), mean(S$a), mean(S$A), sd(S$A), mean(S$l), mean(S$L), mean(S$s), mean(S$f)),"\n", sep = "\t", file = outf, append = T)
  # list with movement info
  movlist <- replicate(pars$Na, data.frame(time = 0, patch = 0), simplify = F)
  # dataframe with parasite dist
  pardb <- as.data.frame(matrix(nrow = pars$total_time%/%tstep+1, ncol = pars$Na+1))
  names(pardb) <- c("time", paste0("A", 1:pars$Na))

  #rates list object
  event_db <- with(pars,
                   list(event_types = c(rep(c("growth","dev_l","death_l",
                                           "death_L","f_decay"),each = N_patches),
                                     rep(c("grazing","death_a",
                                           "dev_a","death_A","immun_gain",
                                           "immun_loss","egg_prod","defecation"), each = Na),
                                     rep("movement",N_patches*Na)),
                     event_indices = c(rep(1:N_patches,5), rep(1:Na,8),rep(1:Na, each = N_patches)),
                     destination = c(rep(NA,N_patches*5+Na*8), rep(1:N_patches, Na)))
  )
  
  # get movement kernel
  movkern <- get_mov_kern(pars$nu, pars$alpha, sqrt(pars$N_patches), sqrt(pars$N_patches))
  
  # calculate initial rates
  event_rates <- get_event_rates0(pars, S, movkern)
  event_times <- get_event_times0(event_rates, pars$Na)
  
  # Main simulation loop
  while (current_time < pars$total_time) {
    if(method == "nrm") {
      # 1. Determine next event (exponential distribution)
      all_times <- unlist(event_times, recursive = T, use.names = F)
      event <- which.min(all_times)
      new_time <- all_times[event]
      
      # 2. Update state variables based on the chosen event
      event_type <- event_db$event_types[event]
      event_index <- event_db$event_indices[event]
      dest <- event_db$destination[event]
      S <- update_state_exact(event_type, event_index, dest, S, pars)
      
      # check for NAs
      nastates <- lapply(S, is.na)|>sapply(any)
      if(any(nastates)) {
        cat("found NA state at ",which(nastates), " at time ", new_time, file = outf, append = T)
        save(S, rates_times, event_type, event_index, dest, prev_rates, file = "errorstate.RData")
        break
      }
      
      # record movement
      if(event_type == 'movement') {
        movlist[[event_index]] <- rbind(movlist[[event_index]], data.frame(time = new_time, patch = dest))
      }
      
      ## 3. Recalculate event rates
      prev_rates <- rates_times$rates
      new_rates <- update_rates_nrm(event_type, event_index, new_time, rates_times, pars, S, movkern)
      event_rates <- new_rates
      
      ## 4. update times
      event_times <- update_times_nrm(event_type, event_index, 
                                            new_rates = new_rates, prev_rates = prev_rates, 
                                            tk = event_times, tm = new_time, dest = dest) 
      # 5. Advance time and record state every tstep minutes
      record_state <- all(floor(new_time)>floor(current_time),floor(new_time)%%tstep==0)
      current_time <- new_time
    } else if(method == "tau") {
      # 1 determine times to next movement for each, and destinations
      movtimes <- apply(rates_times$times$movement, 2, min)
      movdests <- apply(rates_times$times$movement, 2, which.min)
      movorder <- order(movtimes)
      movtimes_ord <- c(movtimes[movorder][1],diff(movtimes[movorder]))
      # register steps
      movlist <- mapply(function(l,t,d) rbind(l, data.frame(time = current_time+t, patch = d)),
                        movlist, movtimes, movdests)

      # update state iteratively, changing the position of the animals
      for (y in 1:pars$Na) {
        ind <- movorder[y]
        dt <- movtimes_ord[y]
        dest <- movdests[ind]
        # update states before movement
        S <- update_state_tau(dt, rates_times, S, pars)
        # update location
        S$animal_locations[ind] <- dest
      }
      # recalculate rates
      rates_times <- get_event_rates0(pars, S, mk)
      # update sim time
      new_time <- current_time+max(movtimes)
      record_state <- all(floor(new_time)>floor(current_time),floor(new_time)%%tstep==0)
      current_time <- new_time
    }
    if(record_state) {
      time_series[ix,] <- c(new_time, mean(S$h), mean(S$a), mean(S$A), sd(S$A), mean(S$l), mean(S$L), mean(S$s), mean(S$f))
      elapsedtime <- as.numeric(difftime(Sys.time(), starttime), units = "hours")
      cat(c(elapsedtime, new_time, mean(S$h), mean(S$a), mean(S$A), sd(S$A), mean(S$l), mean(S$L), mean(S$s), mean(S$f), "\n"), sep = "\t", file = outf, append = T)
      pardb[ix,] <- c(new_time, S$A)
      ix <- ix+1
    }
  }
  return(list(time_series, movlist, pardb))
}

#
#### --- Helper Functions --- ####
#
get_event_rates0 <- function(pars, S, mk) {
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
    grazing <- beta * (h[animal_locations] - h0) * exp(-mu_f * f[animal_locations] *(a+A)*Lambda)
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
    defecation <- sapply(s, \(x) ifelse(x>s0,f_dep*(x-s0),0))
    # Movement rates for each animal
    movement <- t(mk[animal_locations,]*h)
    
    #create output lists
    rates <- list(growth = growth, dev_l = dev_l, death_l = death_l, death_L = death_L, f_decay = f_decay, 
                grazing = grazing, 
                death_a = death_a, dev_a = dev_a, death_A = death_A, immun_gain = immun_gain, immun_loss = immun_loss, egg_prod = egg_prod, 
                defecation = defecation, 
                movement = movement)
    
    # replace negatives for 0s
    rates <- lapply(rates, \(x) replace(x, x<0, 0))
    
    return(rates)
  }
  )
}

## initial event times, for the NRM method
get_event_times0 <- function(r, Na) {
  times <- lapply(r, \(x) 1/x*log(1/runif(length(x))))
  # replace 0s
  times <- lapply(times, \(x) replace(x, x<0, 0))
  # substitute NAs and NaNs with inf
  times$movement <- matrix(times$movement, ncol = Na)
  return(times)
}


# function to get the event rates using the Optimized method (sorting method)
# that does not recalculate every rate, but rather only those affected by the
# latest event
update_rates_nrm <- function(event_type, event_index, tau_mu, rates_times, pars, S, mk) {
  with(c(pars, S, rates_times), {
    # only update based on previous event
    if(event_type == "growth") {
      rates$growth[event_index] <- gamma * h[event_index] * (1 - h[event_index] / h_max)
      # check which animals are in the patch
      animals_in_patch <- animal_locations==event_index
      if(any(animals_in_patch)) {
        rates$grazing[animals_in_patch] <- beta * (h[event_index] - h0) * exp(-mu_f * f[event_index] *(a[animals_in_patch]+A[animals_in_patch])*Lambda)
      }
    } else if (event_type=="dev_l") {
      rates$dev_l[event_index]   <- epsilon * l[event_index]
      rates$death_l[event_index] <- omega * l[event_index]
      rates$death_L[event_index] <- rho*L[event_index]
    } else if (event_type=="death_l") {
      rates$death_l[event_index] <- omega * l[event_index]
      rates$dev_l[event_index]   <- epsilon * l[event_index]
    } else if (event_type=="death_L") {
      rates$death_L[event_index] <- rho * L[event_index]
    } else if (event_type=="f_decay") {
      rates$f_decay[event_index] <- phi * f[event_index]
      animals_in_patch <- animal_locations==event_index
      if(any(animals_in_patch)) {
        rates$grazing[animals_in_patch] <- beta * (h[event_index] - h0) * exp(-mu_f * f[event_index] *(a[animals_in_patch]+A[animals_in_patch])*Lambda)
      }
    } else if (event_type=="grazing") {
      patch <- animal_locations[event_index]
      rates$grazing[event_index] <- beta * (h[patch] - h0) * exp(-mu_f * f[patch] *(a[event_index]+A[event_index])*Lambda)
      rates$growth[patch]  <- gamma * h[patch] * (1 - h[patch] / h_max)
      rates$death_a[event_index] <- zeta * a[event_index]
      rates$dev_a[event_index]   <- chi * a[event_index]
      rates$dev_l[patch]   <- epsilon * l[patch]
      rates$death_l[patch] <- omega * l[patch]
      rates$death_L[patch] <- rho * L[patch]
      rates$defecation[event_index] <- ifelse(s[event_index]>s0, f_dep*(s[event_index]-s0),0)
    } else if (event_type=="death_a") {
      rates$death_a[event_index] <- zeta * a[event_index]
      rates$dev_a[event_index]   <- chi * a[event_index]
      rates$immun_gain[event_index] <- (a[event_index] + A[event_index]) * eta
    } else if (event_type=="dev_a") {
      rates$dev_a[event_index]    <- chi * a[event_index]
      rates$death_a[event_index]  <- zeta * a[event_index]
      rates$death_A[event_index]  <- tau * A[event_index]
      rates$egg_prod[event_index] <- lambda * A[event_index]/2
      rates$immun_gain[event_index] <- (a[event_index] + A[event_index]) * eta
    } else if (event_type=="death_A") {
      rates$death_A[event_index]   <- tau * A[event_index]
      rates$egg_prod[event_index]  <- lambda * A[event_index]/2
      rates$immun_gain[event_index] <- (a[event_index] + A[event_index]) * eta
    } else if (event_type=="immun_loss") {
      rates$immun_loss[event_index] <- sig * r[event_index]
    } else if (event_type=="defecation") {
      rates$defecation[event_index] <- ifelse(s[event_index]>s0, f_dep*(s[event_index]-s0), 0)
      patch <- animal_locations[event_index]
      rates$grazing[event_index] <- beta * (h[patch] - h0) * exp(-mu_f * f[patch] *(a[event_index]+A[event_index])*Lambda)
      rates$dev_l[patch]   <- epsilon * l[patch]
      rates$death_l[patch] <- omega * l[patch]
      rates$f_decay[patch] <- phi * f[patch]
    } else if (event_type=="movement") {
      rates$movement[,event_index] <- mk[animal_locations[event_index],]*h
    }
    return(rates)
  }
  )
}

# function to update the absolute time to the next event, for every possible
# event
update_times_nrm <- function(event_type, event_index, new_rates, prev_rates, tk, tm, dest) {
  new_times <- tk
  changed_rates <- mapply("!=",prev_rates, new_rates)|>sapply(any)
  eix <- which(event_type == names(new_times))
  kevents <- setdiff(which(changed_rates),eix)
  if(length(kevents)>0) {
    new_times[kevents] <- mapply(function(p, n, t) {
      A <- t
      # substitute times of Inf (rate = 0) with a random number 
      i <- which(p == 0 & n != p)
      A[i] <- 1/n[i]*log(1/runif(length(i)))+tm
      # substitute non Inf times using the old rate and new rates
      i <- which(p > 0 & n != p)
      A[i] <- p[i]/n[i]*(t[i]-tm)+tm
      return(list(A))
    }, prev_rates[kevents], new_rates[kevents], tk[kevents])
  }
  if(event_type=='movement') {
    new_times$movement[, event_index] <- 1/new_rates$movement[,event_index]*log(1/runif(nrow(new_times$movement)))+tm
  } else {
    new_times[[eix]][event_index] <- 1/new_rates[[eix]][event_index]*log(1/runif(1))+tm
  }
  return(new_times)
}

# function to calculate the rates at which an individual in patch i moves to
# patch j
mov_rate <- function(cell_curr, h, mk) {
  mkcoef <- mk[cell_curr,]
  rate <- mkcoef*h
  return(rate)
}

# function to update state space based on event type. 
update_state_exact <- function(event_type, event_index, dest = NULL, S, pars) {
  with(c(pars, S),{
    if (event_type == "growth") {
      h[event_index] <- h[event_index] + 1
    } else if (event_type == "dev_l") {
      l[event_index] <- max(l[event_index]-1, 0)
      L[event_index] <- L[event_index]+1
    } else if(event_type == "death_l") {
      l[event_index] <- max(l[event_index]-1, 0)
    } else if(event_type == "death_L") {
      L[event_index] <- max(L[event_index]-1, 0)
    } else if (event_type == "f_decay") {
      f[event_index] <- max(f[event_index] - 1, 0)
    } else if (event_type == "grazing") {
      animal    <- event_index
      patch     <- animal_locations[animal]
      h[patch]  <- h[patch]  - 1   # reduce patch sward height
      s[animal] <- s[animal] + 1  # increase stomach content
      l[patch]  <- max(l[patch]  - B/h[patch]*l[patch], 0) # reduce number of larve in patch
      L[patch]  <- max(L[patch]  - B/h[patch]*L[patch], 0) # reduce number of larve in patch
      a[animal] <- a[animal] + theta*(B/h[patch])*L[patch] # increase number of larvae in host
      r[animal] <- r[animal] + psi*B*L[patch]/h[patch] # update host resistance
    } else if(event_type == "death_a") {
      a[event_index] <- max(a[event_index] - 1, 0)
    } else if(event_type == "dev_a") {
      a[event_index] <- max(a[event_index] - 1, 0)
      A[event_index] <- A[event_index] + 1
    } else if(event_type == "death_A") {
      A[event_index] <- max(A[event_index] - 1, 0)
    } else if(event_type == "immun_gain") {
      r[event_index] <- r[event_index] + 1
    } else if(event_type == "immun_loss") {
      r[event_index] <- max(r[event_index] - 1, 0)
    } else if(event_type == "egg_prod") {
      eg[event_index] <- eg[event_index] + 1
    } else if (event_type == "defecation") {
      animal <- event_index
      patch <- animal_locations[animal]
      # Heaviside function (Theta(s_k - s0))
      if (s[animal] >= s0) {
        eg[animal] <- eg[animal] - s0/s[animal]*eg[animal]
        l[patch]   <- l[patch] + s0/s[animal]*eg[animal]
        s[animal]  <- s[animal] - s0
        f[patch]   <- f[patch] + s0
      }
    } else if (event_type == "movement") {
      animal_locations[event_index] <- dest
    } 
    return(list(h = h, f = f, l = l, L = L, 
                animal_locations = animal_locations,
                r = r, a = a, A = A,
                eg = eg, s = s))
  }
  )
}

# tau leap function estimates the number of (non-movement) events that have
# occurred in a given timespan dt and returns the updated state. The number of
# events is drawn from a Poisson distribution
update_state_tau <- function(dt, rates, S, pars) {
  rates <- rates[1:13]
  # calculate the number of times each event happens in the interval dt
  events_N <- lapply(rates, \(x) rpois(length(x), x*dt))
  occurring <- sapply(events_N,\(x) sum(x)>0)
  patches <- S$animal_locations
  with(pars, {
    # growth
    if(occurring[1]) S$h <- S$h+events_N$growth
    # dev_l
    if(occurring[2]) {
      S$l <- S$l-events_N$dev_l
      S$L <- S$L+events_N$dev_l
    }
    # death_l
    if(occurring[3]) S$l <- S$l-events_N$death_l
    # death_L
    if(occurring[4]) S$L <- S$L-events_N$death_L
    # f_decay
    if(occurring[5]) S$f <- S$f-events_N$f_decay
    # grazing
    if(occurring[6]) {
      S$h[patches] <- S$h[patches]-events_N$grazing
      S$s <- S$s + events_N$grazing
      S$l[patches] <- S$l[patches]-B/S$h[patches]*S$l[patches]*events_N$grazing
      S$L[patches] <- S$L[patches]-B/S$h[patches]*S$L[patches]*events_N$grazing
      S$a <- S$a+theta*(B/S$h[patches])*S$L[patches]*events_N$grazing
      S$r <- S$r+psi*B/S$h[patches]*S$L[patches]*events_N$grazing
    }
    # death_a
    if(occurring[7]) S$a <- S$a-events_N$death_a
    # dev_a
    if(occurring[8]) {
      S$a <- S$a-events_N$dev_a
      S$A <- S$A+events_N$dev_a
    }
    # death_A
    if(occurring[9]) S$A <- S$A - events_N$death_A
    # immun_gain
    if(occurring[10])   S$r <- S$r+events_N$immun_gain
    # immun_loss
    if(occurring[11])   S$r <- S$r-events_N$immun_loss
    # egg_prod
    if(occurring[12])   S$eg <- S$eg+events_N$egg_prod
    # defecation
    if(occurring[13]) {
      def_inds <- which(S$s>s0)
      if(length(def_inds)>0) {
        S$eg[def_inds] <- S$eg[def_inds] - s0/S$s[def_inds]*S$eg[def_inds]*events_N$defecation[def_inds]
        S$l[patches][def_inds] <- S$l[patches][def_inds]+s0/S$s[def_inds]*S$eg[def_inds]*events_N$defecation[def_inds]
        S$s[def_inds] <- S$s[def_inds]-s0*events_N$defecation[def_inds]
        S$f[patches][def_inds] <- S$f[patches][def_inds]+s0*events_N$defecation[def_inds]
      }
    }
    # replace neg values
    S <- lapply(S, \(x) replace(x, list = x<0, values = 0))
    return(list(S, events_N))
  }
  )
}

# function to create the movement kernel, a N_patch x N_patch matrix that has
# the coefficient considering Euclidean distance
get_mov_kern <- function(nu, alpha, rw, cl) {
  dists <- as.matrix(dist(expand.grid(1:rw, 1:cl)))
  # power law using the Euclidean distance, and the alpha value
  F_i <- dists^-alpha
  # set diagonal to 0
  diag(F_i) <- 0
  # standardizing sum
  z_i <- rowSums(F_i)
  mkcoefs <- nu/z_i*F_i
  return(mkcoefs)
}
