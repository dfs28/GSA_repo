#### Genome sequence analysis assignment

#### 1) Write a program which implements a HMM ----
#reading the model and parameters from a file. 
#Use it to simulate 115 emitted values from the following model:

#Set up initial parameters

set.seed(15)

S <- c(0, 1) #Hidden state space

V <- c(1, 2, 3, 4, 5) #Emission space

A <- matrix(c(0.8, 0.1, 0.2, 0.9), nrow = 2) #Transition matrix

Mu_0 <- c(0.5, 0.5) #Initial distribution

B <- matrix(c(0.2, 0, 0.5, 0.1, 0.2, 0.4, 0.1, 0.4, 0, 0.1), nrow = 2) #Emission matrix

#Make function to emit chain

hidden_markov <- function(S, V, A, Mu_0, B, clen) {
  #Function to take hidden markov chain and simulate emission states
  
  #Set up vectors to hold output
  states <- data.frame(chain = rep(0, clen), emission = rep(0, clen))
  
  rownames(A) <- rownames(B) <- S
  
  #Set initial hidden and emission states given initial distribution and emission matrix
  states$chain[1] <- sample(S, 1, prob = Mu_0)
  states$emission[1] <- sample(V, 1, prob = B[as.character(states$chain[1]),])
  
  #Work through states and emit based on transmission probabilities
  for (i in 2:clen) {
    states$chain[i] <- sample(S, 1, prob = A[as.character(states$chain[i - 1]), ])
    states$emission[i] <- sample(V, 1, prob = B[as.character(states$chain[i]),])
  }
  
  return(states)
}

states <- hidden_markov(S, V, A, Mu_0, B, 115)

require(ggplot2)
ggplot(states) + geom_col(aes(x = 1:115, y = emission, fill = as.factor(chain))) +
  theme(legend.title = 'Hidden State') +
  theme_bw() +
  xlab('Position') +
  ylab('Emission State')

hidden.states <- states$chain
states <- states$emission

#### 2) Implement the forward algorithm and calculate the log likelihood given previous data ----

likelihood_forward <- function(states, S, V, A, Mu_0, B) {
  #Function to implement the forward algorithm
  
  #Set up dataframe
  alphas <- matrix(ncol = length(Mu_0), nrow = length(states))
  
  #Set up initialisation states
  alphas[1, ] <- Mu_0*B[, states[1]]
  
  #Set up recursion relationship
  for (n in 2:length(states)){
    
    #Multiply the priors by the transition matrix
    sum_alpha_prior <- alphas[n-1, ] %*% A
    
    #Multiply the hidden state probabilities by the emission probabilities
    alphas[n, ] = sum_alpha_prior * B[, states[n]]
  }
  return(alphas)
}

likelihood_tab <- likelihood_forward(states, S, V, A, Mu_0, B)
sum(likelihood_tab[dim(likelihood_tab)[1],])

#### 3) Calculate GC content in 100bp windows ----

S_cerevisiae_III <- readLines("~/Documents/Masters/Course materials/Genomic sequence analysis/S_cereviciae_III.fasta")
S_cerevisiae_III <- S_cerevisiae_III[-c(1:3)]
S_cerevisiae_III <- paste(S_cerevisiae_III, collapse = '')
S_cerevisiae_III <- strsplit(S_cerevisiae_III, split = '')[[1]]

GC_content <- function(input, loc) {
  #Function to find GC content in any 100bp window
  
  #Get window
  window.vec <- input[seq(loc, loc + 99)]
  #Return GC content as percentage
  length(which(window.vec %in% c('G', 'C')))/length(which(is.na(window.vec) == FALSE))
}

#Apply sliding window over sequence
#GC_window_content <- sapply(seq(1, length(S_cerevisiae_III) - 100), GC_content, input = S_cerevisiae_III)

#Try without sliding window
GC_window_content <- sapply(seq(1, length(S_cerevisiae_III), by = 100), GC_content, input = S_cerevisiae_III)

#Get emission state probabilities for HMM at beginning
bin_prop <- sapply(sort(unique(states)), function(x, y) (length(which(y == x))), y = states)
names(bin_prop) <- sort(unique(states))
bin_prop <- bin_prop/length(states)

#Get GC contents corresponding to each bin
bin_cumulative <- sapply(1:length(bin_prop), function(x, y) (sum(y[1:x])), y = bin_prop)
gc_bins <- c()
for (i in 1:length(bin_cumulative)) {
  gc_bins[i] <- max(sort(GC_window_content)[seq(1, floor(bin_cumulative[i]*length(GC_window_content)))])
}

#Calculate emission sequence under this model
GC_emission_seq <- rep(NA, length(GC_window_content))
GC_emission_seq[which(GC_window_content <= gc_bins[1])] <- 1
for (i in 2:length(gc_bins)) {
  GC_emission_seq[which(GC_window_content <= gc_bins[i] & GC_window_content > gc_bins[i-1])] <- i
}

#Calculate the log likelihood of this sequence
log_likelihood_forward <- function(states, S, V, A, Mu_0, B) {
  #Function to implement the forward algorithm
  
  #Set up dataframe - first col is c_n, rest are alpha_hat
  c_n <- vector(length = length(states))
  alpha_hats <- matrix(ncol = length(Mu_0), nrow = length(states))
  
  #Set up initialisation states
  c_n[1] <- sum(Mu_0*B[, states[1]])
  alpha_hats[1,] <- Mu_0*B[, states[1]]/c_n[1]
  
  #Set up recursion relationship
  for (n in 2:length(states)){
    
    #Calculate c_n
    c_n[n] <- sum((alpha_hats[n-1, ] %*% A) * B[, states[n]])
    
    #Calculate alpha_hat[n] (which is essentially alpha normalised by the sum of the alphas - c_n)
    alpha_hats[n, ] <- ((alpha_hats[n-1, ] %*% A) * B[, states[n]])/c_n[n]

  }
  
  #Completion - sum the logs
  return(list(c_n, alpha_hats))

}

sum(log(log_likelihood_forward(GC_emission_seq, S, V, A, Mu_0, B)[[1]]))

#### 4) Extend to implement Baum-Welch algorithm ----

scaled_backward <- function(states, S, V, A, Mu_0, B) {
  #First build function to calculate scaled backward variable
  
  #Set up dataframe
  beta_hats <- matrix(ncol = length(Mu_0), nrow = length(states))
  
  #Run scaled forward algorithm
  forward <- log_likelihood_forward(states, S, V, A, Mu_0, B)
  c_n <- forward[[1]]
  alpha_hats <- forward[[2]]
  
  #Set up initialisation states
  beta_hats[length(states),] <- 1/c_n[length(c_n)]
  
  #Set up recursion relationship
  for (n in (length(states)-1):1){
    
    #Calculate beta_hats[n]
    beta_hats[n, ] <- ((beta_hats[n+1, ] * B[, states[n+1]]) %*% t(A))/c_n[n]
  }
  
  return(list(c_n, alpha_hats, beta_hats))
  
}

backward <- scaled_backward(states, S, V, A, Mu_0, B)

#Calculate expected values for theta

baum_welch <- function(states, S, V, A, Mu_0, B, threshold, iterations) {
  #Function to calculate the expected values of theta
  
  #First get backward and forward variables
  alpha_beta <- scaled_backward(states, S, V, A, Mu_0, B)
  c_n <- alpha_beta[[1]]
  alpha_hats <- alpha_beta[[2]]
  beta_hats <- alpha_beta[[3]]
  
  #Calculate A prime - first caculate E(n_ij)
  exp_n_ij <- function(alpha_hats, beta_hats, states, A, B, i, j) {
    
    #Define function to apply across  all locations
    func <- function(alpha_hats, A, beta_hats, B, n, i, j) {
      alpha_hats[n, i]*A[i,j]*B[j, states[n+1]]*beta_hats[n+1, j]
    }
    
    #Apply function to all 0 to N-1
    sum(sapply(1:(length(states)-1), func, alpha_hats = alpha_hats, 
               beta_hats = beta_hats, A = A, B = B, i = i, j = j))
    
  }
  
  #Function to make A prime
  A_prime_func <- function(alpha_hats, beta_hats, states, A, B, Mu_0) {
    
    #Make A prime matrix
    A_prime <- matrix(nrow = length(Mu_0), ncol = length(Mu_0))
    
    #Fill A prime matrix
    for (i in 1:length(Mu_0)) {
      for (j in 1:length(Mu_0)) {
        A_prime[i,j] <- exp_n_ij(alpha_hats,  beta_hats, states, A, B, i, j)
      }
    }
    
    #Normalise A prime matrix by row sums and return it - this poss should be col? - previously tried with inverse matrix
    t(apply(A_prime, 1, function(x) (x/sum(x))))
  }
  
  #Next calculate Mu_prime - *****currently not correct/ outputting incorrectly*****
  Mu_prime_func <- function(alpha_hats, beta_hats, states, A, B) {
    
    #Mu_prime <- 
      alpha_hats[1,] * ((B[, states[2]] * beta_hats[2,]) %*% t(A))
    #Mu_prime/sum(Mu_prime)
  }
  
  #Next calculate B_prime
  B_prime_func <- function(alpha_hats, beta_hats, states, A, B) {
    
    #Make B prime matrix
    B_prime <- B
    
    #Function to do non-summed numerator
    B_func <- function(alpha_hats, A, beta_hats, B, n, i, j, k, states) {
      if (states[n] == k) {
        return(alpha_hats[n, i] * A[i, j] * B[j, states[n + 1]] * beta_hats[n + 1, j]) 
      } else {return(0)}
    }
    
    
    #Fill B prime matrix - working through the matrix by row then column
    for (trans_i in 1:dim(B)[1]) {
      for (trans_k in 1:dim(B)[2]) {
        B_prime_ik <- c()
        for(trans_j in 1:dim(B)[1]) {
          
          #Calculate numerator using function before, this time summing
          B_prime_ik[trans_j] <- sum(sapply(1:(length(states)-1), B_func, alpha_hats = alpha_hats, 
                                            beta_hats = beta_hats, A = A, B = B, i = trans_i, j = trans_j, states = states, k = trans_k))
        }
        
        #Normalise by the summed expected value of E(n_ij)
        B_prime[trans_i,trans_k] <- sum(B_prime_ik)/sum(sapply(1:2, exp_n_ij, i = trans_i, 
                                                               alpha_hats = alpha_hats,  beta_hats = beta_hats, 
                                                               states = states, A = A, B = B))
                                                               
      }
      
    }
    
    return(B_prime)
  }
  
  #Now iterate over t until threshold < value set above
  L_t <- sum(log(log_likelihood_forward(states, S, V, A, Mu_0, B)[[1]]))
  iterator <- 1
  pb <- txtProgressBar(min = 0, max = iterations, style = 3)
  delta <- 10000
  
  while (delta > threshold & iterator < iterations) {
    
    #First get backward and forward variables
    alpha_beta <- scaled_backward(states, S, V, A, Mu_0, B)
    c_n <- alpha_beta[[1]]
    alpha_hats <- alpha_beta[[2]]
    beta_hats <- alpha_beta[[3]]
    
    #Get updated parameters
    B_prime <- B_prime_func(alpha_hats, beta_hats, states, A, B)
    
    A_prime <- A_prime_func(alpha_hats, beta_hats, states, A, B, Mu_0)
    
    Mu_prime <- Mu_prime_func(alpha_hats, beta_hats, states, A, B)
    
    L_next <- sum(log(log_likelihood_forward(states, S, V, A_prime, Mu_prime, B_prime)[[1]]))
    
    delta <- abs(L_t - L_next)
    delta
    
    #Now update the parameters (theta_t -> theta_t+1)
    B <- B_prime
    A <- A_prime
    Mu_0 <- Mu_prime
    L_t <- L_next
    
    iterator <- iterator + 1
    setTxtProgressBar(pb, iterator)
  }
  
  close(pb)
  
  return(list(A, B, Mu_0, L_t))
}

updated_theta <- baum_welch(GC_emission_seq, S, V, A, Mu_0, B, 1e-5, 10000)
updated_A <- updated_theta[[1]]
updated_B <- updated_theta[[2]]
updated_Mu <- updated_theta[[3]]

## Now implement Viterbi sequence

viterbi_inference <- function(states, Mu_0, S, V, B, A) {
  
  #Set up initial parameters
  N <- length(states)
  phi <- psi <- matrix(nrow = N, ncol = length(Mu_0))
  phi[1,] <- log(Mu_0) + log(B[,states[1]])
  
  #Run recursion relationship
  for (n in 2:N) {
    
    #Recursion of phi
    to.max <- log(A) + phi[n-1,]
    
    #Do maximisation and add to other term
    phi[n,] <- log(B[,states[n]]) + apply(to.max, 2, max)
    
    #Recursion of psi - do argmax
    psi[n,] <- apply(to.max,  2, function(x) (which(x == max(x))))
  }
  
  ##Do traceback
  i_star <- rep(NA, times = N)
  
  #Get Nth value
  i_star[N] <- which(phi[N,] == max(phi[N,]))
  
  #Trace back along vector
  for (n in (N-1):1) {
    i_star[n] <- psi[n+1, i_star[n+1]]
  }
  
  return(i_star)
}

inferred.states <- viterbi_inference(GC_emission_seq, updated_Mu, S, V, updated_B, updated_A)
ggplot() + geom_col(aes(x = 1:length(GC_emission_seq), y = GC_emission_seq, fill = as.factor(inferred.states))) +
  theme(legend.title = 'Hidden State') +
  theme_minimal() +
  xlab('Position') +
  ylab('Emission State')

#Now try to build the continuous emission distribution


