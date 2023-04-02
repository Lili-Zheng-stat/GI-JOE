##GI-JOE (edge-wise) functions 
#----------------------------------------------------#
library(MTS)
library(rTensor)
#----------------------------------------------------#
edge_testing <- function(X, n_total, Ind, Entryest_Sigma, p, N, node_a, node_b, ...){
  ## Perform GI-JOE (edge-wise) for pairwise measured data
  ## Pairwise measured data typically has huge number of rows, and hence general method for computing the quadruple sample sizes and even indexing can be inefficient.
  ## Therefore, we wrote a separate function for GI-JOE (edge-wise) with pairwise measured data for computational reasons.
  
  ## X: n_total by p matrix; Each row of X represents a sample, with two variables being observed, while the rest are all zero
  ## Ind: list of (p-1) list, Ind[[i]][[j]] tells the row indices of X where pair i,j are observed, for i<j
  ## Entryest_Sigma: entry-wise unbiased estimate of the covariance matrix
  ## N: p by p pairwise sample size matrix
  ## node_a, node_b: tested node pair
  ## Additional input: tuning parameter constant c; if not provided, chosen by stability selection.
  ## return: the test statistic for Theta[node_a,node_b]/Theta[node_a,node_a].
  input <- list(...);
  
  ## find PSD projection of Sigma_hat
  proj_out <- proj_cov(Entryest_Sigma, N, 0.01, 1, matrix(rep(0, p^2), p, p), matrix(rep(1, p^2), p, p), 0.001, 500)
  PSD_Sigma  <- proj_out$projected_S
  
  ## neighborhood lasso for a and debiasing vector for (a,b)
  #optimization parameters for neighborhood regression: 
  eta <- 1;tol <- 0.00001
  #eta is the initial step size and tol is the tolerance in the projected gradient descent algorithm
  lambda_scale <- sqrt(log(p)/apply(N, 1, min))
  #lambda_scale is the scaling vector of the tuning parameter vector
  beta_hat <- rep(0, p)
  theta_hat <- rep(0, p)
  #initialize neighborhood regression estimator for node_a and node_b
  ind <- (1 : p)[-node_a]
  node_b_ind <- which(ind == node_b)
  if(length(input) > 0){
    tuning_c <- input[[1]]
    #set tuning constant as the additional input
    if(length(tuning_c > 1)){
      fit_out1 <- fit_nb(PSD_Sigma, node_a, lambda_scale[-node_a] * tuning_c[1], eta, tol)
      tuning_c1 <- tuning_c[1]
      #neighborhood regression for node_a upon [p]\node_a
      fit_out2 <- fit_nb(PSD_Sigma[-node_a,-node_a], node_b_ind, lambda_scale[-c(node_a,node_b)] * tuning_c[2], eta, tol)
      tuning_c2 <- tuning_c[2]
      #neighborhood regression for node_b upon [p]\{node_a, node_b}
    }else{
      fit_out1 <- fit_nb(PSD_Sigma, node_a, lambda_scale[-node_a] * tuning_c, eta, tol)
      tuning_c1 <- tuning_c
      fit_out2 <- fit_nb(PSD_Sigma[-node_a,-node_a], node_b_ind, lambda_scale[-c(node_a,node_b)] * tuning_c, eta, tol)
      tuning_c2 <- tuning_c
    }
  }else{
    #neighborhood regression with stability tuning when no additional input is provided
    n_subsample <- 20; stb_thrs <- 0.05
    fit_out1 <- est_nb_tuning(X, n_total, p, Ind, n_subsample, stb_thrs, PSD_Sigma, node_a, lambda_scale[-node_a], eta, tol)
    tuning_c1 <- fit_out1$tuning_c
    Ind_del <- find_ind_delete_node(Ind, p, node_a); 
    fit_out2 <- est_nb_tuning(X[, -node_a], n_total, p - 1, Ind_del, n_subsample, stb_thrs, PSD_Sigma[-node_a, -node_a], node_b_ind, lambda_scale[-c(node_a, node_b)], eta, tol)
    tuning_c2 <- fit_out2$tuning_c
  }
  beta_hat[-node_a] <- fit_out1$beta_hat
  #estimated -Theta[,node_a]/Theta[node_a,node_a], but with beta_hat[node_a] = 0
  theta_hat[-c(node_a,node_b)] <- - fit_out2$beta_hat
  theta_hat[node_b] <- 1
  tau_varest <- as.numeric(1 / (t(theta_hat) %*% PSD_Sigma %*% theta_hat))
  tau_debiasing <- as.numeric(1 / (PSD_Sigma[node_b, ] %*% theta_hat))
  theta_hat_debiasing <- theta_hat * tau_debiasing
  theta_hat_varest <- theta_hat * tau_varest
  #estimated Theta^{(a)}_{,node_b} defined in the paper
  
  #debiased neighborhood lasso (node_a, node_b)
  beta_tilde = as.numeric(beta_hat[node_b] - t(theta_hat_debiasing) %*% (Entryest_Sigma %*% beta_hat - Entryest_Sigma[, node_a]))
  #variance estimation
  beta_hat2 <- -beta_hat;
  beta_hat2[node_a] <- 1;
  #estimated Theta[,node_a]/Theta[node_a,node_a]
  var_est <- find_var(PSD_Sigma, N, beta_hat2, theta_hat_varest)
  
  #test statistic
  test <- beta_tilde/sqrt(var_est)
  return(list(beta_hat = beta_hat, theta_hat_debiasing = theta_hat_debiasing, theta_hat_varest = theta_hat_varest, beta_tilde = beta_tilde, tuning_c1 = tuning_c1, tuning_c2 = tuning_c2, var_est = var_est, test = test, PSD_Sigma = PSD_Sigma))
}

est_nb_tuning <- function(X, n_total, p, Ind, n_subsample, stb_thrs, Sigma, a, lambda_scale, eta, tol){
  ## estimate neighborhood with tuning parameter selected based on stability
  #find range of tuning_c
  tuning_c_max <- 1
  out <- fit_nb_zerosamplesz(Sigma, a, lambda_scale * tuning_c_max, eta, tol)
  while(sum(out$beta_hat != 0) > 0){
    tuning_c_max <- tuning_c_max * 2
    out <- fit_nb_zerosamplesz(Sigma, a, lambda_scale * tuning_c_max, eta, tol)
  }
  if(tuning_c_max == 1){
    while(sum(out$beta_hat != 0) == 0){
      tuning_c_max <- tuning_c_max / 2
      out <- fit_nb_zerosamplesz(Sigma, a, lambda_scale * tuning_c_max, eta, tol)
    }
    tuning_c_max <- tuning_c_max * 2
  }
  tuning_c_vec <- exp(log(tuning_c_max)-log(10) / 19 * (0 : 19))
  tuning_out <- stb_tuning(X, n_total, p, Ind, n_subsample, a, tuning_c_vec, eta, tol)
  if(tuning_out$stb_vec[1] > stb_thrs){
    choice_ind <- 1
  }else if(max(tuning_out$stb_vec) <= stb_thrs){
    choice_ind <- length(tuning_c_vec)
  }else{
    choice_ind <- min(which(tuning_out$stb_vec > stb_thrs)) - 1
  }

  out_final <- fit_nb(Sigma, a, lambda_scale * tuning_c_vec[choice_ind], eta, tol)
  out_final$stb_vec <- tuning_out$stb_vec
  out_final$theta_stb <- tuning_out$theta_stb
  out_final$tuning_c_vec <- tuning_c_vec
  out_final$tuning_c <- tuning_c_vec[choice_ind]
  return(out_final)
}

stb_tuning <- function(X, n_total, p, Ind, n_subsample, a, tuning_c_vec, eta, tol){
  #use stability tuning; n_subsample is the number of subsamples taken; prob_subsample is the probability of each sample being included in one subsample;
  nbhd_est <- array(rep(0, length(tuning_c_vec) * n_subsample * (p - 1)), dim = c(length(tuning_c_vec), n_subsample, p - 1))
  for(kk in 1 : n_subsample){
    out_subsample <- calc_Sigma_subsample(X, p, Ind, kk)
    Sigma_subsample <- out_subsample$Sigma
    N_subsample <- out_subsample$N
    lambda_scale_tmp <- sqrt(log(p)/apply(N_subsample,1,min))
    #find PSD projection of Sigma_hat
    PSD_Sigma_subsample  <- proj_cov(Sigma_subsample, N_subsample, 0.01, 1, matrix(rep(0, p^2), p, p), matrix(rep(1, p^2), p, p), 0.001, 500)$projected_S
    for(j in 1 : length(tuning_c_vec)){
      out <- fit_nb_zerosamplesz(PSD_Sigma_subsample, a, lambda_scale_tmp[-a] * tuning_c_vec[j], eta, tol)
      nbhd_est[j, kk, ] <- out$beta_hat != 0
    }
  }
  theta_stb <- apply(nbhd_est, c(1, 3), mean)
  stb_vec <- apply(2 * theta_stb * (1 - theta_stb), 1, mean)
  return(list(theta_stb = theta_stb, stb_vec = stb_vec))
}



calc_Sigma_subsample <- function(X, p, Ind, seed){
  #calculate covariate estimate (entrywise) using data X with subsample index selec_ind; Ind tells the sample indices for each pair (i,j) with i<j in X;
  #For any pair, at least one sample is available
  Sigma <- matrix(rep(0, p^2), p, p)
  dSigma <- rep(0, p)
  N_subsample <- matrix(rep(0, p^2), p, p)
  for(i in 1 : (p - 1)){
    for(j in (i + 1) : p){
      Ind_tmp <- Ind[[i]][[j]]
      n_tmp <- length(Ind_tmp)
      if(n_tmp > 144){
        prob_tmp <- 10 / sqrt(n_tmp)
      }else{
        prob_tmp <- 0.8
      }
      set.seed(seed)
      select_ind_tmp <- as.logical(rbinom(n_tmp, 1, prob_tmp))
      Ind_subsample <- Ind_tmp[select_ind_tmp]
      N_subsample[i, j] <- length(Ind_subsample)
      if(N_subsample[i,j] == 0){
        Sigma[i, j] <- 0
      }else{
        Sigma[i, j] <- mean(X[Ind_subsample, i] * X[Ind_subsample, j])
      }
      dSigma[i] <- dSigma[i] + sum(X[Ind_subsample, i]^2)
      dSigma[j] <- dSigma[j] + sum(X[Ind_subsample, j]^2)
    }
  }
  N_subsample <- N_subsample + t(N_subsample)
  diag(N_subsample) <- apply(N_subsample, 1, sum)
  Sigma <- Sigma + t(Sigma)
  dSigma[is.nan(dSigma)] <- 1
  diag(Sigma) <- dSigma / diag(N_subsample)
  return(list(Sigma = Sigma, N = N_subsample))
}



find_var <- function(Sigma, N, beta, theta){
  p <- dim(Sigma)[[1]]
  var <- 0
  supp_theta <- which(theta != 0);
  supp_beta <- which(beta != 0);
  for(j1 in supp_theta){
    for(j2 in supp_beta){
      if(j1 == j2){
        for(j3 in supp_theta){
          if(j3 == j1){
            for(j4 in supp_beta){
              var <- var + (Sigma[j1, j3] * Sigma[j2, j4] + Sigma[j1, j4] * Sigma[j2, j3]) * 
                N[j1, j4] / N[j1, j2] / N[j3, j4] * theta[j1] * theta[j3] * beta[j2] * beta[j4]
            }
          }else{
            for(j4 in c(j1, j3)){
              var <- var + (Sigma[j1, j3] * Sigma[j2, j4] + Sigma[j1, j4] * Sigma[j2, j3]) * 
                N[j1, j3] / N[j1, j2] / N[j3, j4] * theta[j1] * theta[j3] * beta[j2] * beta[j4]
            }
          }
        }
      }else{
        for(j3 in c(j1, j2)){
          for(j4 in c(j1, j2)){
            var <- var + (Sigma[j1, j3] * Sigma[j2, j4] + Sigma[j1, j4] * Sigma[j2, j3]) * 
              N[j1, j2] / N[j1, j2] / N[j3, j4] * theta[j1] * theta[j3] * beta[j2] * beta[j4]
          }
        }
      }
    }
  }
  return(var)
}

proj_cov <- function(S, N, epsilon, mu, init_B, init_Lambda, tol, maxiter){
  ## projection of sample covariance matrix upon positive semi-definite cone w.r.t. weighted l_infty norm
  # S: sample covariance matrix (entry-wise estimation)
  # N: matrix of pair-wise sample sizes
  # init_B and init_Lambda need to be symmetric to ensure the symmetry of Sigma
  B <- init_B;
  Lambda <- init_Lambda;
  original_loss <- NULL; primal_loss <- NULL; Resid_sz <- NULL
  #dual_loss <- NULL
  for(i in 1:maxiter){
    if(i>1){
      Sigma_prev <- Sigma;
    }
    Sigma <- B + S + mu * Lambda;
    eig_out <- eigen(Sigma);
    Sigma <- eig_out[[2]]%*%diag(pmax(eig_out[[1]],epsilon))%*%t(eig_out[[2]])
    #Resid <- Sigma - B - S
    #inner_loss <- c(inner_loss, primal_loss - sum(Lambda * Resid) + 1 / 2 / mu * norm(Resid, "F"))
    A <- Sigma - mu * Lambda - S; B <- A;
    B[N > 0] <- A[N > 0] - proj_wt_l1_vec(A[N > 0], 1 / sqrt(N[N > 0]), mu / 2);
    primal_loss <- c(primal_loss, max(abs(B * sqrt(N)))/2);
    Resid <- Sigma - B - S
    Resid_sz <- c(Resid_sz, norm(Resid, "F"))
    #dual_loss <- inner_loss[inner_iter]
    Lambda_prev <- Lambda;
    Lambda <- Lambda - Resid/mu;
    original_loss <- c(original_loss,max(abs((Sigma - S)*sqrt(N))))
    if(i>1){
      chg_per <- norm(Sigma - Sigma_prev,"F")/norm(Sigma,"F")
      if(chg_per <= tol){
        break;
      }
    }
  }
  return(list(projected_S=Sigma,resid=B,original_loss=original_loss,Resid_sz = Resid_sz, primal_loss <- primal_loss,iteration_num=i,chg_per=chg_per))
}
#----------------------------------------------------#
proj_wt_l1_vec <- function(v,omega,r){
  ## projection of vector v on weighted l1 ball
  # omega: weight vector
  # r: radius
  if(sum(abs(omega*v))<=r){
    return(v)
  }else{
    z <- v/omega;
    sort_ind <- order(z,decreasing = TRUE)
    a <- cumsum((abs(omega*v))[sort_ind])
    b <- cumsum((omega^2)[sort_ind])
    j <- length(z);
    for(i in 1:length(z)){
      j_new <- sum(z[sort_ind]>(a[j]-r)/b[j])
      if(j_new == j){
        break;
      }else{
        j <- j_new;
      }
    }
    c <- (a[j]-r)/b[j];
    proj_v <- pmax(abs(v) - c*omega, 0 )*sign(v)
    return(proj_v)
  }
}
#----------------------------------------------------#

fit_nb <- function(Sigma,a,lambda,eta,tol,...){
  ## use projected gradient descent to solve neighborhood lasso
  # lambda: the p-1 dimensional penalty parameter
  # Sigma: the estimated p by p sample covariance matrix
  # a: target node
  # eta: intial step size parameter
  # ... can be initial beta (p-1 dimensional), if not supplied, beta_init = 0
  p <- dim(Sigma)[1]
  input <- list(...);
  eta_init=eta;
  if(length(input)>1){
    iter <- input[[2]]
    beta_init <- input[[1]]
  }else{
    iter <- 1000
    if(length(input)>0){
      beta_init <- input[[1]]
    }else{
      beta_init <- rep(0,p-1)
    }
  }
  beta <- beta_init
  kk=0;
  chg_per=Inf;
  A <- Sigma[-a,-a];
  b <- Sigma[-a,a]
  
  # initialize loss 
  loss=t(beta)%*%A%*%beta/2-sum(b*beta)+sum(abs(lambda*beta));
  
  while(kk<iter&chg_per>tol){
    grad <- A%*%beta - b
    # determine initial step size^(-1)
    if(kk > 0){
      denom <- sum((beta-beta_prev)*(grad-grad_prev));
      num <- sum((beta-beta_prev)^2);
      eta = num/denom;
    }
    accept <- F;
    while(!accept){
      grad_upd <- beta-grad*eta;
      # soft-thresholding
      beta_temp <- pmax(abs(grad_upd)-eta*lambda,0)*sign(grad_upd)
      loss_temp <- t(beta_temp)%*%A%*%beta_temp/2 - sum(b*beta_temp) + sum(abs(lambda*beta_temp));
      # check acceptance: non-increase in loss
      if(kk==0){
        accept <- (loss_temp <= loss)|eta<=0.001;
      }else{
        accept <- (loss_temp<=loss[kk+1])|eta<=0.001;
      }
      eta <- eta/2;
    }
  # check convergence
  beta_prev <- beta;
  grad_prev <- grad;
  beta <- beta_temp
  nm_prev <- norm(beta_prev,type='2')
  nm_diff <- norm(beta-beta_prev,type='2')
  if(nm_prev > 0){
    chg_per <- nm_diff/nm_prev;
  }else{
    chg_per <- nm_diff;
  }
  
  kk <- kk+1;
  loss <- c(loss,loss_temp);
  }
 if(chg_per <= tol){
   cvg <- T;
 }else{
   cvg=F;
 }
  return(list(beta_hat = beta,loss_path = loss,final_grad = grad,cvg = cvg))
}
fit_nb_zerosamplesz <- function(Sigma,a,lambda,eta,tol,...){
  #if some node has zero sample size with another node, and lambda has an infinit entry, do not conisder this node in the neighborhood regression.
  p <- dim(Sigma)[1]
  input <- list(...);
  if(length(input)>1){
    iter <- input[[2]]
    beta_init <- input[[1]]
  }else{
    iter <- 1000
    if(length(input)>0){
      beta_init <- input[[1]]
    }else{
      beta_init <- rep(0,p-1)
    }
  }
  est_ind <- ((1:p)[-a])[lambda != Inf]
  out <- fit_nb(Sigma[c(a, est_ind), c(a, est_ind)], 1, lambda[lambda != Inf], eta, tol, beta_init[lambda != Inf], iter)
  beta <- rep(0, p - 1); beta[lambda != Inf] <- out$beta_hat; 
  out$beta_hat <- beta;
  return(out)
}




find_ind_delete_node <- function(Ind, p, node_a){
  #find the new list of indices when node_a is deleted from nodes 1:p
  #the new list has p-1 elements, each including p-1 lists of indices
  Ind_del <- list()
  for(i in 1 : (p - 2)){
      if(i < node_a){
        Ind_del[[i]] <- Ind[[i]][-node_a]
     }else{
        Ind_del[[i]] <- Ind[[i + 1]][-node_a]
     }
  }
  return(Ind_del)
}

