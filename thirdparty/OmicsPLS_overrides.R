# GPL-3, Copyright Said el Bouhaddani
# This file inclues patches applied/developed but not yet publishd to CRAN;
# this files is to be removed once all the fixes make it into the release.
# Some of the fixes and changes authored by M. Krassowski

library(OmicsPLS)


fixed_pow_o2m <- function(X, Y, n, tol = 1e-10, max_iterations = 100) {
  input_checker(X, Y)
  stopifnot(n == round(n))
  #  message("High dimensional problem: switching to power method.\n")
  #  message("initialize Power Method. Stopping crit: sq.err<", tol, " or ", max_iterations, " iter.\n")
  Tt <- NULL
  U <- NULL
  W <- NULL
  C <- NULL
  for (indx in 1:n) {
    Ui = rowSums(Y)
    for (indx2 in 1:max_iterations) {
      tmpp <- Ui
      Wi = crossprod(X, Ui)
      Wi = Wi / c(vnorm(Wi))
      Ti = X %*% Wi
      Ci = crossprod(Y, Ti)
      Ci = Ci / c(vnorm(Ci))
      Ui = Y %*% Ci
      if (mse(tmpp, Ui) < tol) {
        break
      }
    }
    if(ssq(Wi) < 1e-10 || ssq(Ci) < 1e-10){
      Wi <- orth(rep(1,ncol(X)))
      Ci <- orth(rep(1,ncol(Y)))
      for (indx2 in 1:max_iterations) {
        tmpp <- c(Wi, Ci)
        Wi <- orth(t(X) %*% (Y %*% t(Y)) %*% (X %*% Wi))
        Ci <- orth(t(Y) %*% (X %*% t(X)) %*% (Y %*% Ci))
        if (mse(tmpp, c(Wi, Ci)) < tol) {
          message("The initialization of the power method lied in a degenerate space\n")
          message("Initialization changed and power method rerun\n")
          break
        }
      }
    }
    message("Power Method (comp ", indx, ") stopped after ", indx2, " iterations.\n")
    Tt <- cbind(Tt, X %*% Wi)
    U <- cbind(U, Y %*% Ci)
    X <- X - (X %*% Wi) %*% t(Wi)
    Y <- Y - (Y %*% Ci) %*% t(Ci)
    W <- cbind(W, Wi)
    C <- cbind(C, Ci)
  }
  return(list(W = W, C = C, Tt = Tt, U = U))
}


fixed_o2m2 <- function(X, Y, n, nx, ny, stripped = TRUE, tol = 1e-10, max_iterations = 100) {
  
  Xnames = dimnames(X)
  Ynames = dimnames(Y)
  
  X_true <- X
  Y_true <- Y
  
  N <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  
  T_Yosc <- U_Xosc <- matrix(0, N, n)
  W_Yosc <- P_Yosc <- matrix(0, p, n)
  C_Xosc <- P_Xosc <- matrix(0, q, n)
  
  if (nx + ny > 0) {
    # larger principal subspace
    n2 <- n + max(nx, ny)
    
    # if(N<p&N<q){ # When N is smaller than p and q
    #
    # CHANGE
    #
    W_C <- pow_o2m(X, Y, n, tol, max_iterations)
    W <- W_C$W
    C <- W_C$C
    Tt <- W_C$Tt
    U <- W_C$U
    # } cdw = svd(t(Y)%*%X,nu=n2,nv=n2); C=cdw$uW=cdw$v
    
    # Tt = X%*%W;
    
    if (nx > 0) {
      # Orthogonal components in Y
      E_XY <- X - Tt %*% t(W)
      
      udv <- svd(t(E_XY) %*% Tt, nu = nx, nv = 0)
      W_Yosc <- udv$u
      T_Yosc <- X %*% W_Yosc
      P_Yosc <- t(solve(t(T_Yosc) %*% T_Yosc) %*% t(T_Yosc) %*% X)
      X <- X - T_Yosc %*% t(P_Yosc)
      
      # Update T again Tt = X%*%W;
    }
    
    # U = Y%*%C; # 3.2.1. 4
    
    if (ny > 0) {
      # Orthogonal components in Y
      F_XY <- Y - U %*% t(C)
      
      udv <- svd(t(F_XY) %*% U, nu = ny, nv = 0)
      C_Xosc <- udv$u
      U_Xosc <- Y %*% C_Xosc
      P_Xosc <- t(solve(t(U_Xosc) %*% U_Xosc) %*% t(U_Xosc) %*% Y)
      Y <- Y - U_Xosc %*% t(P_Xosc)
      
      # Update U again U = Y%*%C;
    }
  }
  # Re-estimate joint part in n-dimensional subspace if(N<p&N<q){ # When N is smaller than p and q
  W_C <- pow_o2m(X, Y, n, tol, max_iterations)
  W <- W_C$W
  C <- W_C$C
  Tt <- W_C$Tt
  U <- W_C$U
  # } cdw = svd(t(Y)%*%X,nu=n,nv=n); # 3.2.1. 1 C=cdw$u;W=cdw$v Tt = X%*%W; # 3.2.1. 2 U = Y%*%C; #
  # 3.2.1. 4
  
  # Inner relation parameters
  B_U <- solve(t(U) %*% U) %*% t(U) %*% Tt
  B_T <- solve(t(Tt) %*% Tt) %*% t(Tt) %*% U
  
  # Residuals and R2's
  if(stripped){
    E <- Ff <- X_hat <- Y_hat <- as.matrix(0)
  } else {
    E <- X_true - Tt %*% t(W) - T_Yosc %*% t(P_Yosc)
    Ff <- Y_true - U %*% t(C) - U_Xosc %*% t(P_Xosc)
    Y_hat <- Tt %*% B_T %*% t(C)
    X_hat <- U %*% B_U %*% t(W)
  }
  H_TU <- Tt - U %*% B_U
  H_UT <- U - Tt %*% B_T
  
  # R2
  R2Xcorr <- (ssq(Tt)/ssq(X_true))
  R2Ycorr <- (ssq(U)/ssq(Y_true))
  R2X_YO <- (ssq(T_Yosc)/ssq(X_true))
  R2Y_XO <- (ssq(U_Xosc)/ssq(Y_true))
  R2Xhat <- (ssq(U %*% B_U)/ssq(X_true))
  R2Yhat <- (ssq(Tt %*% B_T)/ssq(Y_true))
  R2X <- R2Xcorr + R2X_YO
  R2Y <- R2Ycorr + R2Y_XO
  
  rownames(Tt) <- rownames(T_Yosc) <- rownames(H_TU) <- Xnames[[1]]
  rownames(U) <- rownames(U_Xosc) <- rownames(H_UT) <- Ynames[[1]]
  rownames(W) <- rownames(P_Yosc) <- rownames(W_Yosc) <- Xnames[[2]]
  rownames(C) <- rownames(P_Xosc) <- rownames(C_Xosc) <- Ynames[[2]]
  
  model <- list(Tt = Tt, W. = W, U = U, C. = C, E = E, Ff = Ff, T_Yosc = T_Yosc, P_Yosc. = P_Yosc, W_Yosc = W_Yosc, 
                U_Xosc = U_Xosc, P_Xosc. = P_Xosc, C_Xosc = C_Xosc, B_U = B_U, B_T. = B_T, H_TU = H_TU, H_UT = H_UT, 
                X_hat = X_hat, Y_hat = Y_hat, R2X = R2X, R2Y = R2Y, R2Xcorr = R2Xcorr, R2Ycorr = R2Ycorr, R2X_YO = R2X_YO, 
                R2Y_XO = R2Y_XO, R2Xhat = R2Xhat, R2Yhat = R2Yhat)
  class(model) <- "pre.o2m"
  return(model)
}


fast_svd = function (X, Y, n, nx, ny, stripped = FALSE, tol = 1e-10, max_iterations = 100) 
{
    Xnames = dimnames(X)
    Ynames = dimnames(Y)
    if (stripped) {
        return(o2m_stripped2(X, Y, n, nx, ny, tol, max_iterations))
    }
    X_true <- X
    Y_true <- Y
    N <- nrow(X)
    p <- ncol(X)
    q <- ncol(Y)
    T_Yosc <- U_Xosc <- matrix(0, N, n)
    W_Yosc <- P_Yosc <- matrix(0, p, n)
    C_Xosc <- P_Xosc <- matrix(0, q, n)

    if (nx + ny > 0) {
        n2 <- n + max(nx, ny)
        cdw <- svd(t(Y) %*% X, nu =n2, nv = n2)
        
        if (nx > 0) {
            
            W <- cdw$v
            Tt <- X %*% W
            E_XY <- X - Tt %*% t(W)
            udv <- svd(t(E_XY) %*% Tt, nu = nx, nv = 0)

            W_Yosc <- udv$u
            
            T_Yosc <- X %*% W_Yosc
            P_Yosc <- t(solve(t(T_Yosc) %*% T_Yosc) %*% t(T_Yosc) %*% X)
            X <- X - T_Yosc %*% t(P_Yosc)
        }
        if (ny > 0) {
            C <- cdw$u
            U <- Y %*% C
            F_XY <- Y - U %*% t(C)
            udv <- svd(t(F_XY) %*% U, nu = ny, nv = 0)

            C_Xosc <- udv$u
            U_Xosc <- Y %*% C_Xosc
            P_Xosc <- t(solve(t(U_Xosc) %*% U_Xosc) %*% t(U_Xosc) %*% Y)
            Y <- Y - U_Xosc %*% t(P_Xosc)
        }
    }

    cdw <- svd(t(Y) %*% X, nu = min(n, ncol(X)), nv = min(n, ncol(Y)))
    C <- cdw$u
    W <- cdw$v
    Tt <- X %*% W
    U <- Y %*% C
        
    B_U <- solve(t(U) %*% U) %*% t(U) %*% Tt
    B_T <- solve(t(Tt) %*% Tt) %*% t(Tt) %*% U
    H_TU <- Tt - U %*% B_U
    H_UT <- U - Tt %*% B_T
    R2Xcorr <- (ssq(Tt)/ssq(X_true))
    R2Ycorr <- (ssq(U)/ssq(Y_true))
    R2X_YO <- (ssq(T_Yosc)/ssq(X_true))
    R2Y_XO <- (ssq(U_Xosc)/ssq(Y_true))
    R2Xhat <- (ssq(U %*% B_U)/ssq(X_true))
    R2Yhat <- (ssq(Tt %*% B_T)/ssq(Y_true))
    R2X <- R2Xcorr + R2X_YO
    R2Y <- R2Ycorr + R2Y_XO
    rownames(Tt) <- rownames(T_Yosc) <- rownames(H_TU) <- Xnames[[1]]
    rownames(U) <- rownames(U_Xosc) <- rownames(H_UT) <- Ynames[[1]]
    rownames(W) <- rownames(P_Yosc) <- rownames(W_Yosc) <- Xnames[[2]]
    rownames(C) <- rownames(P_Xosc) <- rownames(C_Xosc) <- Ynames[[2]]

    model <- list(Tt = Tt, W. = W, U = U, C. = C, E = 0, Ff = 0, 
        T_Yosc = T_Yosc, P_Yosc. = P_Yosc, W_Yosc = W_Yosc, U_Xosc = U_Xosc, 
        P_Xosc. = P_Xosc, C_Xosc = C_Xosc, B_U = B_U, B_T. = B_T, 
        H_TU = H_TU, H_UT = H_UT, X_hat = 0, Y_hat = 0, R2X = R2X, 
        R2Y = R2Y, R2Xcorr = R2Xcorr, R2Ycorr = R2Ycorr, R2X_YO = R2X_YO, 
        R2Y_XO = R2Y_XO, R2Xhat = R2Xhat, R2Yhat = R2Yhat)
    class(model) <- "o2m"
    return(model)
}


predict_o2m <- function(object, newdata, XorY = c("X","Y"), ...) {
  XorY = match.arg(XorY)
  Xnames = dimnames(newdata)
  if(!is.matrix(newdata)){
    message("newdata has class ",class(newdata),", trying to convert with as.matrix.",sep="")
    newdata <- as.matrix(newdata)
    dimnames(newdata) <- Xnames
  }
  input_checker(newdata)
  
  pred = switch(XorY, 
                Y = with(object, (newdata - newdata %*% C_Xosc %*% t(C_Xosc)) %*% C. %*% B_U %*% t(W.)), 
                X = with(object, (newdata - newdata %*% W_Yosc %*% t(W_Yosc)) %*% W. %*% B_T. %*% t(C.)))
  
  return(pred)
}


fixed_o2m = function (X, Y, n, nx, ny, stripped = FALSE, p_thresh = 3000, 
    q_thresh = p_thresh, tol = 1e-10, max_iterations = 100) 
{
    tic <- proc.time()
    Xnames = dimnames(X)
    Ynames = dimnames(Y)
    if (!is.matrix(X)) {
        message("X has class ", class(X), ", trying to convert with as.matrix.", 
            sep = "")
        X <- as.matrix(X)
        dimnames(X) <- Xnames
    }
    if (!is.matrix(Y)) {
        message("Y has class ", class(Y), ", trying to convert with as.matrix.", 
            sep = "")
        Y <- as.matrix(Y)
        dimnames(Y) <- Ynames
    }
    input_checker(X, Y)
    ssqX = ssq(X)
    ssqY = ssq(Y)
    if (length(n) > 1 | length(nx) > 1 | length(ny) > 1) 
        stop("Number of components should be scalars, not vectors")
    #if (ncol(X) < n + max(nx, ny) || ncol(Y) < n + max(nx, ny)) 
    #    stop("n + max(nx, ny) =", n + max(nx, ny), " exceed # columns in X or Y")
    
    if (max(ncol(X), ncol(Y)) < n)
        stop("n =", n, " exceed # columns in X or Y")
    if (ncol(X) < nx || ncol(Y) < ny) 
        stop("nx = ", nx, " or ny = ", ny, " exceed # columns in X or Y, respectively")
    
    if (nx != round(abs(nx)) || ny != round(abs(ny))) 
        stop("n, nx and ny should be non-negative integers")
    if (p_thresh != round(abs(p_thresh)) || q_thresh != round(abs(q_thresh))) 
        stop("p_thresh and q_thresh should be non-negative integers")
    if (max_iterations != round(abs(max_iterations))) 
        stop("max_iterations should be a non-negative integer")
    if (tol < 0) 
        stop("tol should be non-negative")
    if (nrow(X) < n + max(nx, ny)) 
        stop("n + max(nx, ny) = ", n + max(nx, ny), " exceed sample size N = ", 
            nrow(X))
    if (nrow(X) == n + max(nx, ny)) 
        warning("n + max(nx, ny) = ", n + max(nx, ny), " equals sample size")
    if (n != round(abs(n)) || n <= 0) {
        stop("n should be a positive integer")
    }
    if (any(abs(colMeans(X)) > 1e-05)) {
        message("Data is not centered, proceeding...")
    }
    highd = FALSE
    if ((ncol(X) > p_thresh && ncol(Y) > q_thresh)) {
        highd = TRUE
        message("Using high dimensional mode with tolerance ", 
            tol, " and max iterations ", max_iterations)
        model = o2m2(X, Y, n, nx, ny, stripped, tol, max_iterations)
    }
    else if (stripped) {
        model = o2m_stripped(X, Y, n, nx, ny)
    }

    else {
        X_true <- X
        Y_true <- Y
        N <- nrow(X)
        p <- ncol(X)
        q <- ncol(Y)
        T_Yosc <- U_Xosc <- matrix(0, N, n)
        W_Yosc <- P_Yosc <- matrix(0, p, n)
        C_Xosc <- P_Xosc <- matrix(0, q, n)
        if (nx + ny > 0) {
            # why?????
            n2 <- n + max(nx, ny)
            cdw <- svd(t(Y) %*% X, nu = n2, nv = n2)
            C <- cdw$u
            W <- cdw$v
            if (nx > 0) {
                Tt <- X %*% W
                E_XY <- X - Tt %*% t(W)
                udv <- svd(t(E_XY) %*% Tt, nu = nx, nv = 0)
                W_Yosc <- udv$u
                T_Yosc <- X %*% W_Yosc
                P_Yosc <- t(solve(t(T_Yosc) %*% T_Yosc) %*% t(T_Yosc) %*% 
                  X)
                X <- X - T_Yosc %*% t(P_Yosc)
                # Tt <- X %*% W
            }
            if (ny > 0) {
                U <- Y %*% C
                F_XY <- Y - U %*% t(C)
                udv <- svd(t(F_XY) %*% U, nu = ny, nv = 0)
                C_Xosc <- udv$u
                U_Xosc <- Y %*% C_Xosc
                P_Xosc <- t(solve(t(U_Xosc) %*% U_Xosc) %*% t(U_Xosc) %*% 
                  Y)
                Y <- Y - U_Xosc %*% t(P_Xosc)
                # U <- Y %*% C
            }
        }
        cdw <- svd(t(Y) %*% X, nu = n, nv = n)
        C <- cdw$u
        W <- cdw$v
        Tt <- X %*% W
        U <- Y %*% C
        B_U <- solve(t(U) %*% U) %*% t(U) %*% Tt
        B_T <- solve(t(Tt) %*% Tt) %*% t(Tt) %*% U
        E <- X_true - Tt %*% t(W) - T_Yosc %*% t(P_Yosc)
        Ff <- Y_true - U %*% t(C) - U_Xosc %*% t(P_Xosc)
        H_TU <- Tt - U %*% B_U
        H_UT <- U - Tt %*% B_T
        Y_hat <- Tt %*% B_T %*% t(C)
        X_hat <- U %*% B_U %*% t(W)
        R2Xcorr <- (ssq(Tt)/ssqX)
        R2Ycorr <- (ssq(U)/ssqY)
        R2X_YO <- (ssq(T_Yosc)/ssqX)
        R2Y_XO <- (ssq(U_Xosc)/ssqY)
        R2Xhat <- (ssq(U %*% B_U)/ssqX)
        R2Yhat <- (ssq(Tt %*% B_T)/ssqY)
        R2X <- R2Xcorr + R2X_YO
        R2Y <- R2Ycorr + R2Y_XO
        rownames(Tt) <- rownames(T_Yosc) <- rownames(E) <- rownames(H_TU) <- Xnames[[1]]
        rownames(U) <- rownames(U_Xosc) <- rownames(Ff) <- rownames(H_UT) <- Ynames[[1]]
        rownames(W) <- rownames(P_Yosc) <- rownames(W_Yosc) <- colnames(E) <- Xnames[[2]]
        rownames(C) <- rownames(P_Xosc) <- rownames(C_Xosc) <- colnames(Ff) <- Ynames[[2]]
        model <- list(Tt = Tt, W. = W, U = U, C. = C, E = E, 
            Ff = Ff, T_Yosc = T_Yosc, P_Yosc. = P_Yosc, W_Yosc = W_Yosc, 
            U_Xosc = U_Xosc, P_Xosc. = P_Xosc, C_Xosc = C_Xosc, 
            B_U = B_U, B_T. = B_T, H_TU = H_TU, H_UT = H_UT, 
            X_hat = X_hat, Y_hat = Y_hat, R2X = R2X, R2Y = R2Y, 
            R2Xcorr = R2Xcorr, R2Ycorr = R2Ycorr, R2X_YO = R2X_YO, 
            R2Y_XO = R2Y_XO, R2Xhat = R2Xhat, R2Yhat = R2Yhat)
        class(model) <- "o2m"
    }
    toc <- proc.time() - tic
    model$flags = c(time = toc[3], list(n = n, nx = nx, ny = ny, 
        stripped = stripped, highd = highd, call = match.call(), 
        ssqX = ssqX, ssqY = ssqY, varXjoint = apply(model$Tt, 
            2, ssq), varYjoint = apply(model$U, 2, ssq), varXorth = apply(model$P_Y, 
            2, ssq) * apply(model$T_Y, 2, ssq), varYorth = apply(model$P_X, 
            2, ssq) * apply(model$U_X, 2, ssq)))
    return(model)
}