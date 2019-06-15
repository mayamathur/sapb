
# modified from package weightr to also provide Hessian

# prepare to assign weightfunct_mm as a new function
#https://stackoverflow.com/questions/53063629/my-modified-version-of-package-function-cant-find-the-packages-other-internal
#environment(weightfunct)<-as.environment("package:weightr")


weightfunct_mm = function (effect, v, steps = c(0.025, 1), mods = NULL, weights = NULL, 
          fe = FALSE, table = FALSE, pval = NULL) 
{
  
  warning("Maya's version of the function")
  
  neglike_unadj <- function(pars) {
    if (fe == FALSE) {
      vc = pars[1]
      beta = pars[2:(npred + 2)]
      mn = XX %*% beta
      eta = sqrt(v + vc)
      b = 1/2 * sum(((effect - mn)/eta)^2)
      c = sum(log(eta))
    }
    else {
      beta = pars[1:(npred + 1)]
      mn = XX %*% beta
      eta = sqrt(v)
      b = 1/2 * sum(((effect - mn)/eta)^2)
      c = sum(log(eta))
    }
    return(b + c)
  }
  neglike_adj <- function(pars) {
    if (fe == FALSE) {
      vc = pars[1]
      beta = pars[2:(npred + 2)]
      if (is.null(weights) == FALSE) {
        w = weights
      }
      else {
        w = c(1, pars[(npred + 3):((nsteps - 2) + (npred + 
                                                     3))])
      }
      contrib = log(w[wt])
      mn = XX %*% beta
      a = sum(contrib)
      eta = sqrt(v + vc)
    }
    else {
      beta = pars[1:(npred + 1)]
      if (is.null(weights) == FALSE) {
        w = weights
      }
      else {
        w = c(1, pars[(npred + 2):((nsteps - 2) + (npred + 
                                                     2))])
      }
      contrib = log(w[wt])
      mn = XX %*% beta
      a = sum(contrib)
      eta = sqrt(v)
    }
    b = 1/2 * sum(((effect - mn)/eta)^2)
    c = sum(log(eta))
    Bij <- matrix(rep(0, number * nsteps), nrow = number, 
                  ncol = nsteps)
    bi = -si * qnorm(steps[1])
    Bij[, 1] = 1 - pnorm((bi - mn)/eta)
    if (nsteps > 2) {
      for (j in 2:(length(steps) - 1)) {
        bi = -si * qnorm(steps[j])
        bilast = -si * qnorm(steps[j - 1])
        Bij[, j] = pnorm((bilast - mn)/eta) - pnorm((bi - 
                                                       mn)/eta)
      }
    }
    bilast = -si * qnorm(steps[length(steps) - 1])
    Bij[, length(steps)] = pnorm((bilast - mn)/eta)
    swbij = 0
    for (j in 1:length(steps)) swbij = swbij + w[j] * Bij[, 
                                                          j]
    d = sum(log(swbij))
    return(-a + b + c + d)
  }
  if (is.null(mods)) {
    npred <- 0
    data <- data.frame(effect, v)
  }
  else {
    if (typeof(mods) == "language") {
      XX <- model.matrix(mods, model.frame(mods, na.action = "na.pass"))
      npred <- dim(XX)[2] - 1
      data <- data.frame(effect, v, XX)
    }
  }
  if (any(is.na(data))) {
    data <- na.omit(data)
    removed <- as.numeric(na.action(data))
  }
  else {
    removed <- NULL
  }
  effect <- data[, 1]
  v <- data[, 2]
  if (npred == 0) {
    XX <- cbind(rep(1, length(effect)))
  }
  else {
    XX <- as.matrix(data[, (3:(npred + 3))])
  }
  if (length(effect) != length(v)) {
    stop("Your vector of effect sizes and your vector of sampling variances are not the same length. Please check your data.")
  }
  if (identical(effect, v)) {
    stop("Your vector of effect sizes is exactly the same as your vector of sampling variances. Please check your data.")
  }
  if (min(v) < 0) {
    stop("Sampling variances cannot be negative. Please check your data.")
  }
  si <- sqrt(v)
  if (is.null(pval)) {
    p <- 1 - pnorm(effect/sqrt(v))
  }
  else {
    p <- pval
  }
  if (max(steps) != 1) {
    steps <- c(steps, 1)
  }
  if (max(steps) > 1) {
    stop("p-value cutpoints cannot be greater than 1.")
  }
  if (min(steps) < 0) {
    stop("p-value cutpoints cannot be negative.")
  }
  if (length(unique(steps)) != length(steps)) {
    stop("Two or more p-value cutpoints are identical.")
  }
  if (is.null(weights)) {
    steps <- sort(steps)
  }
  if (is.null(weights) == FALSE) {
    if (min(weights) < 0) {
      stop("Weights for p-value intervals cannot be negative.")
    }
    if (length(weights) != length(steps)) {
      stop("The number of weights does not match the number of p-value intervals created.")
    }
    new <- cbind(steps, weights)
    steps <- new[order(steps), 1]
    weights <- new[order(steps), 2]
  }
  number <- length(effect)
  nsteps <- length(steps)
  wt <- rep(1, number)
  for (i in 1:number) {
    if (p[i] <= steps[1]) 
      wt[i] = 1
    for (j in 2:nsteps) {
      if (steps[j - 1] < p[i] && p[i] <= steps[j]) 
        wt[i] = j
    }
    if (p[i] > steps[nsteps - 1]) 
      wt[i] = nsteps
  }
  intervaltally <- function(p, steps) {
    p1 <- cut(p, breaks = c(-Inf, steps), labels = steps)
    return(p1)
  }
  pvalues <- as.numeric(table(intervaltally(p, steps)))
  sampletable <- function(p, pvalues, steps) {
    nsteps <- length(steps)
    results <- matrix(nrow = length(pvalues), ncol = 1)
    results[, 1] <- pvalues
    rowlabels <- c(0, length(results[, 1]))
    rowlabels[1] <- paste(c("p-values <", steps[1]), collapse = "")
    for (i in 2:nsteps) {
      rowlabels[i] <- paste(c(steps[i - 1], "< p-values <", 
                              steps[i]), collapse = " ")
    }
    resultsb <- data.frame(results, row.names = c(rowlabels))
    colnames(resultsb) <- c("Frequency")
    return(resultsb)
  }
  if (sum(table(intervaltally(p, steps)) == 0) >= 1) {
    warning("At least one of the p-value intervals contains no effect sizes, leading to estimation problems. Consider re-specifying the cutpoints.")
  }
  if (sum(table(intervaltally(p, steps)) > 0 & table(intervaltally(p, 
                                                                   steps)) <= 3) >= 1) {
    warning("At least one of the p-value intervals contains three or fewer effect sizes, which may lead to estimation problems. Consider re-specifying the cutpoints.")
  }
  if (is.null(mods)) {
    if (fe == FALSE) {
      pars <- c(mean(v)/4, mean(effect), rep(0, (nsteps - 
                                                   1)))
      output_unadj <- optim(par = pars[1:2], fn = neglike_unadj, 
                            lower = c(0, -Inf), method = "L-BFGS-B", hessian = TRUE)
      output_adj <- optim(par = pars, fn = neglike_adj, 
                          lower = c(0, -Inf, rep(0.01, (nsteps - 1))), 
                          method = "L-BFGS-B", hessian = TRUE)
      results <- list(output_unadj, output_adj, steps = steps, 
                      mods = mods, weights = weights, fe = fe, table = table, 
                      effect = effect, v = v, npred = npred, nsteps = nsteps, 
                      p = p, XX = XX, removed = removed)
      if (is.null(weights)) {
        suppressWarnings(results2 <- list(unadj_est = c(output_unadj$par), 
                                          unadj_se = c(sqrt(diag(solve(output_unadj$hessian)))), 
                                          adj_est = c(output_adj$par), adj_se = c(sqrt(diag(solve(output_adj$hessian)))), 
                                          z_unadj = c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))), 
                                          z_adj = c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))), 
                                          p_unadj = c(2 * pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))), 
                                          p_adj = c(2 * pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))), 
                                          ci.lb_unadj = c(output_unadj$par - qnorm(0.975) * 
                                                            sqrt(diag(solve(output_unadj$hessian)))), 
                                          ci.ub_unadj = c(output_unadj$par + qnorm(0.975) * 
                                                            sqrt(diag(solve(output_unadj$hessian)))), 
                                          ci.lb_adj = c(output_adj$par - qnorm(0.975) * 
                                                          sqrt(diag(solve(output_adj$hessian)))), ci.ub_adj = c(output_adj$par + 
                                                                                                                  qnorm(0.975) * sqrt(diag(solve(output_adj$hessian))))))
      }
      if (is.null(weights) == FALSE) {
        results2 <- list(unadj_est = c(output_unadj$par), 
                         unadj_se = c(sqrt(diag(solve(output_unadj$hessian)))), 
                         adj_est = c(output_adj$par), z_unadj = c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))), 
                         p_unadj = c(2 * pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))), 
                         ci.lb_unadj = c(output_unadj$par - qnorm(0.975) * 
                                           sqrt(diag(solve(output_unadj$hessian)))), 
                         ci.ub_unadj = c(output_unadj$par + qnorm(0.975) * 
                                           sqrt(diag(solve(output_unadj$hessian)))))
      }
    }
    if (fe == TRUE) {
      pars <- c(mean(effect), rep(1, (nsteps - 1)))
      output_unadj <- optim(par = pars[1], fn = neglike_unadj, 
                            lower = c(-Inf), method = "L-BFGS-B", hessian = TRUE)
      output_adj <- optim(par = pars, fn = neglike_adj, 
                          lower = c(-Inf, rep(0.01, (nsteps - 1))), method = "L-BFGS-B", 
                          hessian = TRUE)
      results <- list(output_unadj, output_adj, steps = steps, 
                      mods = mods, weights = weights, fe = fe, table = table, 
                      effect = effect, v = v, npred = npred, nsteps = nsteps, 
                      p = p, removed = removed)
      if (is.null(weights)) {
        suppressWarnings(results2 <- list(unadj_est = c(output_unadj$par), 
                                          unadj_se = c(sqrt(diag(solve(output_unadj$hessian)))), 
                                          adj_est = c(output_adj$par), adj_se = c(sqrt(diag(solve(output_adj$hessian)))), 
                                          z_unadj = c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))), 
                                          z_adj = c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))), 
                                          p_unadj = c(2 * pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))), 
                                          p_adj = c(2 * pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))), 
                                          ci.lb_unadj = c(output_unadj$par - qnorm(0.975) * 
                                                            sqrt(diag(solve(output_unadj$hessian)))), 
                                          ci.ub_unadj = c(output_unadj$par + qnorm(0.975) * 
                                                            sqrt(diag(solve(output_unadj$hessian)))), 
                                          ci.lb_adj = c(output_adj$par - qnorm(0.975) * 
                                                          sqrt(diag(solve(output_adj$hessian)))), ci.ub_adj = c(output_adj$par + 
                                                                                                                  qnorm(0.975) * sqrt(diag(solve(output_adj$hessian))))))
      }
      if (is.null(weights) == FALSE) {
        results2 <- list(unadj_est = c(output_unadj$par), 
                         unadj_se = c(sqrt(diag(solve(output_unadj$hessian)))), 
                         adj_est = c(output_adj$par), z_unadj = c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))), 
                         p_unadj = c(2 * pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))), 
                         ci.lb_unadj = c(output_unadj$par - qnorm(0.975) * 
                                           sqrt(diag(solve(output_unadj$hessian)))), 
                         ci.ub_unadj = c(output_unadj$par + qnorm(0.975) * 
                                           sqrt(diag(solve(output_unadj$hessian)))))
      }
    }
  }
  else {
    if (typeof(mods) == "language") {
      if (fe == FALSE) {
        pars <- c(mean(v)/4, mean(effect), rep(0, npred), 
                  rep(1, (nsteps - 1)))
        output_unadj <- optim(par = pars[1:(npred + 2)], 
                              fn = neglike_unadj, lower = c(0, rep(-Inf, 
                                                                   (npred + 1))), method = "L-BFGS-B", hessian = TRUE)
        output_adj <- optim(par = pars, fn = neglike_adj, 
                            lower = c(0, rep(-Inf, (npred + 1)), rep(0.01, 
                                                                     (nsteps - 1))), method = "L-BFGS-B", hessian = TRUE)
        results <- list(output_unadj, output_adj, steps = steps, 
                        mods = mods, weights = weights, fe = fe, table = table, 
                        effect = effect, v = v, npred = npred, nsteps = nsteps, 
                        p = p, XX = XX, removed = removed)
        if (is.null(weights)) {
          suppressWarnings(results2 <- list(unadj_est = c(output_unadj$par), 
                                            unadj_se = c(sqrt(diag(solve(output_unadj$hessian)))), 
                                            adj_est = c(output_adj$par), adj_se = c(sqrt(diag(solve(output_adj$hessian)))), 
                                            z_unadj = c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))), 
                                            z_adj = c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))), 
                                            p_unadj = c(2 * pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))), 
                                            p_adj = c(2 * pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))), 
                                            ci.lb_unadj = c(output_unadj$par - qnorm(0.975) * 
                                                              sqrt(diag(solve(output_unadj$hessian)))), 
                                            ci.ub_unadj = c(output_unadj$par + qnorm(0.975) * 
                                                              sqrt(diag(solve(output_unadj$hessian)))), 
                                            ci.lb_adj = c(output_adj$par - qnorm(0.975) * 
                                                            sqrt(diag(solve(output_adj$hessian)))), 
                                            ci.ub_adj = c(output_adj$par + qnorm(0.975) * 
                                                            sqrt(diag(solve(output_adj$hessian))))))
        }
        if (is.null(weights) == FALSE) {
          suppressWarnings(results2 <- list(unadj_est = c(output_unadj$par), 
                                            unadj_se = c(sqrt(diag(solve(output_unadj$hessian)))), 
                                            adj_est = c(output_adj$par), z_unadj = c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))), 
                                            p_unadj = c(2 * pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))), 
                                            ci.lb_unadj = c(output_unadj$par - qnorm(0.975) * 
                                                              sqrt(diag(solve(output_unadj$hessian)))), 
                                            ci.ub_unadj = c(output_unadj$par + qnorm(0.975) * 
                                                              sqrt(diag(solve(output_unadj$hessian))))))
        }
      }
      if (fe == TRUE) {
        pars <- c(mean(effect), rep(0, npred), rep(1, 
                                                   (nsteps - 1)))
        output_unadj <- optim(par = pars[1:(npred + 1)], 
                              fn = neglike_unadj, lower = c(rep(-Inf, (npred + 
                                                                         1))), method = "L-BFGS-B", hessian = TRUE)
        output_adj <- optim(par = pars, fn = neglike_adj, 
                            lower = c(rep(-Inf, (npred + 1)), rep(0.01, 
                                                                  (nsteps - 1))), method = "L-BFGS-B", hessian = TRUE)
        results <- list(output_unadj, output_adj, steps = steps, 
                        mods = mods, weights = weights, fe = fe, table = table, 
                        effect = effect, v = v, npred = npred, nsteps = nsteps, 
                        p = p, XX = XX, removed = removed)
        if (is.null(weights)) {
          suppressWarnings(results2 <- list(unadj_est = c(output_unadj$par), 
                                            unadj_se = c(sqrt(diag(solve(output_unadj$hessian)))), 
                                            adj_est = c(output_adj$par), adj_se = c(sqrt(diag(solve(output_adj$hessian)))), 
                                            z_unadj = c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))), 
                                            z_adj = c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))), 
                                            p_unadj = c(2 * pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))), 
                                            p_adj = c(2 * pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))), 
                                            ci.lb_unadj = c(output_unadj$par - qnorm(0.975) * 
                                                              sqrt(diag(solve(output_unadj$hessian)))), 
                                            ci.ub_unadj = c(output_unadj$par + qnorm(0.975) * 
                                                              sqrt(diag(solve(output_unadj$hessian)))), 
                                            ci.lb_adj = c(output_adj$par - qnorm(0.975) * 
                                                            sqrt(diag(solve(output_adj$hessian)))), 
                                            ci.ub_adj = c(output_adj$par + qnorm(0.975) * 
                                                            sqrt(diag(solve(output_adj$hessian)))) ) )
        }
        if (is.null(weights) == FALSE) {
          
          # BOOKMARK: HERE, RETURN adj_se = c(sqrt(diag(solve(output_unadj$hessian[1:2,1:2]))))
          # I THINK THIS IS FINE BECAUSE BASICALLY TREATS SES FOR SELECTION PROB AS 0
          
          warning("adj_se=", sqrt(diag(solve(output_unadj$hessian[1:2,1:2]))) )
          
          # ~~~ MM ADDED THE ADJ INFERENCE TO THIS
          suppressWarnings(results2 <- list(unadj_est = c(output_unadj$par), 
                                            unadj_se = c(sqrt(diag(solve(output_unadj$hessian)))), 
                                            
                                            adj_est = c(output_adj$par), z_unadj = c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))), 
                                            p_unadj = c(2 * pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))), 
                                            
                                            ci.lb_unadj = c(output_unadj$par - qnorm(0.975) * 
                                                              sqrt(diag(solve(output_unadj$hessian)))), 
                                            ci.ub_unadj = c(output_unadj$par + qnorm(0.975) * 
                                                              sqrt(diag(solve(output_unadj$hessian)))),
                           
                                            ci.lb_adj = -9 )
                           )
            
                                     
        }
      }
    }
    else {
      stop("Moderators must be entered as \"mods= ~ x1 + x2\"")
    }
  }
  class(results) <- c("weightfunct")
  
  # ~~~ MM CHANGED THIS TO ALSO SPIT OUT THE ADJUSTED HESSIAN
  return( list( results = results, output_adj = output_adj ) )
}


# assign weightfunct_mm 
#unlockBinding("weightfunct", as.environment("package:weightr"))
#assign("weightfunct", weightfunct_mm, as.environment("package:weightr"))
#lockBinding("weightfunct", as.environment("package:weightr"))