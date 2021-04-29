gwasQxE = function(geno,
                   pheno,
                   trait,
                   Env,
                   n.PC,
                   Fixed = NULL,
                   n.core = 1,
                   scale.pheno.by.trial = TRUE)
{
  makeZ = function(row.id, col.id) {
    Z = matrix(0, nrow=length(row.id), ncol=length(col.id))
    rownames(Z) = row.id
    colnames(Z) = col.id
    Z[cbind(1:nrow(Z),match(row.id,col.id))] = 1
    return(Z)
  }
  myAmat = function(X) {
    X = X + 1
    n = nrow(X)
    m = ncol(X)
    v1 = matrix(1, n, 1)
    if (n.core > 1) {
      it = split(1:m, factor(cut(1:m, n.core, labels = FALSE)))
      resit = parallel::mclapply(it, function(markers) { apply(X[, markers], 2, mean) /2 }, mc.cores = n.core)
      p = unlist(resit)
    } else {
      p = apply(X, 2, mean) / 2
    }
    q = 1 - p
    var.A = 2 * mean(p * q)
    Mp = tcrossprod(v1, matrix(p, m, 1))
    W = X - 2 * Mp
    A = tcrossprod(W) / var.A / m
    return(A)
  }
  make.full = function(X) {
    svd.X = svd(X)
    r = max(which(svd.X$d > 1e-08))
    return(as.matrix(svd.X$u[, 1:r]))
  }
  calc.LL = function(ZX) {
    I = ncol(ZX)
    Add = LL1 = LL2 = rep(NA, I)
    Wald = Eff = matrix(NA, nrow = I, ncol = L)
    colnames(Eff) = paste0("effect@", Envs)
    colnames(Wald) = paste0("wald@", Envs)
    for (i in 1:I) {
      Xin = cbind(ZX[, i, drop = FALSE], FIXED)
      soln = gaston::lmm.diago(Y = y,
                               X = Xin,
                               eigenK = eiK,
                               verbose = FALSE)
      Add[i] = soln$BLUP_beta[1] / 2
      LL1[i] = gaston::lmm.diago.profile.likelihood(tau = soln$tau,
                                                    s2 = soln$sigma2,
                                                    Y = y,
                                                    X = Xin,
                                                    eigenK = eiK)[1, 1]
      pxi = Ze
      xi = ZX[, i, drop = FALSE]
      for (j in 1:ncol(pxi)) pxi[, j] = xi * pxi[, j]
      Xin = cbind(pxi, FIXED)
      soln = gaston::lmm.diago(Y = y,
                               X = Xin,
                               eigenK = eiK,
                               verbose = FALSE)
      Eff[i, ] = soln$BLUP_beta[1:L]
      Wald[i, ] = (soln$BLUP_beta[1:L] / sqrt(diag(soln$varbeta)[1:L]))^2
      LL2[i] = gaston::lmm.diago.profile.likelihood(tau = soln$tau,
                                                    s2 = soln$sigma2,
                                                    Y = y,
                                                    X = Xin,
                                                    eigenK = eiK)[1, 1]
    }
    P.ame = pchisq(2 * (LL1 - LL0), 1, lower.tail = FALSE)
    P.all = pchisq(2 * (LL2 - LL0), L, lower.tail = FALSE)
    LLmax = apply(cbind(LL1, LL0), 1, max, na.rm = TRUE)
    P.int = pchisq(2 * (LL2 - LLmax), L - 1, lower.tail = FALSE)
    for (i in 1:ncol(Wald)) Wald[, i] = -log10(pchisq(Wald[, i], 1, lower.tail = FALSE))
    res = data.frame(P.ame = -log10(P.ame),
                     P.all = -log10(P.all),
                     P.int = -log10(P.int),
                     Add = Add)
    res = cbind(res, Eff, Wald)
    return(res)
  }
  print("#################", quote = FALSE)
  print("#### gwasQxE ####", quote = FALSE)
  print("#################", quote = FALSE)
  print("## Trait", quote = FALSE)
  print(trait, quote = FALSE)
  print("## Environment", quote = FALSE)
  print(Env, quote = FALSE)
  X = t(geno[, -c(1:3)])
  colnames(X) = geno[, 1]
  pheno = pheno[!is.na(pheno[[trait]]), ]
  uid = unique(pheno[, 1])
  print("## following samples do not have phenotype", quote = FALSE)
  no.pheno = rownames(X)[!is.element(rownames(X), uid)]
  if (length(no.pheno) != 0) print(no.pheno, quote = FALSE)
  X = X[is.element(rownames(X), uid), ]
  uid = unique(rownames(X))
  print("## following samples do not have genotype", quote = FALSE)
  no.geno = unique(pheno[!is.element(pheno[, 1], uid), 1])
  if (length(no.geno) != 0) print(no.geno, quote = FALSE)
  pheno = pheno[is.element(pheno[, 1], uid), ]
  FIXED = matrix(1, nrow = nrow(pheno), ncol = 1)
  if (!is.null(Fixed)) FIXED = cbind(FIXED, pheno[, is.element(colnames(pheno), Fixed)])
  MAF = apply(X + 1, 2, function(x) { mean(x, na.rm = TRUE) }) / 2
  MAF = apply(rbind(MAF, 1 - MAF), 2, min)
  G = myAmat(X)
  if (n.PC > 0) S = eigen(G)$vectors[, 1:n.PC, drop = FALSE]
  Z = makeZ(pheno[, 1], rownames(X))
  ZX = Z %*% X
  ZGZ = Z %*% G %*% t(Z)
  if (n.PC > 0) FIXED = cbind(FIXED, Z %*% S)
  y = pheno[[trait]]
  trials = unique(pheno$trial)
  if (scale.pheno.by.trial) {
    for (i in trials) y[pheno$trial == i] = scale(y[pheno$trial == i])
  }
  pheno$trial = as.numeric(as.factor(pheno$trial))
  Tc = makeZ(pheno$trial, unique(pheno$trial))
  ZGtZ = tcrossprod(Tc) * ZGZ
  FIXED = cbind(FIXED, Tc)
  FIXED = make.full(FIXED)
  Envs = unique(pheno[[Env]])
  Ze = makeZ(pheno[[Env]], Envs)
  if (Env != "trial") {
    ZGeZ = tcrossprod(Ze) * ZGZ
    soln = gaston::lmm.aireml(Y = y, X = FIXED, K = list(ZGZ, ZGeZ, ZGtZ), verbose = FALSE)
    print("## variance components", quote = FALSE)
    print(paste("sigma_G:", soln$tau[1], "sigma_Ge:", soln$tau[2], "sigma_Gt:", soln$tau[3]), quote = FALSE)
    K = (soln$tau[1] * ZGZ) + (soln$tau[2] * ZGeZ) + (soln$tau[3] * ZGtZ)
    eiK = eigen(K)
  } else {
    soln = gaston::lmm.aireml(Y = y, X = FIXED, K = list(ZGZ, ZGtZ), verbose = FALSE)
    print("## variance components", quote = FALSE)
    print(paste("sigma_G:", soln$tau[1], "sigma_Gt", soln$tau[2]), quote = FALSE)
    K = (soln$tau[1]) * ZGZ + (soln$tau[2] * ZGtZ)
    eiK = eigen(K)
  }
  n = nrow(ZX)
  L = length(Envs)
  ZX = ZX + 1
  soln = gaston::lmm.diago(Y = y,
                           X = FIXED,
                           eigenK = eiK,
                           verbose = FALSE)
  LL0 = gaston::lmm.diago.profile.likelihood(tau = soln$tau,
                                             s2 = soln$sigma2,
                                             Y = y,
                                             X = FIXED,
                                             eigenK = eiK)[1, 1]
  if (n.core > 1) {
    m = ncol(ZX)
    it = split(1:m, factor(cut(1:m, n.core, labels = FALSE)))
    resit = parallel::mclapply(it, function(markers) { calc.LL(ZX[, markers]) }, mc.cores = n.core)
    res = as.data.frame(matrix(NA, nrow = m, ncol = ncol(resit[[1]])))
    colnames(res) = colnames(resit[[1]])
    for (ii in 1:n.core) res[unlist(it[[ii]]), ] = resit[[ii]]
    res = cbind(geno[,1:3], res, MAF)
  } else {
    res = calc.LL(ZX)
    res = cbind(geno[,1:3], res, MAF)
  }
  return(res)
}
