gwasQxE_all = function(geno,
                       pheno,
                       trait,
                       Env,
                       n.PC,
                       Fixed,
                       n.core,
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
    LL1_G = LL1_GT = LL1_GE = LL1_GTE = LL1_GTd = rep(NA, I)
    LLami_G = LLami_GT = LLami_GE = LLami_GTE = LLami_GTd = rep(NA, I)
    LLgge_G = LLgge_GT = LLgge_GE = LLgge_GTE = LLgge_GTd = rep(NA, I)
    for (i in 1:I) {
      Xin = cbind(ZX[, i, drop = FALSE], FIXED)
      soln = gaston::lmm.diago(Y = y,
                               X = Xin,
                               eigenK = eiG,
                               verbose = FALSE)
      LL1_G[i] = gaston::lmm.diago.profile.likelihood(tau = soln$tau,
                                                      s2 = soln$sigma2,
                                                      Y = y,
                                                      X = Xin,
                                                      eigenK = eiG)[1, 1]
      soln = gaston::lmm.diago(Y = y,
                               X = Xin,
                               eigenK = eiGE,
                               verbose = FALSE)
      LL1_GE[i] = gaston::lmm.diago.profile.likelihood(tau = soln$tau,
                                                       s2 = soln$sigma2,
                                                       Y = y,
                                                       X = Xin,
                                                       eigenK = eiGE)[1, 1]
      soln = gaston::lmm.diago(Y = y,
                               X = Xin,
                               eigenK = eiGT,
                               verbose = FALSE)
      LL1_GT[i] = gaston::lmm.diago.profile.likelihood(tau = soln$tau,
                                                       s2 = soln$sigma2,
                                                       Y = y,
                                                       X = Xin,
                                                       eigenK = eiGT)[1, 1]
      soln = gaston::lmm.diago(Y = y,
                               X = Xin,
                               eigenK = eiGTE,
                               verbose = FALSE)
      LL1_GTE[i] = gaston::lmm.diago.profile.likelihood(tau = soln$tau,
                                                        s2 = soln$sigma2,
                                                        Y = y,
                                                        X = Xin,
                                                        eigenK = eiGTE)[1, 1]

      soln = gaston::lmm.diago(Y = y,
                               X = Xin,
                               eigenK = eiGTd,
                               verbose = FALSE)
      LL1_GTd[i] = gaston::lmm.diago.profile.likelihood(tau = soln$tau,
                                                        s2 = soln$sigma2,
                                                        Y = y,
                                                        X = Xin,
                                                        eigenK = eiGTd)[1, 1]
      pxi = Ze
      xi = ZX[, i, drop = FALSE]
      for (j in 1:ncol(pxi)) pxi[, j] = xi * pxi[, j]
      Xin = cbind(xi, pxi, FIXED)
      soln = gaston::lmm.diago(Y = y,
                               X = Xin,
                               eigenK = eiG,
                               verbose = FALSE)
      LLami_G[i] = gaston::lmm.diago.profile.likelihood(tau = soln$tau,
                                                           s2 = soln$sigma2,
                                                           Y = y,
                                                           X = Xin,
                                                           eigenK = eiG)[1, 1]
      soln = gaston::lmm.diago(Y = y,
                               X = Xin,
                               eigenK = eiGE,
                               verbose = FALSE)
      LLami_GE[i] = gaston::lmm.diago.profile.likelihood(tau = soln$tau,
                                                            s2 = soln$sigma2,
                                                            Y = y,
                                                            X = Xin,
                                                            eigenK = eiGE)[1, 1]
      soln = gaston::lmm.diago(Y = y,
                               X = Xin,
                               eigenK = eiGT,
                               verbose = FALSE)
      LLami_GT[i] = gaston::lmm.diago.profile.likelihood(tau = soln$tau,
                                                            s2 = soln$sigma2,
                                                            Y = y,
                                                            X = Xin,
                                                            eigenK = eiGT)[1, 1]
      soln = gaston::lmm.diago(Y = y,
                               X = Xin,
                               eigenK = eiGTE,
                               verbose = FALSE)
      LLami_GTE[i] = gaston::lmm.diago.profile.likelihood(tau = soln$tau,
                                                             s2 = soln$sigma2,
                                                             Y = y,
                                                             X = Xin,
                                                             eigenK = eiGTE)[1, 1]

      soln = gaston::lmm.diago(Y = y,
                               X = Xin,
                               eigenK = eiGTd,
                               verbose = FALSE)
      LLami_GTd[i] = gaston::lmm.diago.profile.likelihood(tau = soln$tau,
                                                             s2 = soln$sigma2,
                                                             Y = y,
                                                             X = Xin,
                                                             eigenK = eiGTd)[1, 1]
      pxi = Ze
      xi = ZX[, i, drop = FALSE]
      for (j in 1:ncol(pxi)) pxi[, j] = xi * pxi[, j]
      Xin = cbind(pxi, FIXED)
      soln = gaston::lmm.diago(Y = y,
                               X = Xin,
                               eigenK = eiG,
                               verbose = FALSE)
      LLgge_G[i] = gaston::lmm.diago.profile.likelihood(tau = soln$tau,
                                                           s2 = soln$sigma2,
                                                           Y = y,
                                                           X = Xin,
                                                           eigenK = eiG)[1, 1]
      soln = gaston::lmm.diago(Y = y,
                               X = Xin,
                               eigenK = eiGE,
                               verbose = FALSE)
      LLgge_GE[i] = gaston::lmm.diago.profile.likelihood(tau = soln$tau,
                                                            s2 = soln$sigma2,
                                                            Y = y,
                                                            X = Xin,
                                                            eigenK = eiGE)[1, 1]
      soln = gaston::lmm.diago(Y = y,
                               X = Xin,
                               eigenK = eiGT,
                               verbose = FALSE)
      LLgge_GT[i] = gaston::lmm.diago.profile.likelihood(tau = soln$tau,
                                                            s2 = soln$sigma2,
                                                            Y = y,
                                                            X = Xin,
                                                            eigenK = eiGT)[1, 1]
      soln = gaston::lmm.diago(Y = y,
                               X = Xin,
                               eigenK = eiGTE,
                               verbose = FALSE)
      LLgge_GTE[i] = gaston::lmm.diago.profile.likelihood(tau = soln$tau,
                                                             s2 = soln$sigma2,
                                                             Y = y,
                                                             X = Xin,
                                                             eigenK = eiGTE)[1, 1]
      soln = gaston::lmm.diago(Y = y,
                               X = Xin,
                               eigenK = eiGTd,
                               verbose = FALSE)
      LLgge_GTd[i] = gaston::lmm.diago.profile.likelihood(tau = soln$tau,
                                                             s2 = soln$sigma2,
                                                             Y = y,
                                                             X = Xin,
                                                             eigenK = eiGTd)[1, 1]
    }

    AME_G = -log10(pchisq(2 * (LL1_G - LL0_G), 1, lower.tail = FALSE))
    AME_GE = -log10(pchisq(2 * (LL1_GE - LL0_GE), 1, lower.tail = FALSE))
    AME_GT = -log10(pchisq(2 * (LL1_GT - LL0_GT), 1, lower.tail = FALSE))
    AME_GTE = -log10(pchisq(2 * (LL1_GTE - LL0_GTE), 1, lower.tail = FALSE))
    AME_GTd = -log10(pchisq(2 * (LL1_GTd - LL0_GTd), 1, lower.tail = FALSE))

    AMI_G = -log10(pchisq(2 * (LLami_G - LL0_G), L + 1, lower.tail = FALSE))
    AMI_GE = -log10(pchisq(2 * (LLami_GE - LL0_GE), L + 1, lower.tail = FALSE))
    AMI_GT = -log10(pchisq(2 * (LLami_GT - LL0_GT), L + 1, lower.tail = FALSE))
    AMI_GTE = -log10(pchisq(2 * (LLami_GTE - LL0_GTE), L + 1, lower.tail = FALSE))
    AMI_GTd = -log10(pchisq(2 * (LLami_GTd - LL0_GTd), L + 1, lower.tail = FALSE))

    LLL_G = apply(cbind(LL1_G, LL0_G), 1, max, na.rm = TRUE)
    AMIqxe_G = -log10(pchisq(2 * (LLami_G - LLL_G), L, lower.tail = FALSE))
    LLL_GE = apply(cbind(LL1_GE, LL0_GE), 1, max, na.rm = TRUE)
    AMIqxe_GE = -log10(pchisq(2 * (LLami_GE - LLL_GE), L, lower.tail = FALSE))
    LLL_GT = apply(cbind(LL1_GT, LL0_GT), 1, max, na.rm = TRUE)
    AMIqxe_GT = -log10(pchisq(2 * (LLami_GT - LLL_GT), L, lower.tail = FALSE))
    LLL_GTE = apply(cbind(LL1_GTE, LL0_GTE), 1, max, na.rm = TRUE)
    AMIqxe_GTE = -log10(pchisq(2 * (LLami_GTE - LLL_GTE), L, lower.tail = FALSE))
    LLL_GTd = apply(cbind(LL1_GTd, LL0_GTd), 1, max, na.rm = TRUE)
    AMIqxe_GTd = -log10(pchisq(2 * (LLami_GTd - LLL_GTd), L, lower.tail = FALSE))

    GGE_G = -log10(pchisq(2 * (LLgge_G - LL0_G), L, lower.tail = FALSE))
    GGE_GE = -log10(pchisq(2 * (LLgge_GE - LL0_GE), L, lower.tail = FALSE))
    GGE_GT = -log10(pchisq(2 * (LLgge_GT - LL0_GT), L, lower.tail = FALSE))
    GGE_GTE = -log10(pchisq(2 * (LLgge_GTE - LL0_GTE), L, lower.tail = FALSE))
    GGE_GTd = -log10(pchisq(2 * (LLgge_GTd - LL0_GTd), L, lower.tail = FALSE))

    LLL_G = apply(cbind(LL1_G, LL0_G), 1, max, na.rm = TRUE)
    GGEqxe_G = -log10(pchisq(2 * (LLgge_G - LLL_G), L - 1, lower.tail = FALSE))
    LLL_GE = apply(cbind(LL1_GE, LL0_GE), 1, max, na.rm = TRUE)
    GGEqxe_GE = -log10(pchisq(2 * (LLgge_GE - LLL_GE), L - 1, lower.tail = FALSE))
    LLL_GT = apply(cbind(LL1_GT, LL0_GT), 1, max, na.rm = TRUE)
    GGEqxe_GT = -log10(pchisq(2 * (LLgge_GT - LLL_GT), L - 1, lower.tail = FALSE))
    LLL_GTE = apply(cbind(LL1_GTE, LL0_GTE), 1, max, na.rm = TRUE)
    GGEqxe_GTE = -log10(pchisq(2 * (LLgge_GTE - LLL_GTE), L - 1, lower.tail = FALSE))
    LLL_GTd = apply(cbind(LL1_GTd, LL0_GTd), 1, max, na.rm = TRUE)
    GGEqxe_GTd = -log10(pchisq(2 * (LLgge_GTd - LLL_GTd), L - 1, lower.tail = FALSE))

    res = data.frame(AME_G = AME_G,
                     AME_GE = AME_GE,
                     AME_GT = AME_GT,
                     AME_GTE = AME_GTE,
                     AME_GTd = AME_GTd,

                     AMI_G = AMI_G,
                     AMI_GE = AMI_GE,
                     AMI_GT = AMI_GT,
                     AMI_GTE = AMI_GTE,
                     AMI_GTd = AMI_GTd,

                     AMIqxe_G = AMIqxe_G,
                     AMIqxe_GE = AMIqxe_GE,
                     AMIqxe_GT = AMIqxe_GT,
                     AMIqxe_GTE = AMIqxe_GTE,
                     AMIqxe_GTd = AMIqxe_GTd,

                     GGE_G = GGE_G,
                     GGE_GE = GGE_GE,
                     GGE_GT = GGE_GT,
                     GGE_GTE = GGE_GTE,
                     GGE_GTd = GGE_GTd,

                     GGEqxe_G = GGEqxe_G,
                     GGEqxe_GE = GGEqxe_GE,
                     GGEqxe_GT = GGEqxe_GT,
                     GGEqxe_GTE = GGEqxe_GTE,
                     GGEqxe_GTd = GGEqxe_GTd)

    return(res)
  }
  print("#####################", quote = FALSE)
  print("#### gwasQxE_all ####", quote = FALSE)
  print("#####################", quote = FALSE)
  print("## Trait", quote = FALSE)
  print(trait, quote = FALSE)
  print("## Environment", quote = FALSE)
  print(Env, quote = FALSE)
  X = t(geno[,-c(1:3)])
  colnames(X) = geno$marker
  pheno = pheno[!is.na(pheno[[trait]]), ]
  uid = unique(pheno[,1])
  print("## following samples do not have phenotype", quote = FALSE)
  no.pheno = rownames(X)[!is.element(rownames(X), uid)]
  if (length(no.pheno) != 0) print(no.pheno, quote = FALSE)
  X = X[is.element(rownames(X), uid), ]
  uid = unique(rownames(X))
  print("## following samples do not have genotype", quote = FALSE)
  no.geno = unique(pheno[!is.element(pheno[,1], uid), 1])
  if (length(no.geno) != 0) print(no.geno, quote = FALSE)
  pheno = pheno[is.element(pheno[,1], uid), ]
  FIXED = matrix(1, nrow = nrow(pheno), ncol = 1)
  if (!is.null(Fixed)) FIXED = cbind(FIXED, pheno[, is.element(colnames(pheno), Fixed)])
  G = myAmat(X)
  if (n.PC > 0) S = eigen(G)$vectors[, 1:n.PC, drop = FALSE]
  Z = makeZ(pheno[,1], rownames(X))
  ZX = Z %*% X
  G = Z %*% G %*% t(Z)
  if (n.PC > 0) FIXED = cbind(FIXED, Z %*% S)
  y = pheno[[trait]]
  trials = unique(pheno$trial)
  if (scale.pheno.by.trial) {
    for (i in trials) y[pheno$trial == i] = scale(y[pheno$trial == i])
  }
  pheno$trial = as.numeric(as.factor(pheno$trial))
  Tc = makeZ(pheno$trial, unique(pheno$trial))
  GT = tcrossprod(Tc) * G
  FIXED = cbind(FIXED, Tc)
  Envs = unique(pheno[[Env]])
  Ze = makeZ(pheno[[Env]], Envs)
  GE = tcrossprod(Ze) * G
  FIXED = make.full(FIXED)
  eiG = eigen(G)
  soln = gaston::lmm.aireml(Y = y, X = FIXED, K = list(G, GE), verbose = FALSE)
  GGE = (soln$tau[1] * G) + (soln$tau[2] * GE)
  eiGE = eigen(GGE)
  soln = gaston::lmm.aireml(Y = y, X = FIXED, K = list(G, GT), verbose = FALSE)
  GGT = (soln$tau[1] * G) + (soln$tau[2] * GT)
  eiGT = eigen(GGT)
  soln = gaston::lmm.aireml(Y = y, X = FIXED, K = list(G, GE, GT), verbose = FALSE)
  GGTE = (soln$tau[1] * G) + (soln$tau[2] * GE) + (soln$tau[3] * GT)
  eiGTE = eigen(GGTE)
  K.list = vector(mode = "list", length = 1 + ncol(Tc))
  K.list[[1]] = G
  for (i in 1:ncol(Tc)) {
    GTd = GT
    GTd[colnames(GT) != i, colnames(GT) != i] = 0
    K.list[[(i + 1)]] = GTd
  }
  soln = gaston::lmm.aireml(Y = y, X = FIXED, K = K.list, verbose = FALSE)
  GTd = soln$tau[1] * K.list[[1]]
  for (i in 2:length(K.list))   GGTd = GTd + soln$tau[i] * K.list[[i]]
  eiGTd = eigen(GGTd)
  n = nrow(ZX)
  L = length(Envs)
  ZX = ZX + 1
  soln = gaston::lmm.diago(Y = y,
                           X = FIXED,
                           eigenK = eiG,
                           verbose = FALSE)
  LL0_G = gaston::lmm.diago.profile.likelihood(tau = soln$tau,
                                               s2 = soln$sigma2,
                                               Y = y,
                                               X = FIXED,
                                               eigenK = eiG)[1, 1]
  soln = gaston::lmm.diago(Y = y,
                           X = FIXED,
                           eigenK = eiGE,
                           verbose = FALSE)
  LL0_GE = gaston::lmm.diago.profile.likelihood(tau = soln$tau,
                                                s2 = soln$sigma2,
                                                Y = y,
                                                X = FIXED,
                                                eigenK = eiGE)[1, 1]
  soln = gaston::lmm.diago(Y = y,
                           X = FIXED,
                           eigenK = eiGT,
                           verbose = FALSE)
  LL0_GT = gaston::lmm.diago.profile.likelihood(tau = soln$tau,
                                                s2 = soln$sigma2,
                                                Y = y,
                                                X = FIXED,
                                                eigenK = eiGT)[1, 1]
  soln = gaston::lmm.diago(Y = y,
                           X = FIXED,
                           eigenK = eiGTE,
                           verbose = FALSE)
  LL0_GTE = gaston::lmm.diago.profile.likelihood(tau = soln$tau,
                                                 s2 = soln$sigma2,
                                                 Y = y,
                                                 X = FIXED,
                                                 eigenK = eiGTE)[1, 1]
  soln = gaston::lmm.diago(Y = y,
                           X = FIXED,
                           eigenK = eiGTd,
                           verbose = FALSE)
  LL0_GTd = gaston::lmm.diago.profile.likelihood(tau = soln$tau,
                                                 s2 = soln$sigma2,
                                                 Y = y,
                                                 X = FIXED,
                                                 eigenK = eiGTd)[1, 1]
  if (n.core > 1) {
    m = ncol(ZX)
    it = split(1:m, factor(cut(1:m, n.core, labels = FALSE)))
    resit = parallel::mclapply(it, function(markers) { calc.LL(ZX[, markers]) }, mc.cores = n.core)
    res = as.data.frame(matrix(NA, nrow = m, ncol = ncol(resit[[1]])))
    colnames(res) = colnames(resit[[1]])
    for (ii in 1:n.core) res[unlist(it[[ii]]), ] = resit[[ii]]
    res = cbind(geno[,1:3], res)
  } else {
    res = calc.LL(ZX)
    res = cbind(geno[,1:3], res)
  }
  return(res)
}
