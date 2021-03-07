simPhenoQxE = function(geno,
                       n.qtl,
                       seed,
                       q.Effs,
                       Env,
                       w_ratio,
                       PVE.qtl,
                       h2,
                       domEff = FALSE)
{
  makeZ = function(row.id, col.id) {
    Z = matrix(0, nrow=length(row.id), ncol=length(col.id))
    rownames(Z) = row.id
    colnames(Z) = col.id
    Z[cbind(1:nrow(Z),match(row.id,col.id))] = 1
    return(Z)
  }
  adj.omega = function(omega, XB, u, PVE, h2)	{
    expected = PVE/h2
    observed = var(XB[,1]) / var(apply(XB, 1, sum) + (u * omega))
    dif = abs(expected - observed)
    return(dif)
  }
  adj.psi = function(psi, g_values, epsilon, h2)	{
    observed = var(g_values) / var(g_values + (epsilon * psi))
    dif = abs(h2 - observed)
    return(dif)
  }
  set.seed(seed)
  chroms = unique(geno$chrom)
  q.chr = sample(chroms, n.qtl)
  simQTLs = c()
  for (i in 1:n.qtl) {
    mrks = geno$marker[geno$chrom == q.chr[i]]
    QTLi = sample(mrks, 1)
    simQTLs = c(simQTLs, QTLi)
  }
  X = t(geno[,-c(1:3)])
  colnames(X) = geno$marker
  Q = as.matrix(X[, is.element(geno$marker, simQTLs), drop = FALSE])
  Qq = c()
  for (i in 1:ncol(Q)) {
    xi = Q[, i, drop = FALSE]
    if (domEff) xi[xi == 0] = 2
    xi = as.numeric(xi %*% t(q.Effs))
    xi = as.numeric(scale(xi) / (1 / sqrt(PVE.qtl)))
    Qq = cbind(Qq, xi)
  }
  colnames(Qq) = simQTLs
  n.trials = length(q.Effs)
  n = nrow(X)
  G = tcrossprod(X) / (ncol(X) - 1)
  Z = rep(rownames(G), n.trials)
  Z = makeZ(Z, rownames(G))
  ZGZ = Z %*% G %*% t(Z)
  Tc = c()
  for (i in 1:n.trials) Tc = c(Tc, rep(i, n))
  Zt = makeZ(Tc, 1:n.trials)
  ZTZ = tcrossprod(Zt)
  ZGtZ = ZGZ * ZTZ
  Ec = c()
  for (i in 1:n.trials) Ec = c(Ec, rep(Env[i], n))
  Ze = makeZ(Ec, unique(Ec))
  ZEZ = tcrossprod(Ze)
  ZGeZ = ZGZ * ZEZ
  mu = numeric(nrow(ZGZ))
  u_g = scale(MASS::mvrnorm(n = 1, mu = mu, Sigma = ZGZ))
  u_gt = scale(MASS::mvrnorm(n = 1, mu = mu, Sigma = ZGtZ))
  u_ge = scale(MASS::mvrnorm(n = 1, mu = mu, Sigma = ZGeZ))
  u = (w_ratio[1] * u_g) + (w_ratio[2] * u_ge) + (w_ratio[3] * u_gt)
  u = as.numeric(scale(u))
  omega = optimize(adj.omega, interval = c(1e-06, 2), Qq, u, PVE.qtl, h2)$minimum
  g_values = apply(Qq, 1, sum) + (u * omega)
  epsilon = as.numeric(scale(rnorm(length(g_values), sd = 1)))
  psi = optimize(adj.psi, interval = c(1e-06, 2), g_values, epsilon, h2)$minimum
  epsilon = epsilon * psi
  y = g_values + epsilon
  PVE.each = as.data.frame(matrix(NA, nrow = n.trials, ncol = 2))
  colnames(PVE.each) = c("PVE.qtl", "h2")
  for (i in 1:n.trials) PVE.each$PVE.qtl[i] = var(Qq[(n*(i-1)+1):(n*i),1]) / var(y[(n*(i-1)+1):(n*i)])
  for (i in 1:n.trials) PVE.each$h2[i] = var(g_values[(n*(i-1)+1):(n*i)]) / var(y[(n*(i-1)+1):(n*i)])
  realized_PVE.qtl = var(Qq[,1]) / var(y)
  realized_h2 = var(g_values) / var(y)
  print(c(paste("PVE.qtl =", realized_PVE.qtl), c(paste("h2 =", realized_h2))))
  pheno = data.frame(line_name = rep(rownames(X), n.trials),
                   env = Ec,
                   trial = Tc,
                   y = y)
  ret = list(pheno = pheno,
             simQTLs = simQTLs,
             realized_PVE.qtl = realized_PVE.qtl,
             realized_h2 = realized_h2,
             PVE.each = PVE.each)
  return(ret)
}

