# 7-Mar-2021
# This R script was written for 
# Exploring Efficient Linear Mixed Models to Detect Quantitative Trait Locus-by-environment Interactions
# correspondence: yame.repos@gmail.com


#####################################################
### comparison of LMMs using simulated phenotypes ###
#####################################################

### libraries ####
if (!is.element("gaston", rownames(installed.packages()))) install.packages("gaston")
if (!is.element("gwasQxE", rownames(installed.packages()))) install.packages("gwasQxE_0.1.0.tar.gz", repos = NULL, type = "source")
require(gaston)
require(gwasQxE)
require(pROC)

### genotype file ####
geno = read.csv("SL96_geno_QxE.csv")

### simulation parameters ####
n.rep = 100; h2 = 0.5; n.qtl = 3; PVE.qtl = 0.10

### simulations ####
dir.name = "simResults"
system(paste0("mkdir ", dir.name))
for (condition in 1:8) {
  if (condition == 1) { Name = "Persistence" ; Env = c(1, 1, 2, 2) ; q.Effs = c(1, 1, 1, 1) ; w_ratio = c(1, 1, 1) ; domEff = FALSE ; method = "add" }
  if (condition == 2) { Name = "Divergence" ; Env = c(1, 1, 2, 2) ; q.Effs = c(0.5, 0.5, 1, 1) ; w_ratio = c(1, 1, 1) ; domEff = FALSE ; method = "add" }
  if (condition == 3) { Name = "Convergence" ; Env = c(1, 1, 2, 2) ; q.Effs = c(0, 0, 1, 1) ; w_ratio = c(1, 1, 1) ; domEff = FALSE ; method = "add" }
  if (condition == 4) { Name = "Crossover" ; Env = c(1, 1, 2, 2) ; q.Effs = c(-1, -1, 1, 1) ; w_ratio = c(1, 1, 1) ; domEff = FALSE ; method = "add" }
  if (condition == 5) { Name = "Convergence.axa" ; Env = c(1, 1, 2, 2) ; q.Effs = c(0, 0, 1, 1) ; w_ratio = c(1, 1, 1); domEff = FALSE ; method = "add" }
  if (condition == 6) { Name = "Convergence.axd" ; Env = c(1, 1, 2, 2) ; q.Effs = c(0, 0, 1, 1) ; w_ratio = c(1, 1, 1); domEff = FALSE ; method = "dom" }
  if (condition == 7) { Name = "Convergence.dxa" ; Env = c(1, 1, 2, 2) ; q.Effs = c(0, 0, 1, 1) ; w_ratio = c(1, 1, 1); domEff = TRUE ; method = "add" }
  if (condition == 8) { Name = "Convergence.dxd" ; Env = c(1, 1, 2, 2) ; q.Effs = c(0, 0, 1, 1) ; w_ratio = c(1, 1, 1); domEff = TRUE ; method = "dom" }
  
  ### GWAS score matrices
  M = matrix(NA, nrow = nrow(geno), ncol = n.rep)
  rownames(M) = geno$marker
  AME_G = AME_GE = AME_GT = AME_GTE = AME_GTd = M
  AMI_G = AMI_GE = AMI_GT = AMI_GTE = AMI_GTd = M
  AMIqxe_G = AMIqxe_GE = AMIqxe_GT = AMIqxe_GTE = AMIqxe_GTd = M
  GGE_G = GGE_GE = GGE_GT = GGE_GTE = GGE_GTd = M
  GGEqxe_G = GGEqxe_GE = GGEqxe_GT = GGEqxe_GTE = GGEqxe_GTd = M
  
  ### simulated QTL matrices
  simQTLs = matrix(NA, nrow = n.rep, ncol = n.qtl)
  
  ### simulations
  for (i in 1:n.rep) {
    seed = i
    check.sim = TRUE
    while(check.sim) {
      sim.pheno = simPhenoQxE(geno, n.qtl, seed, q.Effs, Env, w_ratio, PVE.qtl, h2, domEff)
      ch.PVE = abs(PVE.qtl - sim.pheno$realized_PVE.qtl)
      ch.h2 = abs(h2 - sim.pheno$realized_h2)
      check = ((ch.h2 < 1e-03) & (ch.PVE < 1e-03))
      if (check) check.sim = FALSE
      seed = seed + 999
    }
    simQTLs[i, ] = sim.pheno$simQTLs
    pheno = sim.pheno$pheno
    if (method == "add") res = gwasQxE_all(geno, pheno, trait = "y", Env = "env", n.PC = 2, Fixed = NULL, n.core = 100)
    if (method == "dom") res = gwasQxEd_all(geno, pheno, trait = "y", Env = "env", n.PC = 2, Fixed = NULL, n.core = 100)
    AME_G[,i] = res$AME_G
    AME_GE[,i] = res$AME_GE
    AME_GT[,i] = res$AME_GT
    AME_GTE[,i] = res$AME_GTE
    AME_GTd[,i] = res$AME_GTd
    AMI_G[,i] = res$AMI_G
    AMI_GE[,i] = res$AMI_GE
    AMI_GT[,i] = res$AMI_GT
    AMI_GTE[,i] = res$AMI_GTE
    AMI_GTd[,i] = res$AMI_GTd
    AMIqxe_G[,i] = res$AMIqxe_G
    AMIqxe_GE[,i] = res$AMIqxe_GE
    AMIqxe_GT[,i] = res$AMIqxe_GT
    AMIqxe_GTE[,i] = res$AMIqxe_GTE
    AMIqxe_GTd[,i] = res$AMIqxe_GTd
    GGE_G[,i] = res$GGE_G
    GGE_GE[,i] = res$GGE_GE
    GGE_GT[,i] = res$GGE_GT
    GGE_GTE[,i] = res$GGE_GTE
    GGE_GTd[,i] = res$GGE_GTd
    GGEqxe_G[,i] = res$GGEqxe_G
    GGEqxe_GE[,i] = res$GGEqxe_GE
    GGEqxe_GT[,i] = res$GGEqxe_GT
    GGEqxe_GTE[,i] = res$GGEqxe_GTE
    GGEqxe_GTd[,i] = res$GGEqxe_GTd
  }
  
  ### saving simulation results
  res = list(AME_G, AME_GE, AME_GT, AME_GTE, AME_GTd,
             AMI_G, AMI_GE, AMI_GT, AMI_GTE, AMI_GTd,
             AMIqxe_G, AMIqxe_GE, AMIqxe_GT, AMIqxe_GTE, AMIqxe_GTd,
             GGE_G, GGE_GE, GGE_GT, GGE_GTE, GGE_GTd,
             GGEqxe_G, GGEqxe_GE, GGEqxe_GT, GGEqxe_GTE, GGEqxe_GTd)
  names(res) = c("AME_G", "AME_GE", "AME_GT", "AME_GTE", "AME_GTd",
                 "AMI_G", "AMI_GE", "AMI_GT", "AMI_GTE", "AMI_GTd",
                 "AMIqxe_G", "AMIqxe_GE", "AMIqxe_GT", "AMIqxe_GTE", "AMIqxe_GTd",
                 "GGE_G", "GGE_GE", "GGE_GT", "GGE_GTE", "GGE_GTd",
                 "GGEqxe_G", "GGEqxe_GE", "GGEqxe_GT", "GGEqxe_GTE", "GGEqxe_GTd")
  result.file.name = paste0(dir.name, "/", Name, "_res.RData")
  save(res, file = result.file.name)
  result.file.name = paste0(dir.name, "/", Name, "_ans.RData")
  save(simQTLs, file = result.file.name)
}

### FDR function ####
calcFDR = function(vec, beta = 0.05) {
  vec = unique(vec)
  pvalues = 10^(-1 * vec)
  N.test = length(pvalues)
  smooth.df = 3
  lambda = seq(0.01, 0.99, 0.01)
  pi0 = rep(0, length(lambda))
  for (i in 1:length(lambda)) {
    pi0[i] = mean(pvalues >= lambda[i])/(1 - lambda[i])
  }
  if (length(pi0[is.na(pi0)]) > 1) {
    pi0 = 1
  } else {
    spi0 = smooth.spline(lambda, pi0, df = smooth.df)
    pi0 = predict(spi0, x = max(lambda))$y
  }
  if (pi0 < 0.5) pi0 = 1
  pi0 = min(pi0, 1)
  p.sorted = sort(pvalues)
  qvalues = rep(NA, N.test)
  for (i in 1:N.test) {
    qvalues[i] = (pi0 * N.test * p.sorted[i]) / i
  }
  qvalues[qvalues > 1] = 1
  q.sig = which.min(qvalues <= beta)
  lowest.sig = p.sorted[q.sig]
  highest.notsig = p.sorted[q.sig + 1]
  thr.st = mean(c(lowest.sig, highest.notsig))
  return(-log10(thr.st))
}

### data reshaping ####
simResults.reshape = function(vec) {
  vec = as.character(vec)
  res = matrix(NA, nrow = length(vec), ncol = 2)
  for (i in 1:length(vec)) res[i, ] = unlist(strsplit(vec[i], split = "_")) 
  colnames(res) = c("fixed", "random")
  return(res)
}

### recall, precision, F-measure, and AUC ####
file.names = dir(dir.name)
res.names = file.names[grep("_res", file.names)]
res.names = res.names[-grep("x", res.names)]
res.names = paste(dir.name, res.names, sep = "/")
for (i in res.names) {
  load(i)
  Name = unlist(strsplit(i, "_res"))[1]
  ans.names = paste0(Name, "_ans.RData")
  load(ans.names)
  Name = unlist(strsplit(Name, "/"))[2]
  n.rep = ncol(res[[1]])
  markers = rownames(res[[1]])
  M.sig = M.recall = M.precision = M.auc = matrix(0, nrow = n.rep, ncol = length(res))
  colnames(M.sig) = colnames(M.recall) = colnames(M.precision) = colnames(M.auc) = names(res)
  for (j in 1:length(res)) {
    M = res[[j]]
    ix = unlist(strsplit(names(res)[j], "_"))
    auc = numeric(n.rep)
    for (k in 1:n.rep) {
      if (length(unique(M[, k])) < 100) next
      auc[k] = suppressMessages(pROC::roc(is.element(markers, simQTLs[k, ]), M[, k])$auc)
      M.auc[, j] = auc
    }
    for (k in 1:n.rep) {
      Mk = M[, k]
      if (length(unique(Mk)) < 100) next
      threshold = calcFDR(Mk)
      true.signals = Mk[is.element(rownames(M), simQTLs[k,])]
      M.recall[k, j] = length(true.signals[true.signals >= threshold]) / ncol(simQTLs)
      M.precision[k, j] = length(true.signals[true.signals >= threshold]) / length(Mk[Mk >= threshold])
    }
  }
  
  mean.auc = apply(M.auc, 2, mean, na.rm = TRUE)
  sd.auc = apply(M.auc, 2, sd, na.rm = TRUE)
  mean.auc = reshape2::melt(mean.auc)
  sd.auc = reshape2::melt(sd.auc)
  indicator = rep("AUC", nrow(sd.auc))
  dat.auc = cbind(indicator, simResults.reshape(rownames(mean.auc)), mean.auc, sd.auc)
  colnames(dat.auc)[ncol(dat.auc)] = "sd"
  
  mean.recall = apply(M.recall, 2, mean, na.rm = TRUE)
  mean.recall[is.na(mean.recall)] = 0
  sd.rec = apply(M.recall, 2, sd, na.rm = TRUE)
  mean.recall = reshape2::melt(mean.recall)
  sd.rec = reshape2::melt(sd.rec)
  indicator = rep("Recall", nrow(sd.rec))
  dat.recall = cbind(indicator, simResults.reshape(rownames(mean.recall)), mean.recall, sd.rec)
  colnames(dat.recall)[ncol(dat.recall)] = "sd"
  
  mean.precision = apply(M.precision, 2, mean, na.rm = TRUE)
  mean.precision[is.na(mean.precision)] = 0
  sd.precision = apply(M.precision, 2, sd, na.rm = TRUE)
  mean.precision = reshape2::melt(mean.precision)
  sd.precision = reshape2::melt(sd.precision)
  indicator = rep("Precision", nrow(sd.precision))
  dat.precision = cbind(indicator, simResults.reshape(rownames(mean.precision)), mean.precision, sd.precision)
  colnames(dat.precision)[ncol(dat.precision)] = "sd"
  
  M.fmeasure = (2 * M.recall * M.precision) / (M.recall + M.precision) 
  mean.fmeasure = apply(M.fmeasure, 2, mean, na.rm = TRUE)
  mean.fmeasure[is.na(mean.fmeasure)] = 0
  sd.fmeasure = apply(M.fmeasure, 2, sd, na.rm = TRUE)
  mean.fmeasure = reshape2::melt(mean.fmeasure)
  sd.fmeasure = reshape2::melt(sd.fmeasure)
  indicator = rep("F-measure", nrow(sd.fmeasure))
  dat.fmeasure = cbind(indicator, simResults.reshape(rownames(mean.fmeasure)), mean.fmeasure, sd.fmeasure)
  colnames(dat.fmeasure)[ncol(dat.fmeasure)] = "sd"
  
  dat1 = rbind(dat.recall, dat.precision, dat.fmeasure, dat.auc)
  
  gp = ggplot2::ggplot(data = dat1, ggplot2::aes(x = factor(fixed, levels = c("AME", "AMI", "GGE", "AMIqxe", "GGEqxe")), y = value, fill = factor(random, levels = c("G", "GE", "GT", "GTE", "GTd")))) +
    ggplot2::facet_wrap( ~ factor(indicator, levels = c("Recall", "Precision",   "F-measure", "AUC")), ncol = 1, scales = "free_y") +
    ggplot2::scale_x_discrete(labels = c("Additive Main", "AMMI-type all", "GGE-type all", "AMMI-type int", "GGE-type int")) +
    ggplot2::geom_bar(stat="identity", color="black", position = ggplot2::position_dodge()) + 
    ggplot2::geom_errorbar(ggplot2::aes(ymin = value, ymax = value + sd), position = ggplot2::position_dodge(0.9), width = 0.2) + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12, angle = 45, hjust = 1, colour = "black")) +
    ggplot2::theme(axis.title.y = ggplot2::element_text(size = 14)) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 12, colour = "black")) +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 12, face = "bold")) +
    ggplot2::labs(x = "", y = "") + 
    ggplot2::theme(legend.position = 'none') + 
    ggplot2::theme(legend.text = ggplot2::element_text(size = 12)) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = "rge")) +
    ggplot2::scale_fill_manual(values = c("#FFFFFF", "#E0E0E0", "#BABABA", "#888888", "#000000"))
  gg.name = paste0(Name, "_power.png")
  ggplot2::ggsave(filename = gg.name, plot = gp, dpi = 1200, width = 4, height = 6)
}

### power to detect dominant effect QTL ####
file.names = dir(dir.name)
res.names = file.names[grep("_res", file.names)]
res.names = res.names[grep("x", res.names)]
res.names = paste(dir.name, res.names, sep = "/")
dat.auc = dat.recall = c()
for (i in res.names) {
  load(i)
  Name = unlist(strsplit(i, "_res"))[1]
  ans.names = paste0(Name, "_ans.RData")
  load(ans.names)
  Name = unlist(strsplit(Name, "/"))[2]
  QTL = unlist(strsplit(Name, "x"))[1]
  LMM = unlist(strsplit(Name, "x"))[2]
  n.rep = ncol(res[[1]])
  markers = rownames(res[[1]])
  M.sig = M.recall = M.precision = M.auc = matrix(0, nrow = n.rep, ncol = length(res))
  colnames(M.sig) = colnames(M.recall) = colnames(M.precision) = colnames(M.auc) = names(res)
  for (j in 1:length(res)) {
    M = res[[j]]
    ix = unlist(strsplit(names(res)[j], "_"))
    auc = numeric(n.rep)
    for (k in 1:n.rep) {
      if (length(unique(M[, k])) < 100) next
      auc[k] = suppressMessages(pROC::roc(is.element(markers, simQTLs[k, ]), M[, k])$auc)
      M.auc[, j] = auc
    }
    for (k in 1:n.rep) {
      Mk = M[, k]
      if (length(unique(Mk)) < 100) next
      threshold = calcFDR(Mk)
      true.signals = Mk[is.element(rownames(M), simQTLs[k,])]
      M.recall[k, j] = length(true.signals[true.signals >= threshold]) / ncol(simQTLs)
      M.precision[k, j] = length(true.signals[true.signals >= threshold]) / length(Mk[Mk >= threshold])
    }
  }
  
  mean.auc = apply(M.auc, 2, mean, na.rm = TRUE)
  sd.auc = apply(M.auc, 2, sd, na.rm = TRUE)
  mean.auc = reshape2::melt(mean.auc)
  sd.auc = reshape2::melt(sd.auc)
  indicator = rep("AUC", nrow(sd.auc))
  dat.auc.i = cbind(QTL, LMM, indicator, simResults.reshape(rownames(mean.auc)), mean.auc, sd.auc)
  colnames(dat.auc.i)[ncol(dat.auc.i)] = "sd"
  dat.auc = rbind(dat.auc, dat.auc.i)
  
  mean.recall = apply(M.recall, 2, mean, na.rm = TRUE)
  mean.recall[is.na(mean.recall)] = 0
  sd.rec = apply(M.recall, 2, sd, na.rm = TRUE)
  mean.recall = reshape2::melt(mean.recall)
  sd.rec = reshape2::melt(sd.rec)
  indicator = rep("Recall", nrow(sd.rec))
  dat.recall.i = cbind(QTL, LMM, indicator, simResults.reshape(rownames(mean.recall)), mean.recall, sd.rec)
  colnames(dat.recall.i)[ncol(dat.recall.i)] = "sd"
  dat.recall = rbind(dat.recall, dat.recall.i)
} 

dat.auc = dat.auc[dat.auc$random == "GTE", ]
gp = ggplot2::ggplot(data = dat.auc, ggplot2::aes(x = factor(fixed, levels = c("AME", "AMI", "GGE", "AMIqxe", "GGEqxe")), y = value, fill = LMM)) +
  ggplot2::facet_wrap( ~ QTL, ncol = 2, scales = "free_y") +
  ggplot2::scale_x_discrete(labels = c("Additive Main", "AMMI-type all", "GGE-type all", "AMMI-type int", "GGE-type int")) +
  ggplot2::geom_bar(stat="identity", color="black", position = ggplot2::position_dodge()) + 
  ggplot2::geom_errorbar(ggplot2::aes(ymin = value, ymax = value + sd), position = ggplot2::position_dodge(0.9), width = 0.2) + 
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12, angle = 45, hjust = 1, colour = "black")) +
  ggplot2::theme(axis.title.y = ggplot2::element_text(size = 14)) +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 12, colour = "black")) +
  ggplot2::theme(strip.text = ggplot2::element_text(size = 12, face = "bold")) +
  ggplot2::labs(x = "", y = "AUC") + 
  ggplot2::theme(legend.position = 'none') + 
  ggplot2::theme(legend.text = ggplot2::element_text(size = 12)) +
  ggplot2::guides(fill = ggplot2::guide_legend(title = "rge")) +
  ggplot2::scale_fill_manual(values = c("#FFFFFF", "#000000"))
gg.name = paste0(Name, "Convergence_add.vs.dom_AUC.png")
ggplot2::ggsave(filename = gg.name, plot = gp, dpi = 1200, width = 8, height = 2.6)

dat.recall = dat.recall[dat.recall$random == "GTE", ]
gp = ggplot2::ggplot(data = dat.recall, ggplot2::aes(x = factor(fixed, levels = c("AME", "AMI", "GGE", "AMIqxe", "GGEqxe")), y = value, fill = LMM)) +
  ggplot2::facet_wrap( ~ QTL, ncol = 2, scales = "free_y") +
  ggplot2::scale_x_discrete(labels = c("Additive Main", "AMMI-type all", "GGE-type all", "AMMI-type int", "GGE-type int")) +
  ggplot2::geom_bar(stat="identity", color="black", position = ggplot2::position_dodge()) + 
  ggplot2::geom_errorbar(ggplot2::aes(ymin = value, ymax = value + sd), position = ggplot2::position_dodge(0.9), width = 0.2) + 
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12, angle = 45, hjust = 1, colour = "black")) +
  ggplot2::theme(axis.title.y = ggplot2::element_text(size = 14)) +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 12, colour = "black")) +
  ggplot2::theme(strip.text = ggplot2::element_text(size = 12, face = "bold")) +
  ggplot2::labs(x = "", y = "Recall") + 
  ggplot2::theme(legend.position = 'none') + 
  ggplot2::theme(legend.text = ggplot2::element_text(size = 12)) +
  ggplot2::guides(fill = ggplot2::guide_legend(title = "rge")) +
  ggplot2::scale_fill_manual(values = c("#FFFFFF", "#000000"))
gg.name = paste0(Name, "Convergence_add.vs.dom_Recall.png")
ggplot2::ggsave(filename = gg.name, plot = gp, dpi = 1200, width = 8, height = 2.6)


gp = ggplot2::ggplot(data = dat1, ggplot2::aes(x = factor(FG, levels = c("AME", "AMI", "GGE", "AMIqxe", "GGEqxe")), y = value, fill = factor(mlmm, levels = c("a", "d")))) +
  ggplot2::facet_wrap( ~ mqtl, ncol = 2) +
  ggplot2::scale_x_discrete(labels = c("Additive Main", "AMMI-type QxE", "GGE-type QxE", "iAMMI-type QxE", "iGGE-type QxE")) +
  ggplot2::geom_bar(stat="identity", color="black", position = ggplot2::position_dodge()) + 
  ggplot2::geom_errorbar(ggplot2::aes(ymin = value, ymax = value + sd), position = ggplot2::position_dodge(0.9), width = 0.2) + 
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12, angle = 45, hjust = 1, colour = "black")) +
  ggplot2::theme(axis.title.y = ggplot2::element_text(size = 14)) +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 12, colour = "black")) +
  ggplot2::theme(strip.text = ggplot2::element_text(size = 12, face = "bold")) +
  ggplot2::labs(x = "", y = "mean.auc") + 
  ggplot2::theme(legend.position = 'none') + 
  ggplot2::theme(legend.text = ggplot2::element_text(size = 12)) +
  ggplot2::scale_fill_manual(values = c("#FFFFFF", "#000000"))
gg.name = paste0("gwasQxEd_all_auc.png")
ggplot2::ggsave(filename = gg.name, plot = gp, dpi = 1200, width = 8, height = 2.6)


#### followings are not yet

### genomic inflation factor ####
file.names = dir(dir.name)
res.names = file.names[grep("_res", file.names)]
res.names = res.names[-grep("x", res.names)]
res.names = paste(dir.name, res.names, sep = "/")
dat.lambda = c()
for (i in res.names) {
  load(i)
  Name = unlist(strsplit(i, "_res"))[1]
  Name = unlist(strsplit(Name, "/"))[2]
  n.rep = ncol(res[[1]])
  markers = rownames(res[[1]])
  M.lambda = matrix(NA, nrow = n.rep, ncol = length(res))
  colnames(M.lambda) = names(res)
  for (j in 1:length(res)) {
    M = res[[j]]
    lambda = numeric(n.rep)
    for (k in 1:n.rep) lambda[k] = median(log10(quantile(10^(-1 * M[, k]), 0.5, na.rm = TRUE)) / log10(0.5))
    M.lambda[, j] = lambda
  }
  Lam = reshape2::melt(M.lambda)
  dat.i = cbind(Name, simResults.reshape(Lam$Var2), Lam)
  dat.lambda = rbind(dat.lambda, dat.i)
}  
gp = ggplot2::ggplot(data = dat.lambda, ggplot2::aes(x = factor(fixed, levels = c("AME", "AMI", "GGE", "AMIqxe", "GGEqxe")), y = value, fill = factor(random, levels = c("G", "GE", "GT", "GTE", "GTd")))) +
  ggplot2::geom_boxplot(outlier.colour = NA) +
  ggplot2::scale_x_discrete(labels = c("Additive Main", "AMMI-type all", "GGE-type all", "AMMI-type int", "GGE-type int")) +
  ggplot2::geom_hline(yintercept = 1, color = "red") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12, angle = 45, hjust = 1, colour = "black")) +
  ggplot2::theme(axis.title.y = ggplot2::element_text(size = 14)) +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 12, colour = "black")) +
  ggplot2::theme(strip.text = ggplot2::element_text(size = 12, face = "bold")) +
  ggplot2::labs(x = "", y = expression(lambda[GC]), size = 14) + 
  ggplot2::ylim(0, 5) + 
  ggplot2::guides(fill = ggplot2::guide_legend(title = "rge")) +
  ggplot2::theme(legend.position = 'right') + 
  ggplot2::scale_fill_manual(values = c("#FFFFFF", "#E0E0E0", "#BABABA", "#888888", "#000000"))
ggplot2::ggsave(filename = "Genomic_inflation_factor.png", plot = gp, dpi = 1200, width = 7, height = 4)

QQ = matrix(0, ncol = 3, nrow = 1001)
for (i in res.names) {
  QQi = matrix(0, ncol = 3, nrow = 1001)
  load(i)
  Name = unlist(strsplit(i, "_res"))[1]
  Name = unlist(strsplit(Name, "/"))[2]
  n.rep = ncol(res[[1]])
  markers = rownames(res[[1]])
  ln = names(res)
  res = res[grep("_GTE$", names(res))]
  res = res[c(1, 2, 4)]
  for (j in 1:length(res)) {
    M = res[[j]]
    for (k in 1:ncol(M)) M[,k] = 10^(-1 * sort(M[,k]))
    qt.vec = seq(0, 1, by = 0.001)
    qt.mat = matrix(NA, nrow = length(qt.vec), ncol = ncol(M))
    for (k in 1:ncol(M)) {
      Mk = as.numeric(M[,k])
      Mk = Mk[Mk > 0]
      Mk = Mk[Mk < 1]
      qt.mat[,k] = quantile(Mk, qt.vec)
    }
    QQi[,j] = apply(qt.mat, 1, mean)
  }
  QQ = QQ + QQi
}
QQ = QQ / length(res.names)
jpeg(filename = "QQplot_GTE.random.jpeg", quality = 100, pointsize = 42, width = 960, height = 1080)
col.list = c("green4", "red", "blue")
ep = -log10(ppoints(nrow(QQ)))
op = -log10(QQ[,1])
plot(ep, op, ylim = c(0, 4.5), xlim = c(0, 4), type = "l", col = col.list[1], lwd = 0.5, xlab = expression(Expected ~ ~-log[10](italic(p))),
     ylab = expression(Observed ~ ~-log[10](italic(p))), axes = F)
lines(c(0, 4), c(0, 4), lty = 2, col = "black", lwd = 8)
for (i in 1:3) {
  par(new = TRUE)
  op = -log10(QQ[,i])
  plot(ep, op, ylim = c(0, 4.5), xlim = c(0, 4), type = "l", col = col.list[i], lwd = 8, xlab = "", ylab = "", axes = F)
}
axis(1, lwd = 5)
axis(2, lwd = 5)
dev.off()

QQ = matrix(0, ncol = 5, nrow = 1001)
for (i in res.names) {
  QQi = matrix(0, ncol = 5, nrow = 1001)
  load(i)
  Name = unlist(strsplit(i, "_res"))[1]
  Name = unlist(strsplit(Name, "/"))[2]
  n.rep = ncol(res[[1]])
  markers = rownames(res[[1]])
  ln = names(res)
  res = res[grep("GGE_", names(res))]
  for (j in 1:length(res)) {
    M = res[[j]]
    for (k in 1:ncol(M)) {
      if (length(unique(M[,k])) < 100) next
      M[,k] = 10^(-1 * sort(M[,k]))
    }
    qt.vec = seq(0, 1, by = 0.001)
    qt.mat = matrix(NA, nrow = length(qt.vec), ncol = ncol(M))
    for (k in 1:ncol(M)) {
      Mk = as.numeric(M[,k])
      if (length(unique(Mk)) < 100) next
      Mk = Mk[Mk > 0]
      Mk = Mk[Mk < 1]
      qt.mat[,k] = quantile(Mk, qt.vec)
    }
    QQi[,j] = apply(qt.mat, 1, mean, na.rm = TRUE)
  }
  QQ = QQ + QQi
}
QQ = QQ / length(res.names)
jpeg(filename = "QQplot_GGE.fixed.jpeg", quality = 100, pointsize = 42, width = 960, height = 1080)
col.list = c("green4", "red", "orange", "blue", "skyblue")
ep = -log10(ppoints(nrow(QQ)))
op = -log10(QQ[,1])
plot(ep, op, ylim = c(0, 4.5), xlim = c(0, 4), type = "l", col = col.list[1], lwd = 0.5, xlab = expression(Expected ~ ~-log[10](italic(p))),
     ylab = expression(Observed ~ ~-log[10](italic(p))), axes = F)
lines(c(0, 4), c(0, 4), lty = 2, col = "black", lwd = 8)
for (i in 1:5) {
  par(new = TRUE)
  op = -log10(QQ[,i])
  plot(ep, op, ylim = c(0, 4.5), xlim = c(0, 4), type = "l", col = col.list[i], lwd = 8, xlab = "", ylab = "", axes = F)
}
axis(1, lwd = 5)
axis(2, lwd = 5)
dev.off()

file.names = dir(dir.name)
res.names = file.names[grep("_res", file.names)]
res.names = res.names[grep("x", res.names)]
res.names = paste(dir.name, res.names, sep = "/")
dat.lambda = c()
for (i in res.names) {
  load(i)
  Name = unlist(strsplit(i, "_res"))[1]
  Name = unlist(strsplit(Name, "/"))[2]
  n.rep = ncol(res[[1]])
  markers = rownames(res[[1]])
  M.lambda = matrix(NA, nrow = n.rep, ncol = length(res))
  colnames(M.lambda) = names(res)
  for (j in 1:length(res)) {
    M = res[[j]]
    lambda = numeric(n.rep)
    for (k in 1:n.rep) lambda[k] = median(log10(quantile(10^(-1 * M[, k]), 0.5, na.rm = TRUE)) / log10(0.5))
    M.lambda[, j] = lambda
  }
  Lam = reshape2::melt(M.lambda)
  dat.i = cbind(Name, simResults.reshape(Lam$Var2), Lam)
  dat.lambda = rbind(dat.lambda, dat.i)
}  
gp = ggplot2::ggplot(data = dat.lambda, ggplot2::aes(x = factor(fixed, levels = c("AME", "AMI", "GGE", "AMIqxe", "GGEqxe")), y = value, fill = factor(random, levels = c("G", "GE", "GT", "GTE", "GTd")))) +
  ggplot2::geom_boxplot(outlier.colour = NA) +
  ggplot2::scale_x_discrete(labels = c("Additive Main", "AMMI-type all", "GGE-type all", "AMMI-type int", "GGE-type int")) +
  ggplot2::geom_hline(yintercept = 1, color = "red") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12, angle = 45, hjust = 1, colour = "black")) +
  ggplot2::theme(axis.title.y = ggplot2::element_text(size = 14)) +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 12, colour = "black")) +
  ggplot2::theme(strip.text = ggplot2::element_text(size = 12, face = "bold")) +
  ggplot2::labs(x = "", y = expression(lambda[GC]), size = 14) + 
  ggplot2::ylim(0, 5) + 
  ggplot2::guides(fill = ggplot2::guide_legend(title = "rge")) +
  ggplot2::theme(legend.position = 'right') + 
  ggplot2::scale_fill_manual(values = c("#FFFFFF", "#E0E0E0", "#BABABA", "#888888", "#000000"))
ggplot2::ggsave(filename = "Dominant_Genomic_inflation_factor.png", plot = gp, dpi = 1200, width = 7, height = 4)


