#! /usr/bin/env python
import argparse as ap
import sys

import finemap_1 as fm
import numpy as np
import limix.her as her
import pandas as pd
import scipy.linalg as linalg

from numpy.linalg import multi_dot as mdot
from pandas_plink import read_plink
from scipy import stats
from sklearn import linear_model as lm

mvn = stats.multivariate_normal

exec(open('/project/nmancuso_8/zeyunlu/projects/mefocus/scripts/focus.functions_2pop.py').read())

def main(args):
    argp = ap.ArgumentParser(description="Simulate TWAS using real genotype data",
                             formatter_class=ap.ArgumentDefaultsHelpFormatter)

    argp.add_argument("prefix_pop1", help="Prefix to PLINK-formatted data for population 1")
    argp.add_argument("prefix_pop2", help="Prefix to PLINK-formatted data for population 2")

    argp.add_argument("--ngwas_pop1", default=100000, type=int, help="Sample size for GWAS panel for pop 1")
    argp.add_argument("--ngwas_pop2", default=100000, type=int, help="Sample size for GWAS panel for pop 2")

    argp.add_argument("--nqtl_pop1", default=500, type=int, help="Sample size for eQTL panel for pop 1")
    argp.add_argument("--nqtl_pop2", default=500, type=int, help="Sample size for eQTL panel for pop 2")

    argp.add_argument("--model", choices=["10pct", "1pct", "1snp", "pois2"], default="10pct",
                      help="SNP model for generating gene expression. 10pct = 10%% of SNPs, 1pct = 1%% of SNPs, 1snp = 1 SNP")

    argp.add_argument("--share-model", choices=["shared", "indep"], default="shared",
                      help="The eQTL model across populations. Shared or independent.")

    #argp.add_argument("--rhog", default=0.5, type=float, help="Amount of genetic correlation in eQTL effect sizes if 'shared' is the model")

    argp.add_argument("--eqtl-h2", default=0.1, type=float, help="The narrow-sense heritability of gene expression")

    argp.add_argument("--ngenes", default=10, type=int, help="The total number of genes in the region")

    argp.add_argument("--var-explained_p1", default=0.01, type=float,
                      help="Variance explained in complex trait by gene expression for population 1")
    argp.add_argument("--var-explained_p2", default=0.01, type=float,
                      help="Variance explained in complex trait by gene expression for population 2")
    argp.add_argument("-o", "--output", help="Output prefix")
    argp.add_argument("--seed", type=int, help="Seed for random number generation")
    argp.add_argument("--sim", type=int, help="Simulation index")
    argp.add_argument("--locus", type=int, help="Locus index")
    argp.add_argument("--pop", help="population type")

    args = argp.parse_args(args)
    np.random.seed(int(args.seed))

    LD_pop1, L_pop1, LDSCs_pop1, bim_1, maf_1, p_pop1, LD_pop2, L_pop2, LDSCs_pop2, bim_2, maf_2, p_pop2 = gen_ld(args.prefix_pop1, args.prefix_pop2)

    #b_qtls_pop1, b_qtls_pop2 = sim_beta(args.model, args.share_model, args.rhog, args.eqtl_h2, p_pop1)
    b_qtls_pop1, b_qtls_pop2 = sim_beta(args.model, args.share_model, args.eqtl_h2, p_pop1)


    # simulate GWAS under assumption that expression => downstream trait
    gwas_pop1, alpha_pop1 = sim_gwas(L_pop1, args.ngwas_pop1, b_qtls_pop1, args.var_explained_p1)
    gwas_pop2, alpha_pop2 = sim_gwas(L_pop2, args.ngwas_pop2, b_qtls_pop2, args.var_explained_p2) # alpha is not constrained

    # meta analysis
    gwas_pop1["weights"] = 1 / gwas_pop1.se.values ** 2
    gwas_pop2["weights"] = 1 / gwas_pop2.se.values ** 2
    gwas_meta = gwas_pop1.copy(deep=True)
    gwas_meta.beta = ((gwas_pop1.weights * gwas_pop1.beta) + (gwas_pop2.weights * gwas_pop2.beta)) / (gwas_pop1.weights.values + gwas_pop2.weights.values)
    gwas_meta.se = np.sqrt(1 / (gwas_pop1.weights.values + gwas_pop2.weights.values))

    # sample eQTL reference pop genotypes from MVN approx and perform eQTL scan + fit LASSO
    eqtl_pop1, coef_pop1, LD_qtl_pop1 = sim_eqtl(L_pop1, args.nqtl_pop1, b_qtls_pop1, args.eqtl_h2)
    eqtl_pop2, coef_pop2, LD_qtl_pop2 = sim_eqtl(L_pop2, args.nqtl_pop2, b_qtls_pop2, args.eqtl_h2)

    if np.isclose(np.var(coef_pop1), 0):
        coef_pop1[np.random.choice(int(p_pop1), 2, replace=False)] = 1e-6
    if np.isclose(np.var(coef_pop2), 0):
        coef_pop2[np.random.choice(int(p_pop2), 2, replace=False)] = 1e-6

    # compute TWAS statistics
    score_pop1, within_var_pop1 = compute_twas(gwas_pop1, coef_pop1, LD_qtl_pop1)
    score_pop2, within_var_pop2 = compute_twas(gwas_pop2, coef_pop2, LD_qtl_pop2)
    score_meta, within_var_meta = compute_twas(gwas_meta, coef_pop1, LD_qtl_pop1)
    z_twas_meta, _ = calc_twas(score_meta, within_var_meta)

    output_pop1, min_p_val, mean_chi2, med_chi2, z_twas_pop1, p_twas_pop1 = compute(gwas_pop1, eqtl_pop1, score_pop1, within_var_pop1, bim_1, maf_1, LDSCs_pop1, b_qtls_pop1, alpha_pop1, coef_pop1)
    output_pop2, min_p_val, mean_chi2, med_chi2, z_twas_pop2, p_twas_pop2 = compute(gwas_pop2, eqtl_pop2, score_pop2, within_var_pop2, bim_2, maf_2, LDSCs_pop2, b_qtls_pop2, alpha_pop2, coef_pop2)

    zscores_pop1 = [z_twas_pop1]
    zscores_pop2 = [z_twas_pop2]
    zscores_meta = [z_twas_meta]
    wgt_pop1 = [coef_pop1]
    wgt_pop2 = [coef_pop2]

    for idx in range(args.ngenes - 1):
        #b_qtls_pop1, b_qtls_pop2 = sim_beta(args.model, args.share_model, args.rhog, args.eqtl_h2, p_pop1)
        b_qtls_pop1, b_qtls_pop2 = sim_beta(args.model, args.share_model, args.eqtl_h2, p_pop1)

        eqtl_pop1, coef_pop1, _ = sim_eqtl(L_pop1, args.nqtl_pop1, b_qtls_pop1, args.eqtl_h2, return_ld = False)
        eqtl_pop2, coef_pop2, _ = sim_eqtl(L_pop2, args.nqtl_pop2, b_qtls_pop2, args.eqtl_h2, return_ld = False)

        if np.isclose(np.var(coef_pop1), 0):
            coef_pop1[np.random.choice(int(p_pop1), 2, replace=False)] = 1e-6
        if np.isclose(np.var(coef_pop2), 0):
            coef_pop2[np.random.choice(int(p_pop2), 2, replace=False)] = 1e-6

        wgt_pop1.append(coef_pop1)
        wgt_pop2.append(coef_pop2)

        score_pop1, within_var_pop1 = compute_twas(gwas_pop1, coef_pop1, LD_qtl_pop1)
        score_pop2, within_var_pop2 = compute_twas(gwas_pop2, coef_pop2, LD_qtl_pop2)
        score_meta, within_var_meta = compute_twas(gwas_meta, coef_pop1, LD_qtl_pop1)

        z_twas_pop1, p_twas_pop1 = calc_twas(score_pop1, within_var_pop1)
        z_twas_pop2, p_twas_pop2 = calc_twas(score_pop2, within_var_pop2)
        z_twas_meta, p_twas_meta = calc_twas(score_meta, within_var_meta)

        zscores_pop1.append(z_twas_pop1)
        zscores_pop2.append(z_twas_pop2)
        zscores_meta.append(z_twas_meta)

    wmats = [np.array(wgt_pop1).T, np.array(wgt_pop2).T]
    zscores = [np.array(zscores_pop1), np.array(zscores_pop2)]
    zscores_meta = np.array(zscores_meta)
    ldmats = [LD_qtl_pop1, LD_qtl_pop2]

    prior_prob = 1 / args.ngenes

    # Better estimate of the variance explained by the causal gene?
    V_pop1 = mdot([np.array(wgt_pop1), LD_qtl_pop1, np.array(wgt_pop1).T]) # calculate gene covariance matrix from weight and LD matrices
    V_pop1_scaled = np.diag(1 / np.sqrt(np.diag(V_pop1))) # scale this matrix
    cor_pop1 = mdot([V_pop1_scaled,  V_pop1,  V_pop1_scaled])
    prior_chisq_pop1 = mdot([zscores[0], linalg.pinv(cor_pop1), zscores[0]])

    V_pop2 = mdot([np.array(wgt_pop2), LD_qtl_pop2, np.array(wgt_pop2).T]) # calculate gene covariance matrix from weight and LD matrices
    V_pop2_scaled = np.diag(1 / np.sqrt(np.diag(V_pop2))) # scale this matrix
    cor_pop2 = mdot([V_pop2_scaled,  V_pop2,  V_pop2_scaled])
    prior_chisq_pop2 = mdot([zscores[1], linalg.pinv(cor_pop2), zscores[1]])

    unbiased_prior_chisq_pop1=prior_chisq_pop1 - args.ngenes
    unbiased_prior_chisq_pop2=prior_chisq_pop2 - args.ngenes

    prior_chisq = [unbiased_prior_chisq_pop1, unbiased_prior_chisq_pop2]

    if (unbiased_prior_chisq_pop1 < 0):
        prior_chisq[0] = prior_chisq_pop1
    if (unbiased_prior_chisq_pop2 < 0):
        prior_chisq[1] = prior_chisq_pop2

    pval_p1 = stats.distributions.chi2.sf(zscores[0]**2, df=1)
    pval_p2 = stats.distributions.chi2.sf(zscores[1]**2, df=1)
    pval_meta = stats.distributions.chi2.sf(zscores_meta**2, df=1)

    pips_j, null_res_j = fm.fine_map(zscores, wmats, ldmats, prior_prob=prior_prob, prior_chisq=prior_chisq)
    pips_p1, null_res_p1 = fm.fine_map([zscores[0]], [wmats[0]], [ldmats[0]], prior_prob=prior_prob, prior_chisq=[prior_chisq[0]])
    pips_p2, null_res_p2 = fm.fine_map([zscores[1]], [wmats[1]], [ldmats[1]], prior_prob=prior_prob, prior_chisq=[prior_chisq[1]])
    pips_meta, null_meta = fm.fine_map([zscores_meta], [wmats[0]], [ldmats[0]], prior_prob=prior_prob, prior_chisq=[prior_chisq[0]])

    # indicate true model or not
    if args.var_explained_p1 == 0:
        true_model = [1] + [0] * args.ngenes
    else:
        true_model = [0] + [1] + [0] * (args.ngenes - 1)

    output = pd.DataFrame({"sim": args.sim,
                           "locus": args.locus,
                           "pop": args.pop,
                           "snp.model": args.model,
                           "eqtl.model": args.share_model,
                           "h2g": args.eqtl_h2,
                           "h2ge1": args.var_explained_p1,
                           "h2ge2": args.var_explained_p2,
                           "n1": args.ngwas_pop1,
                           "nge1": args.nqtl_pop1,
                           "n2": args.ngwas_pop2,
                           "nge2": args.nqtl_pop2,
                           "true_model": true_model,
                           "zscores.1": [0] + list(zscores[0]),
                           "p_val.1": [1] + list(pval_p1),
                           "zscores.2": [0] + list(zscores[1]),
                           "p_val.2": [1] + list(pval_p2),
                           "zscores.meta": [0] + list(zscores_meta),
                           "p_val.meta": [1] + list(pval_meta),
                           "ME.pip": [null_res_j] + list(pips_j),
                           "pop1.pip": [null_res_p1] + list(pips_p1),
                           "pop2.pip": [null_res_p2] + list(pips_p2),
                           "meta.pip": [null_meta] + list(pips_meta),
                           })
    output.to_csv("{}.focus.tsv".format(args.output), sep="\t", index=False)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
