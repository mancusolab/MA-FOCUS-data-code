import numpy as np
import scipy.linalg as lin

from itertools import chain, combinations
from numpy.linalg import multi_dot as mdot


def estimate_cor(wmat, ldmat, intercept=False):
    """
    Estimate the sample correlation structure for predicted expression.
    :param wmat: numpy.ndarray eQTL weight matrix for a risk region
    :param ldmat: numpy.ndarray LD matrix for a risk region
    :param intercept: bool to return the intercept variable or not
    :return: tuple (pred_expr correlation, intercept variable; None if intercept=False)
    """
    wcov = mdot([wmat.T, ldmat, wmat])
    scale = np.diag(1 / np.sqrt(np.diag(wcov)))
    wcor = mdot([scale, wcov, scale])

    if intercept:
        inter = mdot([scale, wmat.T, ldmat])
        return wcor, inter
    else:
        return wcor, None


def assoc_test(weights, gwas, ldmat, heterogeneity=False):
    """
    TWAS association test.
    :param weights: numpy.ndarray of eQTL weights
    :param gwas: pyfocus.GWAS object
    :param ldmat: numpy.ndarray LD matrix
    :param heterogeneity:  bool estimate variance from multiplicative random effect
    :return: tuple (beta, se)
    """

    p = ldmat.shape[0]
    assoc = np.dot(weights, gwas.Z)
    if heterogeneity:
        resid = assoc - gwas.Z
        resid_var = mdot([resid, lin.pinvh(ldmat), resid]) / p
    else:
        resid_var = 1

    se = np.sqrt(resid_var * mdot([weights, ldmat, weights]))

    return assoc, se


def get_resid(zscores, swld, wcor):
    """
    Regress out the average pleiotropic signal tagged by TWAS at the region
    :param zscores: numpy.ndarray TWAS zscores
    :param swld: numpy.ndarray intercept variable
    :param wcor: numpy.ndarray predicted expression correlation
    :return: tuple (residual TWAS zscores, intercept z-score)
    """
    m, m = wcor.shape
    m, p = swld.shape

    # create mean factor
    intercept = swld.dot(np.ones(p))

    # estimate under the null for variance components, i.e. V = SW LD SW
    wcor_inv, rank = lin.pinvh(wcor, return_rank=True)

    numer = mdot([intercept.T, wcor_inv, zscores])
    denom = mdot([intercept.T, wcor_inv, intercept])
    alpha = numer / denom
    resid = zscores - intercept * alpha

    s2 = mdot([resid, wcor_inv, resid]) / (rank - 1)
    inter_se = np.sqrt(s2 / denom)
    inter_z = alpha / inter_se

    return resid, inter_z


def bayes_factor(zscores, idx_set, wcor, prior_chisq, prb, use_log=True):
    """
    Compute the Bayes Factor for the evidence that a set of genes explain the observed association signal under the
    correlation structure.
    :param zscores: numpy.ndarray TWAS zscores
    :param idx_set: list the indices representing the causal gene-set
    :param wcor: numpy.ndarray predicted expression correlation
    :param prior_chisq: float prior effect-size variances scaled by GWAS sample size
    :param prb:  float prior probability for a gene to be causal
    :param use_log: bool whether to compute the log Bayes Factor
    :return: float the Bayes Factor (log Bayes Factor if use_log = True)
    """
    idx_set = np.array(idx_set)

    m = len(zscores)

    # only need genes in the causal configuration using FINEMAP BF trick
    nc = len(idx_set)
    cur_chi2 = prior_chisq / nc

    cur_wcor = wcor[idx_set].T[idx_set].T
    cur_zscores = zscores[idx_set]

    # compute SVD for robust estimation
    if nc > 1:
        cur_U, cur_EIG, _ = lin.svd(cur_wcor)
        scaled_chisq = (cur_zscores.T.dot(cur_U)) ** 2
    else:
        cur_U, cur_EIG = 1, cur_wcor
        scaled_chisq = cur_zscores ** 2

    # log BF 
    cur_bf = 0.5 * -np.sum(np.log(1 + cur_chi2 * cur_EIG)) + \
        0.5 * np.sum((cur_chi2 / (1 + cur_chi2 * cur_EIG)) * scaled_chisq)

    # log prior
    cur_prior = nc * np.log(prb) + (m - nc) * np.log(1 - prb)

    if use_log:
        return cur_bf, cur_prior
    else:
        return np.exp(cur_bf), np.exp(cur_prior)


def fine_map(zscores, wmats, ldmats, prior_chisq, max_genes=3, ridge=0.1, prior_prob=1e-3, credible_level=0.9, plot=False):
    p = len(zscores)

    m = len(zscores[0])
    rm = range(m)
    pips = np.zeros(m)

    wcors = []
    swlds = []

    for pdx in range(p):
        wcor, swld = estimate_cor(wmats[pdx], ldmats[pdx])
        wcors.append(wcor)
        swlds.append(swld)

    k = m if max_genes > m else max_genes
    null_res = m * np.log(1 - prior_prob)
    marginal = null_res
    # enumerate all subsets
    for subset in chain.from_iterable(combinations(rm, n) for n in range(1, k + 1)):

        local = 0
        for pdx in range(p):
            local_bf, local_prior = bayes_factor(zscores[pdx], subset, wcors[pdx], prior_chisq[pdx], prior_prob)
            if pdx == 0:
                local += local_bf + local_prior
            else:
                local += local_bf 

        # keep track for marginal likelihood
        marginal = np.logaddexp(local, marginal)

        # marginalize the posterior for marginal-posterior on causals
        for idx in subset:
            if pips[idx] == 0:
                pips[idx] = local
            else:
                pips[idx] = np.logaddexp(pips[idx], local)

    pips = np.exp(pips - marginal)
    null_res = np.exp(null_res - marginal)

    return pips, null_res
