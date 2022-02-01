#! /usr/bin/env python

# For 2 population scripts
# UPDATES - Number of SNPs can be chosen from a truncated Poisson distrubution with mean 2 and minimum 1 - "pois2"

def compute_ld(G):
    G = G.T

    # estimate LD for population from PLINK data
    n, p = [float(x) for x in G.shape]
    p_int = int(p)
    mafs = np.mean(G, axis=0) / 2
    G -= mafs * 2
    G /= np.std(G, axis=0)

    # regularize so that LD is PSD
    LD = np.dot(G.T, G) / n + np.eye(p_int) * 0.1

    # compute cholesky decomp for faster sampling/simulation
    # L = linalg.cholesky(LD, lower=True)
    L = linalg.cholesky(LD, lower=True, check_finite=False)
    mu = 2 * mafs
    adj_mu = linalg.solve_triangular(L, mu, lower=True, check_finite=False)

    # compute LD-scores for reports
    ldscs = np.sum(LD ** 2, axis=0)

    return LD, L, ldscs, mafs, p, adj_mu


def gen_ld(prefix_pop1, prefix_pop2):
    # read in plink data
    bim_1, fam_1, G_1 = read_plink(prefix_pop1, verbose=False)
    bim_2, fam_2, G_2 = read_plink(prefix_pop2, verbose=False)

    snps = pd.merge(bim_1, bim_2, how="inner", on="snp")

    G_1 = G_1[snps.i_x.values, :].compute()
    G_2 = G_2[snps.i_y.values, :].compute()

    bim_1 = bim_1.iloc[snps.i_x.values].reset_index()
    bim_2 = bim_2.iloc[snps.i_y.values].reset_index()

    LD_1, L_1, LDSCs_1, maf_1, p_1, mu1 = compute_ld(G_1)
    LD_2, L_2, LDSCs_2, maf_2, p_2, mu2 = compute_ld(G_2)

    return LD_1, L_1, LDSCs_1, bim_1, maf_1, p_1, mu1, LD_2, L_2, LDSCs_2, bim_2, maf_2, p_2, mu2


def fit_lasso(Z, y, h2g):
    """
    Infer eqtl coefficients using LASSO regression. Uses the PLINK-style coordinate descent algorithm
    that is bootstrapped by the current h2g estimate.

    :param Z:  numpy.ndarray n x p genotype matrix
    :param y: numpy.ndarray gene expression for n individuals
    :param h2g: float the -estimated- h2g from reference panel

    :return: (numpy.ndarray, float, float) tuple of the LASSO coefficients, the r-squared score, and log-likelihood
    """
    n, p = Z.shape

    def _gen_e():
        e = np.random.normal(size=n)
        return np.linalg.norm(Z.T.dot(e), np.inf)

    # PLINK-style LASSO
    lambda_max = np.linalg.norm(Z.T.dot(y), np.inf) / float(n)

    min_tmp = np.median([_gen_e() for _ in range(1000)])
    sige = np.sqrt(1.0 - h2g + (1.0 / float(n)))
    lambda_min = (sige / n) * min_tmp

    # 100 values spaced logarithmically from lambda-min to lambda-max
    alphas = np.exp(np.linspace(np.log(lambda_min), np.log(lambda_max), 100))

    # fit LASSO solution using coordinate descent, updating with consecutively smaller penalties
    lasso = lm.Lasso(fit_intercept=True, warm_start=True)
    for penalty in reversed(alphas):
        lasso.set_params(alpha=penalty)
        lasso.fit(Z, y)

    coef = lasso.coef_

    r2 = lasso.score(Z, y)
    ystar = lasso.predict(Z)
    s2e = sum((y - ystar) ** 2) / (n - 1)

    logl = sum(stats.norm.logpdf(y, loc=ystar, scale=np.sqrt(s2e)))

    return coef, r2, logl


#def sim_beta(model, share_model, rhog, eqtl_h2, n_snps):
def sim_beta(model, share_model, eqtl_h2, n_snps):
    """
    Sample qtl effects under a specified architecture.

    :param model: str the model to simulate under. choices="10pct", "1pct", "1"
    :param share_model: str the sharing model across two populations
    :param rhog: float the genetic correlation when effects are shared
    :param eqtl_h2: float the heritability of gene expression
    :param n_snps: the total number of snps at the region

    :return: numpy.ndarray of causal effects
    """
    # simulate trait architecture
    mapper = {"10pct": 0.1 * n_snps, "1pct": 0.01 * n_snps, "1snp": 1, "pois2": max(1, np.random.poisson(2))}
    n_qtls = int(mapper[model])

    # select which SNPs are causal
    c_qtls = np.random.choice(range(int(n_snps)), n_qtls)
    b_qtls = np.zeros(int(n_snps))
    b_qtls_2 = np.zeros(int(n_snps))

    # sample effects from normal prior
    b_qtls[c_qtls] = np.random.normal(loc=0, scale=np.sqrt(eqtl_h2 / n_qtls), size=n_qtls)

    if share_model != "shared":
        c_qtls = np.random.choice(range(int(n_snps)), n_qtls)
        b_qtls_2[c_qtls] = np.random.normal(loc=0, scale=np.sqrt(eqtl_h2 / n_qtls), size=n_qtls)
    else:
        #b_qtls_2[c_qtls] = np.random.normal(loc=rhog * b_qtls[c_qtls], scale=np.sqrt((1 - rhog**2) * eqtl_h2 / n_qtls), size=n_qtls)
        b_qtls_2[c_qtls] = b_qtls[c_qtls]

    return b_qtls, b_qtls_2


def sim_trait(g, h2g):
    """
    Simulate a complex trait as a function of latent genetic value and env noise.

    :param g: numpy.ndarray of latent genetic values
    :param h2g: float the heritability of the trait in the population

    :return: numpy.ndarray of simulated phenotype
    """
    n = len(g)

    if h2g > 0:
        s2g = np.var(g, ddof=1)
        s2e = s2g * ( (1.0 / h2g ) - 1 )
        e = np.random.normal(0, np.sqrt(s2e), n)
        y = g + e
    else:
        e = np.random.normal(0, 1, n)
        y = e

    # standardize
    y -= np.mean(y)
    y /= np.std(y)

    return y


def sim_geno(L, n, mu):
    """
    Sample genotypes from an MVN approximation.

    :param L: numpy.ndarray lower cholesky factor of the p x p LD matrix for the population
    :param n: int the number of genotypes to sample

    :return: numpy.ndarray n x p centered/scaled genotype matrix
    """
    p, p = L.shape

    Z = L.dot(np.random.normal(loc=mu, size=(n, p)).T).T
    Z -= np.mean(Z, axis=0)
    Z /= np.std(Z, axis=0)

    return Z


def regress(Z, pheno):
    """
    Perform a marginal linear regression for each snp on the phenotype.

    :param Z: numpy.ndarray n x p genotype matrix to regress over
    :param pheno: numpy.ndarray phenotype vector

    :return: pandas.DataFrame containing estimated beta and standard error
    """
    betas = []
    ses = []
    pvals = []
    for snp in Z.T:
        beta, inter, rval, pval, se = stats.linregress(snp, pheno)
        betas.append(beta)
        ses.append(se)
        pvals.append(pval)

    gwas = pd.DataFrame({"beta":betas, "se":ses, "pval":pvals})

    return gwas


def sim_gwas(L, ngwas, b_qtls, var_explained, mu, alpha=None):
    """
    Simulate a GWAS using `ngwas` individuals such that genetics explain `var_explained` of phenotype.

    :param L: numpy.ndarray lower cholesky factor of the p x p LD matrix for the population
    :param ngwas: int the number of GWAS genotypes to sample
    :param b_qtls: numpy.ndarray latent eQTL effects for the causal gene
    :param var_explained: float the amount of phenotypic variance explained by genetic component of gene expression
        :return: (pandas.DataFrame, float) estimated GWAS beta and standard error, causal GE effect
    """
    Z_gwas = sim_geno(L, ngwas, mu)

    # var_explained should only reflect that due to genetics
    gwas_expr = np.dot(Z_gwas, b_qtls)
    if var_explained > 0:
        if alpha is None:
            alpha = np.random.normal(loc=0, scale=np.sqrt(var_explained))
    else:
        if alpha is None:
            alpha = 0

    gwas_expr /= np.std(gwas_expr)
    y = sim_trait(gwas_expr * alpha, var_explained)

    gwas = regress(Z_gwas, y)

    return gwas, alpha


def sim_eqtl(L, nqtl, b_qtls, eqtl_h2, mu, return_ld = True):
    """
    Simulate an eQLT study using `nqtl` individuals.

    :param L: numpy.ndarray lower cholesky factor of the p x p LD matrix for the population
    :param nqtl: int the number of eQTL-panel genotypes to sample
    :param b_qtls: numpy.ndarray latent eQTL effects for the causal gene
    :param eqtl_h2: float the amount of expression variance explained by linear model of SNPs

    :return:  (pandas.DataFrame, numpy.ndarray, numpy.ndarray) DataFrame of eQTL scan, vector of LASSO eQTL coefficients
        and LD estimated from eQTL reference panel.
    """
    Z_qtl = sim_geno(L, nqtl, mu)
    n, p = [float(x) for x in  Z_qtl.shape]

    # GRM and LD
    A = np.dot(Z_qtl, Z_qtl.T) / p
    if return_ld:
        LD_qtl = np.dot(Z_qtl.T, Z_qtl) / n
    else:
        LD_qtl = None

    # simulate gene expression
    gexpr = sim_trait(np.dot(Z_qtl, b_qtls), eqtl_h2)

    # get marginal eQTLs for reporting
    eqtl = regress(Z_qtl, gexpr)

    # fit predictive model using LASSO
    h2g = her.estimate(gexpr, "normal", A, verbose=False)

    # fit LASSO to get predictive weights
    coef, r2, logl = fit_lasso(Z_qtl, gexpr, h2g)

    return (eqtl, coef, LD_qtl)

def sim_eqtl_proxy(L, nqtl, b_qtls, rho_tissue, eqtl_h2, mu, return_ld = True):
    """
    Simulate an eQLT study using `nqtl` individuals for a proxy tissue.

    :param L: numpy.ndarray lower cholesky factor of the p x p LD matrix for the population
    :param nqtl: int the number of eQTL-panel genotypes to sample
    :param b_qtls: numpy.ndarray latent eQTL effects for the causal gene
    :param eqtl_h2: float the amount of expression variance explained by linear model of SNPs

    :return:  (pandas.DataFrame, numpy.ndarray, numpy.ndarray) DataFrame of eQTL scan, vector of LASSO eQTL coefficients
        and LD estimated from eQTL reference panel.
    """
    Z_qtl = sim_geno(L, nqtl, mu)
    n, p = [float(x) for x in  Z_qtl.shape]

    # GRM and LD
    A = np.dot(Z_qtl, Z_qtl.T) / p
    if return_ld:
        LD_qtl = np.dot(Z_qtl.T, Z_qtl) / n
    else:
        LD_qtl = None

    # simulate gene expression
    gexpr = sim_trait(np.dot(Z_qtl, b_qtls), eqtl_h2)

    # simulate gene expression in a proxy tissue
    gexpr_p = np.random.normal(gexpr*rho_tissue, np.sqrt(1-rho_tissue**2))

    # get marginal eQTLs for reporting
    eqtl = regress(Z_qtl, gexpr_p)

    # fit predictive model using LASSO
    h2g = her.estimate(gexpr_p, "normal", A, verbose=False)

    # fit LASSO to get predictive weights
    coef, r2, logl = fit_lasso(Z_qtl, gexpr_p, h2g)

    return (eqtl, coef, LD_qtl)

def compute_twas(gwas, coef, LD):
    """
    Compute the TWAS test statistics.

    :param gwas: pandas.DataFrame containing estimated GWAS beta and standard error
    :param coef: numpy.ndarray LASSO eQTL coefficients
    :param LD:  numpy.ndarray p x p LD matrix

    :return: (float, float) the TWAS score and variance estimates
    """
    # compute Z scores
    Z = gwas.beta.values / gwas.se.values

    # score and variance
    score = np.dot(coef, Z)
    within_var = mdot([coef, LD, coef])

    return score, within_var


def calc_twas(score, within_var):
    if within_var > 0:
        z_twas = score / np.sqrt(within_var)
        p_twas = 2 * stats.norm.sf(np.abs(z_twas))
    else:
        # on underpowered/low-h2g genes LASSO can set all weights to 0 and effectively break the variance estimate
        z_twas = 0
        p_twas = 1

    return z_twas, p_twas


def compute(gwas, eqtl, score, within_var, bim, mafs, ldscs, b_qtls, alpha, coef):
    min_p_val = np.min(gwas.pval.values)
    mean_chi2 = np.mean((gwas.beta.values / gwas.se.values) ** 2)
    med_chi2 = np.median((gwas.beta.values / gwas.se.values) ** 2)

    z_twas, p_twas = calc_twas(score, within_var)

    # output the GWAS, eQTL, and LASSO estimates
    output = bim.drop(columns=["cm", "i"])
    output["maf"] = mafs
    output["ld.score"] = ldscs
    output["gwas.beta"] = gwas.beta
    output["gwas.se"] = gwas.se
    output["gwas.true"] = b_qtls * alpha
    output["eqtl.beta"] = eqtl.beta
    output["eqtl.se"] = eqtl.se
    output["eqtl.true"] = b_qtls
    output["eqtl.lasso"] = coef

    return output, min_p_val, mean_chi2, med_chi2, z_twas, p_twas
