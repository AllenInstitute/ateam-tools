import warnings
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
import statsmodels.api as sm
from statsmodels.stats.anova import anova_lm

def ols_model(data, formula_rhs, feature, anova_type=2, cov_type='HC3'):
    metrics = ['aic', 'bic', 'fvalue', 'f_pvalue', 'llf', 'rsquared', 'rsquared_adj', 'nobs']
    formula = f"{feature} ~ {formula_rhs}"
    res = smf.ols(formula=formula, data=data).fit(cov_type=cov_type)
    fit_dict = {name: getattr(res, name) for name in metrics}

    if anova_type:
        anova = anova_lm(res, typ=anova_type)
        pvals = anova["PR(>F)"].dropna().rename(lambda x: f"pval_{x}")
        fvals = anova["F"].dropna().rename(lambda x: f"fval_{x}")
        eta = anova["sum_sq"].dropna().apply(lambda x: x/(x+res.ssr)).rename(lambda x: f"eta_p_{x}")
        fit_dict.update(pvals.to_dict())
        fit_dict.update(fvals.to_dict())
        fit_dict.update(eta.to_dict())
    return fit_dict, res
    
def fit_models(data, formulas, features, formula_names=None, feature_names=None, cov_type='HC3'):
    feature_names = feature_names or [str(x) for x in features]
    formula_names = formula_names or [str(x) for x in formulas]
    all_fits = []
    for feature, feature_name in zip(features, feature_names):
        for formula, formula_name in zip(formulas, formula_names):
            fit_dict, results = ols_model(data, formula, feature, cov_type=cov_type)
            all_fits.append(dict(fit_dict, model=formula_name, feature=feature_name))
        
    fits_df = pd.DataFrame(all_fits)
    return fits_df

def anova_all(data, feature, f1, f2='cluster', f1_name='gene', f2_name='cluster', cov_type='HC3', weighted=False, quote_f1=True):
    out = {
        "feature":feature,
        f1_name:f1,
        f2_name:f2
    }
    data = data.dropna(subset=[f1, f2])
    if quote_f1:
        f1 = f"Q('{f1}')"
    weights = data["SeuratMappingScore"]**2 if weighted else 1.
    res3 = smf.wls(formula=f"{feature}~{f1}*{f2}", data=data, weights=weights).fit(cov_type=cov_type)
    anova = anova_lm(res3, typ=2)
    out.update({
                    f"p_{f2_name}_c_{f1_name}": anova.loc[f2,"PR(>F)"],
                    f"p_{f1_name}_c_{f2_name}": anova.loc[f1,"PR(>F)"],
                    "p_interaction": anova.loc[f"{f1}:{f2}","PR(>F)"],
    })
    
    res1 = smf.wls(formula=f"{feature}~{f1}", data=data, weights=weights).fit(cov_type=cov_type)
    res2 = smf.wls(formula=f"{feature}~{f2}", data=data, weights=weights).fit(cov_type=cov_type)
    res12 = smf.wls(formula=f"{feature}~{f1}+{f2}", data=data, weights=weights).fit(cov_type=cov_type)
    anova_add = anova_lm(res12, typ=2)
    # anova3 = anova_lm(res12, res3)
    out.update({
        f"p_{f1_name}": res1.f_pvalue[0][0] if not np.isscalar(res1.f_pvalue) else res1.f_pvalue,
        f"p_{f2_name}": res2.f_pvalue[0][0] if not np.isscalar(res2.f_pvalue) else res2.f_pvalue,
        f"p_{f2_name}_c_{f1_name}_add": anova_add.loc[f2,"PR(>F)"],
        f"p_{f1_name}_c_{f2_name}_add": anova_add.loc[f1,"PR(>F)"],
        # "p_interaction_2":anova3.loc[1,"Pr(>F)"],
    })
    return out

def plot_fit(data, feature, formula, x=None, res=None, cluster='cluster', ax=None, legend=False, print_attr=None, print_pvals=True, cv=False, **sns_args):
    if not ax:
        fig, ax = plt.subplots()
#     data = data.dropna(subset=variables+[feature])
    if res is None:
        out_dict, res = ols_model(data, formula, feature)
    if cv:
        cv_dict = ols_cv(data, formula, feature)
        out_dict.update(cv_dict)

    x = x or formula.replace('*','+').split('+')[0].strip()
    # TODO: use hue_order 
    data = data.sort_values(cluster)
    sns_args['s'] = sns_args.get('s', 25)
    sns.scatterplot(data=data, y=feature, x=x, hue=cluster, ax=ax, legend=legend, **sns_args)
    hue = cluster if cluster in formula else None
    c = None if cluster in formula else 'k'
    y_fit = res.fittedvalues.reindex(data.index)
    sns.lineplot(data=data, y=y_fit, x=x, hue=hue, color=c, legend=False, ax=ax)
    ax.set_xlabel(getattr(x, "label", None) or x)
    ax.set_ylabel(getattr(feature, "label", None) or feature)
    # if legend:
        # plt.legend(bbox_to_anchor=(1,1))
    summary = ''
    if print_attr:
        value = out_dict.get(print_attr)
        attr_name = getattr(print_attr, "label", None) or print_attr
        summary = f"{attr_name} = {value:.2g}\n"
    if print_pvals:
        anova = anova_lm(res, typ=2)
        pvals = anova["PR(>F)"].dropna()
        # summary += ", ".join(f"$p_{{{key}}}$={pvals[key]:.2g}" for key in pvals.index)
        summary += ", ".join(f"p_{key}={pvals[key]:.2g}" for key in pvals.index)
    ax.text(0.5, 0.99, summary, transform=ax.transAxes,
        verticalalignment='top', horizontalalignment='center')
    sns.despine()
    return out_dict

def plot_grid(data, plot_fcn, row_args, col_args,  col_titles=None, sharex=False, sharey=False, figsize=None, **plot_args):
    nrows = len(row_args)
    ncols = len(col_args)
    fig, axs = plt.subplots(nrows, ncols, sharex, sharey, figsize=figsize)
    if len(axs.shape)==1:
        axs = np.reshape(axs, (nrows, ncols))
    for i in np.arange(nrows):
        for j in np.arange(ncols):
            ax = axs[i,j]
            plot_fcn(data, row_args[i], col_args[j],  ax=ax, **plot_args)
            if i==0 and col_titles:
                ax.set_title(col_titles[j])
    return fig, axs

def plot_formula_lhs(data, feature, formula, ax=None, **kwargs):
    if not ax:
        fig, ax = plt.subplots()
    x = formula.replace('*','+').split('+')[0].strip()
    return plot_fit(data, feature, formula, x, ax=ax, **kwargs)

def plot_formula_xy(data, formula_pattern, feature_x, feature_y, ax=None, **kwargs):
    if not ax:
        fig, ax = plt.subplots()
    formula = formula_pattern.format(feature=feature_x)
    return plot_fit(data, feature_y, formula, x=feature_x, ax=ax, **kwargs)


# different fit options
import patsy
import sklearn.linear_model as sklm
import sklearn.model_selection as skms
# TODO: make this fully interchangeable, not using res from statsmodels result for plot
def ols_cv(data, formula_rhs, feature):
    formula = f"{feature} ~ {formula_rhs}"
    y, X = patsy.dmatrices(formula, data, return_type='matrix')
    scorers = ["r2", "explained_variance"]
    # TODO: adapt stratified KFold to use groups
    # cv = skms.RepeatedStratifiedKFold(n_splits=5, n_repeats=10)
    cv = skms.RepeatedKFold(n_splits=3, n_repeats=10)
    estimator = sklm.LinearRegression()
    scores = skms.cross_validate(estimator, X, y, scoring=scorers, cv=cv)
    fit_dict = {key: np.mean(val) for key, val in scores.items() if "test" in key}
    std_dict = {key+"_std": np.std(val) for key, val in scores.items() if "test" in key}
    fit_dict.update(std_dict)
    return fit_dict

# not sure if this approach is a valid use of ANOVA
def fit_partial_models(data, features, regressor, base_formulas):
    all_fits = []
    features = [str(x) for x in features]
    base_formulas = [str(x) for x in base_formulas]
    for feature in features:
        for formula in base_formulas:
            fitdata = data.dropna(subset=[feature, regressor]+formula.split('+'))
            _, base = ols_model(fitdata, formula, feature, anova_type=False)
            Y = base.resid
            X = patsy.dmatrix(regressor, fitdata)
            full = sm.OLS(Y,X).fit()
            anova = anova_lm(base, full)
            pvals = anova["Pr(>F)"].dropna().rename(lambda x: f"pval_{x}")
            all_fits.append(dict(pvals, model=formula, feature=feature))
    fits_df = pd.DataFrame(all_fits).pivot(index='feature', columns='model')
    return fits_df