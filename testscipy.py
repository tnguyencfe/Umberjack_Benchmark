import scipy.stats
import numpy as np

exp_freq_distro_arr = [999, 1]
scaled_exp_freq_distro_arr = [55.944, 0.056]
obs_freq_distro_arr = [56, 0]


# chi2 : float
#         The test statistic.
#     p : float
#         The p-value of the test
#     dof : int
#         Degrees of freedom
#     expected : ndarray, same shape as `observed`
#         The expected frequencies, based on the marginal sums of the table.

result = scipy.stats.chi2_contingency(observed=[scaled_exp_freq_distro_arr, obs_freq_distro_arr], lambda_="log-likelihood")
print result


result = scipy.stats.chi2_contingency(observed=[obs_freq_distro_arr, scaled_exp_freq_distro_arr], lambda_="log-likelihood")
print result

result = scipy.stats.chi2_contingency(observed=[obs_freq_distro_arr, exp_freq_distro_arr], lambda_="log-likelihood")
print result

result = scipy.stats.chi2_contingency(observed=[[0.56, 0.01], [0.999, 0.001]], lambda_="log-likelihood")
print result


result = scipy.stats.chi2_contingency(observed=[[56, 1], [999, 1]], lambda_="log-likelihood")
print result

result = scipy.stats.chi2_contingency(observed=[[56, 1], [9990, 10]], lambda_="log-likelihood")
print result

result = scipy.stats.chi2_contingency(observed=[[56, 1], [0.999, 0.001]], lambda_="log-likelihood")
print result

result = scipy.stats.chi2_contingency(observed=[[0.56, 0.1], [999, 1]], lambda_="log-likelihood")
print "here" + str(result)

result = scipy.stats.chi2_contingency(observed=[[56, 0], [999, 1]], lambda_="log-likelihood")
print "zero" + str(result)

result = scipy.stats.chi2_contingency(observed=[[56, 0.1], [999, 1]], lambda_="log-likelihood")
print "dec" + str(result)

result = scipy.stats.chi2_contingency(observed=[[56, 1], [999, 0.1]], lambda_="log-likelihood")
print "dec exp" + str(result)


result = scipy.stats.chi2_contingency(observed=[[56, 0], [999, 1]], lambda_="log-likelihood")
print "Scale1" + str(result)

result = scipy.stats.chi2_contingency(observed=[[560, 0], [999, 1]], lambda_="log-likelihood")
print "Scale2" + str(result)


# https://gist.github.com/brentp/570896
def gtest(f_obs, f_exp=None, ddof=0):
    """
    http://en.wikipedia.org/wiki/G-test
    The G test can test for goodness of fit to a distribution
    Parameters
    ----------
    f_obs : array
        observed frequencies in each category
    f_exp : array, optional
        expected frequencies in each category.  By default the categories are
        assumed to be equally likely.
    ddof : int, optional
        adjustment to the degrees of freedom for the p-value
    Returns
    -------
    chisquare statistic : float
        The chisquare test statistic
    p : float
        The p-value of the test.
    Notes
    -----
    The p-value indicates the probability that the observed distribution is
    drawn from a distribution given frequencies in expected.
    So a low p-value inidcates the distributions are different.
    Examples
    --------
    >>> gtest([9.0, 8.1, 2, 1, 0.1, 20.0], [10, 5.01, 6, 4, 2, 1])
    (117.94955444335938, 8.5298516190930345e-24)
    >>> gtest([1.01, 1.01, 4.01], [1.00, 1.00, 4.00])
    (0.060224734246730804, 0.97033649350189344)
    >>> gtest([2, 1, 6], [4, 3, 2])
    (8.2135343551635742, 0.016460903780063787)
    References
    ----------
    http://en.wikipedia.org/wiki/G-test
    """
    f_obs = np.asarray(f_obs, 'f')
    k = f_obs.shape[0]
    f_exp = np.array([np.sum(f_obs, axis=0) / float(k)] * k, 'f') \
                if f_exp is None \
                else np.asarray(f_exp, 'f')
    g = 2 * np.add.reduce(f_obs * np.log(f_obs / f_exp))
    return g, scipy.stats.chisqprob(g, k - 1 - ddof)


result = gtest(f_obs=scaled_exp_freq_distro_arr, f_exp=obs_freq_distro_arr)
print result


result = gtest(f_obs=obs_freq_distro_arr, f_exp=scaled_exp_freq_distro_arr)
print result

result = gtest(f_obs=obs_freq_distro_arr, f_exp=exp_freq_distro_arr)
print result