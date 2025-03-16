import pandas as pd
import numpy as np
import scipy.stats


# bootstraping is based on https://aegis4048.github.io/non-parametric-confidence-interval-with-bootstrap

def bootstrap_simulation(sample_data, num_realizations = 100):
    sample_data = pd.DataFrame(sample_data)

    n = sample_data.shape[0]
    boot = []
    for i in range(num_realizations):
        real = np.random.choice(sample_data.values.flatten(), size=n)
        boot.append(real)
        
    columns = ['Real ' + str(i + 1) for i in range(num_realizations)]
    
    return pd.DataFrame(boot, index=columns).T

def calc_sum_stats(boot_df):
    boot_df = pd.DataFrame(boot_df)

    sum_stats = boot_df.describe().T[['mean', 'std', 'min', 'max']]
    sum_stats['median'] = boot_df.median()
    sum_stats['skew'] = boot_df.skew()
    sum_stats['kurtosis'] = boot_df.kurtosis()
    sum_stats['IQR'] = boot_df.quantile(0.75) - boot_df.quantile(0.25)
    return sum_stats.T


def calc_bounds(conf_level):
    
    assert (conf_level < 1), "Confidence level must be smaller than 1"
    
    margin = (1 - conf_level) / 2
    upper = conf_level + margin
    lower = margin
    return margin, upper, lower


def calc_confidence_interval(df_sum_stats, conf_level = 0.95): 
    
    margin, upper, lower = calc_bounds(conf_level)
    
    conf_int_df = df_sum_stats.T.describe(percentiles=[lower, 0.5, upper]).iloc[4:7, :].T
    conf_int_df.columns = ['P' + str(round(lower * 100, 1)), 'P50', 'P' + str(round(upper * 100, 1))]
    return conf_int_df 


def print_confidence_interval(conf_df, conf_level):
    print('By {}% chance, the following statistics will fall within the range of:\n'.format(round(conf_level * 100, 1)))
    
    margin, upper, lower = calc_bounds(conf_level)
    
    upper_str = 'P' + str(round(upper * 100, 1))
    lower_str = 'P' + str(round(lower * 100, 1))
    
    for stat in conf_df.T.columns:
        lower_bound = round(conf_df[lower_str].T[stat], 1)
        upper_bound = round(conf_df[upper_str].T[stat], 1)

        mean = round(conf_df['P50'].T[stat], 1)
        print("{0:<10}: {1:>10}  ~ {2:>10} , AVG = {3:>5}".format(stat, lower_bound, upper_bound, mean))



def get_median_ci(df_list, num_realizations, conf_level):
    bootstrap_data = bootstrap_simulation(pd.DataFrame(df_list), num_realizations = num_realizations)
    boot_sum_stats = calc_sum_stats(bootstrap_data)
    conf_int = calc_confidence_interval(boot_sum_stats, conf_level = conf_level)

    ci_1 = round(conf_int.T["median"].iloc[0], 1)
    median = round(conf_int.T["median"].iloc[1], 1)
    ci_2 = round(conf_int.T["median"].iloc[2], 1)

    return ci_1, median, ci_2


def get_mean_ci(df_list, num_realizations, conf_level):
    bootstrap_data = bootstrap_simulation(pd.DataFrame(df_list), num_realizations = num_realizations)
    boot_sum_stats = calc_sum_stats(bootstrap_data)
    conf_int = calc_confidence_interval(boot_sum_stats, conf_level = conf_level)

    ci_1 = round(conf_int.T["mean"].iloc[0], 1)
    mean = round(conf_int.T["mean"].iloc[1], 1)
    ci_2 = round(conf_int.T["mean"].iloc[2], 1)

    return ci_1, mean, ci_2