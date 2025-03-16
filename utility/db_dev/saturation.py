import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.cm import get_cmap
from scipy.optimize import curve_fit

from scipy.stats import rv_discrete

from db_dev.bootstrap import *


def top_cells_index(df, cell_types = "all", condition = "rt_condition", portion = 0.05):
    index_subset_list = []
    
    if cell_types == "none":
        return df.index

    if cell_types == "all":
        df_copy = df.copy()

    else:
        df_copy = df.copy()
        df = df[df["cell_type"].isin(cell_types)]
    
    for cond in list(df[condition].unique()):
        df_condition = df[df[condition] == cond]
    
        for cell_type in list(df["cell_type"].unique()):
            df_cell_type = df_condition[df_condition["cell_type"] == cell_type]
            index_subset = df_cell_type[(df_cell_type["1.0"] >= df_cell_type["1.0"].quantile(portion)) & (df_cell_type["1.0"] <= df_cell_type["1.0"].quantile(1 - portion))].index
            index_subset_list.append(index_subset)
          
    index_subset_list = [item for sublist in index_subset_list for item in sublist]
    
    return index_subset_list


def group_center_index(df, group = "rt_condition", center = "median", index_subset_list = []):
    df = df.iloc[index_subset_list]
    if center == "median":
        return df.groupby(group).median(numeric_only = True).transpose()
    if center == "mean":
        return df.groupby(group).mean(numeric_only = True).transpose()


def growth(X, Ymax, k):    
    return Ymax * (X/(k + X))

@np.errstate(divide = "ignore", invalid = "ignore")
def sequencing_saturation(good_reads, all_reads):
    return 1 - (good_reads/all_reads)

def find_nearest(array, value):
    array = np.asarray(array)
    idx = np.nanargmin((np.abs(array - value)))
    return idx

class Saturation():
    "This is a Saturation class for storing Saturation curve data"

    def __init__ (self, 
                  read = None,
                  umi = None,
                  gene = None,
                  rt_condition = None):
        
        self.umi = umi
        self.read = read
        self.gene = gene
        self.rt_condition = rt_condition
        self.saturation_umi_list = []
        self.saturation_gene_list = []
          
    def fit_umi(self, num_realizations = 100, conf_level = 0.95):    
        popt_exp, pcov_exp = curve_fit(growth,
                  self.read,
                  self.umi,
                  bounds=([0, 0], [50000, 50000]))
        
        self.popt_exp_umi = popt_exp
        self.Ymax_umi = popt_exp[0]
        self.k_umi = popt_exp[1]
        
        self.read_fit_umi = np.arange(0, 50000, 10)
        self.umi_fit = growth(self.read_fit_umi, *popt_exp)
        
        self.ss_fit_umi = sequencing_saturation(self.umi_fit, self.read_fit_umi)
        self.saturation_umi = self.umi_fit[find_nearest(self.ss_fit_umi, 0.8)]

        self.saturation_umi_ci1, self.saturation_umi_median, self.saturation_umi_ci2 = get_median_ci(self.saturation_umi_list, num_realizations, conf_level)
        self.saturation_umi_mean_ci1, self.saturation_umi_mean, self.saturation_umi_mean_ci2 = get_mean_ci(self.saturation_umi_list, num_realizations, conf_level)
      
    def fit_gene(self, num_realizations = 100, conf_level = 0.95):
        popt_exp, pcov_exp = curve_fit(growth,
                  self.read,
                  self.gene,
                  bounds=([0, 0], [50000, 50000]))        
        
        self.popt_exp_gene = popt_exp
        self.Ymax_gene = popt_exp[0]
        self.k_gene = popt_exp[1]
        
        self.read_fit_gene = np.arange(0, 50000, 10)
        self.gene_fit = growth(self.read_fit_gene, *popt_exp)
        
        self.ss_fit_gene = sequencing_saturation(self.gene_fit, self.read_fit_gene)
        self.saturation_gene = self.gene_fit[find_nearest(self.ss_fit_gene, 0.8)]

        self.saturation_gene_ci1, self.saturation_gene_median, self.saturation_gene_ci2 = get_median_ci(self.saturation_gene_list, num_realizations, conf_level)
        self.saturation_gene_mean_ci1, self.saturation_gene_mean, self.saturation_gene_mean_ci2 = get_mean_ci(self.saturation_gene_list, num_realizations, conf_level)
    
### Drawing saturation curves
def draw_saturation_umi(saturation_list, condition_name, n = 0, num_realizations = 100, conf_level = 0.95):
    for (i, saturation_object) in enumerate(saturation_list):
        if condition_name == saturation_object.rt_condition:
            
            condition = saturation_list[i]
            print(i, saturation_object.rt_condition)
        
    colors = get_cmap("tab10").colors
    
    condition.fit_umi(num_realizations, conf_level)
    
    plt.plot(condition.read, 
             condition.umi, 
             marker = "o", 
             linestyle = "None",
             color = colors[n],
             label = condition.rt_condition+"\n"
             +str(round(condition.saturation_umi))+"\n"
             +str(round(condition.saturation_umi_median)) + "[" + str(round(condition.saturation_umi_ci1)) + ";" + str(round(condition.saturation_umi_ci2)) + "]")
    plt.plot(condition.read_fit_umi, growth(condition.read_fit_umi, *condition.popt_exp_umi), lw = 1, color = colors[n])
    plt.xlabel("Reads per cell")
    plt.ylabel("Median UMIs per cell")

    plt.xlim([0,10000])
    plt.ylim(bottom = 0, top = 2500)
    plt.legend(loc = "upper right", bbox_to_anchor=(1.6, 1.025))


def draw_saturation_gene(saturation_list, condition_name, n = 0, num_realizations = 100, conf_level = 0.95):
    for (i, saturation_object) in enumerate(saturation_list):
        if condition_name == saturation_object.rt_condition:
            condition = saturation_list[i]
            print(i, saturation_object.rt_condition)
    colors = get_cmap("tab10").colors
    condition.fit_gene(num_realizations, conf_level)
    
    plt.plot(condition.read, 
             condition.gene, 
             marker = "o", 
             linestyle = "None",
             color = colors[n],
             label = condition.rt_condition+"\n"
             +str(round(condition.saturation_gene))+"\n"
             +str(round(condition.saturation_gene_median)) + "[" + str(round(condition.saturation_gene_ci1)) + ";" + str(round(condition.saturation_gene_ci2)) + "]")
    plt.plot(condition.read_fit_gene, growth(condition.read_fit_gene, *condition.popt_exp_gene), lw = 1, color = colors[n])

    plt.xlabel("Reads per cell")
    plt.ylabel("Median Genes per cell")
    
    plt.xlim([0,5000])
    plt.ylim(bottom = 0, top = 1500)
    plt.legend(loc = "upper right", bbox_to_anchor=(1.6, 1.025))
    

### Saturation analysis functions
def saturation_analsysis_index(read, umi, gene, group, center, cell_types, condition, portion):
    index_subset_list = top_cells_index(read, cell_types = cell_types, condition = condition, portion = portion)

    read = group_center_index(read, group = group, center = center, index_subset_list = index_subset_list)
    umi = group_center_index(umi, group = group, center = center, index_subset_list = index_subset_list)
    gene = group_center_index(gene, group = group, center = center, index_subset_list = index_subset_list)

    saturation_list = []
    
    for condition in read:
        saturation_data = Saturation(read = read[condition].values, 
                                     umi = umi[condition].values, 
                                     gene = gene[condition].values,
                                     rt_condition = condition)
        saturation_list.append(saturation_data)
    
    return saturation_list, index_subset_list, read


# Saturation analysis for single cell saturation objects
def cell_saturation_analsysis_index(read, umi, gene, group, center, cell_types, condition, portion, index_subset_list = []):
    
    saturation_list = []

    read = read.iloc[index_subset_list]
    umi = umi.iloc[index_subset_list]
    gene = gene.iloc[index_subset_list]

    read = read[["0.05", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"]].T
    read = read.set_axis(read.T.index, axis = "columns")

    umi = umi[["0.05", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"]].T
    umi = umi.set_axis(umi.T.index, axis = "columns")

    gene = gene[["0.05", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"]].T
    gene = gene.set_axis(gene.T.index, axis = "columns")

    
    for condition in read:
        saturation_data = CellSaturation(read = read[condition].values, 
                                        umi = umi[condition].values, 
                                        gene = gene[condition].values,
                                        index = condition)
        saturation_list.append(saturation_data)
    
    return saturation_list


# Saturation analysis with CI
def saturation_analysis_ci(read, umi, gene, group, center, cell_types, condition, portion):
    saturation_list_group, index_subset_list, read_grouped = saturation_analsysis_index(read, umi, gene, group, center, cell_types, condition, portion)
    saturation_list_cell = cell_saturation_analsysis_index(read, umi, gene, group, center, cell_types, condition, portion, index_subset_list)
    group_dictionary = get_groups_indexes(read, read_grouped, group)

    objects_merger(saturation_list_cell, saturation_list_group, group_dictionary)

    return saturation_list_group

# Single cell object and grouped cells object merger
def objects_merger(single_cell_objects_list, grouped_cells_objects_list, groups_indexes):
    for group in groups_indexes:
        indexes = groups_indexes[group]
        for single_cell_object in single_cell_objects_list:
            for grouped_cells_object in grouped_cells_objects_list:
                if (single_cell_object.index in indexes) & (grouped_cells_object.rt_condition == group):
                    single_cell_object.fit_umi()
                    grouped_cells_object.saturation_umi_list.append(single_cell_object.saturation_umi)
                    single_cell_object.fit_gene()
                    grouped_cells_object.saturation_gene_list.append(single_cell_object.saturation_gene)


# Get dictionary containing group name and index atributed to that group
def get_groups_indexes(original_data, grouped_data, group):
    dictionary = {key: original_data[original_data[group] == key].index for key in list(grouped_data.columns)}
    return dictionary

# Saturation class for cell
class CellSaturation():
    "This is a Saturation class for storing cell Saturation curve data"

    def __init__ (self, 
                  read = None,
                  umi = None,
                  gene = None,
                  index = None):
        
        self.umi = umi
        self.read = read
        self.gene = gene
        self.index = index


    def fit_umi(self):    
        popt_exp, pcov_exp = curve_fit(growth,
                  self.read,
                  self.umi,
                  bounds=([0, 0], [50000, 50000]))
        
        self.popt_exp_umi = popt_exp
        self.Ymax_umi = popt_exp[0]
        self.k_umi = popt_exp[1]
        
        self.read_fit_umi = np.arange(0, 50000, 10)
        self.umi_fit = growth(self.read_fit_umi, *popt_exp)
        
        self.ss_fit_umi = sequencing_saturation(self.umi_fit, self.read_fit_umi)
        self.saturation_umi = self.umi_fit[find_nearest(self.ss_fit_umi, 0.8)]
        
      
    def fit_gene(self):
        popt_exp, pcov_exp = curve_fit(growth,
                  self.read,
                  self.gene,
                  bounds=([0, 0], [50000, 50000]))        
        
        self.popt_exp_gene = popt_exp
        self.Ymax_gene = popt_exp[0]
        self.k_gene = popt_exp[1]
        
        self.read_fit_gene = np.arange(0, 50000, 10)
        self.gene_fit = growth(self.read_fit_gene, *popt_exp)
        
        self.ss_fit_gene = sequencing_saturation(self.gene_fit, self.read_fit_gene)
        self.saturation_gene = self.gene_fit[find_nearest(self.ss_fit_gene, 0.8)]

