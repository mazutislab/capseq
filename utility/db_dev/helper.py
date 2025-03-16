import os
import numpy as np
import pandas as pd
import time
import datetime

import scipy
import scipy.cluster
from scipy import sparse

from sklearn.decomposition import PCA, TruncatedSVD
import sklearn.cluster
from sklearn.cluster import SpectralClustering

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap

from anndata import AnnData

#import fastcluster
#from adjustText import adjust_text

import pickle

import statsmodels.sandbox.stats.multicomp


# helper functions are based on https://github.com/rapolaszilionis/utility_functions

###################################################################################################

def filter_abund_genes(
                        E,
                        min_counts,
                        min_cells,
                        ):
    
    """Get boolean mask for selecting genes expressed at at least min_counts in at least min_cells.
    Input:
        E - sparse matrix (scipy.sparse)
        min_counts = counts at least
        min_cells = cells at least
        
    Return: boolean mask
    """

    gmask = np.array((E>=min_counts).sum(axis=0))[0]>=min_cells
    print(sum(gmask),"genes passing abundance filter")
    return gmask
    

###################################################################################################

def shuffle_rows(Z,seed=None,sparse=False,nonzeros=None):
    
    """
    input:
            Z - np.array, each row will be shuffled separately
            seed - seed for randomization
            sparse - boolean, default False, specify True if sparse matrix
            nonzeros - np.array with number of nonzero values per row, optional
            
    returns:
            copy of Z with shuffled columns
            
    Note: the sparse version was not particularly fast on with my toy datasets around 2000 x 3000 in size.
    Will try to improve only when faced with a very large dataset. My intuition is that it should be memory-efficient
    at least.
    """
    
    nrows = Z.shape[0]
    ncols = Z.shape[1]
    
    if sparse:
        
        # get the number of non-zero values per row
        if nonzeros is None:
            nonzeros = np.array((Z>0).sum(axis=1).T)[0]

        np.random.seed(seed)
        c = [np.random.choice(np.arange(ncols),nrix,replace=False)+ncols*rowid\
            for rowid, nrix in enumerate(nonzeros)]
        
        # flatten list of lists
        c = [item for sublist in c for item in sublist]
        
        Zshuf = Z.reshape(1,nrows*ncols).tocsr()
        
        Zshuf.indices = np.array(c)
        Zshuf = Zshuf.reshape((nrows,ncols))
        

    else:
        b = range(ncols)
        np.random.seed(seed) #to be able to recreate the exact same results
        c = np.array([np.random.permutation(b) for i in range(nrows)])
        d = np.arange(nrows)*ncols
        e = (c.T+d).T
        Zshuf = Z.flatten()[e.flatten()].reshape(Z.shape)
        
    return Zshuf
    
###################################################################################################

def find_num_pc(Z,n=10,start_pc = 200,sparse=False,
                    svd_solver='randomized',return_all_pca=False,
                    print_progression=True):
    
    """Find the number of non-random principle components.
    Steps:
        get observed eigenvalues for matrix Z
        for n rounds:
            shuffle each column of Z separately the get Zshuffled
            get the eigenvalues of Zshuffled
            record i: the number observed eigenvalues larger than the largest random eigenvalue.
            
        Consider the 5th percentile across i1,i2 ... in the number of non-random principle components
        with 95% confidence.
        
    input:
        Z - mean centered variance stabilized matrix (along the columns)
        n - number of shuffling rounds
        start_pc - guest the max number of PCs to consider, this an attempt to make the code more efficient
        sparse - boolean, default False, if True, will use truncateSVD
        svd_solver - ‘arpack’ or ‘randomized’. If sparse=False, ‘auto’ or ‘full’ are also supported.
            Check the documentation for sklearn.decomposition.PCA and sklearn.decomposition.TruncatedSVD
            for more details.
            On a quick test randomized seem faster than arpack
        return_all_pca - if True, will return a list with pca objects calculated on the random data
        print_progression - verbose or not
        
        
    returns:
        a dictionary with:
        'list_non_rand' - list with numbers of observed eigenv. larger than the largest random eigenv.
        'pca' - pca object from sklearn.decomposition after fitting the observed Zscores
        'num_pc' - integer, number of non-random PCs
        
        """
    
    start=time.time()
    
    nonzeros = None
    if sparse:
        nonzeros = np.array((Z>0).sum(axis=0))[0]
        

    # Get observed eigenvalues:
    insuff_pcs = True
    maxpc = 0
    
    if return_all_pca:
        rnd_pca = []
    
    l = []
    counter=0
    while insuff_pcs:
        # a loop to calculate more observed PCs in case maxpc (e.g. 200) not enough
        # this is an attempt to make the code more efficient
        maxpc+=int(start_pc)
        numpc = min(maxpc,Z.shape[0]-1,Z.shape[1]-1)
        
        # choose PCA method
        if sparse:
            pca = TruncatedSVD(n_components=numpc,algorithm=svd_solver)
            pca_shuff = TruncatedSVD(n_components=numpc,algorithm=svd_solver)
        else:
            pca = PCA(n_components=numpc,svd_solver=svd_solver)
            pca_shuff = PCA(n_components=numpc,svd_solver=svd_solver)
            

        # get observed eigenvalues
        if print_progression:
            print("calculating the first %d observed eigenvalues..."%numpc)
        pca.fit(Z)
        ev_obs = np.msort(pca.explained_variance_)[::-1] #make sure to sort eigenvalues
        
        # shuffle once
        counter+=1
        if print_progression:
            print("calculating the random eigenvalues for %d rounds of shuffling..."%n)
        Zshuff = shuffle_rows(Z.T,seed=counter,sparse=sparse,nonzeros=nonzeros).T
        pca_shuff.fit(Zshuff)
        
        # if more than just the random eigenvalues needed
        if return_all_pca:
            rnd_pca.append(pca_shuff)
        
        ev_rnd = max(pca_shuff.explained_variance_)
        l.append(ev_rnd)
        
        if print_progression:
            print(counter,'\t',sum(ev_obs>ev_rnd),'\t','%.2f min.'%((time.time()-start)/60.))

        insuff_pcs = (ev_obs<(ev_rnd*0.9)).sum()==0 #this sum is 0 if too few observed PCs calculated

    #iterate more
    while n>counter:
        
        # choose PCA method again
        if sparse:
            pca_shuff = TruncatedSVD(n_components=numpc,algorithm=svd_solver)
        else:
            pca_shuff = PCA(n_components=numpc,svd_solver=svd_solver)
        
        
        counter+=1
        Zshuff = shuffle_rows(Z.T,seed=counter,sparse=sparse,nonzeros=nonzeros).T
        pca_shuff.fit(Zshuff)
        
        # if more than just the random eigenvalues needed
        if return_all_pca:
            rnd_pca.append(pca_shuff)
        
        ev_rnd = max(pca_shuff.explained_variance_)
        l.append(ev_rnd)
        
        if print_progression:
            print(counter,'\t',sum(ev_obs>ev_rnd),'\t','%.2f min.'%((time.time()-start)/60.))

    #for each round of shuffling genes, what is the number of obs eigenvalues larger than the largest random eigenvalue
    nrlarger = [(ev_obs>i).sum() for i in l]
    num_pc = int(np.percentile(np.array(nrlarger),0.05))
    
    res = {'list_non_rand':nrlarger,
           'pca':pca,
           'num_pc':num_pc}
    
    if return_all_pca:
        res['rnd_pca'] = rnd_pca
    
    return res


###################################################################################################  

def spec_clust(A, k):
    """
    Spectral clustering
    Input:
    	A - sparse adjacency matrix
    	k - number of clusters to partition into
    Returns:
    	np.array with labels
    From: https://github.com/AllonKleinLab/SPRING_dev/blob/master/data_prep/helper_functions.py
    2018 12 14
    """
    spec = sklearn.cluster.SpectralClustering(n_clusters=k, random_state = 0,
                                              affinity = 'precomputed', assign_labels = 'discretize')
    return spec.fit_predict(A)

###################################################################################################  

    
def now():
    """spring current date and time as filename-friendly string"""
    return datetime.datetime.now().strftime('%y%m%d_%Hh%M')
    
#######################################################################################################
# From Adrian Veres for saving and loading pandas dataframes (modified)
def save_df(obj, filename):
    np.savez_compressed(filename, data=obj.values, index=obj.index.values, columns=obj.columns.values)
    
def load_df(filename,encoding=u'ASCII'):
    """you may want to specify encoding='latin1'
    when loading python 2 pickle with python 3.
    https://stackoverflow.com/questions/28218466/unpickling-a-python-2-object-with-python-3
    """
    with np.load(filename,encoding=encoding) as f:
        obj = pd.DataFrame(**f)
    return obj

#######################################################################################################

def startfig(w = 4, h = 2, rows = 1, columns = 1, wrs = None, hrs = None, frameon = True, return_first_ax = True):

    '''
    for initiating figures, w and h in centimeters
    example of use:
    a,fig,gs = startfig(w=10,h=2.2,rows=1,columns=3,wr=[4,50,1],hrs=None,frameon=True)
    hrs - height ratios
    wrs - width ratios
    frameon - whether first axes with frame
    
    returns:
    if return_first_ax=True
    a,fig,gs
    else
    fig,gs
    '''
    
    ratio = 0.393701 #1 cm in inch
    myfigsize = (w*ratio, h*ratio)
    fig = plt.figure(figsize = (myfigsize))
    gs = mpl.gridspec.GridSpec(rows, columns, width_ratios = wrs, height_ratios = hrs)
    if return_first_ax == True:
        a = fig.add_subplot(gs[0,0], frameon = frameon)
        return a, fig, gs
    else:
        return fig, gs


#######################################################################################################


def find_file_directories(filename, search_path, include_string = None):
    """Find all directories that contain a file in a given directory and its subdirectories.

    Only directories whose paths end with the specified include_string will be included in the results.
    If include_string is None, all directories will be included.
    """
    directories = []
    for root, dirs, files in os.walk(search_path):
        if filename in files:
            if include_string is None or root.endswith(include_string):
                directories.append(root)
    return directories

#######################################################################################################

# def find_files_by_extension(file_extension, search_path, include_string=None):
#     """Find all file paths with the specified extension in a given directory and its subdirectories.

#     Only file paths whose directories' paths end with the specified include_string will be included in the results.
#     If include_string is None, all file paths will be included.
#     """
#     file_paths = []
#     for root, dirs, files in os.walk(search_path):
#         for file in files:
#             if file.endswith(file_extension):
#                 if include_string is None or root.endswith(include_string):
#                     file_paths.append(os.path.join(root, file))
#     return file_paths

#######################################################################################################


def find_files_by_extension(file_extension, search_path, include_string=None):
    """Find all file paths with the specified extension in a given directory and its subdirectories.

    Only file paths whose directories' paths end with the specified include_string will be included in the results.
    If include_string is None, all file paths will be included.
    """
    file_paths = []
    for root, dirs, files in os.walk(search_path):
        for file in files:
            if file.endswith(file_extension):
                if include_string is None or root.endswith(include_string):
                    file_paths.append([os.path.join(root, file), file])
    
    file_paths.sort(key=lambda x: x[0])
    return file_paths


#######################################################################################################

#for saving dictionaries
def save_stuff(stuff,path):
    u"""for saving dictionaries, but probably works with lists and other pickleable objects"""
    import pickle
    with open(path+u'.pickle', u'wb') as handle:
        pickle.dump(stuff, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
#######################################################################################################
        
def load_stuff(path,encoding='ASCII'):
    """for loading object saved using 'save_stuff'.
    I had to use encoding='bytes' to load in python3 
    certain data pickled in python2."""
    import pickle
    with open(path, u'rb') as handle:
        return pickle.load(handle,encoding=encoding)    
    
#######################################################################################################

def oset(a_list):
    """given a list/1d-array, returns an ordered set (list)"""
    seen = set()
    seen_add = seen.add
    return [x for x in a_list if not (x in seen or seen_add(x))]

#######################################################################################################

def centroids(label,adata,E=None,gene_list=None):
    
    """
    Calculate average gene expression level per cell label (e.g. cluster).
    input:
        - label: name of column that stores the label of interest in adata.obs
        - adata: AnnData object OR a cell x feature pandas dataframe with label as one of the columns
        - E and gene_list: optional and only used when adata is not an AnnData object. In that case
        the cells x genes sparse expression matrix E and the gene_list must be specified
        
    returns:
        pandas dataframe, centroids x genes
        
    """
    
    if isinstance(adata,AnnData):
        E = adata.X
        gene_list = adata.var_names
        meta = adata.obs
    else:
        meta = adata
        
        
    labs = meta[label].unique()
    centroids = {}
    for lab in labs:
        msk = (meta[label] == lab).values #use "values" to turn pd.Series into row-label-less np.array,
                                             #sometimes row labels mess up the order

        centroids[lab] = np.array(E[msk,:].mean(axis=0))[0]
    centroids=pd.DataFrame(centroids).T
    centroids.columns = gene_list
    return centroids

#######################################################################################################

def hier_cluster(datatable,hier_clust_rows=True,hier_clust_cols=True,method='ward',metric='sqrt_correlation'):
    
    """
    assumes that the data table is a pandas dataframe, should also work on an numpy array.
    My favorite combinations is sqrt_correlation distance (proportional to euclidean on zscored data) with
    Ward linkage"""

    import sys
    if "fastcluster" not in sys.modules:
        import fastcluster


    data = datatable.copy()
    row_link=np.nan
    col_link=np.nan
    if hier_clust_rows:
        #hierarchically cluster:
        if metric=='sqrt_correlation':
            pdist = scipy.spatial.distance.pdist(data,metric='correlation')**0.5
        else:
            pdist = scipy.spatial.distance.pdist(data,metric=metric)
        row_link = fastcluster.linkage(pdist, method=method)
        row_order = scipy.cluster.hierarchy.leaves_list(row_link)
        try:
            #pandas-style indexing
            data = data.iloc[row_order,:]
        except:
            #numpy-style indexing
            data = data[row_order,:]
        
    if hier_clust_cols:
        #hierarchically cluster:
        if metric=='sqrt_correlation':
            pdist = scipy.spatial.distance.pdist(data.T,metric='correlation')**0.5
        else:
            pdist = scipy.spatial.distance.pdist(data.T,metric=metric)
        col_link = fastcluster.linkage(pdist, method=method)
        col_order = scipy.cluster.hierarchy.leaves_list(col_link)
        try:
            data = data.iloc[:,col_order]
        except:
            data = data[:,col_order]
        
    return {'data':data,'row_link':row_link,'col_link':col_link}

####################################################################################################### 

def yticks_fancy(a,totick,labels_all, color_dict,emptychar = '',fontsize=5, leftshift=0):
    
    """
    utility function originally made for ticking only a subset of selected genes in a genes x observations heatmap.
    example of use: yticks_fancy(a,['Csf1r','Ccr2','','','Arg1','S100a9'],genes_by_cells.index)
    input:
        a - axis with heatmap
        totick - list of yticklabels to display. Use the string defined by
        emptychar to add spacing between groups of genes.
        labels_all - all yticklabels.
        emptychar - string that will be treated as white space
        
    returns: nothing
    
    """

    a.set_yticks([])
    #leftshift = 0
    totick = np.array(totick)
    nr_slots = len(totick)
    tickmask = np.array([i!=emptychar for i in totick])
    totick = totick[tickmask]
    y_right = np.array([pd.Index(labels_all).get_loc(i) for i in totick])
    
    #if genes were not typed in in the correct order, account for that to avoid lines crossing
    tickorder = np.argsort(y_right)
    y_right = y_right[tickorder]
    totick = totick[tickorder]
    y_left = np.linspace(0,len(labels_all),nr_slots)[tickmask]
    for l,r,gene in zip(y_left,y_right,totick):
        a.plot((-0.8-leftshift,-0.5-leftshift),(r,r),lw=0.5,color='0.2')
        a.plot((-1.2-leftshift,-0.8-leftshift),(l,r),lw=0.5,color='0.2')
        a.plot((-1.5-leftshift,-1.2-leftshift),(l,l),lw=0.5,color='0.2')
        a.text(-1.6-(leftshift*1.6),l,gene,ha='right',va='center',fontsize=fontsize, color = color_dict[gene])


####################################################################################################### 

def showspines(an_axes,top=False,right=False,bottom=False,left=False):
    """
    for specifying which spines to make visible in a plot.
    input: 
        an_axes - matplotlib axes object
    returns: nothing

    """
    #after reinstalling conda, top and left switches places...
    [i for i in an_axes.spines.items()][3][1].set_visible(top) #top line
    [i for i in an_axes.spines.items()][1][1].set_visible(right) #right line
    [i for i in an_axes.spines.items()][2][1].set_visible(bottom) #bottom line
    [i for i in an_axes.spines.items()][0][1].set_visible(left) #left line
    an_axes.tick_params(bottom=bottom,right=right,left=left,top=top)


#######################################################################################################

def mwu(cg1,cg2,genes,print_progression=True):
    """perform MWU test for each gene comparing
    cells group 1 (cg1) and cell group 2 (cg2).
    Input:
        - cg1 and cg2: expression matrixes to compared, np.array, cells x genes
        - genes: gene list
        - if print_progression, will print a message every 1000 genes
    returns:
        pd.DataFrame with results, includes FDR calculation
    
    """
    
    # calculate the average per group compared to add to results.
    m1 = cg1.mean(axis=0)
    m2 = cg2.mean(axis=0)
    
    res = []
    counter = 0
    for i in range(cg1.shape[1]):
        counter+=1
        if print_progression:
            if int(counter/1000)==counter/1000.:
                print(counter)
                
        res.append(scipy.stats.mannwhitneyu(cg1[:,i],cg2[:,i],
                                            alternative='two-sided'))

    us = [i[0] for i in res]
    ps = [i[1] for i in res]
    import statsmodels
    cps = statsmodels.sandbox.stats.multicomp.multipletests(ps,method = 'fdr_bh')[1]

    return pd.DataFrame([us,ps,list(cps),list(m1),list(m2)],
                        columns=list(genes),index = ['U_statistic','p','fdr','mean1','mean2']).T

#######################################################################################################