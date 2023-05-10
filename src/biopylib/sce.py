import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns; sns.set(style='whitegrid')
import plotly.express as px
import plotly.graph_objects as go
import scipy
from pathlib import Path,PurePath

# Dimensionality Reduction
from sklearn.decomposition import TruncatedSVD
from sklearn.decomposition import PCA
from sklearn.decomposition import SparsePCA
from sklearn.decomposition import KernelPCA
from sklearn.decomposition import IncrementalPCA
from sklearn.decomposition import TruncatedSVD
from sklearn.decomposition import MiniBatchDictionaryLearning
from sklearn.decomposition import FastICA
from sklearn.manifold import Isomap
from sklearn.manifold import MDS
from sklearn.manifold import LocallyLinearEmbedding
from sklearn.manifold import TSNE
from sklearn.random_projection import GaussianRandomProjection
from sklearn.random_projection import SparseRandomProjection
import umap 
from sklearn import metrics

# Clustering
from sklearn.cluster import KMeans
from sklearn.cluster import DBSCAN
from sklearn.cluster import MeanShift
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import FeatureAgglomeration

'''

CLASS FOR BASIC SCE Related Operations

'''

class SCE:
  
    # constructor
    def __init__(self,file=None,          # set SCE file path
                      title='SCE File',   # set SCE file name
                      rs = 72):           # set random number id 
      
        self.file = Path(file)            # read file string 
        self.title= title                 # case name string
        self.assays = {}                  # dictionary containing main data
        self.reducedDim = {}              # dictionary containing dim red matrices
      
        self.rs = rs                      # random state
        self.dic_dr = None                # reserve for dim red dictionary
      
    ''' Show loaded SCE data info '''
    # print loaded SCE data 
      
    def info(self):
        print('SCE Instance Information:')
        print(f'Loaded File: {self.file}')
        print(f'Case Name: {self.title}')
        print(f'Input Data Shape: {self.__data.X.shape}')
        print(f'Assays ({len(self.assays)}) : {([ i for i in self.assays.keys()])}')
        row_names = (f'rownames({self.__data.X.shape[0]}) :'
                     f' [{self.__data.obs.index[0]} {self.__data.obs.index[1]}'
                     f' ... {self.__data.obs.index[-1]} ]')
        print(row_names)
        col_names = (f'colnames({self.__data.X.shape[1]}) : [{self.__data.var.index[1]}'
                     f'{self.__data.var.index[1]}'
                     f' ... {self.__data.var.index[-1]} ]')
        print(col_names)
      
        coldata_names = (f'colData names({len(self.__data.obs.columns.tolist())})'
                         f': {[i for i in self.__data.obs.columns.tolist()]}')
      
        rowdata_names = (f'rowData names({len(self.__data.var.columns.tolist())})'
                         f': {[i for i in self.__data.var.columns.tolist()]}')
        print(coldata_names)
        print(rowdata_names)
        reddim_names = (f'reducedDimNames ({len(self.reducedDim.keys())}) : '
                        f'{([ i for i in self.reducedDim.keys()])}')
        print(reddim_names)
      
    ''' READ SCE DATASET '''
      
    # Read H5AD format using scanpy
    def read_h5ad(self):
      
        load_id = (f"Loading : {self.file.name.split('.')[0]}"
                   f"| {self.file.name.split('.')[1]} Format")
        print(load_id)
      
        # load H5AD data
        self.__data = sc.read_h5ad(self.file)
      
        # list of genes (columns)
        genes = [t.upper() for t in self.__data.var.index ]
        self.__data.var.index = genes
      
        # store count matrix
        self.assays['counts'] = pd.DataFrame(self.__data.X,
                                columns=self.__data.var.index,
                                index=self.__data.obs.index)
      
        # observations (additional data for each cell)
        self.obs = self.__data.obs 
      
        # variables (additional data for each gene)
        self.var = self.__data.var 
      
        print('[note] h5ad file loaded')
      
    ''' Quality Control '''
      
    def qc(self):
      
        def total_counts():
            ser = self.assays['counts'].sum(axis=1)
            ser.name = 'total_counts'
            return ser
      
        # The number of genes with at least 1 count in a cell. 
        # calculated for all cells
        def n_genes_by_counts():
          
            data = self.assays['counts']
            ser = (data > 0.0).sum(axis=1)
            ser.name = 'n_genes_by_counts'
            return ser
      
        tcounts = total_counts()
        ngbyc = n_genes_by_counts()
        data = tcounts.to_frame().join(ngbyc)
      
        return data
  
    ''' SCE DATA analysis '''
  
    # expressions per cell
    def epc(self):
      
        fig = plt.figure(figsize = (13,6)); c = 0
        c+=1; fig.add_subplot(1,2,c);
        v = np.asarray(self.assays['counts'].sum(axis = 1)).ravel()
        plt.plot(np.sort(v)) 
        plt.title('Expression per cell')
        plt.xlabel('cells sorted');plt.ylabel('Expression Total')
      
        c+=1; fig.add_subplot(1,2,c);
        plt.hist(np.sort(v), bins = 100) 
        plt.title('Expression per cell')
      
        plt.show()
        pd.Series(v).describe()
      
        fig = plt.figure(figsize = (13,6)); c = 0
        c+=1; fig.add_subplot(1,2,c);
        v = np.asarray((self.assays['counts'] != 0).sum(axis = 1)).ravel() 
        plt.plot(np.sort(v)) 
        plt.title('Count genes expressed per cell')
        plt.xlabel('cells sorted')
        plt.ylabel('Count expressed genes')
      
        c+=1; fig.add_subplot(1,2,c);
        plt.hist(np.sort(v), bins = 100) 
        plt.title('Count genes expressed per cell')
      
        plt.show()
        pd.Series(v).describe()
      
    ''' LOGCOUNTS '''
    # For scaling of feature data, log data 
      
    # arguments: assay - str name from 
      
    def logcounts(self,assay='counts'):
      
        if(assay not in list(self.assays.keys())):
            print(f'{assay} not found in assays matrix')
        else:
            print(f'[note] {assay} found')
          
        # column sums (for each feature/gene)
        libsizes = self.assays[assay].sum(axis=1)
        # scale factor for each col
        sizefactors = libsizes/np.mean(libsizes) 
        # scale data w/ log & store 
        logcounts = np.log2((self.assays[assay].T/sizefactors).T + 1) 
        self.assays['logcounts'] = logcounts
        print(f'[note] Finished creating logcounts for matrix {assay}')
      
    ''' CLUSTERING '''
      
    # Show cluster information 
    def clust_info(self):
      
        print('Hyperparameter Information')
        print('dbscan: [eps][min_samples]')
      
    # Main clustering method 
    def clust(self,assay='counts',
                   method='kmeans',
                   nclust=3,  # general
                   eps=1,min_samples=1 # DBSCAN
             ):
      
        X = self.assays[assay]
      
        if(method is 'kmeans'):
          
            note_id1 = (f"Creating cluster on assay: {assay}",
                        f"using {method}, {nclust} components")
            print(note_id1)
            kmeans_model = KMeans(n_clusters=nclust,
                                  random_state=self.rs).fit(X)
          
            labels = kmeans_model.labels_
            score = metrics.silhouette_score(X,
                                             labels, 
                                             metric='euclidean')
          
            print(f'Silhouette Score: {score}')
            self.obs['kmeans'] = labels
          
        # Don't need to specify the number of clusters
        elif(method is 'dbscan'):
          
            print(f'Creating cluster on assay: {assay} using {method}')
            model = DBSCAN(eps=eps,
                           min_samples=min_samples)
            model.fit(X)
            print(model.labels_) # cluster assignments
            self.obs['dbscan'] = model.labels_
          
        elif(method is 'agglomerative'):
          
            print(f'Creating cluster on assay: {assay} using {method}') 
            model = AgglomerativeClustering(n_clusters=nclust)
            model.fit(X)
          
            # cluster assignments
            print(model.labels_)
            self.obs['agglomerative'] = model.labels_
          
        else:
          
            print('Method not found')
          
#         # Don't need to specift the number of clusters
#         elif(method is 'meanshift'):
          
#             model = MeanShift()
#             model.fit(X)
          
#             print(f'labels: {model.labels_}') # cluster assignment
#             print(model.cluster_centers_) # cluster centroids
          
    ''' DIMENSIONALITY REDUCTION '''
    # Set Parameter Dictionary for DR
          
    def param_dr(self,method='PCA'):
      
        self.__param_dr_ID = method
      
        if(method is 'PCA'):
            self.dic_dr = {"n_components":2}
        elif(method is 'SPCA'):
            self.dic_dr = {'n_components':2,"alpha":1,
                            "ridge_alpha":0.01,"random_state":self.rs,
                            "max_iter":1000,"n_jobs":-1,
                            "tol":1e-08,"method":"lars"}
        if(method is 'KPCA'):
            self.dic_dr = {'n_components':2,"kernel":"linear",
                        "gamma":None,"degree":3,"coef0":1,
                        "kernel_params":None,"alpha":1.0,
                        "fit_inverse_transform":False,
                        "eigen_solver":"auto","tol":0,
                        "max_iter":None,"iterated_power":'auto',
                        "remove_zero_eig":False,'random_state':self.rs,
                        "n_nobs":-1}
          
        elif(method is 'IPCA'):
            self.dic_dr = {'n_components':2,'whiten':False,'batch_size':None}
        elif(method is 'TSVD'):
            self.dic_dr = {'n_components':2, 'algorithm':'randomized', 
                           'n_iter':5, 'random_state':self.rs,'tol':0.0}
        elif(method is 'GPR'):
            self.dic_dr = {'n_components':'auto','eps':0.1,'random_state':self.rs}
        elif(method is 'SPR'):
            self.dic_dr = {'n_components':'auto', 'density':'auto', 'eps':0.1,
                           'dense_output':False, 'random_state': self.rs}
        elif(method is 'MBD'):
            self.dic_dr = {'n_components':2, 'alpha':1, "n_iter":1000, 
                           "fit_algorithm":'lars', 'n_jobs':-1, 
                           'batch_size':3, 'shuffle':True, 'dict_init':None, 
                           'transform_algorithm':'omp','transform_n_nonzero_coefs':None,
                           'transform_alpha':None, 'verbose':False,'split_sign':False,
                           'random_state':self.rs, 'positive_code':False,
                           'positive_dict':False,'transform_max_iter':1000}
        elif(method is 'ISO'):
            self.dic_dr = {'n_neighbors':5,'n_components':2,'eigen_solver':'auto',
                           'tol':0, 'max_iter':None,'path_method':'auto', 
                           'neighbors_algorithm':'auto', 'n_jobs':self.rs,
                           'metric':'minkowski','p':2, 'metric_params':None}
        elif(method is 'MDS'):
            self.dic_dr = {'n_components':2,'metric':True,"n_init":4, 
                           'max_iter':300,'verbose':0, 'eps':0.001, 'n_jobs':None,
                           'random_state':self.rs,'dissimilarity':'euclidean'}
        elif(method is 'LLE'):
            self.dic_dr = {'n_neighbors':5, 'n_components':2, 'reg':0.001, 
                           'eigen_solver':'auto','tol':1e-06, 'max_iter':100, 
                           'method':'standard', 'hessian_tol':0.0001,
                           'modified_tol':1e-12, 'neighbors_algorithm':'auto', 
                           'random_state':None, 'n_jobs':self.rs}
        elif(method is 'TSNE'):
            self.dic_dr = {'n_components':2,'perplexity':30.0, 
                           'early_exaggeration':12.0,'learning_rate':'warn', 'n_iter':1000, 
                            'n_iter_without_progress':300, 'min_grad_norm':1e-07,
                            'metric':'euclidean', 'init':'warn', 'verbose':0, 
                            'random_state':self.rs, 'method':'barnes_hut', 'angle':0.5,
                            'n_jobs':self.rs, 'square_distances':'legacy'}
        elif(method is 'FICA'):
            self.dic_dr = {'n_components':2,'algorithm':'parallel','whiten':True, 
                        'fun':'logcosh', 'fun_args':None,'max_iter':200, 'tol':0.0001,
                        'w_init':None, 'random_state':self.rs}
          
        print('> DR Parameters Set')
      
    ''' Dimensional Reduction '''
    # dimensionally reduce the assay data
      
    def dr(self,assay='counts',    # choose which assay to apply DR
                method='PCA',      # dimensionality reduction approach
                n_comp=2,          # reduce the number of dimensions to n_comp
                use_dimRed=None):  
      
        # If we already have parameters set
        if(self.dic_dr is not None):
          
            # if set parameter method is not identical
            if(self.__param_dr_ID is not method):
                self.param_dr(method=method)
                self.dic_dr['n_components'] = n_comp # set basic param
                verb_param_id = (f'param_id != dr method, setting' 
                                 f'default parameters for {method}') 
                print(verb_param_id)
              
        # use default parameters if nothing is preset
        elif(self.dic_dr is None):
            self.param_dr(method=method)
            self.dic_dr['n_components'] = n_comp # set basic param
          
        # Select Unsupervised Learning Model
        if(method is 'PCA'): model = PCA(**self.dic_dr)
        if(method is 'SPCA'): model = SparsePCA(**self.dic_dr)    
        if(method is 'KPCA'): model = KernelPCA(**self.dic_dr)
        if(method is 'IPCA'): model = IncrementalPCA(**self.dic_dr)
        if(method is 'TSVD'): model = TruncatedSVD(**self.dic_dr)
        if(method is 'GRP'):  model = GaussianRandomProjection(**self.dic_dr)    
        if(method is 'SRP'): model = SparseRandomProjection(**self.dic_dr)
        if(method is 'MBD'): model = MiniBatchDictionaryLearning(**self.dic_dr)
        if(method is 'ISO'): model = Isomap(**self.dic_dr)   
        if(method is 'MDS'): model = MDS(**self.dic_dr)
        if(method is 'LLE'): model = LocallyLinearEmbedding(**self.dic_dr)    
        if(method is 'TSNE'): model = TSNE(**self.dic_dr)
        if(method is 'FICA'): model = FastICA(**self.dic_dr)
        if(method is 'UMAP'): model = umap.UMAP(**self.dic_dr)
      
        if(use_dimRed is None): # Fit on Assay Data
      
            X = model.fit_transform(self.assays[assay])
            self.reducedDim[method] = X # Save Data into 
            print(f'Data saved into reducedDim: {method}')
      
        else: # Fit on already dimensionally reduced data
      
            print(f'Fitting on reducedDim data: {use_dimRed}')
            X = model.fit_transform(self.reducedDim[use_dimRed])
            self.reducedDim[f'{use_dimRed}_{method}'] = X # Save Data into 
            print(f'Data saved into reducedDim: {use_dimRed}_{method}')
      
        if(method is 'PCA'):
            self.__pcavar = np.round(model.explained_variance_ratio_.sum()*100,2)
            print(f'PCA - Total Explained Variance: {self.__pcavar} %')
            self.__pcalabels = {
                str(i): f"PC{i+1} ({var:.1f}%)"
                for i, var in enumerate(model.explained_variance_ratio_ * 100)
            }
            print(self.__pcalabels)
          
        print('\nParameters used in DR:')
        print(self.dic_dr,'\n')
      
    ''' Plot Dimensionally Reduced Data '''
      
    def plot_dr(self,method='PCA',components=[0,1],hue=None):
      
        # activate only if we have applied dr method
        if(method in self.reducedDim.keys()):
          
            ttitle = self.title + f' {method}  cells: ' + \
                      str(self.__data.X.shape[0]) + \
                      ' genes: ' + str(self.__data.X.shape[1]) + \
                      ' | geneDim: ' + str(self.reducedDim[method].shape[1]) 
          
            # Plot Dimensions (1,2)
            if(hue is not None):
                fig = px.scatter(self.reducedDim[method],
                                 x=components[0], y=components[1], 
                                 color=self.obs[hue].values.astype('object'))
#                                  color_discrete_sequence=px.colors.qualitative.T10)
                fig.update_layout(template='plotly_white',height=500,width=600,
                                 font=dict(family='sans-serif',size=12),title=ttitle)
#                 fig.update_traces(marker=dict(size=5.5,opacity=0.75,)) # msize
                fig.update_traces(marker_line=dict(width=1, color='DarkSlateGray'))
                if(method is 'PCA'):
                    fig.update_xaxes(title_text=f'{self.__pcalabels[str(components[0])]}')
                    fig.update_yaxes(title_text=f'{self.__pcalabels[str(components[1])]}')
                else:
                    fig.update_xaxes(title_text=f'{method}{str(components[0])}')
                    fig.update_yaxes(title_text=f'{method}{str(components[1])}')
                fig.show()
              
        else:
            print('Data not found in reducedDim')