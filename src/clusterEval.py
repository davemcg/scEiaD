
#%%
import numpy as np
import pandas as pd
import hnswlib
from scipy.sparse import csr_matrix
import igraph as ig
import leidenalg
import time
import os 
import sys
import louvain
import parc
from sklearn import datasets
import parc
import os 
os.chdir('/data/swamyvs/scEiad_subcelltype')
#%%
"""
Need to break into the following pieces:
KNN generation
Pruning 
leiden
louvain
merging
"""
class clusterEval:
    def __init__(self, data, knn, nthreads):
        self.knn = knn
        self.data = data
        self.nthreads = nthreads
        self.random_seed=42
        self.dist_std_local = 3
        self.time_smallpop = 15

    def make_csrmatrix_noselfloop(self, neighbor_array, distance_array, local_pruning):
        # neighbor array not listed in in any order of proximity
        row_list = []
        col_list = []
        weight_list = []
        n_neighbors = neighbor_array.shape[1]
        n_cells = neighbor_array.shape[0]
        rowi = 0
        discard_count = 0
        # locally prune based on (squared) l2 distance
        if local_pruning:
            print('commencing local pruning based on Euclidean distance metric at',
                  self.dist_std_local, 's.dev above mean')
            distance_array = distance_array + 0.1
            for row in neighbor_array:
                distlist = distance_array[rowi, :]
                to_keep = np.where(distlist < np.mean(
                    distlist) + self.dist_std_local * np.std(distlist))[0]  # 0*std
                updated_nn_ind = row[np.ix_(to_keep)]
                updated_nn_weights = distlist[np.ix_(to_keep)]
                discard_count = discard_count + (n_neighbors - len(to_keep))
                for ik in range(len(updated_nn_ind)):
                    if rowi != row[ik]:  # remove self-loops
                        row_list.append(rowi)
                        col_list.append(updated_nn_ind[ik])
                        dist = np.sqrt(updated_nn_weights[ik])
                        weight_list.append(1/(dist+0.1))
                rowi = rowi + 1
        else:  # dont prune based on distance
            row_list.extend(list(np.transpose(
                np.ones((n_neighbors, n_cells)) * range(0, n_cells)).flatten()))
            col_list = neighbor_array.flatten().tolist()
            weight_list = (1. / (distance_array.flatten() + 0.1)).tolist()

        csr_graph = csr_matrix((np.array(weight_list), (np.array(row_list), np.array(col_list))),
                               shape=(n_cells, n_cells))
        return csr_graph

    def buildNeighborGraph(self, nn_space, ef_construction, local_pruning, global_pruning, jac_std_global):
                # ef always should be >K. higher ef, more accurate query
        ef_query = max(100, self.knn + 1)
        num_dims = self.data.shape[1]
        n_elements = self.data.shape[0]
        p = hnswlib.Index(space=nn_space, dim=num_dims)
        p.set_num_threads(self.nthreads)
        if (num_dims > 30) & (n_elements <= 50000):
            p.init_index(max_elements=n_elements, ef_construction=ef_construction,
                             M=48)  # good for scRNA seq where dimensionality is high
        else:
            p.init_index(max_elements=n_elements, ef_construction=ef_construction, M=24 )
        
        p.add_items(self.data)
        p.set_ef(ef_query)
        neighbor_array, distance_array = p.knn_query(self.data, self.knn)
        csr_array = self.make_csrmatrix_noselfloop(neighbor_array, distance_array, local_pruning)
        self.neighbor_array = neighbor_array
        sources, targets = csr_array.nonzero()

        edgelist = list(zip(sources, targets))
        nn_graph = ig.Graph(edgelist, edge_attrs={'weight': csr_array.data.tolist()})
        if global_pruning:
            print("Running global pruning")
            n_elements = self.data.shape[0]
            edgelist_copy = edgelist.copy()
            sim_list = nn_graph.similarity_jaccard(pairs=edgelist_copy)
            sim_list_array = np.asarray(sim_list)
            edge_list_copy_array = np.asarray(edgelist_copy)
            if jac_std_global == 'median':
                threshold = np.median(sim_list)
            else:
                threshold = np.mean(sim_list) - jac_std_global * np.std(sim_list)
            strong_locs = np.where(sim_list_array > threshold)[0]
            new_edgelist = list(edge_list_copy_array[strong_locs])
            sim_list_new = list(sim_list_array[strong_locs])

            graph = ig.Graph(n=n_elements, edges=list(new_edgelist),
                            edge_attrs={'weight': sim_list_new})
            # print('Share of edges kept after Global Pruning %.2f' % (len(strong_locs) / len(sim_list)), '%')
        self.nn_graph = graph
        return


    def run_perturbation(self):
        ## perform pertrubation similar to the sanes bipolar paper
        print("Poortoorb ing")
        graph = self.nn_graph
        ### add randomly add and remove edges to graph
        sub_sample_size = int(len(graph.es) * .05)
        edges_to_remove = np.random.choice(
            list(range(len(graph.es))), size=sub_sample_size, replace=False)
        edges_to_add = [tuple(np.random.choice(list(range(len(graph.vs))),
                                               size=2,
                                               replace=False))
                        for _ in range(sub_sample_size)]
        graph.delete_edges(edges_to_remove)
        graph.add_edges(edges_to_add)
       ### add noise to edge weights
        j = 0
        for i in range(len(graph.es)):
            if graph.es[i]['weight'] is None:
                graph.es[i]['weight'] = np.random.uniform(.6, 1.66)
            else:
                graph.es[i]['weight'] = graph.es[i]['weight'] * \
                    np.random.uniform(.6, 1.66)
        self.nn_graph = graph
        return

    def run_leiden(self, vertex_partition_method, n_iter, resolution, jac_weighted_edges=None):
        self.nn_graph.simplify(combine_edges='sum')
        n_elements = self.data.shape[0]
        partition = leidenalg.find_partition(self.nn_graph, vertex_partition_method,
                                             weights=jac_weighted_edges, n_iterations=n_iter, 
                                             resolution_parameter = resolution, seed=self.random_seed)
        labels = np.asarray(partition.membership)
        return labels
    def run_louvain(self, vertex_partition_method, resolution,  jac_weighted_edges=None):
        self.nn_graph.simplify(combine_edges='sum')
        n_elements = self.data.shape[0]
        partition = louvain.find_partition(self.nn_graph, vertex_partition_method,
                                             weights=jac_weighted_edges, 
                                             resolution_parameter = resolution, seed=self.random_seed)
        labels = np.asarray(partition.membership)
        
        return labels

    def merge_singletons(self, labels, small_pop ):
        n_elements = self.data.shape[0]
        
        labels = np.reshape(labels, ( n_elements, 1))
        dummy, labels = np.unique(
            list(labels.flatten()), return_inverse=True)
        small_pop_exist = True
        time_smallpop_start = time.time()
        while (small_pop_exist == True) & ((time.time() - time_smallpop_start) < self.time_smallpop):
            small_pop_list = []
            small_cluster_list = []
            small_pop_exist = False
            
            for cluster in set(list(labels.flatten())):
                population = len(np.where(labels == cluster)[0])
                if population < small_pop:
                    small_pop_exist = True
                    small_pop_list.append(list(np.where(labels == cluster)[0]))
                    small_cluster_list.append(cluster)
            labels = reassign_small_cluster_cells(labels, small_pop_list, small_cluster_list, self.neighbor_array)
        return labels 


def reassign_small_cluster_cells(labels, small_pop_list, small_cluster_list, neighbor_array):
    for small_cluster in small_pop_list:
        for single_cell in small_cluster:
            old_neighbors = neighbor_array[single_cell]
            group_of_old_neighbors = labels[old_neighbors]
            group_of_old_neighbors = list(group_of_old_neighbors.flatten())
            available_neighbours = set(
                group_of_old_neighbors) - set(small_cluster_list)
            if len(available_neighbours) > 0:
                available_neighbours_list = [value for value in group_of_old_neighbors if
                                             value in list(available_neighbours)]
                best_group = max(available_neighbours_list,
                                 key=available_neighbours_list.count)
            else:
                best_group = max(set(group_of_old_neighbors),
                                 key=group_of_old_neighbors.count)
            labels[single_cell] = best_group
    return labels

        
        
#%%

        
#%%

ex_data = pd.read_csv('testing/amacrin_mm_scvi_dim.csv', index_col=0)
# parc1 = parc.PARC(ex_data,
#              jac_std_global='median',
#              random_seed=1,
#              num_threads=1,
#              small_pop=10,  # disable small pop
#              knn=30,
#              jac_weighted_edges=True,
#              partition_type="RBConfigurationVertexPartition",
#              resolution_parameter=.9)
# parc1.run_PARC()


# %%
clu_obj = clusterEval(data=ex_data, knn=30,  nthreads=1)
clu_obj.buildNeighborGraph(nn_space='l2', ef_construction=150,local_pruning=True, global_pruning=True, jac_std_global ='median')
labels_leiden_clueval = clu_obj.run_leiden(
    vertex_partition_method=leidenalg.RBConfigurationVertexPartition,
    n_iter=5,
    resolution=.9,
    jac_weighted_edges='weight'
)
#%%
# print('running_louvain')
# labels_louvain_clueval = clu_obj.run_louvain(
#     vertex_partition_method=louvain.RBConfigurationVertexPartition,
#     resolution=1.5,
#     jac_weighted_edges='weight'
# )

new_labels = clu_obj.merge_singletons(labels_leiden_clueval, 10)



# %%


# %%
PARC_labels_leiden = labels_leiden_clueval
dummy, PARC_labels_leiden = np.unique(list(PARC_labels_leiden.flatten()), return_inverse=True)
small_pop_list = []
small_cluster_list = []
small_pop_exist = False

# %%
small_pop=10
neighbor_array = clu_obj.neighbor_array
for cluster in set(PARC_labels_leiden):
    population = len(np.where(PARC_labels_leiden == cluster)[0]) # get cluster counts 

    if population < small_pop:  # 10
        small_pop_exist = True

        small_pop_list.append(list(np.where(PARC_labels_leiden == cluster)[0]))
        small_cluster_list.append(cluster)
# %%


def reassign_small_cluster_cells(labels, small_pop_list, small_cluster_list, neighbor_array):
    for small_cluster in small_pop_list:
        for single_cell in small_cluster:
            old_neighbors = neighbor_array[single_cell]
            group_of_old_neighbors = labels[old_neighbors]
            group_of_old_neighbors = list(group_of_old_neighbors.flatten())
            available_neighbours = set(
                group_of_old_neighbors) - set(small_cluster_list)
            if len(available_neighbours) > 0:
                available_neighbours_list = [value for value in group_of_old_neighbors if
                                             value in list(available_neighbours)]
                best_group = max(available_neighbours_list,
                                 key=available_neighbours_list.count)
            else:
                best_group = max(set(group_of_old_neighbors),
                                 key=group_of_old_neighbors.count)
            labels[single_cell] = best_group
    return labels
