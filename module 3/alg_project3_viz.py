"""
Example code for creating and visualizing
cluster of county-based cancer risk data

Note that you must download the file
http://www.codeskulptor.org/#alg_clusters_matplotlib.py
to use the matplotlib version of this code
"""

# Flavor of Python - desktop or CodeSkulptor
DESKTOP = True

import math
import random
import urllib2
import alg_cluster
import project3 as p3
import matplotlib.pyplot as plt

# conditional imports
if DESKTOP:
    #import alg_project3_solution      # desktop project solution
    import project3 as alg_project3_solution
    import alg_clusters_matplotlib
else:
    #import userXX_XXXXXXXX as alg_project3_solution   # CodeSkulptor project solution
    import alg_clusters_simplegui
    import codeskulptor
    codeskulptor.set_timeout(30)


###################################################
# Code to load data tables

# URLs for cancer risk data tables of various sizes
# Numbers indicate number of counties in data table

DIRECTORY = "http://commondatastorage.googleapis.com/codeskulptor-assets/"
DATA_3108_URL = DIRECTORY + "data_clustering/unifiedCancerData_3108.csv"
DATA_896_URL = DIRECTORY + "data_clustering/unifiedCancerData_896.csv"
DATA_290_URL = DIRECTORY + "data_clustering/unifiedCancerData_290.csv"
DATA_111_URL = DIRECTORY + "data_clustering/unifiedCancerData_111.csv"


def load_data_table(data_url):
    """
    Import a table of county-based cancer risk data
    from a csv format file
    """
    data_file = urllib2.urlopen(data_url)
    data = data_file.read()
    data_lines = data.split('\n')
    print "Loaded", len(data_lines), "data points"
    data_tokens = [line.split(',') for line in data_lines]
    return [[tokens[0], float(tokens[1]), float(tokens[2]), int(tokens[3]), float(tokens[4])] 
            for tokens in data_tokens]


############################################################
# Code to create sequential clustering
# Create alphabetical clusters for county data

def sequential_clustering(singleton_list, num_clusters):
    """
    Take a data table and create a list of clusters
    by partitioning the table into clusters based on its ordering
    
    Note that method may return num_clusters or num_clusters + 1 final clusters
    """
    
    cluster_list = []
    cluster_idx = 0
    total_clusters = len(singleton_list)
    cluster_size = float(total_clusters)  / num_clusters
    
    for cluster_idx in range(len(singleton_list)):
        new_cluster = singleton_list[cluster_idx]
        if math.floor(cluster_idx / cluster_size) != \
           math.floor((cluster_idx - 1) / cluster_size):
            cluster_list.append(new_cluster)
        else:
            cluster_list[-1] = cluster_list[-1].merge_clusters(new_cluster)
            
    return cluster_list
                

def compute_distortion(cluster_list, data_table):

    sum = 0
    for each in cluster_list:
        sum += each.cluster_error(data_table)
    return sum



#####################################################################
# Code to load cancer data, compute a clustering and 
# visualize the results
#Q10
# def hierarchical_clustering_list (cluster_list, num_clusters_start, num_clusters_end):
#     """
#     Compute a hierarchical clustering of a set of clusters
#     Note: the function may mutate cluster_list
#
#     Input: List of clusters, integer number of clusters
#     Output: List of clusters whose length is num_clusters
#     """
#
#     dict_clusters ={}
#
#     nodelist = list(cluster_list)
#     clusters = list(nodelist)
#     clusters.sort(key = lambda cluster: cluster.horiz_center())
#
#     count = num_clusters_end
#
#     while len(clusters) > num_clusters_start:
#
#         if len(clusters) == count:
#             dict_clusters[count] = list(clusters)
#             count -= 1
#
#         closest = p3.fast_closest_pair(clusters)
#         clusters[closest[1]].merge_clusters(clusters[closest[2]])
#         clusters.remove(clusters[closest[2]])
#         clusters.sort(key = lambda cluster: cluster.horiz_center())
#
#     return dict_clusters

def run_example(table, method):
    """
    Load a data table, compute a list of clusters and 
    plot a list of clusters

    Set DESKTOP = True/False to use either matplotlib or simplegui
    """
    #data_table = load_data_table(DATA_3108_URL)
    #data_table = load_data_table(DATA_290_URL)
    data_table = load_data_table(table)

    singleton_list = []
    for line in data_table:
        singleton_list.append(alg_cluster.Cluster(set([line[0]]), line[1], line[2], line[3], line[4]))
        
    #cluster_list = sequential_clustering(singleton_list, 15)
    #print "Displaying", len(cluster_list), "sequential clusters"

    cluster_distortion_dict ={}
    start = 20
    end = 6

    count = start

    new_list = list(singleton_list)

    while count >=end:
        if method == 'h_cluster':

            cluster_list = alg_project3_solution.hierarchical_clustering(new_list, count)
            cluster_distortion_dict[count] = compute_distortion(cluster_list, data_table)
            new_list = cluster_list

        elif method == 'k_cluster':

            cluster_list = alg_project3_solution.kmeans_clustering(singleton_list, count, 5)
            cluster_distortion_dict[count] = compute_distortion(cluster_list, data_table)
            #new_list = cluster_list

        count -=1


    #print "Displaying", len(cluster_list), "hierarchical clusters"
    #print "Displaying", len(cluster_list), "hierarchical clusters cluster error"

    #cluster_list = alg_project3_solution.kmeans_clustering(singleton_list, 9, 5)
    #print "Displaying", len(cluster_list), "k-means clusters"
    #print "Displaying", len(cluster_list), "k-means clusters cluster error"

    # draw the clusters using matplotlib or simplegui
    if DESKTOP:
        #alg_clusters_matplotlib.plot_clusters(data_table, cluster_list, False)
        #alg_clusters_matplotlib.plot_clusters(data_table, cluster_list, True)  #add cluster centers
        #print compute_distortion(cluster_list, data_table)
        return cluster_distortion_dict
    else:
        alg_clusters_simplegui.PlotClusters(data_table, cluster_list)   # use toggle in GUI to add cluster centers
    
#print run_example(DATA_111_URL)

table_list = [DATA_111_URL, DATA_290_URL, DATA_896_URL]

distortion_dict_h = run_example(DATA_896_URL, 'h_cluster')
distortion_dict_k = run_example(DATA_896_URL, 'k_cluster')

distortion_list_h = []
distortion_list_k = []
n = list(range(6,21))

for each in n:
    distortion_list_h.append(distortion_dict_h[each]/1000000000)
    distortion_list_k.append(distortion_dict_k[each]/1000000000)

#plot graph
plt.plot(n, distortion_list_h)
plt.plot(n, distortion_list_k)
plt.legend(['hierarchical clustering', 'k-means clustering'], loc='upper right')
plt.title('Distortion by hierarchical clustering \n and k-means clustering of 896 county data set')
plt.ylabel('Distortion [ * 10^(-9) ]')
plt.xlabel('Number of clusters')
plt.show()