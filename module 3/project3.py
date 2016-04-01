"""
project 3
"""

import alg_cluster

def slow_closest_pair(cluster_list):
    """
    from a list of nodes, return distance of the closest pair
    :param cluster_list: list of nodes
    :return: tuple(dist, idx1, idx2), idx1<idx2, dist is the distance between the closest pair
    """
    nodelist = list(cluster_list)
    closest = [float("inf"), -1, -1]
    for idx1 in range(len(nodelist)):
        for idx2 in range(len(nodelist)):
            if idx1 != idx2:
                dist = nodelist[idx1].distance(nodelist[idx2])
                if dist < closest[0]:
                    closest[0] = dist
                    closest[1] = min(idx1, idx2)
                    closest[2] = max(idx1, idx2)

    return tuple(closest)

def fast_closest_pair(cluster_list):
    """
    Compute the distance between the closest pair of clusters in a list (fast)
    Input: cluster_list is list of clusters SORTED such that horizontal positions of their
    centers are in ascending order

    Output: tuple of the form (dist, idx1, idx2) where the centers of the clusters
    cluster_list[idx1] and cluster_list[idx2] have minimum distance dist.
    """

    nodelist = list(cluster_list)

    if len(nodelist) <= 3:
        return slow_closest_pair(nodelist)

    half = len(nodelist)/2
    left = nodelist[0:half]
    right = nodelist[half:]

    left_close = fast_closest_pair(left)
    right_close = fast_closest_pair(right)

    # re-index right_close
    right_close = (right_close[0], right_close[1] + half, right_close[2] + half)

    if left_close[0] < right_close[0]:
        closest = left_close
    else:
        closest = right_close

    mid = 0.5 * (nodelist[half].horiz_center()+nodelist[half-1].horiz_center())
    pair_strip = closest_pair_strip(cluster_list,mid, closest[0])

    if pair_strip[0] < closest[0]:
        closest = pair_strip

    return closest

def closest_pair_strip(cluster_list, horiz_center, half_width):
    """
    :param cluster_list: Cluster objects
    :param horiz_center: horiz_center specifies the horizontal position of the center line for a vertical strip
    :param half_width: half_width specifies the maximal distance of any point in the strip from the center line
    :return: a tuple corresponding to the closest pair of clusters that lie in the specified strip
    """
    closest = [float("inf"), -1, -1]
    nodelist = list(cluster_list)
    s_set=[]

    for each in nodelist:
        if abs(each.horiz_center()-horiz_center)< half_width:
            s_set.append(each)

    s_set.sort(key = lambda cluster: cluster.vert_center())

    if len(nodelist):
        for idx1 in range (0, len(s_set)-1):
            for idx2 in range (idx1+1, len(s_set)):
                dist = s_set[idx1].distance(s_set[idx2])
                if dist < closest[0]:
                    closest[0] = dist
                    closest[1] = min(nodelist.index(s_set[idx1]), nodelist.index(s_set[idx2]))
                    closest[2] = max(nodelist.index(s_set[idx1]), nodelist.index(s_set[idx2]))

    return tuple(closest)

def hierarchical_clustering(cluster_list, num_clusters):
    """
    Compute a hierarchical clustering of a set of clusters
    Note: the function may mutate cluster_list

    Input: List of clusters, integer number of clusters
    Output: List of clusters whose length is num_clusters
    """
    nodelist = list(cluster_list)
    clusters = list(nodelist)
    clusters.sort(key = lambda cluster: cluster.horiz_center())

    while len(clusters) > num_clusters:

        closest = fast_closest_pair(clusters)
        clusters[closest[1]].merge_clusters(clusters[closest[2]])
        clusters.remove(clusters[closest[2]])
        clusters.sort(key = lambda cluster: cluster.horiz_center())

    return clusters

def kmeans_clustering(cluster_list, num_clusters, num_iterations):
    """
    Compute the k-means clustering of a set of clusters
    Note: the function may not mutate cluster_list

    Input: List of clusters, integers number of clusters and number of iterations
    Output: List of clusters whose length is num_clusters
    """

    nodelist = list(cluster_list)
    centers = []
    nodelist_by_pop = list(nodelist)
    nodelist_by_pop.sort(key = lambda cluster: cluster.total_population(), reverse=True)

    # position initial clusters at the location of clusters with largest populations
    for idx in range(num_clusters):
        centers.append(alg_cluster.Cluster(set([]), nodelist_by_pop[idx].horiz_center(), nodelist_by_pop[idx].vert_center(), nodelist_by_pop[idx].total_population(), 0))

    for idx in range(num_iterations):
        # set empty list of clusters with clusters = num_clusters
        results = [alg_cluster.Cluster(set([]), 0, 0, 0, 0) for _ in range(num_clusters)]
        for each in nodelist:
            shortest_dist = float("inf")
            center_position = 0
            for center in centers:
                dist = each.distance(center)
                if dist < shortest_dist:
                    shortest_dist = dist
                    center_position = centers.index(center)

            results[center_position].merge_clusters(each)

        # reset centers
        centers = list(results)

    return results

