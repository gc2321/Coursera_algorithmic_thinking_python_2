import matplotlib.pyplot as plt
import math
import random
import project3 as p3
import timeit
import alg_cluster

def gen_random_clusters(num_clusters):
    """
    generate random clusters
    :param num_clusters:  number of clusters
    :return:  list of clusters
    """
    cluster_list =[]

    for i in range (num_clusters):
        numbers = []
        for j in range(2):
            numbers.append(random.uniform(0, 1)*random.choice([-1,1]))

        cluster_list.append([numbers[0], numbers[1]])

    return cluster_list

#print gen_random_clusters(3)

#Q1 Comparing slow_closest_pair and fast_closest_pair
'''
n = list(range(2, 201))

time_slow=[]
time_fast=[]

for each in n:
    list = gen_random_clusters(each)
    cluster_list =[]

    for node in list:

        cluster_list.append(alg_cluster.Cluster(set([]), node[0], node[1], 0, 0))

    start1 = timeit.default_timer()
    p3.slow_closest_pair(cluster_list)
    stop1 = timeit.default_timer() - start1
    time_slow.append(stop1)

    start2 = timeit.default_timer()
    p3.fast_closest_pair(cluster_list)
    stop2 = timeit.default_timer() - start2
    time_fast.append(stop2)

#plot graph
plt.plot(n,time_slow)
plt.plot(n,time_fast)
plt.legend(['slow_closest_pair', 'fast_closest_pair'], loc='upper left')
plt.title('Comparing efficiency of slow_closest_pair and fast_closest_pair')
plt.ylabel('Time (s)')
plt.xlabel('Number of points')
plt.show()
'''


