
def mystery(list, l, r):

    if l > r:
        return -1

    m = (l+r)/2

    if list[m] == m:
        return m

    if list[m]< m:
        return mystery(list, m+1, r)
    else:
        return mystery(list, l, m-1)

print mystery([-2,1,2,4,7,12,15], 0, 6)