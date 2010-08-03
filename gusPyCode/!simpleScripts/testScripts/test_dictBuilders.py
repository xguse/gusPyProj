from collections import defaultdict
def oldDict(lst):
    d = {}
    for l in lst:
        k = '%s_%s' % (l[0],l[1])
        if k in d:
            d[k]+=1
        else:
            d[k]=0
            d[k]+=1

def tryDict(lst):
    d = {}
    for l in lst:
        k = '%s_%s' % (l[0],l[1])
        try:
            d[k]+=1
        except KeyError:
            d[k]=0
            d[k]+=1

def defDict(lst):
    