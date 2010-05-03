def norm2ReadSets(list1, list2):
    """Takes two lists -> Returns two Normalized Lists.
    Normalizes Read counts in each index based
    on the list with the smallest total (this 
    list is not changed). Indexes in list with 
    larger total are multiplied by 
    (float(sum(smaller))/sum(larger)) to scale the
    read counts to those expected if both total reads
    equaled sum(smaller)."""
    if not len(list1) == len(list2):
        raise Exception, "ERROR: len(list1) != len(list2)."
    
    # ++ float-ify lists ++
    for i in range(len(list1)):
        list1[i] = float(list1[i])
        list2[i] = float(list2[i]) 
    
    totL1,totL2 = sum(list1),sum(list2)
    
    normd1,normd2 = [],[]
    
    if totL1 < totL2:
        for i in range(len(list1)):
            normd1.append(list1[i])
            normd2.append(list2[i]*(totL1/totL2))
    else:
        for i in range(len(list1)):
            normd2.append(list2[i])
            normd1.append(list1[i]*(totL2/totL1))
    
    return normd1,normd2