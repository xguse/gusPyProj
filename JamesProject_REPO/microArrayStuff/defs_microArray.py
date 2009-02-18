

class geneVectors:
    """Class to serve as gene expression data type.  Takes a list of '\t' separated strings in
    the form: probe/gene;expression value 1; ...; expression value N.  May eventually be 
    expanded to incorperate column titles."""
    
    def __init__(self, listOfVectors):

        
        self.vectors = {}
        
        for each in listOfVectors:
            fields = each.split('\t')
            listOfValues = fields[1:]
            for i in range(len(listOfValues)):
                listOfValues[i] = float(listOfValues[i])
            self.vectors[fields[0]] = listOfValues
            

            
    