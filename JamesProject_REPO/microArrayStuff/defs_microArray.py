

class geneVectors:
    
    
    def __init__(self, listOfVectors):
        import math
        from decimal import Decimal
        
        self.vectors = {}
        
        for each in listOfVectors:
            fields = each.split('\t')
            listOfValues = fields[1:]
            for i in range(len(listOfValues)):
                listOfValues[i] = Decimal(listOfValues[i])
            self.vectors[fields[0]] = listOfValues
            
    