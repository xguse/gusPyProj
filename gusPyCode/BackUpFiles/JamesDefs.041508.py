def groupByField (listOfTabbedStrings, fieldGroupedBy): 
    """ WARNING!! listOfTabbedStrings will be destroyed!!
    Example:
    fieldGroupedBy = 0
    listOfTabbedStrings = ['a\t1','a\t2','b\t1']
    result-> [[['a','1'],['a','2']],[['b','1']] ]
    """

    # Convert listOfTabbedStrings to a list of lists with 'lowest' list
    # representing a list of original tab seperated fields
    listOfListsByTab = []
    while listOfTabbedStrings != []:
        listOfListsByTab.append(listOfTabbedStrings.pop(0).split('\t'))
        
    exonList = []

    listOfExonLists = []
    
    while listOfListsByTab != 'end':

        # If there is a fresh and clean exonList:
        #     add the first coding region of the first/next gene to exonList
        if exonList == []:
            exonList.append(listOfListsByTab.pop(0))

        # If the next BioMart record list matches the one(s) in exonList:
        #     add it to exonList
        elif listOfListsByTab[0][fieldGroupedBy] == exonList[0][fieldGroupedBy]:
            exonList.append(listOfListsByTab.pop(0))

            # Check to see if you just popped the last record:
            #    - export last exonList
            #    - cull exonList
            #    - set listOfListsByTab to 'end' to stop the loop
            if listOfListsByTab == []:
                listOfExonLists.append(exonList)
                #print len(listOfExonLists), '\n'
                print exonList[0][0]
                exonList = []
                listOfListsByTab = 'end'
        
        # Otherwise append whole exonList to listOfExonLists and clean exonList for next record group
        else:
            listOfExonLists.append(exonList)
            #print len(listOfExonLists), '\n'
            print exonList[0][0]
            exonList = []
            
    print 'The groupByField function produced ',len(listOfExonLists),' groups.\n\n'
    return listOfExonLists

