from scipy import stats
import csv

def importXLS(filePath):
    reader = csv.reader(open(filePath, 'rU').readlines())
    fileData = []
    for row in reader:
        fileData.append(row)
    return fileData

def collectColumn(tableData,colNum):
    rData = []
    for row in tableData:
        try:
            rData.append(float(row[colNum]))
        except ValueError:
            print "Non-Float found: %s" % (row[colNum])
            
    print rData[:10]
            

f1path = '/Users/biggus/Documents/James/Data/AffyAgapStuff/home-stretch_time-course_expression.csv'
#f2path = ''


f1 = importXLS(f1path)
#f2 = importXLS(f2path)

collectColumn(f1,3)


