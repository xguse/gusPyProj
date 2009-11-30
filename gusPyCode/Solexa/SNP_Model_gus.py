import sys
import math



class ReadsFeature(object):
    def __init__(self, Coordi, Vari_Offset, QualStr):
            
        self.Coordi      = Coordi
        self.Vari_Offset = Vari_Offset
        self.Quals       = []
        
        # initialize self.Quals:
        self._qualStr2probList(QualStr)
    
    def _qualStr2probList(self, QualStr):
        """Internal method to initialize self.Quals as list of probabilities (that
        the nucleotide was miscalled due to machine error) from the QualStr
        
        QualStr => abaa_b`bab`abaaaa`aaaaaaaaaaa^``aY`aaba`
        """
        for char in QualStr:
            self.Quals.append(qual2prob(char))
            
    def qualTest(self,index,thresh=0.001):
        """Decide whether supplied nucleotide index in read passes quality threshold.
        """
        
        pass


class Window(object):
    startPoint=0
    endPoint=39
    numReads=1
    numVariation=0

def lineScan(lineStr):
    #print "Hello!"
    lineList = lineStr.strip('\n').split('\t')  # strip trailing NewLine and split on tabs
    
    Coordi  = lineList[12]
    Dire    = lineList[13]
    Varia   = lineList[14]
    QualStr = lineList[9]  # example => abaa_b`bab`abaaaa`aaaaaaaaaaa^``aY`aaba`
    
    #print line
    #print Coordi,Dire,Varia
    return Coordi,Dire,Varia,QualStr




def FeatureExtract(Coordi,Dire,Varia,QualStr):
    """Converts string vals in call signature to correct python data types
    and constructs the ReadsFeature obj.
    """
    #1st step: to cut the feature into int & letter list
    CutSet=[]
    frag=''
    for i in range(len(Varia)):
        if Varia[i]>='0' and Varia[i]<='9':
            frag=frag+Varia[i]
            if i==len(Varia)-1:
                CutSet.append(int(frag))
                frag=''
        else:
            if frag!='':
                gap=int(frag)
                CutSet.append(gap)
                frag=''
            CutSet.append(Varia[i])
    #return CutSet #debug test: OK!

    #2nd step: to construct the abnormal base pair offset list
    Offset=[]
    index=0;
    for item in CutSet:
        if type(item)==int:
            index=index+item
        else:
            Offset.append(index)
            index=index+1
    #return Offset  #debug test: OK!

    #3rd step: to construct the ReadsFeature data structure
    ##read=ReadsFeature()
    #print Coordi,Coordi.isdigit()
    Coordi=int(Coordi)
    if Dire!="F":
        for index in range(len(Offset)):
            Offset[index]=39-Offset[index]
        Offset.reverse()
    ##read.Vari_Offset=Offset
    read = ReadsFeature(Coordi,Offset,QualStr) # I changed ReadsFeature class to include a constructor to handle the QualStr conversion.
    return read

class SNP_Assess:
    Coordi=0
    NumCover=0
    NumVarCover=0
    VarRatio=0
    def __init__(self, Coor, Cover, Var, Ratio):
        self.Coordi=Coor
        self.NumCover=Cover
        self.NumVarCover=Var
        self.VarRatio=Ratio
    def fout(self,file):
        info=str(self.Coordi)+"  "+str(self.NumCover)+"  "+str(self.NumVarCover)+"  "+str(self.VarRatio)+'\n'
        file.write(info)
    def cout(self):
        print self.Coordi,self.NumCover,self.NumVarCover,self.VarRatio

def BufferUpdate(read,file):
    #if the buffer is empty when entering this function, it 
    #should be at the begining of the scan
    if Buffer==[]:
        for i0 in range(40):
            #-----Gus: THIS IS THE 1st PLACE WHERE YOU WANT TO USE QUALITY SCORE---
            if read.Vari_Offset.count(i0)>0:
                baseTmp1=SNP_Assess(read.Coordi+i0,1,1,0)
            else:
                baseTmp1=SNP_Assess(read.Coordi+i0,1,0,0)
            #----------------------------------------------------------------------
            Buffer.append(baseTmp1)
        return Buffer
    
    CurrentBase=SNP_Assess(0,0,0,0)
    while Buffer!=[] and Buffer[0].Coordi<read.Coordi:
        #>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<
        print Buffer[0].Coordi,Buffer[0].NumCover,Buffer[0].NumVarCover,Buffer[0].VarRatio
        Buffer[0].fout(file)
        CurrentBase=Buffer.pop(0)
    
    #this read's coordinate is beyond the buffer's range so 
    #the buffer is popped until empty
    if Buffer==[]:
        #record the gap
        gapIndex=CurrentBase.Coordi+1
        while gapIndex<read.Coordi:
            #>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<
            print gapIndex,0,0,0
            baseTmp3=SNP_Assess(gapIndex,0,0,0)
            baseTmp3.fout(file)
            gapIndex=gapIndex+1
        #rebuild the buffer
        for i1 in range(40):
            #-----Gus: THIS IS THE 2nd PLACE WHERE YOU WANT TO USE QUALITY SCORE---
            if read.Vari_Offset.count(i1)>0:
                baseTmp2=SNP_Assess(read.Coordi+i1,1,1,0)
            else:
                baseTmp2=SNP_Assess(read.Coordi+i1,1,0,0)
            #----------------------------------------------------------------------
            Buffer.append(baseTmp2)
        return Buffer

    #this read's coordinate is within the buffer's range so
    #the overlap area is updated while the new area is expended
    if Buffer[0].Coordi==read.Coordi:
        for i in range(40):
            if i<len(Buffer):
                #-----Gus: THIS IS THE 3rd PLACE WHERE YOU WANT TO USE QUALITY SCORE---
                Buffer[i].NumCover=Buffer[i].NumCover+1
                if read.Vari_Offset.count(i)>0:
                    Buffer[i].NumVarCover=Buffer[i].NumVarCover+1
                #----------------------------------------------------------------------
            else:
                #-----Gus: THIS IS THE 4th PLACE WHERE YOU WANT TO USE QUALITY SCORE---
                if read.Vari_Offset.count(i)>0:
                    baseAppen=SNP_Assess(read.Coordi+i,1,1,0)
                else:
                    baseAppen=SNP_Assess(read.Coordi+i,1,0,0)
                #----------------------------------------------------------------------
                Buffer.append(baseAppen)
        return Buffer


def BufferClear(Buffer,file):
    while Buffer!=[]:
        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<
        print Buffer[0].Coordi,Buffer[0].NumCover,Buffer[0].NumVarCover,Buffer[0].VarRatio
        Buffer[0].fout(file)
        Buffer.pop(0)
    return Buffer

def qual2prob(qualChar):
    """Takes a single Phred quality score character and returns a the probability
    that the corresponding base was called incorrectly.
    Q_phred = -10 * log10(p)
    p = 10^(Q_phred/-10)
    """
    q = ord(qualChar)-64
    return pow(10,(float(q)/-10))

if __name__=="__main__":
    
    # Soft test of command line input
    usage = 'python %s inFile outFile [qualThresh]\nDefault qualThresh = 0.001' % (sys.argv[0].split('\t')[-1])
    if len(sys.argv) < 3:
        print usage
        exit(1)
        
    file=open(sys.argv[1],'r')
    landscape=open(sys.argv[2],'w')
    
    if sys.argv[3]:
        qualThresh = float(sys.argv[3])
    else:
        qualThresh = 0.001  # Should find out what best defualt value should be here!

    Buffer=[]

    #Coordi=''
    #Dire=''
    #Varia=''
    #read=ReadsFeature()

    line=file.readline()
    while line!='' and line!='\n':
        Coordi,Dire,Varia,QualStr=lineScan(line)
        read=FeatureExtract(Coordi,Dire,Varia,QualStr)
        Buffer=BufferUpdate(read,landscape)
        line=file.readline()
    Buffer=BufferClear(Buffer,landscape)
    file.close()