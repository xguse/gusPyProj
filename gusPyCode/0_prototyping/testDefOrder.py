
def defA(x,y):
    return x+y

def defB(x,y):
    return defA(x,y)

if __name__ == "__main__":
    result = defB(2,4) 
    
    print result