from math import log

def f(x):
	return log((x)**56) + log((1-x)**44)
	
data = {}

for x in range(0,100):
	x1 = x/100.0
	data[f(x1)] = x1

k = max(data.keys())	
	
print 'MLE = %s' % (data[k])