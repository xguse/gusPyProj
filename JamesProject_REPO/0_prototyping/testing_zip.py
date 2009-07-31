import pprint

l = [[['x'],['m1'],['m2'],['m3']],
     [['g1'],[1,2],[3,4],[5,6]],
     [['g2'],[7,8],[9,10],[11,12]]]

zl = zip(*l)

pprint.pprint(l)
pprint.pprint(zl)