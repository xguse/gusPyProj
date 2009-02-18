from TAMO import MotifTools,MotifMetrics,seq


motifs = ['CACWCWC',
          'CGGNGGT',
          'TTGACGTK',
          'TGACGTKT',
          'TTTGACRK',
          'WMAACAAA',
          'WGATAAG',
          'STTTTATAC',
          'GCSCGCS']

for i in range(0,len(motifs)):
    motifs[i] = MotifTools.Motif_from_text(motifs[i])
    
cluster = seq.Fasta.file2dict('/Users/biggus/Documents/James/Data/ClusterDefs/TC-Fastas/TC-8.fas')

cKeys = cluster.keys()
cKeys.sort()

result = ['#\t%s' % ('\t'.join(map(lambda l: l.oneletter, motifs)))]

for k in cKeys:
    hit = [k]
    for m in motifs:
        if m.scan(cluster[k])[0]:
            hit.append('1')
        else:
            hit.append('0')
    result.append('%s' % ('\t'.join(hit)))
    
for r in result:
    print r

