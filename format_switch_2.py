import operator
raw=open("jump/cutoff_1a3/"
   "transition_path_stat")
out=open("jump/cutoff_1a3/"
   "transiton_stat_fastinput",'wb')
m=sorted(map(lambda x:map(int,[t for t
   in x.split()if t!='']),raw.readlines()),
   key=operator.itemgetter(1,2,3))
for a in m: out.write('  '.join(map(str,a)))
raw.close()
out.close()