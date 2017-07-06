raw=open("H-Bond/H_Bond_O_index_during_transition")
out=open("H-Bond/H_Bond_O_index_during_transition_read",'wb')
m=sorted(list(raw.readlines()))
print(list(m[0].split()))
for n in range(len(m)):
    line=[t for t in m[n].split(' ') if t!='']
    for i in range(4,len(line)):
        o=[]
        o.append(line[0])
        o.append(line[1])
        o.append(line[2])
        o.append(line[i].rstrip('\n'))
        out.write('  '.join(o))
        out.write('\n')
raw.close()
out.close()