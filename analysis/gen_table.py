from bioservices.uniprot import UniProt

from collections import defaultdict
import re

blast = open("staph_out.txt", 'r')
table = defaultdict(int)

for line in blast:
    match = re.search("\A> \w", line)
    if(match):
        match = re.search("\w{,5}_\w{,5}\s", line[2:])
        if(match):
            table[match.group()] += 1
blast.close()

print(len(table))
blastOut = open("blast.tsv", 'w')
u = UniProt(verbose=False)
goTable = defaultdict(int)

i = 0
for item in table:
    blastOut.write(item + '\t' + str(table[item]) + '\n')
    i += 1
    if(i % 25 == 0):
        print(i)
    value = u.search(item, columns="go")
    value = value.split(';')
    for val in value:
        match = re.search("GO:\d*", val)
        if(match):
            goTable[match.group()] += (1 * table[item])
blastOut.close()
goTable = sorted(goTable.items(), key=lambda v: v[1], reverse=True)
goOut = open("go.tsv", 'w')

for item in goTable:
    goOut.write(str(item[0]) + '\t' + str(item[1]))
    goOut.write('\n')

goOut.close()