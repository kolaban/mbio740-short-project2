from collections import defaultdict
import re

values = defaultdict(int)

blast = open("Sulfurovum_lithotrophicum.go", 'r')
for line in blast:
    line = line.split(';')
    for val in line:
        match = re.search("GO:\d*", val)
        if(match):
            values[match.group()] += 1
blast.close()

values = sorted(values.items(), key=lambda k_v: k_v[1], reverse=True)
goOut = open("go_sample.tsv", 'w')

for item in values:
    goOut.write(str(item[0]) + '\t' + str(item[1]))
    goOut.write('\n')

goOut.close()