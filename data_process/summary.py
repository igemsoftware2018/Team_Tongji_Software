import re
file = open('reaction.csv')
pattern = 'C\d\d\d\d\d'
pairs_num = 0

index = 0
for line in file.readlines():
    index += 1
    line = ''.join(line.strip().split(',')[2:])
    pairs = re.findall(pattern, line)
    if pairs == []:
                print(index)
    pairs_num += len(pairs)/2


print(pairs_num)
