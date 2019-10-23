name = 'data/pemor.txt'
fp = open(name,'r')
s = "\n".join(fp.readlines())
fp.close()
print(s)