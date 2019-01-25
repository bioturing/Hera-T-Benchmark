from sys import argv

f = open(argv[1], "r")
o = open(argv[2], "w")
r = 0
for l in f:
    r += 1
    l = l[:-1].split(",")
    for c in range(0, len(l) - 1):
        if float(l[c]) != 0:
            o.write("%u %u %u\n" % (c + 1, r, round(float(l[c]))))
