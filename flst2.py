#!/usr/bin/env python2
# mistake started with file list, but really I should be reading a file
import os, sys, re

fn = sys.argv[1]
renli = [] # rename list
with open(fn) as f:
    for ln in f:
        renli.append([nm for nm in ln.strip().split(' ')])

rx1=re.compile(r'(.+)\.jpg')

for i in renli:
    m=rx1.match(i[1])
    if(m):
       fpn=m.group(1) +'.png'
       nfpn=i[0] +'.png'
       # print "os.rename("+fpn+", "+nfpn+")"
       if os.path.isfile(fpn):
           os.rename(fpn, nfpn)
