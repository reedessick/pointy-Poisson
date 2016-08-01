#!/usr/bin/python
usage = "lists2UnionAndIntersection.py list1 list2"
description = "takes the union and intersection of two lists. The union is written to stdout and the intersection is written to stderr"

#-------------------------------------------------

import sys

from optparse import OptionParser

#-------------------------------------------------

parser = OptionParser(usage=usage, description=description)

parser.add_option('-v', '--verbose', default=False, action="store_true")

opts, args = parser.parse_args()

if len(args)!=2:
    raise ValueError('please supply exactly 2 input arguments\n%s'%usage)

#-------------------------------------------------

### read in lists

if opts.verbose:
    print "reading list1 from : %s"%args[0]
file_obj = open(args[0], 'r')
list1 = sorted([_.strip() for _ in file_obj.readlines()])
file_obj.close()
N1 = len(list1)
if opts.verbose:
    print "found %d items"%N1 

if opts.verbose:
    print "reading list2 from : %s"%args[1]
file_obj = open(args[1], 'r')
list2 = sorted([_.strip() for _ in file_obj.readlines()])
file_obj.close()
N2 = len(list2)
if opts.verbose:
    print "found %d items"%N2 

#-------------------------------------------------

### iterate through lists, which we know are sorted
i1 = 0
i2 = 0
u = 0
i = 0

while (i1<N1) and (i2<N2):
    item1 = list1[i1]
    item2 = list2[i2]

    if item1 == item2: ### items match -> interesection
        print >> sys.stdout, item1 ### add to union
        print >> sys.stderr, item1 ### add to intersection
        i1 += 1
        i2 += 1
        u += 1
        i += 1

    elif item1 < item2: ### add item1 to union, increment i1
        print >> sys.stdout, item1 ### add to union
        i1 +=1
        u += 1
    else: ### add item2 to union, increment i2
        print >> sys.stdout, item2 ### add to union
        i2 += 1
        u += 1

### finish up 
while (i1<N1):
    print >> sys.stdout, list1[i1]
    i1 += 1
    u += 1
while (i2<N2):
    print >> sys.stdout, list2[i2]
    i2 += 1
    u += 1

if opts.verbose:
    print "found %d items in the union"%u
    print "found %d items in the intersection"%i
