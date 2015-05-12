#!/usr/bin/python
import csv
import sys
import string

DIVISOR=10000
list_of_lines_1 = list(csv.reader(open(sys.argv[1], "rb"), delimiter=' '))
list_of_lines_2 = list(csv.reader(open(sys.argv[1], "rb"), delimiter=' '))
list_of_nodes = []
lut = {}
fragments_per_link={}
i=1
for nodes in list_of_lines_1:
    if nodes[0] not in list_of_nodes:
        lut["\""+nodes[0]+"\""] = i
        i=i+1
        list_of_nodes.append(nodes[0])
    if nodes[1] not in list_of_nodes:
        lut["\""+nodes[1]+"\""] = i
        i=i+1
        list_of_nodes.append(nodes[1])
print "*Vertices",len(list_of_nodes)
i=1
for nodes in list_of_nodes:
        print i, "\""+nodes+"\""
        i=i+1

for arg_index in range(1,len(sys.argv)):
    list_of_lines_1 = list(csv.reader(open(sys.argv[arg_index], "rb"), delimiter=' '))
    list_of_lines_2 = list(csv.reader(open(sys.argv[arg_index], "rb"), delimiter=' '))
    
    for columns_1 in list_of_lines_1:
        found = False
        for columns_2 in list_of_lines_2:
            #found matching a-b and b-a in different lines
            if columns_2[1] == columns_1[0] and columns_2[0] == columns_1[1]:
                found = True
                #order them (a-b) with "a" being the lower in alphabetic order
                if columns_1[0] < columns_2[0]:
                    key="\""+columns_1[0]+ "\"--\""+columns_1[1]+"\""
                    if key in fragments_per_link:
                        fragments_per_link[key]+=max(int(columns_1[2])+int(columns_2[2]),1)
                    else:
                        fragments_per_link[key]=max(int(columns_1[2])+int(columns_2[2]),1)
        #didn't find b-a match for a-b.
        if not found:
            #however, it might be existing from previous files
            key1="\""+columns_1[0]+ "\"--\""+columns_1[1]+"\""
            key2="\""+columns_1[1]+ "\"--\""+columns_1[0]+"\""
            min_key=min(key1,key2)
            if min_key not in fragments_per_link:
                fragments_per_link[min_key]=max(int(columns_1[2]),1)
            else:
                fragments_per_link[min_key]+=max(int(columns_1[2]),1)

print "*Edges", len(fragments_per_link)
for key in fragments_per_link:
    print lut[key.split("--")[0]], lut[key.split("--")[1]], int(fragments_per_link[key])
