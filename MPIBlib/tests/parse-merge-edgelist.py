#!/usr/bin/python
import csv
import sys
import string
import networkx as nx

DIVISOR=10000
list_of_lines_1 = list(csv.reader(open(sys.argv[1], "rb"), delimiter=' '))
list_of_lines_2 = list(csv.reader(open(sys.argv[1], "rb"), delimiter=' '))
list_of_nodes = []
fragments_per_link={}
for nodes in list_of_lines_1:
    if nodes[0] not in list_of_nodes:
        list_of_nodes.append(nodes[0])
    if nodes[1] not in list_of_nodes:
        list_of_nodes.append(nodes[1])

for nodes in list_of_nodes:
    if -1 != string.find(nodes, "172.16.113."):
        print "\"%s\" 1" % nodes
    if -1 != string.find(nodes,"172.16.1.") or -1 != string.find(nodes,"172.16.2."):
        print "\"%s\" 2" % nodes
    if -1 != string.find(nodes,"172.16.16."):
        print "\"%s\" 3" % nodes
    if -1 != string.find(nodes,"10.69.0.") or -1 != string.find(nodes,"10.69.1."):
        print "\"%s\" 4" % nodes

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
                    #print "\"%s\"--\"%s\"[len=%lf];" % (columns_1[0], columns_1[1], DIVISOR/max(int(columns_1[2]) + int(columns_2[2]),1))
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

G=nx.Graph()
for key in fragments_per_link:
    #print "%s[len=%lf];" % (key,max(1.,DIVISOR/fragments_per_link[key]))
    #print "%s %lf" % (key,max(1.,DIVISOR/fragments_per_link[key]))
    G.add_edge(key.split("--")[0], key.split("--")[1], weight=int(fragments_per_link[key]))
    #print "%s[weight=%lf];" % (key,max(1.,fragments_per_link[key]))
nx.write_weighted_edgelist(G, "edgelist")
