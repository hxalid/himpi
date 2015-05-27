#! /usr/bin/env python
import time
import os

timestr = time.strftime("%Y%m%d-%H%M%S")

P_MAX = 8
P_MIN = 4
BIN_CMD = ' ./test_hbcast'

#My own bcast implementations
#CMND = 'mpirun -n %s  -hostfile my_nodes_graphene --mca btl_tcp_if_exclude eth0  --mca btl openib,sm,self --mca pml ^cm  --mca plm_rsh_agent "ssh -q -o StrictHostKeyChecking=no" '
#ARGS = ' -m 256 -M 16777216 -r 1 -i 30 -s 0 -a %s'
#CMD = CMND + BIN_CMD + ARGS


#for p in [128, 64, 32, 16, 8]:
#    for alg in [7, 6]:#, 5, 3, 2, 1, 0]:
#        outFileName = " > hbcast_logs/hbcast_test_p%s_alg%s.log_%s" % (p, alg, timestr)
#        cmdStr = CMD % (p, alg) + outFileName
#        print "Executing: %s" % (cmdStr) 
#        os.system(cmdStr)


###########################################################

# Disable IB, use Ethernet, multicore, bind to core

CMND = 'mpirun -n %s  -hostfile my_nodes_graphene_p%s --bind-to-core  --mca coll_tuned_use_dynamic_rules 1 --mca coll_tuned_bcast_algorithm %s  --mca btl_tcp_if_exclude ib0,lo,myri0 --mca btl self,sm,tcp --mca pml ^cm  --mca plm_rsh_agent "ssh -q -o StrictHostKeyChecking=no" '
ARGS = ' -m %s -M %s -r 1 -i 30 -s 0 -a 4 -d 2'
CMD = CMND + BIN_CMD + ARGS


for p in [512, 256, 128, 64, 32, 16]:
    for alg in [6, 5, 4, 3, 2, 1, 0]:
        for msg in [2**x for x in range(10, 25)]:
            outFileName = " > hbcast_logs_03012015/hbcast_tuned_eth_multicore_p%s_alg%s_msg%s.log_%s" % (p, alg, msg, timestr)
            cmdStr = CMD % (p, p, alg, msg, msg) + outFileName
#        print "Executing: %s" % (cmdStr)
#        os.system(cmdStr)


#Disable Ethernet, use IB, multicore, bind to core

CMND = 'mpirun -n %s  -hostfile my_nodes_graphene_p%s --bind-to-core  --mca coll_tuned_use_dynamic_rules 1 --mca coll_tuned_bcast_algorithm %s  --mca coll_tuned_bcast_algorithm_segmentsize %s  --mca btl_tcp_if_exclude eth0  --mca btl openib,sm,self --mca pml ^cm  --mca plm_rsh_agent "ssh -q -o StrictHostKeyChecking=no" '
ARGS = ' -m %s -M %s -r 1 -i 15 -s 0 -a 4 -d 2'
CMD = CMND + BIN_CMD + ARGS


for p in [512, 256, 128, 64, 32, 16]:
    for alg in [3, 2, 0]:
	for msg in [2**x for x in range(10, 25)]:
	    for seg in [0, 512, 1024, 16384, 32768, 65536]:
                outFileName = " > hbcast_logs_15022015/hbcast_multicore_ib_p%s_alg%s_msg%s_seg%s.log_%s" % (p, alg, msg, seg, timestr)
                cmdStr = CMD % (p, p/4, alg, seg, msg, msg) + outFileName
                print "Executing: %s" % (cmdStr)
		sys.exit(1)
                os.system(cmdStr)






# Disable IB, use Ethernet

CMND = 'mpirun -n %s  -hostfile my_nodes_graphene  --mca coll_tuned_use_dynamic_rules 1 --mca coll_tuned_bcast_algorithm %s  --mca btl_tcp_if_exclude ib0,lo,myri0 --mca btl self,sm,tcp --mca pml ^cm  --mca plm_rsh_agent "ssh -q -o StrictHostKeyChecking=no" '
ARGS = ' -m %s -M %s -r 1 -i 30 -s 0 -a 4 -d 2'
CMD = CMND + BIN_CMD + ARGS


for p in [128, 64, 32, 16, 8]:
    for alg in [3, 2, 0]:
	for msg in [2**x for x in range(10, 25)]:
            outFileName = " > hbcast_logs_03012015/hbcast_tuned_eth_p%s_alg%s_msg%s.log_%s" % (p, alg, msg, timestr)
            cmdStr = CMD % (p, alg, msg, msg) + outFileName
#        print "Executing: %s" % (cmdStr)
#        os.system(cmdStr)



#Disable Ethernet, use IB

CMND = 'mpirun -n %s  -hostfile my_nodes_graphene --mca coll_tuned_use_dynamic_rules 1 --mca coll_tuned_bcast_algorithm %s --mca coll_tuned_bcast_algorithm_segmentsize %s  --mca btl_tcp_if_exclude eth0  --mca btl openib,sm,self --mca pml ^cm  --mca plm_rsh_agent "ssh -q -o StrictHostKeyChecking=no" '
ARGS = ' -m %s -M %s -r 1 -i 15 -s 0 -a 4 -d 2'
CMD = CMND + BIN_CMD + ARGS


for p in [128, 64, 32, 16, 8]:
    for alg in [3, 2, 0]:
	for msg in [2**x for x in range(10, 25)]:
	    for seg in [0, 512, 1024, 16384, 32768, 65536]:
                outFileName = " > hbcast_logs_15022015/hbcast_ib_p%s_alg%s_msg%s.log_%s" % (p, alg, msg, timestr)
                cmdStr = CMD % (p, alg, seg, msg, msg) + outFileName
                print "Executing: %s" % (cmdStr)
                os.system(cmdStr)











