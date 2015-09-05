import sys
import os

CONFIG_PATH = './silgetsin/'
ALG_COUNT = 7
ALG_ID = '7'

class HMPI_Runner(object):
    
    def __init__(self, ):
        self.alg = 0
        
    @staticmethod    
    def build_command(num_procs, host_file, exe_name, dynamic_rules_filename, exe_args):
        output_file = dynamic_rules_filename[len(CONFIG_PATH):]
        mca_params_dynamics = '--mca coll_tuned_use_dynamic_rules 1 --mca coll_tuned_dynamic_rules_filename %s' % (dynamic_rules_filename)
        #mca_params_network = '--mca btl_tcp_if_exclude eth0  --mca btl openib,sm,self --mca pml ^cm  --mca plm_rsh_agent "ssh -q -o StrictHostKeyChecking=no" '
        return "mpirun -n %s -hostfile %s %s %s %s &> log_%s " % (num_procs, host_file, mca_params_dynamics, exe_name, exe_args, output_file)
        #return "%s --prefix %s -n %s -hostfile %s %s %s %s %s &> log_%s " % (MPI_EXE, MPI_PREFIX, num_procs, host_file, mca_params_dynamics, mca_params_network, exe_name, exe_args, output_file)
    
    @staticmethod    
    def execute(num_procs, host_file, exe_name, dynamic_rules_filename, exe_args):
        command = HMPI_Runner.build_command(num_procs, host_file, exe_name, dynamic_rules_filename, exe_args)
        print(command)
        os.system(command)

    @staticmethod
    def create_config_file(file_name, coll_ids, comm_sizes, comm_msg_sizes):
        with open(file_name, 'w') as fout:
            num_colls = len(coll_ids)
            fout.write(str(num_colls) + '\n')
            for nc in coll_ids:
                fout.write(nc + '\n')
                num_comm_size = len(comm_sizes)
                fout.write(str(num_comm_size) + '\n')
                for cs in comm_sizes:
                    fout.write(str(cs) + '\n')
                    num_msg_sizes = len(comm_msg_sizes.get(str(cs)))
                    fout.write(str(num_msg_sizes) + '\n')
                    msg_list = comm_msg_sizes.get(str(cs))
                    for m in msg_list:
                        fout.write(m + '\n')
            
            



if __name__ == '__main__':
    p = 128
    groups = [2, 4, 8]
    pgs = [16, 32, 64, p]
    comm_sizes = groups+pgs
        
    alg_count = ALG_COUNT
    alg_id = ALG_ID
    coll_ids = [alg_id] 
    msg_sizes= [0]
    msg_sizes.extend([2**m for m in range(0, 17, 2)])
            
    for alg_in in range(0, alg_count):
        for alg_out in range(0, alg_count):
            file_name = '%sompi_conf_p%s_alg_in%s_alg_out%s' % (CONFIG_PATH, p, alg_in, alg_out)
            
            comm_msgs = {}
            msg_sizes_g = []
            msg_sizes_pg = []   
                                    
            for m in msg_sizes:
                msg_sizes_pg.append('%s %s %s %s' % (m, alg_in, 0, 0))
                msg_sizes_g.append('%s %s %s %s' % (m, alg_out, 0, 0))
            
            for g in groups:
                comm_msgs[str(g)] = msg_sizes_g
            for pg in pgs:
                comm_msgs[str(pg)] = msg_sizes_pg
            
            
            HMPI_Runner.create_config_file(file_name, coll_ids, comm_sizes, comm_msgs)
         #   HMPI_Runner.execute(16, 'hostlarim', 'test_hbcast', file_name, '-i 1 -m 1024 -M 1024')



# [message_size alg topo segmentation] 
#    msg_sizes_p0 = ['0 1 0 0', '1024 0 0 0', '8192 6 0 0', '16384 5 0 0', '32768 6 0 0', '262144 3 0 0', '524288 4 0 0']
#    msg_sizes_p16 = ['0 1 0 0', '1024 6 0 0', '8192 6 0 0', '16384 5 0 0', '32768 2 0 0', '262144 3 0 0', '524288 4 0 0']
#    HMPI_Runner.create_config_file('fayildi-bu', [alg_id], ['0', '16'], {'0':msg_sizes_p0, '16':msg_sizes_p16})
  
#MPI_EXE = '/home/khasanov/software/my_openmpi-1.8.4/bin/mpirun'
#MPI_PREFIX = '/home/khasanov/software/my_openmpi-1.8.4/'                  
                    
                    
                    
                    