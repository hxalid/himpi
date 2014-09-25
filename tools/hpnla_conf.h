#ifndef HPNLA_CONF_H_
#define HPNLA_CONF_H_

#include <stdio.h>
#include <mpi.h>


#ifdef __cplusplus
extern "C" {
#endif
  /*!
    \defgroup hpnla_conf_group Configuration of MPI process
    Provides configuration settings specific to each MPI process. 
    - Hardware to run on.
    - File names to read and write to.
    - Distribution to use.

    Allows the distribution outputted by the partitioner to be associated 
    with a specific device on a host. Devices are addressed by their hostname and rank_intra (MPI rank within a node). 
    The device_type parameter can be used to cause the execution of the kernel or routine to be done on a specific device, 
    eg. cpu or gpu, however since it is device_type char* more specific device types such as gpu1, gpu2, fpga etc. may be used.
   
    This data is stored in hpnla_device_info structure and input and output
    to conf_file. It is also used in file names of FPMs and in the
    partition.dist file.
   
    Example conf_file
   
    \code
#hostname                      rank_intra  bind  device_type  suboptions
adonis-1.grenoble.grid5000.fr  0           na     gpu           <some optional params for gpu on adonis-1: device_id_start number_devices memory_alloc_mode algorithm is_overlap>
adonis-1.grenoble.grid5000.fr  1           na     cpu           <some optional params for first cpu on adonis-1>
adonis-1.grenoble.grid5000.fr  2           na     cpu           
edel-1.grenoble.grid5000.fr    0           na     cpu           OMP_NUM_THREADS=8 <If implemented, this MPI process could use all cores on host>
    \endcode
   
    Example FPM filename:
    \code
adonis-1.grenoble.grid5000.fr.0.gpu.fpm
    \endcode
   
    Example partition.dist file:
    \code
#size 4
#D 1400
#host                          rank_intra  d    time       i  x       complexity  speed
adonis-1.grenoble.grid5000.fr  0           800  1.000e-02  0  0       800         80000
adonis-1.grenoble.grid5000.fr  1           100  1.000e-02  1  800     100         10000
adonis-1.grenoble.grid5000.fr  2           100  1.000e-02  2  900     100         10000
edel-1.grenoble.grid5000.fr    0           400  1.000e-02  3  1000    400         40000
    \endcode
   
    \{
   */

  /*! Maximum lenght of the device string */
  #define HPNLA_DEVICE_MAX_STRING 12

  /*! Maximum lenght of the string for process/thread binding */
  #define HPNLA_BIND_MAX_STRING 256

  /*! Maximum lenght of the string for suboptions */
  #define HPNLA_SUBOPTION_MAX_STRING 256

  /*! Default device */
  #define HPNLA_DEVICE_DEFAULT cpu

  /*! Provides info about the device the MPI process is running on. */
  typedef struct hpnla_process_conf {
    /*! Hostname */
    char* hostname;
    /*! Rank within node */
    int rank_intra;
    /*! string giving CPU core or device number to bind to */
    char* bind;
    /*! Device type, string (eg. cpu, gpu) */
    char* device_type;
    /*! Suboptions */
    char* subopts;
    // todo: remove the rest
    /*! Number of process-specific parameters (optional, parsed by user coded routine) */
    int argc;
    /*! Array of process-specific parameters (optional, parsed by user coded routine) */
    char** argv;
  } hpnla_process_conf;

  
  /*! Frees memory allocated in hpnla_process_conf */
  void hpnla_conf_free(hpnla_process_conf conf);

  /*! Reads a file and returns the configuration struct for that MPI process.
   * \param comm MPI communicator (should be MPI_COMM_WORLD).
   * \param filename filename string (default: ./conf_file).
   * */
  hpnla_process_conf hpnla_get_conf(MPI_Comm comm, char* filename);

  /*! Reads a machine file and returns an array of configuration structs. This function does not use MPI communicator
   * and hence is suitable for serial programmes such as \ref partitioner. conf's are ordered as per machine file.
   * \param filename (default: conf_file).
   * \param size on return will be the size of that array and equal to the size of comm
   * */
  hpnla_process_conf* hpnla_get_conf_all(char* filename, int* size);

  /*! Reads a file and returns an array of configuration structs in the same order as processes MPI rank. 
   * \param comm MPI communicator (should be MPI_COMM_WORLD).
   * \param filename (default: conf_file).
   * \param size on return will be the size of that array and equal to the size of comm
   * */
  hpnla_process_conf* hpnla_get_conf_all_sorted(MPI_Comm comm, char* filename, int* size);

  /*! Prints a default usable conf_file for all of the processes in the current mpirun. 
   * This file can then be customised further by the user.
   * \param comm MPI communicator (should be MPI_COMM_WORLD).
   * \param root root processor to print
   * \param file File pointer to output file (default: conf_file).
   * \param default_device_type prints all device_type's to be this (eg "cpu").
   * \param subopts Suboptions passed to kernal
   * */
void hpnla_print_conf(MPI_Comm comm, int root, FILE* file, char* default_device_type, char* subopts);

  /*!
   * Creates multiple communicators, each consisting of the processes running on the same processor.
   * \note Do not forget to free the result communicators.
   * \param comm MPI communicator
   * \param comm_intra intra communicator (non-MPI_COMM_NULL at all processes)
   * \return error status
   */
  int hpnla_comm_intra(MPI_Comm comm, MPI_Comm* comm_intra);

  /*!
   * Creates a single communicator consisting of the processes running on different processors
   * \note Do not forget to free the result communicator.
   * \param comm MPI communicator
   * \param comm_inter inter communicator (non-MPI_COMM_NULL only on the first process running on the processor)
   * \return error status
   */
  int hpnla_comm_inter(MPI_Comm comm, MPI_Comm* comm_inter);

  /*!
   * Gathers hostnames from all processors in comm to root process. Returns to root process a pointer to array of size: MPI_MAX_PROCESSOR_NAME * size.
   * \param root root processor to hold result
   * \param comm MPI communicator
   */
  char** hpnla_gather_hostnames(int root, MPI_Comm comm);

  /*!
   * Same as hpnla_gather_hostnames but all processors will have result.
   * \param comm MPI communicator
   */
  char** hpnla_allgather_hostnames(MPI_Comm comm);

  /*!
   * \}
   */
#ifdef __cplusplus
}
#endif

#endif /*HPNLA_CONF_H_*/

