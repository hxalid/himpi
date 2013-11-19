//todo: remove unnecessary
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_heapsort.h>
#include "hpnla_conf.h"
#include "hpnla_debug.h"

void hpnla_conf_free(hpnla_process_conf conf) {
  free(conf.hostname);
  free(conf.bind);
  free(conf.device_type);
  free(conf.subopts);
//todo: remove the rest
  int i;
  for(i = 0; i < conf.argc; i++){
    free(conf.argv[i]);
  }
  free(conf.argv);
}

//local function
int match(const char *s, char p){
  int c = 1;
  while (*s != '\0') {
    if (strncmp(s++, &p, 1)) continue;
    c++;
  }
  return c;
}

//local function
int get_list_param(const char* input, int* argc, char*** argv){
  *argv = NULL;
  if(input==NULL){
    *argc = 0;
    return -1;
  }
  if(strlen(input) < 2) {
    *argc = 0;
    return -1;
  }
  char *line = strdup(input);
  int size = match(line, ' ');
  *argv = (char**)malloc(sizeof(char*) * size);
  char *pch = strtok (line," \t\n");
  int i = 0;
  while (pch != NULL)
  {
    if(strcmp(pch,"") != 0 ) {
      (*argv)[i] = strdup(pch);
      i++;
    }
    pch = strtok (NULL, " \t\n");
  }
  *argc = i;
  free(line); line = NULL;
  return 0;
}

//local function
//todo: is this really needed?
int* get_index(MPI_Comm comm, MPI_Comm comm_intra){
  int rank_intra = 0;
  int size;
  int* index;

  MPI_Comm_rank(comm_intra, &rank_intra);
  MPI_Comm_size(comm, &size);
  index = (int*)malloc(sizeof(int) * size);
  MPI_Allgather(&rank_intra, 1, MPI_INT, index, 1, MPI_INT, comm);
  return index;
}

hpnla_process_conf* hpnla_get_conf_all(char* filename, int* size){
  if (filename == NULL){
    fprintf(stderr, "Error filename null %s\n", __func__);
    MPI_Abort(MPI_COMM_WORLD, 200);
  }
  FILE* stream;
  stream = fopen(filename, "r");
  if (stream == NULL) {
    debug_print(0, (size_t) 0,(size_t) 0,
        "Try to open the configuration file %s\n", filename);
    perror("fopen");
    MPI_Abort(MPI_COMM_WORLD, 201);
  }
  char* line = NULL;
  size_t number = 0;
  int err;
  int n = 0;
  hpnla_process_conf* confs = NULL;
  while ((err = getline(&line, &number, stream)) != -1) {
    // skip comment headers
    if (line[0] == '#'|| line[0] == '\n')
      continue;
    confs = realloc(confs, sizeof(hpnla_process_conf) * (n + 1));
    confs[n].hostname = malloc(sizeof(char) * MPI_MAX_PROCESSOR_NAME);
    confs[n].bind = malloc(sizeof(char) * FUPERMOD_BIND_MAX_STRING);
    confs[n].device_type = malloc(sizeof(char) * FUPERMOD_DEVICE_MAX_STRING);
    confs[n].subopts = malloc(sizeof(char) * FUPERMOD_SUBOPTION_MAX_STRING);
    int pos = 0;
    int err = sscanf(line, "%s %d %s %s %s %n", confs[n].hostname, &confs[n].rank_intra, confs[n].bind, confs[n].device_type, confs[n].subopts, &pos);
    if (err < 4) {
      fprintf(stderr, "Error reading line%d: \"%s\" err:%d in %s\n", n, line, err, __func__);
      MPI_Abort(MPI_COMM_WORLD, 204);
    }
    if (err < 5) { //No subopts.
      confs[n].subopts[0] = '\0';
      confs[n].argc = 0;
      confs[n].argv = NULL;
    } else {
      get_list_param(&line[pos], &(confs[n].argc), &(confs[n].argv));
    }
    n++;
  }
  free(line);
  *size = n;
  return confs;
}

hpnla_process_conf* hpnla_get_conf_all_sorted(MPI_Comm comm, char* filename, int* size){
  //get all confs from file.
  hpnla_process_conf* confs = hpnla_get_conf_all(filename, size);

  int mpi_size;
  MPI_Comm_size(comm, &mpi_size);
  if (*size != mpi_size) {
    fprintf(stderr, "Error:conf file size != mpi_size %d != %d\n", *size, mpi_size);
    MPI_Abort(comm, 205);
  }
  
  int rank_intra;
  MPI_Comm comm_node;
  hpnla_comm_intra(comm, &comm_node);
  MPI_Comm_rank(comm_node, &rank_intra);
  MPI_Comm_free(&comm_node);

  char** names = hpnla_allgather_hostnames(comm);
  int* rank_intras = (int*)malloc(sizeof(int) * *size);
  MPI_Allgather(&rank_intra, 1, MPI_INT, rank_intras, 1, MPI_INT, comm);
  hpnla_process_conf* confs_sorted = (hpnla_process_conf*)malloc(sizeof(hpnla_process_conf) * *size);
  
  int i;
  for(i = 0; i < *size; i++){
    int j;
    for(j = 0; j < *size; j++){
      if (confs[j].rank_intra == rank_intras[i] && (strcmp(confs[j].hostname, names[i]) == 0)) {
        confs_sorted[i] = confs[j];
        break;
      }
    }
    if (j == *size) {
      fprintf(stderr, "Error: No configuration found for %s %d\n", names[i], rank_intras[i]);
      MPI_Abort(comm, 202);
    }
  }
  return confs_sorted;
}

hpnla_process_conf hpnla_get_conf(MPI_Comm comm, char* filename){
  //get all confs from file.
  int size;
  hpnla_process_conf* confs = hpnla_get_conf_all(filename, &size);

  //get hostname
  char name[MPI_MAX_PROCESSOR_NAME];
  int len;
  MPI_Get_processor_name(name, &len);

  int rank_intra;
  MPI_Comm comm_node;
  hpnla_comm_intra(comm, &comm_node);
  MPI_Comm_rank(comm_node, &rank_intra);
  MPI_Comm_free(&comm_node);
  hpnla_process_conf* conf = NULL;
  int i;
  for(i = 0; i < size; i++){
    if (confs[i].rank_intra == rank_intra && (strcmp(confs[i].hostname, name) == 0)) {
      conf = &(confs[i]);
      break;
    } else {
      hpnla_conf_free(confs[i]);
    }
  }
  if (conf == NULL) {
    debug_print(0, (size_t) 0, (size_t) 0,
        "No configuration found for %s %d\n", name, rank_intra );
    MPI_Abort(comm, 202);
  }
  //free(confs); confs = NULL;
  return *conf;
}

void hpnla_print_conf(MPI_Comm comm, int root, FILE* file, char* default_device_type, char* subopts){
  //TODO
  if(file == NULL) file = stdout;
  if(default_device_type == NULL) default_device_type = "cpu";

  MPI_Comm comm_node;
  hpnla_comm_intra(comm, &comm_node);

  char** names = hpnla_allgather_hostnames(comm);

  int* index = get_index( comm,  comm_node);
  MPI_Comm_free(&comm_node);
  int size;
  MPI_Comm_size(comm, &size);
  int i=0;
  int rank;
  MPI_Comm_rank(comm, &rank);
  if(rank == root){
    fprintf(file, "#hostname\trank_intra\tbind\tdevice_type\tsuboptions\n");
    for(i = 0; i < size; i++){
      fprintf(file, "%s\t%d\t%s\t%s\t%s\n", names[i], index[i], "na", default_device_type, subopts);
    }
  }
  free(names[0]);
  free(names);
  free(index);
  return;
}


int hpnla_comm_intra(MPI_Comm comm, MPI_Comm* comm_intra) {
  char name[MPI_MAX_PROCESSOR_NAME];
  int len;
  int size;
  char** names = hpnla_allgather_hostnames(comm);
  int color = -1;
  int i = 0;

  MPI_Get_processor_name(name, &len);
  MPI_Comm_size(comm, &size);
  while (i < size){
    if (strcmp(name, names[i]) == 0) {
      break;
    }
    i++;
  }
  color = i;
  free(names);
  // split the communicator
  int rank;
  MPI_Comm_rank(comm, &rank);
  return MPI_Comm_split(comm, color, rank, comm_intra);
}

int hpnla_comm_inter(MPI_Comm comm, MPI_Comm* comm_inter) {
  // collect processor names
  char name[MPI_MAX_PROCESSOR_NAME];
  int len;
  MPI_Get_processor_name(name, &len);
  int size;
  MPI_Comm_size(comm, &size);
  char* names = (char*)malloc(sizeof(char) * MPI_MAX_PROCESSOR_NAME * size);
  MPI_Allgather(name, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, names, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, comm);
  // sort processor names
  size_t* indices = (size_t*)malloc(sizeof(size_t) * size);
  gsl_heapsort_index(indices, names, size, sizeof(char) * MPI_MAX_PROCESSOR_NAME, (int(*)(const void*, const void*))strcmp);
  // colouring (0, key) for the first process on the processor, MPI_UNDEFINED elsewhere
  int color = MPI_UNDEFINED;
  int key = 0;
  int rank;
  MPI_Comm_rank(comm, &rank);
  int i = 0;
  // colour until the current process is found
  do
    if (i == 0 || strcmp(&names[MPI_MAX_PROCESSOR_NAME * indices[i - 1]], &names[MPI_MAX_PROCESSOR_NAME * indices[i]])) {
      if (indices[i] == (size_t) rank) {
        // the first process on the processor
        color = 0;
      } else {
        // the next processor
        key++;
      }
    }
  while (indices[i++] != (size_t) rank);
  free(names);
  free(indices);
  // split the communicator
  return MPI_Comm_split(comm, color, key, comm_inter);
}

char** hpnla_gather_hostnames(int root, MPI_Comm comm) {
  char** friendly_names = NULL;
  char* hostnames = NULL;
  char hostname[MPI_MAX_PROCESSOR_NAME];
  int len;
  MPI_Get_processor_name(hostname, &len);
  int rank;
  MPI_Comm_rank(comm, &rank);
  int size;
  if (rank == root){
    MPI_Comm_size(comm, &size);
    hostnames = (char*)malloc(sizeof(char) * MPI_MAX_PROCESSOR_NAME * size);
  }
  MPI_Gather(hostname, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, hostnames, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, root, comm);
  if (rank == root){
    int i;
    for(i = 0; i < size;i++) friendly_names[i] = &hostnames[MPI_MAX_PROCESSOR_NAME * i];
  }
  return friendly_names;
}

char** hpnla_allgather_hostnames(MPI_Comm comm) {
  char name[MPI_MAX_PROCESSOR_NAME];
  int len;
  int size;
  int i;
  char** friendly_names;
  char*  names;

  MPI_Get_processor_name(name, &len);
  MPI_Comm_size(comm, &size);
  friendly_names = malloc(sizeof(char*) * size);
  names          = malloc(sizeof(char) * MPI_MAX_PROCESSOR_NAME * size);

  MPI_Allgather(name, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, names, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, comm);

  for( i = 0; i < size;i++) friendly_names[i] = &names[MPI_MAX_PROCESSOR_NAME * i];
  return friendly_names;
}

