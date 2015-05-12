/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "mpi.h"

using namespace std;

#define BUF_SIZE 32768


int main(int argc,char *argv[])
{
	int    n, myid, numprocs, i;
	double startwtime = 0.0, endwtime;
	int    namelen;
	char   processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI::Init(argc,argv);

	int num_procs = MPI::COMM_WORLD.Get_size();
	myid = MPI::COMM_WORLD.Get_rank();
	MPI::Get_processor_name(processor_name,namelen);


	//read input file
	if (argc < 2) {
		MPI::COMM_WORLD.Abort(-1);
	}
	ostringstream oss;
	int length;
	ifstream ifs;
	ofstream ofs;
	if (myid == 0) {

		ifs.open(argv[1], ios_base::binary);
		if (!ifs) {
			cerr << "can't open input file";
			MPI::COMM_WORLD.Abort(-1);
		}
		ifs.seekg(0,ios::end);
		length = ifs.tellg();
		ifs.seekg (0, ios::beg);
	}
	else{
		oss << argv[1] << processor_name;
		ofs.open(oss.str().c_str(),ios_base::binary);
		if (!ofs) {
			cerr << "can't open output file";
			MPI::COMM_WORLD.Abort(-1);
		}
	}

	MPI::COMM_WORLD.Bcast(&length, 1, MPI_INT, 0);

	char* buffer = new char[BUF_SIZE];
	for (int i=0; i<length/BUF_SIZE; i++){
		if (myid == 0)
			ifs.read(buffer,BUF_SIZE);
		MPI::COMM_WORLD.Bcast(buffer, BUF_SIZE, MPI_CHAR, 0);
		if (myid != 0) {
			ofs.write(buffer, BUF_SIZE);
		}
	}
	ifs.read(buffer,BUF_SIZE);
	if (myid == 0)
		ifs.read(buffer,BUF_SIZE);
	MPI::COMM_WORLD.Bcast(buffer, BUF_SIZE, MPI_CHAR, 0);
	if (myid != 0) {
		ofs.write(buffer, BUF_SIZE);
	}
	

	if (myid == 0)
		ifs.close();

	if (myid != 0)
		ofs.close();

	MPI::Finalize();
	return 0;
}
