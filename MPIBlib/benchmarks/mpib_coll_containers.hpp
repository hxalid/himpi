#ifndef MPIB_COLL_CONTAINERS_HPP_
#define MPIB_COLL_CONTAINERS_HPP_

#include "mpib_coll_benchmarks.h"
//#include <malloc.h>
#include <stdlib.h>

#include "mpib_defs.h"

/*!
 * \defgroup coll_containers Containers for collective communication operations
 * Data structures describing collective communications to be measured,
 * which are used as an argument of collective benchmark functions.
 * Collective container allocators usually have an argument that points to the implementation of
 * the collective communication operation, for example \ref MPIB_Scatter_container_alloc has an argument
 * \c MPIB_Scatter, pointer to a scatter implementation. This provides three-level extension:
 * 
 * \dotfile benchmarks/coll_containers.dot
 * 
 * The same container can be used in different benchmarks. Container can call different algorithms
 * of the same collective operation.
 * \{
 */

/*!
 * Base container for Scatter(v)/Gather(v), Bcast. Containes a buffer, which can be of
 * the same (Bcast) or different (Scatter(v)/Gather(v)) sizes at processes.
 */
class MPIB_buffer_container: public MPIB_coll_container {
protected:
	char* buffer;

public:
	MPIB_buffer_container() {
		MPIB_coll_container::operation = "";
		MPIB_coll_container::finalize = finalize;
	}

	static int finalize(void* _this, MPI_Comm comm, int root) {
		MPIB_buffer_container* container = (MPIB_buffer_container*) _this;
		free(container->buffer);
		return 0;
	}
};

/*! Base container for Scatter/Gather */
class MPIB_SG_container: public MPIB_buffer_container {
public:
	MPIB_SG_container() {
		MPIB_coll_container::initialize = initialize;
	}

	static int initialize(void* _this, MPI_Comm comm, int root, int M) {
		MPIB_SG_container* container = (MPIB_SG_container*) _this;
		int rank;
		MPI_Comm_rank(comm, &rank);
		int size;
		MPI_Comm_size(comm, &size);
		container->buffer =
				rank == root ?
						(char*) malloc(sizeof(char) * M * size) :
						(char*) malloc(sizeof(char) * M);
		return 0;
	}
};

/*! Scatter container */
class MPIB_Scatter_container: public MPIB_SG_container {
private:
	MPIB_Scatter scatter;

public:
	MPIB_Scatter_container(MPIB_Scatter scatter) {
		MPIB_coll_container::operation = "Scatter";
		MPIB_coll_container::execute = execute;
		this->scatter = scatter;
	}

	static int execute(void* _this, MPI_Comm comm, int root, int M) {
		MPIB_Scatter_container* container = (MPIB_Scatter_container*) _this;
		return container->scatter(container->buffer, M, MPI_CHAR,
				container->buffer, M, MPI_CHAR, root, comm);
	}
};

/*! Gather container */
class MPIB_Gather_container: public MPIB_SG_container {
private:
	MPIB_Gather gather;

public:
	MPIB_Gather_container(MPIB_Gather gather) {
		MPIB_coll_container::operation = "Gather";
		MPIB_coll_container::execute = execute;
		this->gather = gather;
	}

	static int execute(void* _this, MPI_Comm comm, int root, int M) {
		MPIB_Gather_container* container = (MPIB_Gather_container*) _this;
		return container->gather(container->buffer, M, MPI_CHAR,
				container->buffer, M, MPI_CHAR, root, comm);
	}
};

/*! Bcast container */
class MPIB_Bcast_container: public MPIB_buffer_container {
private:
	MPIB_Bcast bcast;

public:
	MPIB_Bcast_container(MPIB_Bcast bcast) {
		MPIB_coll_container::operation = "Bcast";
		MPIB_coll_container::initialize = initialize;
		MPIB_coll_container::execute = execute;
		this->bcast = bcast;
	}

	static int initialize(void* _this, MPI_Comm comm, int root, int M) {
		MPIB_Bcast_container* container = (MPIB_Bcast_container*) _this;
		container->buffer = (char*) malloc(sizeof(char) * M);
		return 0;
	}

	static int execute(void* _this, MPI_Comm comm, int root, int M) {
		MPIB_Bcast_container* container = (MPIB_Bcast_container*) _this;
		return container->bcast(container->buffer, M, MPI_CHAR, root, comm);
	}
};

/*! Reduce container */
class MPIB_Reduce_container: public MPIB_coll_container {
private:
	MPIB_Reduce reduce;
	char* sendbuf;
	char* recvbuf;

public:
	MPIB_Reduce_container(MPIB_Reduce reduce) {
		MPIB_coll_container::operation = "Reduce";
		MPIB_coll_container::initialize = initialize;
		MPIB_coll_container::execute = execute;
		MPIB_coll_container::finalize = finalize;
		this->reduce = reduce;
	}

	static int initialize(void* _this, MPI_Comm comm, int root, int M) {
		MPIB_Reduce_container* container = (MPIB_Reduce_container*) _this;
		container->sendbuf = (char*) malloc(sizeof(char) * M);
		int rank;
		MPI_Comm_rank(comm, &rank);
		container->recvbuf =
				rank == root ? (char*) malloc(sizeof(char) * M) : NULL;
		return 0;
	}

	static int execute(void* _this, MPI_Comm comm, int root, int M) {
		MPIB_Reduce_container* container = (MPIB_Reduce_container*) _this;
		return container->reduce(container->sendbuf, container->recvbuf, M,
				MPI_CHAR, MPI_MAX, root, comm);
	}

	static int finalize(void* _this, MPI_Comm comm, int root) {
		MPIB_Reduce_container* container = (MPIB_Reduce_container*) _this;
		free(container->sendbuf);
		free(container->recvbuf);
		return 0;
	}
};

/*! Container for MPI_Comm_dup-MPI_Comm_free */
class MPIB_Comm_dup_free_container: public MPIB_coll_container {
public:
	MPIB_Comm_dup_free_container() {
		MPIB_coll_container::operation = "MPI_Comm_dup-MPI_Comm_free";
		MPIB_coll_container::execute = execute;
	}

	static int execute(void* _this, MPI_Comm comm, int root, int M) {
		MPI_Comm newcomm;
		MPI_Comm_dup(comm, &newcomm);
		MPI_Comm_free(&newcomm);
		return 0;
	}
};

/*! Base container for Scatterv/Gatherv */
class MPIB_SGv_container: public MPIB_buffer_container {
protected:
	const double* factors;
	int count;
	int* counts;
	int* displs;

public:
	MPIB_SGv_container(const double* _factors) :
			factors(_factors) {
		MPIB_coll_container::initialize = initialize;
		MPIB_coll_container::finalize = finalize;
	}

	static int initialize(void* _this, MPI_Comm comm, int root, int M) {
		MPIB_SGv_container* container = (MPIB_SGv_container*) _this;
		int count = 0;
		int size;
		MPI_Comm_size(comm, &size);
		container->counts = (int*) malloc(sizeof(int) * size);
		container->displs = (int*) malloc(sizeof(int) * size);
		for (int i = 0; i < size; i++) {
			count += container->counts[i] = M * container->factors[i];
			container->displs[i] =
					i == 0 ?
							0 :
							container->displs[i - 1] + container->counts[i - 1];
		}
		int rank;
		MPI_Comm_rank(comm, &rank);
		container->count = container->counts[rank];
		container->buffer = (char*) malloc(
				sizeof(char) * (rank == root ? count : container->count));
		return 0;
	}

	static int finalize(void* _this, MPI_Comm comm, int root) {
		MPIB_buffer_container::finalize(_this, comm, root);
		MPIB_SGv_container* container = (MPIB_SGv_container*) _this;
		free(container->counts);
		free(container->displs);
		return 0;
	}
};

/*! Scatterv container */
class MPIB_Scatterv_container: public MPIB_SGv_container {
private:
	MPIB_Scatterv scatterv;

public:
	MPIB_Scatterv_container(MPIB_Scatterv scatterv, const double* factors) :
			MPIB_SGv_container(factors) {
		MPIB_coll_container::operation = "Scatterv";
		MPIB_coll_container::execute = execute;
		this->scatterv = scatterv;
	}

	static int execute(void* _this, MPI_Comm comm, int root, int M) {
		MPIB_Scatterv_container* container = (MPIB_Scatterv_container*) _this;
		return container->scatterv(container->buffer, container->counts,
				container->displs, MPI_CHAR, container->buffer,
				container->count, MPI_CHAR, root, comm);
	}
};

/*! Gatherv container */
class MPIB_Gatherv_container: public MPIB_SGv_container {
private:
	MPIB_Gatherv gatherv;

public:
	MPIB_Gatherv_container(MPIB_Gatherv gatherv, const double* factors) :
			MPIB_SGv_container(factors) {
		MPIB_coll_container::operation = "Gatherv";
		MPIB_coll_container::execute = execute;
		this->gatherv = gatherv;
	}

	static int execute(void* _this, MPI_Comm comm, int root, int M) {
		MPIB_Gatherv_container* container = (MPIB_Gatherv_container*) _this;
		return container->gatherv(container->buffer, container->count, MPI_CHAR,
				container->buffer, container->counts, container->displs,
				MPI_CHAR, root, comm);
	}
};

/*! Alltoall container */
class MPIB_Alltoall_container: public MPIB_coll_container {
private:
	MPIB_Alltoall alltoall;
	char* sendbuf;
	char* recvbuf;

public:
	MPIB_Alltoall_container(MPIB_Alltoall alltoall) {
		MPIB_coll_container::operation = "Alltoall";
		MPIB_coll_container::initialize = initialize;
		MPIB_coll_container::execute = execute;
		MPIB_coll_container::finalize = finalize;
		this->alltoall = alltoall;
	}

	static int initialize(void* _this, MPI_Comm comm, int root, int M) {
		MPIB_Alltoall_container* container = (MPIB_Alltoall_container*) _this;
		int size;
		MPI_Comm_size(comm, &size);
		container->sendbuf = (char*) malloc(sizeof(char) * M * size);
		container->recvbuf = (char*) malloc(sizeof(char) * M * size);
		return 0;
	}

	static int execute(void* _this, MPI_Comm comm, int root, int M) {
		MPIB_Alltoall_container* container = (MPIB_Alltoall_container*) _this;
		return container->alltoall(container->sendbuf, M, MPI_CHAR,
				container->recvbuf, M, MPI_CHAR, comm);
	}

	static int finalize(void* _this, MPI_Comm comm, int root) {
		MPIB_Alltoall_container* container = (MPIB_Alltoall_container*) _this;
		free(container->sendbuf);
		free(container->recvbuf);
		return 0;
	}
};

/*
 *  Hierarchical MPI containers
 */

/*! HBcast container */
class MPIB_HBcast_container: public MPIB_buffer_container {
private:
	MPIB_HBcast hbcast;
	int num_groups;
	int num_levels;
	int alg_in;
	int alg_out;
public:
	MPIB_HBcast_container(MPIB_HBcast hbcast, int num_groups, int num_levels,
			int alg_in, int alg_out) {
		MPIB_coll_container::operation = "HBcast";
		MPIB_coll_container::initialize = initialize;
		MPIB_coll_container::execute = execute;
		this->hbcast = hbcast;
		this->num_groups = num_groups; //TODO: check +1
		this->num_levels = num_levels;
		this->alg_in = alg_in;
		this->alg_out = alg_out;
	}

	static int initialize(void* _this, MPI_Comm comm, int root, int M) {
		MPIB_HBcast_container* container = (MPIB_HBcast_container*) _this;
		container->buffer = (char*) malloc(sizeof(char) * M);
		return 0;
	}

	static int execute(void* _this, MPI_Comm comm, int root, int M) {
		MPIB_HBcast_container* container = (MPIB_HBcast_container*) _this;
		return container->hbcast(container->buffer, M, MPI_CHAR, root, comm,
				container->num_groups, container->num_levels, container->alg_in,
				container->alg_out); //TODO
	}
};

/*! HReduce container */
class MPIB_HReduce_container: public MPIB_coll_container {
private:
	MPIB_HReduce hreduce;
	char* sendbuf;
	char* recvbuf;
	int num_groups;
	int num_levels;
	int alg_in;
	int alg_out;
public:
	MPIB_HReduce_container(MPIB_HReduce hreduce, int num_groups, int num_levels,
			int alg_in, int alg_out) {
		MPIB_coll_container::operation = "HReduce";
		MPIB_coll_container::initialize = initialize;
		MPIB_coll_container::execute = execute;
		MPIB_coll_container::finalize = finalize;
		this->hreduce = hreduce;
		this->num_groups = num_groups; //TODO
		this->num_levels = num_levels;
		this->alg_in = alg_in;
		this->alg_out = alg_out;
	}

	static int initialize(void* _this, MPI_Comm comm, int root, int M) {
		MPIB_HReduce_container* container = (MPIB_HReduce_container*) _this;
		container->sendbuf = (char*) malloc(sizeof(char) * M);
		int rank;
		MPI_Comm_rank(comm, &rank);
		container->recvbuf =
				rank == root ? (char*) malloc(sizeof(char) * M) : NULL;
		return 0;
	}

	static int execute(void* _this, MPI_Comm comm, int root, int M) {
		MPIB_HReduce_container* container = (MPIB_HReduce_container*) _this;
		return container->hreduce(container->sendbuf, container->recvbuf, M,
				MPI_CHAR, MPI_MAX, root, comm, container->num_groups,
				container->num_levels, container->alg_in, container->alg_out); //TODO
	}

	static int finalize(void* _this, MPI_Comm comm, int root) {
		MPIB_HReduce_container* container = (MPIB_HReduce_container*) _this;
		free(container->sendbuf);
		free(container->recvbuf);
		return 0;
	}
};


class MPIB_HAllreduce_container: public MPIB_coll_container {
private:
	MPIB_HAllreduce hallreduce;
	char* sendbuf;
	char* recvbuf;
	int num_groups;
	int num_levels;
	int alg_in;
	int alg_out;
public:
	MPIB_HAllreduce_container(MPIB_HAllreduce hallreduce, int num_groups, int num_levels,
			int alg_in, int alg_out) {
		MPIB_coll_container::operation = "HAllreduce";
		MPIB_coll_container::initialize = initialize;
		MPIB_coll_container::execute = execute;
		MPIB_coll_container::finalize = finalize;
		this->hallreduce = hallreduce;
		this->num_groups = num_groups; //TODO
		this->num_levels = num_levels;
		this->alg_in = alg_in;
		this->alg_out = alg_out;
	}

	static int initialize(void* _this, MPI_Comm comm, int root, int M) {
		MPIB_HAllreduce_container* container = (MPIB_HAllreduce_container*) _this;
		container->sendbuf = (char*) malloc(sizeof(char) * M);
		container->recvbuf = (char*) malloc(sizeof(char) * M);
		return 0;
	}

	static int execute(void* _this, MPI_Comm comm, int root, int M) {
		MPIB_HAllreduce_container* container = (MPIB_HAllreduce_container*) _this;
		return container->hallreduce(container->sendbuf, container->recvbuf, M,
				MPI_CHAR, MPI_MAX, comm, container->num_groups,
				container->num_levels, container->alg_in, container->alg_out); //TODO
	}

	static int finalize(void* _this, MPI_Comm comm, int root) {
		MPIB_HAllreduce_container* container = (MPIB_HAllreduce_container*) _this;
		free(container->sendbuf);
		free(container->recvbuf);
		return 0;
	}
};



/*! HScatter container */
class MPIB_HScatter_container: public MPIB_SG_container {
private:
	MPIB_HScatter hscatter;
	int num_groups;
	int num_levels;
	int alg_in;
	int alg_out;
public:
	MPIB_HScatter_container(MPIB_HScatter hscatter, int num_groups,
			int num_levels, int alg_in, int alg_out) {
		MPIB_coll_container::operation = "HScatter";
		MPIB_coll_container::initialize = initialize;
		MPIB_coll_container::execute = execute;
		this->hscatter = hscatter;
		this->num_groups = num_groups; // check if I need to add 1
		this->num_levels = num_levels;
		this->alg_in = alg_in;
		this->alg_out = alg_out;
	}

	static int execute(void* _this, MPI_Comm comm, int root, int M) {
		MPIB_HScatter_container* container = (MPIB_HScatter_container*) _this;
		return container->hscatter(container->buffer, M, MPI_CHAR,
				container->buffer, M, MPI_CHAR, root, comm,
				container->num_groups, container->num_levels, container->alg_in,
				container->alg_out);
	}
};



/*! HScatter container */
class MPIB_HGather_container: public MPIB_SG_container {
private:
	MPIB_HGather hgather;
	int num_groups;
	int num_levels;
	int alg_in;
	int alg_out;
public:
	MPIB_HGather_container(MPIB_HGather hgather, int num_groups,
			int num_levels, int alg_in, int alg_out) {
		MPIB_coll_container::operation = "HGather";
		MPIB_coll_container::initialize = initialize;
		MPIB_coll_container::execute = execute;
		this->hgather = hgather;
		this->num_groups = num_groups; //TODO: check if I need +1
		this->num_levels = num_levels;
		this->alg_in = alg_in;
		this->alg_out = alg_out;
	}

	static int execute(void* _this, MPI_Comm comm, int root, int M) {
		MPIB_HGather_container* container = (MPIB_HGather_container*) _this;
		return container->hgather(container->buffer, M, MPI_CHAR,
				container->buffer, M, MPI_CHAR, root, comm,
				container->num_groups, container->num_levels, container->alg_in,
				container->alg_out);
	}
};





/*!
 * \}
 */

#endif /*MPIB_COLL_CONTAINERS_HPP_*/
