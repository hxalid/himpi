#ifndef DOXYMAIN_H_
#define DOXYMAIN_H_

/*!
 * \mainpage Introduction
 *  The Hierarchical MPI (HiMPI) addresses the problem of reduction of the communication cost of traditional message-passing data-parallel applications on large-scale distributed-memory computing systems. The approach it employs is a traditional methodology widely used for dealing with the complexity of coordination and management of a large number of actors, namely, the hierarchical approach. According to this technique, thousands or millions of actors are structured, and instead of interacting with a large number of peers, they coordinate their activities with one superior and a small number of peers and inferiors. This way the overhead of interaction is significantly reduced.

The hierarchy is achieved by organizing the MPI processes into logical groups. The HiMPI library is layered on top of the existing MPI implementations, provides automatic estimation and selections of the optimal parameters in the hierarchical algorithms and is portable to any parallel platforms that supports MPI. Currently, HiMPI supports MPI broadcast, reduce, allreduce, scatter and gather operations.
 *
 *
 *
 *
 *
 *
 * \section authors Authors
 * Khalid Hasanov, Alexey Lastovetsky
 * 
 * Heterogeneous Computing Laboratory\n
 * School of Computer Science and Informatics, University College Dublin\n
 * Belfield, Dublin 4, Ireland\n
 * http://hcl.ucd.ie
 * 
 * \latexonly
 * \bibliographystyle{plain}
 * \bibliography{mpib}
 * \endlatexonly
 * 
 * \page installation Installation
 * \include README
 * 
 * \page design The software design
 *  Under construction
 * 
 */

#endif /*DOXYMAIN_H_*/
