/*
 * SliceTraversals.h
 *
 *  Created on: Feb 27, 2020
 *      Author: tchipevn
 */

#ifndef SLICE1D_H_
#define SLICE1D_H_

#include <vector>
#include <omp.h>

class Slice1D {
public:
	Slice1D() {
		int numThreads;

		#pragma omp parallel
		{
			#pragma omp master
			numThreads = omp_get_num_threads();
		}
		_locks.resize(numThreads+1, nullptr);

		#pragma omp parallel
		{
			int myId = omp_get_thread_num();
			_locks[myId+1] = new omp_lock_t();
			omp_init_lock(_locks[myId+1]);
		}
		_locks[0] = new omp_lock_t();
		omp_init_lock(_locks[0]);
	}

	~Slice1D() {
		#pragma omp parallel
		{
			int myId = omp_get_thread_num();
			omp_destroy_lock(_locks[myId+1]);
			delete _locks[myId+1];
		}
		omp_destroy_lock(_locks[0]);
		delete _locks[0];
	}

	enum LockType {
		MY_LOCK = 0,
		NEXT_LOCK = 1
	};

	void acquireLock(int myId, LockType l) {
		omp_set_lock(_locks[myId + l]);
	}

	void releaseLock(int myId, LockType l) {
		omp_unset_lock(_locks[myId + l]);
	}

	/** note: this spans the space [0, numCells) */
	static int getStart(int numCells, int numThreads, int threadId) {
		return numCells * threadId / numThreads;
	}

	/** note: this spans the space [0, numCells) */
	static int getEnd(int numCells, int numThreads, int threadId) {
		return numCells * (threadId + 1) / numThreads;
	}

	/** the maximal number of threads ensuring two slices per thread */
	static int getMaxThreads(int numGridpoints) {
		return numGridpoints / 2;
	}

	/** the maximal usable number of threads */
	static int getActualThreads(int numGridpoints, int ompNumThreads) {
		return std::min(ompNumThreads, getMaxThreads(numGridpoints));
	}


private:
	std::vector<omp_lock_t *> _locks;

};

#endif /* SLICE1D_H_ */
