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
	Slice1D(int numThreads, int numGridpoints) :
		_numThreads(numThreads), _numGridpoints(numGridpoints) {
		_locks.resize(_numThreads+1, nullptr);

		for (int i = 0; i < _numThreads + 1; ++i) {
			_locks[i] = new omp_lock_t();
			omp_init_lock(_locks[i]);
		}
	}

	~Slice1D() {
		for (int i = 0; i < _numThreads + 1; ++i) {
			omp_destroy_lock(_locks[i]);
			delete _locks[i];
		}
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
	int getStart(int threadId) {
		return _numGridpoints * threadId / _numThreads;
	}

	/** note: this spans the space [0, numCells) */
	int getEnd(int threadId) const {
		return _numGridpoints * (threadId + 1) / _numThreads;
	}

	/** the maximal number of threads ensuring two slices per thread */
	int getMaxThreads() const {
		return _numGridpoints / 2;
	}

	/** the maximal usable number of threads */
	int getActualThreads() const {
		return std::min(_numThreads, getMaxThreads());
	}


private:
	int _numThreads;
	int _numGridpoints;
	std::vector<omp_lock_t *> _locks;

};

#endif /* SLICE1D_H_ */
