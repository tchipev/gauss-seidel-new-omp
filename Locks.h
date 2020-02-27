/*
 * SliceTraversals.h
 *
 *  Created on: Feb 27, 2020
 *      Author: tchipevn
 */

#ifndef LOCKS_H_
#define LOCKS_H_

#include <vector>
#include <omp.h>

class Locks {
public:
	Locks() {
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

	~Locks() {
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

private:

	std::vector<omp_lock_t *> _locks;

};

#endif /* LOCKS_H_ */
