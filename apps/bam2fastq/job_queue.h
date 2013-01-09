// ==========================================================================
//                                BAM 2 FASTQ
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Simple, thread-safe data structure for jobs (intervals on chromsomes).
// ==========================================================================

#ifndef SANDBOX_BAMBI_APPS_BAM2FASTQ_JOB_QUEUE_H_
#define SANDBOX_BAMBI_APPS_BAM2FASTQ_JOB_QUEUE_H_

#include <queue>

#include <omp.h>

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class ConverterJob
// ----------------------------------------------------------------------------

// Data structure that encapsulates one conversion job.  A conversion job consists of an interval on a chromsome.  Any
// overlapping window computation is handled by the conversion threads themselves.
//
// A rId value of -1 leads marks an invalid job.

struct ConverterJob
{
    // Id/index of the reference id.
    int rId;
    // Begin position of the job on the reference/chromosome.
    int beginPos;
    // End position of the job on the reference/chromsome.
    int endPos;

    ConverterJob() : rId(-1), beginPos(0), endPos(0)
    {}

    ConverterJob(int rId, int beginPos, int endPos) : rId(rId), beginPos(beginPos), endPos(endPos)
    {}
};

// ----------------------------------------------------------------------------
// Class JobQueue
// ----------------------------------------------------------------------------

// A data structure for the thread-safe accessing of threads.

class JobQueue
{
public:
    // The actual job queue.
    std::queue<ConverterJob> _queue;

    // The OpenMP lock for thread safety.
    mutable omp_lock_t _lock;

    // Constructor, initialize lock.
    JobQueue()
    {
        omp_init_lock(&_lock);
    }

    // Destructor, free lock again.
    ~JobQueue()
    {
        omp_destroy_lock(&_lock);
    }

    // Push one lock.
    void push(ConverterJob const & job)
    {
        // Lock data structure.
        omp_set_lock(&_lock);

        _queue.push(job);

        // Unlock data structure.
        omp_unset_lock(&_lock);
    }

    // Push a range of locks.
    template <typename TIter>
    void push(TIter itBegin, TIter itEnd)
    {
        // Lock data structure.
        omp_set_lock(&_lock);

        for (TIter it = itBegin; it != itEnd; ++it)
            _queue.push(*it);

        // Unlock data structure.
        omp_unset_lock(&_lock);
    }

    // Returns whether the queue is empty.
    bool empty() const
    {
        bool res = _queue.empty();
        return res;
    }

    // Returns the number of elements in the queue.
    size_t size() const
    {
        return _queue.size();
    }

    // Get the next job.  Return true if there was a job left and false otherwise.
    bool pop(ConverterJob & job)
    {
        // Initialize return value with false (could not get job), and initially mark job as invalid.
        bool res = false;
        job.rId = -1;

        // Lock data structure.
        omp_set_lock(&_lock);

        if (!_queue.empty())
        {
            job = _queue.front();
            _queue.pop();
            res = true;
        }

        // Unlock data structure.
        omp_unset_lock(&_lock);

        return res;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef SANDBOX_BAMBI_APPS_BAM2FASTQ_JOB_QUEUE_H_
