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
// Main work code for the BAM to FASTQ conversion.
// ==========================================================================

#ifndef SANDBOX_BAMBI_APPS_BAM2FASTQ_CONVERTER_THREAD_H_
#define SANDBOX_BAMBI_APPS_BAM2FASTQ_CONVERTER_THREAD_H_

#include "job_queue.h"

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class ConversionOptions
// ----------------------------------------------------------------------------

// Data structure that encapsulates only the options required for the conversion.

struct ConversionOptions
{
    // The path to the input file.
    seqan::CharString inputPath;

    // The maximal template length.
    int maxTemplateLength;

    ConversionOptions() : maxTemplateLength(0)
    {}
};

// ----------------------------------------------------------------------------
// Class ConverterThread
// ----------------------------------------------------------------------------

class ConverterThread
{
public:
    // The id of the thread.  -1 for invalid.
    int _threadId;
    // The conversion options.
    ConversionOptions _options;
    // The job queue to use.  NULL if not set.
    JobQueue * _queue;
    // The current job to process.
    ConverterJob _job;
    
    ConverterThread() : _threadId(-1), _queue()
    {}

    ConverterThread(int threadId, ConversionOptions const & options, JobQueue & queue) :
            _threadId(threadId), _options(options), _queue(&queue)
    {}

    void convertOrphans()
    {}

    void convertMapped(ConverterJob const & job)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef SANDBOX_BAMBI_APPS_BAM2FASTQ_CONVERTER_THREAD_H_
