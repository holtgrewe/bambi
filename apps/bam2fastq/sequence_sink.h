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
// Thread-safe writing to sequence files.
// ==========================================================================

#ifndef SANDBOX_BAMBI_APPS_BAM2FASTQ_SEQUENCE_SINK_H_
#define SANDBOX_BAMBI_APPS_BAM2FASTQ_SEQUENCE_SINK_H_

#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include <omp.h>


// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

typedef seqan::StringSet<seqan::CharString, seqan::Owner<seqan::ConcatDirect<> > > TCharStringSet;
typedef seqan::StringSet<seqan::Dna5String, seqan::Owner<seqan::ConcatDirect<> > > TDna5StringSet;

class SequenceSink
{
public:
    // The names of the left/right/single output files.
    seqan::CharString pathLeft, pathRight, pathSingle;

    // Sequence stream objects for writing out the sequences.
    seqan::SequenceStream outLeft, outRight, outSingle;

    // Whether or not we have already opened the sequence streams.
    bool isOpenLeft, isOpenRight, isOpenSingle;

    // Whether or not to write out paired-end reads in an interleaved fashion.
    bool interleaved;
    // Whether or not to compress the output files.
    bool gzip;

    // The OpenMP lock for thread safety.
    omp_lock_t _lock;

    SequenceSink(bool interleaved, bool gzip, char const * pathLeft, char const * pathRight, char const * pathSingle) :
            interleaved(false), gzip(gzip), pathLeft(pathLeft), pathRight(pathRight), pathSingle(pathSingle),
            isOpenLeft(false), isOpenRight(false), isOpenSingle(false)
    {
        omp_init_lock(&_lock);
    }

    ~SequenceSink()
    {
        omp_destroy_lock(&_lock);
    }

    int write(TCharStringSet const & leftIds, TDna5StringSet const & leftSeqs, TCharStringSet const & leftQuals,
              TCharStringSet const & rightIds, TDna5StringSet const & rightSeqs, TCharStringSet const & rightQuals,
              TCharStringSet const & singleIds, TDna5StringSet const & singleSeqs, TCharStringSet const & singleQuals)
    {
        SEQAN_ASSERT_EQ(length(leftIds), length(rightIds));

        // Lock data structure.
        RaiiLock lockGuard(_lock);

        // Write out to single files.
        if (!empty(singleIds))
        {
            if (!isOpenSingle && _open(outSingle, isOpenSingle, pathSingle) != 0)
                return 1;
            if (writeAll(outSingle, singleIds, singleSeqs, singleQuals) != 0)
            {
                SEQAN_OMP_PRAGMA(critical (io))
                {
                    std::cerr << "ERROR: Could not write to single output file!\n";
                }
                return 1;
            }
        }

        // Nothing to be done for paired reads.  Early exit.
        if (empty(leftIds))
            return 0;
        
        // Write out mates.
        if (!interleaved)
        {
            if (!isOpenLeft && _open(outLeft, isOpenLeft, pathLeft) != 0)
                return 1;
            if (writeAll(outLeft, leftIds, leftSeqs, leftQuals) != 0)
            {
                SEQAN_OMP_PRAGMA(critical (io))
                {
                    std::cerr << "ERROR: Could not write to left output file!\n";
                }
                return 1;
            }
            if (!isOpenRight && _open(outRight, isOpenRight, pathRight) != 0)
                return 1;
            if (writeAll(outRight, rightIds, rightSeqs, rightQuals) != 0)
            {
                SEQAN_OMP_PRAGMA(critical (io))
                {
                    std::cerr << "ERROR: Could not write to right output file!\n";
                }
                return 1;
            }
        }
        else
        {
            if (!isOpenLeft && _open(outLeft, isOpenLeft, pathLeft) != 0)
                return 1;
            for (unsigned i = 0; i < length(leftIds); ++i)
            {
                if (writeRecord(outLeft, leftIds[i], leftSeqs[i], leftQuals[i]) != 0)
                {
                    SEQAN_OMP_PRAGMA(critical (io))
                    {
                        std::cerr << "ERROR: Could not write to left output file!\n";
                    }
                    return 1;
                }
                if (writeRecord(outLeft, rightIds[i], rightSeqs[i], rightQuals[i]) != 0)
                {
                    SEQAN_OMP_PRAGMA(critical (io))
                    {
                        std::cerr << "ERROR: Could not write to right output file!\n";
                    }
                    return 1;
                }
            }
        }

        return 0;
    }

    int _open(seqan::SequenceStream & outStream, bool & flag, seqan::CharString const & path)
    {
        seqan::SequenceStream::FileType fileType = gzip ? seqan::SequenceStream::GZ : seqan::SequenceStream::PLAIN_TEXT;
        open(outStream, toCString(path), seqan::SequenceStream::WRITE, seqan::SequenceStream::FASTQ, fileType);
        if (!isGood(outStream))
        {
            SEQAN_OMP_PRAGMA(critical (io))
            {
                std::cerr << "ERROR: Could not open " << path << " for writing!\n";
            }
            return 1;
        }
        flag = true;
        return 0;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef SANDBOX_BAMBI_APPS_BAM2FASTQ_SEQUENCE_SINK_H_
