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

#include <memory>

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
// Class ConversionStats
// ----------------------------------------------------------------------------

// Statistics on conversion.
struct ConversionStats
{
    int numOrphans;
    int numSingletonOrphans;
    int numPairedOrphans;

    ConversionStats() : numOrphans(0), numSingletonOrphans(0), numPairedOrphans(0)
    {}

    void countSingletonOrphan(int num = 1)
    {
        numSingletonOrphans += num;
    }

    void countPairedOrphan(int num = 1)
    {
        numPairedOrphans += num;
    }
};

// ----------------------------------------------------------------------------
// Class ConverterThread
// ----------------------------------------------------------------------------

class ConverterThread
{
public:
    // The state of the converter thread.
    enum {
        START,       // At beginning, nothing done yet.
        FIRST_STEP,  // In the first step, writing to temporary sequence files.
        SECOND_STEP  // In the second step, filling output FASTQ files.
    } _state;

    // The id of the thread.  -1 for invalid.
    int _threadId;
    // The conversion options.
    ConversionOptions _options;
    // The job queue to use.  NULL if not set.
    JobQueue * _queue;
    // The current job to process.
    ConverterJob _job;

    // The conversion statistics.
    ConversionStats _stats;

    // State for the printing of dots (every 100 packages).
    int _dotCounter;
    bool _dot;
    
    // Path to the per-thread mate FASTQ file.  Mates will be written interleaved, the left ones come first.
    seqan::CharString _matesFastqPath;
    // Path to the pre-thread singleton FASTQ file.
    seqan::CharString _singletonFastqPath;
    // Path to the pile files for left/right mates that are discordant.
    seqan::CharString _leftPileFastqPath, _rightPileFastqPath;

    // SequenceStreams for mates/singleton FASTQ files.
    std::shared_ptr<seqan::SequenceStream> _matesSeqStream, _singletonSeqStream;
    // SequenceStreams for left/right pile FASTQ files.
    std::shared_ptr<seqan::SequenceStream> _leftPileSeqStream, _rightPileSeqStream;

    // Stream for reading BAM file.
    std::shared_ptr<seqan::BamStream> _bamStream;
    // Index on BAM file.
    std::shared_ptr<seqan::BamIndex<seqan::Bai> > _baiIndex;

    ConverterThread() : _state(START), _threadId(-1), _queue(), _dotCounter(0), _dot(false)
    {}

    ConverterThread(int threadId, ConversionOptions const & options, JobQueue & queue) :
            _state(START), _threadId(threadId), _options(options), _queue(&queue), _dotCounter(0), _dot(false)
    {}

    // TODO(holtgrew): Delete temporary files on destruction.

    // Convert orphans from BAM file and append to the sequence streams.
    void convertOrphans();

    // Convert mapped files from the given job and append to the sequence streams.
    void convertMapped(ConverterJob const & job);

    // Begin first step, opening input and output files.
    //
    // Return status code for indicating success/failure.
    int _startFirstStep();

    // Update state for dot printing.
    void _updateDot()
    {
        if (_dotCounter > 100)
            _dotCounter = 0;
        _dotCounter += 1;
        if (_dotCounter > 100)
            _dot = true;
    }

    // Returns whether or not enough jobs were processed to print a dot and reset the state to "no dot" before
    // returning.
    bool dot()
    {
        bool res = _dot;
        _dot = false;
        return res;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function ConverterThread::convertMapped()
// ----------------------------------------------------------------------------

void ConverterThread::convertMapped(ConverterJob const & job)
{
    if (_state == START)
        _startFirstStep();
    _updateDot();
}

// ----------------------------------------------------------------------------
// Function ConverterThread::convertOrphans()
// ----------------------------------------------------------------------------

void ConverterThread::convertOrphans()
{
    if (_state == START)
        _startFirstStep();

    // Get shortcuts to streams and indices, also reduces overhead through dereferencing the shared_ptr.
    seqan::BamStream & bamStream = *_bamStream;
    seqan::BamIndex<seqan::Bai> & baiIndex = *_baiIndex;
    seqan::SequenceStream & matesStream = *_matesSeqStream;
    seqan::SequenceStream & singletonStream = *_singletonSeqStream;

    if (!jumpToOrphans(bamStream, baiIndex))
        return;  // TODO(holtgrew): Indicate failure.

    seqan::Dna5String leftSeq, rightSeq;

    seqan::BamAlignmentRecord leftRecord, rightRecord;
    while (!atEnd(bamStream))
    {
        if (readRecord(leftRecord, bamStream) != 0)
            return; // TODO(holtgrew): Indicate failure.
        leftSeq = leftRecord.seq;
        if (hasFlagRC(leftRecord))
            reverseComplement(leftSeq);
        SEQAN_ASSERT(hasFlagUnmapped(leftRecord));

        // Write out singletons immediately and continue.
        if (!hasFlagMultiple(leftRecord))
        {
            _stats.countSingletonOrphan();
            if (writeRecord(singletonStream, leftRecord.qName, leftSeq, leftRecord.qual) != 0)
                return;  // TODO(holtgrew): Indicate failure.
            continue;
        }

        // Read second mate if is paired.
        if (readRecord(rightRecord, bamStream) != 0)
            return;  // TODO(holtgrew): Indicate failure.
        rightSeq = rightRecord.seq;
        if (hasFlagRC(rightRecord))
            reverseComplement(rightSeq);
        SEQAN_ASSERT(hasFlagUnmapped(rightRecord));

        SEQAN_CHECK(leftRecord.qName == rightRecord.qName, "Adjacent orphan paired records must have the same QNAME!");

        // Write out first and right mate in the correct order.
        _stats.countPairedOrphan(2);
        if (hasFlagFirst(leftRecord))
        {
            if (writeRecord(matesStream, leftRecord.qName, leftSeq, leftRecord.qual) != 0)
                return;  // TODO(holtgrew): Indicate failure.
            if (writeRecord(matesStream, rightRecord.qName, rightSeq, rightRecord.qual) != 0)
                return;  // TODO(holtgrew): Indicate failure.
        }
        else
        {
            if (writeRecord(matesStream, rightRecord.qName, rightSeq, rightRecord.qual) != 0)
                return;  // TODO(holtgrew): Indicate failure.
            if (writeRecord(matesStream, leftRecord.qName, leftSeq, leftRecord.qual) != 0)
                return;  // TODO(holtgrew): Indicate failure.
        }
    }
}

// ----------------------------------------------------------------------------
// Function ConverterThread::_startFirstStep()
// ----------------------------------------------------------------------------

int ConverterThread::_startFirstStep()
{
    if (_state != START)
        return 1;  // Invalid state.

    // Open BAM file.
    _bamStream.reset(new seqan::BamStream());
    open(*_bamStream, toCString(_options.inputPath));
    if (!isGood(*_bamStream))
        return 1;

    // Open BAI file.
    seqan::CharString baiPath = _options.inputPath;
    append(baiPath, ".bai");
    _baiIndex.reset(new seqan::BamIndex<seqan::Bai>());
    if (read(*_baiIndex, toCString(baiPath)) != 0)
        return 1;

    // Generate temporary file name.
    _matesFastqPath = SEQAN_TEMP_FILENAME();
    append(_matesFastqPath, "_m.fq");  // Mate Pairs
    _singletonFastqPath = SEQAN_TEMP_FILENAME();
    append(_singletonFastqPath, "_s.fq");  // Singletons
    _leftPileFastqPath = SEQAN_TEMP_FILENAME();
    append(_leftPileFastqPath, "_lp.fq");  // Left Pile
    _rightPileFastqPath = SEQAN_TEMP_FILENAME();
    append(_rightPileFastqPath, "_rp.fq");  // Right Pile

    // Open temporary files.
    _matesSeqStream.reset(new seqan::SequenceStream());
    open(*_matesSeqStream, toCString(_matesFastqPath), seqan::SequenceStream::WRITE);
    if (!isGood(*_matesSeqStream))
        return 1;

    _singletonSeqStream.reset(new seqan::SequenceStream());
    open(*_singletonSeqStream, toCString(_singletonFastqPath), seqan::SequenceStream::WRITE);
    if (!isGood(*_singletonSeqStream))
        return 1;

    _leftPileSeqStream.reset(new seqan::SequenceStream());
    open(*_leftPileSeqStream, toCString(_singletonFastqPath), seqan::SequenceStream::WRITE);
    if (!isGood(*_leftPileSeqStream))
        return 1;

    _rightPileSeqStream.reset(new seqan::SequenceStream());
    open(*_rightPileSeqStream, toCString(_singletonFastqPath), seqan::SequenceStream::WRITE);
    if (!isGood(*_rightPileSeqStream))
        return 1;

    _state = FIRST_STEP;

    return 0;
}

#endif  // #ifndef SANDBOX_BAMBI_APPS_BAM2FASTQ_CONVERTER_THREAD_H_
