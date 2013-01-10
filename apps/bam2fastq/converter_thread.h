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

#include <seqan/sequence.h>
#include <seqan/store.h>  // For NameStoreCache.

#include "job_queue.h"
#include "bam_scanner_cache.h"

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
    int numMapped;
    int numSingletonMapped;
    int numPairedMapped;

    ConversionStats() :
            numOrphans(0), numSingletonOrphans(0), numPairedOrphans(0), numMapped(0), numSingletonMapped(0),
            numPairedMapped(0)
    {}

    void countSingletonOrphan(int num = 1)
    {
        numOrphans += num;
        numSingletonOrphans += num;
    }

    void countPairedOrphan(int num = 1)
    {
        numOrphans += num;
        numPairedOrphans += num;
    }

    void countSingletonMapped(int num = 1)
    {
        numMapped += num;
        numSingletonMapped += num;
    }

    void countPairedMapped(int num = 1)
    {
        numMapped += num;
        numPairedMapped += num;
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

    // The left-mapping mate cache data structure for BAM scanning.
    BamScannerCache scannerCache;

    ConverterThread() :
            _state(START), _threadId(-1), _queue(), _dotCounter(0), _dot(false)
    {}

    ConverterThread(int threadId, ConversionOptions const & options, JobQueue & queue) :
            _state(START), _threadId(threadId), _options(options), _queue(&queue), _dotCounter(0), _dot(false)
    {}

    ~ConverterThread()
    {
        cleanup();
    }

    // Remove temporary files again.
    void cleanup()
    {
        remove(toCString(_matesFastqPath));
        remove(toCString(_singletonFastqPath));
        remove(toCString(_leftPileFastqPath));
        remove(toCString(_rightPileFastqPath));
    }

    // Convert orphans from BAM file and append to the sequence streams.
    void convertOrphans();

    // Convert mapped files from the given job and append to the sequence streams.
    void convertMapped(ConverterJob const & job);

    // Write resulting pairs interleaved to the given SequenceStream opened for output.
    int writeResult(seqan::SequenceStream & peOut, seqan::SequenceStream & singletonOut);

    // Begin first step, opening input and output files.
    //
    // Return status code for indicating success/failure.
    int _startFirstStep();

    // Begin second step, re-opening output files for reading.
    int _startSecondStep();

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

    // Get shortcuts to streams and indices, also reduces overhead through dereferencing the shared_ptr.
    seqan::BamStream & bamStream = *_bamStream;
    seqan::BamIndex<seqan::Bai> & baiIndex = *_baiIndex;
    seqan::SequenceStream & matesStream = *_matesSeqStream;
    seqan::SequenceStream & singletonStream = *_singletonSeqStream;
    seqan::SequenceStream & leftPileStream = *_leftPileSeqStream;
    seqan::SequenceStream & rightPileStream = *_rightPileSeqStream;

    seqan::BamAlignmentRecord record;
    seqan::Dna5String seq, otherSeq;
    record.rId = -1;

    int mtl = _options.maxTemplateLength;

    bool DEBUG = false;
    bool FD = false;

    if (FD || DEBUG)
    {
        SEQAN_OMP_PRAGMA(critical (io))
        {
            std::cerr << "--------------------- NEW TILE ---------------------\n"
                      << "  begin pos - mtl\t" << job.beginPos - mtl << "\n"
                      << "  begin pos      \t" << job.beginPos << "\n"
                      << "  end pos        \t" << job.endPos << "\n";
        }
    }

    // Jump to the first record processed on the tile.  Return from function if a record is found that is right of the
    // tile and none on the tile.  A record to be processed on the tile is part of a pair that is completely on the tile
    // or that spans the tile and the one before.
    jumpToPos(bamStream, job.rId, std::max(job.beginPos - mtl, 0), baiIndex);
    if (atEnd(bamStream))
        return;  // Done, no error.
    while (!atEnd(bamStream))
    {
        if (readRecord(record, bamStream) != 0)
            return;  // TODO(holtgrew): Indicate error.
        if (record.rId > job.rId || (record.rId == job.rId && record.pos >= job.endPos))
            return;  // Done, no record in tile in BAM file.

        if (hasFlagSecondary(record))
            continue;  // Skip, we only want the primary mapping location.
        if (record.rId < job.rId)
            continue;  // Skip if on contig left of the current one.
        // If we reach here then the record is on the same contig.

        // Stop if the whole pair is on the tile.
        if (record.pos >= job.beginPos && record.pos < job.endPos &&
            record.pNext >= job.beginPos && record.pNext < job.endPos)
            break;
        // Stop if the pair spans the tile to the left and does not have a too long insert size.
        if (record.pos < job.beginPos && record.pNext >= job.beginPos &&
            abs(record.pos - record.pNext) <= _options.maxTemplateLength)
            break;
        // Stop if the read is placed on this tile and the pair needs special handling below because of too long
        // template length.
        if (record.pos >= job.beginPos &&
            abs(record.pos - record.pNext) > _options.maxTemplateLength)
            break;

        // Otherwise, we have to look at the next.
    }
    if (FD)
    {
        SEQAN_OMP_PRAGMA(critical (io))
        {
            std::cerr << "FIRST ON TILE IS\n";
            write2(std::cerr, record, bamStream.bamIOContext, seqan::Sam());
        }
    }
    // The first BAM record on the tile is now stored in record or it is the first after.

    // If this was already the last record, then put it onto the pile or into the singleton stream.
    if (atEnd(bamStream))
    {
        // Restore original sequence and quality strings.
        seq = record.seq;
        if (hasFlagRC(record))
        {
            reverseComplement(seq);
            reverse(record.qual);
        }

        if (!hasFlagMultiple(record))
        {
            if (DEBUG)
            {
                SEQAN_OMP_PRAGMA(critical (io))
                {
                    write2(std::cerr, record, bamStream.bamIOContext, seqan::Sam());
                    std::cerr << " => singleton stream\n";
                }
            }
            _stats.countSingletonMapped();
            if (writeRecord(singletonStream, record.qName, seq, record.qual) != 0)
                return;  // TODO(holtgrew): Indicate failure.
        }
        else
        {
            if (hasFlagFirst(record))
            {
                if (DEBUG)
                {
                    SEQAN_OMP_PRAGMA(critical (io))
                    {
                        write2(std::cerr, record, bamStream.bamIOContext, seqan::Sam());
                        std::cerr << " => left pile stream\n";
                    }
                }
                _stats.countPairedMapped();
                if (writeRecord(leftPileStream, record.qName, seq, record.qual) != 0)
                    return;  // TODO(holtgrew): Indicate failure.
            }
            else
            {
                if (DEBUG)
                {
                    SEQAN_OMP_PRAGMA(critical (io))
                    {
                        write2(std::cerr, record, bamStream.bamIOContext, seqan::Sam());
                        std::cerr << " => right pile stream\n";
                    }
                }
                _stats.countPairedMapped();
                if (writeRecord(rightPileStream, record.qName, seq, record.qual) != 0)
                    return;  // TODO(holtgrew): Indicate failure.
            }
        }
    }

    // Iterate over the BAM stream until we are right of the tile.
    while (!atEnd(bamStream))
    {
        if (record.pos >= job.endPos || record.rId != job.rId)
            break;  // Break out if not on this tile.

        bool skip = false;
        if (hasFlagSecondary(record))
            skip = true;  // Skip secondary mappings.
        if (record.pos < job.beginPos && record.pNext < job.beginPos)
            skip = true;  // Skip if whole pair left of tile.
        if (record.pNext >= job.endPos)
            skip = true;  // Skip if pair hangs over tile to the right.
        if (record.pos >= job.beginPos &&
            abs(record.pos - record.pNext) > _options.maxTemplateLength)
            skip = false;  // Do NOT skip if needs special treatment.
        if (skip)
        {
            // Skip secondary mappings.
            if (readRecord(record, bamStream) != 0)
                return;  // TODO(holtgrew): Indicate error.
            continue;
        }

        if (FD)
        {
            if (record.qName == "chr17_3060049_3060217_0:0:0_1:0:0_2ffad")
            {
                write2(std::cerr, record, bamStream.bamIOContext, seqan::Sam());
                std::cerr << " => BEING PROCESSED\n";
            }
        }

        // Restore original sequence and quality strings.
        seq = record.seq;
        if (hasFlagRC(record))
        {
            reverseComplement(seq);
            reverse(record.qual);
        }

        // Handle the case of singleton records first.
        if (!hasFlagMultiple(record))
        {
            if (DEBUG)
            {
                SEQAN_OMP_PRAGMA(critical (io))
                {
                    write2(std::cerr, record, bamStream.bamIOContext, seqan::Sam());
                    std::cerr << " => singleton stream\n";
                }
            }
            _stats.countSingletonMapped();
            if (writeRecord(singletonStream, record.qName, seq, record.qual) != 0)
                return;  // TODO(holtgrew): Indicate failure.
            // Read next record and go to beginning of loop.
            if (readRecord(record, bamStream) != 0)
                return;  // TODO(holtgrew): Indicate error.
            continue;
        }

        // Handle the case where record is handled through the piles.  This happens if the record does not map within
        // the maximum template length or on different contigs.
        if (record.rId != record.rNextId || abs(record.pos - record.pNext) > _options.maxTemplateLength)
        {
            if (hasFlagFirst(record))
            {
                if (DEBUG)
                {
                    SEQAN_OMP_PRAGMA(critical (io))
                    {
                        write2(std::cerr, record, bamStream.bamIOContext, seqan::Sam());
                        std::cerr << " => left pile stream\n";
                    }
                }
                _stats.countPairedMapped();
                if (writeRecord(leftPileStream, record.qName, seq, record.qual) != 0)
                    return;  // TODO(holtgrew): Indicate failure.
            }
            else
            {
                if (DEBUG)
                {
                    SEQAN_OMP_PRAGMA(critical (io))
                    {
                        write2(std::cerr, record, bamStream.bamIOContext, seqan::Sam());
                        std::cerr << " => right pile stream\n";
                    }
                }
                _stats.countPairedMapped();
                if (writeRecord(rightPileStream, record.qName, seq, record.qual) != 0)
                    return;  // TODO(holtgrew): Indicate failure.
            }
            // Read next record and go to beginning of loop.
            if (readRecord(record, bamStream) != 0)
                return;  // TODO(holtgrew): Indicate error.
            continue;
        }

        // If we reach here then the record is a mate in a pair that has a sufficiently short template length and both
        // mates map on the same contig.  If this is the left mapping mate then we simply insert record into the
        // BamScannerCache.  If this is the right mapping mate then we have to be able to find the left mapping mate in
        // the BamScannerCache.  In the case of mapping a the same position, we put record into the cache if there is no
        // such record in the cache already and we must be the first.
        SEQAN_ASSERT_EQ(record.rId, record.rNextId);
        if (record.pos < record.pNext || (record.pos == record.pNext && !containsMate(scannerCache, record)))
        {
            if (DEBUG)
            {
                SEQAN_OMP_PRAGMA(critical (io))
                {
                    write2(std::cerr, record, bamStream.bamIOContext, seqan::Sam());
                    std::cerr << " => into cache\n";
                }
            }
            insertRecord(scannerCache, record);
            // Read next record and go to beginning of loop.
            if (readRecord(record, bamStream) != 0)
                return;  // TODO(holtgrew): Indicate error.
            continue;
        }
        else
        {
            typedef seqan::Iterator<BamScannerCache, seqan::Standard>::Type TIter;
            TIter otherIt = findMate(scannerCache, record);
            if (otherIt == end(scannerCache, seqan::Standard()))
            {
                if (DEBUG)
                {
                    SEQAN_OMP_PRAGMA(critical (io))
                    {
                        std::cerr << "MATE COULD NOT BE FOUND FOR:\n";
                        write2(std::cerr, record, bamStream.bamIOContext, seqan::Sam());
                        std::cerr << "job.beginPos = " << job.beginPos << "\n"
                                  << "job.beginPos - mtl = " << job.beginPos - mtl << "\n"
                                  << "job.endPos = " << job.endPos << "\n";
                    }
                }
                SEQAN_FAIL("Could not find mate for record but it must be there!");
            }
            seqan::BamAlignmentRecord const & otherRecord = otherIt->second;
            // Restore original sequence for the other record.
            otherSeq = otherRecord.seq;
            if (hasFlagRC(otherRecord))
            {
                reverseComplement(otherSeq);
                reverse(otherRecord.qual);
            }

            // Write out sequences into the result stream.
            _stats.countPairedMapped(2);
            if (hasFlagFirst(record))
            {
                if (DEBUG)
                {
                    SEQAN_OMP_PRAGMA(critical (io))
                    {
                        write2(std::cerr, record, bamStream.bamIOContext, seqan::Sam());
                        std::cerr << " => mates stream (as left)\n";
                        write2(std::cerr, otherRecord, bamStream.bamIOContext, seqan::Sam());
                        std::cerr << " => mates stream (as right)\n";
                    }
                }
                if (writeRecord(matesStream, record.qName, seq, record.qual) != 0)
                    return;  // TODO(holtgrew): Indicate failure!
                if (writeRecord(matesStream, otherRecord.qName, otherSeq, otherRecord.qual) != 0)
                    return;  // TODO(holtgrew): Indicate failure!
            }
            else
            {
                if (DEBUG)
                {
                    SEQAN_OMP_PRAGMA(critical (io))
                    {
                        write2(std::cerr, otherRecord, bamStream.bamIOContext, seqan::Sam());
                        std::cerr << " => mates stream (as left)\n";
                        write2(std::cerr, record, bamStream.bamIOContext, seqan::Sam());
                        std::cerr << " => mates stream (as right)\n";
                    }
                }
                if (writeRecord(matesStream, otherRecord.qName, otherSeq, otherRecord.qual) != 0)
                    return;  // TODO(holtgrew): Indicate failure!
                if (writeRecord(matesStream, record.qName, seq, record.qual) != 0)
                    return;  // TODO(holtgrew): Indicate failure!
            }
            
            // Remove other entry from cache again.
            unsigned oldSize = scannerCache._map.size();
            erase(scannerCache, otherIt);
            SEQAN_CHECK(scannerCache._map.size() + 1 == oldSize, "Must decrease size!");
            SEQAN_CHECK(!containsMate(scannerCache, record), "Must not be there any more!");
        }

        // Read next record.
        if (readRecord(record, bamStream) != 0)
            return;  // TODO(holtgrew): Indicate error.
    }

    if (!empty(scannerCache))
        dumpCache(scannerCache);
    SEQAN_CHECK(empty(scannerCache), "Invalid pairing, probably some POS/PNEXT, RID/RNEXT field is wrong!");

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

    // Used for storing the sequence below.  Declare here to allow reusing.
    seqan::Dna5String leftSeq, rightSeq;
    seqan::BamAlignmentRecord leftRecord, rightRecord;

    // Scan over orphans in BAM file.
    while (!atEnd(bamStream))
    {
        if (readRecord(leftRecord, bamStream) != 0)
            return; // TODO(holtgrew): Indicate failure.
        leftSeq = leftRecord.seq;
        if (hasFlagRC(leftRecord))
        {
            reverseComplement(leftSeq);
            reverse(leftRecord.qual);
        }
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
        {
            reverseComplement(rightSeq);
            reverse(rightRecord.qual);
        }
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
    SEQAN_OMP_PRAGMA(critical(temp_filename))
    {
        _leftPileFastqPath = SEQAN_TEMP_FILENAME();
        _rightPileFastqPath = SEQAN_TEMP_FILENAME();
        _matesFastqPath = SEQAN_TEMP_FILENAME();
        _singletonFastqPath = SEQAN_TEMP_FILENAME();
    }

    append(_matesFastqPath, "_m.fq");  // Mate Pairs
    append(_singletonFastqPath, "_s.fq");  // Singletons
    append(_leftPileFastqPath, "_lp.fq");  // Left Pile
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
    open(*_leftPileSeqStream, toCString(_leftPileFastqPath), seqan::SequenceStream::WRITE);
    if (!isGood(*_leftPileSeqStream))
        return 1;

    _rightPileSeqStream.reset(new seqan::SequenceStream());
    open(*_rightPileSeqStream, toCString(_rightPileFastqPath), seqan::SequenceStream::WRITE);
    if (!isGood(*_rightPileSeqStream))
        return 1;

    _state = FIRST_STEP;

    return 0;
}

// ----------------------------------------------------------------------------
// Function ConverterThread::_startSecondStep()
// ----------------------------------------------------------------------------

int ConverterThread::_startSecondStep()
{
    if (_state != FIRST_STEP)
        return 1;  // Invalid state.

    // Re-open temporary files for reading.
    _matesSeqStream.reset(new seqan::SequenceStream());
    open(*_matesSeqStream, toCString(_matesFastqPath));
    if (!isGood(*_matesSeqStream))
        return 1;

    _singletonSeqStream.reset(new seqan::SequenceStream());
    open(*_singletonSeqStream, toCString(_singletonFastqPath));
    if (!isGood(*_singletonSeqStream))
        return 1;

    _leftPileSeqStream.reset(new seqan::SequenceStream());
    open(*_leftPileSeqStream, toCString(_leftPileFastqPath));
    if (!isGood(*_leftPileSeqStream))
        return 1;

    _rightPileSeqStream.reset(new seqan::SequenceStream());
    open(*_rightPileSeqStream, toCString(_rightPileFastqPath));
    if (!isGood(*_rightPileSeqStream))
        return 1;

    _state = SECOND_STEP;

    return 0;
}

// ----------------------------------------------------------------------------
// Function ConverterThread::writeResult()
// ----------------------------------------------------------------------------

int ConverterThread::writeResult(seqan::SequenceStream & matesOut, seqan::SequenceStream & singletonOut)
{
    // Start second step if necessary.
    if (_state != SECOND_STEP && _startSecondStep() != 0)
    {
        std::cerr << "ERROR: Could not start second step for thread " << omp_get_thread_num() << "\n";
        return 1;
    }

    seqan::SequenceStream & matesStream = *_matesSeqStream;
    seqan::SequenceStream & singletonStream = *_singletonSeqStream;

    // Buffers used below.
    seqan::CharString id, seq, qual;
    
    // Write out all singletons.
    while (!atEnd(singletonStream))
    {
        if (readRecord(id, seq, qual, singletonStream) != 0)
        {
            std::cerr << "ERROR: Could not load singleton from temporary file for thread " << omp_get_thread_num() << "\n";
            return 1;
        }
        if (writeRecord(singletonOut, id, seq, qual) != 0)
        {
            std::cerr << "ERROR: Could not write singleton record for thread " << omp_get_thread_num() << "\n";
            return 1;
        }
    }

    // Write out all PE reads.
    while (!atEnd(matesStream))
    {
        if (readRecord(id, seq, qual, matesStream) != 0)
        {
            std::cerr << "ERROR: Could not load paired from temporary file for thread " << omp_get_thread_num() << "\n";
            return 1;
        }
        if (writeRecord(matesOut, id, seq, qual) != 0)
        {
            std::cerr << "ERROR: Could not write paired record for thread " << omp_get_thread_num() << "\n";
            return 1;
        }
    }

    return 0;
}

#endif  // #ifndef SANDBOX_BAMBI_APPS_BAM2FASTQ_CONVERTER_THREAD_H_
