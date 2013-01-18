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
#include "sequence_sink.h"

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

    // The chunk length.
    unsigned chunkLength;

    // The verbosity.
    int verbosity;

    ConversionOptions() : maxTemplateLength(0), chunkLength(0), verbosity(1)
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
// Class ResultBuffer
// ----------------------------------------------------------------------------

// A buffer of left/right/single sequence ids/seqs/quals.
//
// It allows easy pushing of records into it and will classify the records accordingly.

class ResultBuffer
{
public:
    // The number of sequences in the buffer.
    unsigned numSeqs;

    // Buffer for left/right ids, quals, and sequences.
    TCharStringSet leftIds, leftQuals, rightIds, rightQuals, singleIds, singleQuals;
    TDna5StringSet leftSeqs, rightSeqs, singleSeqs;

    ResultBuffer() : numSeqs(0)
    {}

    // Write to sink.
    int writeToSink(SequenceSink & sink)
    {
        return sink.write(leftIds, leftSeqs, leftQuals, rightIds, rightSeqs, rightQuals,
                          singleIds, singleSeqs, singleQuals);
    }

    void clear()
    {
        numSeqs = 0;

        seqan::clear(leftIds);
        seqan::clear(leftSeqs);
        seqan::clear(leftQuals);
        seqan::clear(rightIds);
        seqan::clear(rightSeqs);
        seqan::clear(rightQuals);
        seqan::clear(singleIds);
        seqan::clear(singleSeqs);
        seqan::clear(singleQuals);
    }

    void insert(seqan::BamAlignmentRecord const & record)
    {
        numSeqs += 1;

        if (!hasFlagMultiple(record))
        {
            appendValue(singleIds, record.qName);
            appendValue(singleSeqs, record.seq);
            appendValue(singleQuals, record.qual);
            if (hasFlagRC(record))
                _applyRC(back(singleSeqs), back(singleQuals));
            return;
        }

        if (hasFlagFirst(record))
        {
            appendValue(leftIds, record.qName);
            appendValue(leftSeqs, record.seq);
            appendValue(leftQuals, record.qual);
            if (hasFlagRC(record))
                _applyRC(back(leftSeqs), back(leftQuals));
        }
        else
        {
            appendValue(rightIds, record.qName);
            appendValue(rightSeqs, record.seq);
            appendValue(rightQuals, record.qual);
            if (hasFlagRC(record))
                _applyRC(back(rightSeqs), back(rightQuals));
        }
    }

    void insert(seqan::CharString const & leftId, seqan::Dna5String const & leftSeq,
                seqan::CharString const & leftQual, seqan::CharString const & rightId,
                seqan::Dna5String const & rightSeq, seqan::CharString const & rightQual)
    {
        appendValue(leftIds, leftId);
        appendValue(leftSeqs, leftSeq);
        appendValue(leftQuals, leftQual);
        appendValue(rightIds, rightId);
        appendValue(rightSeqs, rightSeq);
        appendValue(rightQuals, rightQual);
    }

    // Apply reverse-complement to seq and reverse the qualities infix.
    template <typename TSeqInfix, typename TQualInfix>
    void _applyRC(TSeqInfix seq, TQualInfix quals) const
    {
        reverseComplement(seq);
        reverse(quals);
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
        FIRST_STEP   // In the first step, writing to temporary sequence files.
    } _state;

    // The id of the thread.  -1 for invalid.
    int _threadId;
    // The SequenceSink to write to.
    SequenceSink * _sink;
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

    // Stream for reading BAM file.
    std::shared_ptr<seqan::BamStream> _bamStream;
    // Index on BAM file.
    std::shared_ptr<seqan::BamIndex<seqan::Bai> > _baiIndex;

    // Buffer for ids, sequences and quality strings for left/right mates.
    ResultBuffer buffer;
    // The "pile" for reads having a too large template length.
    ResultBuffer pile;

    // The left-mapping mate cache data structure for BAM scanning.
    BamScannerCache scannerCache;

    ConverterThread() :
            _state(START), _threadId(-1), _sink(), _queue(), _dotCounter(0), _dot(false)
    {}

    ConverterThread(SequenceSink & sink, int threadId, ConversionOptions const & options, JobQueue & queue) :
            _state(START), _threadId(threadId), _sink(&sink), _options(options), _queue(&queue),
            _dotCounter(0), _dot(false)
    {}

    // Convert orphans from BAM file and append to the sequence streams.
    void convertOrphans();

    // Convert mapped files from the given job and append to the sequence streams.
    void convertMapped(ConverterJob const & job);

    // Begin first step, opening input and output files.
    //
    // Return status code for indicating success/failure.
    int _start();

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
    double tileStart = sysTime();
    if (_state == START)
        _start();

    // Get shortcuts to streams and indices, also reduces overhead through dereferencing the shared_ptr.
    seqan::BamStream & bamStream = *_bamStream;
    seqan::BamIndex<seqan::Bai> & baiIndex = *_baiIndex;

    seqan::BamAlignmentRecord record;
    seqan::Dna5String seq, otherSeq;
    record.rId = -1;

    int mtl = _options.maxTemplateLength;

    // Jump to the first record processed on the tile.  Return from function if a record is found that is right of the
    // tile and none on the tile.  A record to be processed on the tile is part of a pair that is completely on the tile
    // or that spans the tile and the one before.
    jumpToRegion(bamStream, job.rId, std::max(job.beginPos - mtl, 0), job.endPos, baiIndex);
    if (atEnd(bamStream))
        return;  // Done, no error.
    while (!atEnd(bamStream))
    {
        // Read record.
        if (readRecord(record, bamStream) != 0)
            return;  // TODO(holtgrew): Indicate error.

        // Skip record if not necessary.
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

        buffer.insert(record);
        if (!hasFlagMultiple(record))
            _stats.countSingletonMapped();
        else
            _stats.countPairedMapped();
    }

    // Iterate over the BAM stream until we are right of the tile.
    double writeTime = 0;
    while (true)  // We will break out before readRecord() on atEnd(bamStream).
    {
        // Write out chunk if it is over the size.
        if (buffer.numSeqs >= _options.chunkLength)
        {
            if (_options.verbosity >= 2)
            {
                SEQAN_OMP_PRAGMA(critical (io))
                {
                    std::cerr << "thread " << omp_get_thread_num() << " BEGINS writing to sink.\n";
                }
            }
            double writeStart = sysTime();
            if (buffer.writeToSink(*_sink) != 0)
                return;  // TOOD(holtgrew): Indicate failure.
            buffer.clear();
            writeTime += sysTime() - writeStart;
            if (_options.verbosity >= 2)
            {
                SEQAN_OMP_PRAGMA(critical (io))
                {
                    std::cerr << "thread " << omp_get_thread_num() << " STOPPED writing to sink.\n";
                }
            }
        }

        // Break if record right of this tile.  Skip if we are not to process this record for this tile.
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
            // Read next record.
            if (atEnd(bamStream))
                break;  // Stop if there are no more records.
            if (readRecord(record, bamStream) != 0)
                return;  // TODO(holtgrew): Indicate error.
            continue;
        }

        // Handle the case of singleton records first.
        if (!hasFlagMultiple(record))
        {
            buffer.insert(record);
            _stats.countSingletonMapped();

            // Read next record and go to beginning of loop.
            if (atEnd(bamStream))
                break;  // Stop if there are no more records.
            if (readRecord(record, bamStream) != 0)
                return;  // TODO(holtgrew): Indicate error.
            continue;
        }

        // Handle the case where record is handled through the piles.  This happens if the record does not map within
        // the maximum template length or on different contigs.
        if (record.rId != record.rNextId || abs(record.pos - record.pNext) > _options.maxTemplateLength)
        {
            pile.insert(record);
            _stats.countPairedMapped();
            // Read next record and go to beginning of loop.
            if (atEnd(bamStream))
                break;  // Stop if there are no more records.
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
            insertRecord(scannerCache, record);
            // Read next record and go to beginning of loop.
            if (atEnd(bamStream))
                break;  // Stop if there are no more records.
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
                SEQAN_OMP_PRAGMA(critical (io))
                {
                    std::cerr << "WARNING: Mate could not be found for:\n";
                    write2(std::cerr, record, bamStream.bamIOContext, seqan::Sam());
                    if (_options.verbosity >= 2)
                    {
                        std::cerr << "job.beginPos = " << job.beginPos << "\n"
                                  << "job.beginPos - mtl = " << job.beginPos - mtl << "\n"
                                  << "job.endPos = " << job.endPos << "\n";
                    }
                }
                //SEQAN_FAIL("Could not find mate for record but it must be there!");
            }
            else
            {
                seqan::BamAlignmentRecord const & otherRecord = otherIt->second;

                // Write out sequences into the result stream.
                _stats.countPairedMapped(2);
                buffer.insert(record);
                buffer.insert(otherRecord);

                // Remove other entry from cache again.
                unsigned oldSize = scannerCache._map.size();
                erase(scannerCache, otherIt);
                SEQAN_CHECK(scannerCache._map.size() + 1 == oldSize, "Must decrease size!");
                SEQAN_CHECK(!containsMate(scannerCache, record), "Must not be there any more!");
            }
        }

        // Read next record.
        if (atEnd(bamStream))
            break;  // Stop if there are no more records.
        if (readRecord(record, bamStream) != 0)
            return;  // TODO(holtgrew): Indicate error.
    }

    if (!empty(scannerCache))
    {
        SEQAN_OMP_PRAGMA(critical (io))
        {
            std::cerr << "WARNING: Did not find mate for some records in the cache.\n"
                      << "The reason could be erroneous data in the BAM file or a bug in the program.\n"
                      << "\n"
                      << "The problematic records are as follows:\n"
                      << "<begin>\n";
            for (auto it = scannerCache._map.begin(); it != scannerCache._map.end(); ++it)
                write2(std::cerr, it->second, bamStream.bamIOContext, seqan::Sam());
            std::cerr << "<end>\n";
        }
    }
    clear(scannerCache);

    double writeStart = sysTime();
    if (_options.verbosity >= 2)
    {
        SEQAN_OMP_PRAGMA(critical (io))
        {
            std::cerr << "thread " << omp_get_thread_num() << " BEGINS writing to sink.\n";
        }
    }
    if (buffer.writeToSink(*_sink) != 0)
        return;  // TOOD(holtgrew): Indicate failure.
    buffer.clear();
    writeTime += sysTime() - writeStart;
    if (_options.verbosity >= 2)
    {
        SEQAN_OMP_PRAGMA(critical (io))
        {
            std::cerr << "thread " << omp_get_thread_num() << " STOPPED writing to sink.\n";
            double now = sysTime();
            fprintf(stderr, "THREAD %d PROCESSED TILE (%d, %d, %d).\tTOTAL: %f s\tWRITING: %f s\n", omp_get_thread_num(), job.rId, job.beginPos, job.endPos, now - tileStart, writeTime);
        }
    }

    _updateDot();
}

// ----------------------------------------------------------------------------
// Function ConverterThread::convertOrphans()
// ----------------------------------------------------------------------------

void ConverterThread::convertOrphans()
{
    if (_state == START)
        _start();

    if (_options.verbosity >= 2)
    {
        SEQAN_OMP_PRAGMA(critical (io))
        {
            std::cerr << "thread " << omp_get_thread_num() << " begins converting orphans.\n";
        }
    }

    // Get shortcuts to streams and indices, also reduces overhead through dereferencing the shared_ptr.
    seqan::BamStream & bamStream = *_bamStream;
    seqan::BamIndex<seqan::Bai> & baiIndex = *_baiIndex;

    if (!jumpToOrphans(bamStream, baiIndex))
        return;  // TODO(holtgrew): Indicate failure.

    // Used for storing the sequence below.  Declare here to allow reusing.
    seqan::BamAlignmentRecord record;

    // Scan over orphans in BAM file.
    while (!atEnd(bamStream))
    {
        // Write out chunk if it is over the size.
        if (buffer.numSeqs >= _options.chunkLength)
        {
            if (_options.verbosity >= 2)
            {
                SEQAN_OMP_PRAGMA(critical (io))
                {
                    std::cerr << "thread " << omp_get_thread_num() << " BEGINS writing to sink.\n";
                }
            }
            if (buffer.writeToSink(*_sink) != 0)
                return;  // TOOD(holtgrew): Indicate failure.
            buffer.clear();
            if (_options.verbosity >= 2)
            {
                SEQAN_OMP_PRAGMA(critical (io))
                {
                    std::cerr << "thread " << omp_get_thread_num() << " STOPPED writing to sink.\n";
                }
            }
        }

        // Read first/single record.
        if (readRecord(record, bamStream) != 0)
            return; // TODO(holtgrew): Indicate failure.
        SEQAN_ASSERT(hasFlagUnmapped(record));

        // Append singletons to buffer and read next.
        if (!hasFlagMultiple(record))
        {
            buffer.insert(record);
            continue;
        }

        // Insert first mate.
        SEQAN_ASSERT(hasFlagMultiple(record));
        buffer.insert(record);

        // Read second mate.
        if (readRecord(record, bamStream) != 0)
            return; // TODO(holtgrew): Indicate failure.
        SEQAN_ASSERT(hasFlagUnmapped(record));
        SEQAN_ASSERT(hasFlagMultiple(record));
        buffer.insert(record);

        SEQAN_ASSERT_EQ(length(buffer.leftIds), length(buffer.rightIds));
    }

    // Write remaining records to output file.
    if (_options.verbosity >= 2)
    {
        SEQAN_OMP_PRAGMA(critical (io))
        {
            std::cerr << "thread " << omp_get_thread_num() << " BEGINS writing to sink.\n";
        }
    }
    if (buffer.writeToSink(*_sink) != 0)
        return;  // TOOD(holtgrew): Indicate failure.
    buffer.clear();
    if (_options.verbosity >= 2)
    {
        SEQAN_OMP_PRAGMA(critical (io))
        {
            std::cerr << "thread " << omp_get_thread_num() << " STOPPED writing to sink.\n";
        }
    }

    if (_options.verbosity >= 2)
    {
        SEQAN_OMP_PRAGMA(critical (io))
        {
            std::cerr << "thread " << omp_get_thread_num() << " is done converting orphans.\n";
        }
    }
}

// ----------------------------------------------------------------------------
// Function ConverterThread::_start()
// ----------------------------------------------------------------------------

int ConverterThread::_start()
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

    _state = FIRST_STEP;

    return 0;
}

#endif  // #ifndef SANDBOX_BAMBI_APPS_BAM2FASTQ_CONVERTER_THREAD_H_
