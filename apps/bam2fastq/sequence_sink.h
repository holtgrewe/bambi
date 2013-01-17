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

#include <zlib.h>

#include <omp.h>

// TODO(holtgrew): This could be made somewhat cleaner.

// ============================================================================
// Forwards
// ============================================================================

bool compressData(seqan::CharString & outBuffer, seqan::CharString const & inBuffer);

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

    // Sequence stream objects for writing out the sequences when not using any compression.
    seqan::SequenceStream outLeft, outRight, outSingle;

    // The memory mapped string objects to use when writing out with compression.
    seqan::String<char, seqan::MMap<> > mmapLeft, mmapRight, mmapSingle;

    // The buffers to use.
    seqan::String<seqan::CharString> textBuffersL, textBuffersR, textBuffersS;
    seqan::String<seqan::CharString> compressedBuffersL, compressedBuffersR, compressedBuffersS;

    // Whether or not we have already opened the sequence streams.
    bool isOpenLeft, isOpenRight, isOpenSingle;

    // Whether or not to write out paired-end reads in an interleaved fashion.
    bool interleaved;
    // Whether or not to compress the output files.
    bool gzip;
    // Number of threads using this sink.
    int numThreads;
    // Verbosity level.
    int verbosity;

    // The names of the left/right/single output files.
    seqan::CharString pathLeft, pathRight, pathSingle;

    // Whether or not we have already opened the sequence streams.
    bool isOpenLeft, isOpenRight, isOpenSingle;

    // The OpenMP lock for thread safety.
    omp_lock_t _lock;

    SequenceSink(bool interleaved, bool gzip, int numThreads, char const * pathLeft, char const * pathRight,
                 char const * pathSingle, int verbosity) :
            pathLeft(pathLeft), pathRight(pathRight), pathSingle(pathSingle), isOpenLeft(false),
            isOpenRight(false), isOpenSingle(false), interleaved(interleaved), gzip(gzip), numThreads(numThreads),
            verbosity(verbosity)
    {
        omp_init_lock(&_lock);
        resize(textBuffersL, numThreads);
        resize(textBuffersR, numThreads);
        resize(textBuffersS, numThreads);
        resize(compressedBuffersL, numThreads);
        resize(compressedBuffersR, numThreads);
        resize(compressedBuffersS, numThreads);
    }

    ~SequenceSink()
    {
        omp_destroy_lock(&_lock);
    }

    int write(TCharStringSet const & leftIds, TDna5StringSet const & leftSeqs, TCharStringSet const & leftQuals,
              TCharStringSet const & rightIds, TDna5StringSet const & rightSeqs, TCharStringSet const & rightQuals,
              TCharStringSet const & singleIds, TDna5StringSet const & singleSeqs, TCharStringSet const & singleQuals)
    {
        if (gzip)
            return _writeCompressed(leftIds, leftSeqs, leftQuals, rightIds, rightSeqs, rightQuals,
                                    singleIds, singleSeqs, singleQuals);
        else  // plain text
            return _writeUncompressed(leftIds, leftSeqs, leftQuals, rightIds, rightSeqs, rightQuals,
                                      singleIds, singleSeqs, singleQuals);
    }

    // Write out compressed data.  The compression is done in parallel and only for writing out the gzip block do we
    // obtain the lock.
    int _writeCompressed(
            TCharStringSet const & leftIds, TDna5StringSet const & leftSeqs, TCharStringSet const & leftQuals,
            TCharStringSet const & rightIds, TDna5StringSet const & rightSeqs, TCharStringSet const & rightQuals,
            TCharStringSet const & singleIds, TDna5StringSet const & singleSeqs, TCharStringSet const & singleQuals)
    {
        SEQAN_ASSERT_EQ(length(leftIds), length(rightIds));
        int tid = omp_get_thread_num();

        // Write single reads to FASTQ buffer.
        for (unsigned i = 0; i < length(singleIds); ++i)
            if (writeRecord(textBuffersS[tid], singleIds[i], singleSeqs[i], singleQuals[i], seqan::Fastq()) != 0)
                return 1;  // Error writing.

        // Write paired reads to FASTQ buffer(s).
        seqan::CharString * bufferL = &textBuffersL[tid];
        seqan::CharString * bufferR = interleaved ? &textBuffersL[tid] : &textBuffersR[tid];
        for (unsigned i = 0; i < length(leftIds); ++i)
        {
            if (writeRecord(*bufferL, leftIds[i], leftSeqs[i], leftQuals[i], seqan::Fastq()) != 0)
                return 1;  // Error writing.
            if (writeRecord(*bufferR, rightIds[i], rightSeqs[i], rightQuals[i], seqan::Fastq()) != 0)
                return 1;  // Error writing.
        }

        // Compress single buffer.
        if (!empty(textBuffersS[tid]))
        {
            if (!compressData(compressedBuffersS[tid], textBuffersS[tid]))
            {
                SEQAN_OMP_PRAGMA(critical (io))
                {
                    std::cerr << "Compression of single reads failed for thread " << tid << "\n";
                }
                return 1;
            }
            if (verbosity >= 2)
            {
                SEQAN_OMP_PRAGMA(critical (io))
                {
                    std::cerr << "Thread " << tid << " compressed " << length(textBuffersS[tid])
                              << " bytes to " << length(compressedBuffersS[tid]) << "\n";
                }
            }
        }

        // Compress paired buffer(s).
        if (!empty(textBuffersL[tid]))
        {
            if (!compressData(compressedBuffersL[tid], textBuffersL[tid]))
            {
                SEQAN_OMP_PRAGMA(critical (io))
                {
                    std::cerr << "Compression of left reads failed for thread " << tid << "\n";
                }
                return 1;
            }
            if (verbosity >= 2)
            {
                SEQAN_OMP_PRAGMA(critical (io))
                {
                    std::cerr << "Thread " << tid << " compressed " << length(textBuffersL[tid])
                              << " bytes to " << length(compressedBuffersL[tid]) << "\n";
                }
            }
        }
        if (!empty(textBuffersR[tid]))
        {
            if (!compressData(compressedBuffersR[tid], textBuffersR[tid]))
            {
                SEQAN_OMP_PRAGMA(critical (io))
                {
                    std::cerr << "Compression of right reads failed for thread " << tid << "\n";
                }
                return 1;
            }
            if (verbosity >= 2)
            {
                SEQAN_OMP_PRAGMA(critical (io))
                {
                    std::cerr << "Thread " << tid << " compressed " << length(textBuffersR[tid])
                              << " bytes to " << length(compressedBuffersR[tid]) << "\n";
                }
            }
        }

        // Obtain lock and write out compressed data.
        {
            RaiiLock lockGuard(_lock);
            if (!empty(compressedBuffersS[tid]))
            {
                if (verbosity >= 2)
                {
                    SEQAN_OMP_PRAGMA(critical (io))
                    {
                        std::cerr << "Thread " << tid << " writes " << length(compressedBuffersS[tid])
                                  << " bytes to " << pathSingle << "\n";
                    }
                }
                if (!isOpenSingle && _openCompressed(mmapSingle, isOpenSingle, pathSingle) != 0)
                    return 1;  // Could not open single output file.
                append(mmapSingle, compressedBuffersS[tid]);
                // append(mmapSingle, textBuffersS[tid]);
            }

            if (!empty(compressedBuffersL[tid]))
            {
                if (verbosity >= 2)
                {
                    SEQAN_OMP_PRAGMA(critical (io))
                    {
                        std::cerr << "Thread " << tid << " writes " << length(compressedBuffersL[tid])
                                  << " bytes to " << pathLeft << "\n";
                    }
                }
                if (!isOpenLeft && _openCompressed(mmapLeft, isOpenLeft, pathLeft) != 0)
                    return 1;  // Could not open left output file.
                append(mmapLeft, compressedBuffersL[tid]);
                // unsigned oldLength = length(mmapLeft);
                // resize(mmapLeft, oldLength + length(compressedBuffersL[tid]));
                // std::copy(begin(compressedBuffersL[tid], seqan::Standard()), end(compressedBuffersL[tid], seqan::Standard()),
                //           begin(mmapLeft, seqan::Standard()) + oldLength);
                // close(mmapLeft);
                // append(mmapLeft, textBuffersL[tid]);
            }

            if (!empty(compressedBuffersR[tid]))
            {
                if (verbosity >= 2)
                {
                    SEQAN_OMP_PRAGMA(critical (io))
                    {
                        std::cerr << "Thread " << tid << " writes " << length(compressedBuffersR[tid])
                                  << " bytes to " << pathRight << "\n";
                    }
                }
                if (!isOpenRight && _openCompressed(mmapRight, isOpenRight, pathRight) != 0)
                    return 1;  // Could not open right output file.
                append(mmapRight, compressedBuffersR[tid]);
                // append(mmapRight, textBuffersR[tid]);
            }
        }
        

        // Clear buffers again.
        clear(textBuffersL[tid]);
        clear(textBuffersR[tid]);
        clear(textBuffersS[tid]);
        clear(compressedBuffersL[tid]);
        clear(compressedBuffersR[tid]);
        clear(compressedBuffersS[tid]);
        return 0;
    }

    // Writing out uncompressed data.
    int _writeUncompressed(
            TCharStringSet const & leftIds, TDna5StringSet const & leftSeqs, TCharStringSet const & leftQuals,
            TCharStringSet const & rightIds, TDna5StringSet const & rightSeqs, TCharStringSet const & rightQuals,
            TCharStringSet const & singleIds, TDna5StringSet const & singleSeqs, TCharStringSet const & singleQuals)
    {
        SEQAN_ASSERT_EQ(length(leftIds), length(rightIds));

        // Lock data structure.
        RaiiLock lockGuard(_lock);

        // Write out to single files.
        if (!empty(singleIds))
        {
            if (!isOpenSingle && _openUncompressed(outSingle, isOpenSingle, pathSingle) != 0)
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
            if (!isOpenLeft && _openUncompressed(outLeft, isOpenLeft, pathLeft) != 0)
                return 1;
            if (writeAll(outLeft, leftIds, leftSeqs, leftQuals) != 0)
            {
                SEQAN_OMP_PRAGMA(critical (io))
                {
                    std::cerr << "ERROR: Could not write to left output file!\n";
                }
                return 1;
            }
            if (!isOpenRight && _openUncompressed(outRight, isOpenRight, pathRight) != 0)
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
            if (!isOpenLeft && _openUncompressed(outLeft, isOpenLeft, pathLeft) != 0)
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

    int _openCompressed(seqan::String<char, seqan::MMap<> > & out, bool & flag, seqan::CharString const & path)
    {
        if (!open(out, toCString(path)))//, seqan::OPEN_WRONLY | seqan::OPEN_CREATE))
        {
            SEQAN_OMP_PRAGMA(critical (io))
            {
                std::cerr << "ERROR: Could not open " << path << " for writing!\n";
            }
            return 1;
        }
        clear(out);
        flag = true;
        return 0;
    }

    int _openUncompressed(seqan::SequenceStream & outStream, bool & flag, seqan::CharString const & path)
    {
        open(outStream, toCString(path), seqan::SequenceStream::WRITE, seqan::SequenceStream::FASTQ);
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

// ----------------------------------------------------------------------------
// Function compressData()
// ----------------------------------------------------------------------------

// Compress the data into a gzip block.  This consists of a gzip header, the deflated data, and a footer.

inline
bool compressData(seqan::CharString & outBuffer, seqan::CharString const & inBuffer)
{
    // We will try different values for outLen in case compression yields more data then there is in the output buffer.
    resize(outBuffer, 1.5 * length(inBuffer));

    unsigned const HEADER_LENGTH = 10;
    unsigned const FOOTER_LENGTH = 8;

    // The following flags are taken from the BAM file format.
    const int GZIP_WINDOW_BITS = -15; // no zlib header
    const int Z_DEFAULT_MEM_LEVEL = 8;

    // Write header.
    outBuffer[0] = '\037';  // ID1
    outBuffer[1] = '\213';  // ID2
    outBuffer[2] = 8;   // CM
    outBuffer[3] = 0;   // FLG
    outBuffer[4] = 0;   // MTIME
    outBuffer[5] = 0;   // MTIME
    outBuffer[6] = 0;   // MTIME
    outBuffer[7] = 0;   // MTIME
    outBuffer[8] = 0;   // XFL
    outBuffer[9] = '\xff';   // OS

    // Initialize the z_stream object we will use for compression.
    z_stream zStream;
    zStream.zalloc = Z_NULL;
    zStream.zfree = Z_NULL;
    zStream.opaque = Z_NULL;
    zStream.total_in = 0;
    zStream.avail_in = length(inBuffer);
    zStream.total_out = 0;
    zStream.avail_out = length(outBuffer) - HEADER_LENGTH - FOOTER_LENGTH;
    zStream.next_in = (Bytef*)&inBuffer[0];
    zStream.next_out = (Bytef*)(&outBuffer[0] + HEADER_LENGTH);

    int err = deflateInit2(&zStream, Z_DEFAULT_COMPRESSION, Z_DEFLATED, GZIP_WINDOW_BITS, Z_DEFAULT_MEM_LEVEL,
                           Z_DEFAULT_STRATEGY);
    if (err != Z_OK)
        return false;  // Error.
    while ((err = deflate(&zStream, Z_FINISH)) == Z_OK)
    {
        if ((err = deflateEnd(&zStream)) != Z_OK)
            return false;  // An error occured.
        if (length(outBuffer) > 10 * length(inBuffer))
            return false;  // Should never happen, though.
        resize(outBuffer, 1.5 * length(outBuffer));
        zStream.avail_out = length(outBuffer) - HEADER_LENGTH - FOOTER_LENGTH;
        zStream.next_out = (Bytef*)(&outBuffer[0] + HEADER_LENGTH + zStream.total_out);
    }
    if (err != Z_STREAM_END)
        return false;  // Could not finish compression.

    // Resize outBuffer to total bytes written.
    resize(outBuffer, zStream.total_out + HEADER_LENGTH + FOOTER_LENGTH);
    Bytef * ptr = (Bytef*)(&outBuffer[0] + zStream.total_out + HEADER_LENGTH);
    __int32 crc = crc32(0L, Z_NULL, 0);
    crc = crc32(crc, (Bytef*)&inBuffer[0], length(inBuffer));
    memcpy(ptr, (Bytef*)((char*)&crc), 4);
    ptr += 4;
    memcpy(ptr, (Bytef*)((char*)&zStream.total_in), 4);
    // Finish compression.
    if (deflateEnd(&zStream) != Z_OK)
        return false;

    return true;  // Success!
}

#endif  // #ifndef SANDBOX_BAMBI_APPS_BAM2FASTQ_SEQUENCE_SINK_H_
