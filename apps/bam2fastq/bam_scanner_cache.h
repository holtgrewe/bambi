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
// A cache to be used for storing BAM records until finding the right mate
// when scanning over a SAM or BAM file.
// ==========================================================================

#ifndef SANDBOX_BAMBI_APPS_BAM2FASTQ_BAM_SCANNER_CACHE_H_
#define SANDBOX_BAMBI_APPS_BAM2FASTQ_BAM_SCANNER_CACHE_H_

#include <set>

#include <seqan/bam_io.h>

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class BamScannerCache
// ----------------------------------------------------------------------------

// Data structure for storing first-see mates.

class BamScannerCache
{
public:
    // The Key is a pair of (genomic pos, name) where genomic pos is a pair of (rId, pos).
    typedef std::pair<std::pair<int, int>, std::string> TKey;
    // A mapping from the key type to the BamAlignmentRecord at this position.
    typedef std::map<TKey, seqan::BamAlignmentRecord> TMap;

    TMap _map;

    BamScannerCache()
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

// Returns iterator for the cache.

namespace seqan {

template <>
struct Iterator<BamScannerCache>
{
    typedef typename BamScannerCache::TMap::const_iterator Type;
};

template <>
struct Iterator<BamScannerCache const>
{
    typedef BamScannerCache const TCache_;
    typedef typename TCache_::TMap::const_iterator Type;
};

}  // namespace seqan

// ============================================================================
// Functions
// ============================================================================

void dumpCache(BamScannerCache const & cache)
{
    for (auto it = cache._map.begin(); it != cache._map.end(); ++it)
        std::cerr << it->first.first.first << "\t" << it->first.first.second << "\t" << it->first.second << "\n";
}

// ----------------------------------------------------------------------------
// Function insertRecord()
// ----------------------------------------------------------------------------

void insertRecord(BamScannerCache & cache, seqan::BamAlignmentRecord const & record)
{
    BamScannerCache::TKey key(std::make_pair(record.rId, record.pos), toCString(record.qName));
    cache._map[key] = record;
}

// ----------------------------------------------------------------------------
// Function findMate()
// ----------------------------------------------------------------------------

seqan::Iterator<BamScannerCache>::Type
findMate(BamScannerCache & cache, seqan::BamAlignmentRecord const & record)
{
    BamScannerCache::TKey key(std::make_pair(record.rNextId, record.pNext), toCString(record.qName));
    return cache._map.find(key);
}

// ----------------------------------------------------------------------------
// Function end()
// ----------------------------------------------------------------------------

seqan::Iterator<BamScannerCache>::Type
end(BamScannerCache & cache, seqan::Standard const & /*tag*/)
{
    return cache._map.end();
}

// ----------------------------------------------------------------------------
// Function erase()
// ----------------------------------------------------------------------------

void erase(BamScannerCache & cache, seqan::Iterator<BamScannerCache const>::Type it)
{
    cache._map.erase(it);
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

bool empty(BamScannerCache const & cache)
{
    return cache._map.empty();
}

// ----------------------------------------------------------------------------
// Function containsMate()
// ----------------------------------------------------------------------------

bool containsMate(BamScannerCache const & cache, seqan::BamAlignmentRecord const & record)
{
    BamScannerCache::TKey key(std::make_pair(record.rNextId, record.pNext), toCString(record.qName));
    return (bool)cache._map.count(key);
}

#endif  // #ifndef SANDBOX_BAMBI_APPS_BAM2FASTQ_BAM_SCANNER_CACHE_H_
