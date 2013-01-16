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
// Parallel conversion of BAM files (with an existing BAI index) to FASTQ.
// ==========================================================================

// TODO(holtgrew): One assumption here is that one thread can saturate the disk bandwidth when writing uncompressed file.

#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include <seqan/parallel.h>
#include <seqan/sequence.h>

#include "job_queue.h"
#include "converter_thread.h"
#include "sequence_sink.h"

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class Options
// --------------------------------------------------------------------------

// This struct stores the options from the command line.

struct Options
{
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // The number of threads to use.
    int numThreads;
    // The length of a tile.
    int tileLength;
    // The maximal expected template length.  Everything that has a greater tlen will be handled with the "pile"
    // handling, i.e. integrated in a second step.
    int maxTemplateLength;
    // Chunk length for buffering input records.
    unsigned chunkLength;

    // Path to input BAM file.
    seqan::CharString inputPath;
    // Path to output FASTQ file for left/right paired and singleton reads.  outputPathLeft is used in case
    // interleavedOut is true.
    seqan::CharString outputPathLeft, outputPathRight, outputPathSingle;

    // Whether or not to write interleavedly to outputPathLeft.
    bool interleavedOut;
    // Whether or not to gzip the output.
    bool gzipOutput;
    // The string to use for joining "1", "2", and "S" to the input in the output filename.
    seqan::CharString outSep;

    Options() :
            verbosity(1), numThreads(1), tileLength(0), maxTemplateLength(0),
            interleavedOut(false), gzipOutput(false)
    {}
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function yesNo()
// --------------------------------------------------------------------------

char const * yesNo(bool b)
{
    return b ? "YES" : "NO";
}

// --------------------------------------------------------------------------
// Function printOptions()
// --------------------------------------------------------------------------

void printOptions(std::ostream & out, Options const & options)
{
    out << "__OPTIONS____________________________________________________________________\n"
        << "\n"
        << "VERBOSITY                  \t" << options.verbosity << "\n"
        << "\n"
        << "NUM THREADS                \t" << options.numThreads << "\n"
        << "\n"
        << "TILE LENGTH                \t" << options.tileLength << "\n"
        << "MAX TPL LENGTH             \t" << options.maxTemplateLength << "\n"
        << "CHUNK LENGTH               \t" << options.chunkLength << "\n"
        << "\n"
        << "WRITE INTERLEAVED          \t" << yesNo(options.interleavedOut) << "\n"
        << "GZIP OUTPUT                \t" << yesNo(options.gzipOutput) << "\n"
        << "FILE NAME SEGMENT SEPARATOR\t\"" << options.outSep << "\"\n"
        << "\n"
        << "INPUT PATH                 \t" << options.inputPath << "\n"
        << "OUTPUT PATH LEFT           \t" << options.outputPathLeft << "\n"
        << "OUTPUT PATH RIGHT          \t" << options.outputPathRight << "\n"
        << "OUTPUT PATH SINGLE         \t" << options.outputPathSingle << "\n";
}

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(Options & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("bam2fastq");
    // Set short description, version, and date.
    setShortDescription(parser, "Convert BAM to FASTQ");
    setVersion(parser, "0.1");
    setDate(parser, "January 2013");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fB-i\\fP \\fIIN.bam\\fP");
    addDescription(parser,
                   "Convert the input BAM file into a FASTQ file.  The BAI index for this file. "
                   "must already exist.");
    addDescription(parser,
                   "In the first step, each thread scans over its part of the genome and writes out "
                   "all singletons and mate pairs that have a template length below a threshold. "
                   "In a second step, all other mate pairs are processed.  See section Parameter "
                   "Overview below.  See the section Output Files for details on the "
                   "output file name creation.");

    addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

    // Tuning Options

    addSection(parser, "Tuning Options");
    addOption(parser, seqan::ArgParseOption("nt", "num-threads", "Number of threads to use.",
                                            seqan::ArgParseOption::INTEGER, "THREADS"));
    setMinValue(parser, "num-threads", "1");

    addOption(parser, seqan::ArgParseOption("", "chunk-length", "Number of records to read before writing out.",
                                            seqan::ArgParseOption::INTEGER, "NUM"));
    setMinValue(parser, "chunk-length", "1");
    setDefaultValue(parser, "chunk-length", "320000");

    addOption(parser, seqan::ArgParseOption("", "tile-length", "Tile length.",
                                            seqan::ArgParseOption::INTEGER, "LEN"));
    setMinValue(parser, "tile-length", "10000");
    setDefaultValue(parser, "tile-length", "1000000");

    addOption(parser, seqan::ArgParseOption("l", "max-template-length", "Maximal expected template length.",
                                            seqan::ArgParseOption::INTEGER, "LEN"));
    setMinValue(parser, "max-template-length", "100");
    setDefaultValue(parser, "max-template-length", "1000");

    // Input / Output
    addSection(parser, "Input / Output");
    addOption(parser, seqan::ArgParseOption("i", "in-file",
                                            "Input BAM File.  A corresponding BAI file must already exist.",
                                            seqan::ArgParseOption::INPUTFILE, "IN"));
    setRequired(parser, "in-file");

    addOption(parser, seqan::ArgParseOption("o", "out-file",
                                            "Output FASTQ file base name (use OUT to get OUT_1.fastq, OUT_2.fastq, "
                                            "and OUT_S.fastq).  Defaults to IN if input is IN.bam.",
                                            seqan::ArgParseOption::STRING, "OUT"));

    addOption(parser, seqan::ArgParseOption("oL", "out-file-left", "Output path for left paired reads.",
                                            seqan::ArgParseOption::OUTPUTFILE, "OUT"));

    addOption(parser, seqan::ArgParseOption("oR", "out-file-right", "Output path for right paired reads.",
                                            seqan::ArgParseOption::OUTPUTFILE, "OUT"));

    addOption(parser, seqan::ArgParseOption("oS", "out-file-single", "Output path for single reads.",
                                            seqan::ArgParseOption::OUTPUTFILE, "OUT"));

    addOption(parser, seqan::ArgParseOption("z", "gzip-output", "Compress the output files using gzip."));

    addOption(parser, seqan::ArgParseOption("s", "out-name-sep",
                                            "Separator to use in output filename.  If the output file name is "
                                            "automatically determined, it is the \"\\fB_\\fP\" in the output "
                                            "file names NAME_1.fastq, NAME_2.fastq, NAME_S.fastq",
                                            seqan::ArgParseOption::STRING, "SEP"));
    setDefaultValue(parser, "out-name-sep", "_");

    addOption(parser, seqan::ArgParseOption("", "interleaved",
                                            "Write output interleaved. Only one file for paired reads is generated "
                                            "in this case."));

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBbam2fastq\\fP \\fB-t\\fP 8 \\fB-i\\fP IN.bam \\fB-o\\fP OUT.fq",
                "Convert \\fIIN.bam\\fP to \\fIOUT.fq\\fP using 8 threads.");

    // Add Parameter Overview Section.
    addTextSection(parser, "Parameter Overview");
    addText(parser,
            "The program will split the genome is tiles of the same length (\\fB--tile-length\\fP). "
            "Each tile is processed atomically by one thread.  Mate pairs that have a TLEN less than "
            "\\fB--max-template-length\\fP are directly processed in the first step.  Set this "
            "parameter such that most of the pairs in your data are below this range.");
    addText(parser,
            "A greater value for \\fB--max-template-length\\fP allows more pairs to be directly processed "
            "but also increases the memory requirement as more reads in a pair with large insert sizes "
            "are kept in memory.");
    addText(parser,
            "The value for \\fB--tile-length\\fP determines the granularity.  It has to be greater than "
            "\\fB--max-template-length\\fP and is probably best left in he mega base order of magnitude.");

    // Add Output File Section.
    addTextSection(parser, "Output Files");
    addText(parser,
            "The program will write out the paired reads as (possibly split into left and right reads) "
            "and the singleton reads separately.  By default, given the input file \\fIFILE.bam\\fP, "
            "the output paths \\fIFILE_1.fastq\\fP, \\fIFILE_2.fastq\\fP, \\fIFILE_S.fastq\\fP are "
            "used.  The paired and singleton files are only created when there are any such reads.");
    addText(parser,
            "The parameters \\fB--out-name-sep\\fP, \\fB--interleaved\\fP and \\fB-z\\fP/\\fB--gzip-out\\fB "
            "also influence the automatically generated out file names.  You can enforce the file names by "
            "setting \\fB-oL\\fP/\\fB--out-file-left\\fP, \\fB-oR\\fP/\\fB--out-file-right\\fP, and "
            "\\fB-oS\\fP/\\fB--out-file-single\\fP.");
    

    // Add Author Section.
    addTextSection(parser, "Author");
    addText(parser, "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;

    getOptionValue(options.numThreads, parser, "num-threads");

    getOptionValue(options.chunkLength, parser, "chunk-length");
    getOptionValue(options.tileLength, parser, "tile-length");
    getOptionValue(options.maxTemplateLength, parser, "max-template-length");

    options.interleavedOut = isSet(parser, "interleaved");
    options.gzipOutput = isSet(parser, "gzip-output");
    getOptionValue(options.outSep, parser, "out-name-sep");

    getOptionValue(options.inputPath, parser, "in-file");

    // Get out file basename from arguments or from --in-file.
    seqan::CharString outFileBasename;
    getOptionValue(outFileBasename, parser, "out-file");
    if (empty(outFileBasename))
    {
        outFileBasename = options.inputPath;
        // Find rightmost '.'.
        for (unsigned i = 0; i < length(outFileBasename); ++i)
        {
            if (outFileBasename[length(outFileBasename) - i - 1] == '.')
            {
                resize(outFileBasename, length(outFileBasename) - i - 1);
                break;
            }
        }
    }

    // Generate out file left/right/single paths if not set.
    if (isSet(parser, "out-file-left"))
    {
        getOptionValue(options.outputPathLeft, parser, "out-file-left");
    }
    else
    {
        options.outputPathLeft = outFileBasename;
        if (!options.interleavedOut)
        {
            append(options.outputPathLeft, options.outSep);
            append(options.outputPathLeft, "1.fastq");
        }
        if (options.gzipOutput)
            append(options.outputPathLeft, ".gz");
    }
    if (isSet(parser, "out-file-right"))
    {
        getOptionValue(options.outputPathRight, parser, "out-file-right");
    }
    else
    {
        options.outputPathRight = outFileBasename;
        if (!options.interleavedOut)
        {
            append(options.outputPathRight, options.outSep);
            append(options.outputPathRight, "2.fastq");
        }
        if (options.gzipOutput)
            append(options.outputPathRight, ".gz");
    }
    if (isSet(parser, "out-file-single"))
    {
        getOptionValue(options.outputPathSingle, parser, "out-file-single");
    }
    else
    {
        options.outputPathSingle = outFileBasename;
        if (!options.interleavedOut)
        {
            append(options.outputPathSingle, options.outSep);
            append(options.outputPathSingle, "S.fastq");
        }
        if (options.gzipOutput)
            append(options.outputPathSingle, ".gz");
    }

    return seqan::ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function processPiles()
// --------------------------------------------------------------------------

// Process pile members of threads.
int processPiles(SequenceSink & sink, seqan::String<ConverterThread> & threads, Options const & options)
{
    // Build mapping from left read name to (thread id, thread-local seq id).
    typedef std::map<seqan::CharString, std::pair<unsigned, unsigned> > TMap;
    TMap leftIdMap;

    for (unsigned i = 0; i < length(threads); ++i)
        for (unsigned j = 0; j < length(threads[i].pile.leftIds); ++j)
            leftIdMap[threads[i].pile.leftIds[j]] = std::make_pair(i, j);

    // Process right reads, chunk-wise for using with sequencSink.
    ResultBuffer buffer;
    for (unsigned i = 0; i < length(threads); ++i)
    {
        ResultBuffer const & pile = threads[i].pile;
        for (unsigned j = 0; j < length(pile.rightIds); ++j)
        {
            TMap::const_iterator it = leftIdMap.find(pile.rightIds[j]);
            if (it == leftIdMap.end())
                continue;  // Skip if not in left pile.
            int seqId = it->second.second;
            ResultBuffer const & otherPile = threads[it->second.first].pile;
            buffer.insert(
                    otherPile.leftIds[seqId], otherPile.leftSeqs[seqId], otherPile.leftQuals[seqId],
                    pile.rightIds[j], pile.rightSeqs[j], pile.rightQuals[j]);

        }
        if (buffer.writeToSink(sink) != 0)
            return 1;
        buffer.clear();
    }

    return 0;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    double startTime = 0, overallStartTime = sysTime();

    // Parse the command line.
    seqan::ArgumentParser parser;
    Options options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Print program header.
    std::cerr << "BAM 2 FASTQ\n"
              << "===========\n\n";

    // Print the command line arguments back to the user.
    printOptions(std::cerr, options);

    // Precomputation step.
    std::cerr << "\n"
              << "__PRECOMPUTATION_____________________________________________________________\n"
              << "\n";
    // Check that we can open the BAM file for reading.
    std::cerr << "Opening BAM File " << options.inputPath << " ..." << std::flush;
    seqan::BamStream bamStream(toCString(options.inputPath));
    if (!isGood(bamStream))
    {
        std::cerr << " ERROR\n"
                  << "\n"
                  << "Could not open file " << options.inputPath << " for reading.\n";
        return 1;
    }
    std::cerr << " OK\n";

    // Open BAI file for reading.
    seqan::CharString baiFilename = options.inputPath;
    append(baiFilename, ".bai");
    std::cerr << "Loading BAI Index " << baiFilename << " ..." << std::flush;
    seqan::BamIndex<seqan::Bai> baiIndex;
    if (read(baiIndex, toCString(baiFilename)) != 0)
    {
        std::cerr << " ERROR\n"
                  << "\n"
                  << "Could not read index from file " << baiFilename << "\n";
        return 1;
    }
    std::cerr << " OK\n";
    // Tile the genome.
    std::cerr << "Computing splitters ..." << std::flush;
    seqan::StringSet<seqan::String<int> > splitters;
    resize(splitters, length(bamStream.header.sequenceInfos));
    for (unsigned seqId = 0; seqId < length(bamStream.header.sequenceInfos); ++seqId)
    {
        int chromLength = bamStream.header.sequenceInfos[seqId].i2;
        int numTiles = chromLength / options.tileLength + 1;
        seqan::computeSplitters(splitters[seqId], chromLength, numTiles);
    }
    JobQueue jobQueue;
    for (unsigned seqId = 0; seqId < length(splitters); ++seqId)
    {
        for (unsigned i = 0; i + 1 < length(splitters[seqId]); ++i)
            jobQueue.push(ConverterJob(seqId, splitters[seqId][i], splitters[seqId][i + 1]));
    }
    std::cerr << " OK\n";
    std::cerr << "  Partitioned genome into " << jobQueue.size() << " tiles.\n";

    // Perform parallel conversion.
    std::cerr << "\n"
              << "__CONVERSION_________________________________________________________________\n"
              << "\n";

    // Building conversion options.
    ConversionOptions convOptions;
    convOptions.inputPath = options.inputPath;
    convOptions.maxTemplateLength = options.maxTemplateLength;
    convOptions.chunkLength = options.chunkLength;

    // The SequenceSink allows thread-safe writing of blocks of FASTQ data.
    SequenceSink sink(options.interleavedOut, options.gzipOutput, toCString(options.outputPathLeft),
                      toCString(options.outputPathRight), toCString(options.outputPathSingle));

    std::cerr << "Creating threads ..." << std::flush;
    seqan::String<ConverterThread> threads;
    for (int i = 0; i < options.numThreads; ++i)
    {
        ConverterThread thread(sink, i, convOptions, jobQueue);
        appendValue(threads, thread);
    }
    std::cerr << " OK\n";

    startTime = sysTime();
    std::cerr << "Concordant Conversion ... " << std::flush;

    // Fire up threads.  The master does the orphan conversion first and then grabs jobs as all other threads.
    SEQAN_OMP_PRAGMA(parallel num_threads(options.numThreads))
    {
        // Convert options by master.
        SEQAN_OMP_PRAGMA(master)
        {
            threads[omp_get_thread_num()].convertOrphans();
            SEQAN_OMP_PRAGMA(critical (io))
            {
                std::cerr << "[orphans done]";
            }
        }

        bool stop = false;
        while (!stop)
        {
            ConverterJob job;
            if (!(stop = !jobQueue.pop(job)))
            {
                threads[omp_get_thread_num()].convertMapped(job);
            }
            SEQAN_OMP_PRAGMA(critical (io))
            {
                if (threads[omp_get_thread_num()].dot())
                    std::cerr << ".";
                if (options.verbosity >= 2)
                    std::cerr << "[thread " << omp_get_thread_num() << ", job done, rId=" << job.rId << ", beginPos=" << job.beginPos << ", endPos=" << job.endPos << "]";
            }
        }
    }

    std::cerr << " OK\n";
    std::cerr << "  Took " << (sysTime() - startTime) << " s\n";

    startTime = sysTime();
    std::cerr << "Remaining Conversion ... " << std::flush;

    if (processPiles(sink, threads, options) != 0)
        return 1;
    // // Process piles.
    // //
    // // Read in all left pile sequences.
    // typedef seqan::StringSet<seqan::CharString> TCharStringSet;
    // TCharStringSet ids, seqs, quals;
    // for (unsigned i = 0; i < length(threads); ++i)
    // {
    //     seqan::CharString id, seq, qual;
    //     seqan::SequenceStream & leftPileStream = *threads[i]._leftPileSeqStream;
    //     while (!atEnd(leftPileStream))
    //     {
    //         if (readRecord(id, seq, qual, leftPileStream) != 0)
    //         {
    //             std::cerr << "\nERROR: Could not write record from left pile of " << i << "\n";
    //             return 1;
    //         }
    //         appendValue(ids, id);
    //         appendValue(seqs, seq);
    //         appendValue(quals, qual);
    //     }
    // }

    // // Build cache over the pile sequence names.
    // seqan::NameStoreCache<TCharStringSet> idsCache(ids);

    // // Read in all right pile sequences.
    // for (unsigned i = 0; i < length(threads); ++i)
    // {
    //     seqan::CharString id, seq, qual;
    //     seqan::SequenceStream & rightPileStream = *threads[i]._rightPileSeqStream;
    //     while (!atEnd(rightPileStream))
    //     {
    //         if (readRecord(id, seq, qual, rightPileStream) != 0)
    //         {
    //             std::cerr << "\nERROR: Could not write record from left pile of " << i << "\n";
    //             return 1;
    //         }
    //         unsigned idx = 0;
    //         if (!getIdByName(ids, id, idx, idsCache))
    //         {
    //             std::cerr << "ERROR: Could not find read " << id << " in left pile.\n";
    //             return 1;
    //         }
    //         if (writeRecord(peOut, ids[idx], seqs[idx], quals[idx]) != 0)
    //         {
    //             std::cerr << "\nERROR: Could not write to PE output file\n";
    //             return 1;
    //         }
    //         if (writeRecord(peOut, id, seq, qual) != 0)
    //         {
    //             std::cerr << "\nERROR: Could not write to PE output file\n";
    //             return 1;
    //         }
    //     }
    // }
    // std::cerr << " OK\n";
    // std::cerr << "  Took " << (sysTime() - startTime) << " s\n";

    std::cerr << "\nDone converting BAM to FASTQ\n";
    // Sum up statistics on orphans.
    int numOrphans = 0, numSingletonOrphans = 0, numPairedOrphans = 0;
    int numMapped = 0, numSingletonMapped = 0, numPairedMapped = 0;
    for (unsigned i = 0; i < length(threads); ++i)
    {
        numOrphans += threads[i]._stats.numOrphans;
        numSingletonOrphans += threads[i]._stats.numSingletonOrphans;
        numPairedOrphans += threads[i]._stats.numPairedOrphans;
        numMapped += threads[i]._stats.numMapped;
        numSingletonMapped += threads[i]._stats.numSingletonMapped;
        numPairedMapped += threads[i]._stats.numPairedMapped;
    }
    std::cerr << "  Converted orphans:      " << numOrphans << "\n"
              << "              singletons: " << numSingletonOrphans << "\n"
              << "              paired:     " << numPairedOrphans << "\n"
              << "            mapped reads: " << numMapped << "\n"
              << "              singletons: " << numSingletonMapped << "\n"
              << "              paired:     " << numPairedMapped << "\n"
              << "\n"
              << "Took " << sysTime() - overallStartTime << " s\n";

    return 0;
}
