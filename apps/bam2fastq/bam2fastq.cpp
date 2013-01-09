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

    // Path to input BAM file.
    seqan::CharString inputPath;
    // Path to output FASTQ file.
    seqan::CharString outputPath;

    Options() :
            verbosity(1), numThreads(1), tileLength(0), maxTemplateLength(0)
    {}
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function printOptions()
// --------------------------------------------------------------------------

void printOptions(std::ostream & out, Options const & options)
{
    out << "__OPTIONS____________________________________________________________________\n"
        << "\n"
        << "VERBOSITY\t" << options.verbosity << "\n"
        << "\n"
        << "NUM THREADS\t" << options.numThreads << "\n"
        << "\n"
        << "TILE LENGTH\t" << options.tileLength << "\n"
        << "MAX TPL LENGTH\t" << options.maxTemplateLength << "\n"
        << "\n"
        << "INPUT PATH\t" << options.inputPath << "\n"
        << "OUTPUT PATH\t" << options.outputPath << "\n";
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
    setVersion(parser, "0.0");
    setDate(parser, "January 2013");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fB-i\\fP \\fIIN.bam\\fP \\fB-o\\fP \\fIOUT.fq\\fP");
    addDescription(parser,
                   "Convert the input BAM file into a FASTQ file.  The BAI index for this file. "
                   "must already exist.");
    addDescription(parser,
                   "In the first step, each thread scans over its part of the genome and writes out "
                   "all singletons and mate pairs that have a template length below a threshold. "
                   "In a second step, all other mate pairs are processed.  See section Parameter "
                   "Overview below.");

    addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

    // Tuning Options

    addSection(parser, "Tuning Options");
    addOption(parser, seqan::ArgParseOption("nt", "num-threads", "Number of threads to use.",
                                            seqan::ArgParseOption::INTEGER, "THREADS"));
    setMinValue(parser, "num-threads", "1");

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
                                            "Output FASTQ file.",
                                            seqan::ArgParseOption::OUTPUTFILE, "OUT"));
    setRequired(parser, "out-file");

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

    getOptionValue(options.tileLength, parser, "tile-length");
    getOptionValue(options.maxTemplateLength, parser, "max-template-length");

    getOptionValue(options.inputPath, parser, "in-file");
    getOptionValue(options.outputPath, parser, "out-file");

    return seqan::ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
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

    std::cerr << "Creating threads ..." << std::flush;
    seqan::String<ConverterThread> threads;
    for (int i = 0; i < options.numThreads; ++i)
    {
        ConverterThread thread(i, convOptions, jobQueue);
        appendValue(threads, thread);
    }
    std::cerr << " OK\n";

    std::cerr << "Conversion ..." << std::flush;

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
            }
        }
    }

    std::cerr << " OK\n";

    std::cerr << "Joining temporary FASTQ files ..." << std::flush;
    // Open output files.
    unsigned dotPos = length(options.outputPath);  // Rightmost dot position.
    for (unsigned i = length(options.outputPath) - 1; i > 0; --i)
        if (options.outputPath[i] == '.')
        {
            dotPos = i;
            break;
        }
    seqan::CharString pePath = options.outputPath;
    seqan::CharString singletonPath = options.outputPath;
    insert(singletonPath, dotPos, ".singletons");
    seqan::SequenceStream peOut(toCString(pePath), seqan::SequenceStream::WRITE);
    if (!isGood(peOut))
    {
        std::cerr << "\nERROR: Could not open " << pePath << " for writing!\n";
        return 1;
    }
    seqan::SequenceStream singletonOut(toCString(singletonPath), seqan::SequenceStream::WRITE);
    if (!isGood(singletonOut))
    {
        std::cerr << "\nERROR: Could not open " << singletonPath << " for writing!\n";
        return 1;
    }
    // Write out data.    
    for (unsigned i = 0; i < length(threads); ++i)
    {
        if (threads[i].writeResult(peOut, singletonOut) != 0)
            return 1;
    }
    std::cerr << " OK\n";

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
              << "              paired:     " << numPairedMapped << "\n";

    return 0;
}
