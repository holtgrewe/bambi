# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox/bambi/apps/bam2fastq/CMakeLists.txt
#
# CMakeLists.txt file for bam2fastq.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_sandbox_bambi_apps_bam2fastq)
message (STATUS "Configuring sandbox/bambi/apps/bam2fastq")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES ZLIB OpenMP)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your application.
add_executable (bam2fastq
                bam2fastq.cpp
                bam_scanner_cache.h
                converter_thread.h
                job_queue.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (bam2fastq ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (bam2fastq)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install bam2fastq in ${PREFIX}/bin directory
install (TARGETS bam2fastq
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/bam2fastq for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (bam2fastq)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:bam2fastq")
  set (CPACK_PACKAGE_NAME "bam2fastq")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "bam2fastq")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")

  seqan_configure_cpack_app (bam2fastq "bam2fastq")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:bam2fastq")

