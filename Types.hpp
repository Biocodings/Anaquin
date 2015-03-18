#ifndef AS_TYPES_HPP
#define AS_TYPES_HPP

#include <string>
#include <vector>
#include <sstream>

typedef std::string ChromoID;
typedef std::string GeneID;

typedef long long Reads;

typedef float FPKM;

// Number of lines in a file (most likely a large file)
typedef long long Lines;

// The amount added to a sample
typedef long Amount;

typedef std::string GeneID;
typedef std::string TranscriptID;

typedef float Percentage;

// Eg: 388488 from the first matching base
typedef long long BasePair;

#endif