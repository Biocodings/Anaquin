#ifndef AS_TYPES_HPP
#define AS_TYPES_HPP

#include <string>

typedef std::string GeneID;
typedef std::string ChromoID;
typedef std::string FeatureName;

typedef long long Reads;

typedef float FPKM;
typedef float Expression;

// Number of lines in a file (most likely a large file)
typedef long long Lines;

typedef float Fold;

typedef std::string GeneID;
typedef std::string IsoformID;

typedef float Percentage;

// Eg: 388488 from the first matching base
typedef long long BasePair;

#endif