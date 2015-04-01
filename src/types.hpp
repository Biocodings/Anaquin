#ifndef GI_TYPES_HPP
#define GI_TYPES_HPP

#include <string>

typedef double Percentage;

typedef std::string GeneID;
typedef std::string ChromoID;
typedef std::string VariantID;
typedef std::string FeatureName;

typedef std::string Sequence;

typedef long long Reads;

typedef float FPKM;
typedef float Concentration;

// Number of lines in a file (most likely a large file)
typedef long long Lines;

typedef float Fold;

typedef std::string GeneID;
typedef std::string IsoformID;

// Eg: 388488 from the first matching base
typedef long long BasePair;

#endif