#ifndef AS_PARSER_C_TRACKING_HPP
#define AS_PARSER_C_TRACKING_HPP

#include "Types.hpp"
#include <functional>

enum CTrackingStatus
{
    OK,
    HIData
};

struct CTracking
{
    GeneID geneID;

    FPKM fpkm;
    FPKM lFPKM;
    FPKM uFPKM;

    CTrackingStatus status;
};

struct ParserCTracking
{
    static bool parse(const std::string &file, std::function<void (const CTracking &)>);
};

#endif