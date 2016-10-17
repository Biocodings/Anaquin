#ifndef R_WRITER_HPP
#define R_WRITER_HPP

#include <set>
#include <math.h>
#include <numeric>
#include "stats/linear.hpp"
#include <boost/format.hpp>

// Defined in main.cpp
extern std::string date();

// Defined in main.cpp
extern std::string __full_command__;

namespace Anaquin
{
    class MappingStats;
    
    struct RWriter
    {
        static Scripts createLogistic(const FileName    &,
                                      const std::string &,
                                      const std::string &,
                                      const std::string &,
                                      const std::string &,
                                      const std::string &,
                                      bool showLOQ);
        
        static Scripts createMultiLinear(const FileName    &,
                                         const Path        &,
                                         const std::string &,
                                         const std::string &,
                                         const std::string &,
                                         const std::string &,
                                         const std::string &,
                                         const std::string &,
                                         bool showLOQ,
                                         bool shouldLog,
                                         const std::string &extra = "");

        static Scripts createFold(const FileName    &,
                                  const Path        &,
                                  const std::string &,
                                  const std::string &,
                                  const std::string &,
                                  const std::string &,
                                  const std::string &,
                                  bool shouldLog,
                                  const std::string &extra = "");

        static Scripts createRLinear(const FileName    &,
                                     const Path        &,
                                     const std::string &,
                                     const std::string &,
                                     const std::string &,
                                     const std::string &,
                                     const std::string &,
                                     const std::string &,
                                     bool showLOQ);

        static Scripts createScript(const FileName &, const Scripts &);
        static Scripts createScript(const FileName &, const Scripts &, const std::string &);
    };
}

#endif
