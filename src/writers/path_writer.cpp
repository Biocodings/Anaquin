#include "writers/path_writer.hpp"
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <boost/format.hpp>

using namespace Spike;

void PathWriter::open(const std::string &file)
{
    if (!path.empty())
    {
        system((boost::format("mkdir -p %1%") % path).str().c_str());
    }

    _o = std::shared_ptr<std::ofstream>(new std::ofstream(!path.empty() ? path + "/" + file : file));

    if (!_o->good())
    {
        throw std::runtime_error("Failed to open: " + file);
    }
}

void PathWriter::write(const std::string &line)
{
    *(_o) << line << std::endl;
}