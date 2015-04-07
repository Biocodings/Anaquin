#include <fstream>
#include "writers/writer.hpp"

using namespace Spike;

struct Spike::WriterImpl
{
    std::shared_ptr<std::ofstream> o;
};

Writer::Writer(const std::string &file)
{
    _impl = std::shared_ptr<Spike::WriterImpl>(new Spike::WriterImpl());
    _impl->o = std::shared_ptr<std::ofstream>(new std::ofstream(file));
    
    if (!_impl->o->good())
    {
        throw std::runtime_error("ds.sd.ds");
    }
}

void Writer::write(const std::string &line)
{
    *(_impl->o) << line << std::endl;
}