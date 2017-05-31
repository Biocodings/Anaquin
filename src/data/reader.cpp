#include <fstream>
#include <sstream>
#include <iostream>
#include <assert.h>
#include "data/reader.hpp"
#include <boost/algorithm/string.hpp>

using namespace Anaquin;

struct Anaquin::ReaderInternal
{
    Line line;
    
    // Defined only for file input
    std::string file;
    
    // Implementation for file
    std::shared_ptr<std::ifstream> f;

    // Implementation for memory
    std::shared_ptr<std::stringstream> s;
};

Reader::Reader(const Reader &r)
{
    _imp = new ReaderInternal();
    _imp->line = r._imp->line;
    _imp->f = r._imp->f;
    _imp->s = r._imp->s;

    // Make sure we start off from the default state
    reset();
}

Reader::Reader(const std::string &file, DataMode mode)
{
    if (file.empty())
    {
        throw std::runtime_error("Empty file name");
    }
    
    _imp = new ReaderInternal();
    _imp->file = file;

    if (mode == DataMode::File)
    {
        const auto f = std::shared_ptr<std::ifstream>(new std::ifstream(file));

        if (!f->good())
        {
            throw InvalidFileError(file);
        }
        else if (f->peek() == std::ifstream::traits_type::eof())
        {
            throw InvalidFileError(file);
        }

        _imp->f = f;
    }
    else
    {
        if (file.empty())
        {
            throw InvalidFileError(file);
        }

        _imp->s = std::shared_ptr<std::stringstream>(new std::stringstream(file));
    }
}

Reader::~Reader()
{
    delete _imp;
}

Line Reader::lastLine() const
{
    return _imp->line;
}

void Reader::reset()
{
    if (_imp->f)
    {
        _imp->f->clear();
        _imp->f->seekg(0, std::ios::beg);
    }
    
    if (_imp->s)
    {
        _imp->s->clear();
        _imp->s->seekg(0, std::ios::beg);
    }
}

std::string Reader::src() const
{
    return _imp->file;
}

bool Reader::nextLine(std::string &line) const
{
    retry:

    if ((_imp->f && std::getline(*_imp->f, line)) || (_imp->s && std::getline(*_imp->s, line)))
    {
        if (line.empty())
        {
            goto retry;
        }

        boost::trim(line);
        return true;
    }
    else
    {
        return false;
    }
}

bool Reader::nextTokens(std::vector<std::string> &toks, const std::string &c) const
{
    _imp->line.clear();

    if (nextLine(_imp->line))
    {
        toks.clear();
        boost::split(toks, _imp->line, boost::is_any_of(c));
        return true;
    }
    else
    {
        return false;
    }
}
