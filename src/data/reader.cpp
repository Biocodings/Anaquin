#include <fstream>
#include <sstream>
#include <iostream>
#include <assert.h>
#include "data/reader.hpp"
#include <boost/algorithm/string.hpp>

using namespace Spike;

struct Spike::ReaderInternal
{
    std::string t;
    std::shared_ptr<std::ifstream> f;
    std::shared_ptr<std::stringstream> s;
};

Reader::Reader(const std::string &file, DataMode mode)
{
    _imp = new ReaderInternal();

    if (mode == DataMode::File)
    {
        const auto f = std::shared_ptr<std::ifstream>(new std::ifstream(file));

        if (!f->good())
        {
            throw InvalidFileError(file);
        }
        else if (f->peek() == std::ifstream::traits_type::eof())
        {
            throw EmptyFileError(file);
        }

        _imp->f = f;
    }
    else
    {
        if (file.empty())
        {
            throw EmptyFileError(file);
        }

        _imp->s = std::shared_ptr<std::stringstream>(new std::stringstream(file));
    }
}

Reader::~Reader()
{
    delete _imp;
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

        return true;
    }
    else
    {
        return false;
    }
}

bool Reader::nextTokens(std::vector<std::string> &tokens, const std::string &c) const
{
    _imp->t.clear();

    if (nextLine(_imp->t))
    {
        tokens.clear();
        boost::split(tokens, _imp->t, boost::is_any_of(c));
        return true;
    }
    else
    {
        return false;
    }
}