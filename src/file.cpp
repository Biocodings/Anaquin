#include <assert.h>
#include <fstream>
#include <iostream>
#include "file.hpp"
#include <boost/algorithm/string.hpp>

struct FileInternal
{
    std::string t;
    std::shared_ptr<std::ifstream> f;
};

File::File(const std::string &file)
{
    const auto f = std::shared_ptr<std::ifstream>(new std::ifstream(file));
    
    if (!f->good())
    {
        std::cerr << "Failed to load " << file << std::endl;
        throw std::runtime_error("Failed to load " + file);
    }
   
    _imp = new FileInternal();
    _imp->f = f;
}

File::~File()
{
    delete _imp;
}

bool File::nextLine(std::string &line) const
{
    if (std::getline(*_imp->f, line))
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool File::nextTokens(std::vector<std::string> &tokens, const std::string &c) const
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