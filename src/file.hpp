#ifndef GI_FILE_HPP
#define GI_FILE_HPP

#include <vector>
#include <string>

namespace Spike
{
    struct FileInternal;
    
    class File
    {
        public:
            File(const std::string &file);
            ~File();
        
            // Returns the next line in the file
            bool nextLine(std::string &) const;
        
            // Returns the next line and parse it into tokens
            bool nextTokens(std::vector<std::string> &, const std::string &c) const;
        
        private:
            FileInternal *_imp;
    };
}

#endif