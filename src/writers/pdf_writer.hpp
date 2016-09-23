#ifndef PDF_WRITER_HPP
#define PDF_WRITER_HPP

#include "tools/script.hpp"

extern Anaquin::Scripts ReportScript();

namespace Anaquin
{
    class PDFWriter
    {
        public:

            inline void open(const FileName &output)
            {
                _output = output;
            }

            inline void addTitle(const std::string &title)
            {
                _title = title;
            }
        
            inline void addFile(const FileName &file, const std::string &title = "")
            {
                _files.push_back(file);
                _subTits.push_back(title == "" ? file : title);
            }

            inline void create(const Path &path)
            {
                assert(!_title.empty());
                
                /*
                 * Eg: <output> <title> <path> <files> <titles
                 */
                
                const auto cmd = ((boost::format("%1% %2% %3% %4% %5%")) % _output
                                                                         % _title
                                                                         % path
                                                                         % concat(_files, "")
                                                                         % concat(_subTits, "")).str();
                //Script::run(ReportScript(), "python", cmd);
            }

        private:

            std::string _title;

            // Where the report is written to
            FileName _output;
        
            // Files to be added to the report
            std::vector<FileName> _files;
        
            // Titles for each of the file
            std::vector<std::string> _subTits;
    };
}

#endif