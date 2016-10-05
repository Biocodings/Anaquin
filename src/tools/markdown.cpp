#include <sstream>
#include <iostream>
#include <boost/format.hpp>
#include "tools/markdown.hpp"

using namespace Anaquin;

Scripts MarkDown::RItem::generate() const
{
    std::stringstream ss;
    
    ss << "\n## ";
    ss << _title;
    ss << "\n\n";
    ss << "```{r results='hide', message=FALSE, warning=FALSE, echo=FALSE, fig.align='center', fig.height=49}\n";
    ss << _txt;
    ss << "\n```\n\n" << std::endl;

    return ss.str();
}

Scripts MarkDown::TextItem::generate() const
{
    std::stringstream ss;
    
    ss << "\n## ";
    ss << _title;
    ss << "\n\n";
    ss << "```{ eval=FALSE}\n";
    ss << _txt;
    ss << "\n```\n\n" << std::endl;
    
    return ss.str();
}

Scripts MarkDown::Section::generate() const
{
    std::stringstream ss;
    
    ss << "\n# " << _title << std::endl;
    
    for (auto &item : _items)
    {
        ss << item->generate() << "\\pagebreak" << std::endl;
    }
    
    return ss.str();
}

Scripts MarkDown::generate(const std::string &title) const
{
    std::stringstream ss;
    
    const auto header = "---\n"
                        "title: 'Anaquin: %s'\n"
                        "header-includes: \\usepackage{graphicx}\n"
                        "output:\n"
                        "    pdf_document:\n"
                        "        keep_tex: true\n"
                        "        toc: yes\n"
                        "        toc_depth: 2\n"
                        "---\n\n\\pagebreak";

    ss << (boost::format(header) % title).str() << std::endl;
    
    for (auto &sect : _sections)
    {
        ss << sect.generate() << std::endl;
        ss << "\\pagebreak" << std::endl;
    }
    
    return ss.str();
}
