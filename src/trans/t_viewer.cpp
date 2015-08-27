#include "data/script.hpp"
#include "trans/t_viewer.hpp"

using namespace Anaquin;

void TViewer::generate(const std::string &file, const Options &options)
{
    Script::viewer("Trans ABCD scripts/accepted_hits.bam");
}