#include <Rcpp.h>
#include "data/reader.hpp"
#include "parsers/parser_csv.hpp"

using namespace Rcpp;
using namespace Anaquin;

// [[Rcpp::export]]
List TestTKExpress()
{
  std::string s;
  
  ParserCSV::parse(Reader("/Users/tedwong/Sources/QA/results.csv"), [&](const ParserCSV::Data &d, const ParserProgress &p)
  {
    s = d[0];
  }, ",");
  

    return List::create(Rcpp::Named("Value") = s);
}