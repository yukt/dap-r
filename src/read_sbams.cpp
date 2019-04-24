#include <Rcpp.h>
#include "parser.h"
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List read_sbams(const char* file, bool normalize = true)
{
  parser pars;
  ifstream infile(file);
  string line;
  
  while(getline(infile,line))
    pars.process_line(line);
  int numgroup = pars.pheno_vec.size();
  
  if(normalize)
  {
    for(int i=0;i<numgroup;i++)
    {
      pars.regress_cov(pars.pheno_vec[i], pars.covar_vec[i], pars.geno_vec[i]);
    }
  }
  
  List data;
  if(numgroup==1)
  {
    int p = pars.geno_vec[0].size();
    CharacterVector geno_names(p);
    for(int i = 0; i<p; i++)
      geno_names[i]=pars.geno_map[i];
    List geno = wrap(pars.geno_vec[0]);
    geno.attr("names") = geno_names;
    data = List::create(Named("pheno") = wrap(pars.pheno_vec[0]),
                        Named("geno")  = geno,
                        Named("controlled") = wrap(pars.covar_vec[0]),
                        Named("file") = wrap(file));
  }
  // currently do not support multiple groups
  
  data.attr("class") = "sbams";
  return data;
}


