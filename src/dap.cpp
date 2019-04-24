#include <Rcpp.h>
#include "controller.h"
#include "result_parser.h"
#include <stdlib.h>
#include <stdio.h>
#include <cstring>

using namespace Rcpp;
using namespace std;

//' DAP
//'
//' @param arg a list
//' @return Nothing.
//' @export
// [[Rcpp::export]]
List dap(List arg) {
  char grid_file[128];
  char data_file[128];
  char zval_file[128];
  char est_file[128];
  char ld_file[128];
  char out_file[128];
  char log_file[128];
  char gene_name[64];

  char prior_file[128];

  memset(gene_name,0,64);
  memset(grid_file,0,128);
  memset(out_file,0,128);
  memset(log_file,0,128);
  memset(data_file,0,128);
  memset(ld_file,0,128);
  memset(zval_file,0,128);
  memset(est_file,0,128);

  memset(prior_file,0,128);



  int ld_format = 1; // for correlation matrix


  double abf_option = -1;

  int msize =-1;

  int output_all = 0;

  double pes = 1.0;
  double pi1 = -1;
  double lambda = 0.5;
  double ld_control = -1;

  int size_limit = -1;

  int sample_size = -1;
  double syy = -1;

  double snp_select_thresh = -1;
  double size_select_thresh = -1;

  // alternative non-fm running options
  int run_scan = 0;
  int extract_ss = 0;

  int thread = 1;

  vector< string > mystrings =  arg.attr("names");

  for(int i=0; i<arg.size(); i++)
  {
    // required data files and additional info
    if(mystrings[i]=="data")
    {
      strcpy(data_file, arg[i]);
      continue;
    }
    if(mystrings[i]=="grid")
    {
      strcpy(grid_file, arg[i]);
      continue;
    }
    if(mystrings[i]=="d_zval")
    {
      strcpy(zval_file, arg[i]);
      continue;
    }
    if(mystrings[i]=="data_ld")
    {
      strcpy(ld_file, arg[i]);
      continue;
    }
    if(mystrings[i]=="data_ld2")
    {
      strcpy(ld_file, arg[i]);
      ld_format = 2;
      continue;
    }
    if(mystrings[i]=="d_est")
    {
      strcpy(est_file, arg[i]);
      continue;
    }
    if(mystrings[i]=="d_n")
    {
      sample_size = atoi(arg[i]);
      continue;
    }
    if(mystrings[i]=="d_syy")
    {
      syy = atof(arg[i]);
      continue;
    }

    // prior file
    if(mystrings[i]=="prior")
    {
      strcpy(prior_file, arg[i]);
      continue;
    }

    // output file
    if(mystrings[i]=="output")
    {
      strcpy(out_file, arg[i]);
      continue;
    }
    // if(mystrings[i]=="logfile")
    // {
    //   strcpy(log_file, arg[i]);
    //   continue;
    // }

    stop("Unknow option ", mystrings[i]);
  }

  controller con;

  if(strlen(data_file)!=0){
    con.initialize(data_file,grid_file);
  }else if(strlen(zval_file)!=0 && strlen(ld_file)!=0 ){
    con.initialize(zval_file, ld_file, grid_file,sample_size, ld_format);
  }else if(strlen(ld_file)!=0 && strlen(est_file)!=0 && sample_size >0 && syy>0){
    con.initialize(est_file, ld_file, grid_file, sample_size, syy,ld_format);
  }else{
    stop("No suitable input data specified! \n");
  }

  con.set_outfile(out_file, log_file);
  con.set_gene(gene_name);
  con.set_abf_option(abf_option);
  con.set_thread(thread);
  con.set_size_limit(size_limit);

  if(ld_control>=0)
    con.set_ld_control(ld_control);

  if(msize>=1){
    con.set_max_size(msize);
  }

  if(strlen(prior_file)==0){
    if(pi1 != -1){
      if(0<pi1 && pi1 <1){
        con.set_prior(pi1);
      }else{
        fprintf(stderr, "Warning: pi1 specification is outside the range, ignored...\n");
        con.set_prior_exp(pes);
      }
    }else{
      // default
      con.set_prior_exp(pes);
    }
  }else
    con.set_prior(prior_file);

  if(output_all == 1)
    con.set_output_all();

  if(snp_select_thresh>=0)
    con.set_snp_select_thresh(snp_select_thresh);

  if(size_select_thresh >=0)
    con.set_size_select_thresh(size_select_thresh);

  con.run_option = 0;

  if(run_scan){
    con.run_option = 1;
  }

  if(extract_ss==1){
    con.run_option =2;
  }
  if(extract_ss==2){
    con.run_option =3;
  }

  con.print_dap_config();
  con.run();

  result_parser result(out_file);
  DataFrame model_summary = DataFrame::create( Named("rank") = wrap(result.model_rank),
                                               Named("size") = wrap(result.model_size),
                                               Named("posterior") = wrap(result.model_posterior),
                                               Named("score") = wrap(result.model_score),
                                               Named("configure") = wrap(result.model_configure));
  DataFrame SNP_summary = DataFrame::create( Named("rank") = wrap(result.snp_rank),
                                             Named("name") = wrap(result.snp_name),
                                             Named("PIP") = wrap(result.snp_pip),
                                             Named("score") = wrap(result.snp_score),
                                             Named("cluster") = wrap(result.snp_cluster));
  DataFrame cluster_summary = DataFrame::create( Named("cluster") = wrap(result.cluster),
                                                 Named("member SNP") = wrap(result.member),
                                                 Named("average PIP") = wrap(result.cluster_pip),
                                                 Named("average r2") = wrap(result.cluster_r2));
  List result_list = List::create(Named("models") = model_summary,
                                  Named("SNPs") = SNP_summary,
                                  Named("clusters") = cluster_summary,
                                  Named("N") = wrap(result.N),
                                  Named("size.mean") = wrap(result.model_size_mean),
                                  Named("size.sd") = wrap(result.model_size_sd),
                                  Named("logNC") = wrap(result.LogNC),
                                  Named("log10NC") = wrap(result.Log10NC),
                                  Named("PIP.min") = wrap(result.min_pip));
  result_list.attr("class") = "dap";
  return result_list;
}


