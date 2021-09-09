#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

// This is a Rcpp function to create design matrix.
// Input: A data.frame in .maf format. All columns should be as character and EA column should be as numeric.
//        Two vectors of genes and patients' ID as character.
// Output: Design matrix as numeric
// Author: Saeid Parvandeh, November 2018

// [[Rcpp::export]]
NumericMatrix DesignMatrix(DataFrame x, CharacterVector genes, CharacterVector patients) {
  NumericVector   ea    = x["ACTION"];
  CharacterVector gene  = x["GENE"];
  CharacterVector id    = x["Tumor_Sample_Barcode"];
  NumericMatrix   dMatrix(patients.size(), genes.size());
  std::fill( dMatrix.begin(), dMatrix.end(), NumericVector::get_na() ) ;
  for (int g = 0; g < genes.size(); g++){
    Rcpp::NumericVector idx;
    for(int i = 0; i < gene.size(); i++){
      if (gene[i] == genes[g]) {
        idx.push_back(i);
      }
    }
    NumericVector   sub_ea    = ea[idx];
    CharacterVector sub_gene  = gene[idx];
    CharacterVector sub_id    = id[idx];
    for (int p = 0; p < patients.size(); p++){
      Rcpp::NumericVector EAs;
      for (int j = 0; j < sub_id.size(); j++){
        if ((sub_id[j] == patients[p]) & (sub_gene[j] == genes[g])){
          EAs.push_back(sub_ea[j]);
        }
      }
      if (Rf_length(EAs)>=1){
        dMatrix(p, g) = Rcpp::max(EAs);
      }
    }
  }
  return(dMatrix);
}