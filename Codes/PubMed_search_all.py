# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 12:07:47 2017

@author: teng-kuei
"""
from Bio import Entrez
import pandas as pd
import os
def search(query):
    Entrez.email = 'tkhsu00@gmail.com'
    handle = Entrez.esearch(db='pubmed', 
                            sort='relevance', 
                            retmax='100000',
                            retmode='xml', 
                            term=query)
    results = Entrez.read(handle)
    return results
    
def fetch_details(id_list):
    ids = ','.join(id_list)
    Entrez.email = 'tkhsu00@gmail.com'
    handle = Entrez.efetch(db='pubmed',
                           retmode='xml',
                           id=ids)
    results = Entrez.read(handle)
    return results
    
def chunks(l, n):
    n = max(1, n)
    return [l[i:i+n] for i in range(0, len(l), n)]    
    
if __name__ == '__main__':
    df = pd.read_csv("path_to_the_directory",sep=",")    
    #df = pd.read_csv("/cedar/hsu/EA_gene_Feb2017/COAD.tsv",sep="\t")    
    genes = df["Gene"].tolist()
    for count, gene in enumerate(genes[::-1]):
        fn = "path_to_the_directory"+gene+".tsv"
        if not os.path.isfile(fn):
            with open(fn,"w") as f:
                f.write("Gene\tNo.literature_w_Cancer\tNo.literature\n")
            
                print (count+1,gene)
                
                num_paper, num_paper_ca = 0,0
                
                results_ca = search("\""+gene+"\" AND (\"gene\" or \"protein\") AND \"CANCER\"")
                id_list_ca = results_ca['IdList']
                if len(id_list_ca)>0:
                    papers_ca = fetch_details(id_list_ca)
                    num_paper_ca = len(papers_ca['PubmedArticle'])            
                
                results = search("\""+gene+"\" AND (\"gene\" or \"protein\")")
                id_list = results['IdList']
                if len(id_list)>0:
                    papers = fetch_details(id_list)
                    num_paper = len(papers['PubmedArticle'])
                
                f.write("{}\t{}\t{}\n".format(gene,num_paper_ca,num_paper))
