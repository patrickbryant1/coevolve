#include <unordered_map>
#include <stdexcept> 
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cctype>
#include <math.h>
#include <cmath>
#include <array>
#include <map>

std::string split(std::string str, const char delim, int pos){
    int count = 0, pointer = 0;
    std::string token, prev;
    std::stringstream ss(str);

    while (std::getline (ss, token, delim)) {count+=1;}
    ss.str(""); ss.clear(); ss.str(str);

    std::vector< std::string > splitted;

    splitted.resize(count);

    while (std::getline (ss, token, delim)) {splitted[pointer]=token; pointer+=1;}

    if (pos >= 0){
	if (pos >= count) throw std::runtime_error("Out of splitting range!");
        token=splitted[pos]; 
    }
    else {
	if (pos < -count) throw std::runtime_error("Out of splitting range!");
        token=splitted[count+pos];
    }

    return token;
}

double MI(std::vector< std::vector < int > > parsedmsa,
          std::vector< std::array < double, 21 > > msafreqs,
	  int i, int j) {

    int x, y;
    double mi = 0, Pij;
    unsigned int seq, count = 0;
    std::array< int, 2 > couple;
    std::map < std::array< int, 2 >, double > jointfreqs;

    std::array < double, 21 > Pi = msafreqs[i], Pj = msafreqs[j];

	// for each seq in parsed msa, count the co-occurrences of aa 	
	// in positions i & j, if none of those are gaps	

    for (seq=0; seq<parsedmsa.size(); seq++) {
	if (parsedmsa[seq][i]!=21 && parsedmsa[seq][j]!=21) {
    	    couple[0]=parsedmsa[seq][i];
	    couple[1]=parsedmsa[seq][j];
	    jointfreqs[couple] ++;
        }
    }

	// normalize each count on the total number of aminoacid	
        // combinations occurrences					

//    for (auto coupleidx=jointfreqs.begin(); coupleidx!=jointfreqs.end(); coupleidx++) {
//        count += coupleidx->second;
//    }
//    for (auto coupleidx=jointfreqs.begin(); coupleidx!=jointfreqs.end(); coupleidx++) {
//        jointfreqs[coupleidx->first] = coupleidx->second/count;
//    }

	// given the joint frequencies of co-occurrence (Pij) and the marginal
        // frequencies (Pi[x] & Pj[y]) for each pair of residues, calculates the 
        // mutual information

    for (auto coupleidx=jointfreqs.begin(); coupleidx!=jointfreqs.end(); coupleidx++) {
	Pij = coupleidx->second;
        x = coupleidx->first[0];
	y = coupleidx->first[1];
        mi += Pij*log(Pij/(Pi[x]*Pj[y]));
     }
    return mi;
}

int main (int argc, char** argv){

    if (argc != 3) {
    std::cout << "2 command line arguments required (in this order): " << std::endl;
    std::cout << "	path of input MSA file in a3m format;" << std::endl;
    std::cout << "      path to save output" << std::endl;
    return 0;
    }

    std::string msapath = argv[1], outpath = argv[2];

    double mi, maxmi;
    unsigned int i = 0, j = 0, minlength = 0, msalength = 0, gapcount = 0, hitcount = 0, count = 0;
    std::string msaline, apcline, msahit = "", residue, name, id1, id2;
    std::ifstream msafile;
    std::ofstream outfile;
    std::vector< std::vector < int > > parsedmsa;

    std::vector< std::vector < double > > mi_matrix;
    std::vector< std::array < double, 21 > > msafreqs;

    std::unordered_map< char, int > aamap;

    aamap['A'] = 1, aamap['C'] = 2, aamap['D'] = 3, aamap['E'] = 4, aamap['F'] = 5;
    aamap['G'] = 6, aamap['H'] = 7, aamap['I'] = 8, aamap['K'] = 9, aamap['L'] = 10;
    aamap['M'] = 11, aamap['N'] = 12, aamap['P'] = 13, aamap['Q'] = 14, aamap['R'] = 15;
    aamap['S'] = 16, aamap['T'] = 17, aamap['V'] = 18, aamap['W'] = 19, aamap['Y'] = 20;
    aamap['U'] = 21, aamap['Z'] = 21, aamap['X'] = 21, aamap['J'] = 21, aamap['B'] = 21; 
    aamap['O'] = 21, aamap['-'] = 21;

    name = split(msapath, '/', -1);
    name = split(name, '.', 0);

    msafile.open(msapath);
    if(!msafile) throw std::runtime_error("Unable to open msa file!");
    outfile.open(outpath+name+".txt");
    if(!outfile) throw std::runtime_error("Unable to open output file!");


	// Check MSA parameters: number of hits (hitcount) and alignment length (minlength)

    while (std::getline(msafile, msaline)) {
        if (msaline.at(0) == '>') {
	    if (minlength == 0 || msalength < minlength) {minlength = msalength;}
	    msalength = 0; hitcount ++;
	}
        else {msalength += msaline.length();}
    }

	// Reset input file read and initialize datastructure 

    msafile.clear();
    msafile.seekg(0, std::ios::beg);
    parsedmsa.resize(hitcount, std::vector< int > (minlength, 0) );
    msafreqs.resize(minlength, std::array< double, 21 > {0});

	// MSA parsing, removing lowcase residues and seqs with more than 90% gaps

    j = 0;
    while (std::getline(msafile, msaline)) {
        if (msaline.at(0) != '>') {
            for (i=0; i<msaline.length(); i++) {if (msaline.at(i) == '-') {gapcount ++;}}
            msahit += msaline;
        }
        if (msaline.at(0) == '>' && msahit != "") {
            if (gapcount/minlength < 0.9) {           
                for (i=0; i<msahit.length(); i++) {
                    if (aamap.find(msahit.at(i)) != aamap.end()) {parsedmsa[j][i] = aamap.at(msahit.at(i));}
                }
                j++;
            }
            msahit = ""; gapcount = 0;
        }
        //std::cout << msaline << std::endl;
    }

	// adds the last MSA hit, which is complete reaching EOF (out of the while loop)

    if (gapcount/minlength < 0.9) {
        for (i=0; i<msahit.length(); i++) { 
            if (aamap.find(msahit.at(i)) != aamap.end()) {parsedmsa[j][i] = aamap.at(msahit.at(i));}
        }
    }
    j++; parsedmsa.resize(j);
     
	// Counting of aminoacids occurrence for each position in the MSA

    for (j=0; j<parsedmsa.size(); j++) {
        for (i=0; i<minlength; i++) {msafreqs[i][parsedmsa[j][i]]++;}
    }

	// Computing frequencies for the obtained counts

//    for (i=0; i<msafreqs.size(); i++) {
//        count = 0;
//	for (j=0; j<20; j++) {count += msafreqs[i][j];}
//	for (j=0; j<20; j++) {msafreqs[i][j] = msafreqs[i][j]/count;}
//    }


	// MI calculation

    mi_matrix.resize(minlength, std::vector < double > (minlength, 0.0));

    maxmi = 0;
    for (i=0; i<msalength; i++) {
        for (j=i+1; j<msalength; j++) {
            mi = MI(parsedmsa, msafreqs, i, j); 
            mi_matrix[i][j] = mi;
	    mi_matrix[j][i] = mi;
            //std::cout<< i << " "<< j << " " << mi << std::endl;
            if (std::abs(mi)>maxmi) {maxmi = std::abs(mi);}
        }
    }
    
    msafile.close();

	// MI matrix save 
 
    for (i=0; i<msalength; i++) {
        for (j=0; j<msalength; j++) { outfile << mi_matrix[i][j] << " ";}
        outfile << "\n";
        }
    outfile.close();

    std::cout << name << " " << maxmi << std::endl;
    std::cout << "Done!" << std::endl;

    return 0;

}
