#include "stdafx.h"
#include <iostream>
#include <fstream>
#include "./../kmc_api/kmer_api.h"
#include "./../kmc_api/kmc_file.h"
#include "nc_utils.h"

int main(int argc, char* argv[]) {	
	CKMCFile file;
	//file.OpenForRA("/seq/schatz/tbenavi/hifi/DMEL/KMC/DMEL"); //assume success
	file.OpenForRA(argv[1]);
	//currently assuming that there is only a single sequence
	//and that it is on only one line of the file
	std::ifstream input_file;
	input_file.open(argv[2]);
	std::string line;
	getline(input_file, line);
	getline(input_file, line);
	//remove "\n" at end of string
	//line.erase(line.length()-1);
	std::vector<uint32_t> v;
	//file.GetCountersForRead("TAAGTACCGTTTAGTTTTAACCACTCCCAAGCGGCGCA", v);
	//file.GetCountersForRead(argv[1],v);
	file.GetCountersForRead(line, v);
	
	for (auto c : v) {
		std::cout << c << "\n";
	}
	//std::cout << '\n';
	return 0;
}
