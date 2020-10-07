#include "stdafx.h"
#include <iostream>
#include "./../kmc_api/kmer_api.h"
#include "./../kmc_api/kmc_file.h"
#include "nc_utils.h"

int main(int argc, char* argv[]) {	
	CKMCFile file;
	//file.OpenForRA("/seq/schatz/tbenavi/hifi/DMEL/KMC/DMEL"); //assume success
	file.OpenForRA(argv[1]);
	std::vector<uint32_t> v;
	//file.GetCountersForRead("TAAGTACCGTTTAGTTTTAACCACTCCCAAGCGGCGCA", v);
	//file.GetCountersForRead(argv[1],v);
	file.GetCountersForRead(argv[2],v);
	
	for (auto c : v) {
		std::cout << c << " ";
	}
	std::cout << '\n';
	return 0;
}
