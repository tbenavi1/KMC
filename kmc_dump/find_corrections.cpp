#include "stdafx.h"
#include <iostream>
#include "./../kmc_api/kmer_api.h"
#include "./../kmc_api/kmc_file.h"
#include "nc_utils.h"
#include <queue>
#include <set>
#include <numeric>
#include <map>

class CKmerAPI_Derived: public CKmerAPI
{
public:
	CKmerAPI_Derived(uint kmer_length = 0) : CKmerAPI(kmer_length)
	{
	}
	std::vector<CKmerAPI_Derived> CandidateKmers()
	{
		std::vector<CKmerAPI_Derived> candidates;
		std::string original_str;
		this->to_string(original_str);
		std::map<char, std::string> replacements;
		replacements['A']="CGT";
		replacements['C']="AGT";
		replacements['G']="ACT";
		replacements['T']="ACG";
		//Add kmers one SNP away
		uint i = kmer_length/2;
		for (uint j=0; j<3; j++)
		{
			std::string forward_edited_str = original_str;
			std::string reverse_edited_str;
			forward_edited_str.replace(i,1,replacements[original_str[i]].substr(j,1));
			CKmerAPI_Derived kmer_object_new(kmer_length);
			kmer_object_new.from_string(forward_edited_str);
			kmer_object_new.reverse();
			kmer_object_new.to_string(reverse_edited_str);
			if (forward_edited_str < reverse_edited_str)
			{
				kmer_object_new.reverse();
			}
			candidates.push_back(kmer_object_new);
		}
		//Add kmers two SNPs away
		//for (uint i=0; i<kmer_length; i++)
		for (uint i2=0; i2<kmer_length; i2++)
		{
			if (i2 == i) continue;
			for (uint j=0; j<3; j++)
			{
				for (uint j2=0; j2<3; j2++)
				{
					std::string forward_edited_str = original_str;
					std::string reverse_edited_str;
					forward_edited_str.replace(i,1,replacements[original_str[i]].substr(j,1));
					forward_edited_str.replace(i2,1,replacements[original_str[i2]].substr(j2,1));
					CKmerAPI_Derived kmer_object_new(kmer_length);
					kmer_object_new.from_string(forward_edited_str);
					kmer_object_new.reverse();
					kmer_object_new.to_string(reverse_edited_str);
					if (forward_edited_str < reverse_edited_str)
					{
						kmer_object_new.reverse();
					}
					candidates.push_back(kmer_object_new);
				}
			}
		}
		return candidates;
	}
//	std::vector<CKmerAPI_Derived> CandidateKmers()
//	{
//		std::vector<CKmerAPI_Derived> candidates;
//		std::string original_str;
//		this->to_string(original_str);
//		std::map<char, std::string> replacements;
//		replacements['A']="CGT";
//		replacements['C']="AGT";
//		replacements['G']="ACT";
//		replacements['T']="ACG";
//		//Add kmers one SNP away
//		for (uint i=0; i<kmer_length; i++)
//		{
//			for (uint j=0; j<3; j++)
//			{
//				std::string forward_edited_str = original_str;
//				std::string reverse_edited_str;
//				forward_edited_str.replace(i,1,replacements[original_str[i]].substr(j,1));
//				CKmerAPI_Derived kmer_object_new(kmer_length);
//				kmer_object_new.from_string(forward_edited_str);
//				kmer_object_new.reverse();
//				kmer_object_new.to_string(reverse_edited_str);
//				if (forward_edited_str < reverse_edited_str)
//				{
//					kmer_object_new.reverse();
//				}
//				candidates.push_back(kmer_object_new);
//			}
//		}
//		//Add kmers two SNPs away
//		for (uint i=0; i<kmer_length; i++)
//		{
//			for (uint i2=0; i2<i; i2++)
//			{
//				for (uint j=0; j<3; j++)
//				{
//					for (uint j2=0; j2<3; j2++)
//					{
//						std::string forward_edited_str = original_str;
//						std::string reverse_edited_str;
//						forward_edited_str.replace(i,1,replacements[original_str[i]].substr(j,1));
//						forward_edited_str.replace(i2,1,replacements[original_str[i2]].substr(j2,1));
//						CKmerAPI_Derived kmer_object_new(kmer_length);
//						kmer_object_new.from_string(forward_edited_str);
//						kmer_object_new.reverse();
//						kmer_object_new.to_string(reverse_edited_str);
//						if (forward_edited_str < reverse_edited_str)
//						{
//							kmer_object_new.reverse();
//						}
//						candidates.push_back(kmer_object_new);
//					}
//				}
//			}
//		}
//		return candidates;
//	}
};

struct family_member
{
	CKmerAPI_Derived kmer_object;
	uint64 counter;
};

std::vector<family_member> get_corrections(family_member my_family_member, CKMCFile *database)
{
	uint64 candidate_counter;
	std::vector<family_member> corrections;
	CKmerAPI_Derived kmer_object = my_family_member.kmer_object;
	//For each candidate kmer within the edit distance from the error kmer
	for (CKmerAPI_Derived kmer_candidate_object : kmer_object.CandidateKmers())
	{
		//If the candidate kmer is in the database of genomic kmers (and set its counter to candidate_counter)
		if (database->CheckKmer(kmer_candidate_object, candidate_counter))
		{
			family_member candidate_family_member = {kmer_candidate_object, candidate_counter};
			corrections.push_back(candidate_family_member);
		}
	}
	return corrections;
}

int _tmain(int argc, char* argv[])
{
	CKMCFile kmer_data_base_listing;
	CKMCFile kmer_data_base_RA;
	int32 i;
	uint32 min_count_to_set = 0;
	uint32 max_count_to_set = 0;
	std::string input_file_name;
	std::string input_file_name2;
	std::string output_file_name;

	FILE * out_file;
	//------------------------------------------------------------
	// Parse input parameters
	//------------------------------------------------------------
	if(argc < 3)
	{
		return EXIT_FAILURE;
	}

	for(i = 1; i < argc; ++i)
	{
		if(argv[i][0] == '-')
		{	
			if(strncmp(argv[i], "-ci", 3) == 0)
				min_count_to_set = atoi(&argv[i][3]);
			else if(strncmp(argv[i], "-cx", 3) == 0)
					max_count_to_set = atoi(&argv[i][3]);
		}
		else
			break;
	}

	if(argc - i < 2)
	{ 
		return EXIT_FAILURE;
	}

	input_file_name = std::string(argv[i++]);
	input_file_name2 = std::string(argv[i++]);
	output_file_name = std::string(argv[i]);

	if((out_file = fopen (output_file_name.c_str(),"wb")) == NULL)
	{
		return EXIT_FAILURE;
	}

	setvbuf(out_file, NULL ,_IOFBF, 1 << 24);

	//------------------------------------------------------------------------------
	// Open kmer database for listing and print kmers within min_count and max_count
	//------------------------------------------------------------------------------

	if (!kmer_data_base_listing.OpenForListing(input_file_name))
	{
		return EXIT_FAILURE ;
	}
	else if (!kmer_data_base_RA.OpenForRA(input_file_name2))
	{
		return EXIT_FAILURE ;
	}
	else
	{
		uint32 _kmer_length;
		uint32 _mode;
		uint32 _counter_size;
		uint32 _lut_prefix_length;
		uint32 _signature_len;
		uint32 _min_count;
		uint64 _max_count;
		uint64 _total_kmers;

		kmer_data_base_listing.Info(_kmer_length, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);

		char str[1024];
		uint32 counter_len;
		
		CKmerAPI_Derived err_object(_kmer_length);
		
		if(min_count_to_set)
		if (!(kmer_data_base_listing.SetMinCount(min_count_to_set)))
				return EXIT_FAILURE;
		if(max_count_to_set)
		if (!(kmer_data_base_listing.SetMaxCount(max_count_to_set)))
				return EXIT_FAILURE;
		if(min_count_to_set)
		if (!(kmer_data_base_RA.SetMinCount(min_count_to_set)))
				return EXIT_FAILURE;
		if(max_count_to_set)
		if (!(kmer_data_base_RA.SetMaxCount(max_count_to_set)))
				return EXIT_FAILURE;

		if (_mode) //quake compatible mode
		{
			float counter;
			while (kmer_data_base_listing.ReadNextKmer(err_object, counter))
			{
				err_object.to_string(str);
				str[_kmer_length] = '\t';
				counter_len = CNumericConversions::Double2PChar(counter, 6, (uchar*)str + _kmer_length + 1);				
				str[_kmer_length + 1 + counter_len] = '\n';
				fwrite(str, 1, _kmer_length + counter_len + 2, out_file);			
			}
		}
		else
		{
			uint64 err_counter;
			std::vector<family_member> corrections;
			while (kmer_data_base_listing.ReadNextKmer(err_object, err_counter))
			{
				family_member my_family_member = {err_object, err_counter};
				corrections = get_corrections(my_family_member, &kmer_data_base_RA);
				if (corrections.size() >= 1)
				{
					std::string err_str;
					err_object.to_string(err_str);
					CKmerAPI_Derived correction;
					std::sort(corrections.begin(), corrections.end(), [](family_member f1, family_member f2) {return f1.counter > f2.counter;});
					correction = (corrections.at(0)).kmer_object;
					std::string corr_str;
					correction.to_string(corr_str);
					std::string line;
					line.append(err_str + "\t" + corr_str + "\n");
					const char *c = line.c_str();
					fprintf(out_file, "%s", c);
				}
			}
		}
		fclose(out_file);
		kmer_data_base_listing.Close();
		kmer_data_base_RA.Close();
	}

	return EXIT_SUCCESS; 
}
