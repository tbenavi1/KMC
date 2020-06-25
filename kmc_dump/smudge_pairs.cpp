#include "stdafx.h"
#include <iostream>
#include "./../kmc_api/kmer_api.h"
#include "./../kmc_api/kmc_file.h"
#include "nc_utils.h"
#include <queue>
#include <set>
#include <numeric>


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
		std::string nucleotides = "ACGT";
		i = kmer_length/2
		for (uint j=0; j<nucleotides.length(); j++)
		{
			if (original_str[i] != nucleotides[j]) //changed < to !=
			{
				std::string forward_edited_str = original_str;
				std::string reverse_edited_str;
				forward_edited_str.replace(i,1,nucleotides.substr(j,1));
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
		return candidates;
	}
};

struct family_member
{
	CKmerAPI_Derived kmer_object;
	uint64 counter;
};

bool operator <(const family_member& x, const family_member& y)
{
	return x.kmer_object < y.kmer_object;
}


std::vector<family_member> get_family(family_member my_family_member, CKMCFile *database)
{
	uint64 candidate_counter;
	std::vector<family_member> family;
	std::set<CKmerAPI_Derived> family_objects;
	std::queue<CKmerAPI_Derived> Q;
	CKmerAPI_Derived kmer_object = my_family_member.kmer_object;
	family_objects.insert(kmer_object);
	family.push_back(my_family_member);
	Q.push(kmer_object);
	std::string my_family_member_kmer;
	my_family_member.kmer_object.to_string(my_family_member_kmer);
	while (!Q.empty())
	{
		kmer_object = Q.front();
		Q.pop();
		//For each candidate kmer
		for (CKmerAPI_Derived kmer_candidate_object : kmer_object.CandidateKmers())
		{
			std::string kmer_candidate;
			kmer_candidate_object.to_string(kmer_candidate);
			//If the candidate kmer is later in the alphabet than the original family member and it is not yet in family and it is in the database
			//Actually, the check for candidate kmer later in the alphabet than original is superfluous. If we correctly track visited members we will never run into this case
			//(kmer_candidate.compare(my_family_member_kmer) > 0) && 
			if ((family_objects.find(kmer_candidate_object) == family_objects.end()) && database->CheckKmer(kmer_candidate_object, candidate_counter))
			{
				//add candidate kmer to the family, and to the queue
				family_objects.insert(kmer_candidate_object);
				family_member candidate_family_member = {kmer_candidate_object, candidate_counter};
				family.push_back(candidate_family_member);
				Q.push(kmer_candidate_object);
			}
		}
	}
	return family;
}

int _tmain(int argc, char* argv[])
{
	CKMCFile kmer_data_base_listing;
	CKMCFile kmer_data_base_RA;
	int32 i;
	uint32 min_count_to_set = 0;
	uint32 max_count_to_set = 0;
	std::string input_file_name;
	std::string output_file_name;
	std::string output_file_name2;

	FILE * out_file;
	FILE * out_file2;
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
	output_file_name = std::string(argv[i++]);
	output_file_name2 = std::string(argv[i]);

	if((out_file = fopen (output_file_name.c_str(),"wb")) == NULL)
	{
		return EXIT_FAILURE;
	}

	if((out_file2 = fopen (output_file_name2.c_str(),"wb")) == NULL)
	{
		return EXIT_FAILURE;
	}


	setvbuf(out_file, NULL ,_IOFBF, 1 << 24);
	setvbuf(out_file2, NULL ,_IOFBF, 1 << 24);

	//------------------------------------------------------------------------------
	// Open kmer database for listing and print kmers within min_count and max_count
	//------------------------------------------------------------------------------

	if (!kmer_data_base_listing.OpenForListing(input_file_name))
	{
		return EXIT_FAILURE ;
	}
	else if (!kmer_data_base_RA.OpenForRA(input_file_name))
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
		
		CKmerAPI_Derived kmer_object(_kmer_length);
		
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
			while (kmer_data_base_listing.ReadNextKmer(kmer_object, counter))
			{
				kmer_object.to_string(str);
				str[_kmer_length] = '\t';
				counter_len = CNumericConversions::Double2PChar(counter, 6, (uchar*)str + _kmer_length + 1);				
				str[_kmer_length + 1 + counter_len] = '\n';
				fwrite(str, 1, _kmer_length + counter_len + 2, out_file);			
			}
		}
		else
		{
			uint64 counter;
			std::set<family_member> visited_members;
			std::vector<family_member> family;
			while (kmer_data_base_listing.ReadNextKmer(kmer_object, counter))
			{
				//std::string kmer_test;
				//kmer_object.to_string(kmer_test);
				//std::cout << kmer_test << "\n";
				family_member my_family_member = {kmer_object, counter};
				std::set<family_member>::iterator it = visited_members.find(my_family_member);
				//If this kmer has not been visited
				if (it == visited_members.end())
				{
					//get the family (connected component) of this kmer
					family = get_family(my_family_member, &kmer_data_base_RA);
					std::string kmer_str;
					kmer_object.to_string(kmer_str);
					std::cout << family.size() << "\n";
					if (family.size() >= 2) // changed from > to >=
					{
						//add family members to visited, except the original which we shouldn't see again
						visited_members.insert(family.begin()+1, family.end()); //added +1
					}
					//if ((family.size() == 2) && (20 <= (family.at(0)).counter) && ((family.at(0)).counter <= 50) && (20 <= (family.at(1)).counter) && ((family.at(1)).counter <= 50))
					if ((family.size() == 2))
					{
						std::sort(family.begin(), family.end(), [](family_member f1, family_member f2) {return f1.counter < f2.counter;});
						std::string line1;
						std::string line2;
						line1.append(std::to_string((family.at(0)).counter) + "\t" + std::to_string((family.at(1)).counter));
						std::string kmer_str1;
						std::string kmer_str2;
						(family.at(0)).kmer_object.to_string(kmer_str1);
						(family.at(1)).kmer_object.to_string(kmer_str2);
						line2.append(kmer_str1 + "\t" + kmer_str2);
						line1.append("\n");
						line2.append("\n");
						const char *c = line1.c_str();
						const char *c2 = line2.c_str();
						fprintf(out_file, "%s", c);
						fprintf(out_file2, "%s", c2);
					}
				}
				//if the kmer has been visited
				else
				{
					visited_members.erase(it);
				}
			}
		}
		fclose(out_file);
		kmer_data_base_listing.Close();
		kmer_data_base_RA.Close();
	}

	return EXIT_SUCCESS; 
}
