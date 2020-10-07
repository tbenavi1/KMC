#include "stdafx.h"
#include <iostream>
#include "./../kmc_api/kmer_api.h"
#include "./../kmc_api/kmc_file.h"
#include "nc_utils.h"
#include <list>
#include <vector>
#include <math.h>
#include <fstream>
#include <algorithm>
#include <exception>

int get_type(uint32_t& coverage, int& error_threshold, int& het_threshold, int& unique_threshold)
{
	if (coverage <= error_threshold)
	{
		return 0;
	}
	else if (coverage <= het_threshold)
	{
		return 1;
	}
	else if (coverage <= unique_threshold)
	{
		return 2;
	}
	else
	{
		return 3;
	}
}

std::vector<std::string> get_adjacent(CKMCFile& file, std::string& kmer, int& error_threshold, int& het_threshold, int& unique_threshold)
{
	std::vector<uint32_t> v;
	std::vector<std::string> adjacent_kmers;
	//for all possible nucleotide extensions from kmer
	for (char const &c: "ACGT") {
		std::string adjacent_kmer = kmer.substr(1)+c;
		file.GetCountersForRead(adjacent_kmer, v);
		int current_type = get_type(v[0], error_threshold, het_threshold, unique_threshold);
		//if the adjacent kmer is not an error
		if (current_type > 0) {
		  adjacent_kmers.push_back(adjacent_kmer);
		}
	}
	return adjacent_kmers;
}

std::vector<std::string> get_paths(CKMCFile& file, int& error_threshold, int& het_threshold, int& unique_threshold, std::string& before_first_error_kmer, std::string& after_last_error_kmer, int& min_distance_of_path, int& max_distance_of_path, int& max_nodes_to_search, int& k, bool& queue_broken)
{
	//This function finds paths starting from before_first_error_kmer and ending at
	//after_last_error_kmer where each kmer in the path is a nonerror kmer.
	//We follow all paths, but cut the depth of any path to max_distance_of_path.
	//We also ensure a minimum depth equal to min_distance_of_path.
	std::list<std::string> queue;
	queue.push_back(before_first_error_kmer);
	//Initialize edited_paths to store all the edited paths that are found.
	std::vector<std::string> edited_paths;
	//We use i as a counter for how many nodes have been visited in the search.
	//If we haven't finished the search within max_nodes_to_search nodes, 
	//we break the search and don't edit this error block.
	//This drastically speeds up the run time for those few regions that can't be edited.
	//Thankfully, it doesn't seem to impact effectiveness, since most error blocks
	//can be edited before this threshold.
	int i = 0;
	//This flag keeps track of whether we had to stop the search early.
	queue_broken = false;
	while(!queue.empty())
	{
		i++;
		std::string current_path = queue.front();
		std::string current_kmer = current_path.substr(current_path.length()-k);
		queue.pop_front();
		int current_depth = current_path.length()-k;
		//If we have to terminate search early
		if (i > max_nodes_to_search)
		{
			//std::cout << "queue broken" << '\n';
			queue_broken = true;
			break;
		}
		//If the depth of this node hasn't exceed the max distance of the path
		if (current_depth <= max_distance_of_path)
		{
			//Extend the path by one nucleotide, keep the ones that are not error kmers
			std::vector<std::string> adjacent_kmers = get_adjacent(file, current_kmer, error_threshold, het_threshold, unique_threshold);
			for (auto adjacent_kmer : adjacent_kmers)
			{
				std::string edited_path = current_path + adjacent_kmer.back();
				//std::cout << edited_path << '\n';
				//If we have found a path of nonerror kmers which bridges the error block
				//and doesn't terminate too early (i.e. before min_distance_of_path)
				if ((adjacent_kmer == after_last_error_kmer) && (current_depth + 1 >= min_distance_of_path))
				{
					edited_path.erase(edited_path.end()-k, edited_path.end());
					edited_paths.push_back(edited_path.substr(1));
				}
				//Else we haven't found a path yet
				else
				{
					queue.push_back(edited_path);
				}
			}
		}
	}
	return edited_paths;
}

std::string extend_left_unique(std::string& beginning_portion, int& i, int& k, int& error_threshold, int& het_threshold, int& unique_threshold, CKMCFile& file, int& read_number)
{
	std::string extension;
	//The beginning portion should be of length i+k
	std::string kmer = beginning_portion.substr(i);
	std::vector<uint32_t> v;
	int j = 0;
	while (j < i)
	{
		bool found_path = false;
		std::string new_kmer;
		for (char const &c: "ACGT")
		{
			std::string adjacent_kmer = c+kmer.substr(0, k-1);
			file.GetCountersForRead(adjacent_kmer, v);
			int current_type = get_type(v[0], error_threshold, het_threshold, unique_threshold);
			//if the adjacent kmer is not an error
			if (current_type > 0)
			{
				//we terminate at a branching path
				if (found_path)
				{
					j--;
					//remove the previous character we added
					extension = extension.substr(1);
					extension = beginning_portion.substr(0, i-j) + extension;
					//std::cout << "we terminate at a branching path" << '\n';
					return extension;
				}
				else
				{
					found_path = true;
					j++;
					new_kmer = adjacent_kmer;
					extension = c+extension;
				}
			}
		}
		//we terminate at no paths
		if (!found_path)
		{
			extension = beginning_portion.substr(0, i-j) + extension;
			//std::cout << "we terminate at no paths" << '\n';
			return extension;
		}
		kmer = new_kmer;
	}
	//we terminate at the end of the read
	//std::cout << "we terminate at the end of the read" << '\n';
	return extension;
}


std::string extend_right_unique(std::string& ending_portion, int& i, int& k, int& error_threshold, int& het_threshold, int& unique_threshold, CKMCFile& file, int& read_number)
{
	std::string extension;
	//The ending portion should be of length k+(read_length-first_error_idx-k+1)
	std::string kmer = ending_portion.substr(0, k);
	std::vector<uint32_t> v;
	int j = 0;
	while (j < i)
	{
		bool found_path = false;
		std::string new_kmer;
		for (char const &c: "ACGT")
		{
			std::string adjacent_kmer = kmer.substr(1)+c;
			file.GetCountersForRead(adjacent_kmer, v);
			int current_type = get_type(v[0], error_threshold, het_threshold, unique_threshold);
			//if the adjacent kmer is not an error
			if (current_type > 0)
			{
				//we terminate at a branching path
				if (found_path)
				{
					j--;
					//remove the previous character we added
					//extension = extension.substr(1);
					extension.pop_back();
					//extension = beginning_portion.substr(0, i-j) + extension;
					//std::cout << extension << '\n';
					//std::cout << ending_portion.substr(j+k) << '\n';
					extension = extension + ending_portion.substr(j+k);
					//std::cout << "we terminate at a branching path" << '\n';
					return extension;
				}
				else
				{
					found_path = true;
					j++;
					new_kmer = adjacent_kmer;
					//extension = c+extension;
					extension = extension+c;
				}
			}
		}
		//we terminate at no paths
		if (!found_path)
		{
			//extension = beginning_portion.substr(0, i-j) + extension;
			//std::cout << extension << '\n';
			//std::cout << ending_portion.substr(j+k) << '\n';
			extension = extension + ending_portion.substr(j+k);
			//std::cout << "we terminate at no paths" << '\n';
			return extension;
		}
		kmer = new_kmer;
	}
	//we terminate at the end of the read
	//std::cout << "we terminate at the end of the read" << '\n';
	return extension;
}

bool IsErr (uint32_t coverage)
{
	int error_threshold = 9;
	return (coverage <= error_threshold);
}

bool IsRep (uint32_t coverage)
{
	int unique_threshold = 72;
	return (coverage > unique_threshold);
}

std::string remove_err (std::vector<uint32_t>& v, std::string& read, int& read_number, CKMCFile& file, std::ofstream& erredits_output_file, std::ofstream& errpaths_output_file)
{
	//initialize variables
	std::string edited_read;
	int k = 21;
	int error_threshold = 9; //8 for dmel, 9 for atha
	int het_threshold = 30; //21 for dmel, 30 for atha
	int unique_threshold = 72; //60 for dmel, 72 for atha

	//iterate over counts to edit errors
	int previous_type = -1;
	int first_nonerror_idx;
	int last_nonerror_idx;
	int first_error_idx;
	int last_error_idx;
	std::string before_first_error_kmer;
	std::string after_last_error_kmer;
	for (int i = 0; i < v.size(); i++)
	{
		int current_type = get_type(v[i], error_threshold, het_threshold, unique_threshold);
		//if kmer is error
		if (current_type == 0)
		{
			//if this is the first kmer of the read
			if (previous_type == -1)
			{
				first_error_idx = i;
				last_nonerror_idx = i-1;
			}
			//if previous kmer was an error, we are continuing the error block
			if (previous_type == 0)
			{
				;
			}
			//if previous kmer was not an error, we are leaving the nonerror block
			if (previous_type > 0)
			{
				//get kmer that is right before the first error kmer of error block
				first_error_idx = i;
				last_nonerror_idx = i-1;
				before_first_error_kmer = read.substr(i-1, k);
				std::string nonerror_portion = read.substr(first_nonerror_idx, last_nonerror_idx - first_nonerror_idx + 1);
				edited_read += nonerror_portion;
			}
			previous_type = current_type;
		}
		//else if kmer is nonerror
		else if (current_type > 0) {
			//if this is the first kmer of the read
			if (previous_type == -1)
			{
				first_nonerror_idx = i;
				last_error_idx = i-1;
			}
			//if previous kmer was error, and we are at the beginning of the read
			if (previous_type == 0 && before_first_error_kmer.empty())
			{
				//std::cout << "We are at the beginning of read " << read_number << '\n';
				//The very beginning of the read is an error portion, we are currently not editing
				last_error_idx = i-1;
				first_nonerror_idx = i;
				after_last_error_kmer = read.substr(i, k); // added
				std::string original_error_portion = read.substr(0, i);
				std::string original_error_block = original_error_portion + after_last_error_kmer;
				erredits_output_file << ">read" << read_number << "_firsterrorkmer" << first_error_idx << "_lasterrorkmer" << last_error_idx << "_original" << '\n';
				erredits_output_file << " " << original_error_portion << " " << after_last_error_kmer << '\n';
				std::vector<uint32_t> w;
				file.GetCountersForRead(original_error_block, w);
				for (int j=0; j < w.size(); j++)
				{
					erredits_output_file << w.at(j) << " ";
				}
				erredits_output_file << '\n';
				std::string beginning_portion = read.substr(0, i+k);
				std::string edited_error_portion = extend_left_unique(beginning_portion, i, k, error_threshold, het_threshold, unique_threshold, file, read_number);
				std::string edited_error_block = edited_error_portion + after_last_error_kmer;
				erredits_output_file << ">read" << read_number << "_firsterrorkmer" << first_error_idx << "_lasterrorkmer" << last_error_idx << "_edited" << '\n';
				erredits_output_file << " " << edited_error_portion << " " << after_last_error_kmer << '\n';
				file.GetCountersForRead(edited_error_block, w);
				for (int j=0; j < w.size(); j++)
				{
					erredits_output_file << w.at(j) << " ";
				}
				erredits_output_file << '\n';
				edited_read += edited_error_portion;
			}
			//if previous kmer was error, we have left the error block
			if (previous_type == 0 && !before_first_error_kmer.empty())
			{
				int number_of_error_kmers = i - first_error_idx;
				//If the position of after_last_error_kmer overlaps before_first_error_kmer
				//we keep progressing as if nothing has happened, waiting to find another non_error kmer
				if (number_of_error_kmers < k)
				{
					//std::cout << "number of error kmers is less than k. kmer: " << i << '\n';
					current_type = previous_type;
					continue;
				}
				//get kmer that is right after the last error kmer of block
				last_error_idx = i-1;
				first_nonerror_idx = i;
				after_last_error_kmer = read.substr(i, k);
				int min_distance_of_path = k;
				int max_distance_of_path = ceil(1.2 * number_of_error_kmers);
				int max_nodes_to_search = 1000;
				bool queue_broken;
				std::vector<std::string> edited_paths = get_paths(file, error_threshold, het_threshold, unique_threshold, before_first_error_kmer, after_last_error_kmer, min_distance_of_path, max_distance_of_path, max_nodes_to_search, k, queue_broken);
				//Only edit if we fully traveled all of the paths (with depth between min and max distance)
				//and there exists only exactly one path that bridges the gap
				if (!queue_broken && edited_paths.size() == 1)
				{
					std::string original_error_portion = read.substr(first_error_idx+k-1, last_error_idx - first_error_idx + 2 - k);
					std::string original_error_block = before_first_error_kmer + original_error_portion + after_last_error_kmer;
					erredits_output_file << ">read" << read_number << "_firsterrorkmer" << first_error_idx << "_lasterrorkmer" << last_error_idx << "_original" << '\n';
					erredits_output_file << before_first_error_kmer << " " << original_error_portion << " " << after_last_error_kmer << '\n';
					std::vector<uint32_t> w;
					file.GetCountersForRead(original_error_block, w);
					for (int j=0; j < w.size(); j++)
					{
						erredits_output_file << w.at(j) << " ";
					}
					erredits_output_file << '\n';
					std::string edited_error_portion = edited_paths[0].substr(k-1);
					std::string edited_error_block = before_first_error_kmer + edited_error_portion + after_last_error_kmer;
					erredits_output_file << ">read" << read_number << "_firsterrorkmer" << first_error_idx << "_lasterrorkmer" << last_error_idx << "_edited" << '\n';
					erredits_output_file << before_first_error_kmer << " " << edited_error_portion << " " << after_last_error_kmer << '\n';
					//std::vector<uint32_t> w;
					file.GetCountersForRead(edited_error_block, w);
					for (int j=0; j < w.size(); j++)
					{
						erredits_output_file << w.at(j) << " ";
					}
					erredits_output_file << '\n';
					edited_read += edited_paths[0];
				}
				//Either we didn't fully travel all of the paths (with depth between min and max distance)
				//or there does not exist only exactly one path that bridges that gap
				else
				{
					std::string uneditable_error_portion = read.substr(first_error_idx, last_error_idx - first_error_idx + 1);
					//output_file << uneditable_error_portion;
					edited_read += uneditable_error_portion;
					//Let's keep track of all the possible paths, just for bookkeeping. 
					if ((!queue_broken) && (edited_paths.size() > 1)) //only record those paths where the search completed and there are multiple
					{
						std::string portion;
						errpaths_output_file << ">read" << read_number << "_firsterrorkmer" << first_error_idx << "_lasterrorkmer" << last_error_idx << "_original" << '\n';
						if (uneditable_error_portion.length() >= k-1)
						{
							errpaths_output_file << before_first_error_kmer << " " << uneditable_error_portion.substr(k-1) << " " << after_last_error_kmer << '\n';
							portion = before_first_error_kmer + uneditable_error_portion.substr(k-1) + after_last_error_kmer;
						}
						else
						{
							errpaths_output_file << before_first_error_kmer.front() << " " << uneditable_error_portion << " " << after_last_error_kmer << '\n';
							portion = before_first_error_kmer.front() + uneditable_error_portion + after_last_error_kmer;
						}
						std::vector<uint32_t> w;
						file.GetCountersForRead(portion, w);
						for (int j=0; j < w.size(); j++)
						{
							errpaths_output_file << w.at(j) << " ";
						}
						errpaths_output_file << '\n';
						for (int l = 0; l < edited_paths.size(); l++)
						{
							std::string edited_path = edited_paths[l];
							errpaths_output_file << ">read" << read_number << "_firsterrorkmer" << first_error_idx << "_lasterrorkmer" << last_error_idx << "_path" << l << '\n';
							if (edited_path.length() >= k-1)
							{
								errpaths_output_file << before_first_error_kmer << " " << edited_path.substr(k-1) << " " << after_last_error_kmer << '\n'; //I hope this is correct, we'll check
								portion = before_first_error_kmer + edited_path.substr(k-1) + after_last_error_kmer;
							}
							else 
							{
								errpaths_output_file << before_first_error_kmer.front() << " " << edited_path << " " << after_last_error_kmer << '\n';
								portion = before_first_error_kmer.front() + edited_path + after_last_error_kmer;
							}
							std::vector<uint32_t> w;
							file.GetCountersForRead(portion, w);
							for (int j=0; j < w.size(); j++)
							{
								errpaths_output_file << w.at(j) << " ";
							}
							errpaths_output_file << '\n';
						}
					}
				}
			}
			//if previous kmer is nonerror, we are continuing a non error block
			if (previous_type > 0)
			{
				;
			}
			previous_type = current_type;
		}
	}
	//We have reached the end of the read, let's make sure we have added the last bit of the read
	if (previous_type == 0)
	{
		//std::cout << "We are at the end of read " << read_number << '\n';
		//We have "left" the error portion of the read
		last_error_idx = v.size()-1;
		first_nonerror_idx = v.size();
		std::string original_error_portion = read.substr(first_error_idx+k-1);
		std::string original_error_block = before_first_error_kmer + original_error_portion;
		erredits_output_file << ">read" << read_number << "_firsterrorkmer" << first_error_idx << "_lasterrorkmer" << last_error_idx << "_original" << '\n';
		erredits_output_file << before_first_error_kmer << " " << original_error_portion << " " << '\n';
		std::vector<uint32_t> w;
		file.GetCountersForRead(original_error_block, w);
		for (int j=0; j < w.size(); j++)
		{
			erredits_output_file << w.at(j) << " ";
		}
		erredits_output_file << '\n';
		std::string ending_portion = read.substr(first_error_idx-1);
		int original_error_portion_length = v.size()-first_error_idx;
		std::string edited_error_portion = extend_right_unique(ending_portion, original_error_portion_length, k, error_threshold, het_threshold, unique_threshold, file, read_number);
		std::string edited_error_block = before_first_error_kmer + edited_error_portion;
		erredits_output_file << ">read" << read_number << "_firsterrorkmer" << first_error_idx << "_lasterrorkmer" << last_error_idx << "_edited" << '\n';
		erredits_output_file << before_first_error_kmer << " " << edited_error_portion << " " << '\n';
		file.GetCountersForRead(edited_error_block, w);
		for (int j=0; j < w.size(); j++)
		{
			erredits_output_file << w.at(j) << " ";
		}
		erredits_output_file << '\n';
		edited_read += edited_error_portion;
	}
	if (previous_type > 0)
	{
		//We have "left" the nonerror portion of the read
		first_error_idx = v.size();
		last_nonerror_idx = v.size()-1;
		std::string nonerror_portion = read.substr(first_nonerror_idx, last_nonerror_idx - first_nonerror_idx + k);
		//output_file << nonerror_portion;
		edited_read += nonerror_portion;
	}
	//output_file << '\n';
	edited_read += '\n';
	return edited_read;
}

//std::string remove_het ()
//{
//
//}

int main(int argc, char* argv[])
{
	//load KMC database
	CKMCFile file;
	file.OpenForRA(argv[1]);
	
	//initialize file streams
	std::string line;
	std::ifstream input_file(argv[2]);
	std::ofstream err_output_file(argv[3]);
	std::ofstream erredits_output_file(argv[4]);
	std::ofstream errpaths_output_file(argv[5]);
	std::ofstream het_output_file(argv[6]);
	std::ofstream hetedits_output_file(argv[7]);
	std::ofstream hetpaths_output_file(argv[8]);

	//load reads
	int line_num = -1;
	//currently assuming fasta format
	while (getline(input_file, line))
	{
		line_num++;
		if (line_num % 2 == 0)
		{
			//write read header to err output file and het output file
			err_output_file << line << '\n';
			het_output_file << line << '\n';
		}
		else
		{
			std::string read = line;
			int read_number = (line_num+1)/2;
			if (read_number%10000==0)
			{
				std::cout << read_number << '\n';
			}
			
			//get counters of kmers in read
			std::vector<uint32_t> v;
			file.GetCountersForRead(read, v);
			int num_kmers = v.size();
			
			//check if too many error kmers
			//if greater than 25% of kmers in the read are errors
			//discard the read to not waste time
			int num_err = std::count_if(v.begin(), v.end(), IsErr);
			double err_fraction = static_cast<double>(num_err)/num_kmers;
			if (err_fraction >= 0.25)
			{
				//std::cout << "read number " << read_number << " has " << err_fraction << " percent error kmers, discarding." << '\n';
				err_output_file << "N" << '\n'; //hopefully this is fine
				het_output_file << "N" << '\n';
				continue;
			}

			//remove errors from the read to get edited read
			std::string edited_read = remove_err(v, read, read_number, file, erredits_output_file, errpaths_output_file);

			//write edited read to err_output_file
			err_output_file << edited_read;

			//get counters of kmers in edited read
			//std::vector<uint32_t> v;
			//file.GetCountersForRead(edited_read, v);
			//int num_kmers = v.size();

			//check if too many repetitive kmers
			//if greater than 25% of kmers in the edited read are repetitive
			//discard the read to not waste time
			//int num_rep = std::count_if(v.begin(), v.end(), IsRep);
			//double rep_fraction = static_cast<double>(num_rep)/num_kmers;
			//if (rep_fraction >= 0.25)
			//{
			//	std::cout << "read number " << read_number << " has " << rep_fraction << " percent repetitive kmers, discarding." << '\n';
			//	het_output_file << "N" << '\n';
			//	continue;
			//}

			//"remove" het from the read to get smoothed read
			//std::string smoothed_read = remove_het();

			//write smoothed read to het_output_file
			//het_output_file << smoothed_read;
		}
  }
  input_file.close();
  err_output_file.close();
  erredits_output_file.close();
  errpaths_output_file.close();
  het_output_file.close();
  hetedits_output_file.close();
	hetpaths_output_file.close();
  return 0;
}
