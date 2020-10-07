#include <utility>                   // for std::pair
#include "stdafx.h"
#include <iostream>
#include "./../kmc_api/kmer_api.h"
#include "./../kmc_api/kmc_file.h"
#include "nc_utils.h"
#include <list>
#include <unordered_set>
#include <map>
#include <vector>
#include <math.h>
#include <fstream>
#include <numeric>

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
			std::cout << "queue broken" << '\n';
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

int main(int argc, char* argv[])
{
	//load KMC database
	CKMCFile file;
	file.OpenForRA(argv[1]);
	
	//initialize file streams
	std::string line;
	std::ifstream input_file(argv[2]);
	std::ofstream output_file(argv[3]);
	std::ofstream paths_output_file(argv[4]);
	//std::ofstream edits_output_file(argv[5]);
	std::ofstream editlengths_output_file(argv[5]);

	//load reads
	int line_num = -1;
	//currently assuming fasta format
	while (getline(input_file, line))
	{
		line_num++;
		std::cout << line_num << '\n';
		if (line_num % 2 == 0)
		{
			//write read header to output file
			output_file << line << '\n';
		}
		else
		{
			std::string read = line;
				
			//get counters of kmers in read
			std::vector<uint32_t> v;
			file.GetCountersForRead(read, v);
			//check whether the read has too many errors to edit
			int num_error = std::count(v.begin(), v.end(), 1) + std::count(v.begin(), v.end(), 2) + std::count(v.begin(), v.end(), 3) + std::count(v.begin(), v.end(), 4) + std::count(v.begin(), v.end(), 5) + std::count(v.begin(), v.end(), 6) + std::count(v.begin(), v.end(), 7) + std::count(v.begin(), v.end(), 8);
			int num_kmers = v.size();
			double error_fraction = static_cast<double>(num_error)/num_kmers;
			//if greater than 25% of kmers in the read are errors, discard the read to not waste time
			if (error_fraction >= 0.25)
			{
				int read_number = (line_num+1)/2;
				std::cout << "read number " << read_number << " has " << error_fraction << " percent error kmers, discarding." << '\n';
				output_file << "N" << '\n'; //hopefully this is fine
				continue;
			}
			
			//initialize variables
			int k = 21;
			int error_threshold = 8;
			int het_threshold = 21;
			int unique_threshold = 60;

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
				//std::cout << i << '\n';
				int current_type = get_type(v[i], error_threshold, het_threshold, unique_threshold);
				//if kmer is error
				//if (current_type == 0)
				//{
				//	current_type=previous_type;
				//	//don't edit previous_type_idx
				//}
				if ((current_type == 0) || (current_type == 1) || (current_type == 3))
				{
					//if this is the first kmer of the read
					if (previous_type == -1)
					{
						first_error_idx = i;
						last_nonerror_idx = i-1;
						//previous_type = current_type;
						//previous_type_idx = i;
					}
					//if previous kmer was an error, we are continuing the error block
					if ((previous_type == 0) || (previous_type == 1) || (previous_type == 3))
					{
						;
					}
					//if previous kmer was not an error, we are leaving the nonerror block
					//if (previous_type > 0)
					if (previous_type == 2) //changed from !=1 to > 1 so doesn't overlap == -1, actually
					{
						//get kmer that is right before the first error kmer of error block
						first_error_idx = i;
						last_nonerror_idx = i-1;
						before_first_error_kmer = read.substr(i-1, k);
						std::string nonerror_portion = read.substr(first_nonerror_idx, last_nonerror_idx - first_nonerror_idx + 1);
						output_file << nonerror_portion;
					}
					previous_type = current_type;
				}
				//else if kmer is nonerror
				//else if (current_type > 0)
				//else if (current_type != 1)
				if (current_type == 2)
				{
					//if this is the first kmer of the read
					if (previous_type == -1)
					{
						first_nonerror_idx = i;
						last_error_idx = i-1;
						//previous_type = current_type;
					}
					//if previous kmer was error, and we are at the beginning of the read
					//if (previous_type == 0 && before_first_error_kmer.empty())
					if (((previous_type == 0) || (previous_type == 1) || (previous_type == 3)) && before_first_error_kmer.empty())
					{
						//The very beginning of the read is an error portion, we are currently not editing
						last_error_idx = i-1;
						first_nonerror_idx = i;
						std::string uneditable_error_portion = read.substr(first_error_idx, last_error_idx - first_error_idx + 1);
						output_file << uneditable_error_portion;
						//previous_type = current_type;
					}
					//if previous kmer was error, we have left the error block
					//if (previous_type == 0 && !before_first_error_kmer.empty())
					if (((previous_type == 0) || (previous_type == 1) || (previous_type == 3)) && !before_first_error_kmer.empty())
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
						//and there exists only exactly one that bridges the gap
						if (queue_broken || edited_paths.size() != 2) //changed from 1 to 2
						{
							std::string uneditable_error_portion = read.substr(first_error_idx, last_error_idx - first_error_idx + 1);
							output_file << uneditable_error_portion;
							//Let's keep track of all the possible paths, just for bookkeeping. 
							if ((!queue_broken) && (edited_paths.size() > 1)) //only record those paths where the search completed and there are multiple
							{
								std::string portion;
								paths_output_file << ">read" << (line_num+1)/2 << "_firsterrorkmer" << first_error_idx << "_original" << '\n';
								if (uneditable_error_portion.length() >= k-1)
								{
									paths_output_file << before_first_error_kmer << " " << uneditable_error_portion.substr(k-1) << " " << after_last_error_kmer << '\n';
									portion = before_first_error_kmer + uneditable_error_portion.substr(k-1) + after_last_error_kmer;
								}
								else
								{
									paths_output_file << before_first_error_kmer.front() << " " << uneditable_error_portion << " " << after_last_error_kmer << '\n';
									portion = before_first_error_kmer.front() + uneditable_error_portion + after_last_error_kmer;
								}
								std::vector<uint32_t> w;
								file.GetCountersForRead(portion, w);
								for (int j=0; j < w.size(); j++)
								{
									paths_output_file << w.at(j) << " ";
								}
								paths_output_file << '\n';
								for (int l = 0; l < edited_paths.size(); l++)
                {
                  std::string edited_path = edited_paths[l];
                  paths_output_file << ">read" << (line_num+1)/2 << "_firsterrorkmer" << first_error_idx << "_path" << l << '\n';
									if (edited_path.length() >= k-1)
									{
                  	paths_output_file << before_first_error_kmer << " " << edited_path.substr(k-1) << " " << after_last_error_kmer << '\n'; //I hope this is correct, we'll check
										portion = before_first_error_kmer + edited_path.substr(k-1) + after_last_error_kmer;
									}
									else 
									{
										paths_output_file << before_first_error_kmer.front() << " " << edited_path << " " << after_last_error_kmer << '\n';
										portion = before_first_error_kmer.front() + edited_path + after_last_error_kmer;
									}
									std::vector<uint32_t> w;
									file.GetCountersForRead(portion, w);
									for (int j=0; j < w.size(); j++)
									{
										paths_output_file << w.at(j) << " ";
									}
									paths_output_file << '\n';
                }
							}
						} else {
							//std::string uneditable_error_portion = read.substr(first_error_idx, last_error_idx - first_error_idx + 1);
							//std::string edited_path = edited_paths[0];
							//edits_output_file << ">read" << (line_num+1)/2 << "_firsterrorkmer" << first_error_idx << "_original" << '\n';
							//edits_output_file << before_first_error_kmer << " " << uneditable_error_portion.substr(k-1) << " " << after_last_error_kmer << '\n';
							//edits_output_file << ">read" << (line_num+1)/2 << "_firsterrorkmer" << first_error_idx << "_path0" << '\n';
							//edits_output_file << before_first_error_kmer << " " << edited_path.substr(k-1) << " " << after_last_error_kmer << '\n';
							int original_length = read.substr(first_error_idx, last_error_idx - first_error_idx + 1).substr(k-1).length();
							//std::string edited_path = edited_paths[0].substr(k-1);
							int edited_length = edited_paths[0].substr(k-1).length();
							int length_difference = edited_length - original_length;
							editlengths_output_file << length_difference << '\n';
							std::vector<uint32_t>w0;
							std::vector<uint32_t>w1;
							file.GetCountersForRead(edited_paths[0], w0);
							file.GetCountersForRead(edited_paths[1], w1);
							float average0 = std::accumulate(w0.begin(), w0.end(), 0.0)/w0.size();
							float average1 = std::accumulate(w1.begin(), w1.end(), 0.0)/w1.size();
							std::cout << "Average0: " << average0 << '\n';
							std::cout << "Average1: " << average1 << '\n';
							if (average0 > average1)
							{
								output_file << edited_paths[0];
							}
							else
							{
								output_file << edited_paths[1];
							}
						}
					//previous_type = current_type;
					}
					//if previous kmer is nonerror, we are continuing a non error block
					//if (previous_type > 0)
					if (previous_type == 2)
					{
						;
					}
					previous_type = current_type;
				}
			}
			//We have reached the end of the read, let's make sure we have added the last bit of the read
			//if (previous_type == 0)
			//if (previous_type == 1)
			if ((previous_type == 0) || (previous_type == 1) || (previous_type == 3))
			{
				//We have "left" the error portion of the read
				last_error_idx = v.size()-1;
				first_nonerror_idx = v.size();
				std::string uneditable_error_portion = read.substr(first_error_idx, last_error_idx - first_error_idx + k);
				output_file << uneditable_error_portion;
			}
			//if (previous_type > 0)
			if (previous_type == 2)
			{
				//We have "left" the nonerror portion of the read
				first_error_idx = v.size();
				last_nonerror_idx = v.size()-1;
				std::string nonerror_portion = read.substr(first_nonerror_idx, last_nonerror_idx - first_nonerror_idx + k);
				output_file << nonerror_portion;
			}
			output_file << '\n';
    }
  }
	input_file.close();
	output_file.close();
  return 0;
}
