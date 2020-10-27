#include "./../kmc_api/kmc_file.h"
#include "./../kmc_api/kmer_api.h"
#include "nc_utils.h"
#include "stdafx.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <list>
#include <math.h>
#include <numeric>
#include <stdlib.h>
#include <vector>

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

std::vector<std::string> get_adjacent(CKMCFile& file, std::string& kmer, int& error_threshold, int& het_threshold, int& unique_threshold, bool& going_right)
{
	std::vector<uint32_t> v;
	std::vector<std::string> adjacent_kmers;
	//for all possible nucleotide extensions from kmer
	for (char const &c: "ACGT")
	{
		std::string adjacent_kmer;
		if (going_right)
		{
			adjacent_kmer = kmer.substr(1)+c;
		}
		else
		{
			adjacent_kmer = c+kmer.substr(0, kmer.length()-1);
		}
		file.GetCountersForRead(adjacent_kmer, v);
		int current_type = get_type(v[0], error_threshold, het_threshold, unique_threshold);
		//if the adjacent kmer is not an error
		if (current_type > 0)
		{
		  adjacent_kmers.push_back(adjacent_kmer);
		}
	}
	return adjacent_kmers;
}

bool is_left_anchor(std::string& previous_kmer, int& previous_count, int& k, CKMCFile& file, int& error_threshold, int& het_threshold, int& unique_threshold)
{
	//if ((previous_kmer.length() != k) || (previous_count > unique_threshold))
	if (previous_kmer.length() != k)
	{
		return false;
	}
	bool going_right = true;
	std::vector<std::string> adjacent_kmers = get_adjacent(file, previous_kmer, error_threshold, het_threshold, unique_threshold, going_right);
	if (adjacent_kmers.size() == 2)
	{
		std::vector<uint32_t> v1;
		file.GetCountersForRead(adjacent_kmers[0], v1);
		int count1 = v1[0];
		std::vector<uint32_t> v2;
		file.GetCountersForRead(adjacent_kmers[1], v2);
		int count2 = v2[0];
		//If the sum of the coverages of the two branches is within 3 of the homozygous portion
		if ((previous_count - 3 <= count1 + count2) && (count1 + count2 <= previous_count + 3))
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	else
	{
		return false;
	}
}

bool is_right_anchor(std::string& current_kmer, int& current_count, int& k, CKMCFile& file, int& error_threshold, int& het_threshold, int& unique_threshold)
{
	//if ((current_kmer.length() != k) || (current_count > unique_threshold))
	if (current_kmer.length() != k)
	{
		return false;
	}
	bool going_right = false;
	std::vector<std::string> adjacent_kmers = get_adjacent(file, current_kmer, error_threshold, het_threshold, unique_threshold, going_right);
	if (adjacent_kmers.size() == 2)
	{
		std::vector<uint32_t> v1;
		file.GetCountersForRead(adjacent_kmers[0], v1);
		int count1 = v1[0];
		std::vector<uint32_t> v2;
		file.GetCountersForRead(adjacent_kmers[1], v2);
		int count2 = v2[0];
		//If the sum of the coverages of the two branches is within 3 of the homozygous portion
		if ((current_count - 3 <= count1 + count2) && (count1 + count2 <= current_count + 3))
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	else
	{
		return false;
	}
}

int get_type_het(int& previous_type, std::string& previous_kmer, std::string& current_kmer, int& previous_count, int& current_count, int& k, CKMCFile& file, int& error_threshold, int& het_threshold, int& unique_threshold, std::string& anchor_found, std::string& previous_anchor_found)
{
	//If previously at the beginning of the read
	if (previous_type == -1)
	{
		//if previous kmer is a left anchor
		if ((current_count < (previous_count - error_threshold)) && is_left_anchor(previous_kmer, previous_count, k, file, error_threshold, het_threshold, unique_threshold))
		{
			//we have left the hom region
			std::cout << "previous was left anchor" << '\n';
			//if (anchor_found != "none")
			//{
			//	previous_anchor_found = anchor_found;
			//}
			anchor_found == "left";
			return 1;
		}
		//if current kmer is a right anchor
		if ((current_count > (previous_count + error_threshold)) && is_right_anchor(current_kmer, current_count, k, file, error_threshold, het_threshold, unique_threshold))
		{
			//we have left the nonhom region
			std::cout << "current is right anchor" <<'\n';
			//if (anchor_found != "none")
			//{
			//previous_anchor_found = anchor_found;
			//}
			anchor_found = "right";
			return 2;
		}
		//else we continue the beginning of the read
		//return -1;
		//To keep it consistent with previous get_type, we need type to be -1 only before we start the read
		//So, in this case we will do what was done before, using the coverage and coverage thresholds to
		//determine whether the kmer is homozygous or not
		if ((het_threshold < current_count) && (current_count <= unique_threshold))
		{
			if (anchor_found != "none")
			{
			previous_anchor_found = anchor_found;
			}
			anchor_found = "none"; // is this right?
			return 2;
		}
		else
		{
			if (anchor_found != "none")
			{
			previous_anchor_found = anchor_found;
			}
			anchor_found = "none"; // is this right?
			return 1;
		}
	}
	//If previously in a nonhomozygous region
	if (previous_type == 1)
	{
		//If current kmer is a right anchor
		if ((current_count > (previous_count + error_threshold)) && is_right_anchor(current_kmer, current_count, k, file, error_threshold, het_threshold, unique_threshold))
		{
			//we have left the nonhom region
			std::cout << "current is right anchor" << '\n';
			//if (anchor_found != "none")
			//{
			//previous_anchor_found = anchor_found;
			//}
			anchor_found = "right";
			return 2;
		}
		//If previous kmer is a left anchor
		//if ((current_count < (previous_count - error_threshold)) && is_left_anchor(previous_kmer, previous_count, k, file, error_threshold, het_threshold, unique_threshold))
		//If current kmer is a left anchor
		else if (is_left_anchor(current_kmer, current_count, k, file, error_threshold, het_threshold, unique_threshold))
		{
			//this is a weird case where we find two left anchors in a row before a right anchor
			std::cout << "we found two left anchors in a row" << '\n';
			//if (anchor_found != "none")
			//{
			//previous_anchor_found = anchor_found;
			//}
			anchor_found = "left";
			return 2;
		}
		//If current kmer is a right anchor
		//if ((current_count > (previous_count + error_threshold)) && is_right_anchor(current_kmer, current_count, k, file, error_threshold, het_threshold, unique_threshold))
		//{
		//	//we have left the nonhom region
		//	std::cout << "current is right anchor" << '\n';
		//	anchor_found = "right";
		//	return 2;
		//}
		//else current kmer is not a right anchor
		else
		{
			//we continue the nonhom region
			if (anchor_found != "none")
			{
			previous_anchor_found = anchor_found;
			}
			anchor_found = "none";
			return 1;
		}
	}
	//If previously in a homozygous region
	if (previous_type == 2)
	{
		//If current kmer is a right anchor
		if ((current_count > (previous_count + error_threshold)) && is_right_anchor(current_kmer, current_count, k, file, error_threshold, het_threshold, unique_threshold))
		{
			//this is a weird case where we find two right anchors in a row before a left anchor
			std::cout << "we found two right anchors in a row" << '\n';
		}
		//If current kmer is a left anchor
		if (is_left_anchor(current_kmer, current_count, k, file, error_threshold, het_threshold, unique_threshold))
		{
			std::cout << "current is left anchor" << '\n';
			//if (anchor_found != "none")
			//{
			//previous_anchor_found = anchor_found;
			//}
			anchor_found = "left";
			return 2; //do I need this here?
		}
		//If previous kmer is a left anchor
		else if ((current_count < (previous_count - error_threshold)) && is_left_anchor(previous_kmer, previous_count, k, file, error_threshold, het_threshold, unique_threshold)) // do I need the else here?
		{
			//we have left the hom region
			//std::cout << "previous was left anchor" << '\n';
			//if (anchor_found != "none")
			//{
			//previous_anchor_found = anchor_found;
			//}
			anchor_found = "left";
			return 1;
		}
		//else previous kmer is not a left anchor (and current kmer is not a left anchor)
		else
		{
			//we continue the hom region
			if (anchor_found != "none")
			{
			previous_anchor_found = anchor_found;
			}
			anchor_found = "none";
			return 2;
		}
	}
}

std::vector<std::string> get_paths(CKMCFile& file, int& error_threshold, int& het_threshold, int& unique_threshold, std::string& left_anchor_kmer, std::string& right_anchor_kmer, int& min_distance_of_path, int& max_distance_of_path, int& max_nodes_to_search, int& k, bool& queue_broken)
{
	std::string starting_anchor_kmer;
	bool going_right;
	if (left_anchor_kmer.empty())
	{
		//We are at the beginning of the read.
		//This function will find paths starting from right_anchor_kmer continuing left to the
		//beginning of the read where each kmer in the path is a nonerror kmer.
		starting_anchor_kmer = right_anchor_kmer;
		going_right = false;
	}
	else if (right_anchor_kmer.empty())
	{
		//We are at the end of the read.
		//This function will find paths starting from left_anchor_kmer continuing right to the
		//end of the read where each kmer in the path is a nonerror kmer.
		starting_anchor_kmer = left_anchor_kmer;
		going_right = true;
	}
	else
	{
		//We are in the middle of the read.
		//This function will find paths starting from left_anchor_kmer continuing right until
		//right_anchor_kmer where each kmer in the path is a nonerror kmer.
		starting_anchor_kmer = left_anchor_kmer;
		going_right = true;
	}
	//In every case, we follow all paths, but cut the depth of any path to max_distance_of_path.
	//We also ensure a minimum depth equal to min_distance_of_path.
	std::list<std::string> queue;
	queue.push_back(starting_anchor_kmer);
	//Initialize paths to store all the paths that are found.
	std::vector<std::string> paths;
	//We use i as a counter for how many nodes have been visited in the search.
	//If we haven't finished the search within max_nodes_to_search nodes, 
	//we break the search.
	//This drastically speeds up the run time for some regions.
	//Thankfully, it doesn't seem to impact effectiveness, since most searches
	//complete before this threshold.
	int i = 0;
	//This flag keeps track of whether we had to stop the search early.
	queue_broken = false;
	while(!queue.empty())
	{
		i++;
		std::string current_path = queue.front();
		std::string current_kmer;
		if (going_right)
		{
			current_kmer = current_path.substr(current_path.length()-k);
		}
		else
		{
			current_kmer = current_path.substr(0, k);
		}
		queue.pop_front();
		int current_depth = current_path.length()-k;
		//If we have to terminate search early
		if (i > max_nodes_to_search)
		{
			//std::cout << "queue broken" << '\n';
			queue_broken = true;
			break;
		}
		//If the depth of this node hasn't exceeded the max distance of the path
		if (current_depth <= max_distance_of_path)
		{
			//Extend the path by one nucleotide, keep the ones that are not error kmers
			std::vector<std::string> adjacent_kmers = get_adjacent(file, current_kmer, error_threshold, het_threshold, unique_threshold, going_right);
			for (auto adjacent_kmer : adjacent_kmers)
			{
				std::string path;
				if (going_right)
				{
					path = current_path + adjacent_kmer.back();
				}
				else
				{
					path = adjacent_kmer.front() + current_path;
				}
				bool end_condition;
				//If we are in the middle of the read, we end when we have found a path
				//of nonerror kmers which bridges the anchor kmers
				//and doesn't terminate too early (i.e. before min_distance_of_path)
				if (!left_anchor_kmer.empty() && !right_anchor_kmer.empty())
				{
					end_condition = ((adjacent_kmer == right_anchor_kmer) && (current_depth + 1 >= min_distance_of_path));
				}
				//If we are at either end of the read, we end when we have found a path
				//of nonerror kmers which continues until the end of the read
				//and doesn't terminate too early (i.e. before min_distance_of_path)
				else
				{
					end_condition = ((current_depth + 1 == max_distance_of_path) && (current_depth + 1 >= min_distance_of_path));
				}
				if (end_condition)
				{
					if (!left_anchor_kmer.empty())
					{
						//remove left_anchor_kmer from path
						path.erase(path.begin(), path.begin()+k);
					}
					if (!right_anchor_kmer.empty())
					{
						//remove right_anchor_kmer from path
						path.erase(path.end()-k, path.end());
					}
					paths.push_back(path);
				}
				//Else we haven't found a path yet
				else
				{
					queue.push_back(path);
				}
			}
		}
	}
	return paths;
}

//bool IsErr (uint32_t coverage)
//{
//	int error_threshold = 7; //9 for atha, 5 for simulated_hifi, 4 for amaranth, 8 for dmel, 12 for msie, 7 for msyl
//	return (coverage <= error_threshold);
//}

//bool IsRep (uint32_t coverage)
//{
//	int unique_threshold = 60; //72 for atha, 50 for simulated_hifi, 60 for amaranth, 60 for dmel, 100 for msie, 60 for msyl
//	return (coverage > unique_threshold);
//}

void write_error_paths(bool& queue_broken, std::vector<std::string>& edited_error_portions, std::ofstream& erredits_output_file, std::ofstream& errpaths_output_file, std::string& edited_read, int& read_number, int& first_error_idx, int& last_error_idx, std::string& before_first_error_kmer, std::string& original_error_portion, std::string& after_last_error_kmer, CKMCFile& file)
{
	std::ofstream* errwrite_output_file;
	//we finished the search and presumably we have found one homozygous path or two heterozygous paths
	if ((!queue_broken) && ((edited_error_portions.size() == 1) or (edited_error_portions.size() == 2)))
	{
		errwrite_output_file = &erredits_output_file;
		if (!before_first_error_kmer.empty())
		{
			edited_read += before_first_error_kmer.substr(1) + edited_error_portions[0];
		}
		else
		{
			edited_read += edited_error_portions[0];
		}
	}
	//there are no paths, or there are more than two paths, or the search wasn't finished.
	//We are currently not editing.
	else
	{
		errwrite_output_file = &errpaths_output_file;
		if (!before_first_error_kmer.empty())
		{
			edited_read += before_first_error_kmer.substr(1) + original_error_portion;
		}
		else
		{
			edited_read += original_error_portion;
		}
	}
	std::string original_error_block = before_first_error_kmer + original_error_portion + after_last_error_kmer;
	*errwrite_output_file << ">read" << read_number << "_firsterrorkmer" << first_error_idx << "_lasterrorkmer" << last_error_idx << "_original" << '\n';
	*errwrite_output_file << before_first_error_kmer << " " << original_error_portion << " " << after_last_error_kmer << '\n';
	std::vector<uint32_t> w;
	file.GetCountersForRead(original_error_block, w);
	for (int j=0; j < w.size(); j++)
	{
		*errwrite_output_file << w.at(j) << " ";
	}
	*errwrite_output_file << '\n';
	for (int l = 0; l < edited_error_portions.size(); l++)
	{
		std::string edited_error_portion = edited_error_portions[l];
		std::string edited_error_block = before_first_error_kmer + edited_error_portion + after_last_error_kmer;
		*errwrite_output_file << ">read" << read_number << "_firsterrorkmer" << first_error_idx << "_lasterrorkmer" << last_error_idx << "_edited" << l << '\n';
		*errwrite_output_file << before_first_error_kmer << " " << edited_error_portion << " " << after_last_error_kmer << '\n';
		file.GetCountersForRead(edited_error_block, w);
		for (int j=0; j < w.size(); j++)
		{
			*errwrite_output_file << w.at(j) << " ";
		}
		*errwrite_output_file << '\n';
	}
}

void write_nonhom_paths(bool& queue_broken, std::vector<std::string>& smoothed_nonhom_portions, std::ofstream& hetedits_output_file, std::ofstream& hetpaths_output_file, std::string& smoothed_read, int& read_number, int& first_nonhom_idx, int& last_nonhom_idx, std::string& before_first_nonhom_kmer, std::string& original_nonhom_portion, std::string& after_last_nonhom_kmer, CKMCFile& file)
{	
	std::string left_part;
	if (!before_first_nonhom_kmer.empty())
	{
		left_part = before_first_nonhom_kmer.substr(1);
	}
	else
	{
		left_part = before_first_nonhom_kmer;
	}
	std::string right_part;
	if (!after_last_nonhom_kmer.empty())
	{
		right_part = after_last_nonhom_kmer.substr(0, after_last_nonhom_kmer.length()-1);
	}
	else
	{
		right_part = after_last_nonhom_kmer;
	}
	std::ofstream* hetwrite_output_file;
	//we finished the search and presumably we have found two heterozygous paths
	if ((!queue_broken) && (smoothed_nonhom_portions.size() == 2))
	{
		hetwrite_output_file = &hetedits_output_file;
		std::string smoothed_nonhom_portion0 = smoothed_nonhom_portions[0];
		std::string smoothed_nonhom_portion1 = smoothed_nonhom_portions[1];
		std::vector<uint32_t> w;
		file.GetCountersForRead(left_part + smoothed_nonhom_portion0 + right_part, w);
		float average0 = std::accumulate(w.begin(), w.end(), 0.0)/w.size();
		file.GetCountersForRead(left_part + smoothed_nonhom_portion1 + right_part, w);
		float average1 = std::accumulate(w.begin(), w.end(), 0.0)/w.size();
		std::string primary_smoothed_nonhom_portion;
		std::string alternate_smoothed_nonhom_portion;
		if (average0 > average1)
		{
			primary_smoothed_nonhom_portion = smoothed_nonhom_portion0;
			alternate_smoothed_nonhom_portion = smoothed_nonhom_portion1;
		}
		else
		{
			primary_smoothed_nonhom_portion = smoothed_nonhom_portion1;
			alternate_smoothed_nonhom_portion = smoothed_nonhom_portion0;
		}
		smoothed_read += left_part + primary_smoothed_nonhom_portion;
		std::cout << left_part + primary_smoothed_nonhom_portion << '\n';
		smoothed_nonhom_portions[0] = primary_smoothed_nonhom_portion;
		smoothed_nonhom_portions[1] = alternate_smoothed_nonhom_portion;
	}
	//there is not exactly two paths, or the search wasn't finished.
	//We are currently not smoothing.
	else
	{
		hetwrite_output_file = &hetpaths_output_file;
		smoothed_read += left_part + original_nonhom_portion;
		std::cout << left_part + original_nonhom_portion << '\n';
	}
	std::string original_nonhom_block = before_first_nonhom_kmer + original_nonhom_portion + after_last_nonhom_kmer;
	*hetwrite_output_file << ">read" << read_number << "_firstnonhomkmer" << first_nonhom_idx << "_lastnonhomkmer" << last_nonhom_idx << "_original" << '\n';
	*hetwrite_output_file << before_first_nonhom_kmer << " " << original_nonhom_portion << " " << after_last_nonhom_kmer << '\n';
	std::vector<uint32_t> w;
	file.GetCountersForRead(original_nonhom_block, w);
	for (int j=0; j < w.size(); j++)
	{
		*hetwrite_output_file << w.at(j) << " ";
	}
	*hetwrite_output_file << '\n';
	for (int l = 0; l < smoothed_nonhom_portions.size(); l++)
	{
		std::string smoothed_nonhom_portion = smoothed_nonhom_portions[l];
		std::string smoothed_nonhom_block = before_first_nonhom_kmer + smoothed_nonhom_portion + after_last_nonhom_kmer;
		*hetwrite_output_file << ">read" << read_number << "_firstnonhomkmer" << first_nonhom_idx << "_lastnonhomkmer" << last_nonhom_idx << "_edited" << l << '\n';
		*hetwrite_output_file << before_first_nonhom_kmer << " " << smoothed_nonhom_portion << " " << after_last_nonhom_kmer << '\n';
		file.GetCountersForRead(smoothed_nonhom_block, w);
		for (int j=0; j < w.size(); j++)
		{
			*hetwrite_output_file << w.at(j) << " ";
		}
		*hetwrite_output_file << '\n';
	}
}

std::string remove_err (std::vector<uint32_t>& v, std::string& read, int& read_number, CKMCFile& file, std::ofstream& erredits_output_file, std::ofstream& errpaths_output_file, int& error_threshold, int& het_threshold, int& unique_threshold, int& max_nodes_to_search, double& distance_multiplier)
{
	//initialize variables
	std::string edited_read;
	int k = 21;
	//int error_threshold = 7; //8 for dmel, 9 for atha, 5 for simulated_hifi, 4 for amaranth, 8 for dmel, 12 for msie, 7 for msyl
	//int het_threshold = 29; //21 for dmel, 30 for atha, 20 for simulated_hifi at 0.5% error, 29 for amaranth, 29 for dmel, 48 for msie, 29 for msyl
	//int unique_threshold = 60; //60 for dmel, 72 for atha, 50 for simulated_hifi, 60 for amaranth, 60 for dmel, 100 for msie, 60 for msyl

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
		//if kmer is nonerror
		if (current_type > 0)
		{
			//if this is the first kmer of the read
			if (previous_type == -1)
			{
				first_nonerror_idx = i;
				last_error_idx = i-1;
			}
			//if previous kmer was error, and we are at the beginning of the read
			if (previous_type == 0 && before_first_error_kmer.empty())
			{
				//The very beginning of the read is an error portion
				after_last_error_kmer = read.substr(i, k);
				int min_distance_of_path = 0;
				int max_distance_of_path = i;
				//int max_nodes_to_search = 1000;
				bool queue_broken = false;
				std::vector<std::string> edited_error_portions = get_paths(file, error_threshold, het_threshold, unique_threshold, before_first_error_kmer, after_last_error_kmer, min_distance_of_path, max_distance_of_path, max_nodes_to_search, k, queue_broken);
				last_error_idx = i-1;
				first_nonerror_idx = i;
				std::string original_error_portion = read.substr(0, i);
				write_error_paths(queue_broken, edited_error_portions, erredits_output_file, errpaths_output_file, edited_read, read_number, first_error_idx, last_error_idx, before_first_error_kmer, original_error_portion, after_last_error_kmer, file);
			}
			//if previous kmer was error, we have left the error block
			if (previous_type == 0 && !before_first_error_kmer.empty())
			{
				int number_of_error_kmers = i - first_error_idx;
				//If the position of after_last_error_kmer overlaps before_first_error_kmer
				//we keep progressing as if nothing has happened, waiting to find another non_error kmer
				if (number_of_error_kmers < k)
				{
					current_type = previous_type;
					continue;
				}
				//get kmer that is right after the last error kmer of block
				after_last_error_kmer = read.substr(i, k);
				int min_distance_of_path = k;
				int max_distance_of_path = ceil(distance_multiplier * number_of_error_kmers);
				//int max_nodes_to_search = 1000;
				bool queue_broken;
				std::vector<std::string> edited_error_portions = get_paths(file, error_threshold, het_threshold, unique_threshold, before_first_error_kmer, after_last_error_kmer, min_distance_of_path, max_distance_of_path, max_nodes_to_search, k, queue_broken);
				last_error_idx = i-1;
				first_nonerror_idx = i;
				std::string original_error_portion = read.substr(first_error_idx+k-1, last_error_idx - first_error_idx + 2 - k);
				write_error_paths(queue_broken, edited_error_portions, erredits_output_file, errpaths_output_file, edited_read, read_number, first_error_idx, last_error_idx, before_first_error_kmer, original_error_portion, after_last_error_kmer, file);
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
		//We have "left" the error portion of the read
		after_last_error_kmer = "";
		int min_distance_of_path = 0;
		int max_distance_of_path = v.size()-first_error_idx;
		//int max_nodes_to_search = 1000;
		bool queue_broken = false;
		std::vector<std::string> edited_error_portions = get_paths(file, error_threshold, het_threshold, unique_threshold, before_first_error_kmer, after_last_error_kmer, min_distance_of_path, max_distance_of_path, max_nodes_to_search, k, queue_broken);
		last_error_idx = v.size()-1;
		first_nonerror_idx = v.size();
		std::string original_error_portion = read.substr(first_error_idx+k-1);
		write_error_paths(queue_broken, edited_error_portions, erredits_output_file, errpaths_output_file, edited_read, read_number, first_error_idx, last_error_idx, before_first_error_kmer, original_error_portion, after_last_error_kmer, file);
	}
	if (previous_type > 0)
	{
		//We have "left" the nonerror portion of the read
		first_error_idx = v.size();
		last_nonerror_idx = v.size()-1;
		std::string nonerror_portion = read.substr(first_nonerror_idx, last_nonerror_idx - first_nonerror_idx + k);
		edited_read += nonerror_portion;
	}
	edited_read += '\n';
	return edited_read;
}

std::string smooth_het (std::vector<uint32_t>& v, std::string& read, int& read_number, CKMCFile& file, std::ofstream& hetedits_output_file, std::ofstream& hetpaths_output_file, int& error_threshold, int& het_threshold, int& unique_threshold, int& max_nodes_to_search, double& distance_multiplier)
{
	//initialize variables
	std::string smoothed_read;
	int k = 21;
	//int error_threshold = 7; //8 for dmel, 9 for atha, 5 for simulated_hifi, 4 for amaranth, 8 for dmel, 12 for msie, 7 for msyl
	//int het_threshold = 29; //21 for dmel, 30 for atha, 20 for simulated_hifi at 0.5% error, 29 for amaranth, 29 for dmel, 48 for msie, 29 for msyl
	//int unique_threshold = 60; //60 for dmel, 72 for atha, 50 for simulated_hifi, 60 for amaranth, 60 for dmel, 100 for msie, 60 for msyl

	//iterate over counts to smoothe het
	int previous_type = -1;
	int first_hom_idx;
	int last_hom_idx;
	int first_nonhom_idx;
	int last_nonhom_idx;
	std::string before_first_nonhom_kmer;
	std::string after_last_nonhom_kmer;
	for (int i = 0; i < v.size(); i++)
	{
		std::cout << i << '\t' << v[i] << '\n';
		//int current_type = get_type(v[i], error_threshold, het_threshold, unique_threshold);
		std::string previous_kmer;
		int previous_count;
		if (i == 0)
		{
			previous_kmer = "";
			previous_count = 0;
		}
		else
		{
			previous_kmer = read.substr(i-1, k);
			std::vector<uint32_t> previous_count_vector;
			file.GetCountersForRead(previous_kmer, previous_count_vector);
			previous_count = previous_count_vector[0];
		}
		std::string current_kmer = read.substr(i, k);
		//std::cout << current_kmer << '\n';
		std::vector<uint32_t> current_count_vector;
		file.GetCountersForRead(current_kmer, current_count_vector);
		int current_count = current_count_vector[0];
		std::string anchor_found;
		std::string previous_anchor_found;
		int current_type = get_type_het(previous_type, previous_kmer, current_kmer, previous_count, current_count, k, file, error_threshold, het_threshold, unique_threshold, anchor_found, previous_anchor_found);
		//int current_type;
		//if (i==0)
		//{
		//	current_type = 2;
		//}
		//else
		//{
		//	current_type = 1;
		//}
		//if kmer is nonhom
		if ((current_type == 0) || (current_type == 1) || (current_type == 3))
		{
			//if this is the first kmer of the read
			if (previous_type == -1)
			{
				first_nonhom_idx = i;
				last_hom_idx = i-1;
			}
			//if previous kmer was nonhom, we are continuing the nonhom block
			if ((previous_type == 0) || (previous_type == 1) || (previous_type == 3))
			{
				;
			}
			//if previous kmer was hom, we are leaving the hom block
			if (previous_type == 2)
			{
				//get kmer that is right before the first nonhom kmer of nonhom block
				first_nonhom_idx = i;
				last_hom_idx = i-1;
				before_first_nonhom_kmer = read.substr(i-1, k);
				std::string hom_portion = read.substr(first_hom_idx, last_hom_idx - first_hom_idx + 1);
				smoothed_read += hom_portion;
				std::cout << hom_portion << '\n';
			}
			previous_type = current_type;
		}
		//if kmer is hom
		if (current_type == 2)
		{
			//if this is the first kmer of the read
			if (previous_type == -1)
			{
				first_hom_idx = i;
				last_nonhom_idx = i-1;
			}
			//if previous kmer was nonhom, and we are at the beginning of the read
			if (((previous_type == 0) || (previous_type == 1) || (previous_type == 3)) && before_first_nonhom_kmer.empty())
			{
				//The very beginning of the read is an nonhom portion
				after_last_nonhom_kmer = read.substr(i, k);
				int min_distance_of_path = 0;
				int max_distance_of_path = i;
				//int max_nodes_to_search = 1000;
				bool queue_broken = false;
				std::vector<std::string> smoothed_nonhom_portions = get_paths(file, error_threshold, het_threshold, unique_threshold, before_first_nonhom_kmer, after_last_nonhom_kmer, min_distance_of_path, max_distance_of_path, max_nodes_to_search, k, queue_broken);
				last_nonhom_idx = i-1;
				first_hom_idx = i;
				std::string original_nonhom_portion = read.substr(0, i);
				write_nonhom_paths(queue_broken, smoothed_nonhom_portions, hetedits_output_file, hetpaths_output_file, smoothed_read, read_number, first_nonhom_idx, last_nonhom_idx, before_first_nonhom_kmer, original_nonhom_portion, after_last_nonhom_kmer, file);
			}
			//if previous kmer was nonhom, we have left the nonhom block
			if (((previous_type == 0) || (previous_type == 1) || (previous_type == 3)) && !before_first_nonhom_kmer.empty())
			{
				if (anchor_found == "left")
				{
					after_last_nonhom_kmer = read.substr(i, k);
					std::string hom_portion;
					if (previous_anchor_found == "left")
					{
						hom_portion = read.substr(first_hom_idx+1, i-first_hom_idx-1);
					}
					else
					{
						std::cout << "OMG!" << '\n';
						hom_portion = read.substr(first_hom_idx+1, i-first_hom_idx-1);
					}
					std::cout << hom_portion << '\n';
					last_nonhom_idx = i-1;
					first_hom_idx = i;
					smoothed_read += hom_portion;
				}
				else
				{
				int number_of_nonhom_kmers = i - first_nonhom_idx;
				//If the position of after_last_nonhom_kmer overlaps before_first_nonhom_kmer
				//we keep progressing as if nothing has happened, waiting to find another hom kmer
				if ((anchor_found == "right") && (number_of_nonhom_kmers < k))
				{
					std::cout << "left anchor and right anchor overlap" << '\n';
					current_type = previous_type;
					continue;
				}
				//get kmer that is right after the last nonhom kmer of block
				after_last_nonhom_kmer = read.substr(i, k);
				int min_distance_of_path = k;
				int max_distance_of_path = ceil(distance_multiplier * number_of_nonhom_kmers);
				//int max_nodes_to_search = 1000;
				bool queue_broken;
				std::vector<std::string> smoothed_nonhom_portions = get_paths(file, error_threshold, het_threshold, unique_threshold, before_first_nonhom_kmer, after_last_nonhom_kmer, min_distance_of_path, max_distance_of_path, max_nodes_to_search, k, queue_broken);
				last_nonhom_idx = i-1;
				first_hom_idx = i;
				std::string original_nonhom_portion = read.substr(first_nonhom_idx+k-1, last_nonhom_idx - first_nonhom_idx + 2 - k);
				write_nonhom_paths(queue_broken, smoothed_nonhom_portions, hetedits_output_file, hetpaths_output_file, smoothed_read, read_number, first_nonhom_idx, last_nonhom_idx, before_first_nonhom_kmer, original_nonhom_portion, after_last_nonhom_kmer, file);
				}
			}
			//if previous kmer is hom, we are continuing a hom block
			if (previous_type == 2)
			{
				if (anchor_found == "left")
				{
					std::string hom_portion = read.substr(first_hom_idx, i-first_hom_idx);
					smoothed_read += hom_portion;
					std::cout << hom_portion << '\n';
					std::cout << "this is new" << '\n';
					first_hom_idx=i;
					last_nonhom_idx=i-1;
				}
			}
			previous_type = current_type;
		}
	}
	//We have reached the end of the read, let's make sure we have added the last bit of the read
	if ((previous_type == 0) || (previous_type == 1) || (previous_type == 3))
	{
		//We have "left" the nonhom portion of the read
		//If we have a homozygous on left with which to anchor
		if (first_nonhom_idx > 0)
		{
			after_last_nonhom_kmer = "";
			int min_distance_of_path = 0;
			int max_distance_of_path = v.size()-first_nonhom_idx;
			//max_distance_of_path = 50;
			//int max_nodes_to_search = 1000;
			bool queue_broken = false;
			std::vector<std::string> smoothed_nonhom_portions = get_paths(file, error_threshold, het_threshold, unique_threshold, before_first_nonhom_kmer, after_last_nonhom_kmer, min_distance_of_path, max_distance_of_path, max_nodes_to_search, k, queue_broken);
			last_nonhom_idx = v.size()-1;
			first_hom_idx = v.size();
			std::string original_nonhom_portion = read.substr(first_nonhom_idx+k-1);
			write_nonhom_paths(queue_broken, smoothed_nonhom_portions, hetedits_output_file, hetpaths_output_file, smoothed_read, read_number, first_nonhom_idx, last_nonhom_idx, before_first_nonhom_kmer, original_nonhom_portion, after_last_nonhom_kmer, file);
		}
		else
		{
			smoothed_read += read;
		}
	}
	if (previous_type == 2)
	{
		//We have "left" the hom portion of the read
		first_nonhom_idx = v.size();
		last_hom_idx = v.size()-1;
		std::string hom_portion = read.substr(first_hom_idx, last_hom_idx - first_hom_idx + k);
		smoothed_read += hom_portion;
		std::cout << hom_portion << '\n';
	}
	smoothed_read += '\n';
	return smoothed_read;
}

std::string getFileExt (const std::string &s)
{
	size_t i = s.find_last_of('.');
	if (i != std::string::npos)
	{
		return(s.substr(i+1, s.length() - i));
	}
	return("");
}

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
	int error_threshold = atoi(argv[9]);
	int het_threshold = atoi(argv[10]);
	int unique_threshold = atoi(argv[11]);
	double allowed_err_fraction = std::stod(argv[12]);
	double allowed_rep_fraction = std::stod(argv[13]);
	int max_nodes_to_search = atoi(argv[14]);
	double distance_multiplier = std::stod(argv[15]);

	//load reads
	int line_num = -1;
	
	//is the input fasta or fastq?
	//note: the output will be fasta format, since quality values
	//will not match once the read is edited and smoothed
	int num_lines_per_read;
	if ((getFileExt(argv[2]) == "fasta") || (getFileExt(argv[2]) == "fa"))
	{
		num_lines_per_read = 2;
	}
	if ((getFileExt(argv[2]) == "fastq") || (getFileExt(argv[2]) == "fq"))
	{
		num_lines_per_read = 4;
	}
	
	std::string header;
	while (getline(input_file, line))
	{
		line_num++;
		//if on line with the header
		if (line_num % num_lines_per_read == 0)
		{
			//save header, but make sure it is a fasta header, not a fastq header
			header = ">" + line.substr(1);
		}
		//if on line with the read
		if (line_num % num_lines_per_read == 1)
		{
			std::string read = line;
			int read_number = (line_num+num_lines_per_read-1)/num_lines_per_read;
			if (read_number%10000==0)
			{
				std::cout << read_number << '\n';
			}
			
			//get counters of kmers in read
			std::vector<uint32_t> v;
			file.GetCountersForRead(read, v);
			int num_kmers = v.size();
			
			//check if too many error kmers
			//if greater than allowed_err_fraction of kmers in the read are errors
			//discard the read to not waste time
			//int num_err = std::count_if(v.begin(), v.end(), IsErr);
			int num_err = 0;
			for (int i = 0; i < v.size(); i++)
			{
				if (v[i] <= error_threshold)
				{
					num_err++;
				}
			}
			double err_fraction = static_cast<double>(num_err)/num_kmers;
			if (err_fraction >= allowed_err_fraction)
			{
				//std::cout << "read number " << read_number << " has " << err_fraction << " percent error kmers, discarding." << '\n';
				continue;
			}

			//remove errors from the read to get edited read
			std::string edited_read = remove_err(v, read, read_number, file, erredits_output_file, errpaths_output_file, error_threshold, het_threshold, unique_threshold, max_nodes_to_search, distance_multiplier);

			//write header and edited read to err_output_file
			err_output_file << header << '\n';
			err_output_file << edited_read;
			edited_read.pop_back();

			//get counters of kmers in edited read
			file.GetCountersForRead(edited_read, v);
			num_kmers = v.size();

			//check if too many repetitive kmers
			//if greater than allowed_rep_fraction of kmers in the edited read are repetitive
			//discard the read to not waste time
			//int num_rep = std::count_if(v.begin(), v.end(), IsRep);
			int num_rep = 0;
			for (int i = 0; i < v.size(); i++)
			{
				if (v[i] > unique_threshold)
				{
					num_rep++;
				}
			}
			double rep_fraction = static_cast<double>(num_rep)/num_kmers;
			if (rep_fraction >= allowed_rep_fraction)
			{
				//std::cout << "read number " << read_number << " has " << rep_fraction << " percent repetitive kmers, discarding." << '\n';
				het_output_file << header << '\n';
				het_output_file << edited_read << '\n';
				continue;
			}

			//smoothe het from the edited read to get smoothed read
			std::string smoothed_read = smooth_het(v, edited_read, read_number, file, hetedits_output_file, hetpaths_output_file, error_threshold, het_threshold, unique_threshold, max_nodes_to_search, distance_multiplier);

			//write header and smoothed read to het_output_file
			het_output_file << header << '\n';
			het_output_file << smoothed_read;
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
