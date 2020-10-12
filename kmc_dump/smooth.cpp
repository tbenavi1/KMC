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

std::vector<std::string> get_paths(CKMCFile& file, int& error_threshold, int& het_threshold, int& unique_threshold, std::string& left_anchor_kmer, std::string& right_anchor_kmer, int& min_distance_of_path, int& max_distance_of_path, int& max_nodes_to_search, int& k, bool& queue_broken)
{
	//This function finds paths starting from left_anchor_kmer and ending at
	//right_anchor_kmer where each kmer in the path is a nonerror kmer.
	//We follow all paths, but cut the depth of any path to max_distance_of_path.
	//We also ensure a minimum depth equal to min_distance_of_path.
	std::list<std::string> queue;
	queue.push_back(left_anchor_kmer);
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
		//If the depth of this node hasn't exceeded the max distance of the path
		if (current_depth <= max_distance_of_path)
		{
			//Extend the path by one nucleotide, keep the ones that are not error kmers
			bool going_right = true;
			std::vector<std::string> adjacent_kmers = get_adjacent(file, current_kmer, error_threshold, het_threshold, unique_threshold, going_right);
			for (auto adjacent_kmer : adjacent_kmers)
			{
				std::string path = current_path + adjacent_kmer.back();
				//If we have found a path of nonerror kmers which bridges the anchor kmers
				//and doesn't terminate too early (i.e. before min_distance_of_path)
				if ((adjacent_kmer == right_anchor_kmer) && (current_depth + 1 >= min_distance_of_path))
				{
					path.erase(path.end()-k, path.end());
					path.erase(path.begin(), path.begin()+k);
					//paths.push_back(path.substr(1));
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

std::vector<std::string> get_paths_het(CKMCFile& file, int& error_threshold, int& het_threshold, int& unique_threshold, std::string& left_anchor_kmer, std::string& right_anchor_kmer, int& min_distance_of_path, int& max_distance_of_path, int& max_nodes_to_search, int& k, bool& queue_broken)
{
	//This function finds paths starting from left_anchor_kmer and ending at
	//right_anchor_kmer where each kmer in the path is a nonerror kmer.
	//We follow all paths, but cut the depth of any path to max_distance_of_path.
	//We also ensure a minimum depth equal to min_distance_of_path.
	std::list<std::string> queue;
	queue.push_back(left_anchor_kmer);
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
		//If the depth of this node hasn't exceeded the max distance of the path
		if (current_depth <= max_distance_of_path)
		{
			//Extend the path by one nucleotide, keep the ones that are not error kmers
			bool going_right = true;
			std::vector<std::string> adjacent_kmers = get_adjacent(file, current_kmer, error_threshold, het_threshold, unique_threshold, going_right);
			for (auto adjacent_kmer : adjacent_kmers)
			{
				std::string path = current_path + adjacent_kmer.back();
				//If we have found a path of nonerror kmers which bridges the anchor kmers
				//and doesn't terminate too early (i.e. before min_distance_of_path)
				if ((adjacent_kmer == right_anchor_kmer) && (current_depth + 1 >= min_distance_of_path))
				{
					path.erase(path.end()-k, path.end());
					//path.erase(path.begin(), path.begin()+k);
					paths.push_back(path.substr(1));
					//paths.push_back(path);
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

std::vector<std::string> get_paths_left(CKMCFile& file, int& error_threshold, int& het_threshold, int& unique_threshold, std::string& right_anchor_kmer, int& min_distance_of_path, int& max_distance_of_path, int& max_nodes_to_search, int& k, bool& queue_broken)
{
	//This function finds paths starting from right_anchor_kmer and
	//going to the left where each kmer in the path is a nonerror kmer.
	//We follow all paths, but cut the depth of any path to max_distance_of_path.
	//We also ensure a minimum depth equal to min_distance_of_path.
	std::list<std::string> queue;
	queue.push_back(right_anchor_kmer);
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
		//std::string current_kmer = current_path.substr(current_path.length()-k);
		std::string current_kmer = current_path.substr(0, k);
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
			bool going_right = false;
			std::vector<std::string> adjacent_kmers = get_adjacent(file, current_kmer, error_threshold, het_threshold, unique_threshold, going_right);
			for (auto adjacent_kmer : adjacent_kmers)
			{
				//std::string path = current_path + adjacent_kmer.back();
				std::string path = adjacent_kmer.front() + current_path;
				//std::cout << path << '\n';
				//If we have found a path of nonerror kmers which bridges the anchor kmers
				//and doesn't terminate too early (i.e. before min_distance_of_path)
				//if ((adjacent_kmer == right_anchor_kmer) && (current_depth + 1 >= min_distance_of_path))
				if ((current_depth + 1 == max_distance_of_path) && (current_depth + 1 >= min_distance_of_path))
				{
					path.erase(path.end()-k, path.end());
					//path.erase(path.begin(), path.begin()+k);
					//paths.push_back(path.substr(1));
					//path.pop_back();
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

std::vector<std::string> get_paths_right(CKMCFile& file, int& error_threshold, int& het_threshold, int& unique_threshold, std::string& left_anchor_kmer, int& min_distance_of_path, int& max_distance_of_path, int& max_nodes_to_search, int& k, bool& queue_broken)
{
	//This function finds paths starting from left_anchor_kmer and
	//going to the right where each kmer in the path is a nonerror kmer.
	//We follow all paths, but cut the depth of any path to max_distance_of_path.
	//We also ensure a minimum depth equal to min_distance_of_path.
	std::list<std::string> queue;
	queue.push_back(left_anchor_kmer);
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
		std::string current_kmer = current_path.substr(current_path.length()-k);
		//std::string current_kmer = current_path.substr(0, k);
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
			bool going_right = true;
			std::vector<std::string> adjacent_kmers = get_adjacent(file, current_kmer, error_threshold, het_threshold, unique_threshold, going_right);
			for (auto adjacent_kmer : adjacent_kmers)
			{
				std::string path = current_path + adjacent_kmer.back();
				//std::string path = adjacent_kmer.front() + current_path;
				//std::cout << path << '\n';
				//If we have found a path of nonerror kmers which bridges the anchor kmers
				//and doesn't terminate too early (i.e. before min_distance_of_path)
				//if ((adjacent_kmer == right_anchor_kmer) && (current_depth + 1 >= min_distance_of_path))
				if ((current_depth + 1 == max_distance_of_path) && (current_depth + 1 >= min_distance_of_path))
				{
					//path.erase(path.end()-k, path.end());
					path.erase(path.begin(), path.begin()+k);
					//paths.push_back(path.substr(1));
					//path.pop_back();
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
			return extension;
		}
		kmer = new_kmer;
	}
	//we terminate at the end of the read
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
					extension.pop_back();
					extension = extension + ending_portion.substr(j+k);
					return extension;
				}
				else
				{
					found_path = true;
					j++;
					new_kmer = adjacent_kmer;
					extension = extension+c;
				}
			}
		}
		//we terminate at no paths
		if (!found_path)
		{
			extension = extension + ending_portion.substr(j+k);
			return extension;
		}
		kmer = new_kmer;
	}
	//we terminate at the end of the read
	return extension;
}

//std::vector<std::string> extend_left_het(std::string& beginning_portion, int& i, int& k, int& error_threshold, int& het_threshold, int& unique_threshold, CKMCFile& file, int& read_number)
//{
//	std::vector<std::string> extensions;
//	//The beginning portion should be of length i+k
//	std::string kmer1 = beginning_portion.substr(i);
//	std::string kmer2 = beginning_portion.substr(i);
//	std::vector<uint32_t> v1;
//	std::vector<uint32_t> v2;
//	int j = 0;
//	bool found_branch = false;
//	while (j < i)
//	{
//		//we are before the branch
//		if (!found_branch)
//		{
//			//bool found_path = false;
//			int num_paths_found = 0;
//			//std::string new_kmer;
//			std::string new_kmer1;
//			//std::string new_kmer2;
//			for (char const &c: "ACGT")
//			{
//				std::string adjacent_kmer1 = c+kmer1.substr(0, k-1);
//				//std::string adjacent_kmer2 = c+kmer.substr(0, k-1);
//				file.GetCountersForRead(adjacent_kmer1, v1);
//				//file.GetCountersForRead(adjacent_kmer2, v2);
//				int current_type1 = get_type(v1[0], error_threshold, het_threshold, unique_threshold);
//				//int current_type2 = get_type(v2[0], error_threshold, het_threshold, unique_threshold);
//				//if the adjacent kmer is not an error
//				if (current_type1 > 0)
//				{
//					num_paths_found++;
//					if (num_paths_found == 1)
//					{
//						j++;
//						new_kmer1 = adjacent_kmer1;
//						new_kmer2 = adjacent_kmer1;
//						//extension = c+extension;
//						extensions[0] = c + extensions[0];
//						extensions[1] = c + extensions[1];
//					}
//					if (num_paths_found == 2)
//					{
//						//j--;
//						//remove the previous character we added
//						//extension = extension.substr(1);
//						//extension = beginning_portion.substr(0, i-j) + extension;
//						new_kmer2 = adjacent_kmer1;
//						extensions[1] = extensions[1].substr(1);
//						extensions[1] = c + extensions[1];
//						found_branch = true;
//						//extensions[1] = beginning_portion.substr(0, i-j) + extensions[1];
//						//std::cout << "we terminate at a branching path" << '\n';
//						//return extension;
//					}
//					//We terminate since there are too many paths
//					//return the original sequence, with correct length
//					if (num_paths_found > 2)
//					{
//						extensions[0] = beginning_portion.substr(0, i);
//						extensions[1] = beginning_portion.substr(0, i);
//						return extensions;
//					}
//				}
//			}
//			//we terminate at no paths
//			if (num_paths_found == 0)
//			{
//				//extension = beginning_portion.substr(0, i-j) + extension;
//				extensions[0] = beginning_portion.substr(0, i);
//				extensions[1] = beginning_portion.substr(0, i);
//				//std::cout << "we terminate at no paths" << '\n';
//				return extensions;
//			}
//			kmer1 = new_kmer1;
//			kmer2 = new_kmer2;
//		}
//		//we are after the branch
//		else
//		{
//			if (j > 1)
//			{
//				std::cout << "the branch was not the first position, j: " << j << '\n';
//			}
//			int num_paths_found = 0;
//			
//		}
//	}
//	//we terminate at the end of the read
//	//std::cout << "we terminate at the end of the read" << '\n';
//	return extensions;
//}
//
//std::vector<std::string> extend_right_het(std::string& ending_portion, int& i, int& k, int& error_threshold, int& het_threshold, int& unique_threshold, CKMCFile& file, int& read_number)
//{
//	std::vector<std::string> extensions;
//	//The ending portion should be of length k+(read_length-first_error_idx-k+1)
//	std::string kmer = ending_portion.substr(0, k);
//	std::vector<uint32_t> v;
//	int j = 0;
//	while (j < i)
//	{
//		bool found_path = false;
//		std::string new_kmer;
//		for (char const &c: "ACGT")
//		{
//			std::string adjacent_kmer = kmer.substr(1)+c;
//			file.GetCountersForRead(adjacent_kmer, v);
//			int current_type = get_type(v[0], error_threshold, het_threshold, unique_threshold);
//			//if the adjacent kmer is not an error
//			if (current_type > 0)
//			{
//				//we terminate at a branching path
//				if (found_path)
//				{
//					j--;
//					//remove the previous character we added
//					//extension = extension.substr(1);
//					extension.pop_back();
//					//extension = beginning_portion.substr(0, i-j) + extension;
//					//std::cout << extension << '\n';
//					//std::cout << ending_portion.substr(j+k) << '\n';
//					extension = extension + ending_portion.substr(j+k);
//					//std::cout << "we terminate at a branching path" << '\n';
//					return extension;
//				}
//				else
//				{
//					found_path = true;
//					j++;
//					new_kmer = adjacent_kmer;
//					//extension = c+extension;
//					extension = extension+c;
//				}
//			}
//		}
//		//we terminate at no paths
//		if (!found_path)
//		{
//			//extension = beginning_portion.substr(0, i-j) + extension;
//			//std::cout << extension << '\n';
//			//std::cout << ending_portion.substr(j+k) << '\n';
//			extension = extension + ending_portion.substr(j+k);
//			//std::cout << "we terminate at no paths" << '\n';
//			return extension;
//		}
//}

std::vector<std::string> extend_left_het(std::string& beginning_portion, int& i, int& k, int& error_threshold, int& het_threshold, int& unique_threshold, CKMCFile& file, int& read_number)
{
	std::vector<std::string> extensions;
	extensions.push_back(beginning_portion.substr(0, i));
	extensions.push_back(beginning_portion.substr(0, i));
	return extensions;
}

std::vector<std::string> extend_right_het(std::string& ending_portion, int& i, int& k, int& error_threshold, int& het_threshold, int& unique_threshold, CKMCFile& file, int& read_number)
{
	std::vector<std::string> extensions;
	extensions.push_back(ending_portion.substr(k));
	extensions.push_back(ending_portion.substr(k));
	return extensions;
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

void write_error_paths(bool& queue_broken, std::vector<std::string>& edited_error_portions, std::ofstream& erredits_output_file, std::ofstream& errpaths_output_file, std::string& edited_read, int& read_number, int& first_error_idx, int& last_error_idx, std::string& before_first_error_kmer, std::string& original_error_portion, std::string& after_last_error_kmer, CKMCFile& file)
{
	std::ofstream* errwrite_output_file;
	//we finished the search and presumably we have found one homozygous path or two heterozygous paths
	if ((!queue_broken) && ((edited_error_portions.size() == 1) or (edited_error_portions.size() == 2)))
	{
		//std::cout << "hi" << '\n';
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
		//std::cout << "read " << read_number << "does not have 1 or two error paths" << '\n';
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
				//last_error_idx = i-1;
				//first_nonerror_idx = i;
				//after_last_error_kmer = read.substr(i, k);
				//std::string original_error_portion = read.substr(0, i);
				//std::string original_error_block = original_error_portion + after_last_error_kmer;
				//std::string right_anchor_kmer = read.substr(i, k);
				//before_first_error_kmer = "";
				after_last_error_kmer = read.substr(i, k);
				int min_distance_of_path = 0;
				int max_distance_of_path = i;
				int max_nodes_to_search = 1000;
				bool queue_broken = false;
				std::vector<std::string> edited_error_portions = get_paths_left(file, error_threshold, het_threshold, unique_threshold, after_last_error_kmer, min_distance_of_path, max_distance_of_path, max_nodes_to_search, k, queue_broken);
				last_error_idx = i-1;
				first_nonerror_idx = i;
				//after_last_error_kmer = read.substr(i, k);
				std::string original_error_portion = read.substr(0, i);
				//std::string original_error_block = original_error_portion + after_last_error_kmer;
				write_error_paths(queue_broken, edited_error_portions, erredits_output_file, errpaths_output_file, edited_read, read_number, first_error_idx, last_error_idx, before_first_error_kmer, original_error_portion, after_last_error_kmer, file);
//				std::ofstream* err_output_file;
//				//we finished the search and presumably we have found one homozygous path or two heterozygous paths
//				if ((!queue_broken) && ((edited_error_portions.size() == 1) or (edited_error_portions.size() == 2)))
//				{
//					//std::cout << "hi" << '\n';
//					err_output_file = &erredits_output_file;
//					edited_read += edited_error_portions[0];
//				}
//				//there are no paths, or there are more than two paths, or the search wasn't finished.
//				//We are currently not editing.
//				else
//				{
//					//std::cout << "read " << read_number << "does not have 1 or two error paths" << '\n';
//					err_output_file = &errpaths_output_file;
//					edited_read += original_error_portion;
//				}
//				*err_output_file << ">read" << read_number << "_firsterrorkmer" << first_error_idx << "_lasterrorkmer" << last_error_idx << "_original" << '\n';
//				*err_output_file << " " << original_error_portion << " " << after_last_error_kmer << '\n';
//				std::vector<uint32_t> w;
//				file.GetCountersForRead(original_error_block, w);
//				for (int j=0; j < w.size(); j++)
//				{
//					*err_output_file << w.at(j) << " ";
//				}
//				*err_output_file << '\n';
//				for (int l = 0; l < edited_error_portions.size(); l++)
//				{
//					std::string edited_error_portion = edited_error_portions[l];
//					std::string edited_error_block = edited_error_portion + after_last_error_kmer;
//					*err_output_file << ">read" << read_number << "_firsterrorkmer" << first_error_idx << "_lasterrorkmer" << last_error_idx << "_edited" << l << '\n';
//					*err_output_file << " " << edited_error_portion << " " << after_last_error_kmer << '\n';
//					file.GetCountersForRead(edited_error_block, w);
//					for (int j=0; j < w.size(); j++)
//					{
//						*err_output_file << w.at(j) << " ";
//					}
//					*err_output_file << '\n';
//				}
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
				int max_distance_of_path = ceil(1.2 * number_of_error_kmers);
				int max_nodes_to_search = 1000;
				bool queue_broken;
				std::vector<std::string> edited_paths = get_paths(file, error_threshold, het_threshold, unique_threshold, before_first_error_kmer, after_last_error_kmer, min_distance_of_path, max_distance_of_path, max_nodes_to_search, k, queue_broken);
				last_error_idx = i-1;
				first_nonerror_idx = i;
				std::string original_error_portion = read.substr(first_error_idx+k-1, last_error_idx - first_error_idx + 2 - k);
				write_error_paths(queue_broken, edited_paths, erredits_output_file, errpaths_output_file, edited_read, read_number, first_error_idx, last_error_idx, before_first_error_kmer, original_error_portion, after_last_error_kmer, file);
//				//Only edit if we fully traveled all of the paths (with depth between min and max distance)
//				//and there exists only exactly one path that bridges the gap
//				if (!queue_broken && edited_paths.size() == 1)
//				{
//					std::string original_error_portion = read.substr(first_error_idx+k-1, last_error_idx - first_error_idx + 2 - k);
//					std::string original_error_block = before_first_error_kmer + original_error_portion + after_last_error_kmer;
//					erredits_output_file << ">read" << read_number << "_firsterrorkmer" << first_error_idx << "_lasterrorkmer" << last_error_idx << "_original" << '\n';
//					erredits_output_file << before_first_error_kmer << " " << original_error_portion << " " << after_last_error_kmer << '\n';
//					std::vector<uint32_t> w;
//					file.GetCountersForRead(original_error_block, w);
//					for (int j=0; j < w.size(); j++)
//					{
//						erredits_output_file << w.at(j) << " ";
//					}
//					erredits_output_file << '\n';
//					std::string edited_error_portion = edited_paths[0].substr(k-1);
//					std::string edited_error_block = before_first_error_kmer + edited_error_portion + after_last_error_kmer;
//					erredits_output_file << ">read" << read_number << "_firsterrorkmer" << first_error_idx << "_lasterrorkmer" << last_error_idx << "_edited" << '\n';
//					erredits_output_file << before_first_error_kmer << " " << edited_error_portion << " " << after_last_error_kmer << '\n';
//					file.GetCountersForRead(edited_error_block, w);
//					for (int j=0; j < w.size(); j++)
//					{
//						erredits_output_file << w.at(j) << " ";
//					}
//					erredits_output_file << '\n';
//					edited_read += edited_paths[0];
//				}
//				//Either we didn't fully travel all of the paths (with depth between min and max distance)
//				//or there does not exist only exactly one path that bridges that gap
//				else
//				{
//					std::string uneditable_error_portion = read.substr(first_error_idx, last_error_idx - first_error_idx + 1);
//					edited_read += uneditable_error_portion;
//					//Let's keep track of all the possible paths, just for bookkeeping.
//					//We only record those paths where the search completed and there are multiple.
//					if ((!queue_broken) && (edited_paths.size() > 1))
//					{
//						std::string portion;
//						errpaths_output_file << ">read" << read_number << "_firsterrorkmer" << first_error_idx << "_lasterrorkmer" << last_error_idx << "_original" << '\n';
//						if (uneditable_error_portion.length() >= k-1)
//						{
//							errpaths_output_file << before_first_error_kmer << " " << uneditable_error_portion.substr(k-1) << " " << after_last_error_kmer << '\n';
//							portion = before_first_error_kmer + uneditable_error_portion.substr(k-1) + after_last_error_kmer;
//						}
//						else
//						{
//							errpaths_output_file << before_first_error_kmer.front() << " " << uneditable_error_portion << " " << after_last_error_kmer << '\n';
//							portion = before_first_error_kmer.front() + uneditable_error_portion + after_last_error_kmer;
//						}
//						std::vector<uint32_t> w;
//						file.GetCountersForRead(portion, w);
//						for (int j=0; j < w.size(); j++)
//						{
//							errpaths_output_file << w.at(j) << " ";
//						}
//						errpaths_output_file << '\n';
//						for (int l = 0; l < edited_paths.size(); l++)
//						{
//							std::string edited_path = edited_paths[l];
//							errpaths_output_file << ">read" << read_number << "_firsterrorkmer" << first_error_idx << "_lasterrorkmer" << last_error_idx << "_path" << l << '\n';
//							if (edited_path.length() >= k-1)
//							{
//								errpaths_output_file << before_first_error_kmer << " " << edited_path.substr(k-1) << " " << after_last_error_kmer << '\n';
//								portion = before_first_error_kmer + edited_path.substr(k-1) + after_last_error_kmer;
//							}
//							else 
//							{
//								errpaths_output_file << before_first_error_kmer.front() << " " << edited_path << " " << after_last_error_kmer << '\n';
//								portion = before_first_error_kmer.front() + edited_path + after_last_error_kmer;
//							}
//							std::vector<uint32_t> w;
//							file.GetCountersForRead(portion, w);
//							for (int j=0; j < w.size(); j++)
//							{
//								errpaths_output_file << w.at(j) << " ";
//							}
//							errpaths_output_file << '\n';
//						}
//					}
//				}
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
		//std::cout << "hi" << '\n';
		//std::string left_anchor_kmer = read.substr(first_error_idx-1, k);
		after_last_error_kmer = "";
		int min_distance_of_path = 0;
		int max_distance_of_path = v.size()-first_error_idx;
		int max_nodes_to_search = 1000;
		bool queue_broken = false;
		//std::cout << "A" << '\n';
		std::vector<std::string> edited_error_portions = get_paths_right(file, error_threshold, het_threshold, unique_threshold, before_first_error_kmer, min_distance_of_path, max_distance_of_path, max_nodes_to_search, k, queue_broken);
		//std::cout << "B" << '\n';
		last_error_idx = v.size()-1;
		first_nonerror_idx = v.size();
		std::string original_error_portion = read.substr(first_error_idx+k-1);
		//std::cout << "C" << '\n';
		//std::string original_error_block = before_first_error_kmer + original_error_portion;
		write_error_paths(queue_broken, edited_error_portions, erredits_output_file, errpaths_output_file, edited_read, read_number, first_error_idx, last_error_idx, before_first_error_kmer, original_error_portion, after_last_error_kmer, file);
		//std::cout << "D" << '\n';
//		erredits_output_file << ">read" << read_number << "_firsterrorkmer" << first_error_idx << "_lasterrorkmer" << last_error_idx << "_original" << '\n';
//		erredits_output_file << before_first_error_kmer << " " << original_error_portion << " " << '\n';
//		std::vector<uint32_t> w;
//		file.GetCountersForRead(original_error_block, w);
//		for (int j=0; j < w.size(); j++)
//		{
//			erredits_output_file << w.at(j) << " ";
//		}
//		erredits_output_file << '\n';
//		std::string ending_portion = read.substr(first_error_idx-1);
//		int original_error_portion_length = v.size()-first_error_idx;
//		std::string edited_error_portion = extend_right_unique(ending_portion, original_error_portion_length, k, error_threshold, het_threshold, unique_threshold, file, read_number);
//		std::string edited_error_block = before_first_error_kmer + edited_error_portion;
//		erredits_output_file << ">read" << read_number << "_firsterrorkmer" << first_error_idx << "_lasterrorkmer" << last_error_idx << "_edited" << '\n';
//		erredits_output_file << before_first_error_kmer << " " << edited_error_portion << " " << '\n';
//		file.GetCountersForRead(edited_error_block, w);
//		for (int j=0; j < w.size(); j++)
//		{
//			erredits_output_file << w.at(j) << " ";
//		}
//		erredits_output_file << '\n';
//		edited_read += edited_error_portion;
	}
	if (previous_type > 0)
	{
		//std::cout << "hey" << '\n';
		//We have "left" the nonerror portion of the read
		first_error_idx = v.size();
		last_nonerror_idx = v.size()-1;
		std::string nonerror_portion = read.substr(first_nonerror_idx, last_nonerror_idx - first_nonerror_idx + k);
		edited_read += nonerror_portion;
	}
	edited_read += '\n';
	return edited_read;
}

std::string smooth_het (std::vector<uint32_t>& v, std::string& read, int& read_number, CKMCFile& file, std::ofstream& hetedits_output_file, std::ofstream& hetpaths_output_file)
{
	//initialize variables
	std::string smoothed_read;
	int k = 21;
	int error_threshold = 9; //8 for dmel, 9 for atha
	int het_threshold = 30; //21 for dmel, 30 for atha
	int unique_threshold = 72; //60 for dmel, 72 for atha

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
		//std::cout << i << '\n';
		int current_type = get_type(v[i], error_threshold, het_threshold, unique_threshold);
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
				last_nonhom_idx = i-1;
				first_hom_idx = i;
				after_last_nonhom_kmer = read.substr(i, k);
				std::string original_nonhom_portion = read.substr(0, i);
				std::string original_nonhom_block = original_nonhom_portion + after_last_nonhom_kmer;
				hetedits_output_file << ">read" << read_number << "_firstnonhomkmer" << first_nonhom_idx << "_lastnonhomkmer" << last_nonhom_idx << "_original" << '\n';
				hetedits_output_file << " " << original_nonhom_portion << " " << after_last_nonhom_kmer << '\n';
				std::vector<uint32_t> w;
				file.GetCountersForRead(original_nonhom_block, w);
				for (int j=0; j < w.size(); j++)
				{
					hetedits_output_file << w.at(j) << " ";
				}
				hetedits_output_file << '\n';
				std::string beginning_portion = read.substr(0, i+k);
				std::vector<std::string> smoothed_nonhom_portions = extend_left_het(beginning_portion, i, k, error_threshold, het_threshold, unique_threshold, file, read_number);
				std::string smoothed_nonhom_portion0 = smoothed_nonhom_portions[0];
				std::string smoothed_nonhom_portion1 = smoothed_nonhom_portions[1];
				file.GetCountersForRead(smoothed_nonhom_portions[0]+after_last_nonhom_kmer.substr(0, k-1), w);
				float average0 = std::accumulate(w.begin(), w.end(), 0.0)/w.size();
				file.GetCountersForRead(smoothed_nonhom_portions[1]+after_last_nonhom_kmer.substr(0, k-1), w);
				float average1 = std::accumulate(w.begin(), w.end(), 0.0)/w.size();
				std::string primary_smoothed_nonhom_portion;
				std::string alternate_smoothed_nonhom_portion;
				if (average0 > average1)
				{
					primary_smoothed_nonhom_portion = smoothed_nonhom_portion0;
					alternate_smoothed_nonhom_portion = smoothed_nonhom_portion1;
					smoothed_read += smoothed_nonhom_portions[0];
				}
				else
				{
					primary_smoothed_nonhom_portion = smoothed_nonhom_portion1;
					alternate_smoothed_nonhom_portion = smoothed_nonhom_portion0;
					smoothed_read += smoothed_nonhom_portions[1];
				}
				std::string primary_smoothed_nonhom_block = primary_smoothed_nonhom_portion + after_last_nonhom_kmer;
				std::string alternate_smoothed_nonhom_block = alternate_smoothed_nonhom_portion + after_last_nonhom_kmer;
				hetedits_output_file << ">read" << read_number << "_firstnonhomkmer" << first_nonhom_idx << "_lastnonhomkmer" << last_nonhom_idx << "_primary" << '\n';
				hetedits_output_file << " " << primary_smoothed_nonhom_portion << " " << after_last_nonhom_kmer << '\n';
				file.GetCountersForRead(primary_smoothed_nonhom_block, w);
				for (int j=0; j < w.size(); j++)
				{
					hetedits_output_file << w.at(j) << " ";
				}
				hetedits_output_file << '\n';
				hetedits_output_file << ">read" << read_number << "_firstnonhomkmer" << first_nonhom_idx << "_lastnonhomkmer" << last_nonhom_idx << "_alternate" << '\n';
				hetedits_output_file << " " << alternate_smoothed_nonhom_portion << " " << after_last_nonhom_kmer << '\n';
				file.GetCountersForRead(alternate_smoothed_nonhom_block, w);
				for (int j=0; j < w.size(); j++)
				{
					hetedits_output_file << w.at(j) << " ";
				}
				hetedits_output_file << '\n';
			}
			//if previous kmer was nonhom, we have left the nonhom block
			if (((previous_type == 0) || (previous_type == 1) || (previous_type == 3)) && !before_first_nonhom_kmer.empty())
			{
				int number_of_nonhom_kmers = i - first_nonhom_idx;
				//If the position of after_last_nonhom_kmer overlaps before_first_nonhom_kmer
				//we keep progressing as if nothing has happened, waiting to find another hom kmer
				if (number_of_nonhom_kmers < k)
				{
					current_type = previous_type;
					continue;
				}
				//get kmer that is right after the last nonhom kmer of block
				last_nonhom_idx = i-1;
				first_hom_idx = i;
				after_last_nonhom_kmer = read.substr(i, k);
				int min_distance_of_path = k;
        int max_distance_of_path = ceil(1.2 * number_of_nonhom_kmers);
				int max_nodes_to_search = 1000;
				bool queue_broken;
				std::vector<std::string> smoothed_paths = get_paths_het(file, error_threshold, het_threshold, unique_threshold, before_first_nonhom_kmer, after_last_nonhom_kmer, min_distance_of_path, max_distance_of_path, max_nodes_to_search, k, queue_broken);
				//Only smoothe if we fully traveled all of the paths (with depth between min and max distance)
				//and there exists only exactly two paths that bridge the gap
				if (!queue_broken && smoothed_paths.size() == 2)
				{
					std::string original_nonhom_portion = read.substr(first_nonhom_idx+k-1, last_nonhom_idx - first_nonhom_idx + 2 - k);
					std::string original_nonhom_block = before_first_nonhom_kmer + original_nonhom_portion + after_last_nonhom_kmer;
					hetedits_output_file << ">read" << read_number << "_firstnonhomkmer" << first_nonhom_idx << "_lastnonhomkmer" << last_nonhom_idx << "_original" << '\n';
					hetedits_output_file << before_first_nonhom_kmer << " " << original_nonhom_portion << " " << after_last_nonhom_kmer << '\n';
					std::vector<uint32_t> w;
					file.GetCountersForRead(original_nonhom_block, w);
					for (int j=0; j < w.size(); j++)
					{
						hetedits_output_file << w.at(j) << " ";
					}
					hetedits_output_file << '\n';
					std::string smoothed_nonhom_portion0 = smoothed_paths[0].substr(k-1);
					std::string smoothed_nonhom_portion1 = smoothed_paths[1].substr(k-1);
					file.GetCountersForRead(smoothed_paths[0]+after_last_nonhom_kmer.substr(0, k-1), w);
					float average0 = std::accumulate(w.begin(), w.end(), 0.0)/w.size();
					file.GetCountersForRead(smoothed_paths[1]+after_last_nonhom_kmer.substr(0, k-1), w);
					float average1 = std::accumulate(w.begin(), w.end(), 0.0)/w.size();
					std::string primary_smoothed_nonhom_portion;
					std::string alternate_smoothed_nonhom_portion;
					if (average0 > average1)
					{
						primary_smoothed_nonhom_portion = smoothed_nonhom_portion0;
						alternate_smoothed_nonhom_portion = smoothed_nonhom_portion1;
						smoothed_read += smoothed_paths[0];
					}
					else
					{
						primary_smoothed_nonhom_portion = smoothed_nonhom_portion1;
						alternate_smoothed_nonhom_portion = smoothed_nonhom_portion0;
						smoothed_read += smoothed_paths[1];
					}
					std::string primary_smoothed_nonhom_block = before_first_nonhom_kmer + primary_smoothed_nonhom_portion + after_last_nonhom_kmer;
					std::string alternate_smoothed_nonhom_block = before_first_nonhom_kmer + alternate_smoothed_nonhom_portion + after_last_nonhom_kmer;
					hetedits_output_file << ">read" << read_number << "_firstnonhomkmer" << first_nonhom_idx << "_lastnonhomkmer" << last_nonhom_idx << "_primary" << '\n';
					hetedits_output_file << before_first_nonhom_kmer << " " << primary_smoothed_nonhom_portion << " " << after_last_nonhom_kmer << '\n';
					file.GetCountersForRead(primary_smoothed_nonhom_block, w);
					for (int j=0; j < w.size(); j++)
					{
						hetedits_output_file << w.at(j) << " ";
					}
					hetedits_output_file << '\n';
					hetedits_output_file << ">read" << read_number << "_firstnonhomkmer" << first_nonhom_idx << "_lastnonhomkmer" << last_nonhom_idx << "_alternate" << '\n';
					hetedits_output_file << before_first_nonhom_kmer << " " << alternate_smoothed_nonhom_portion << " " << after_last_nonhom_kmer << '\n';
					file.GetCountersForRead(alternate_smoothed_nonhom_block, w);
					for (int j=0; j < w.size(); j++)
					{
						hetedits_output_file << w.at(j) << " ";
					}
					hetedits_output_file << '\n';
				}
				//Either we didn't fully travel all of the paths (with depth between min and max distance)
				//or there does not exist only exactly two paths that bridge the gap
				else
				{
					std::string unsmoothable_nonhom_portion = read.substr(first_nonhom_idx, last_nonhom_idx - first_nonhom_idx + 1);
					smoothed_read += unsmoothable_nonhom_portion;
					//Let's keep track of all the possible paths, just for bookkeeping.
					//We only record those paths where the search completed and there are more than 2 
					if ((!queue_broken) && (smoothed_paths.size() > 2))
					{
						std::string portion;
						hetpaths_output_file << ">read" << read_number << "_firstnonhomkmer" << first_nonhom_idx << "_lastnonhomkmer" << last_nonhom_idx << "_original" << '\n';
						if (unsmoothable_nonhom_portion.length() >= k-1)
						{
							hetpaths_output_file << before_first_nonhom_kmer << " " << unsmoothable_nonhom_portion.substr(k-1) << " " << after_last_nonhom_kmer << '\n';
							portion = before_first_nonhom_kmer + unsmoothable_nonhom_portion.substr(k-1) + after_last_nonhom_kmer;
						}
						else
						{
							hetpaths_output_file << before_first_nonhom_kmer.front() << " " << unsmoothable_nonhom_portion << " " << after_last_nonhom_kmer << '\n';
							portion = before_first_nonhom_kmer.front() + unsmoothable_nonhom_portion + after_last_nonhom_kmer;
						}
						std::vector<uint32_t> w;
						file.GetCountersForRead(portion, w);
						for (int j=0; j < w.size(); j++)
						{
							hetpaths_output_file << w.at(j) << " ";
						}
						hetpaths_output_file << '\n';
						for (int l = 0; l < smoothed_paths.size(); l++)
            {
              std::string smoothed_path = smoothed_paths[l];
							hetpaths_output_file << ">read" << read_number << "_firstnonhomkmer" << first_nonhom_idx << "_lastnonhomkmer" << last_nonhom_idx << "_path" << l << '\n';
							if (smoothed_path.length() >= k-1)
							{
								hetpaths_output_file << before_first_nonhom_kmer << " " << smoothed_path.substr(k-1) << " " << after_last_nonhom_kmer << '\n';
								portion = before_first_nonhom_kmer + smoothed_path.substr(k-1) + after_last_nonhom_kmer;
							}
							else 
							{
								hetpaths_output_file << before_first_nonhom_kmer.front() << " " << smoothed_path << " " << after_last_nonhom_kmer << '\n';
								portion = before_first_nonhom_kmer.front() + smoothed_path + after_last_nonhom_kmer;
							}
							std::vector<uint32_t> w;
							file.GetCountersForRead(portion, w);
							for (int j=0; j < w.size(); j++)
							{
								hetpaths_output_file << w.at(j) << " ";
							}
							hetpaths_output_file << '\n';
            }
					}
				}
			}
			//if previous kmer is hom, we are continuing a hom block
			if (previous_type == 2)
			{
				;
			}
			previous_type = current_type;
		}
	}
	//We have reached the end of the read, let's make sure we have added the last bit of the read
	if ((previous_type == 0) || (previous_type == 1) || (previous_type == 3))
	{
		//We have "left" the nonhom portion of the read
		last_nonhom_idx = v.size()-1;
		first_hom_idx = v.size();
		//If we have a homozygous kmer on left with which to anchor
		if (first_nonhom_idx > 0)
		{
			std::string original_nonhom_portion = read.substr(first_nonhom_idx+k-1);
			std::string original_nonhom_block = before_first_nonhom_kmer + original_nonhom_portion;
			hetedits_output_file << ">read" << read_number << "_firstnonhomkmer" << first_nonhom_idx << "_lastnonhomkmer" << last_nonhom_idx << "_original" << '\n';
			hetedits_output_file << before_first_nonhom_kmer << " " << original_nonhom_portion << " " << '\n';
			std::vector<uint32_t> w;
			file.GetCountersForRead(original_nonhom_block, w);
			for (int j=0; j < w.size(); j++)
			{
				hetedits_output_file << w.at(j) << " ";
			}
			hetedits_output_file << '\n';
			std::string ending_portion = read.substr(first_nonhom_idx-1);
			int original_nonhom_portion_length = v.size()-first_nonhom_idx;
			std::vector<std::string> smoothed_nonhom_portions = extend_right_het(ending_portion, original_nonhom_portion_length, k, error_threshold, het_threshold, unique_threshold, file, read_number);
			std::string smoothed_nonhom_portion0 = smoothed_nonhom_portions[0];
			std::string smoothed_nonhom_portion1 = smoothed_nonhom_portions[1];
			file.GetCountersForRead(before_first_nonhom_kmer.substr(1) + smoothed_nonhom_portion0, w);
			float average0 = std::accumulate(w.begin(), w.end(), 0.0)/w.size();
			file.GetCountersForRead(before_first_nonhom_kmer.substr(1) + smoothed_nonhom_portion1, w);
			float average1 = std::accumulate(w.begin(), w.end(), 0.0)/w.size();
			std::string primary_smoothed_nonhom_portion;
			std::string alternate_smoothed_nonhom_portion;
			if (average0 > average1)
			{
				primary_smoothed_nonhom_portion = smoothed_nonhom_portion0;
				alternate_smoothed_nonhom_portion = smoothed_nonhom_portion1;
				smoothed_read += smoothed_nonhom_portions[0];
			}
			else
			{
				primary_smoothed_nonhom_portion = smoothed_nonhom_portion1;
				alternate_smoothed_nonhom_portion = smoothed_nonhom_portion0;
				smoothed_read += smoothed_nonhom_portions[1];
			}
			std::string primary_smoothed_nonhom_block = before_first_nonhom_kmer + primary_smoothed_nonhom_portion;
			std::string alternate_smoothed_nonhom_block = before_first_nonhom_kmer + alternate_smoothed_nonhom_portion;
			hetedits_output_file << ">read" << read_number << "_firstnonhomkmer" << first_nonhom_idx << "_lastnonhomkmer" << last_nonhom_idx << "_primary" << '\n';
			hetedits_output_file << before_first_nonhom_kmer << " " << primary_smoothed_nonhom_portion << " " << '\n';
			file.GetCountersForRead(primary_smoothed_nonhom_block, w);
			for (int j=0; j < w.size(); j++)
			{
				hetedits_output_file << w.at(j) << " ";
			}
			hetedits_output_file << '\n';
			hetedits_output_file << ">read" << read_number << "_firstnonhomkmer" << first_nonhom_idx << "_lastnonhomkmer" << last_nonhom_idx << "_alternate" << '\n';
			hetedits_output_file << before_first_nonhom_kmer << " " << alternate_smoothed_nonhom_portion << " " << '\n';
			file.GetCountersForRead(alternate_smoothed_nonhom_block, w);
			for (int j=0; j < w.size(); j++)
			{
				hetedits_output_file << w.at(j) << " ";
			}
			hetedits_output_file << '\n';
		}
		//else, there is no homozygous kmer in the read
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
	}
	smoothed_read += '\n';
	return smoothed_read;
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

	//load reads
	int line_num = -1;
	//currently assuming fasta format
	while (getline(input_file, line))
	{
		line_num++;
		//if (line_num % 2 == 0)
		//if (line_num % 4 != 1)
		if (line_num % 4 == 0)
		{
			//write read header to err output file and het output file
			err_output_file << line << '\n';
			het_output_file << line << '\n';
		}
		//else
		if (line_num % 4 == 1)
		{
			std::string read = line;
			//std::cout << "size of read: " << read.size() << '\n';
			//int read_number = (line_num+1)/2;
			int read_number = (line_num+3)/4;
			//if (read_number%10000==0)
			//{
				std::cout << read_number << '\n';
			//}
			
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
				err_output_file << "N" << '\n';
				het_output_file << "N" << '\n';
				continue;
			}

			//remove errors from the read to get edited read
			std::string edited_read = remove_err(v, read, read_number, file, erredits_output_file, errpaths_output_file);
			//std::cout << "edited read size " << edited_read.size() << '\n';

			//write edited read to err_output_file
			err_output_file << edited_read;
			edited_read.pop_back();

			//get counters of kmers in edited read
			file.GetCountersForRead(edited_read, v);
			num_kmers = v.size();

			//check if too many repetitive kmers
			//if greater than 25% of kmers in the edited read are repetitive
			//discard the read to not waste time
			int num_rep = std::count_if(v.begin(), v.end(), IsRep);
			double rep_fraction = static_cast<double>(num_rep)/num_kmers;
			if (rep_fraction >= 0.25)
			{
				//std::cout << "read number " << read_number << " has " << rep_fraction << " percent repetitive kmers, discarding." << '\n';
				het_output_file << "N" << '\n';
				continue;
			}
			//smoothe het from the edited read to get smoothed read
			std::string smoothed_read = smooth_het(v, edited_read, read_number, file, hetedits_output_file, hetpaths_output_file);
			//write smoothed read to het_output_file
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
