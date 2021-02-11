#include "./../kmc_api/kmc_file.h"
#include "./../kmc_api/kmer_api.h"
#include "nc_utils.h"
#include "stdafx.h"
#include <algorithm>
#include <bits/stdc++.h>
#include <ctime>
#include <fstream>
#include <iostream>
#include <list>
#include <math.h>
#include <mutex>
#include <numeric>
#include <stdlib.h>
#include <thread>
#include <unistd.h>
#include <vector>

//smooth kmcdb reads.fastq errremoved.fasta erredits.fasta errpaths.fasta hetremoved.fasta hetedits.fasta hetpath1.fasta hetpath2.fasta hetpath3.fasta hetpath4.fasta hetpath5.fasta hetpath6.fasta error_threshold het_threshold unique_threshold anchor_threshold allowed_err_fraction allowed_rep_fraction max_nodes_to_search distance_multiplier strict > counts.txt

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
		if ((current_type > 0)) //&& (current_type < 3)) //Added the 2nd condition 
		{
		  adjacent_kmers.push_back(adjacent_kmer);
		}
	}
	return adjacent_kmers;
}

bool is_left_anchor(std::string& previous_kmer, int& previous_count, int& k, CKMCFile& file, int& error_threshold, int& het_threshold, int& unique_threshold, int& anchor_threshold)
{
	if ((previous_kmer.length() != k) || (previous_count > unique_threshold)) // > unique_threshold or anchor_threshold
	//if (previous_kmer.length() != k)
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

bool is_right_anchor(std::string& current_kmer, int& current_count, int& k, CKMCFile& file, int& error_threshold, int& het_threshold, int& unique_threshold, int& anchor_threshold)
{
	if ((current_kmer.length() != k) || (current_count > unique_threshold)) // > unique_threshold or anchor_threshold?
	//if (current_kmer.length() != k)
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

int get_type_het(int& previous_type, std::string& previous_kmer, std::string& current_kmer, int& previous_count, int& current_count, int& k, CKMCFile& file, int& error_threshold, int& het_threshold, int& unique_threshold, int& anchor_threshold, std::string& anchor_found)
{
	//If "previously" at the beginning of the read
	if (previous_type == -1)
	{
		//if current kmer is a right anchor
		if ((current_count > (previous_count + error_threshold)) && is_right_anchor(current_kmer, current_count, k, file, error_threshold, het_threshold, unique_threshold, anchor_threshold))
		{
			//std::cout << "current is right anchor" << '\n';
			anchor_found = "right";
			return 2;
		}
		//if current kmer is a left anchor
		else if (is_left_anchor(current_kmer, current_count, k, file, error_threshold, het_threshold, unique_threshold, anchor_threshold))
		{
			//std::cout << "current is left anchor" << '\n';
			anchor_found = "left";
			return 2;
		}
		//To keep it consistent with previous get_type, we need type to be -1 only before we start the read
		//So, in this case we will do what was done before, using the coverage and coverage thresholds to
		//determine whether the kmer is homozygous or not
		else if ((het_threshold < current_count) && (current_count <= unique_threshold))
		{
			anchor_found = "none";
			return 2;
		}
		else
		{
			anchor_found = "none";
			return 1;
		}
	}
	//If previously in a nonhomozygous region
	else if (previous_type == 1)
	{
		//If current kmer is a right anchor
		if ((current_count > (previous_count + error_threshold)) && is_right_anchor(current_kmer, current_count, k, file, error_threshold, het_threshold, unique_threshold, anchor_threshold))
		{
			//we have left the nonhom region
			//std::cout << "current is right anchor" << '\n';
			anchor_found = "right";
			return 2;
		}
		//If current kmer is a left anchor
		else if (is_left_anchor(current_kmer, current_count, k, file, error_threshold, het_threshold, unique_threshold, anchor_threshold))
		{
			//this is a weird case where we find two left anchors in a row before a right anchor
			std::cout << "we found two left anchors in a row" << '\n';
			anchor_found = "left";
			return 2;
		}
		//else current kmer is not an anchor
		else
		{
			//we continue the nonhom region
			anchor_found = "none";
			return 1;
		}
	}
	//If previously in a homozygous region
	//else if (previous_type == 2)
	else
	{
		//If current kmer is a right anchor
		//if ((current_count > (previous_count + error_threshold)) && is_right_anchor(current_kmer, current_count, k, file, error_threshold, het_threshold, unique_threshold, anchor_threshold))
		//{
			//this is a weird case where we find two right anchors in a row before a left anchor
			//we don't need to set anchor_found or return 2 because we really care
			//if current kmer is a left anchor, ending the homozygous region
			//std::cout << "we found two right anchors in a row" << '\n';
		//}
		//If current kmer is a left anchor
		if (is_left_anchor(current_kmer, current_count, k, file, error_threshold, het_threshold, unique_threshold, anchor_threshold))
		{
			//std::cout << "current is left anchor" << '\n';
			anchor_found = "left";
			return 2;
		}
		//If previous kmer is a left anchor
		else if ((current_count < (previous_count - error_threshold)) && is_left_anchor(previous_kmer, previous_count, k, file, error_threshold, het_threshold, unique_threshold, anchor_threshold))
		{
			//we have left the hom region
			//std::cout << "previous was left anchor" << '\n';
			anchor_found = "left";
			return 1;
		}
		//else previous kmer is not a left anchor (and current kmer is not a left anchor)
		else
		{
			//we continue the hom region
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
					//added this case for if the right and left anchors overlap
					if (!left_anchor_kmer.empty() && !right_anchor_kmer.empty() && (path.size() < 2*k))
					{
						int total_overlaps = 2*k - path.size();
						path.clear();
						for (int number_overlaps = 0; number_overlaps < total_overlaps; number_overlaps++)
						{
							path += "-";
						}
						paths.push_back(path);
					}
					else
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

int minDis(std::string s1, std::string s2, int n, int m, std::vector<std::vector<int>> &dp)
{
	// If any string is empty,
	// return the remaining characters of other string
	if (n==0)
	{
		return m;
	}
	if (m==0)
	{
		return n;
	}
	// To check if the recursive tree
	// for given n & m has already been executed	
	if (dp[n][m]!=-1)
	{
		return dp[n][m];
	}
	// If characters are equal, execute
	// recursive function for n-1, m-1
	
	if (s1[n-1] == s2[m-1])
	{
		if (dp[n-1][m-1] == -1)
		{
			return dp[n][m] = minDis(s1, s2, n-1, m-1, dp);
		}
		else
		{
			return dp[n][m] = dp[n-1][m-1];
		}
	}
	// If characters are nt equal, we need to
	// find the minimum cost out of all 3 operations.
	else
	{
		int m1, m2, m3; // temp variables
		
		if (dp[n-1][m]!=-1)
		{
			m1 = dp[n-1][m];
		}
		else
		{
			m1 = minDis(s1, s2, n-1, m, dp);
		}
		
		if (dp[n][m-1]!=-1)
		{
			m2 = dp[n][m-1];
		}
		else
		{
			m2 = minDis(s1, s2, n, m-1, dp);
		}
		
		if (dp[n-1][m-1]!=-1)
		{
			m3 = dp[n-1][m-1];
		}
		else
		{
			m3 = minDis(s1, s2, n-1, m-1, dp);
		}
		return dp[n][m] = 1 + std::min(m1, std::min(m2, m3));
	}
}

void write_error_paths(bool& queue_broken, std::vector<std::string>& edited_error_portions, std::ofstream& errpaths_queuecomplete0_numpaths0_output_file, std::ofstream& errpaths_queuecomplete0_numpaths1to2_output_file, std::ofstream& errpaths_queuecomplete0_numpaths3plus_output_file, std::ofstream& errpaths_queuecomplete1_numpaths0_output_file, std::ofstream& errpaths_queuecomplete1_numpaths1to2_output_file, std::ofstream& errpaths_queuecomplete1_numpaths3plus_output_file, std::ofstream& erredits_output_file, std::string& edited_read, int& read_number, int& first_error_idx, int& last_error_idx, std::string& before_first_error_kmer, std::string& original_error_portion, std::string& after_last_error_kmer, CKMCFile& file, std::mutex& outputfileMutex)
{
	std::string left_part;
	if (!before_first_error_kmer.empty())
	{
		left_part = before_first_error_kmer.substr(1);
	}
	else
	{
		left_part = before_first_error_kmer;
	}
	std::string right_part;
	if (!after_last_error_kmer.empty())
	{
		right_part = after_last_error_kmer.substr(0, after_last_error_kmer.length()-1);
	}
	else
	{
		right_part = after_last_error_kmer;
	}
	std::size_t found = original_error_portion.find("-");
	std::string original;
	if (found!=std::string::npos)
	{
		original = left_part + right_part.substr(original_error_portion.size());
	}
	else
	{
		original = left_part + original_error_portion + right_part;
	}
	auto is_less_than = [&](std::string edited_error_portion1, std::string edited_error_portion2)
	{
		int n = original.length();
		std::size_t found1 = edited_error_portion1.find("-");
		std::string portion1;
		if (found1!=std::string::npos)
		{
			portion1 = left_part + right_part.substr(edited_error_portion1.size());
		}
		else
		{
			portion1 = left_part + edited_error_portion1 + right_part;
		}
		int m = portion1.length();
		std::vector<std::vector<int>> dp(n+1, std::vector<int>(m+1, -1));
		int dist1 = minDis(original, portion1, n, m, dp);
		std::size_t found2 = edited_error_portion2.find("-");
		std::string portion2;
		if (found2!=std::string::npos)
		{
			portion2 = left_part + right_part.substr(edited_error_portion2.size());
		}
		else
		{
			portion2 = left_part + edited_error_portion2 + right_part;
		}
		m = portion2.length();
		std::vector<std::vector<int>> dp2(n+1, std::vector<int>(m+1, -1));
		int dist2 = minDis(original, portion2, n, m, dp2);
		return dist1 < dist2;
	};
	std::sort(edited_error_portions.begin(), edited_error_portions.end(), is_less_than);
	std::ofstream* errwrite_output_file;
	bool was_edited = false;
	//we finished the search and presumably we have found one homozygous path or two heterozygous paths
	//if ((!queue_broken) && ((edited_error_portions.size() == 1) or (edited_error_portions.size() == 2)))
	if (edited_error_portions.size() >= 1)
	{
		//errwrite_output_file = &erredits_output_file;
		found = edited_error_portions[0].find("-");
		if (found!=std::string::npos)
		{
			edited_read += left_part.substr(0, left_part.size() - edited_error_portions[0].size());
		}
		else
		{
			edited_read += left_part + edited_error_portions[0];
		}
		if (edited_error_portions[0] != original_error_portion)
		{
			was_edited = true;
		}
		//if (!before_first_error_kmer.empty())
		//{
		//	edited_read += before_first_error_kmer.substr(1) + edited_error_portions[0];
		//}
		//else
		//{
		//	edited_read += edited_error_portions[0];
		//}
	}
	//there are no paths, or there are more than two paths, or the search wasn't finished.
	//We are currently not editing.
	else
	{
		//errwrite_output_file = &errpaths_output_file;
		found = original_error_portion.find("-");
		if (found!=std::string::npos)
		{
			edited_read += left_part.substr(0, left_part.size() - original_error_portion.size());
		}
		else
		{
			edited_read += left_part + original_error_portion;
		}
		//if (!before_first_error_kmer.empty())
		//{
		//	edited_read += before_first_error_kmer.substr(1) + original_error_portion;
		//}
		//else
		//{
		//	edited_read += original_error_portion;
		//}
	}
	if ((!queue_broken) && ((edited_error_portions.size() == 1) or (edited_error_portions.size() == 2)))
	{
		errwrite_output_file = &errpaths_queuecomplete1_numpaths1to2_output_file;
	}
	if ((!queue_broken) && (edited_error_portions.size() < 1))
	{
		errwrite_output_file = &errpaths_queuecomplete1_numpaths0_output_file;
	}
	if ((!queue_broken) && (edited_error_portions.size() > 2))
	{
		errwrite_output_file = &errpaths_queuecomplete1_numpaths3plus_output_file;
	}
	if ((queue_broken) && ((edited_error_portions.size() == 1) or (edited_error_portions.size() == 2)))
	{
		errwrite_output_file = &errpaths_queuecomplete0_numpaths1to2_output_file;
	}
	if ((queue_broken) && (edited_error_portions.size() < 1))
	{
		errwrite_output_file = &errpaths_queuecomplete0_numpaths0_output_file;
	}
	if ((queue_broken) && (edited_error_portions.size() > 2))
	{
		errwrite_output_file = &errpaths_queuecomplete0_numpaths3plus_output_file;
	}
	outputfileMutex.lock();
	std::string original_error_block;
	found = original_error_portion.find("-");
	if (found!=std::string::npos)
	{
		original_error_block = before_first_error_kmer + after_last_error_kmer.substr(original_error_portion.size());
	}
	else
	{
		original_error_block = before_first_error_kmer + original_error_portion + after_last_error_kmer;
	}
	*errwrite_output_file << ">read" << read_number << "_firsterrorkmer" << first_error_idx << "_lasterrorkmer" << last_error_idx << "_original" << '\n';
	*errwrite_output_file << before_first_error_kmer << " " << original_error_portion << " " << after_last_error_kmer << '\n';
	if (was_edited)
	{
		erredits_output_file << ">read" << read_number << "_firsterrorkmer" << first_error_idx << "_lasterrorkmer" << last_error_idx << "_original" << '\n';
		erredits_output_file << before_first_error_kmer << " " << original_error_portion << " " << after_last_error_kmer << '\n';
	}
	std::vector<uint32_t> w;
	file.GetCountersForRead(original_error_block, w);
	for (int j=0; j < w.size(); j++)
	{
		*errwrite_output_file << w.at(j) << " ";
		if (was_edited)
		{
			erredits_output_file << w.at(j) << " ";
		}
	}
	*errwrite_output_file << '\n';
	if (was_edited)
	{
		erredits_output_file << '\n';
	}
	for (int l = 0; l < edited_error_portions.size(); l++)
	{
		std::string edited_error_portion = edited_error_portions[l];
		std::string edited_error_block;
		found = edited_error_portion.find("-");
		if (found!=std::string::npos)
		{
			edited_error_block = before_first_error_kmer + after_last_error_kmer.substr(edited_error_portion.size());
		}
		else
		{
			edited_error_block = before_first_error_kmer + edited_error_portion + after_last_error_kmer;
		}
		*errwrite_output_file << ">read" << read_number << "_firsterrorkmer" << first_error_idx << "_lasterrorkmer" << last_error_idx << "_path" << l << '\n';
		*errwrite_output_file << before_first_error_kmer << " " << edited_error_portion << " " << after_last_error_kmer << '\n';
		if ((l==0) && (was_edited))
		{
			erredits_output_file << ">read" << read_number << "_firsterrorkmer" << first_error_idx << "_lasterrorkmer" << last_error_idx << "_edited" << '\n';
			erredits_output_file << before_first_error_kmer << " " << edited_error_portion << " " << after_last_error_kmer << '\n';
		}
		file.GetCountersForRead(edited_error_block, w);
		for (int j=0; j < w.size(); j++)
		{
			*errwrite_output_file << w.at(j) << " ";
			if ((l==0) && (was_edited))
			{
				erredits_output_file << w.at(j) << " ";
			}
		}
		*errwrite_output_file << '\n';
		if ((l==0) && (was_edited))
		{
			erredits_output_file << '\n';
		}
	}
	outputfileMutex.unlock();
}

void write_nonhom_paths(bool& queue_broken, std::vector<std::string>& smoothed_nonhom_portions, std::ofstream& hetpaths_queuecomplete0_numpaths0to1_output_file, std::ofstream& hetpaths_queuecomplete0_numpaths2_output_file, std::ofstream& hetpaths_queuecomplete0_numpaths3plus_output_file, std::ofstream& hetpaths_queuecomplete1_numpaths0to1_output_file, std::ofstream& hetpaths_queuecomplete1_numpaths2_output_file, std::ofstream& hetpaths_queuecomplete1_numpaths3plus_output_file, std::ofstream& hetedits_output_file, std::string& smoothed_read, int& read_number, int& first_nonhom_idx, int& last_nonhom_idx, std::string& before_first_nonhom_kmer, std::string& original_nonhom_portion, std::string& after_last_nonhom_kmer, CKMCFile& file, int& strict, std::mutex& outputfileMutex)
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
	auto is_greater_than = [&](std::string smoothed_nonhom_portion1, std::string smoothed_nonhom_portion2)
	{
		std::vector<uint32_t> w;
		std::size_t found1 = smoothed_nonhom_portion1.find("-");
		//added this case to account for when anchors overlap
		if (found1!=std::string::npos)
		{
			file.GetCountersForRead(left_part + right_part.substr(smoothed_nonhom_portion1.size()), w);
		}
		else
		{
			file.GetCountersForRead(left_part + smoothed_nonhom_portion1 + right_part, w);
		}
		float average1 = std::accumulate(w.begin(), w.end(), 0.0)/w.size();
		std::size_t found2 = smoothed_nonhom_portion2.find("-");
		//added this case to account for when anchors overlap
		if (found2!=std::string::npos)
		{
			file.GetCountersForRead(left_part + right_part.substr(smoothed_nonhom_portion2.size()), w);
		}
		else
		{
			file.GetCountersForRead(left_part + smoothed_nonhom_portion2 + right_part, w);
		}
		float average2 = std::accumulate(w.begin(), w.end(), 0.0)/w.size();
		return average1 > average2;
	};
	std::sort(smoothed_nonhom_portions.begin(), smoothed_nonhom_portions.end(), is_greater_than);
	std::ofstream* hetwrite_output_file;
	//we finished the search and presumably we have found two heterozygous paths
	//actually, let's relax the !queue_broken requirement, jk restoring this requirement
	//and when strict==1 we only smooth when there are exactly two paths
	//and when strict==0 we smooth when there are at least two paths
	//if ((!queue_broken) && (smoothed_nonhom_portions.size() == 2))
	bool end_condition;
	bool was_smoothed = false;
	std::size_t found;
	if (strict==1)
	{
		end_condition = ((!queue_broken) && (smoothed_nonhom_portions.size() == 2));
	}
	else
	{
		end_condition = ((!queue_broken) && (smoothed_nonhom_portions.size() >= 2));
	}
	if (end_condition)
	{
		found = smoothed_nonhom_portions[0].find("-");
		if (found!=std::string::npos)
		{
			smoothed_read += left_part.substr(0, left_part.size() - smoothed_nonhom_portions[0].size());
		}
		else
		{
			smoothed_read += left_part + smoothed_nonhom_portions[0];
		}
		//We have checked that the condition for smoothing was met
		//Now we check whether the path chosen actually differs from the original
		//If so, then there was a true smooth and we add this to hetedits_output_file
		if (smoothed_nonhom_portions[0] != original_nonhom_portion)
		{
			was_smoothed = true;
		}
	}
	//there is not exactly two paths, or the search wasn't finished.
	//We are currently not smoothing.
	else
	{
		found = original_nonhom_portion.find("-");
		if (found!=std::string::npos)
		{
			smoothed_read += left_part.substr(0, left_part.size() - original_nonhom_portion.size());
		}
		else
		{
			smoothed_read += left_part + original_nonhom_portion;
		}
		//std::cout << left_part + original_nonhom_portion << '\n';
	}
	if ((!queue_broken) && (smoothed_nonhom_portions.size() == 2))
	{
		hetwrite_output_file = &hetpaths_queuecomplete1_numpaths2_output_file;
	}
	if ((!queue_broken) && (smoothed_nonhom_portions.size() < 2))
	{
		hetwrite_output_file = &hetpaths_queuecomplete1_numpaths0to1_output_file;
	}
	if ((!queue_broken) && (smoothed_nonhom_portions.size() > 2))
	{
		hetwrite_output_file = &hetpaths_queuecomplete1_numpaths3plus_output_file;
	}
	if ((queue_broken) && (smoothed_nonhom_portions.size() == 2))
	{
		hetwrite_output_file = &hetpaths_queuecomplete0_numpaths2_output_file;
	}
	if ((queue_broken) && (smoothed_nonhom_portions.size() < 2))
	{
		hetwrite_output_file = &hetpaths_queuecomplete0_numpaths0to1_output_file;
	}
	if ((queue_broken) && (smoothed_nonhom_portions.size() > 2))
	{
		hetwrite_output_file = &hetpaths_queuecomplete0_numpaths3plus_output_file;
	}
	outputfileMutex.lock();
	std::string original_nonhom_block;
	found = original_nonhom_portion.find("-");
	if (found!=std::string::npos)
	{
		original_nonhom_block = before_first_nonhom_kmer + after_last_nonhom_kmer.substr(original_nonhom_portion.size());
	}
	else
	{
		original_nonhom_block = before_first_nonhom_kmer + original_nonhom_portion + after_last_nonhom_kmer;
	}
	*hetwrite_output_file << ">read" << read_number << "_firstnonhomkmer" << first_nonhom_idx << "_lastnonhomkmer" << last_nonhom_idx << "_original" << '\n';
	*hetwrite_output_file << before_first_nonhom_kmer << " " << original_nonhom_portion << " " << after_last_nonhom_kmer << '\n';
	if (was_smoothed)
	{
		hetedits_output_file << ">read" << read_number << "_firstnonhomkmer" << first_nonhom_idx << "_lastnonhomkmer" << last_nonhom_idx << "_original" << '\n';
		hetedits_output_file << before_first_nonhom_kmer << " " << original_nonhom_portion << " " << after_last_nonhom_kmer << '\n';
	}
	std::vector<uint32_t> w;
	file.GetCountersForRead(original_nonhom_block, w);
	for (int j=0; j < w.size(); j++)
	{
		*hetwrite_output_file << w.at(j) << " ";
		if (was_smoothed)
		{
			hetedits_output_file << w.at(j) << " ";
		}
	}
	*hetwrite_output_file << '\n';
	if (was_smoothed)
	{
		hetedits_output_file << '\n';
	}
	for (int l = 0; l < smoothed_nonhom_portions.size(); l++)
	{
		std::string smoothed_nonhom_portion = smoothed_nonhom_portions[l];
		std::string smoothed_nonhom_block;
		found = smoothed_nonhom_portion.find("-");
		if (found!=std::string::npos)
		{
			smoothed_nonhom_block = before_first_nonhom_kmer + after_last_nonhom_kmer.substr(smoothed_nonhom_portion.size());
		}
		else
		{
			smoothed_nonhom_block = before_first_nonhom_kmer + smoothed_nonhom_portion + after_last_nonhom_kmer;
		}
		*hetwrite_output_file << ">read" << read_number << "_firstnonhomkmer" << first_nonhom_idx << "_lastnonhomkmer" << last_nonhom_idx << "_path" << l << '\n';
		*hetwrite_output_file << before_first_nonhom_kmer << " " << smoothed_nonhom_portion << " " << after_last_nonhom_kmer << '\n';
		if ((l==0) && (was_smoothed))
		{
			hetedits_output_file << ">read" << read_number << "_firstnonhomkmer" << first_nonhom_idx << "_lastnonhomkmer" << last_nonhom_idx << "_smoothed" << '\n';
			hetedits_output_file << before_first_nonhom_kmer << " " << smoothed_nonhom_portion << " " << after_last_nonhom_kmer << '\n';
		}
		file.GetCountersForRead(smoothed_nonhom_block, w);
		for (int j=0; j < w.size(); j++)
		{
			*hetwrite_output_file << w.at(j) << " ";
			if ((l==0) && (was_smoothed))
			{
				hetedits_output_file << w.at(j) << " ";
			}
		}
		*hetwrite_output_file << '\n';
		if ((l==0) && (was_smoothed))
		{
			hetedits_output_file << '\n';
		}
	}
	outputfileMutex.unlock();
}

std::string remove_err (std::vector<uint32_t>& v, std::string& read, int& read_number, CKMCFile& file, std::ofstream& errpaths_queuecomplete0_numpaths0_output_file, std::ofstream& errpaths_queuecomplete0_numpaths1to2_output_file, std::ofstream& errpaths_queuecomplete0_numpaths3plus_output_file, std::ofstream& errpaths_queuecomplete1_numpaths0_output_file, std::ofstream& errpaths_queuecomplete1_numpaths1to2_output_file, std::ofstream& errpaths_queuecomplete1_numpaths3plus_output_file, std::ofstream& erredits_output_file, int& error_threshold, int& het_threshold, int& unique_threshold, int& max_nodes_to_search, double& distance_multiplier, int& k, std::mutex& outputfileMutex)
{
	//initialize variables
	std::string edited_read;
	//int k = 21;

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
				bool queue_broken = false;
				std::vector<std::string> edited_error_portions = get_paths(file, error_threshold, het_threshold, unique_threshold, before_first_error_kmer, after_last_error_kmer, min_distance_of_path, max_distance_of_path, max_nodes_to_search, k, queue_broken);
				last_error_idx = i-1;
				first_nonerror_idx = i;
				std::string original_error_portion = read.substr(0, i);
				write_error_paths(queue_broken, edited_error_portions, errpaths_queuecomplete0_numpaths0_output_file, errpaths_queuecomplete0_numpaths1to2_output_file, errpaths_queuecomplete0_numpaths3plus_output_file, errpaths_queuecomplete1_numpaths0_output_file, errpaths_queuecomplete1_numpaths1to2_output_file, errpaths_queuecomplete1_numpaths3plus_output_file, erredits_output_file, edited_read, read_number, first_error_idx, last_error_idx, before_first_error_kmer, original_error_portion, after_last_error_kmer, file, outputfileMutex);
			}
			//if previous kmer was error, we have left the error block
			if (previous_type == 0 && !before_first_error_kmer.empty())
			{
				int number_of_error_kmers = i - first_error_idx;
				//If the position of after_last_error_kmer overlaps before_first_error_kmer
				//we keep progressing as if nothing has happened, waiting to find another non_error kmer
				//if (number_of_error_kmers < k)
				//{
				//	current_type = previous_type;
				//	continue;
				//}
				//get kmer that is right after the last error kmer of block
				after_last_error_kmer = read.substr(i, k);
				//int min_distance_of_path = k;
				int min_distance_of_path = std::min(number_of_error_kmers, k);
				int max_distance_of_path = ceil(distance_multiplier * number_of_error_kmers);
				bool queue_broken;
				std::vector<std::string> edited_error_portions = get_paths(file, error_threshold, het_threshold, unique_threshold, before_first_error_kmer, after_last_error_kmer, min_distance_of_path, max_distance_of_path, max_nodes_to_search, k, queue_broken);
				last_error_idx = i-1;
				first_nonerror_idx = i;
				std::string original_error_portion;
				//If the after_last_error_kmer overlaps before_first_error_kmer
				if (last_error_idx - first_error_idx + 2 - k < 0)
				{
					for (int number_overlaps=0; number_overlaps < first_error_idx - last_error_idx + k - 2; number_overlaps++)
					{
						original_error_portion += "-";
					}
				}
				else
				{
					original_error_portion = read.substr(first_error_idx+k-1, last_error_idx - first_error_idx + 2 - k);
				}
				write_error_paths(queue_broken, edited_error_portions, errpaths_queuecomplete0_numpaths0_output_file, errpaths_queuecomplete0_numpaths1to2_output_file, errpaths_queuecomplete0_numpaths3plus_output_file, errpaths_queuecomplete1_numpaths0_output_file, errpaths_queuecomplete1_numpaths1to2_output_file, errpaths_queuecomplete1_numpaths3plus_output_file, erredits_output_file, edited_read, read_number, first_error_idx, last_error_idx, before_first_error_kmer, original_error_portion, after_last_error_kmer, file, outputfileMutex);
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
		bool queue_broken = false;
		std::vector<std::string> edited_error_portions = get_paths(file, error_threshold, het_threshold, unique_threshold, before_first_error_kmer, after_last_error_kmer, min_distance_of_path, max_distance_of_path, max_nodes_to_search, k, queue_broken);
		last_error_idx = v.size()-1;
		first_nonerror_idx = v.size();
		std::string original_error_portion = read.substr(first_error_idx+k-1);
		write_error_paths(queue_broken, edited_error_portions, errpaths_queuecomplete0_numpaths0_output_file, errpaths_queuecomplete0_numpaths1to2_output_file, errpaths_queuecomplete0_numpaths3plus_output_file, errpaths_queuecomplete1_numpaths0_output_file, errpaths_queuecomplete1_numpaths1to2_output_file, errpaths_queuecomplete1_numpaths3plus_output_file, erredits_output_file, edited_read, read_number, first_error_idx, last_error_idx, before_first_error_kmer, original_error_portion, after_last_error_kmer, file, outputfileMutex);
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

std::string smooth_het (std::vector<uint32_t>& v, std::string& read, int& read_number, CKMCFile& file, std::ofstream& hetpaths_queuecomplete0_numpaths0to1_output_file, std::ofstream& hetpaths_queuecomplete0_numpaths2_output_file, std::ofstream& hetpaths_queuecomplete0_numpaths3plus_output_file, std::ofstream& hetpaths_queuecomplete1_numpaths0to1_output_file, std::ofstream& hetpaths_queuecomplete1_numpaths2_output_file, std::ofstream& hetpaths_queuecomplete1_numpaths3plus_output_file, std::ofstream& hetedits_output_file, int& error_threshold, int& het_threshold, int& unique_threshold, int& anchor_threshold, int& max_nodes_to_search, double& distance_multiplier, int& strict, int& k, std::mutex& outputfileMutex)
{
	//initialize variables
	std::string smoothed_read;
	//int k = 21;

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
		//std::cout << i << '\t' << v[i] << '\n';
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
		std::vector<uint32_t> current_count_vector;
		file.GetCountersForRead(current_kmer, current_count_vector);
		int current_count = current_count_vector[0];
		std::string anchor_found;
		int current_type = get_type_het(previous_type, previous_kmer, current_kmer, previous_count, current_count, k, file, error_threshold, het_threshold, unique_threshold, anchor_threshold, anchor_found);
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
				//std::cout << hom_portion << '\n';
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
				bool queue_broken = false;
				std::vector<std::string> smoothed_nonhom_portions = get_paths(file, error_threshold, het_threshold, unique_threshold, before_first_nonhom_kmer, after_last_nonhom_kmer, min_distance_of_path, max_distance_of_path, max_nodes_to_search, k, queue_broken);
				last_nonhom_idx = i-1;
				first_hom_idx = i;
				std::string original_nonhom_portion = read.substr(0, i);
				write_nonhom_paths(queue_broken, smoothed_nonhom_portions, hetpaths_queuecomplete0_numpaths0to1_output_file, hetpaths_queuecomplete0_numpaths2_output_file, hetpaths_queuecomplete0_numpaths3plus_output_file, hetpaths_queuecomplete1_numpaths0to1_output_file, hetpaths_queuecomplete1_numpaths2_output_file, hetpaths_queuecomplete1_numpaths3plus_output_file, hetedits_output_file, smoothed_read, read_number, first_nonhom_idx, last_nonhom_idx, before_first_nonhom_kmer, original_nonhom_portion, after_last_nonhom_kmer, file, strict, outputfileMutex);
			}
			//if previous kmer was nonhom, we have left the nonhom block
			if (((previous_type == 0) || (previous_type == 1) || (previous_type == 3)) && !before_first_nonhom_kmer.empty())
			{
				if (anchor_found == "left")
				{
					after_last_nonhom_kmer = read.substr(i, k);
					std::string hom_portion = read.substr(first_hom_idx+1, i-first_hom_idx-1);
					//std::cout << hom_portion << '\n';
					last_nonhom_idx = i-1;
					first_hom_idx = i;
					smoothed_read += hom_portion;
				}
				else
				{
				int number_of_nonhom_kmers = i - first_nonhom_idx;
				//If the position of after_last_nonhom_kmer overlaps before_first_nonhom_kmer
				//we keep progressing as if nothing has happened, waiting to find another hom kmer
				//if ((anchor_found == "right") && (number_of_nonhom_kmers < k))
				//{
				//	//std::cout << "left anchor and right anchor overlap" << '\n';
				//	current_type = previous_type;
				//	continue;
				//}
				//get kmer that is right after the last nonhom kmer of block
				after_last_nonhom_kmer = read.substr(i, k);
				//int min_distance_of_path = k;
				int min_distance_of_path = std::min(number_of_nonhom_kmers, k);
				int max_distance_of_path = ceil(distance_multiplier * number_of_nonhom_kmers);
				bool queue_broken;
				std::vector<std::string> smoothed_nonhom_portions = get_paths(file, error_threshold, het_threshold, unique_threshold, before_first_nonhom_kmer, after_last_nonhom_kmer, min_distance_of_path, max_distance_of_path, max_nodes_to_search, k, queue_broken);
				last_nonhom_idx = i-1;
				first_hom_idx = i;
				std::string original_nonhom_portion;
				//If the after_last_nonhom_kmer overlaps before_first_nonhom_kmer
				if (last_nonhom_idx - first_nonhom_idx + 2 - k < 0)
				{
					for (int number_overlaps=0; number_overlaps < first_nonhom_idx - last_nonhom_idx + k - 2; number_overlaps++)
					{
						original_nonhom_portion += "-";
					}
				}
				else
				{
					original_nonhom_portion = read.substr(first_nonhom_idx+k-1, last_nonhom_idx - first_nonhom_idx + 2 - k);
				}
				write_nonhom_paths(queue_broken, smoothed_nonhom_portions, hetpaths_queuecomplete0_numpaths0to1_output_file, hetpaths_queuecomplete0_numpaths2_output_file, hetpaths_queuecomplete0_numpaths3plus_output_file, hetpaths_queuecomplete1_numpaths0to1_output_file, hetpaths_queuecomplete1_numpaths2_output_file, hetpaths_queuecomplete1_numpaths3plus_output_file, hetedits_output_file, smoothed_read, read_number, first_nonhom_idx, last_nonhom_idx, before_first_nonhom_kmer, original_nonhom_portion, after_last_nonhom_kmer, file, strict, outputfileMutex);
				}
			}
			//if previous kmer is hom, we are continuing a hom block
			if (previous_type == 2)
			{
				if (anchor_found == "left")
				{
					std::string hom_portion = read.substr(first_hom_idx, i-first_hom_idx);
					smoothed_read += hom_portion;
					//std::cout << hom_portion << '\n';
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
			bool queue_broken = false;
			std::vector<std::string> smoothed_nonhom_portions = get_paths(file, error_threshold, het_threshold, unique_threshold, before_first_nonhom_kmer, after_last_nonhom_kmer, min_distance_of_path, max_distance_of_path, max_nodes_to_search, k, queue_broken);
			last_nonhom_idx = v.size()-1;
			first_hom_idx = v.size();
			std::string original_nonhom_portion = read.substr(first_nonhom_idx+k-1);
			write_nonhom_paths(queue_broken, smoothed_nonhom_portions, hetpaths_queuecomplete0_numpaths0to1_output_file, hetpaths_queuecomplete0_numpaths2_output_file, hetpaths_queuecomplete0_numpaths3plus_output_file, hetpaths_queuecomplete1_numpaths0to1_output_file, hetpaths_queuecomplete1_numpaths2_output_file, hetpaths_queuecomplete1_numpaths3plus_output_file, hetedits_output_file, smoothed_read, read_number, first_nonhom_idx, last_nonhom_idx, before_first_nonhom_kmer, original_nonhom_portion, after_last_nonhom_kmer, file, strict, outputfileMutex);
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
		//std::cout << hom_portion << '\n';
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

void processRead(std::ifstream& input_file, int& line_num, int& num_lines_per_read, CKMCFile& file, int& error_threshold, double& allowed_err_fraction, std::ofstream& err_output_file, std::ofstream& errpaths_queuecomplete0_numpaths0_output_file, std::ofstream& errpaths_queuecomplete0_numpaths1to2_output_file, std::ofstream& errpaths_queuecomplete0_numpaths3plus_output_file, std::ofstream& errpaths_queuecomplete1_numpaths0_output_file, std::ofstream& errpaths_queuecomplete1_numpaths1to2_output_file, std::ofstream& errpaths_queuecomplete1_numpaths3plus_output_file, std::ofstream& erredits_output_file, int& het_threshold, int& unique_threshold, int& max_nodes_to_search, double& distance_multiplier, int& k, int& polish, double& allowed_rep_fraction, std::ofstream& het_output_file, std::ofstream& hetpaths_queuecomplete0_numpaths0to1_output_file, std::ofstream& hetpaths_queuecomplete0_numpaths2_output_file, std::ofstream& hetpaths_queuecomplete0_numpaths3plus_output_file, std::ofstream& hetpaths_queuecomplete1_numpaths0to1_output_file, std::ofstream& hetpaths_queuecomplete1_numpaths2_output_file, std::ofstream& hetpaths_queuecomplete1_numpaths3plus_output_file, std::ofstream& hetedits_output_file, int& anchor_threshold, int& strict, std::mutex& inputfileMutex, std::mutex& outputfileMutex, std::mutex& coutMutex, int& verbose)
{
	std::string line;
	int thread_line_num;
	std::string header;
	inputfileMutex.lock();
	while (getline(input_file, line))
	{
		line_num++;
		thread_line_num = line_num;
		//if on line with the header
		if (thread_line_num % num_lines_per_read == 0)
		{
			//save header, but make sure it is a fasta header, not a fastq header
			header = ">" + line.substr(1);
		}
		//if on line with the read
		if (thread_line_num % num_lines_per_read == 1)
		{
			inputfileMutex.unlock();
			std::string read = line;
			int read_number = (thread_line_num+num_lines_per_read-1)/num_lines_per_read - 1;
			std::time_t thread_result = std::time(nullptr);
			if (verbose == 1)
			{
				coutMutex.lock();
				std::cout << "Analyzing read/contig/scaffold number " << read_number << " at " << std::asctime(std::localtime(&thread_result)) << '\n';
				coutMutex.unlock();
			}
			if (read_number%10000==0)
			{
				//std::cout << read_number << '\n';
			}
			
			//get counters of kmers in read
			std::vector<uint32_t> v;
			file.GetCountersForRead(read, v);
			int num_kmers = v.size();
			
			//check if too many error kmers
			//if greater than allowed_err_fraction of kmers in the read are errors
			//discard the read to not waste time
			int num_err = 0;
			for (int i = 0; i < v.size(); i++)
			{
				//std::cout << v[i] << ' ';
				if (v[i] <= error_threshold)
				{
					num_err++;
				}
			}
			//std::cout << '\n';
			//continue;
			double err_fraction = static_cast<double>(num_err)/num_kmers;
			if (err_fraction > allowed_err_fraction)
			{
				//std::cout << "read number " << read_number << " has " << err_fraction << " percent error kmers, discarding." << '\n';
				continue;
			}

			//remove errors from the read to get edited read
			std::string edited_read = remove_err(v, read, read_number, file, errpaths_queuecomplete0_numpaths0_output_file, errpaths_queuecomplete0_numpaths1to2_output_file, errpaths_queuecomplete0_numpaths3plus_output_file, errpaths_queuecomplete1_numpaths0_output_file, errpaths_queuecomplete1_numpaths1to2_output_file, errpaths_queuecomplete1_numpaths3plus_output_file, erredits_output_file, error_threshold, het_threshold, unique_threshold, max_nodes_to_search, distance_multiplier, k, outputfileMutex);

			//write header and edited read to err_output_file
			outputfileMutex.lock();
			err_output_file << header << '\n';
			err_output_file << edited_read;
			outputfileMutex.unlock();
			edited_read.pop_back();
			if (polish == 1)
			{
				continue;
			}

			//get counters of kmers in edited read
			file.GetCountersForRead(edited_read, v);
			num_kmers = v.size();

			//check if too many repetitive kmers
			//if greater than allowed_rep_fraction of kmers in the edited read are repetitive
			//discard the read to not waste time
			int num_rep = 0;
			for (int i = 0; i < v.size(); i++)
			{
				//std::cout << v[i] << ' ';
				if (v[i] > unique_threshold)
				{
					num_rep++;
				}
			}
			//std::cout << '\n';
			//continue;
			double rep_fraction = static_cast<double>(num_rep)/num_kmers;
			if (rep_fraction > allowed_rep_fraction)
			{
				//std::cout << "read number " << read_number << " has " << rep_fraction << " percent repetitive kmers, discarding." << '\n';
				outputfileMutex.lock();
				het_output_file << header << '\n';
				het_output_file << edited_read << '\n';
				outputfileMutex.unlock();
				continue;
			}

			//smoothe het from the edited read to get smoothed read
			std::string smoothed_read = smooth_het(v, edited_read, read_number, file, hetpaths_queuecomplete0_numpaths0to1_output_file, hetpaths_queuecomplete0_numpaths2_output_file, hetpaths_queuecomplete0_numpaths3plus_output_file, hetpaths_queuecomplete1_numpaths0to1_output_file, hetpaths_queuecomplete1_numpaths2_output_file, hetpaths_queuecomplete1_numpaths3plus_output_file, hetedits_output_file, error_threshold, het_threshold, unique_threshold, anchor_threshold, max_nodes_to_search, distance_multiplier, strict, k, outputfileMutex);

			//write header and smoothed read to het_output_file
			outputfileMutex.lock();
			het_output_file << header << '\n';
			het_output_file << smoothed_read;
			outputfileMutex.unlock();
			inputfileMutex.lock();
		}
	}
	inputfileMutex.unlock();
}

int main(int argc, char* argv[])
{
	//parse arguments
	int c;
	std::ifstream input_file;
	int num_lines_per_read;
	std::string kmcdb;
	std::string outdir;
	std::ofstream err_output_file;
	std::ofstream errpaths_queuecomplete0_numpaths0_output_file;
	std::ofstream errpaths_queuecomplete0_numpaths1to2_output_file;
	std::ofstream errpaths_queuecomplete0_numpaths3plus_output_file;
	std::ofstream errpaths_queuecomplete1_numpaths0_output_file;
	std::ofstream errpaths_queuecomplete1_numpaths1to2_output_file;
	std::ofstream errpaths_queuecomplete1_numpaths3plus_output_file;
	std::ofstream erredits_output_file;
	std::ofstream het_output_file;
	std::ofstream hetpaths_queuecomplete0_numpaths0to1_output_file;
	std::ofstream hetpaths_queuecomplete0_numpaths2_output_file;
	std::ofstream hetpaths_queuecomplete0_numpaths3plus_output_file;
	std::ofstream hetpaths_queuecomplete1_numpaths0to1_output_file;
	std::ofstream hetpaths_queuecomplete1_numpaths2_output_file;
	std::ofstream hetpaths_queuecomplete1_numpaths3plus_output_file;
	std::ofstream hetedits_output_file;
	int k = 0;
	double l = 0;
	int error_threshold = 0;
	int het_threshold = 0;
	int unique_threshold = 0;
	int anchor_threshold = 0;
	double allowed_err_fraction = 1;
	double allowed_rep_fraction = 1;
	int max_nodes_to_search = 1000;
	double distance_multiplier = 1.2;
	int strict = 1;
	int polish = 0;
	int num_threads = 1;
	int verbose = 0;
	
	while ((c = getopt(argc, argv, "hi:j:o:k:l:g:m:a:u:t:pn:d:e:r:s:v")) != -1)
	{
		switch (c)
		{
			case 'h':
				fprintf(stderr, "Usage: %s -i input.fa/fq -j kmcdb -o outdir -k kmersize -l lambda [-t num_threads (default 1)] [-p (run in polish mode, i.e. run only error correction and not het smoothing)] [-n max_nodes_to_search (default 1000)] [-d distance_multiplier (default 1.2)] [-e allowed_err_fraction (default 1.0)] [-r allowed_rep_fraction (default 1.0)] [-s strict (0 or 1, default 1)] [-v (run in verbose mode, i.e. print time stamps for analyses)]\n", argv[0]);
				exit(EXIT_FAILURE);
			case 'i':
				//is the input fasta or fastq?
				//note: the output will be fasta format, since quality values
				//will not match once the read is edited and smoothed
				if ((getFileExt(optarg) == "fasta") || (getFileExt(optarg) == "fa"))
				{
					num_lines_per_read = 2;
				}
				else if ((getFileExt(optarg) == "fastq") || (getFileExt(optarg) == "fq"))
				{
					num_lines_per_read = 4;
				}
				else
				{
					fprintf(stderr, "Input filename must end in .fa .fasta .fq or .fastq.\n");
					exit(EXIT_FAILURE);
				}
				input_file.open(optarg);
				if (!input_file.is_open())
				{
					fprintf(stderr, "Please ensure %s exists.\n", optarg);
					exit(EXIT_FAILURE);
				}
				break;
			case 'j':
				kmcdb = optarg;
				break;
			case 'o':
				outdir = optarg;
				err_output_file.open(outdir+"/errremoved.fasta");
				errpaths_queuecomplete0_numpaths0_output_file.open(outdir+"/errpaths_queuecomplete0_numpaths0.fasta");
				errpaths_queuecomplete0_numpaths1to2_output_file.open(outdir+"/errpaths_queuecomplete0_numpaths1to2.fasta");
				errpaths_queuecomplete0_numpaths3plus_output_file.open(outdir+"/errpaths_queuecomplete0_numpaths3plus.fasta");
				errpaths_queuecomplete1_numpaths0_output_file.open(outdir+"/errpaths_queuecomplete1_numpaths0.fasta");
				errpaths_queuecomplete1_numpaths1to2_output_file.open(outdir+"/errpaths_queuecomplete1_numpaths1to2.fasta");
				errpaths_queuecomplete1_numpaths3plus_output_file.open(outdir+"/errpaths_queuecomplete1_numpaths3plus.fasta");
				erredits_output_file.open(outdir+"/erredits.fasta");
				het_output_file.open(outdir+"/hetremoved.fasta");
				hetpaths_queuecomplete0_numpaths0to1_output_file.open(outdir+"/hetpaths_queuecomplete0_numpaths0to1.fasta");
				hetpaths_queuecomplete0_numpaths2_output_file.open(outdir+"/hetpaths_queuecomplete0_numpaths2.fasta");
				hetpaths_queuecomplete0_numpaths3plus_output_file.open(outdir+"/hetpaths_queuecomplete0_numpaths3plus.fasta");
				hetpaths_queuecomplete1_numpaths0to1_output_file.open(outdir+"/hetpaths_queuecomplete1_numpaths0to1.fasta");
				hetpaths_queuecomplete1_numpaths2_output_file.open(outdir+"/hetpaths_queuecomplete1_numpaths2.fasta");
				hetpaths_queuecomplete1_numpaths3plus_output_file.open(outdir+"/hetpaths_queuecomplete1_numpaths3plus.fasta");
				hetedits_output_file.open(outdir+"/hetedits.fasta");
				break;
			case 'k':
				k = atoi(optarg);
				break;
			case 'l':
				//only set values based on lambda if not already set by hidden parameters g, m, a, or u
				l = std::stod(optarg);
				if (error_threshold == 0)
				{
					error_threshold = round(ceil(0.33 * l));
				}
				if (het_threshold == 0)
				{
					het_threshold = round(ceil(1.5 * l));
				}
				if (unique_threshold == 0)
				{
					unique_threshold = round(ceil(3.5 * l));
				}
				if (anchor_threshold == 0)
				{
					anchor_threshold = round(ceil(2.5 * l));
				}
				break;
			case 'g':
				error_threshold = atoi(optarg);
				break;
			case 'm':
				het_threshold = atoi(optarg);
				break;
			case 'a':
				anchor_threshold = atoi(optarg);
				break;
			case 'u':
				unique_threshold = atoi(optarg);
				break;
			case 't':
				num_threads = atoi(optarg);
				break;
			case 'p':
				polish = 1;
				break;
			case 'n':
				max_nodes_to_search = atoi(optarg);
				break;
			case 'd':
				distance_multiplier = std::stod(optarg);
				break;
			case 'e':
				allowed_err_fraction = std::stod(optarg);
				break;
			case 'r':
				allowed_rep_fraction = std::stod(optarg);
				break;
			case 's':
				strict = atoi(optarg);
				break;
			case 'v':
				verbose = 1;
				break;
			case '?':
				fprintf(stderr, "Option -%c is invalid or requires an argument.\n", optopt);
			default:
				fprintf(stderr, "Usage: %s -i input.fa/fq -j kmcdb -o outdir -k kmersize -l lambda [-t num_threads (default 1)] [-p (run in polish mode, i.e. run only error correction and not het smoothing)] [-n max_nodes_to_search (default 1000)] [-d distance_multiplier (default 1.2)] [-e allowed_err_fraction (default 1.0)] [-r allowed_rep_fraction (default 1.0)] [-s strict (0 or 1, default 1)] [-v (run in verbose mode, i.e. print time stamps for analyses)]\n", argv[0]);
				exit(EXIT_FAILURE);
		}
	}
	//check that required arguments are given
	if (!input_file.is_open())
	{
		fprintf(stderr, "Please provide input file with -i argument.\n");
		exit(EXIT_FAILURE);
	}
	if (kmcdb.empty())
	{
		fprintf(stderr, "Please provide kmcdb with -j argument.\n");
		exit(EXIT_FAILURE);
	}
	if (outdir.empty())
	{
		fprintf(stderr, "Please provide output directory with -o argument.\n");
		exit(EXIT_FAILURE);
	}
	if (k==0)
	{
		fprintf(stderr, "Please provide kmer size with -k argument.\n");
		exit(EXIT_FAILURE);
	}
	if ((error_threshold == 0) || (het_threshold == 0) || (anchor_threshold == 0) || (unique_threshold == 0))
	{
		fprintf(stderr, "Please provide average kmer coverage with -l argument. (Or specify the thresholds manually).\n");
		exit(EXIT_FAILURE);
	}
	//std::cout << error_threshold << " " << het_threshold << " " << anchor_threshold << " " << unique_threshold << "\n";
	//load KMC database
	CKMCFile file;
	std::time_t result = std::time(nullptr);
	std::cout << "Loading KMC database at " << std::asctime(std::localtime(&result)) << '\n';
	file.OpenForRA(kmcdb);
	result = std::time(nullptr);
	std::cout << "KMC database loaded at " << std::asctime(std::localtime(&result)) << '\n';
	
	//global variables
	int line_num = -1;
	
	//make locks for the threads
	std::mutex inputfileMutex;
	std::mutex outputfileMutex;
	std::mutex coutMutex;
	
	//process reads using threads
	std::vector<std::thread> threads;
	for (int i = 0; i < num_threads; i++)
	{
		threads.push_back(std::thread(processRead, std::ref(input_file), std::ref(line_num), std::ref(num_lines_per_read), std::ref(file), std::ref(error_threshold), std::ref(allowed_err_fraction), std::ref(err_output_file), std::ref(errpaths_queuecomplete0_numpaths0_output_file), std::ref(errpaths_queuecomplete0_numpaths1to2_output_file), std::ref(errpaths_queuecomplete0_numpaths3plus_output_file), std::ref(errpaths_queuecomplete1_numpaths0_output_file), std::ref(errpaths_queuecomplete1_numpaths1to2_output_file), std::ref(errpaths_queuecomplete1_numpaths3plus_output_file), std::ref(erredits_output_file), std::ref(het_threshold), std::ref(unique_threshold), std::ref(max_nodes_to_search), std::ref(distance_multiplier), std::ref(k), std::ref(polish), std::ref(allowed_rep_fraction), std::ref(het_output_file), std::ref(hetpaths_queuecomplete0_numpaths0to1_output_file), std::ref(hetpaths_queuecomplete0_numpaths2_output_file), std::ref(hetpaths_queuecomplete0_numpaths3plus_output_file), std::ref(hetpaths_queuecomplete1_numpaths0to1_output_file), std::ref(hetpaths_queuecomplete1_numpaths2_output_file), std::ref(hetpaths_queuecomplete1_numpaths3plus_output_file), std::ref(hetedits_output_file), std::ref(anchor_threshold), std::ref(strict), std::ref(inputfileMutex), std::ref(outputfileMutex), std::ref(coutMutex), std::ref(verbose)));
	}
	for (auto &th : threads)
	{
		th.join();
	}
	
	//close files
	input_file.close();
	err_output_file.close();
	errpaths_queuecomplete0_numpaths0_output_file.close();
	errpaths_queuecomplete0_numpaths1to2_output_file.close();
	errpaths_queuecomplete0_numpaths3plus_output_file.close();
	errpaths_queuecomplete1_numpaths0_output_file.close();
	errpaths_queuecomplete1_numpaths1to2_output_file.close();
	errpaths_queuecomplete1_numpaths3plus_output_file.close();
	erredits_output_file.close();
	het_output_file.close();
	hetpaths_queuecomplete0_numpaths0to1_output_file.close();
	hetpaths_queuecomplete0_numpaths2_output_file.close();
	hetpaths_queuecomplete0_numpaths3plus_output_file.close();
	hetpaths_queuecomplete1_numpaths0to1_output_file.close();
	hetpaths_queuecomplete1_numpaths2_output_file.close();
	hetpaths_queuecomplete1_numpaths3plus_output_file.close();
	hetedits_output_file.close();
	return 0;
}
