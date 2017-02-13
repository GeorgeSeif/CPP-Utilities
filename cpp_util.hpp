#pragma once

#pragma warning(disable : 4996) //_CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <limits>  
#include <math.h>
#include <vector>
#include <numeric> 
#include <algorithm>
#include <iterator>
#include <chrono>  
#include <ctime>   
#include <sstream> 
#include <iomanip> 
#include <string> 
#include <utility>  
#include <complex>
#include <iterator>
#include <list>

using namespace std::complex_literals;

namespace constants
{
	static const double PI = 3.14159265358979323; // The constant PI
	static const double e = 2.71828182845904523; // Euler's constant
	static const double light = 299792458; // Meters per second

	static const char single_space = ' ';
	static const char comma = ',';
	static const char back_slash = '/';
	static const char period = '.';
	static const char equals_sign = '=';
}

namespace cpp_help
{
	// Logger function
	inline void LOG(std::string given_string)
	{
		// Get the current time
		time_t t = time(NULL);

		// Convert time_t to tm
		struct tm* my_time = std::localtime(&t);

		std::cout << "[ " << my_time->tm_mday << "/" << my_time->tm_mon << "/" << 1900 + my_time->tm_year << " "
			<< my_time->tm_hour << ":" << my_time->tm_min << ":" << my_time->tm_sec << " ]  "
			<< given_string << std::endl;
	}

	// Print an array
	template <class T>
	T print_arr(T arr)
	{
		for (int i = 0; i < arr.size(); i++)
		{
			std::cout << arr[i] << std::endl;
		}

		return static_cast<T>(NULL);
	}

}

namespace arr_utils
{
	// Insertion sort sorting algorithm
	template <class T>
	T insertion_sort(T arr)
	{
		for (int j = 1; j < arr.size(); j++)
		{
			// Make sure the key has the same type as the array
			typename T::value_type key = arr[j];

			// Insert arr[j] into the sorted sequence arr[1...j-1]
			int i = j - 1;
			while (i >= 0 && arr[i] > key)
			{
				arr[i + 1] = arr[i];
				i = i - 1;
			}

			arr[i + 1] = key;
		}

		return (arr);
	}

	// Merge to sorted arrays
	// arr is an array and p, q, and r are indices into the array such that p <= q < r.
	template <class T>
	T merge(T arr, int p, int q, int r)
	{
		int n_1 = q - p + 1;
		int n_2 = r - q;

		// First seperate the single array into two sub arrays 
		T left_arr(n_1 + 1);
		T right_arr(n_2 + 1);

		for (int i = 0; i < n_1; i++)
		{
			left_arr[i] = arr[p + i];
		}

		for (int j = 0; j < n_2; j++)
		{
			right_arr[j] = arr[q + 1 + j];
		}

		// Using extra sentinal value to avoid having to check whether either array is empty in each basic step
		left_arr[n_1] = std::numeric_limits<T::value_type>::max();
		right_arr[n_2] = std::numeric_limits<T::value_type>::max();

		// Merge the two arrays together
		int i = 0;
		int j = 0;

		for (int k = p; k <= r; k++)
		{
			if (left_arr[i] <= right_arr[j])
			{
				arr[k] = left_arr[i];
				i++;
			}
			else
			{
				arr[k] = right_arr[j];
				j++;
			}
		}

		return arr;
	}

	// Merge sort sorting algorithm
	template <class T>
	T merge_sort(T arr, int p, int r)
	{
		if (p < r)
		{
			int q = floor((p + r) / 2);
			arr = merge_sort(arr, p, q);
			arr = merge_sort(arr, q + 1, r);
			arr = merge(arr, p, q, r);
		}

		return (arr);
	}

	// Compute the average value in an array
	template <class T>
	double array_avg(T arr)
	{
		typename T::value_type avg = std::accumulate(arr.begin(), arr.end(), 0.0) / arr.size();

		return (avg);
	}

	// Compute the median value in an array
	template <class T>
	double array_median(T arr)
	{

		if (arr.empty())
		{
			std::cout << "The given array/vector is empty" << std::endl;
			return 0;
		}
		else
		{
			// First sort the vector
			std::sort(arr.begin(), arr.end());

			if (arr.size() % 2 == 0) // If the array has an ODD number of elements
			{
				return (arr[arr.size() / 2 - 1] + arr[arr.size() / 2]) / 2;
			}
			else // If the array has an EVEN number of elements
			{
				return arr[arr.size() / 2];
			}
		}

	}

	// Compute the standard deviation an array
	template <class T>
	double array_standard_deviation(T arr)
	{
		double mean = array_avg(arr);

		double sum = 0;
		for (int i = 0; i < arr.size(); i++)
		{
			sum += pow((arr[i] - mean), 2);
		}

		double standard_deviation = sqrt(sum / arr.size());

		return standard_deviation;
	}

	// Find the intersection between two arrays
	template <class T>
	T array_intersection(T arr_1, T arr_2)
	{
		std::sort(arr_1.begin(), arr_1.end());
		std::sort(arr_2.begin(), arr_2.end());

		T arr_intersection;

		std::set_intersection(arr_1.begin(), arr_1.end(),
			arr_2.begin(), arr_2.end(),
			std::back_inserter(arr_intersection));

		return arr_intersection;
	}
}

namespace math_utils
{
	// Find the roots of a quadratic equation
	inline std::pair<std::complex<float>, std::complex<float>> solve_quadratic(float a, float b, float c)
	{
		// The roots
		std::complex<float> x1 = 0;
		std::complex<float> x2 = 0;

		// First compute the determinent
		float determinant = b*b - 4 * a*c;

		if (determinant > 0)
		{
			x1 = (-b + sqrt(determinant)) / (a * 2);
			x2 = (-b - sqrt(determinant)) / (a * 2);
		}

		else if (determinant == 0)
		{
			x1 = (-b + sqrt(determinant)) / (a * 2);
			x2 = (-b + sqrt(determinant)) / (a * 2);
		}
		else if (determinant < 0)
		{
			x1.real(-b / (a * 2));
			x1.imag(sqrt(-determinant) / (a * 2));
			x2.real(-b / (a * 2));
			x2.imag(-sqrt(-determinant) / (a * 2));
		}

		return std::make_pair(x1, x2);
	}

	// Find the magnitude of a complex number
	inline float complex_mag(std::complex<float> complex_num)
	{
		return sqrt(pow(complex_num.real(), 2) + pow(complex_num.imag(), 2));
	}

	// Find the phase in degrees of a complex number 
	inline float complex_phase(std::complex<float> complex_num)
	{
		return (atan(complex_num.imag() / (complex_num.real())) * 180 / constants::PI);
	}

	// Compute the factorial
	inline int factorial(int num)
	{
		int fact = std::tgamma(num + 1);
		return fact;
	}

	// Compute number of possible permutations
	inline int permutations(int set_size, int subset_size)
	{
		int perm = factorial(set_size) / factorial(set_size - subset_size);
		return perm;
	}

	// Compute number of possible combinations
	inline int combinations(int set_size, int subset_size)
	{
		int comb = factorial(set_size) / (factorial(subset_size) * factorial(set_size - subset_size));
		return comb;
	}

	// Euclidean distance between 2D points
	inline float euclidean_dist_2d(std::pair<float, float> point_1, std::pair<float, float> point_2)
	{
		float x_dist = pow(point_1.first - point_2.first, 2);
		float y_dist = pow(point_1.second - point_2.second, 2);
		float dist = sqrt(x_dist + y_dist);
		return dist;
	}

	// RANSAC line fitting algorithm for a vector of 2D points
	inline std::pair<float, float> ransac_2d(std::vector<std::pair<float, float>> data, int num_iters=10, float thresh_dist=5, float inlier_ratio=0.5)
	{
		// data: a 2xn dataset with #n data points
		// num: the minimum number of points. 
		// num_iters : the number of iterations
		// thresh_dist : the threshold of the distances between points and the fitting line
		// inlier_ratio : the threshold of the number of inliers

		int best_inlier_num = 0;
		float best_slope = 0;
		float best_y_intercept = 0;

		for (int i = 0; i < num_iters; i++)
		{
			// Randomly select two points
			std::random_shuffle(std::begin(data), std::end(data));
			std::pair<float, float> point_1 = data.at(0);
			std::pair<float, float> point_2 = data.at(1);

			// Based on the randomly selected points, compute a temporary line model
			float temp_slope = (point_2.second - point_1.second) / (point_2.first - point_1.first);
			float temp_y_intercept = point_2.second - (temp_slope * point_2.first);

			// Compute the distance between all points and the temporary fitting line
			std::vector<float> distances;
			for (int j = 0; j < data.size(); j++)
			{
				std::pair<float, float> actual_pos = data.at(j);
				std::pair<float, float> fitted_pos = std::make_pair(data.at(0).first, data.at(0).first * temp_slope + temp_y_intercept);
				float curr_dist = euclidean_dist_2d(actual_pos, fitted_pos);
				distances.push_back(curr_dist);
			}

			// Compute the inliers with distances smaller than the threshold
			int inliers_count = 0;
			for (int j = 0; j < distances.size(); j++)
			{
				if (distances.at(j) <= thresh_dist)
				{
					inliers_count += 1;
				}
			}

			// Update the number of inliers and fitting model if better model is found 
			if (inliers_count >= round(inlier_ratio * data.size()) && inliers_count > best_inlier_num)
			{
				best_inlier_num = inliers_count;
				best_slope = temp_slope;
				best_y_intercept = temp_y_intercept;
			}

		}

		std::pair<float, float> best_line_params = std::make_pair(best_slope, best_y_intercept);

		return best_line_params;
	}

}

namespace str_utils
{
	// Converts all string characters to lower case
	inline void to_lower_case(std::string &given_string)
	{
		for (auto it = std::begin(given_string); it != std::end(given_string); it++)
		{
			if (*it >= 'A' && *it <= 'Z')
				*it += ('a' - 'A');
		}
	}
	
	// Converts all string characters to upper case
	inline void to_upper_case(std::string &given_string)
	{
		for (auto it = std::begin(given_string); it != std::end(given_string); it++)
		{
			if (*it >= 'a' && *it <= 'z')
				*it -= ('a' - 'A');
		}
	}

	// Outputs puts doubles/floats with "n" significant digits
	template <typename T>
	std::string to_string_with_precision(const T a_value, const int n = 6)
	{
		std::ostringstream out;
		out << std::setprecision(n) << a_value;
		return out.str();
	}

	// Delete all spaces within a string
	inline void delete_spaces(std::string &given_string)
	{
		given_string.erase(remove_if(given_string.begin(), given_string.end(), isspace), given_string.end());
	}

	// Trim all spaces on the left side of the string
	inline void trim_left(std::string &given_string)
	{
		// Count the number of spaces we need to erase
		int trim_count = 0;
		for (int i = 0; i < given_string.length(); i++)
		{
			if (isspace(given_string[i]))
			{
				trim_count++;
			}
			else if (!isspace(given_string[i]))
			{
				break;
			}
		}
		
		// Trim all of the spaces on the left
		if (trim_count > 0)
		{
			given_string.erase(0, trim_count);
		}
	}
	
	// Trim all spaces on the right side of the string
	inline void trim_right(std::string &given_string)
	{
		// Count the number of spaces we need to erase
		int trim_count = 0;
		for (int i = given_string.length() - 1; i >= 0; i--)
		{
			if (isspace(given_string[i]))
			{
				trim_count++;
			}
			else if (!isspace(given_string[i]))
			{
				break;
			}
		}

		// Trim all of the spaces on the end
		if (trim_count > 0)
		{
			int trim_start_index = given_string.length() - trim_count;
			given_string.erase(trim_start_index, given_string.length() - 1);
		}
	}

	// Trim all spaces on the two ends (left and right) of the string
	inline void trim_both_ends(std::string &given_string)
	{
		// Trim the left and right sides
		str_utils::trim_left(given_string);
		str_utils::trim_right(given_string);
	}

	// Extract all seperate word tokens in the string into a vector of strings
	inline std::vector<std::string> split_into_words(std::string given_string)
	{
		// Use a vector of strings to store the words
		// Use two iterators to track the start and end of each word
		std::vector<std::string> words;
		auto start_it = std::begin(given_string);
		auto end_it = std::begin(given_string);

		// Collect all of the words
		for (auto it = std::begin(given_string); it != std::end(given_string); it++)
		{
			end_it = it;
			if (isspace(*it))
			{
				// If we hit the next space, grab the word we just passed
				words.push_back(std::string(start_it, (end_it)));

				// Now that we have the word, jump over all of the spaces to the start of the next word
				for (auto it_2 = it; it_2 != std::end(given_string); it_2++)
				{
					if (isspace(*it_2))
					{
						start_it = it_2 + 1;
						end_it = it_2 + 1;
						break;
					}
				}
			}
			// If we've reached the end of the big string, extract the last word right away because we know there won't be anymore spaces
			else if(it == std::end(given_string) - 1) 
			{
				words.push_back(std::string(start_it, std::next(end_it)));
				break;
			}
		}

		return words;
	}

	// Count the number of occurences of a substring in a string
	inline int count_substrs(std::string big_string, std::string sub_string)
	{
		int count = 0;
		size_t nPos = big_string.find(sub_string, 0); // First occurrence

		// Count the rest of the occurences
		while (nPos != std::string::npos)
		{
			count++;
			nPos = big_string.find(sub_string, nPos + 1);
		}

		int count_int = static_cast<int>(count);

		return count_int;
	}

	// Get all of the substring indices from a string
	inline std::vector<int> get_substr_indices(std::string big_string, std::string sub_string)
	{
		std::vector<int> substr_indices;
		size_t nPos = big_string.find(sub_string, 0); // First occurrence
													  
		while (nPos != std::string::npos)
		{
			substr_indices.push_back(static_cast<int>(nPos));
			nPos = big_string.find(sub_string, nPos + 1);
		}

		return substr_indices;
	}
}