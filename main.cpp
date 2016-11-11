/*
 * Copyright(C) 2016, Blake C. Lucas, Ph.D. (img.science@gmail.com)
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
#include "UberKDTree.h"
#include <iomanip>
using namespace uber;
template<class T, int N> void ClosestPointBruteForce(const std::vector<vec<T, N>>& data, const vec<T, N>& query, vec<T, N>& result) {
	double minDistance = std::numeric_limits<double>::max();
	for (vec<T, N> pt : data) {
		double d = uber::distance(pt, query);
		if (d < minDistance) {
			result = pt;
			minDistance = d;
		}
	}
}
int main(int argc, char *argv[]) {
	using namespace uber;
	const int SPLIT_POLICIES = 3;
	const int DIMENSION_POLICIES = 3;
	int counter = 0;
	size_t index;
	Vector3f trainData;
	Vector3f testData;
	Vector2f resultData;
	std::cout << "--> Reading data ..." << std::endl;
	ReadPointData("sample_data.csv", trainData);
	ReadPointData("query_data.csv", testData);
	std::cout<<"--> Testing serialization ..."<<std::endl;
	{
		KDTree tree;
		tree.build(trainData);
		WriteKDTreeToFile("tree.json", tree);
	}
	std::cout<<"--> Testing deserialization ..."<<std::endl;
	{
		KDTree tree;
		ReadKDTreeFromFile("tree.json", tree);
		counter = 0;
		for (float3 query : testData) {
			float3 result1, result2;
			ClosestPointBruteForce(trainData, query, result1);
			tree.closestPoint(query, result2, index);
			resultData.push_back(float2(uber::distance(query, result2), static_cast<float>(index)));
			double delta = uber::distance(result1, result2);
			if (delta > 1E-10f) {
				std::cout << "Query: " << query << " Kd-Tree Result: " << result2 << " Brute Force Result: " << result1 << " :: " << delta << std::endl;
			}
			if (counter % 10 == 0)
				std::cout << "Testing " << std::setw(4) << 100 * counter / (float) testData.size() << "% ..." << std::endl;
			counter++;
		}
		WritePointData("result_data_serialization.csv", resultData);
		resultData.clear();
	}
	std::vector<std::shared_ptr<KDTree>> kdTrees;
		for (int i = 0; i < SPLIT_POLICIES; i++) {
			for (int j = 0; j < DIMENSION_POLICIES; j++) {
				auto sPolicy = static_cast<KDTree::SplitPolicy>(i);
				auto dPolicy = static_cast<KDTree::DimensionPolicy>(j);
				std::shared_ptr<KDTree> tree(new KDTree(sPolicy, dPolicy));
				std::cout << "--> Building tree for " << sPolicy << "/" << dPolicy << " ..." << std::endl;
				tree->build(trainData);
				kdTrees.push_back(tree);
			}
		}
	std::cout << "--> Done Building Trees" << std::endl;
	counter = 0;
	for (float3 query : testData) {
		float3 result1, result2;
		ClosestPointBruteForce(trainData, query, result1);
		for (auto tree : kdTrees) {
			tree->closestPoint(query, result2, index);
			//Result set includes results for all variants of kd-tree
			resultData.push_back(float2(uber::distance(query, result2), static_cast<float>(index)));
			double delta = uber::distance(result1, result2);
			if (delta > 1E-10f) {
				std::cout << "Query: " << query << " Kd-Tree Result: " << result2 << " Brute Force Result: " << result1 << " :: " << delta << std::endl;
			}
		}
		if (counter % 10 == 0)
			std::cout << "--> Testing " << std::setw(4) << 100 * counter / (float) testData.size() << "% ..." << std::endl;
		counter++;
	}
	WritePointData("result_data_all.csv", resultData);
	std::cout << "--> Done!" << std::endl;
}

