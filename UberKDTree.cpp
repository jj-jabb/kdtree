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
#include "CSVReader.h"
#include <fstream>
#include <queue>
#include <cereal/archives/json.hpp>
namespace uber {
	void WriteKDTreeToFile(const std::string& file, const KDTree& tree) {
		std::ofstream os(file);
		cereal::JSONOutputArchive archive(os);
		archive(cereal::make_nvp("kdtree", tree));
	}
	void ReadKDTreeFromFile(const std::string& file, KDTree& tree) {
		std::ifstream os(file);
		cereal::JSONInputArchive archive(os);
		archive(cereal::make_nvp("kdtree", tree));

	}
	template<class T, int N> void KDTree::buildInternal(const std::vector<vec<T, N>>& data) {
		vec<T, N> mean(T(0));
		vec<T, N> var(T(0));
		vec<T, N> minVal(T(0));
		vec<T, N> maxVal(T(0));
		vec<T, N> minPt, maxPt;
		int splitDimension = 0;
		size_t index = 0;
		T splitPosition = T(0);
		std::queue<std::shared_ptr<KDTreeNodeLeaf<T, N>>>nodes;
		nodes.push(std::shared_ptr<KDTreeNodeLeaf<T, N>>(new KDTreeNodeLeaf<T, N>(data)));
		root = std::dynamic_pointer_cast<KDTreeNode>(nodes.front());
		while (!nodes.empty()) {
			std::shared_ptr<KDTreeNodeLeaf<T, N>> node = nodes.front();
			nodes.pop();
			if (node->size() > MIN_NODE_SIZE) {
				switch (dimPolicy) {
				case DimensionPolicy::NumericalOrder:
					splitDimension = node->depth % N;
					break;
				case DimensionPolicy::LargestVariance: {
					meanAndVariance(node->data, mean, var);
					T maxVar(T(0));
					for (int n = 0; n < N; n++) {
						if (var[n] > maxVar) {
							maxVar = var[n];
							splitDimension = n;
						}
					}
				}
					break;
				case DimensionPolicy::LargestRange: {
					range(node->data, minVal, maxVal);
					vec<T, N> diff = maxVal - minVal;
					T maxDiff(T(0));
					for (int n = 0; n < N; n++) {
						if (diff[n] > maxDiff) {
							maxDiff = diff[n];
							splitDimension = n;
						}
					}
				}
					break;
				}
				switch (splitPolicy) {
				case SplitPolicy::Mean:
					if (dimPolicy != DimensionPolicy::LargestVariance) {
						meanAndVariance(node->data, mean, var);
					}
					splitPosition = mean[splitDimension];
					break;
				case SplitPolicy::Median:
					approximateMedian(node->data, splitPosition, splitDimension);
					break;
				case SplitPolicy::Midpoint:
					if (dimPolicy != DimensionPolicy::LargestRange) {
						range(node->data, minVal, maxVal);
					}
					splitPosition = T(0.5) * (maxVal[splitDimension] + minVal[splitDimension]);
					break;
				}
				std::shared_ptr<KDTreeNodeLeaf<T, N>> leftNode = std::shared_ptr<KDTreeNodeLeaf<T, N>>(new KDTreeNodeLeaf<T, N>());
				std::shared_ptr<KDTreeNodeLeaf<T, N>> rightNode = std::shared_ptr<KDTreeNodeLeaf<T, N>>(new KDTreeNodeLeaf<T, N>());
				index = 0;
				for (vec<T, N> val : node->data) {
					if (val[splitDimension] < splitPosition) {
						leftNode->data.push_back(val);
						leftNode->indexes.push_back(node->indexes[index]);
					} else {
						rightNode->data.push_back(val);
						rightNode->indexes.push_back(node->indexes[index]);

					}
					index++;
				}
				leftNode->depth = node->depth + 1;
				rightNode->depth = node->depth + 1;
				//Construct bounding box for each KD Node
				minPt = node->bounds.min();
				maxPt = node->bounds.max();
				maxPt[splitDimension] = splitPosition;
				leftNode->bounds = box<T, N>(minPt, maxPt - minPt);
				minPt = node->bounds.min();
				maxPt = node->bounds.max();
				minPt[splitDimension] = splitPosition;
				rightNode->bounds = box<T, N>(minPt, maxPt - minPt);
				node->leftChild = leftNode;
				node->rightChild = rightNode;
				node->clear();
				nodes.push(leftNode);
				nodes.push(rightNode);

			}
		}
	}
	template<class T, int N> bool KDTree::closestPointInternal(const vec<T, N>& query, vec<T, N>& result, size_t& index) const {
		if (root.get() == nullptr)
			throw new std::runtime_error("KD-Tree has not been initialized.");

		double minDistance = std::numeric_limits<T>::max();
		std::priority_queue<KDTreeDistance, std::vector<KDTreeDistance>, KDTreeDistance> queue;
		queue.push(KDTreeDistance(query, dynamic_cast<KDTreeNodeLeaf<T, N>*>(root.get())));
		size_t counter = 0;
		bool found = false;
		while (queue.size() > 0) {
			KDTreeDistance kd = queue.top();
			queue.pop();
			if (kd.value <= minDistance || kd.node->isLeaf()) {
				if (kd.node->isLeaf()) {
					auto leaf = dynamic_cast<KDTreeNodeLeaf<T, N>*>(kd.node);
					counter = 0;
					for (vec<T, N> pt : leaf->data) {
						double d = distance(pt, query);
						if (d <= minDistance) {
							result = pt;
							index = leaf->indexes[counter];
							minDistance = d;
							found = true;
						}
						counter++;
					}
				} else {
					queue.push(KDTreeDistance(query, dynamic_cast<KDTreeNodeLeaf<T, N>*>(kd.node->leftChild.get())));
					queue.push(KDTreeDistance(query, dynamic_cast<KDTreeNodeLeaf<T, N>*>(kd.node->rightChild.get())));
				}
			}
		}
		return found;
	}
	void KDTree::build(const std::vector<float2>& data) {
		buildInternal(data);
	}
	void KDTree::build(const std::vector<float3>& data) {
		buildInternal(data);
	}
	void KDTree::build(const std::vector<double2>& data) {
		buildInternal(data);
	}
	void KDTree::build(const std::vector<double3>& data) {
		buildInternal(data);
	}

	bool KDTree::closestPoint(const float2& query, float2& result, size_t& index) const {
		return closestPointInternal(query, result, index);
	}
	bool KDTree::closestPoint(const float3& query, float3& result, size_t& index) const {
		return closestPointInternal(query, result, index);
	}
	bool KDTree::closestPoint(const double2& query, double2& result, size_t& index) const {
		return closestPointInternal(query, result, index);
	}
	bool KDTree::closestPoint(const double3& query, double3& result, size_t& index) const {
		return closestPointInternal(query, result, index);
	}
	template<class T, int N> bool ReadPointDataInternal(const std::string& file, std::vector<vec<T, N>>& data) {
		std::ifstream str(file);
		if (str.is_open()) {
			for (CSVIterator iter(str); iter != CSVIterator(); iter++) {
				CSVRow row = *iter;
				if (row.size() == N) {
					vec<T, N> val;
					for (int n = 0; n < N; n++) {
						val[n] = T(std::atof(row[n].c_str()));
					}
					data.push_back(val);
				}
			}
			return true;
		}
		return false;
	}

	bool ReadPointData(const std::string& file, Vector2f& data) {
		return ReadPointDataInternal(file, data);
	}
	bool ReadPointData(const std::string& file, Vector3f& data) {
		return ReadPointDataInternal(file, data);
	}
	bool ReadPointData(const std::string& file, Vector2d& data) {
		return ReadPointDataInternal(file, data);
	}
	bool ReadPointData(const std::string& file, Vector3d& data) {
		return ReadPointDataInternal(file, data);
	}

	template<class T, int N> bool WritePointDataInternal(const std::string& file, const std::vector<vec<T, N>>& data) {
		std::ofstream out(file.c_str());
		if (out.is_open()) {
			for (vec<T, N> val : data) {
				for (int n = 0; n < N; n++) {
					if (n < N - 1) {
						out << val[n] << ",";
					} else {
						out << val[n] << "\n";
					}
				}
			}
			return true;
		}
		return false;
	}
	bool WritePointData(const std::string& file, const Vector2f& data) {
		return WritePointDataInternal(file, data);
	}
	bool WritePointData(const std::string& file, const Vector3f& data) {
		return WritePointDataInternal(file, data);
	}
	bool WritePointData(const std::string& file, const Vector2d& data) {
		return WritePointDataInternal(file, data);
	}
	bool WritePointData(const std::string& file, const Vector3d& data) {
		return WritePointDataInternal(file, data);
	}
}
