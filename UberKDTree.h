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
#ifndef INCLUDE_UBERKDTREE_H_
#define INCLUDE_UBERKDTREE_H_
#include <cereal/types/memory.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/cereal.hpp>
#include <vector>
#include <memory>
#include "UberMath.h"
namespace uber {
	struct KDTreeNode {
		std::shared_ptr<KDTreeNode> leftChild;
		std::shared_ptr<KDTreeNode> rightChild;
		int depth;
		virtual bool isLeaf() const {
			return false;
		}
		virtual size_t size() const {
			return 0;
		}
		virtual void clear() {
		}
		KDTreeNode(const std::shared_ptr<KDTreeNode>& leftChild = nullptr, const std::shared_ptr<KDTreeNode>& rightChild = nullptr, int splitDimension = -1,
				double splitPosition = 0) :
				leftChild(leftChild), rightChild(rightChild), depth(0){
		}

		template<class Archive> void serialize(Archive & archive) {
			archive(CEREAL_NVP(depth), CEREAL_NVP(leftChild), CEREAL_NVP(rightChild));
		}

		virtual inline ~KDTreeNode() {
		}
	};
	template<class T, int N> struct KDTreeNodeLeaf: KDTreeNode {
		box<T, N> bounds;
		std::vector<vec<T, N>> data;
		std::vector<size_t> indexes;
		KDTreeNodeLeaf(const std::vector<vec<T, N>>& data = std::vector<vec<T, N>>()) :
				data(data) {
			indexes.resize(data.size());
			size_t counter = 0;
			vec<T, N> minPt(T(1E30));
			vec<T, N> maxPt(T(-1E30));
			for (const vec<T, N>& pt : data) {
				minPt = uber::min(pt, minPt);
				maxPt = uber::max(pt, maxPt);
			}
			bounds.position = minPt;
			bounds.dimensions = maxPt - minPt;
			for (size_t& idx : indexes) {
				idx = counter++;
			}
		}
		virtual void clear() override {
			data.clear();
			indexes.clear();
		}
		virtual bool isLeaf() const override {
			return (data.size() > 0);
		}
		virtual size_t size() const override {
			return data.size();
		}
		template<class Archive> void serialize(Archive & archive) {
			archive(cereal::base_class<KDTreeNode>(this), CEREAL_NVP(bounds), CEREAL_NVP(data), CEREAL_NVP(indexes));
		}
	};
	struct KDTreeDistance {
		double value;
		KDTreeNode* node;
		KDTreeDistance() :
				value(0), node(nullptr) {

		}
		template<class T, int N> KDTreeDistance(const vec<T, N>& pt, KDTreeNodeLeaf<T, N>* node) :
				node(node) {
			if (node->bounds.contains(pt)) {
				value=-1;
			} else {
				vec<T, N> minPoint = node->bounds.min();
				vec<T, N> maxPoint = node->bounds.max();
				vec<T, N> ref;
				for (int n = 0; n < N; n++) {
					ref[n] = uber::clamp(pt[n], minPoint[n], maxPoint[n]);
				}
				value=distance(pt, ref);
			}

		}
		bool operator()(const KDTreeDistance& a, const KDTreeDistance& b) {
			return (a.value > b.value);
		}
	};

	class KDTree {
	public:
		enum class SplitPolicy {
			Mean = 0, Median = 1, Midpoint = 2
		};
		enum class DimensionPolicy {
			NumericalOrder = 0, LargestVariance = 1, LargestRange = 2
		};
	private:
		static const int MIN_NODE_SIZE = 16;
		static const int HISTOGRAM_BINS = 256;
		SplitPolicy splitPolicy;
		DimensionPolicy dimPolicy;
		std::shared_ptr<KDTreeNode> root;

		template<class T, int N> void meanAndVariance(const std::vector<vec<T, N>>& data, vec<T, N>& mean, vec<T, N>& var) {
			vec<T, N> sum(T(0));
			vec<T, N> sumSqr(T(0));
			if (data.size() == 0) {
				mean = vec<T, N>(T(0));
				var = vec<T, N>(T(0));
				return;
			}
			for (const vec<T, N>& val : data) {
				sum += val;
				sumSqr += val;
			}
			mean = sum / T(data.size());
			var = sumSqr / T(data.size()) - mean * mean;
		}
		template<class T, int N> void range(const std::vector<vec<T, N>>& data, vec<T, N>& minValue, vec<T, N>& maxValue) {
			if (data.size() == 0) {
				minValue = vec<T, N>(T(0));
				maxValue = vec<T, N>(T(0));
				return;
			}
			minValue = data.front();
			maxValue = data.front();
			for (const vec<T, N>& val : data) {
				minValue = uber::min(minValue, val);
				maxValue = uber::max(maxValue, val);
			}
		}
		template<class T, int N> void approximateMedian(const std::vector<vec<T, N>>& data, T& median, int dim) {
			static const float ZERO_TOLERANCE = 1E-6f;
			std::vector<int> bins(HISTOGRAM_BINS + 1, 0);
			if (data.size() == 0) {
				median = T(0);
				return;
			}
			vec<T, N> minVal, maxVal;
			range(data, minVal, maxVal);
			T binSize = (maxVal[dim] - minVal[dim]) / T(HISTOGRAM_BINS);
			binSize = std::max(binSize, T(ZERO_TOLERANCE));
			int samples = (int) data.size();
			for (const vec<T, N>& val : data) {
				T b = std::floor((val[dim] - minVal[dim]) / binSize + T(0.5));//Round to nearest bin
				bins[std::max(std::min(int(b), HISTOGRAM_BINS), 0)]++;
			}
			int sum = 0;
			int pos = 0;
			for (pos = 0; pos <= HISTOGRAM_BINS; pos++) {
				if (2 * sum < samples) {
					sum += bins[pos];
				} else {
					break;
				}
			}
			median = pos * binSize + minVal[dim];//Good guess for median without the overhead of sorting.

		}
		template<class T, int N> void buildInternal(const std::vector<vec<T, N>>& data);
		template<class T, int N> bool closestPointInternal(const vec<T, N>& query, vec<T, N>& result, size_t& index) const;
	public:
		KDTree(const SplitPolicy& splitPolicy = SplitPolicy::Median, const DimensionPolicy& dimPolicy = DimensionPolicy::LargestVariance) :
				splitPolicy(splitPolicy), dimPolicy(dimPolicy) {
		}
		template<class Archive> void serialize(Archive & archive) {
			archive(CEREAL_NVP(root));
		}
		inline ~KDTree() {
		}

		void build(const std::vector<float2>& data);
		void build(const std::vector<float3>& data);
		void build(const std::vector<double2>& data);
		void build(const std::vector<double3>& data);
		bool closestPoint(const float2& query, float2& result, size_t& index) const;
		bool closestPoint(const float3& query, float3& result, size_t& index) const;
		bool closestPoint(const double2& query, double2& result, size_t& index) const;
		bool closestPoint(const double3& query, double3& result, size_t& index) const;
	};
	template<class C, class R> std::basic_ostream<C, R> & operator <<(std::basic_ostream<C, R> & ss, const KDTree::SplitPolicy& type) {
		switch (type) {
		case KDTree::SplitPolicy::Mean:
			return ss << "Mean";
		case KDTree::SplitPolicy::Median:
			return ss << "Median";
		case KDTree::SplitPolicy::Midpoint:
			return ss << "Midpoint";
		}
		return ss;
	}
	template<class C, class R> std::basic_ostream<C, R> & operator <<(std::basic_ostream<C, R> & ss, const KDTree::DimensionPolicy& type) {
		switch (type) {
		case KDTree::DimensionPolicy::NumericalOrder:
			return ss << "Mean";
		case KDTree::DimensionPolicy::LargestVariance:
			return ss << "Largest Variance";
		case KDTree::DimensionPolicy::LargestRange:
			return ss << "Largest Range";
		}
		return ss;
	}
	typedef std::vector<float2> Vector2f;
	typedef std::vector<float3> Vector3f;
	typedef std::vector<double2> Vector2d;
	typedef std::vector<double3> Vector3d;
	bool ReadPointData(const std::string& file, Vector2f& data);
	bool ReadPointData(const std::string& file, Vector3f& data);
	bool ReadPointData(const std::string& file, Vector2d& data);
	bool ReadPointData(const std::string& file, Vector3d& data);

	bool WritePointData(const std::string& file, const Vector2f& data);
	bool WritePointData(const std::string& file, const Vector3f& data);
	bool WritePointData(const std::string& file, const Vector2d& data);
	bool WritePointData(const std::string& file, const Vector3d& data);

	void WriteKDTreeToFile(const std::string& file, const KDTree& tree);
	void ReadKDTreeFromFile(const std::string& file, KDTree& tree);
	typedef KDTreeNodeLeaf<float, 2> KDTreeNodeLeaf2f;
	typedef KDTreeNodeLeaf<float, 3> KDTreeNodeLeaf3f;
	typedef KDTreeNodeLeaf<double, 2> KDTreeNodeLeaf2d;
	typedef KDTreeNodeLeaf<double, 3> KDTreeNodeLeaf3d;
}
CEREAL_REGISTER_TYPE(uber::KDTreeNodeLeaf2f);
CEREAL_REGISTER_TYPE(uber::KDTreeNodeLeaf3f);
CEREAL_REGISTER_TYPE(uber::KDTreeNodeLeaf2d);
CEREAL_REGISTER_TYPE(uber::KDTreeNodeLeaf3d);

#endif
