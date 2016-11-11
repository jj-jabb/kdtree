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

#ifndef INCLUDE_CSVREADER_H_
#define INCLUDE_CSVREADER_H_
#include <iterator>
#include <iostream>

#include <vector>
#include <string>
namespace uber {
//Borrowed from http://stackoverflow.com/questions/1120140/how-can-i-read-and-parse-csv-files-in-c
	class CSVRow {
	public:
		inline std::string const& operator[](std::size_t index) const {
			return data[index];
		}
		inline std::size_t size() const {
			return data.size();
		}
		void readNextRow(std::istream& str);
	private:
		std::vector<std::string> data;
	};
	class CSVIterator {
	public:
		typedef std::input_iterator_tag iterator_category;
		typedef CSVRow value_type;
		typedef std::size_t difference_type;
		typedef CSVRow* pointer;
		typedef CSVRow& reference;
		CSVIterator(std::istream& str) :
				str(str.good() ? &str : nullptr) {
			++(*this);
		}
		CSVIterator() :
				str(nullptr) {
		}
		CSVIterator& operator++();
		CSVIterator operator++(int);
		CSVRow const& operator*() const {
			return row;
		}
		CSVRow const* operator->() const {
			return &row;
		}
		inline bool operator==(CSVIterator const& rhs) {
			return ((this == &rhs) || ((this->str == nullptr) && (rhs.str == nullptr)));
		}
		inline bool operator!=(CSVIterator const& rhs) {
			return !((*this) == rhs);
		}
	private:
		std::istream* str;
		CSVRow row;
	};
	std::istream& operator>>(std::istream& str, CSVRow& data);
}
#endif
