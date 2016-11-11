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
#include <fstream>
#include <sstream>
#include "CSVReader.h"

//Borrowed from http://stackoverflow.com/questions/1120140/how-can-i-read-and-parse-csv-files-in-c
namespace uber {
	CSVIterator& CSVIterator::operator++() {
		if (str) {
			if (!((*str) >> row)) {
				str = NULL;
			}
		}
		return *this;
	}
	CSVIterator CSVIterator::operator++(int) {
		CSVIterator tmp(*this);
		++(*this);
		return tmp;
	}

	void CSVRow::readNextRow(std::istream& str) {
		std::string line;
		std::getline(str, line);

		std::stringstream lineStream(line);
		std::string cell;

		data.clear();
		while (std::getline(lineStream, cell, ',')) {
			data.push_back(cell);
		}
	}
	std::istream& operator>>(std::istream& str, CSVRow& data) {
		data.readNextRow(str);
		return str;
	}
}
