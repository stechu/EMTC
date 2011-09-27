 /*
 * heapbitset.h
 *
 *  Heap based bitset
 *  @author: Shumo Chu
 */

#ifndef HEAPBITSET_H
#define HEAPBITSET_H

#include <iostream>

class HeapBitset {
public:
	class BitProxy {
	public:
		BitProxy(char* byte, unsigned int index) :
			byte_(byte), index_(index) {

		}

		BitProxy& operator=(const BitProxy& rhs) {
			byte_ = rhs.byte_;
			index_ = rhs.index_;
			return *this;
		}

		BitProxy& operator=(bool bit) {
			if (bit && (*byte_ & (1 << index_)) == 0) {
				*byte_ ^= (1 << index_);
			} else if (!bit && *byte_ & (1 << index_)) {
				*byte_ ^= (1 << index_);
			}

			return *this;
		}

		operator bool() const {
			return *byte_ & (1 << index_);
		}

	private:
		char* byte_;
		unsigned int index_;
	};

	HeapBitset(unsigned int bits) :
		bits_(bits) {
		bitset_ = new char[(bits_ / 8)+1];
	}

	~HeapBitset() {
		delete[] bitset_;
	}

	const BitProxy operator[](unsigned int index) const {
		return BitProxy(&bitset_[index / 8], index % 8);
	}

	BitProxy operator[](unsigned int index) {
		return BitProxy(&bitset_[index / 8], index % 8);
	}

	void reset(){
		memset(bitset_,0,(bits_ / 8)+1);
	}

private:
	char* bitset_;
	unsigned int bits_;
};
#endif
