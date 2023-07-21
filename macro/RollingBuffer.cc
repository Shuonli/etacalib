#include "RollingBuffer.h"
#include <vector>

RollingBuffer::RollingBuffer(int maxLength) : maxLength(maxLength) {}

void RollingBuffer::AddFront(const std::vector<TLorentzVector>& vec) {
    buffer.insert(buffer.begin(), vec);
    if (buffer.size() > maxLength) {
        PopBack();
    }
}

int RollingBuffer::GetSize() const {
    return buffer.size();
}

void RollingBuffer::Flush() {
    buffer.clear();
}

std::vector<TLorentzVector> RollingBuffer::GetElement(int index) const {
    if (index >= 0 && index < buffer.size()) {
        return buffer[index];
    }
    return std::vector<TLorentzVector>(); // Return an empty vector if index is out of bounds
}

void RollingBuffer::PopBack() {
    buffer.pop_back();
}
