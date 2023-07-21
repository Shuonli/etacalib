#ifndef ROLLING_BUFFER_H
#define ROLLING_BUFFER_H

#include <vector>
#include <TLorentzVector.h>

class RollingBuffer {
public:
    RollingBuffer(int maxLength = 1000);

    void AddFront(const std::vector<TLorentzVector>& vec);
    int GetSize() const;
    void Flush();
    std::vector<TLorentzVector> GetElement(int index) const;

private:
    void PopBack();

    int maxLength;
    std::vector<std::vector<TLorentzVector>> buffer;
};

#endif // ROLLING_BUFFER_H