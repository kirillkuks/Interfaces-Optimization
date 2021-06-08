#pragma once
#include "IVector.h"
#include "IMultiIndex.h"
#include <cstddef>

class ICompactControlBlock {
public:
    virtual RC get(IMultiIndex* const& currentIndex, IMultiIndex const* const& bypassOrder) const = 0;
    virtual RC get(IMultiIndex const* const& currentIndex, IVector* const& val) const = 0;

    virtual ~ICompactControlBlock() = 0;

private:
    ICompactControlBlock(ICompactControlBlock const &);
    ICompactControlBlock& operator=(ICompactControlBlock const &);

protected:
    ICompactControlBlock() = default;
};
