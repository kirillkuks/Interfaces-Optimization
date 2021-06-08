#include "IMultiIndex.h"
#include <cstring>

namespace {
    class MultiIndexImpl : public IMultiIndex {
    public:
        MultiIndexImpl(size_t dim);

        IMultiIndex* clone() const override;

        size_t getDim() const override;
        size_t const* getData() const override;
        RC setData(size_t dim, size_t const* const& ptr_data) override;

        RC getAxisIndex(size_t axisIndex, size_t& val) const override;
        RC setAxisIndex(size_t axisIndex, size_t val) override;

        RC incAxisIndex(size_t axisIndex, size_t val) override;

        ~MultiIndexImpl();

    private:
        size_t n_dim;
    };
}

IMultiIndex::~IMultiIndex() {}

IMultiIndex* IMultiIndex::createMultiIndex(size_t dim, size_t* indecis) {
    if(!indecis) {
        return nullptr;
    }

    std::uint8_t* ptr = new std::uint8_t [sizeof(MultiIndexImpl) + sizeof(size_t) * dim];
    if(!ptr) {
        return nullptr;
    }

    IMultiIndex* multy = new(ptr) MultiIndexImpl(dim);
    std::memcpy(ptr + sizeof(MultiIndexImpl), indecis, dim * sizeof(size_t));
    return multy;
}

MultiIndexImpl::MultiIndexImpl(size_t dim) : n_dim{ dim } {}

MultiIndexImpl::~MultiIndexImpl() {}

IMultiIndex* MultiIndexImpl::clone() const {
    size_t size = sizeof(MultiIndexImpl) + n_dim * sizeof(size_t);
    std::uint8_t* data = new std::uint8_t[size];
    if(!data) {
        return nullptr;
    }
    std::memcpy(data, this, size);
    return  (IMultiIndex*)data;
}

size_t MultiIndexImpl::getDim() const {
    return n_dim;
}

size_t const* MultiIndexImpl::getData() const {
    return (size_t*)((std::uint8_t*)(this) + sizeof(MultiIndexImpl));
}

RC MultiIndexImpl::setData(size_t dim, size_t const* const& ptr_data) {
    if(!ptr_data) {
        return RC::NULLPTR_ERROR;
    }
    if(n_dim != dim) {
        return RC::MISMATCHING_DIMENSIONS;
    }
    std::memcpy((std::uint8_t*)(this) + sizeof(MultiIndexImpl), ptr_data, dim * sizeof(size_t));
    return RC::SUCCESS;
}

RC MultiIndexImpl::getAxisIndex(size_t axisIndex, size_t& val) const {
    if(axisIndex >= n_dim) {
        return RC::INDEX_OUT_OF_BOUND;
    }
    val = *(size_t*)((std::uint8_t*)(this) + sizeof(MultiIndexImpl) + axisIndex * sizeof(size_t));
    return RC::SUCCESS;
}

RC MultiIndexImpl::setAxisIndex(size_t axisIndex, size_t val) {
    if(axisIndex >= n_dim) {
        return RC::INDEX_OUT_OF_BOUND;
    }
    *(size_t*)((std::uint8_t*)(this) + sizeof(MultiIndexImpl) + axisIndex * sizeof(size_t)) = val;
    return RC::SUCCESS;
}

RC MultiIndexImpl::incAxisIndex(size_t axisIndex, size_t val) {
    ++*(size_t*)((std::uint8_t*)(this) + sizeof(MultiIndexImpl) + axisIndex * sizeof(size_t));
    return RC::SUCCESS;
}
