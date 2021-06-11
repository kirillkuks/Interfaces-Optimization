#include <memory>
#include <cmath>
#include <cstring>
#include "ICompact.h"
#include "ICompactControlBlock.h"

namespace {
    class CompactImpl : public ICompact {
    public:
        CompactImpl(size_t dim);

        ICompact* clone() const override;

        bool isInside(IVector const* const& vec) const override;
        RC getVectorCopy(IMultiIndex const* index, IVector*& val) const override;
        RC getVectorCoords(IMultiIndex const* index, IVector* const& val) const override;

        RC getLeftBoundary(IVector*& vec) const override;
        RC getRightBoundary(IVector*& vec) const override;

        size_t getDim() const override;
        IMultiIndex* getGrid() const override;

        class IteratorImpl : public IIterator {
        public:
            static IIterator* createIterator(IVector const* vector, IMultiIndex const* const& index, IMultiIndex const* const& bypassOrder, std::shared_ptr<ICompactControlBlock> const& controlBlock);

            IIterator* getNext() override;
            IIterator* clone() const override;

            RC next() override;

            bool isValid() const override;

            RC getVectorCopy(IVector*& val) const override;
            RC getVectorCoords(IVector* const& val) const override;

            ~IteratorImpl();

            static ILogger* pLogger;

        private:
            IteratorImpl(IVector* vector, IMultiIndex* index, IMultiIndex* bypassOrder, std::shared_ptr<ICompactControlBlock> const& controlBlock, size_t dim);

            IVector* vector;
            IMultiIndex* index;
            IMultiIndex* bypassOrder;

            std::shared_ptr<ICompactControlBlock> controlBlock;

            size_t dim;
        };

        class CompactControlBlockImpl : public ICompactControlBlock {
        public:
            CompactControlBlockImpl(ICompact const* const& compact);

            RC get(IMultiIndex* const& currentIndex, IMultiIndex const* const& bypassOrder) const override;
            RC get(IMultiIndex const* const& currentIndex, IVector* const& val) const override;

            static size_t axisToIncreament(IMultiIndex const* const& currentIndex, IMultiIndex const* const& bypassOrder, IMultiIndex const* grid);

            ~CompactControlBlockImpl();

        private:
            ICompact const* compact;

        };

        IIterator* getIterator(IMultiIndex const* const& index, IMultiIndex const* const& bypassOrder) const override;
        IIterator* getBegin(IMultiIndex const* const& bypassOrder) const override;
        IIterator* getEnd(IMultiIndex const* const& bypassOrder) const override;

        static double* axisIncrement(IVector const* op1, IVector const* op2, IMultiIndex const* grid);

        size_t sizeAllocated() const;
        double const* getAxisIncrement() const;
        RC ptrVectorInCompact(IMultiIndex const* index, double*& data) const;

        static RC leastVectorComponents(IVector const* op1, IVector const* op2, IVector*& res);
        static RC biggestVectorComponents(IVector const* op1, IVector const* op2, IVector*& res);

        static bool isValid(IVector const*& op1, IVector const*& op2, IMultiIndex const* const& grid);

        ~CompactImpl();

        static ILogger* pLogger;

    private:
        size_t dim;

        std::shared_ptr<CompactControlBlockImpl> controlBlock;
    };
}

ILogger* CompactImpl::pLogger = nullptr;
ILogger* CompactImpl::IteratorImpl::pLogger = nullptr;

inline ICompact::IIterator::~IIterator() {}

RC ICompact::IIterator::setLogger(ILogger* const logger) {
    if(!logger) {
        return RC::NULLPTR_ERROR;
    }
    CompactImpl::IteratorImpl::pLogger = logger;
    return RC::SUCCESS;
}

ILogger* ICompact::IIterator::getLogger() {
    return CompactImpl::IteratorImpl::pLogger;
}

CompactImpl::IteratorImpl::IteratorImpl(IVector* vector, IMultiIndex* index, IMultiIndex* bypassOrder, std::shared_ptr<ICompactControlBlock> const& controlBlock, size_t dim)
    : vector{ vector }, index{ index }, bypassOrder{ bypassOrder }, controlBlock{ controlBlock }, dim{ dim } {}

CompactImpl::IteratorImpl::~IteratorImpl() {
    delete index; index = nullptr;
    delete bypassOrder; bypassOrder = nullptr;
    delete vector; vector = nullptr;
}

ICompact::IIterator* CompactImpl::IteratorImpl::createIterator(IVector const* vector, IMultiIndex const* const& index, IMultiIndex const* const& bypassOrder, std::shared_ptr<ICompactControlBlock> const& controlBlock) {
    if(!vector || !index || !bypassOrder) {
        if(CompactImpl::IteratorImpl::pLogger) {
            CompactImpl::IteratorImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    size_t dim = index->getDim();
    if(dim != bypassOrder->getDim()) {
        if(CompactImpl::IteratorImpl::pLogger) {
            CompactImpl::IteratorImpl::pLogger->warning(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    IMultiIndex* indexCopy = index->clone();
    if(!indexCopy) {
        if(CompactImpl::IteratorImpl::pLogger) {
            CompactImpl::IteratorImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    IMultiIndex* bypassOrderCopy = bypassOrder->clone();
    if(!bypassOrderCopy) {
        delete indexCopy; indexCopy = nullptr;
        if(CompactImpl::IteratorImpl::pLogger) {
            CompactImpl::IteratorImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    IVector* vec = vector->clone();
    if(!vec) {
        if(CompactImpl::IteratorImpl::pLogger) {
            CompactImpl::IteratorImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        delete indexCopy; indexCopy = nullptr;
        delete bypassOrderCopy; bypassOrderCopy = nullptr;
        return nullptr;
    }

    IIterator* it = new IteratorImpl(vec, indexCopy, bypassOrderCopy, controlBlock, dim);
    if(!it) {
        delete indexCopy; indexCopy = nullptr;
        delete bypassOrderCopy; bypassOrderCopy = nullptr;
        delete vec; vec = nullptr;
        if(CompactImpl::IteratorImpl::pLogger) {
            CompactImpl::IteratorImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    return it;

}

ICompact::IIterator* CompactImpl::IteratorImpl::getNext() {
    IIterator* it = clone();
    if(!it) {
        if(CompactImpl::IteratorImpl::pLogger) {
            CompactImpl::IteratorImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    RC rc = it->next();
    if(rc != RC::SUCCESS) {
        delete it; it = nullptr;
        if(CompactImpl::IteratorImpl::pLogger) {
            CompactImpl::IteratorImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    return it;
}

ICompact::IIterator* CompactImpl::IteratorImpl::clone() const {
    IIterator* it = createIterator(vector, index, bypassOrder, controlBlock);
    if(!it) {
        if(CompactImpl::IteratorImpl::pLogger) {
            CompactImpl::IteratorImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    return it;
}

RC CompactImpl::IteratorImpl::next() {
    RC rc = controlBlock->get(index, bypassOrder);
    if(rc != RC::SUCCESS) {
        delete index; index = nullptr;
        index = nullptr;
        if(CompactImpl::IteratorImpl::pLogger) {
            CompactImpl::IteratorImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return rc;
    }
    rc = controlBlock->get(index, vector);
    if(rc != RC::SUCCESS) {
        delete index; index = nullptr;
        index = nullptr;
        if(CompactImpl::IteratorImpl::pLogger) {
            CompactImpl::IteratorImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return rc;
    }
    return RC::SUCCESS;
}

bool CompactImpl::IteratorImpl::isValid() const {
    return index;
}

RC CompactImpl::IteratorImpl::getVectorCopy(IVector*& val) const {
    if(!vector) {
        if(CompactImpl::IteratorImpl::pLogger) {
            CompactImpl::IteratorImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::NULLPTR_ERROR;
    }
    val = vector->clone();
    if(!val) {
        if(CompactImpl::IteratorImpl::pLogger) {
            CompactImpl::IteratorImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::ALLOCATION_ERROR;
    }
    return RC::SUCCESS;
}

RC CompactImpl::IteratorImpl::getVectorCoords(IVector* const& val) const {
    if(!val || !vector) {
        if(CompactImpl::IteratorImpl::pLogger) {
            CompactImpl::IteratorImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::NULLPTR_ERROR;
    }
    if(val->getDim() != vector->getDim()) {
        if(CompactImpl::IteratorImpl::pLogger) {
            CompactImpl::IteratorImpl::pLogger->warning(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
        }
        return RC::MISMATCHING_DIMENSIONS;
    }

    RC rc = IVector::copyInstance(val, vector);
    if(rc != RC::SUCCESS) {
        if(CompactImpl::IteratorImpl::pLogger) {
            CompactImpl::IteratorImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return rc;
    }
    return RC::SUCCESS;
}

ICompactControlBlock::~ICompactControlBlock() {}

CompactImpl::CompactControlBlockImpl::CompactControlBlockImpl(ICompact const* const& compact) : compact{ compact } {}

CompactImpl::CompactControlBlockImpl::~CompactControlBlockImpl() {}

RC CompactImpl::CompactControlBlockImpl::get(IMultiIndex* const& currentIndex, IMultiIndex const* const& bypassOrder) const {
    if(!compact) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(RC::SOURCE_COMPACT_DESTROYED, __FILE__, __func__, __LINE__);
        }
        return RC::SOURCE_COMPACT_DESTROYED;
    }
    size_t dim = currentIndex->getDim();
    if(dim != bypassOrder->getDim()) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->warning(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
        }
        return RC::MISMATCHING_DIMENSIONS;
    }
    IMultiIndex* grid = compact->getGrid();
    if(!grid) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::ALLOCATION_ERROR;
    }

    size_t ind = axisToIncreament(currentIndex, bypassOrder, grid);
    if(ind == dim) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(RC::INDEX_OUT_OF_BOUND, __FILE__, __func__, __LINE__);
        }
        return RC::INDEX_OUT_OF_BOUND;
    }

    for(size_t i = 0; i < ind; ++i) {
        size_t id;
        bypassOrder->getAxisIndex(i, id);
        currentIndex->setAxisIndex(id, 0);
    }
    size_t id;
    bypassOrder->getAxisIndex(ind, id);
    currentIndex->incAxisIndex(id, 1);

    delete grid; grid = nullptr;
    return RC::SUCCESS;
}

RC CompactImpl::CompactControlBlockImpl::get(IMultiIndex const* const& currentIndex, IVector* const& val) const {
    RC rc = compact->getVectorCoords(currentIndex, val);

    if(rc != RC::SUCCESS) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return rc;
    }
    return RC::SUCCESS;
}

size_t CompactImpl::CompactControlBlockImpl::axisToIncreament(IMultiIndex const* const& currentIndex, IMultiIndex const* const& bypassOrder, IMultiIndex const* grid) {
    size_t dim = grid->getDim();
    size_t const* ptrGrid = grid->getData();
    size_t const* ptrInd = currentIndex->getData();
    size_t const* ptrBypass = bypassOrder->getData();

    for(size_t i = 0; i < dim; ++i) {
        if(ptrInd[ptrBypass[i]] != ptrGrid[ptrBypass[i]]) {
            return i;
        }
    }

    return dim;
}

ICompact* ICompact::createCompact(IVector const* vec1, IVector const* vec2, IMultiIndex const* nodeQuantities) {
    if(!vec1 || !vec2 || !nodeQuantities) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    size_t dim = vec1->getDim();
    if(dim != vec2->getDim() || dim != nodeQuantities->getDim()) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->warning(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    IVector const* ptrVec1Copy = vec1;
    IVector const* ptrVec2Copy = vec2;

    if(!CompactImpl::isValid(ptrVec1Copy, ptrVec2Copy, nodeQuantities)) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(RC::INVALID_ARGUMENT, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    size_t allocated_size = 2 * sizeof(double) * dim + sizeof(size_t) * dim + sizeof(double) * dim;

    std::uint8_t* data = new std::uint8_t[sizeof(CompactImpl) + allocated_size];
    if(!data) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    double* ptrIncrementData = CompactImpl::axisIncrement(ptrVec1Copy, ptrVec2Copy, nodeQuantities);
    if(!ptrIncrementData) {
        delete[] data; data = nullptr;
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    ICompact* compact = new(data) CompactImpl(dim);

    size_t vector_size = sizeof(double) * dim;
    std::memcpy(data + sizeof(CompactImpl), ptrVec1Copy->getData(), vector_size);
    std::memcpy(data + sizeof(CompactImpl) + vector_size, ptrVec2Copy->getData(), vector_size);
    std::memcpy(data + sizeof(CompactImpl) + 2 * vector_size, nodeQuantities->getData(), sizeof(size_t) * dim);
    std::memcpy(data + sizeof(CompactImpl) + 2 * vector_size + sizeof(size_t) * dim, ptrIncrementData, vector_size);

    return compact;
}

ICompact* ICompact::createIntersection(ICompact const* op1, ICompact const* op2, IMultiIndex const* const grid, double tol) {
    if(!op1 || !op2 || !grid) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(RC::INVALID_ARGUMENT, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    size_t dim = op1->getDim();
    if(dim != op2->getDim() || dim != grid->getDim()) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->warning(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    if(std::isnan(tol)) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(RC::INVALID_ARGUMENT, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    IVector* left1, * left2;
    RC rc = op1->getLeftBoundary(left1);
    if(rc != RC::SUCCESS) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    rc = op2->getLeftBoundary(left2);
    if(rc != RC::SUCCESS) {
        delete left1; left1 = nullptr;
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    IVector* intersectionLeft;
    rc = CompactImpl::biggestVectorComponents(left1, left2, intersectionLeft);
    delete left1; left1 = nullptr;
    delete left2; left2 = nullptr;
    if(rc != RC::SUCCESS) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    IVector* right1, * right2;
    rc = op1->getRightBoundary(right1);
    if(rc != RC::SUCCESS) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    rc = op2->getRightBoundary(right2);
    if(rc != RC::SUCCESS) {
        delete right1; right1 = nullptr;
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    IVector* intersectionRight;
    rc = CompactImpl::leastVectorComponents(right1, right2, intersectionRight);
    delete right1; right1 = nullptr;
    delete right2; right2 = nullptr;
    if(rc != RC::SUCCESS) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    double const* intersectionLeftData = intersectionLeft->getData();
    double const* intersectionRightData = intersectionRight->getData();
    for(size_t i = 0; i < dim; ++i) {
        if(intersectionLeftData[i] - intersectionRightData[i] > tol) {
            delete intersectionLeft; intersectionLeft = nullptr;
            delete intersectionRight; intersectionRight = nullptr;
            return nullptr;
        }
        if(std::abs(intersectionRightData[i] - intersectionLeftData[i]) < tol) {
            rc = intersectionLeft->setCord(i, intersectionRightData[i]);
            if(rc != RC::SUCCESS) {
                delete intersectionLeft; intersectionLeft = nullptr;
                delete intersectionRight; intersectionRight = nullptr;
                if(CompactImpl::pLogger) {
                    CompactImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
                }
                return nullptr;
            }
        }
    }

    ICompact* intersection = ICompact::createCompact(intersectionLeft, intersectionRight, grid);
    delete intersectionLeft; intersectionLeft = nullptr;
    delete intersectionRight; intersectionRight = nullptr;
    if(!intersection) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    return intersection;
}

ICompact* ICompact::createCompactSpan(ICompact const* op1, ICompact const* op2, IMultiIndex const* const grid) {
    if(!op1 || !op2 || !grid) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    size_t dim = op1->getDim();
    if(dim != op2->getDim() || dim != grid->getDim()) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->warning(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    IVector* left1, * left2;
    RC rc = op1->getLeftBoundary(left1);
    if(rc != RC::SUCCESS) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    rc = op2->getLeftBoundary(left2);
    if(rc != RC::SUCCESS) {
        delete left1; left1 = nullptr;
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    IVector* spanLeft;
    rc = CompactImpl::leastVectorComponents(left1, left2, spanLeft);
    delete left1; left1 = nullptr;
    delete left2; left2 = nullptr;
    if(rc != RC::SUCCESS) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    IVector* right1, * right2;
    rc = op1->getRightBoundary(right1);
    if(rc != RC::SUCCESS) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    rc = op2->getRightBoundary(right2);
    if(rc != RC::SUCCESS) {
        delete right1; right1 = nullptr;
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    IVector* spanRight;
    rc = CompactImpl::biggestVectorComponents(right1, right2, spanRight);
    delete right1; right1 = nullptr;
    delete right2; right2 = nullptr;
    if(rc != RC::SUCCESS) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    ICompact* span = createCompact(spanLeft, spanRight, grid);
    delete spanLeft; spanLeft = nullptr;
    delete spanRight; spanRight = nullptr;
    if(!span) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    return span;
}

RC ICompact::setLogger(ILogger* const logger) {
    if(!logger) {
        return RC::NULLPTR_ERROR;
    }
    CompactImpl::pLogger = logger;
    return RC::SUCCESS;
}

ILogger* ICompact::getLogger() {
    return CompactImpl::pLogger;
}

inline ICompact::~ICompact() {}

CompactImpl::CompactImpl(size_t dim) : dim{ dim }, controlBlock{ new CompactControlBlockImpl(this) } {}

CompactImpl::~CompactImpl() {}

double* CompactImpl::axisIncrement(IVector const* op1, IVector const* op2, IMultiIndex const* grid) {
    size_t dim = op1->getDim();
    double const* ptrData1 = op1->getData();
    double const* ptrData2 = op2->getData();
    size_t const* ptrGrid = grid->getData();
    double* data = new double[dim];
    if(!data) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    for(size_t i = 0; i < dim; ++i) {
        data[i] = (ptrData2[i] - ptrData1[i]) / ptrGrid[i];
    }
    return data;
}

RC CompactImpl::leastVectorComponents(IVector const* op1, IVector const* op2, IVector*& res) {
    size_t dim = op1->getDim();
    double const* ptrData1 = op1->getData();
    double const* ptrData2 = op2->getData();
    double* data = new double[dim];
    if(!data) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::ALLOCATION_ERROR;
    }

    for(size_t i = 0; i < dim; ++i) {
        data[i] = std::min(ptrData1[i], ptrData2[i]);
    }

    res = IVector::createVector(dim, data);
    if(!res) {
        res = nullptr;
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::ALLOCATION_ERROR;
    }
    return RC::SUCCESS;
}

RC CompactImpl::biggestVectorComponents(IVector const* op1, IVector const* op2, IVector*& res) {
    size_t dim = op1->getDim();
    double const* ptrData1 = op1->getData();
    double const* ptrData2 = op2->getData();
    double* data = new double[dim];
    if(!data) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::ALLOCATION_ERROR;
    }

    for(size_t i = 0; i < dim; ++i) {
        data[i] = std::max(ptrData1[i], ptrData2[i]);
    }

    res = IVector::createVector(dim, data);
    if(!res) {
        res = nullptr;
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::ALLOCATION_ERROR;
    }
    return RC::SUCCESS;
}

bool CompactImpl::isValid(IVector const*& op1, IVector const*& op2, IMultiIndex const* const& grid) {
    size_t dim = op1->getDim();

    size_t const* gridData = grid->getData();
    for(size_t i = 0; i < dim; ++i) {
        if(gridData[i] == 0) {
            return false;
        }
    }

    double const* data1 = op1->getData();
    double const* data2 = op2->getData();
    bool ind = data1[0] <= data2[0];

    for(size_t i = 1; i < dim; ++i) {
        if((data1[i] <= data2[i]) != ind) {
            return false;
        }
    }

    if(!ind) {
        std::swap(op1, op2);
    }
    return true;
}

ICompact* CompactImpl::clone() const {
    size_t size = sizeAllocated();

    std::uint8_t* data = new std::uint8_t[size];
    if(!data) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    ICompact* compact = new(data) CompactImpl(dim);
    size_t instanceSize = sizeof(CompactImpl);
    std::memcpy((std::uint8_t*)data + instanceSize, (std::uint8_t*)this + instanceSize, 2 * sizeof(double) * dim + sizeof(size_t) * dim + sizeof(double) * dim);

    return compact;
}

bool CompactImpl::isInside(IVector const* const& vec) const {
    if(!vec) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return false;
    }
    if(dim != vec->getDim()) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->warning(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
        }
        return false;
    }

    double* left = (double*)((std::uint8_t*)this + sizeof(CompactImpl));
    double* right = (double*)((std::uint8_t*)this + sizeof(CompactImpl) + sizeof(double) * dim);
    double const* ptrData = vec->getData();

    for(size_t i = 0; i < dim; ++i) {
        if(left[i] > ptrData[i] || ptrData[i] > right[i]) {
            return false;
        }
    }

    return true;
}

RC CompactImpl::getVectorCopy(IMultiIndex const* index, IVector*& val) const {
    if(!index) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::NULLPTR_ERROR;
    }
    if(dim != index->getDim()) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->warning(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
        }
        return RC::MISMATCHING_DIMENSIONS;
    }

    double* data;
    RC rc = ptrVectorInCompact(index, data);
    if(rc != RC::SUCCESS) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return rc;
    }

    val = IVector::createVector(dim, data);
    delete[] data; data = nullptr;
    if(!val) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::ALLOCATION_ERROR;
    }
    return RC::SUCCESS;
}

RC CompactImpl::getVectorCoords(IMultiIndex const* index, IVector* const& val) const {
    if(!index && !val) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::NULLPTR_ERROR;
    }
    if(dim != index->getDim() || dim != val->getDim()) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->warning(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
        }
        return RC::MISMATCHING_DIMENSIONS;
    }

    double* data;
    RC rc = ptrVectorInCompact(index, data);
    if(rc != RC::SUCCESS) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return rc;
    }

    rc = val->setData(dim, data);
    if(rc != RC::SUCCESS) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return rc;
    }

    return RC::SUCCESS;
}

RC CompactImpl::ptrVectorInCompact(IMultiIndex const* index, double*& data) const {
    size_t const* ptrInd = index->getData();
    size_t const* ptrGrid = (size_t*)((std::uint8_t*)this + sizeof(CompactImpl) + 2 * sizeof(double) * dim);

    for(size_t i = 0; i < dim; ++i) {
        if(ptrInd[i] > ptrGrid[i]) {
            if(CompactImpl::pLogger) {
                CompactImpl::pLogger->sever(RC::INDEX_OUT_OF_BOUND, __FILE__, __func__, __LINE__);
            }
            return RC::INDEX_OUT_OF_BOUND;
        }
    }

    data = new double[dim];
    if(!data) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::ALLOCATION_ERROR;
    }
    double const* ptrInc = getAxisIncrement();
    std::memcpy(data, (std::uint8_t*)this + sizeof(CompactImpl), sizeof(double) * dim);
    for(size_t i = 0; i < dim; ++i) {
        data[i] += ptrInd[i] * ptrInc[i];
    }
    return RC::SUCCESS;
}

RC CompactImpl::getLeftBoundary(IVector*& vec) const {
    double* vector_data = (double*)((std::uint8_t*)this + sizeof(CompactImpl));
    vec = IVector::createVector(dim, vector_data);
    if(!vec) {
        vec = nullptr;
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::ALLOCATION_ERROR;
    }
    return RC::SUCCESS;
}

RC CompactImpl::getRightBoundary(IVector*& vec) const {
    double* vector_data = (double*)((std::uint8_t*)this + sizeof(CompactImpl) + sizeof(double) * dim);
    vec = IVector::createVector(dim, vector_data);
    if(!vec) {
        vec = nullptr;
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::ALLOCATION_ERROR;
    }
    return RC::SUCCESS;
}

size_t CompactImpl::getDim() const {
    return dim;
}

IMultiIndex* CompactImpl::getGrid() const {
    size_t* index_data = (size_t*)((std::uint8_t*)this + sizeof(CompactImpl) + 2 * sizeof(double) * dim);
    IMultiIndex* index = IMultiIndex::createMultiIndex(dim, index_data);
    if(!index) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    return index;
}

ICompact::IIterator* CompactImpl::getIterator(IMultiIndex const* const& index, IMultiIndex const* const& bypassOrder) const {
    IVector* vec;
    RC rc = getVectorCopy(index, vec);
    if(rc != RC::SUCCESS) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    IIterator* it = IteratorImpl::createIterator(vec, index, bypassOrder, controlBlock);
    delete vec; vec = nullptr;
    if(!it) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    return it;
}

ICompact::IIterator* CompactImpl::getBegin(IMultiIndex const* const& bypassOrder) const {
    size_t* index_data = new size_t[dim];
    if(!index_data) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    for(size_t i = 0; i < dim; ++i) {
        index_data[i] = 0;
    }
    IMultiIndex* multInd = IMultiIndex::createMultiIndex(dim, index_data);
    delete[] index_data; index_data = nullptr;
    if(!multInd) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    return getIterator(multInd, bypassOrder);
}

ICompact::IIterator* CompactImpl::getEnd(IMultiIndex const* const& bypassOrder) const {
    size_t* index_data = (size_t*)((std::uint8_t*)this + sizeof(CompactImpl) + 2 * sizeof(double) * dim);
    IMultiIndex* multInd = IMultiIndex::createMultiIndex(dim, index_data);
    if(!multInd) {
        if(CompactImpl::pLogger) {
            CompactImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    return getIterator(multInd, bypassOrder);
}

size_t CompactImpl::sizeAllocated() const {
    return sizeof(CompactImpl) + 2 * sizeof(double) * dim + sizeof(size_t) * dim + sizeof(double) * dim;
}

double const* CompactImpl::getAxisIncrement() const {
    return (double*)((std::uint8_t*)this + sizeof(CompactImpl) + 2 * sizeof(double) * dim + sizeof(size_t) * dim);
}
