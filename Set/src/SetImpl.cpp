#include <memory>
#include <cstring>
#include "ISet.h"
#include "ISetControlBlock.h"

namespace {
    class SetImpl : public ISet {
    public:
        SetImpl();

        ISet* clone() const override;

        size_t getDim() const override;
        size_t getSize() const override;

        RC getCopy(size_t index, IVector *& val) const override;
        RC findFirstAndCopy(IVector const * const& pat, IVector::NORM n, double tol, IVector *& val) const override;

        RC getCoords(size_t index, IVector * const& val) const override;
        RC findFirstAndCopyCoords(IVector const * const& pat, IVector::NORM n, double tol, IVector * const& val) const override;
        RC findFirst(IVector const* const& pat, IVector::NORM n, double tol) const override;

        RC insert(IVector const* const& val, IVector::NORM n, double tol) override;

        RC remove(size_t index) override;
        RC remove(IVector const* const& pat, IVector::NORM n, double tol) override;

        ~SetImpl();

        class IteratorImpl : public IIterator {
        public:
            static IIterator* createIterator(IVector* vector, size_t index, std::shared_ptr<ISetControlBlock> const& controlBlock);

            IIterator* getNext(size_t indexInc = 1) const override;
            IIterator* getPrevious(size_t indexInc = 1) const override;
            IIterator* clone() const override;

            RC next(size_t indexInc = 1) override;
            RC previous(size_t indexInc = 1) override;

            bool isValid() const override;

            RC makeBegin() override;
            RC makeEnd() override;

            RC getVectorCopy(IVector*& val) const override;
            RC getVectorCoords(IVector* const& val) const override;

            ~IteratorImpl();

            static ILogger* pLogger;

        private:
            IteratorImpl(IVector* vec, size_t index, std::shared_ptr<ISetControlBlock> const& controlBlock);

            IVector* vector;
            size_t index;

            std::shared_ptr<ISetControlBlock> controlBlock;
        };

        class SetControlBlockImpl : public ISetControlBlock {
        public:
            SetControlBlockImpl(SetImpl const* const& set);

            RC getNext(IVector* const& vec, size_t& index, size_t indexInc = 1) const override;
            RC getPrevious(IVector* const& vec, size_t& index, size_t indexInc = 1) const override;

            RC getBegin(IVector* const& vec, size_t& index) const override;
            RC getEnd(IVector* const& vec, size_t& index) const override;

            ~SetControlBlockImpl();

        private:
            RC indexInSet(size_t index, size_t& ind) const;

            SetImpl const* set;
        };

        IIterator *getIterator(size_t index) const override;
        IIterator *getBegin() const override;
        IIterator *getEnd() const override;

        static ILogger* pLogger;

    private:
        size_t size;
        size_t dim;
        double* data;
        size_t memSize;

        size_t curIndex;
        size_t* indexes;

        std::shared_ptr<SetControlBlockImpl> controlBlock;
    };
}

ILogger* SetImpl::pLogger = nullptr;
ILogger* SetImpl::IteratorImpl::pLogger = nullptr;

///////////////////// IIterator

ISet::IIterator::~IIterator() {}

RC ISet::IIterator::setLogger(ILogger* const logger) {
    if(!logger) {
        return RC::NULLPTR_ERROR;
    }
    SetImpl::IteratorImpl::pLogger = logger;
    return RC::SUCCESS;
}

ILogger* ISet::IIterator::getLogger() {
    return SetImpl::IteratorImpl::pLogger;
}

///////////////////// IteratorImpl

SetImpl::IteratorImpl::IteratorImpl(IVector* vec, size_t index, std::shared_ptr<ISetControlBlock> const& controlBlock) : vector{ vec }, index{ index }, controlBlock{ controlBlock } {}

SetImpl::IteratorImpl::~IteratorImpl() {}

ISet::IIterator* SetImpl::IteratorImpl::createIterator(IVector* vector, size_t index, std::shared_ptr<ISetControlBlock> const& controlBlock) {
    if(!vector) {
        if(SetImpl::IteratorImpl::pLogger) {
            SetImpl::IteratorImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    IVector* vec = vector->clone();
    if(!vec) {
        if(SetImpl::IteratorImpl::pLogger) {
            SetImpl::IteratorImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    IIterator* it = new IteratorImpl(vec, index, controlBlock);
    if(!it) {
        delete vec; vec = nullptr;
        if(SetImpl::IteratorImpl::pLogger) {
            SetImpl::IteratorImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    return it;
}

ISet::IIterator* SetImpl::IteratorImpl::getNext(size_t indexInc) const {
    IIterator* it = clone();
    if(!it) {
        if(SetImpl::IteratorImpl::pLogger) {
            SetImpl::IteratorImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    RC rc = it->next();
    if(rc != RC::SUCCESS) {
        delete it; it = nullptr;
        if(SetImpl::IteratorImpl::pLogger) {
            SetImpl::IteratorImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    return it;
}

ISet::IIterator* SetImpl::IteratorImpl::getPrevious(size_t indexInc) const {
    IIterator* it = clone();
    if(!it) {
        if(SetImpl::IteratorImpl::pLogger) {
            SetImpl::IteratorImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    RC rc = it->previous();
    if(rc != RC::SUCCESS) {
        delete it; it = nullptr;
        if(SetImpl::IteratorImpl::pLogger) {
            SetImpl::IteratorImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    return it;
}

ISet::IIterator* SetImpl::IteratorImpl::clone() const {
    return createIterator(vector, index, controlBlock);
}

RC SetImpl::IteratorImpl::next(size_t indexInc) {
    RC rc = controlBlock->getNext(vector, index, indexInc);
    if(rc != RC::SUCCESS) {
        delete vector;
        vector = nullptr;
        if(SetImpl::IteratorImpl::pLogger) {
            SetImpl::IteratorImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return rc;
    }
    return RC::SUCCESS;
}

RC SetImpl::IteratorImpl::previous(size_t indexInc) {
    RC rc = controlBlock->getPrevious(vector, index, indexInc);
    if(rc != RC::SUCCESS) {
        delete vector;
        vector = nullptr;
        if(SetImpl::IteratorImpl::pLogger) {
            SetImpl::IteratorImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return rc;
    }
    return RC::SUCCESS;
}

bool SetImpl::IteratorImpl::isValid() const {
    return vector;
}

RC SetImpl::IteratorImpl::makeBegin() {
    RC rc = controlBlock->getBegin(vector, index);
    if(rc != RC::SUCCESS) {
        delete vector;
        vector = nullptr;
        if(SetImpl::IteratorImpl::pLogger) {
            SetImpl::IteratorImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return rc;
    }
    return RC::SUCCESS;
}

RC SetImpl::IteratorImpl::makeEnd() {
    RC rc = controlBlock->getEnd(vector, index);
    if(rc != RC::SUCCESS) {
        delete vector;
        vector = nullptr;
        if(SetImpl::IteratorImpl::pLogger) {
            SetImpl::IteratorImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return rc;
    }
    return RC::SUCCESS;
}

RC SetImpl::IteratorImpl::getVectorCopy(IVector*& val) const {
    if(!vector) {
        if(SetImpl::IteratorImpl::pLogger) {
            SetImpl::IteratorImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::NULLPTR_ERROR;
    }
    val = vector->clone();
    if(!val) {
        if(SetImpl::IteratorImpl::pLogger) {
            SetImpl::IteratorImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::ALLOCATION_ERROR;
    }
    return RC::SUCCESS;
}

RC SetImpl::IteratorImpl::getVectorCoords(IVector* const& val) const {
    if(!vector || !val) {
        if(SetImpl::IteratorImpl::pLogger) {
            SetImpl::IteratorImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::NULLPTR_ERROR;
    }
    return IVector::copyInstance(val, vector);
}

///////////////////// ISetControlBlock

ISetControlBlock::~ISetControlBlock() {}

///////////////////// SetControlBlockImpl

SetImpl::SetControlBlockImpl::SetControlBlockImpl(SetImpl const* const& set) : set{ set } {}

SetImpl::SetControlBlockImpl::~SetControlBlockImpl() {}

RC SetImpl::SetControlBlockImpl::indexInSet(size_t index, size_t& ind) const {
    size_t i = 0;
    for(; i < set->size && set->indexes[i] < index; ++i) {}
    if(i == set->size) {
        return RC::INDEX_OUT_OF_BOUND;
    }
    ind = i;
    return RC::SUCCESS;
}

RC SetImpl::SetControlBlockImpl::getNext(IVector* const& vec, size_t& index, size_t indexInc) const {
    if(!vec) {
        if(pLogger) {
            pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::NULLPTR_ERROR;
    }
    size_t ind = 0;
    RC rc = indexInSet(index, ind);
    if(rc == RC::INDEX_OUT_OF_BOUND) {
        if(pLogger) {
            pLogger->warning(RC::INDEX_OUT_OF_BOUND, __FILE__, __func__, __LINE__);
        }
        return rc;
    }
    ind += (set->indexes[ind] == index ? indexInc : indexInc - 1);
    if(ind >= set->size) {
        return RC::INDEX_OUT_OF_BOUND;
    }
    index = set->indexes[ind];
    double* ptrData = set->data + set->dim * ind;
    vec->setData(set->dim, ptrData);
    return RC::SUCCESS;
}

RC SetImpl::SetControlBlockImpl::getPrevious(IVector* const& vec, size_t& index, size_t indexInc) const {
    if(!vec) {
        return RC::NULLPTR_ERROR;
    }
    size_t ind = 0;
    RC rc = indexInSet(index, ind);
    if(rc == RC::INDEX_OUT_OF_BOUND) {
        if(pLogger) {
            pLogger->warning(RC::INDEX_OUT_OF_BOUND, __FILE__, __func__, __LINE__);
        }
        return RC::INDEX_OUT_OF_BOUND;
    }
    if((int)ind - (int)indexInc < 0) {
        return RC::INDEX_OUT_OF_BOUND;
    }
    ind -= indexInc;
    index = set->indexes[ind];
    double* ptrData = set->data + set->dim * ind;
    vec->setData(set->dim, ptrData);
    return RC::SUCCESS;
}

RC SetImpl::SetControlBlockImpl::getBegin(IVector* const& vec, size_t& index) const {
    if(!vec) {
        if(pLogger) {
            pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::NULLPTR_ERROR;
    }
    index = set->indexes[0];
    double* ptrData = set->data;
    vec->setData(set->dim, ptrData);
    return RC::SUCCESS;
}

RC SetImpl::SetControlBlockImpl::getEnd(IVector* const& vec, size_t& index) const {
    if(!vec) {
        if(pLogger) {
            pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::NULLPTR_ERROR;
    }
    index = set->indexes[set->size - 1];
    double* ptrData = set->data + set->dim * (set->size - 1);
    vec->setData(set->dim, ptrData);
    return RC::SUCCESS;
}

///////////////////// ISet

ISet* ISet::createSet() {
    return new SetImpl();
}

RC ISet::setLogger(ILogger* const logger) {
    if(!logger) {
        return RC::NULLPTR_ERROR;
    }
    SetImpl::pLogger = logger;
    return RC::SUCCESS;
}

ILogger* ISet::getLogger() {
    return SetImpl::pLogger;
}

ISet* ISet::makeIntersection(ISet const* const& op1, ISet const* const& op2, IVector::NORM n, double tol) {
    if(!op1 || !op2) {
        if(SetImpl::pLogger) {
            SetImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    if(op1->getDim() != op2->getDim()) {
        if(SetImpl::pLogger) {
            SetImpl::pLogger->warning(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    size_t size1 = op1->getSize();
    size_t size2 = op2->getSize();
    if(size1 == 0 || size2 == 0) {
        return createSet();
    }

    ISet* set = createSet();
    if(!set) {
        if(SetImpl::pLogger) {
            SetImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    IVector* vec1, * vec2;
    RC rc = op1->getCopy(0, vec1);
    if(rc != RC::SUCCESS) {
        delete set; set = nullptr;
        if(SetImpl::pLogger) {
            SetImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    rc = op2->getCopy(0, vec2);
    if(rc != RC::SUCCESS) {
        delete set; set = nullptr;
        delete vec1; vec1 = nullptr;
        if(SetImpl::pLogger) {
            SetImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    for(size_t i = 0; i < size1; ++i) {
        rc = op1->getCoords(i, vec1);
        if(rc != RC::SUCCESS) {
            delete set; set = nullptr;
            delete vec1; vec1 = nullptr;
            delete vec2; vec2 = nullptr;
            if(SetImpl::pLogger) {
                SetImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
            }
            return nullptr;
        }
        for(size_t j = 0; j < size2; ++j) {
            rc = op2->getCoords(j, vec2);
            if(rc != RC::SUCCESS) {
                delete set; set = nullptr;
                delete vec1; vec1 = nullptr;
                delete vec2; vec2 = nullptr;
                if(SetImpl::pLogger) {
                    SetImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
                }
                return nullptr;
            }
            if(IVector::equals(vec1, vec2, n, tol)) {
                set->insert(vec1, n, tol);
            }
        }
    }

    delete vec1; vec1 = nullptr;
    delete vec2; vec2 = nullptr;
    return set;
}

ISet* ISet::makeUnion(ISet const* const& op1, ISet const* const& op2, IVector::NORM n, double tol) {
    if(!op1 || !op2) {
        if(SetImpl::pLogger) {
            SetImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    if(op1->getDim() != op2->getDim()) {
        if(SetImpl::pLogger) {
            SetImpl::pLogger->warning(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    size_t size2 = op2->getSize();

    ISet* set = op1->clone();
    if(!set) {
        if(SetImpl::pLogger) {
            SetImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    IVector* vec;
    if(size2 > 0) {
        RC rc = op2->getCopy(0, vec);
        if(rc != RC::SUCCESS) {
            delete set; set = nullptr;
            delete vec; vec = nullptr;
            if(SetImpl::pLogger) {
                SetImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
            }
            return nullptr;
        }
        rc = set->insert(vec, n, tol);
        if(rc != RC::SUCCESS) {
            delete set; set = nullptr;
            delete vec; vec = nullptr;
            if(SetImpl::pLogger) {
                SetImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
            }
            return nullptr;
        }
        for(size_t i = 1; i < size2; ++i) {
            rc = op2->getCoords(i, vec);
            if(rc != RC::SUCCESS) {
                delete set; set = nullptr;
                delete vec; vec = nullptr;
                if(SetImpl::pLogger) {
                    SetImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
                }
                return nullptr;
            }
            set->insert(vec, n, tol);
            if(rc != RC::SUCCESS) {
                delete set; set = nullptr;
                delete vec; vec = nullptr;
                if(SetImpl::pLogger) {
                    SetImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
                }
                return nullptr;
            }
        }

        delete vec; vec = nullptr;
    }
    return set;
}

ISet* ISet::sub(ISet const* const& op1, ISet const* const& op2, IVector::NORM n, double tol) {
    if(!op1 || !op2) {
        if(SetImpl::pLogger) {
            SetImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    if(op1->getDim() != op2->getDim()) {
        if(SetImpl::pLogger) {
            SetImpl::pLogger->warning(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    ISet* set = op1->clone();
    if(!set) {
        if(SetImpl::pLogger) {
            SetImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    size_t size = op2->getSize();
    if(size > 0) {
        IVector* vec;
        RC rc = op2->getCopy(0, vec);
        if(rc != RC::SUCCESS) {
            delete set; set = nullptr;
            delete vec; vec = nullptr;
            if(SetImpl::pLogger) {
                SetImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
            }
            return nullptr;
        }
        rc = set->remove(vec, n, tol);
        if(rc != RC::SUCCESS && rc != RC::VECTOR_NOT_FOUND) {
            delete set; set = nullptr;
            delete vec; vec = nullptr;
            if(SetImpl::pLogger) {
                SetImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
            }
            return nullptr;
        }
        for(size_t i = 1; i < size; ++i) {
            rc = op2->getCoords(i, vec);
            if(rc != RC::SUCCESS) {
                delete set; set = nullptr;
                delete vec; vec = nullptr;
                if(SetImpl::pLogger) {
                    SetImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
                }
                return nullptr;
            }
            rc = set->remove(vec, n, tol);
            if(rc != RC::SUCCESS && rc != RC::VECTOR_NOT_FOUND) {
                delete set; set = nullptr;
                delete vec; vec = nullptr;
                if(SetImpl::pLogger) {
                    SetImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
                }
                return nullptr;
            }
        }
        delete vec; vec = nullptr;
    }

    return set;
}

ISet* ISet::symSub(ISet const* const& op1, ISet const* const& op2, IVector::NORM n, double tol) {
    if(!op1 || !op2) {
        if(SetImpl::pLogger) {
            SetImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    ISet* setSub1 = sub(op1, op2, n, tol);
    if(!setSub1) {
        return nullptr;
    }
    ISet* setSub2 = sub(op2, op1, n, tol);
    if(!setSub2) {
        delete setSub1; setSub1 = nullptr;
        return nullptr;
    }

    ISet* set = makeUnion(setSub1, setSub2, n, tol);

    delete setSub1; setSub1 = nullptr;
    delete setSub2; setSub2 = nullptr;

    return set;
}

bool ISet::equals(ISet const* const& op1, ISet const* const& op2, IVector::NORM n, double tol) {
    return subSet(op1, op2, n, tol) && subSet(op2, op1, n, tol);
}

bool ISet::subSet(ISet const* const& op1, ISet const* const& op2, IVector::NORM n, double tol) {
    if(!op1 || !op2) {
        if(SetImpl::pLogger) {
            SetImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return false;
    }
    if(op1->getDim() != op2->getDim()) {
        if(SetImpl::pLogger) {
            SetImpl::pLogger->warning(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
        }
        return false;
    }
    size_t size = op1->getSize();
    if(size == 0) {
        return true;
    }
    if(op2->getSize() == 0) {
        return false;
    }

    IVector* vec;
    RC rc = op1->getCopy(0, vec);
    if(rc != RC::SUCCESS) {
        if(SetImpl::pLogger) {
            SetImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return false;
    }
    IVector* vec1;
    rc = op2->getCopy(0, vec1);
    if(rc != RC::SUCCESS) {
        delete vec; vec = nullptr;
        if(SetImpl::pLogger) {
            SetImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return false;
    }

    for(size_t i = 0; i < size; ++i) {
        rc = op1->getCoords(i, vec);
        if(rc != RC::SUCCESS) {
            delete vec; vec = nullptr;
            delete vec1; vec1 = nullptr;
            if(SetImpl::pLogger) {
                SetImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
            }
            return false;
        }
        rc = op2->findFirstAndCopyCoords(vec, n, tol, vec1);
        if(rc == RC::VECTOR_NOT_FOUND) {
            delete vec; vec = nullptr;
            delete vec1; vec1 = nullptr;
            return false;
        }
        if(rc != RC::SUCCESS) {
            delete vec; vec = nullptr;
            delete vec1; vec1 = nullptr;
            if(SetImpl::pLogger) {
                SetImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
            }
            return false;
        }
    }

    delete vec; vec = nullptr;
    delete vec1; vec1 = nullptr;
    return true;
}

ISet::~ISet() {}

///////////////////// SetImpl

SetImpl::SetImpl() : size{ 0 }, dim{ 0 }, data{ nullptr }, memSize{ 0 }, curIndex{ 0 }, indexes{ nullptr }, controlBlock{ new SetControlBlockImpl(this) } {}

ISet* SetImpl::clone() const {
    SetImpl* newSet = new SetImpl();
    if(!newSet) {
        if(SetImpl::pLogger) {
            SetImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    newSet->size = size;
    newSet->dim = dim;
    newSet->data = new double[size * dim];
    if(!newSet->data) {
        delete newSet; newSet = nullptr;
        if(SetImpl::pLogger) {
            SetImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    newSet->indexes = new size_t[size];
    if(!newSet->indexes) {
        delete[] newSet->data; newSet->data = nullptr;
        delete newSet; newSet = nullptr;
        if(SetImpl::pLogger) {
            SetImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    std::memcpy(newSet->data, data, size * dim * sizeof(double));
    for(size_t i = 0; i < size; ++i) {
        newSet->indexes[i] = i;
    }
    newSet->curIndex = size;
    newSet->memSize = size;

    return newSet;
}

size_t SetImpl::getDim() const {
    return dim;
}

size_t SetImpl::getSize() const {
    return size;
}

RC SetImpl::getCopy(size_t index, IVector *& val) const {
    if(index >= size) {
        if(SetImpl::pLogger) {
            SetImpl::pLogger->warning(RC::INDEX_OUT_OF_BOUND, __FILE__, __func__, __LINE__);
        }
        return RC::INDEX_OUT_OF_BOUND;
    }
    IVector* vec = IVector::createVector(dim, data + index * dim);
    if(!vec) {
        if(SetImpl::pLogger) {
            SetImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::ALLOCATION_ERROR;
    }
    val = vec;
    return RC::SUCCESS;
}

RC SetImpl::findFirstAndCopy(IVector const* const& pat, IVector::NORM n, double tol, IVector *& val) const {
    if(!pat) {
        if(SetImpl::pLogger) {
            SetImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::NULLPTR_ERROR;
    }
    if(size == 0) {
        if(SetImpl::pLogger) {
            SetImpl::pLogger->info(RC::VECTOR_NOT_FOUND, __FILE__, __func__, __LINE__);
        }
        return RC::VECTOR_NOT_FOUND;
    }
    double* ptrData = data;
    IVector* vec = IVector::createVector(dim, ptrData);
    size_t ind = 0;
    do {
        if(IVector::equals(vec, pat, n, tol)) {
            return getCopy(ind, val);
        }
        ++ind;
        ptrData += dim;
        vec->setData(dim, ptrData);
    } while(ind < size);
    delete vec; vec = nullptr;
    if(SetImpl::pLogger) {
        SetImpl::pLogger->info(RC::VECTOR_NOT_FOUND, __FILE__, __func__, __LINE__);
    }
    return RC::VECTOR_NOT_FOUND;
}

RC SetImpl::getCoords(size_t index, IVector* const& val) const {
    if(index >= size) {
        if(SetImpl::pLogger) {
            SetImpl::pLogger->warning(RC::VECTOR_NOT_FOUND, __FILE__, __func__, __LINE__);
        }
        return RC::INDEX_OUT_OF_BOUND;
    }
    return val->setData(dim, data + dim * index);
}

RC SetImpl::findFirstAndCopyCoords(IVector const* const& pat, IVector::NORM n, double tol, IVector* const& val) const {
    if(!pat) {
        if(SetImpl::pLogger) {
            SetImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::NULLPTR_ERROR;
    }
    if(size == 0) {
        return RC::VECTOR_NOT_FOUND;
    }
    double* ptrData = data;
    IVector* vec = IVector::createVector(dim, ptrData);
    size_t ind = 0;
    do {
        if(IVector::equals(vec, pat, n, tol)) {
            return getCoords(ind, val);
        }
        ++ind;
        ptrData += dim;
        vec->setData(dim, ptrData);
    } while(ind < size);
    delete vec; vec = nullptr;
    return RC::VECTOR_NOT_FOUND;
}

RC SetImpl::findFirst(IVector const* const& pat, IVector::NORM n, double tol) const {
    if(!pat) {
        if(SetImpl::pLogger) {
            SetImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::NULLPTR_ERROR;
    }
    if(size == 0) {
        return RC::VECTOR_NOT_FOUND;
    }
    double* ptrData = data;
    IVector* vec = IVector::createVector(dim, ptrData);
    size_t ind = 0;
    do {
        if(IVector::equals(vec, pat, n, tol)) {
            return RC::SUCCESS;
        }
        ++ind;
        ptrData += dim;
        vec->setData(dim, ptrData);
    } while(ind < size);
    delete vec; vec = nullptr;
    return RC::VECTOR_NOT_FOUND;
}

RC SetImpl::insert(IVector const* const& val, IVector::NORM n, double tol) {
    if(!val) {
        if(SetImpl::pLogger) {
            SetImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::NULLPTR_ERROR;
    }
    if(dim) {
        if(dim != val->getDim()) {
            if(SetImpl::pLogger) {
                SetImpl::pLogger->warning(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
            }
            return RC::MISMATCHING_DIMENSIONS;
        }
        IVector* vec = val->clone();
        RC rc = findFirstAndCopyCoords(val, n, tol, vec);
        if(rc != RC::VECTOR_NOT_FOUND) {
            delete vec; vec = nullptr;
            return rc;
        }
        delete vec; vec = nullptr;
        if(memSize > size) {
            std::memcpy(data + size * dim, val->getData(), dim * sizeof(double));
            indexes[size] = curIndex++;
            ++size;
        } else {
            double* data_ = new double[(size + 1) * dim];
            if(!data_) {
                if(SetImpl::pLogger) {
                    SetImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
                }
                return RC::ALLOCATION_ERROR;
            }
            size_t* indexes_ = new size_t[size + 1];
            if(!indexes_) {
                delete[] data_; data_ = nullptr;
                if(SetImpl::pLogger) {
                    SetImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
                }
                return RC::ALLOCATION_ERROR;
            }
            std::memcpy(data_, data, size * dim * sizeof(double));
            std::swap(data, data_);
            std::memcpy(data + size * dim, val->getData(), dim * sizeof(double));
            std::memcpy(indexes_, indexes, size * sizeof(size_t));
            std::swap(indexes, indexes_);
            indexes[size] = curIndex++;
            ++size;
            ++memSize;
            delete[] data_; data_ = nullptr;
            delete[] indexes_; indexes_ = nullptr;
        }
    }
    else {
        dim = val->getDim();
        memSize = size = 1;
        data = new double[size * dim];
        if(!data) {
            if(SetImpl::pLogger) {
                SetImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
            }
            return RC::ALLOCATION_ERROR;
        }
        indexes = new size_t[size];
        if(!indexes) {
            if(SetImpl::pLogger) {
                SetImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
            }
            return RC::ALLOCATION_ERROR;
        }
        indexes[size - 1] = curIndex++;

        std::memcpy(data, val->getData(), dim * sizeof(double));
    }
    return RC::SUCCESS;
}

RC SetImpl::remove(size_t index) {
    if(index >= size) {
        if(SetImpl::pLogger) {
            SetImpl::pLogger->warning(RC::INDEX_OUT_OF_BOUND, __FILE__, __func__, __LINE__);
        }
        return RC::INDEX_OUT_OF_BOUND;
    }
    double* offset = data + index * dim;
    std::memmove(offset, offset + dim, (size - index - 1) * dim * sizeof(double));
    std::memmove(indexes + index, indexes + index + 1, (size - index - 1) * sizeof(size_t));
    --size;
    return RC::SUCCESS;
}

RC SetImpl::remove(IVector const* const& pat, IVector::NORM n, double tol) {
    if(!pat) {
        if(SetImpl::pLogger) {
            SetImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::NULLPTR_ERROR;
    }
    if(size == 0) {
        if(SetImpl::pLogger) {
            SetImpl::pLogger->info(RC::VECTOR_NOT_FOUND, __FILE__, __func__, __LINE__);
        }
        return RC::VECTOR_NOT_FOUND;
    }
    double* ptrData = data;
    IVector* vec = IVector::createVector(dim, ptrData);
    size_t ind = 0;
    do {
        if(IVector::equals(vec, pat, n, tol)) {
            remove(ind);
            return RC::SUCCESS;
        }
        ++ind;
        ptrData += dim;
        vec->setData(dim, ptrData);
    } while(ind < size);

    delete vec; vec = nullptr;
    return RC::VECTOR_NOT_FOUND;
}

ISet::IIterator* SetImpl::getIterator(size_t index) const {
    IVector* vec;
    RC rc = getCopy(index, vec);
    if(rc != RC::SUCCESS) {
        if(SetImpl::pLogger) {
            SetImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    IIterator* it = IteratorImpl::createIterator(vec, indexes[index], controlBlock);
    delete vec; vec = nullptr;
    if(!it) {
        if(SetImpl::pLogger) {
            SetImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    return it;
}

ISet::IIterator* SetImpl::getBegin() const {
    if(size > 0) {
        return getIterator(0);
    }
    return nullptr;
}

ISet::IIterator* SetImpl::getEnd() const {
    if(size > 0) {
        return getIterator(size - 1);
    }
    return nullptr;
}

SetImpl::~SetImpl() {
    delete[] data; data = nullptr;
    delete[] indexes; indexes = nullptr;
}
