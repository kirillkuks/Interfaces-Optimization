#include <cstring>
#include <cmath>
#include "IVector.h"

namespace {
    class VectorImpl : public IVector {
    public:
        VectorImpl(size_t dim);

        IVector* clone() const override;
        double const* getData() const override;
        RC setData(size_t dim, double const* const& ptr_data) override;

        RC getCord(size_t index, double& val) const override;
        RC setCord(size_t index, double val) override;
        RC scale(double multiplier) override;
        size_t getDim() const override;

        RC inc(IVector const* const& op) override;
        RC dec(IVector const* const& op) override;

        double norm(NORM n) const override;

        RC applyFunction(const std::function<double(double)>& fun) override;
        RC foreach(const std::function<void(double)>& fun) const override;

        size_t sizeAllocated() const override;

        ~VectorImpl();

        static bool isValid(size_t dim, double const* ptr_data);
        RC cordFunction(IVector const* const op, std::function<double(double, double)> func);

        static ILogger* pLogger;

    private:
        size_t n_dim;
    };
}; //end namespace anonymous

ILogger* VectorImpl::pLogger = nullptr;

IVector::~IVector() {}

IVector* IVector::createVector(size_t dim, double const* const& ptr_data) {
    if(ptr_data && VectorImpl::isValid(dim, ptr_data)) {
        std::uint8_t* ptr = new std::uint8_t[sizeof(VectorImpl) + dim * sizeof(double)];
        IVector* vec = new(ptr) VectorImpl(dim);
        std::memcpy(ptr + sizeof(VectorImpl), ptr_data, dim * sizeof(double));
        return vec;
    }
    if(VectorImpl::pLogger) {
        VectorImpl::pLogger->sever(RC::INVALID_ARGUMENT, __FILE__, __func__, __LINE__);
    }
    return nullptr;
}

RC IVector::copyInstance(IVector* dest, IVector const* const& src) {
    if(src) {
        if(dest) {
            if(src->getDim() == dest->getDim()) {
                size_t size = src->sizeAllocated();
                if(size == dest->sizeAllocated()) {
                    if(size_t(std::abs((uint8_t*)dest - (uint8_t*)src)) > size) {
                        std::memcpy((std::uint8_t*)dest, (std::uint8_t*)src, size);
                        return RC::SUCCESS;
                    }
                    if(VectorImpl::pLogger) {
                        VectorImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
                    }
                }
                if(VectorImpl::pLogger) {
                    VectorImpl::pLogger->sever(RC::INVALID_ARGUMENT, __FILE__, __func__, __LINE__);
                }
                return RC::ALLOCATION_ERROR;
            }
            if(VectorImpl::pLogger) {
                VectorImpl::pLogger->sever(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
            }
            return RC::MISMATCHING_DIMENSIONS;
        }
        if(VectorImpl::pLogger) {
            VectorImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::NULLPTR_ERROR;
    }
    if(VectorImpl::pLogger) {
        VectorImpl::pLogger->sever(RC::INVALID_ARGUMENT, __FILE__, __func__, __LINE__);
    }
    return RC::INVALID_ARGUMENT;

}

RC IVector::moveInstance(IVector* const dest, IVector*& src) {
    if(src) {
        if(dest) {
            if(src->getDim() == dest->getDim()) {
                size_t size = src->sizeAllocated();
                if(size == dest->sizeAllocated()) {
                    std::memmove((std::uint8_t*)dest, (std::uint8_t*)src, size);
                    delete src; src = nullptr;
                    return RC::SUCCESS;
                }
                if(VectorImpl::pLogger) {
                    VectorImpl::pLogger->sever(RC::INVALID_ARGUMENT, __FILE__, __func__, __LINE__);
                }
                return RC::INVALID_ARGUMENT;
            }
            if(VectorImpl::pLogger) {
                VectorImpl::pLogger->sever(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
            }
            return RC::MISMATCHING_DIMENSIONS;
        }
        if(VectorImpl::pLogger) {
            VectorImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::NULLPTR_ERROR;
    }
    if(VectorImpl::pLogger) {
        VectorImpl::pLogger->sever(RC::INVALID_ARGUMENT, __FILE__, __func__, __LINE__);
    }
    return RC::INVALID_ARGUMENT;
}

RC IVector::setLogger(ILogger* const logger) {
    if(logger) {
        VectorImpl::pLogger = logger;
        return RC::SUCCESS;
    }
    if(VectorImpl::pLogger) {
        VectorImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
    }
    return RC::NULLPTR_ERROR;
}

IVector* IVector::add(IVector const* const& op1, IVector const* const& op2) {
    if(op1 && op2) {
        size_t size = op1->getDim();
        if(VectorImpl::isValid(size, op1->getData())) {
            if(size == op2->getDim()) {
                IVector* vec = op1->clone();
                if(vec) {
                    RC rc = vec->inc(op2);
                    if(rc == RC::SUCCESS) {
                        return vec;
                    }
                    delete vec; vec = nullptr;
                    if(VectorImpl::pLogger) {
                        VectorImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
                    }
                    return nullptr;
                }
                if(VectorImpl::pLogger) {
                    VectorImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
                }
                return nullptr;
            }
            if(VectorImpl::pLogger) {
                VectorImpl::pLogger->warning(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
            }
            return nullptr;
        }
        if(VectorImpl::pLogger) {
            VectorImpl::pLogger->sever(RC::INVALID_ARGUMENT, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    if(VectorImpl::pLogger) {
        VectorImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
    }
    return nullptr;
}

IVector* IVector::sub(IVector const* const& op1, IVector const* const& op2) {
    if(op1 && op2) {
        size_t size = op1->getDim();
        if(VectorImpl::isValid(size, op1->getData())) {
            if(size == op2->getDim()) {
                IVector* vec = op1->clone();
                if(vec) {
                    RC rc = vec->dec(op2);
                    if(rc == RC::SUCCESS) {
                        return vec;
                    }
                    delete vec; vec = nullptr;
                    if(VectorImpl::pLogger) {
                        VectorImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
                    }
                    return nullptr;
                }
                if(VectorImpl::pLogger) {
                    VectorImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
                }
                return nullptr;
            }

            if(VectorImpl::pLogger) {
                VectorImpl::pLogger->warning(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
            }
            return nullptr;
        }
        if(VectorImpl::pLogger) {
            VectorImpl::pLogger->sever(RC::INVALID_ARGUMENT, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    if(VectorImpl::pLogger) {
        VectorImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
    }
    return nullptr;
}

double IVector::dot(IVector const* const& op1, IVector const* const& op2) {
    if(op1 && op2) {
        size_t size = op1->getDim();

        if(size == op2->getDim()) {
            double dot_product = 0;
            double const* it1 = op1->getData();
            double const* it2 = op2->getData();
            for(size_t i = 0; i < size; ++i, ++it1, ++it2) {
                dot_product += *it1 + *it2;
            }
            return dot_product;
        }
        if(VectorImpl::pLogger) {
            VectorImpl::pLogger->warning(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
        }
        return NAN;
    }
    if(VectorImpl::pLogger) {
        VectorImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
    }
    return NAN;
}

bool IVector::equals(IVector const* const& op1, IVector const* const& op2, NORM n, double tol) {
    IVector* vec = IVector::sub(op1, op2);
    if(vec) {
        double norm = vec->norm(n);
        delete vec; vec = nullptr;
        return norm < tol;
    }
    return false;
}

VectorImpl::VectorImpl(size_t dim) : n_dim{ dim } {}

bool VectorImpl::isValid(size_t dim, double const* ptr_data) {
    double inf = INFINITY;
    for(size_t i = 0; i < dim; ++i) {
        double x = ptr_data[i];
        if(std::isnan(x) || x == inf || x == -inf) {
            return false;
        }
    }
    return true;
}

IVector* VectorImpl::clone() const {
    size_t size = sizeAllocated();

    std::uint8_t* vec = new std::uint8_t[size];
    if(vec) {
        std::memcpy(vec, (std::uint8_t*)this, size);
        return (IVector*)vec;
    }
    if(VectorImpl::pLogger) {
        VectorImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
    }
    return nullptr;
}

double const* VectorImpl::getData() const {
    return (double*)((std::uint8_t*)(this) + sizeof(VectorImpl));
}

RC VectorImpl::setData(size_t dim, double const* const& ptr_data) {
    if(ptr_data) {
        if(n_dim == dim)  {
            std::memcpy((uint8_t*)(this) + sizeof(VectorImpl), ptr_data, dim * sizeof(double));
            return RC::SUCCESS;
        }
        if(VectorImpl::pLogger) {
            VectorImpl::pLogger->warning(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
        }
        return RC::MISMATCHING_DIMENSIONS;
    }
    if(VectorImpl::pLogger) {
        VectorImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
    }
    return RC::NULLPTR_ERROR;
}

RC VectorImpl::getCord(size_t index, double& val) const {
    if(index < n_dim) {
        val = *(double*)((std::uint8_t*)(this) + sizeof(VectorImpl) + index * sizeof(double));
        return RC::SUCCESS;
    }
    if(VectorImpl::pLogger) {
        VectorImpl::pLogger->sever(RC::INDEX_OUT_OF_BOUND, __FILE__, __func__, __LINE__);
    }
    return RC::INDEX_OUT_OF_BOUND;
}

RC VectorImpl::setCord(size_t index, double val) {
    if(index < n_dim) {
        *(double*)((std::uint8_t*)(this) + sizeof(VectorImpl) + index * sizeof(double)) = val;
        return RC::SUCCESS;
    }
    if(VectorImpl::pLogger) {
        VectorImpl::pLogger->sever(RC::INVALID_ARGUMENT, __FILE__, __func__, __LINE__);
    }
    return RC::INDEX_OUT_OF_BOUND;
}

RC VectorImpl::scale(double multiplier) {
    if(multiplier < INFINITY) {
        size_t size = getDim();
        for(size_t i = 0; i < size; ++i) {
            *(double*)((std::uint8_t*)(this) + sizeof(VectorImpl) + i * sizeof(double)) *= multiplier;
        }
        return RC::SUCCESS;
    }
    if(VectorImpl::pLogger) {
        VectorImpl::pLogger->warning(RC::INFINITY_OVERFLOW, __FILE__, __func__, __LINE__);
    }
    return RC::INFINITY_OVERFLOW;
}

size_t VectorImpl::getDim() const {
    return n_dim;
}

RC VectorImpl::cordFunction(IVector const* const op, std::function<double(double, double)> func) {
    if(op) {
        if(isValid(op->getDim(), op->getData())) {
            size_t size = getDim();
            if(size == op->getDim()) {
                double y;
                double* it = (double*)((uint8_t*)this + sizeof(VectorImpl));
                for(size_t i = 0; i < size; ++i, ++it) {
                    if(op->getCord(i, y) == RC::SUCCESS) {
                        *it = func(*it, y);
                    }
                }
                return RC::SUCCESS;
            }
            return RC::INDEX_OUT_OF_BOUND;
        }
        else {
            return RC::INVALID_ARGUMENT;
        }
    }
    else {
        return RC::NULLPTR_ERROR;
    }
}

RC VectorImpl::inc(IVector const* const& op) {
    RC rc = cordFunction(op, [](double x, double y) -> double { return x + y; });
    if(rc != RC::SUCCESS) {
        if(VectorImpl::pLogger) {
            VectorImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
    }
    return rc;
}

RC VectorImpl::dec(IVector const* const& op) {
    RC rc = cordFunction(op, [](double x, double y) -> double { return x - y; });
    if(rc != RC::SUCCESS) {
        if(VectorImpl::pLogger) {
            VectorImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
    }
    return rc;
}

double VectorImpl::norm(NORM n) const {
    double norm = 0;
    size_t size = getDim();

    switch(n) {
    case NORM::FIRST:
        for(size_t i = 0; i < size; ++i) {
            norm += std::abs(*(double*)((std::uint8_t*)(this) + sizeof(VectorImpl) + i * sizeof(double)));
        }
        return norm;
    case NORM::SECOND:
        for(size_t i = 0; i < size; ++i) {
            double cord = *(double*)((std::uint8_t*)(this) + sizeof(VectorImpl) + i * sizeof(double));
            norm += cord * cord;
        }
        return std::sqrt(norm);
    case NORM::CHEBYSHEV:
        for(size_t i = 0; i < size; ++i) {
            double cord = std::abs(*(double*)((std::uint8_t*)(this) + sizeof(VectorImpl) + i * sizeof(double)));
            if(norm < cord) {
                norm = cord;
            }
        }
        return norm;
    default:
        return 0;
    }
    return 0;
}

RC VectorImpl::applyFunction(std::function<double(double)> const& func) {
    size_t size = getDim();
    double* it = (double*)((uint8_t*)this + sizeof(VectorImpl));
    for(size_t i = 0; i < size; ++i, ++it) {
        *it = func(*it);
    }
    return RC::SUCCESS;
}

RC VectorImpl::foreach(std::function<void(double)> const& func) const {
    size_t size = getDim();
    double *it = (double*)((uint8_t*)this +  sizeof(VectorImpl));
    for(size_t i = 0; i < size; ++i, ++it) {
        func(*it);
    }
    return RC::SUCCESS;
}

size_t VectorImpl::sizeAllocated() const {
    return sizeof(VectorImpl) + n_dim * sizeof(double);
}

VectorImpl::~VectorImpl() {}
