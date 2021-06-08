#include <new>
#include <cstring>
#include <cmath>
#include "IProblem.h"

namespace {
    class ProblemImpl : public IProblem {
    public:
        ProblemImpl();
        ProblemImpl(ICompact* const& params, ICompact* const& args);

        IProblem* clone() const override;

        bool isValidParams(IVector const* const& params) const override;
        bool isValidArgs(IVector const* const& args) const override;

        RC setParams(IVector const* const& params) override;
        RC setArgs(IVector const* const& args) override;

        RC setParamsDomain(ICompact const* const& params) override;
        RC setArgsDomain(ICompact const* const& args, ILogger* logger = nullptr) override;

        double evalByArgs(IVector const* const& args) const override;
        double evalByParams(IVector const* const& params) const override;

        ~ProblemImpl();
    private:
        double const* getParams() const;
        double const* getArgs() const;
        double func(double const* params, double const* args) const;

    public:
        static ILogger* pLogger;

    private:
        ICompact* paramsCompact;
        ICompact* argsCompact;

        bool isParamsInit;
        bool isArgsInit;

    public:
        static size_t const dimArgs = 2;
        static size_t const dimParams = 3;
    };
}

ILogger* ProblemImpl::pLogger = nullptr;

IProblem* IProblem::createProblem() {
    std::uint8_t* data = new std::uint8_t[sizeof(ProblemImpl) + sizeof(double) * (ProblemImpl::dimArgs + ProblemImpl::dimParams)];
    if(!data) {
        if(ProblemImpl::pLogger) {
            ProblemImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    IProblem* problem = new(data) ProblemImpl();
    return problem;
}

RC IProblem::setLogger(ILogger* const logger) {
    if(!logger) {
        return RC::NULLPTR_ERROR;
    }
    ProblemImpl::pLogger = logger;
    return RC::SUCCESS;
}

ILogger* IProblem::getLogger() {
    return ProblemImpl::pLogger;
}

IProblem::~IProblem() {}

ProblemImpl::ProblemImpl() : paramsCompact{nullptr}, argsCompact{ nullptr }, isParamsInit{ false }, isArgsInit{ false } {}

ProblemImpl::ProblemImpl(ICompact* const& params, ICompact* const& args) : paramsCompact{ params }, argsCompact{ args }, isParamsInit{ false }, isArgsInit{ false } {}

ProblemImpl::~ProblemImpl() {
    delete paramsCompact;
    delete argsCompact;
}

IProblem* ProblemImpl::clone() const {
    size_t size = sizeof(ProblemImpl) + sizeof(double) * (dimParams + dimArgs);
    std::uint8_t* data = new std::uint8_t[size];
    if(!data) {
        if(ProblemImpl::pLogger) {
            ProblemImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    std::memcpy(data, this, size);

    ICompact* paramsCopy = nullptr;
    if(paramsCompact) {
        paramsCopy = paramsCompact->clone();
        if(!paramsCopy) {
            delete[] data;
            if(ProblemImpl::pLogger) {
                ProblemImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
            }
            return nullptr;
        }
    }
    ICompact* argsCopy = nullptr;
    if(argsCompact) {
        argsCopy = argsCompact->clone();
        if(!argsCopy) {
            delete[] data;
            delete paramsCopy;
            if(ProblemImpl::pLogger) {
                ProblemImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
            }
            return nullptr;
        }
    }

    ProblemImpl* problem = (ProblemImpl*)data;
    problem->paramsCompact = paramsCopy;
    problem->argsCompact = argsCopy;

    return problem;
}

bool ProblemImpl::isValidParams(IVector const* const& params) const {
    if(!paramsCompact) {
        return false;
    }
    if(!params) {
        if(ProblemImpl::pLogger) {
            ProblemImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return false;
    }
    return paramsCompact->isInside(params);
}

bool ProblemImpl::isValidArgs(IVector const* const& args) const {
    if(!argsCompact) {
        return false;
    }
    if(!args) {
        if(ProblemImpl::pLogger) {
            ProblemImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return false;
    }
    return argsCompact->isInside(args);
}

RC ProblemImpl::setParams(IVector const* const& params) {
    if(!paramsCompact) {
        if(ProblemImpl::pLogger) {
            ProblemImpl::pLogger->sever(RC::NO_PARAMS_SET, __FILE__, __func__, __LINE__);
        }
        return RC::NO_PARAMS_SET;
    }
    if(!params) {
        if(ProblemImpl::pLogger) {
            ProblemImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::NULLPTR_ERROR;
    }
    if(!isValidParams(params)) {
        if(params->getDim() != dimParams) {
            if(ProblemImpl::pLogger) {
                ProblemImpl::pLogger->warning(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
            }
            return RC::MISMATCHING_DIMENSIONS;
        }
        if(ProblemImpl::pLogger) {
            ProblemImpl::pLogger->sever(RC::INVALID_ARGUMENT, __FILE__, __func__, __LINE__);
        }
        return RC::INVALID_ARGUMENT;
    }

    std::memcpy((std::uint8_t*)this + sizeof(ProblemImpl), params->getData(), sizeof(double) * dimParams);
    isParamsInit = true;

    return RC::SUCCESS;
}

RC ProblemImpl::setArgs(IVector const* const& args) {
    if(!argsCompact) {
        if(ProblemImpl::pLogger) {
            ProblemImpl::pLogger->sever(RC::NO_ARGS_SET, __FILE__, __func__, __LINE__);
        }
        return RC::NO_ARGS_SET;
    }
    if(!args) {
        if(ProblemImpl::pLogger) {
            ProblemImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::NULLPTR_ERROR;
    }
    if(!isValidArgs(args)) {
        if(args->getDim() != dimArgs) {
            if(ProblemImpl::pLogger) {
                ProblemImpl::pLogger->warning(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
            }
            return RC::MISMATCHING_DIMENSIONS;
        }
        if(ProblemImpl::pLogger) {
            ProblemImpl::pLogger->sever(RC::INVALID_ARGUMENT, __FILE__, __func__, __LINE__);
        }
        return RC::INVALID_ARGUMENT;
    }

    std::memcpy((std::uint8_t*)this + sizeof(ProblemImpl) + sizeof(double) * dimParams, args->getData(), sizeof(double) * dimArgs);
    isArgsInit = true;

    return RC::SUCCESS;
}

RC ProblemImpl::setParamsDomain(ICompact const* const& params) {
    if(!params) {
        if(ProblemImpl::pLogger) {
            ProblemImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::NULLPTR_ERROR;
    }
    if(ProblemImpl::dimParams != params->getDim()) {
        if(ProblemImpl::pLogger) {
            ProblemImpl::pLogger->warning(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
        }
        return RC::MISMATCHING_DIMENSIONS;
    }

    ICompact* compactCopy = params->clone();
    if(!compactCopy) {
        if(ProblemImpl::pLogger) {
            ProblemImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::ALLOCATION_ERROR;
    }
    if(paramsCompact) {
        delete paramsCompact;
    }
    paramsCompact = compactCopy;
    return RC::SUCCESS;
}

RC ProblemImpl::setArgsDomain(ICompact const* const& args, ILogger* logger) {
    if(logger) {
        ProblemImpl::pLogger = logger;
    }

    if(!args) {
        if(ProblemImpl::pLogger) {
            ProblemImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::NULLPTR_ERROR;
    }
    if(dimArgs != args->getDim()) {
        if(ProblemImpl::pLogger) {
            ProblemImpl::pLogger->warning(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
        }
        return RC::MISMATCHING_DIMENSIONS;
    }

    ICompact* compactCopy = args->clone();
    if(!compactCopy) {
        if(ProblemImpl::pLogger) {
            ProblemImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::ALLOCATION_ERROR;
    }
    if(argsCompact) {
        delete argsCompact;
    }
    argsCompact = compactCopy;
    return RC::SUCCESS;
}

double ProblemImpl::evalByArgs(IVector const* const& args) const {
    if(!args) {
        if(ProblemImpl::pLogger) {
            ProblemImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return NAN;
    }
    if(!isValidArgs(args)) {
        if(ProblemImpl::pLogger) {
            ProblemImpl::pLogger->sever(RC::INVALID_ARGUMENT, __FILE__, __func__, __LINE__);
        }
        return NAN;
    }
    if(!isParamsInit) {
        if(ProblemImpl::pLogger) {
            ProblemImpl::pLogger->sever(RC::NO_PARAMS_SET, __FILE__, __func__, __LINE__);
        }
        return NAN;
    }

    return func(getParams(), args->getData());
}

double ProblemImpl::evalByParams(IVector const* const& params) const {
    if(!params) {
        if(ProblemImpl::pLogger) {
            ProblemImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return NAN;
    }
    if(!isValidParams(params)) {
        if(ProblemImpl::pLogger) {
            ProblemImpl::pLogger->sever(RC::INVALID_ARGUMENT, __FILE__, __func__, __LINE__);
        }
        return NAN;
    }
    if(!isArgsInit) {
        if(ProblemImpl::pLogger) {
            ProblemImpl::pLogger->sever(RC::NO_ARGS_SET, __FILE__, __func__, __LINE__);
        }
        return NAN;
    }

    return func(params->getData(), getArgs());
}

double const* ProblemImpl::getParams() const {
    return (double*)((std::uint8_t*)this + sizeof(ProblemImpl));
}

double const* ProblemImpl::getArgs() const {
    return (double*)((std::uint8_t*)this + sizeof(ProblemImpl) + sizeof(double) * dimParams);
}

double ProblemImpl::func(double const* params, double const* args) const {
    double a = params[0], b = params[1], c = params[2];
    double x1 = args[0], x2 = args[1];

    // положительность подкорреного выражения гарантируется проверками до вызова
    return a * x1 + x2 + 4 * std::sqrt(1 + b * x1 * x1 + c * x2 * x2);
}
