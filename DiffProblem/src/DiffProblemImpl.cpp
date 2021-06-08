#include <cstring>
#include <cmath>
#include "IDiffProblem.h"

namespace {
    class DiffProblemImpl : public IDiffProblem {
    public:
        DiffProblemImpl();
        DiffProblemImpl(ICompact* const& params, ICompact* const& args);

        IDiffProblem* clone() const override;

        bool isValidParams(IVector const* const& params) const override;
        bool isValidArgs(IVector const* const& args) const override;

        RC setParams(IVector const* const& params) override;
        RC setArgs(IVector const* const& args) override;

        RC setParamsDomain(ICompact const* const& params) override;
        RC setArgsDomain(ICompact const* const& args, ILogger* logger = nullptr) override;

        double evalByArgs(IVector const* const& params) const override;
        double evalByParams(IVector const* const& params) const override;

        double evalDerivativeByArgs(IVector const* const& args, IMultiIndex const* const& index) const override;
        double evalDerivativeByParams(IVector const* const& params, IMultiIndex const* const& index) const override;

        RC evalGradientByArgs(IVector const* const& args, IVector* const& val) const override;
        RC evalGradientByParams(IVector const* const& params, IVector* const& val) const override;

        ~DiffProblemImpl();

    private:
        double const* getParams() const;
        double const* getArgs() const;
        double func(double const* params, double const* args) const;

        double derivativeArgs(double const* args, size_t ind) const;
        double derivativeParams(double const* params, size_t ind) const;

        /**
        *@brief Вычисляет производную высшего порядка по одной переменной (по аргументам)
        *@param args Точка, в которой считается производная
        *@param ind Индекс переменной, по которой считается производная
        *@param order Порядок производной
        *@return Значение производной в точке args
        */
        double higherOrderDerivativesArgs(double const* const& args, size_t const& ind, size_t order) const;
        /**
        *@brief Вычисляет производную высшего порядка по одной переменной (по параметрам)
        *@param params Точка, в которой считается производная
        *@param ind Индекс переменной, по которой считается производная
        *@param order Порядок производной
        *@return Значение производной в точке params
        */
        double higherOrderDerivativesParams(double const* const& params, size_t const& ind, size_t order) const;

        /**
        *@brief Вычисляет смешанную производную по агрументам
        *
        *@param args Точка, в которой считается производная
        *@param order Порядок дифференцирования в смешанной производной
        *@param curInd Индекс переменной, по которой дифференцируем при текущем вызове
        *@param derNum Число раз, сколько уже раз было дифференцирование по переменной с индексом curInd
        *@return Значение производной в точке args
        */
        double mixedDerivarivesArgs(double const* const& args, size_t const* const& order, size_t curInd, size_t derNum = 0) const;
        /**
        *@brief Вычисляет смешанную производную по параметрам
        *
        *@param params Точка, в которой считается производная
        *@param order Порядок дифференцирования в смешанной производной
        *@param curInd Индекс переменной, по которой дифференцируем при текущем вызове
        *@param derNum Число раз, сколько уже раз было дифференцирование по переменной с индексом curInd
        *@return Значение производной в точке params
        */
        double mixedDerivarivesParams(double const* const& params, size_t const* const& order, size_t curInd, size_t  derNum = 0) const;

    public:
        static ILogger* pLogger;

    private:
        ICompact* paramsCompact;
        ICompact* argsCompact;

        bool isParamsInit;
        bool isArgsInit;

        static double constexpr h = 1e-4;

    public:
        static size_t const dimArgs = 2;
        static size_t const dimParams = 3;

    };

}

ILogger* DiffProblemImpl::pLogger = nullptr;

IDiffProblem* IDiffProblem::createDiffProblem() {
    std::uint8_t* data = new std::uint8_t[sizeof(DiffProblemImpl) + sizeof(double) * (DiffProblemImpl::dimArgs + DiffProblemImpl::dimParams)];
    if(!data) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    IDiffProblem* problem = new(data) DiffProblemImpl();
    return problem;
}

RC IDiffProblem::setLogger(ILogger* const logger) {
    if(!logger) {
        return RC::NULLPTR_ERROR;
    }
    DiffProblemImpl::pLogger = logger;
    return RC::SUCCESS;
}

ILogger* IDiffProblem::getLogger() {
    return DiffProblemImpl::pLogger;
}

IDiffProblem::~IDiffProblem() {}

DiffProblemImpl::DiffProblemImpl() : paramsCompact{ nullptr }, argsCompact{ nullptr }, isParamsInit{ false }, isArgsInit{ false } {}

DiffProblemImpl::DiffProblemImpl(ICompact* const& params, ICompact* const& args) : paramsCompact{ params }, argsCompact{ args }, isParamsInit{ false }, isArgsInit{ false } {}

DiffProblemImpl::~DiffProblemImpl() {
    delete paramsCompact; paramsCompact = nullptr;
    delete argsCompact; argsCompact = nullptr;
}

IDiffProblem* DiffProblemImpl::clone() const {
    size_t size = sizeof(DiffProblemImpl) + sizeof(double) * (dimParams + dimArgs);
    std::uint8_t* data = new std::uint8_t[size];
    if(!data) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    std::memcpy(data, this, size);

    ICompact* paramsCopy = paramsCompact->clone();
    if(!paramsCopy) {
        delete[] data; data = nullptr;
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    ICompact* argsCopy = argsCompact->clone();
    if(!argsCopy) {
        delete[] data; data = nullptr;
        delete paramsCopy; paramsCopy = nullptr;
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    DiffProblemImpl* problem = (DiffProblemImpl*)data;
    problem->paramsCompact = paramsCopy;
    problem->argsCompact = argsCopy;

    return problem;
}

bool DiffProblemImpl::isValidParams(IVector const* const& params) const {
    if(!params) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return false;
    }
    return paramsCompact->isInside(params);
}

bool DiffProblemImpl::isValidArgs(IVector const* const& args) const {
    if(!args) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return false;
    }
    return argsCompact->isInside(args);
}

RC DiffProblemImpl::setParams(IVector const* const& params) {
    if(!paramsCompact) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::NO_PARAMS_SET, __FILE__, __func__, __LINE__);
        }
        return RC::NO_PARAMS_SET;
    }
    if(!params) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::NULLPTR_ERROR;
    }
    if(!isValidParams(params)) {
        if(params->getDim() != dimParams) {
            if(DiffProblemImpl::pLogger) {
                DiffProblemImpl::pLogger->warning(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
            }
            return RC::MISMATCHING_DIMENSIONS;
        }
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::INVALID_ARGUMENT, __FILE__, __func__, __LINE__);
        }
        return RC::INVALID_ARGUMENT;
    }

    std::memcpy((std::uint8_t*)this + sizeof(DiffProblemImpl), params->getData(), sizeof(double) * dimParams);
    isParamsInit = true;

    return RC::SUCCESS;
}

RC DiffProblemImpl::setArgs(IVector const* const& args) {
    if(!argsCompact) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::NO_ARGS_SET, __FILE__, __func__, __LINE__);
        }
        return RC::NO_ARGS_SET;
    }
    if(!args) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::NULLPTR_ERROR;
    }
    if(!isValidArgs(args)) {
        if(args->getDim() != dimArgs) {
            if(DiffProblemImpl::pLogger) {
                DiffProblemImpl::pLogger->warning(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
            }
            return RC::MISMATCHING_DIMENSIONS;
        }
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::INVALID_ARGUMENT, __FILE__, __func__, __LINE__);
        }
        return RC::INVALID_ARGUMENT;
    }

    std::memcpy((std::uint8_t*)this + sizeof(DiffProblemImpl) + sizeof(double) * dimParams, args->getData(), sizeof(double) * dimArgs);
    isArgsInit = true;

    return RC::SUCCESS;
}

RC DiffProblemImpl::setParamsDomain(ICompact const* const& params) {
    if(!params) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::NULLPTR_ERROR;
    }
    if(dimParams != params->getDim()) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->warning(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
        }
        return RC::MISMATCHING_DIMENSIONS;
    }

    ICompact* compactCopy = params->clone();
    if(!compactCopy) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::ALLOCATION_ERROR;
    }
    if(paramsCompact) {
        delete paramsCompact; paramsCompact = nullptr;
    }
    paramsCompact = compactCopy;
    return RC::SUCCESS;
}

RC DiffProblemImpl::setArgsDomain(ICompact const* const& args, ILogger* logger) {
    if(logger) {
        DiffProblemImpl::pLogger = logger;
    }

    if(!args) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::NULLPTR_ERROR;
    }
    if(dimArgs != args->getDim()) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->warning(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
        }
        return RC::MISMATCHING_DIMENSIONS;
    }

    ICompact* compactCopy = args->clone();
    if(!compactCopy) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::ALLOCATION_ERROR;
    }
    if(argsCompact) {
        delete argsCompact; argsCompact = nullptr;
    }
    argsCompact = compactCopy;
    return RC::SUCCESS;
}

double DiffProblemImpl::evalByArgs(IVector const* const& args) const {
    if(!args) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return NAN;
    }
    if(!isValidArgs(args)) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::INVALID_ARGUMENT, __FILE__, __func__, __LINE__);
        }
        return NAN;
    }
    if(!isParamsInit) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::NO_PARAMS_SET, __FILE__, __func__, __LINE__);
        }
        return NAN;
    }

    return func(getParams(), args->getData());
}

double DiffProblemImpl::evalByParams(IVector const* const& params) const {
    if(!params) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return NAN;
    }
    if(!isValidParams(params)) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::INVALID_ARGUMENT, __FILE__, __func__, __LINE__);
        }
        return NAN;
    }
    if(!isArgsInit) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::NO_ARGS_SET, __FILE__, __func__, __LINE__);
        }
        return NAN;
    }

    return func(params->getData(), getArgs());
}

double DiffProblemImpl::evalDerivativeByArgs(IVector const* const& args, IMultiIndex const* const& index) const {
    if(!args || !index) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return NAN;
    }
    if(dimArgs != args->getDim() || dimArgs != index->getDim()) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->warning(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
        }
        return NAN;
    }
    if(!isValidArgs(args)) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::INVALID_ARGUMENT, __FILE__, __func__, __LINE__);
        }
        return NAN;
    }
    if(!isParamsInit) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::NO_PARAMS_SET, __FILE__, __func__, __LINE__);
        }
        return NAN;
    }

    return mixedDerivarivesArgs(args->getData(), index->getData(), dimArgs - 1);
}

double DiffProblemImpl::evalDerivativeByParams(IVector const* const& params, IMultiIndex const* const& index) const {
    if(!params || !index) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return NAN;
    }
    if(dimParams != params->getDim() || dimParams != index->getDim()) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->warning(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
        }
        return NAN;
    }
    if(!isValidParams(params)) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::INVALID_ARGUMENT, __FILE__, __func__, __LINE__);
        }
        return NAN;
    }
    if(!isArgsInit) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::NO_ARGS_SET, __FILE__, __func__, __LINE__);
        }
        return NAN;
    }

    return mixedDerivarivesParams(params->getData(), index->getData(), dimParams - 1);
}

RC DiffProblemImpl::evalGradientByArgs(IVector const* const& args, IVector* const& val) const {
    if(!args || !val) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::NULLPTR_ERROR;
    }
    if(!isParamsInit) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::NO_PARAMS_SET, __FILE__, __func__, __LINE__);
        }
        return RC::NO_PARAMS_SET;
    }
    if(dimArgs != args->getDim() || dimArgs != val->getDim()) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->warning(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
        }
        return RC::MISMATCHING_DIMENSIONS;
    }
    if(!isValidArgs(args)) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::INVALID_ARGUMENT, __FILE__, __func__, __LINE__);
        }
        return RC::INVALID_ARGUMENT;
    }

    double const* argsData = args->getData();

    for(size_t i = 0; i < dimArgs; ++i) {
        double derValue = derivativeArgs(argsData, i);
        if(std::isnan(derValue)) {
            if(DiffProblemImpl::pLogger) {
                DiffProblemImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
            }
           return RC::ALLOCATION_ERROR;
        }
        RC rc = val->setCord(i, derValue);
        if(rc != RC::SUCCESS) {
            if(DiffProblemImpl::pLogger) {
                DiffProblemImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
            }
            return rc;
        }
    }

    return RC::SUCCESS;
}

RC DiffProblemImpl::evalGradientByParams(IVector const* const& params, IVector* const& val) const {
    if(!params || !val) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::NULLPTR_ERROR;
    }
    if(!isArgsInit) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::NO_ARGS_SET, __FILE__, __func__, __LINE__);
        }
        return RC::NO_ARGS_SET;
    }
    if(dimParams != params->getDim() || dimParams != val->getDim()) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->warning(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
        }
        return RC::MISMATCHING_DIMENSIONS;
    }
    if(!isValidParams(params)) {
        if(DiffProblemImpl::pLogger) {
            DiffProblemImpl::pLogger->sever(RC::INVALID_ARGUMENT, __FILE__, __func__, __LINE__);
        }
        return RC::INVALID_ARGUMENT;
    }
    double const* paramsData = params->getData();
    for(size_t i = 0; i < dimParams; ++i) {
        double derValue = derivativeParams(paramsData, i);
        if(std::isnan(derValue)) {
            if(DiffProblemImpl::pLogger) {
                DiffProblemImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
            }
            return RC::ALLOCATION_ERROR;
        }
        RC rc = val->setCord(i, derValue);
        if(rc != RC::SUCCESS) {
            if(DiffProblemImpl::pLogger) {
                DiffProblemImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
            }
            return rc;
        }
    }
    return RC::SUCCESS;
}

double const* DiffProblemImpl::getParams() const {
    return (double*)((std::uint8_t*)this + sizeof(DiffProblemImpl));
}

double const* DiffProblemImpl::getArgs() const {
    return (double*)((std::uint8_t*)this + sizeof(DiffProblemImpl) + sizeof(double) * dimParams);
}

double DiffProblemImpl::func(double const* params, double const* args) const {
    double a = params[0], b = params[1], c = params[2];
    double x1 = args[0], x2 = args[1];

    // положительность подкорреного выражения гарантируется проверками до вызова
    return a * x1 + x2 + 4 * std::sqrt(1 + b * x1 * x1 + c * x2 * x2);
}

double DiffProblemImpl::derivativeArgs(double const* args, size_t ind) const {
    double* incArgs = new double[dimArgs];
    if(!incArgs) {
        return NAN;
    }

    std::memcpy(incArgs, args, sizeof(double) * dimArgs);
    incArgs[ind] += h;
    double const* params = getParams();

    double res = (func(params, incArgs) - func(params, args)) / h;

    delete[] incArgs; incArgs = nullptr;

    return res;
}

double DiffProblemImpl::derivativeParams(double const* params, size_t ind) const {
    double* incParams = new double[dimParams];
    if(!incParams) {
        return NAN;
    }

    std::memcpy(incParams, params, sizeof(double) * dimParams);
    incParams[ind] += h;
    double const* args = getArgs();

    double res = (func(incParams, args) - func(params, args)) / h;

    delete[] incParams; incParams = nullptr;

    return res;
}

double DiffProblemImpl::higherOrderDerivativesArgs(double const* const& args, size_t const& ind, size_t order) const {
    if(order == 1) {
        return derivativeArgs(args, ind);
    }

    double* incArgs = new double[dimArgs];
    if(!incArgs) {
        return NAN;
    }

    std::memcpy(incArgs, args, sizeof(double) * dimArgs);
    incArgs[ind] += h;

    double res = (higherOrderDerivativesArgs(incArgs, ind, order - 1) - higherOrderDerivativesArgs(args, ind, order - 1)) / h;

    delete[] incArgs; incArgs = nullptr;

    return res;
}

double DiffProblemImpl::higherOrderDerivativesParams(double const* const& params, size_t const& ind, size_t order) const {
    if(order == 1) {
        return derivativeParams(params, ind);
    }

    double* incParams = new double[dimParams];
    if(!incParams) {
        return NAN;
    }

    std::memcpy(incParams, params, sizeof(double) * dimParams);
    incParams[ind] += h;

    double res = (higherOrderDerivativesParams(incParams, ind, order - 1) - higherOrderDerivativesParams(params, ind, order - 1)) / h;

    delete[] incParams; incParams = nullptr;

    return res;
}

double DiffProblemImpl::mixedDerivarivesArgs(double const* const& args, size_t const* const& order, size_t curInd, size_t derNum) const {
    if(curInd == 0) {
        if(order[curInd] == 0) {
            return func(getParams(), args);
        }
        else {
            return higherOrderDerivativesArgs(args, curInd, order[curInd]);
        }
    }

    if(derNum == order[curInd]) {
        return mixedDerivarivesArgs(args, order, curInd - 1);
    }

    double* incArgs = new double[dimArgs];
    if(!incArgs) {
        return NAN;
    }

    std::memcpy(incArgs, args, sizeof(double) * dimArgs);
    incArgs[curInd] += h;

    double res = (mixedDerivarivesArgs(incArgs, order, curInd, derNum + 1) - mixedDerivarivesArgs(args, order, curInd, derNum + 1)) / h;

    delete[] incArgs; incArgs = nullptr;

    return res;
}

double DiffProblemImpl::mixedDerivarivesParams(double const* const& params, size_t const* const& order, size_t curInd, size_t derNum) const {
    if(curInd == 0) {
        if(order[curInd] == 0) {
            return func(params, getArgs());
        }
        else {
            return higherOrderDerivativesParams(params, curInd, order[curInd]);
        }
    }

    if(derNum == order[curInd]) {
        return mixedDerivarivesParams(params, order, curInd - 1);
    }

    double* incParams = new double[dimParams];
    if(!incParams) {
        return NAN;
    }

    std::memcpy(incParams, params, sizeof(double) * dimParams);
    incParams[curInd] += h;

    double res = (mixedDerivarivesParams(incParams, order, curInd, derNum + 1) - mixedDerivarivesParams(params, order, curInd, derNum + 1)) / h;

    delete[] incParams; incParams = nullptr;

    return res;
}
