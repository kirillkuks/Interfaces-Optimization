#include <cmath>
#include <functional>
#include "ISolver.h"

namespace {
    class SolverImpl : public ISolver {
    using Gradient = RC (IDiffProblem::*)(IVector const* const& x, IVector* const& val) const;
    using Evaluate = double (IDiffProblem::*)(IVector const* const& x) const;
    enum class SolverBy {
        ARGS,
        PARAMS,
    };

    public:
        SolverImpl();

        ISolver* clone() const override;

        RC setProblem(IDiffProblem const* const& problem) override;

        bool isValidArgsDomain(ICompact const* const& args) const override;
        bool isValidParamsDomain(ICompact const* const& params) const override;
        RC setArgsDomain(ICompact const* const& args) override;
        RC setParamsDomain(ICompact const* const& params) override;

        RC solveByArgs(IVector const* const& initArg, IVector const* const& solverParams) override;
        RC solveByParams(IVector const* const& initParam, IVector const* const& solverParams) override;
        RC getSolution(IVector*& solution) const override;

        ~SolverImpl();

        static ILogger* pLogger;

    private:
        bool isValidSolverParams(IVector const* const& solverParams) const;
        RC getInfoForSolver(SolverBy solverBy, ICompact*& compact, Evaluate& evaluate, Gradient& gradient);
        RC gradientMethod(IVector const* const& initArgs, IVector const* const& solverParams, SolverBy solverBy);
        double calculateAlpha(Evaluate const& evaluate, IVector const* const& grad, double eps) const;
        RC solutionToDomain(ICompact const* const& compact);

    private:
        SolverImpl(IDiffProblem* const& problem, ICompact* const& argsDomain, ICompact* const& paramsDomain, ISet* const& gradSeq, IVector* const& solut);

        IDiffProblem* problem;
        ICompact* argsDomain;
        ICompact* paramsDomain;
        ISet* gradSeq;
        IVector* solution;

        static size_t const solverParamsSize = 1;

    };

}

ILogger* SolverImpl::pLogger = nullptr;

ISolver* ISolver::createSolver() {
    SolverImpl* solver = new SolverImpl();
    if(!solver) {
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }
    return solver;
}

ISolver::~ISolver() {}

SolverImpl::SolverImpl() : problem{ nullptr }, argsDomain{ nullptr }, paramsDomain{ nullptr }, gradSeq{ nullptr }, solution{ nullptr } {}

SolverImpl::SolverImpl(IDiffProblem* const& problem, ICompact* const& argsDomain, ICompact* const& paramsDomain, ISet* const& gradSeq, IVector* const& solution)
    : problem{ problem }, argsDomain{ argsDomain }, paramsDomain{ paramsDomain }, gradSeq{ gradSeq }, solution{ solution } {}

SolverImpl::~SolverImpl() {
    delete problem; problem = nullptr;
    delete argsDomain; argsDomain = nullptr;
    delete paramsDomain; paramsDomain = nullptr;
    delete gradSeq; gradSeq = nullptr;
    delete solution; solution = nullptr;
}

ISolver* SolverImpl::clone() const {
    IDiffProblem* prob = nullptr;
    if(problem) {
        prob = problem->clone();
        if(!prob) {
            if(SolverImpl::pLogger) {
                SolverImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
            }
            return nullptr;
        }
    }
    ICompact* args = nullptr;
    if(argsDomain) {
        args = argsDomain->clone();
        if(!args) {
            delete prob; prob = nullptr;
            if(SolverImpl::pLogger) {
                SolverImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
            }
            return nullptr;
        }
    }
    ICompact* params = nullptr;
    if(paramsDomain) {
        params = paramsDomain->clone();
        if(!params) {
            delete prob; prob = nullptr;
            delete args; args = nullptr;
            if(SolverImpl::pLogger) {
                SolverImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
            }
            return nullptr;
        }
    }
    ISet* seq = nullptr;
    if(gradSeq) {
        seq = gradSeq->clone();
        if(!seq) {
            delete prob; prob = nullptr;
            delete args; args = nullptr;
            delete params; params = nullptr;
            if(SolverImpl::pLogger) {
                SolverImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
            }
            return nullptr;
        }
    }
    IVector* solut = nullptr;
    if(solution) {
        solut = solut->clone();
        if(!solut) {
            delete prob; prob = nullptr;
            delete args; args = nullptr;
            delete params; params = nullptr;
            delete seq; seq = nullptr;
            if(SolverImpl::pLogger) {
                SolverImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
            }
            return nullptr;
        }
    }

    SolverImpl* solver = new SolverImpl(prob, args, params, seq, solut);
    if(!solver) {
        delete prob; prob = nullptr;
        delete args; args = nullptr;
        delete params; params = nullptr;
        delete seq; seq = nullptr;
        delete solut; solut = nullptr;
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return nullptr;
    }

    return solver;
}

RC SolverImpl::setProblem(IDiffProblem const* const& prob) {
    if(!prob) {
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::NULLPTR_ERROR;
    }
    if(problem) {
        delete problem; problem = nullptr;
    }

    problem = prob->clone();
    if(!problem) {
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::ALLOCATION_ERROR;
    }
    return RC::SUCCESS;
}

bool SolverImpl::isValidArgsDomain(ICompact const* const& args) const {
    if(!problem) {
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return false;
    }

    IVector* left, * right;
    RC rc = args->getLeftBoundary(left);
    if(rc != RC::SUCCESS) {
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return false;
    }
    rc = args->getRightBoundary(right);
    if(rc != RC::SUCCESS) {
        delete left; left = nullptr;
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return false;
    }

    bool res = problem->isValidArgs(left);
    delete left; left = nullptr;
    if(!res) {
        delete right; right = nullptr;
        return false;
    }
    res = problem->isValidArgs(right);
    delete right; right = nullptr;

    return res;
}

bool SolverImpl::isValidParamsDomain(ICompact const* const& params) const {
    if(!problem) {
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return false;
    }

    IVector* left, * right;
    RC rc = params->getLeftBoundary(left);
    if(rc != RC::SUCCESS) {
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return false;
    }
    rc = params->getRightBoundary(right);
    if(rc != RC::SUCCESS) {
        delete left; left = nullptr;
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return false;
    }

    bool res = problem->isValidParams(left);
    delete left; left = nullptr;
    if(!res) {
        delete right; right = nullptr;
        return false;
    }
    res = problem->isValidParams(right);
    delete right; right = nullptr;

    return res;
}

RC SolverImpl::setArgsDomain(ICompact const* const& args) {
    if(!args) {
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::NULLPTR_ERROR;
    }
    if(!isValidArgsDomain(args)) {
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->sever(RC::INVALID_ARGUMENT, __FILE__, __func__, __LINE__);
        }
        return RC::INVALID_ARGUMENT;
    }

    ICompact* argsCopy = args->clone();
    if(!argsCopy) {
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::ALLOCATION_ERROR;
    }
    if(argsDomain) {
        delete argsDomain; argsDomain = nullptr;
    }
    argsDomain = argsCopy;

    return RC::SUCCESS;
}

RC SolverImpl::setParamsDomain(ICompact const* const& params) {
    if(!params) {
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::NULLPTR_ERROR;
    }
    if(!isValidParamsDomain(params)) {
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->sever(RC::INVALID_ARGUMENT, __FILE__, __func__, __LINE__);
        }
        return RC::INVALID_ARGUMENT;
    }

    ICompact* paramsCopy = params->clone();
    if(!paramsCopy) {
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::ALLOCATION_ERROR;
    }
    if(paramsDomain) {
        delete paramsDomain; paramsDomain = nullptr;
    }
    paramsDomain = paramsCopy;

    return RC::SUCCESS;
}

RC SolverImpl::solveByArgs(IVector const* const& initArg, IVector const* const& solverParams) {
    if(!initArg || !solverParams) {
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::NULLPTR_ERROR;
    }
    if(!problem) {
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->sever(RC::NO_PROBLEM_SET, __FILE__, __func__, __LINE__);
        }
        return RC::NO_PROBLEM_SET;
    }
    if(!argsDomain) {
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->sever(RC::NO_ARGS_SET, __FILE__, __func__, __LINE__);
        }
        return RC::NO_ARGS_SET;
    }
    if(argsDomain->getDim() != initArg->getDim()) {
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->warning(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
        }
        return RC::MISMATCHING_DIMENSIONS;
    }
    if(solverParams->getDim() != solverParamsSize) {
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->sever(RC::INVALID_ARGUMENT, __FILE__, __func__, __LINE__);
        }
        return RC::INVALID_ARGUMENT;
    }
    if(!argsDomain->isInside(initArg) || !isValidSolverParams(solverParams)) {
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->sever(RC::INVALID_ARGUMENT, __FILE__, __func__, __LINE__);
        }
        return RC::INVALID_ARGUMENT;
    }

    if(gradSeq) {
        delete gradSeq; gradSeq = nullptr;
    }
    gradSeq = ISet::createSet();
    if(!gradSeq) {
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::ALLOCATION_ERROR;
    }

    RC rc = gradientMethod(initArg, solverParams, SolverBy::ARGS);
    if(rc != RC::SUCCESS) {
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return rc;
    }

    return RC::SUCCESS;
}

RC SolverImpl::solveByParams(IVector const* const& initParam, IVector const* const& solverParams) {
    if(!initParam || !solverParams) {
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::NULLPTR_ERROR;
    }
    if(!problem) {
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->sever(RC::NO_PROBLEM_SET, __FILE__, __func__, __LINE__);
        }
        return RC::NO_PROBLEM_SET;
    }
    if(!paramsDomain) {
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->sever(RC::NO_PARAMS_SET, __FILE__, __func__, __LINE__);
        }
        return RC::NO_PARAMS_SET;
    }
    if(paramsDomain->getDim() != initParam->getDim()) {
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->warning(RC::MISMATCHING_DIMENSIONS, __FILE__, __func__, __LINE__);
        }
        return RC::MISMATCHING_DIMENSIONS;
    }
    if(solverParams->getDim() != solverParamsSize) {
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->sever(RC::INVALID_ARGUMENT, __FILE__, __func__, __LINE__);
        }
        return RC::INVALID_ARGUMENT;
    }
    if(!paramsDomain->isInside(initParam) || !isValidSolverParams(solverParams)) {
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->sever(RC::INVALID_ARGUMENT, __FILE__, __func__, __LINE__);
        }
        return RC::INVALID_ARGUMENT;
    }

    if(gradSeq) {
        delete gradSeq; gradSeq = nullptr;
    }
    gradSeq = ISet::createSet();
    if(!gradSeq) {
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::ALLOCATION_ERROR;
    }

    RC rc = gradientMethod(initParam, solverParams, SolverBy::PARAMS);
    if(rc != RC::SUCCESS) {
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->sever(rc, __FILE__, __func__, __LINE__);
        }
        return rc;
    }

    return RC::SUCCESS;
}

RC SolverImpl::getSolution(IVector*& solut) const {
    if(!solution) {
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->sever(RC::NULLPTR_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::NULLPTR_ERROR;
    }
    solut = solution->clone();
    if(!solut) {
        if(SolverImpl::pLogger) {
            SolverImpl::pLogger->sever(RC::ALLOCATION_ERROR, __FILE__, __func__, __LINE__);
        }
        return RC::ALLOCATION_ERROR;
    }
    return RC::SUCCESS;
}

bool SolverImpl::isValidSolverParams(IVector const* const& solverParams) const {
    double const* params = solverParams->getData();
    return params[0] > 0;
}

RC SolverImpl::gradientMethod(IVector const* const& initVal, IVector const* const& solverParams, SolverBy solverBy) {
    double const* solParams = solverParams->getData();
    double eps = solParams[0];
    ICompact* compact;
    Evaluate evaluate;
    Gradient gradient;
    RC rc = getInfoForSolver(solverBy, compact, evaluate, gradient);
    if(rc != RC::SUCCESS) {
        return rc;
    }

    if(solution) {
        delete solution; solution = nullptr;
    }
    solution = initVal->clone();
    delete solution; solution = nullptr;
    solution = initVal->clone();
    if(!solution) {
        return RC::ALLOCATION_ERROR;
    }
    IVector* grad = initVal->clone();

    if(!grad) {
        delete solution; solution = nullptr;
        return RC::ALLOCATION_ERROR;
    }

    rc = (problem->*gradient)(solution, grad);
    if(rc != RC::SUCCESS) {
        delete grad; grad = nullptr;
        delete solution; solution = nullptr;
        return rc;
    }
    double norm = grad->norm(IVector::NORM::SECOND);
    if(std::isnan(norm)) {
        delete grad; grad = nullptr;
        delete solution; solution = nullptr;
        return RC::INVALID_ARGUMENT;
    }

    rc = gradSeq->insert(solution, IVector::NORM::SECOND, 0);
    if(rc != RC::SUCCESS) {
        delete grad; grad = nullptr;
        return rc;
    }

    while(norm >= eps) {
        double alpha_k = calculateAlpha(evaluate, grad, eps);
        if(std::isnan(alpha_k)) {
            return RC::INVALID_ARGUMENT;
        }

        IVector* inc = grad->clone();
        if(!inc) {
            delete grad; grad = nullptr;
            delete solution; solution = nullptr;
            return RC::ALLOCATION_ERROR;
        }
        rc = inc->scale(-alpha_k);
        if(rc != RC::SUCCESS) {
            delete inc; inc = nullptr;
            delete grad; grad = nullptr;
            delete solution; solution = nullptr;
            return rc;
        }

        IVector* xk_1 = IVector::add(solution, inc);
        delete inc; inc = nullptr;
        if(!xk_1) {
            delete grad; grad = nullptr;
            delete solution; solution = nullptr;
            return RC::ALLOCATION_ERROR;
        }

        delete solution; solution = nullptr;
        solution = xk_1->clone();
        delete solution; solution = nullptr;
        solution = xk_1->clone();
        delete xk_1;

        rc = gradSeq->insert(solution, IVector::NORM::SECOND, 0);
        if(rc != RC::SUCCESS) {
            delete grad; grad = nullptr;
            return rc;
        }

        if(!compact->isInside(solution)) {
            rc = solutionToDomain(compact);
            if(rc != RC::SUCCESS) {
                delete grad; grad = nullptr;
                delete solution; solution = nullptr;
                return rc;
            }
            break;
        }

        rc = (problem->*gradient)(solution, grad);
        if(rc != RC::SUCCESS) {
            delete grad; grad = nullptr;
            delete solution; solution = nullptr;
            return rc;
        }
        norm = grad->norm(IVector::NORM::SECOND);
        if(std::isnan(norm)) {
            delete grad; grad = nullptr;
            delete solution; solution = nullptr;
            return RC::INVALID_ARGUMENT;
        }

    }

    delete grad; grad = nullptr;

    return RC::SUCCESS;
}

double f(std::function<double(double)>const& func, double tol, double a = 0, double b = 1) {
    while(b - a > tol) {
        double delta = (b - a) / 1000;
        double x1 = (a + b) / 2 - delta;
        double x2 = (a + b) / 2 + delta;

        if(func(x2) > func(x1)) {
            b = x2;
        }
        else {
            a = x1;
        }
    }

    return b;
}

double SolverImpl::calculateAlpha(Evaluate const& evaluate, IVector const* const& grad, double eps) const {
    double a_k = f([&](double ak) -> double {
        IVector* inc = grad->clone();
        inc->scale(-ak);
        IVector* xk_1 = IVector::add(solution, inc);

        double res = (problem->*evaluate)(xk_1);
        delete xk_1; xk_1 = nullptr;

        delete inc; inc = nullptr;

        return res;
    }, eps);

    return a_k;
}

RC SolverImpl::getInfoForSolver(SolverBy solverBy, ICompact*& compact, Evaluate& evaluate, Gradient& gradient) {
    switch(solverBy) {
    case SolverBy::ARGS:
        compact = argsDomain;
        evaluate = IDiffProblem::evalByArgs;
        gradient = IDiffProblem::evalGradientByArgs;
        return RC::SUCCESS;

    case SolverBy::PARAMS:
        compact = paramsDomain;
        evaluate = IDiffProblem::evalByParams;
        gradient = IDiffProblem::evalGradientByParams;
        return RC::SUCCESS;
    }

    return RC::INVALID_ARGUMENT;
}

RC SolverImpl::solutionToDomain(ICompact const* const& compact) {
    size_t dim = solution->getDim();
    double const* solutionData = solution->getData();
    IVector* left, * right;
    RC rc = compact->getLeftBoundary(left);
    if(rc != RC::SUCCESS) {
        return rc;
    }
    rc = compact->getRightBoundary(right);
    if(rc != RC::SUCCESS) {
        delete left; left = nullptr;
        return rc;
    }

    double const* leftData = left->getData();
    double const* rightData = right->getData();

    for(size_t i = 0; i < dim; ++i) {
        if(solutionData[i] < leftData[i]) {
            rc = solution->setCord(i, leftData[i]);
            if(rc != RC::SUCCESS) {
                delete left; left = nullptr;
                delete right; right = nullptr;
                return rc;
            }
        }
        else if(solutionData[i] > rightData[i]) {
            rc = solution->setCord(i, rightData[i]);
            if(rc != RC::SUCCESS) {
                delete left; left = nullptr;
                delete right; right = nullptr;
                return rc;
            }
        }
    }

    delete left; left = nullptr;
    delete right; right = nullptr;

    return RC::SUCCESS;
}
