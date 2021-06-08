#include <iostream>
#include <windows.h>
#include "IDiffProblem.h"
#include "ISolver.h"
#include "IBroker.h"
#include "IBrokerGetter.h"

int main(int argc, char* argv[]) {
    if(argc != 3) {
        return 1;
    }

    IBrokerGetter* brokerGetter = IBrokerGetter::createGetter();
    IDiffProblem* problem = nullptr;
    if(brokerGetter->getProblem(argv[1], problem) != RC::SUCCESS) {
        std::cout << "!!!\n";
    }
    ISolver* solver = nullptr;
    if(brokerGetter->getSolver(argv[2], solver) != RC::SUCCESS) {
        std::cout << "!!!\n";
    }
    if(!problem) {
        std::cout << "nullptr problem\n";
    }
    if(!solver) {
        std::cout << "nullptr solver\n";
    }

    size_t dimParams = 3, dimArgs = 2;

    size_t* paramsGridData = new size_t[dimParams];
    size_t* argsGridData = new size_t[dimArgs];
    paramsGridData[0] = 20; paramsGridData[1] = 20; paramsGridData[2] = 20;
    argsGridData[0] = 10; argsGridData[1] = 10;

    IMultiIndex* paramsGrid = IMultiIndex::createMultiIndex(dimParams, paramsGridData);
    IMultiIndex* argsGrid = IMultiIndex::createMultiIndex(dimArgs, argsGridData);

    delete[] paramsGridData; paramsGridData = nullptr;
    delete[] argsGridData; argsGridData = nullptr;

    double* paramsLeftData = new double[dimParams];
    double* paramsRightData = new double[dimParams];
    paramsLeftData[0] = 0; paramsLeftData[1] = 0; paramsLeftData[2] = 0;
    paramsRightData[0] = 10; paramsRightData[1] = 10; paramsRightData[2] = 10;

    IVector* paramsLeft = IVector::createVector(dimParams, paramsLeftData);
    IVector* paramsRight = IVector::createVector(dimParams, paramsRightData);

    delete[] paramsLeftData; paramsLeftData = nullptr;
    delete[] paramsRightData; paramsRightData = nullptr;

    ICompact* paramsCompact = ICompact::createCompact(paramsLeft, paramsRight, paramsGrid);
    if(!paramsCompact) {
        std::cout << __LINE__ << "!!!\n";
    }

    delete paramsLeft; paramsLeft = nullptr;
    delete paramsRight; paramsRight = nullptr;

    double* argsLeftData = new double[dimArgs];
    double* argsRightData = new double[dimArgs];
    argsLeftData[0] = -10; argsLeftData[1] = -10;
    argsRightData[0] = 10; argsRightData[1] = 10;

    IVector* argsLeft = IVector::createVector(dimArgs, argsLeftData);
    IVector* argsRight = IVector::createVector(dimArgs, argsRightData);

    delete[] argsLeftData; argsLeftData = nullptr;
    delete[] argsRightData; argsRightData = nullptr;

    ICompact* argsCompact = ICompact::createCompact(argsLeft, argsRight, argsGrid);
    if(!argsCompact) {
        std::cout << __LINE__ << "!!!\n";
    }

    delete argsLeft; argsLeft = nullptr;
    delete argsRight; argsRight = nullptr;

    //IDiffProblem* problem = IDiffProblem::createDiffProblem();
    if(!problem) {
        std::cout << __LINE__ << "!!!\n";
    }
    RC rc = problem->setArgsDomain(argsCompact);
    if(rc != RC::SUCCESS) {
        std::cout << __LINE__ << "!!!\n";
    }
    rc = problem->setParamsDomain(paramsCompact);
    if(rc != RC::SUCCESS) {
        std::cout << __LINE__ << "!!!\n";
    }

    delete argsCompact; argsCompact = nullptr;
    delete paramsCompact; paramsCompact = nullptr;

    double* paramsData = new double[dimParams];
    double* argsData = new double[dimArgs];
    paramsData[0] = 2; paramsData[1] = 4; paramsData[2] = 5;
    argsData[0] = 0.5; argsData[1] = 0.5;

    IVector* params = IVector::createVector(dimParams, paramsData);
    IVector* args = IVector::createVector(dimArgs, argsData);

    delete[] paramsData; paramsData = nullptr;
    delete[] argsData; argsData = nullptr;

    std::cout << "Is valid args: " << problem->isValidArgs(args) << std::endl;
    std::cout << "Is valid params: " << problem->isValidParams(params) << std::endl;

    rc = problem->setArgs(args);
    if(rc != RC::SUCCESS) {
        std::cout << __LINE__ << "!!!\n";
    }

    std::cout << "Eval by args: " << problem->evalByArgs(args) << std::endl;
    std::cout << "Eval by params: " << problem->evalByParams(params) << std::endl;

    rc = problem->setParams(params);
    if(rc != RC::SUCCESS) {
        std::cout << __LINE__ << "!!!\n";
    }

    std::cout << "Eval by args: " << problem->evalByArgs(args) << std::endl;
    std::cout << "Eval by params: " << problem->evalByParams(params) << std::endl;

    //ISolver* solver = ISolver::createSolver();
    rc = solver->setProblem(problem);
    if(rc != RC::SUCCESS) {
        std::cout << __LINE__ << "!!!\n";
    }

    IDiffProblem* problemCopy = problem->clone();
    delete problemCopy;
    problemCopy = problem->clone();

    delete problem; problem = nullptr;

    std::cout << "Eval by args(Copy): " << problemCopy->evalByArgs(args) << std::endl;
    std::cout << "Eval by params(Copy): " << problemCopy->evalByParams(params) << std::endl;

    double* leftSolverArgsData = new double[dimArgs];
    double* rightSolverArgsData = new double[dimArgs];
    leftSolverArgsData[0] = -1; leftSolverArgsData[1] = -5;
    rightSolverArgsData[0] = 7; rightSolverArgsData[1] = 5;

    IVector* leftSolverArgs = IVector::createVector(dimArgs, leftSolverArgsData);
    IVector* rightSolverArg = IVector::createVector(dimArgs, rightSolverArgsData);

    ICompact* solverArgsDomain = ICompact::createCompact(leftSolverArgs, rightSolverArg, argsGrid);
    if(!solverArgsDomain) {
        std::cout << __LINE__ << "!!!\n";
    }

    delete[] leftSolverArgsData; leftSolverArgs = nullptr;
    delete[] rightSolverArgsData; rightSolverArgsData = nullptr;

    delete leftSolverArgs; leftSolverArgs = nullptr;
    delete rightSolverArg; rightSolverArg = nullptr;

    double* leftSolverParamsData = new double[dimParams];
    double* rightSolverParamsData = new double[dimParams];
    leftSolverParamsData[0] = 0; leftSolverParamsData[1] = 0; leftSolverParamsData[2] = 0;
    rightSolverParamsData[0] = 7; rightSolverParamsData[1] = 6; rightSolverParamsData[2] = 7;

    IVector* leftSolverParams = IVector::createVector(dimParams, leftSolverParamsData);
    IVector* rightSolverParams = IVector::createVector(dimParams, rightSolverParamsData);

    delete[] leftSolverParamsData; leftSolverParamsData = nullptr;
    delete[] rightSolverParamsData; rightSolverParamsData = nullptr;

    ICompact* solverParamsDomain = ICompact::createCompact(leftSolverParams, rightSolverParams, paramsGrid);

    delete argsGrid;
    delete paramsGrid;

    delete leftSolverParams; leftSolverParams = nullptr;
    delete rightSolverParams; rightSolverParams = nullptr;

    rc = solver->setArgsDomain(solverArgsDomain);
    if(rc != RC::SUCCESS) {
        std::cout << __LINE__ << "!!!\n";
    }
    rc = solver->setParamsDomain(solverParamsDomain);
    if(rc != RC::SUCCESS) {
        std::cout << __LINE__ << "!!!\n";
    }

    delete solverArgsDomain; solverArgsDomain = nullptr;
    delete solverParamsDomain; solverParamsDomain = nullptr;

    size_t dimSolverParams = 1;
    double* solverParamsData = new double[dimSolverParams];
    solverParamsData[0] = 1E-4;
    IVector* solverParams = IVector::createVector(dimSolverParams, solverParamsData);

    delete[] solverParamsData; solverParamsData = nullptr;

    rc = solver->solveByArgs(args, solverParams);
    if(rc != RC::SUCCESS) {
        std::cout << __LINE__ << "!!!\n";
    }
    IVector* optimal;
    rc = solver->getSolution(optimal);
    if(rc != RC::SUCCESS) {
        std::cout << __LINE__ << "!!!\n";
    }
    std::cout << "Optimal by args:\n";
    double const* optimalData = optimal->getData();
    size_t dimOptimal = optimal->getDim();

    for(size_t i = 0; i < dimOptimal; ++i) {
        std::cout << optimalData[i] << " ";
    }
    std::cout << std::endl;

    delete optimal; optimal = nullptr;

    rc = solver->solveByParams(params, solverParams);
    if(rc != RC::SUCCESS) {
        std::cout << __LINE__ << "!!!\n";
    }
    rc = solver->getSolution(optimal);
    if(rc != RC::SUCCESS) {
        std::cout << __LINE__ << "!!!\n";
    }
    std::cout << "Optimal by params\n";
    optimalData = optimal->getData();
    dimOptimal = optimal->getDim();

    for(size_t i = 0; i < dimOptimal; ++i) {
        std::cout << optimalData[i] << " ";
    }
    std::cout << std::endl;

    delete optimal; optimal = nullptr;

    delete args; args = nullptr;

    delete params; params = nullptr;

    delete problemCopy; problemCopy = nullptr;

    delete solverParams; solverParams = nullptr;

    delete brokerGetter;
    return 0;
}
