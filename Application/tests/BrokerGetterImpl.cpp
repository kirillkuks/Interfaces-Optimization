#include <vector>
#include <iostream>
#include "IBrokerGetter.h"

namespace {
    class BrokerGetterImpl : public IBrokerGetter {
    public:
        BrokerGetterImpl();

        RC getProblem(char const* const& dllName, IDiffProblem*& problem) override;
        RC getSolver(char const* const& dllName, ISolver*& solver) override;

        ~BrokerGetterImpl();

    private:
         void* getBrokerData(char const* const& dllName, IBroker::INTERFACE_IMPL impl);

    private:
        std::vector<LIB> libs;

    };
}

IBrokerGetter* IBrokerGetter::createGetter() {
    return new BrokerGetterImpl();
}

IBrokerGetter::~IBrokerGetter() {}

BrokerGetterImpl::BrokerGetterImpl() : libs{ std::vector<LIB>() } {}

RC BrokerGetterImpl::getProblem(char const* const& dllName, IDiffProblem*& problem) {
    void* problemData = getBrokerData(dllName, IBroker::INTERFACE_IMPL::IPROBLEM);
    if(!problemData) {
        return RC::UNKNOWN;
    }
    problem = reinterpret_cast<IDiffProblem*>(problemData);
    return RC::SUCCESS;
}

RC BrokerGetterImpl::getSolver(char const* const& dllName, ISolver*& solver) {
    void* solverData = getBrokerData(dllName, IBroker::INTERFACE_IMPL::ISOLVER);
    if(!solverData) {
        return RC::UNKNOWN;
    }
    solver = reinterpret_cast<ISolver*>(solverData);
    return RC::SUCCESS;
}

void* BrokerGetterImpl::getBrokerData(char const* const& dllName, IBroker::INTERFACE_IMPL impl) {
    LIB hLib = LOAD_LIBRARY(dllName);
    if(!hLib) {
        std::cout << "Load library error\n";
        return nullptr;
    }

    ptr_getBroker getBroker = (ptr_getBroker)GET_BROKER(hLib);
    if(!getBroker) {
        std::cout << "Import function error\n";
        CLOSE_LIBRARY(hLib);
        return nullptr;
    }

    libs.push_back(hLib);

    IBroker* broker = reinterpret_cast<IBroker*>((*getBroker)());
    void* data = nullptr;

    if(broker->canCastTo(impl)) {
        data = broker->getInterfaceImpl(impl);
    }

    broker->release();

    return data;
}

BrokerGetterImpl::~BrokerGetterImpl() {
    for(auto& hLib : libs) {
        CLOSE_LIBRARY(hLib);
    }
}
