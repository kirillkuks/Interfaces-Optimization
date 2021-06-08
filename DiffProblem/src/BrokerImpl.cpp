#include "IBroker.h"
#include "IDiffProblem.h"
#include "Interfacedllexport.h"

namespace {
    class BrokerImpl : public IBroker {
    public:
        static BrokerImpl* instance() {
            if(!m_instance) {
                m_instance = new (std::nothrow) BrokerImpl();
            }
            return m_instance;
        }

        bool canCastTo(INTERFACE_IMPL impl) const override {
            return impl == INTERFACE_IMPL::IPROBLEM;
        }

        void* getInterfaceImpl(INTERFACE_IMPL impl) const override {
            return canCastTo(impl) ? IDiffProblem::createDiffProblem() : NULL;
        }

        void release() override {
            delete m_instance;
            m_instance = nullptr;
        }

    private:
        static BrokerImpl* m_instance;

    };

    BrokerImpl* BrokerImpl::m_instance = nullptr;

}

IBroker::~IBroker() {}

extern "C" {
    LIB_EXPORT void* getBroker() {
        return (void*)BrokerImpl::instance();
    }
}
