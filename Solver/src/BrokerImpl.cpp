#include "IBroker.h"
#include "ISolver.h"
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

        bool canCastTo(INTERFACE_IMPL impl) const {
            return impl == INTERFACE_IMPL::ISOLVER;
        }

        void* getInterfaceImpl(INTERFACE_IMPL impl) const {
            return canCastTo(impl) ? ISolver::createSolver() : nullptr;
        }

        void release() {
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
