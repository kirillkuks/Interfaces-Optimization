#include "ILogger.h"

namespace {
    class LoggerImpl : public ILogger {
    public:
        LoggerImpl();
        LoggerImpl(char const* const& filename, bool rewrite_if_exist);

        RC log(RC code, Level level, char const* const& srcfile, char const* const& function, int line) override;
        RC log(RC code, Level level) override;

        ~LoggerImpl();

    private:
        std::ostream* os;
    };
}

ILogger* ILogger::createLogger() {
    return new LoggerImpl();
}

ILogger* ILogger::createLogger(char const* const& filename, bool rewrite_if_exist) {
    return new LoggerImpl(filename, rewrite_if_exist);
}

ILogger::~ILogger() {}

LoggerImpl::LoggerImpl() : os{ &std::cout } {}

LoggerImpl::LoggerImpl(char const* const& filename, bool rewrite_if_exist) : os{ new std::ofstream(filename) } {}

LoggerImpl::~LoggerImpl() {
    if(os != &std::cout) {
        delete os;
    }
}

RC LoggerImpl::log(RC code, Level level, char const* const& srcfile, char const* const& function, int line) {
    *os << "At file: " << srcfile << "; At line: " << line << std::endl;
    return RC::SUCCESS;
}

RC LoggerImpl::log(RC code, Level level) {
    return RC::SUCCESS;
}

