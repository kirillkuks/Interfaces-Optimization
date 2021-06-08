#include <windows.h>
#include "IDiffProblem.h"
#include "ISolver.h"
#include "IBroker.h"
#include "RC.h"

typedef void* (*ptr_getBroker)();

#ifndef _WIN32
    #define LIB void*
    #define LOAD_LIBRARY(lib_name) dlopen(lib_name, RTLD_LAZY)
    #define GET_BROKER(lib) dlsym(lib, "getBroker")
    #define CLOSE_LIBRARY(lib) dlclose(lib)
#else
    #define LIB HMODULE
    #define LOAD_LIBRARY(lib_name) LoadLibrary(lib_name)
    #define GET_BROKER(lib) GetProcAddress(lib, "getBroker")
    #define CLOSE_LIBRARY(lib) FreeLibrary(lib)
#endif

class IBrokerGetter {
public:
    static IBrokerGetter* createGetter();

    virtual RC getProblem(char const* const& dllName, IDiffProblem*& problem) = 0;
    virtual RC getSolver(char const* const& dllName, ISolver*& solver) = 0;

    virtual ~IBrokerGetter() = 0;

protected:
    IBrokerGetter() = default;

private:
    IBrokerGetter(IBrokerGetter const&) = delete;
    IBrokerGetter& operator=(IBrokerGetter const&) = delete;
};
