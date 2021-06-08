#pragma once
#include <cstddef>
#include "Interfacedllexport.h"
#include "ILogger.h"

class LIB_EXPORT IMultiIndex {
public:
    static IMultiIndex* createMultiIndex(size_t dim, size_t* indecis);

    virtual IMultiIndex* clone() const = 0;

    virtual size_t getDim() const = 0;
    virtual size_t const* getData() const = 0;
    virtual RC setData(size_t dim, size_t const* const& ptr_data) = 0;

    static RC SetLogger(ILogger* const pLogger);

    virtual RC getAxisIndex(size_t axisIndex, size_t& val) const = 0;
    virtual RC setAxisIndex(size_t axisIndex, size_t val) = 0;

    virtual RC incAxisIndex(size_t axisIndex, size_t val) = 0;

    virtual ~IMultiIndex() = 0;

private:
    IMultiIndex(IMultiIndex const&) = delete;
    IMultiIndex& operator=(IMultiIndex const&) = delete;

protected:
    IMultiIndex() = default;

};
