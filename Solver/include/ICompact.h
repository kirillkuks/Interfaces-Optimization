#pragma once
#include "Interfacedllexport.h"
#include "ILogger.h"
#include "IVector.h"
#include "RC.h"
#include "IMultiIndex.h"

class LIB_EXPORT ICompact {
public:
    static ICompact* createCompact(IVector const * vec1, IVector const * vec2, IMultiIndex const *nodeQuantities);

    virtual ICompact *clone() const = 0;

    static RC setLogger(ILogger* const logger);
    static ILogger* getLogger();

    virtual bool isInside(IVector const * const&vec) const = 0;
    virtual RC getVectorCopy(IMultiIndex const *index, IVector *& val) const = 0;
    virtual RC getVectorCoords(IMultiIndex const *index, IVector * const& val) const = 0;

    virtual RC getLeftBoundary(IVector *& vec) const = 0;
    virtual RC getRightBoundary(IVector *& vec) const = 0;
    virtual size_t getDim() const = 0;
    virtual IMultiIndex* getGrid() const = 0;

    static ICompact* createIntersection(ICompact const *op1, ICompact const *op2, IMultiIndex const* const grid, double tol);
    static ICompact* createCompactSpan(ICompact const *op1, ICompact const *op2, IMultiIndex const* const grid);

    class IIterator {
    public:
        virtual IIterator * getNext() = 0;
        virtual IIterator * clone() const = 0;

        static RC setLogger(ILogger * const pLogger);
        static ILogger* getLogger();

        virtual RC next() = 0;

        virtual bool isValid() const = 0;

        virtual RC getVectorCopy(IVector *& val) const = 0;
        virtual RC getVectorCoords(IVector * const& val) const = 0;

        virtual ~IIterator() = 0;

    private:
        IIterator(const IIterator&) = delete;
        IIterator& operator=(const IIterator&) = delete;

    protected:
        IIterator() = default;
    };

    virtual IIterator* getIterator(IMultiIndex const * const&index, IMultiIndex const * const &bypassOrder) const = 0;
    virtual IIterator* getBegin(IMultiIndex const * const &bypassOrder) const = 0;
    virtual IIterator* getEnd(IMultiIndex const * const &bypassOrder) const = 0;

    virtual ~ICompact() = 0;

private:
    ICompact(const ICompact& compact) = delete;
    ICompact& operator=(const ICompact& compact) = delete;

protected:
    ICompact() = default;
};
