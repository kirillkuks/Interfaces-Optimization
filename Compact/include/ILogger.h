#pragma once
#include <string>
#include <iostream>
#include <ostream>
#include <fstream>
#include "RC.h"

class __declspec(dllexport) ILogger {
public:
    enum class Level {
        SEVER,
        WARNING,
        INFO
    };

    static ILogger* createLogger();

    static ILogger* createLogger(char const* const& filename, bool rewrite_if_exist = true);

    virtual RC log(RC code, Level level, char const* const& srcfilem, char const* const& function, int line) = 0;
    virtual RC log(RC code, Level level) = 0;

    virtual RC sever(RC code, char const* const& srcfile, char const* const& function, int line) {
        return log(code, Level::SEVER, srcfile, function, line);
    }
    virtual RC sever(RC code) {
        return log(code, Level::SEVER);
    }

    virtual RC warning(RC code, char const* const& srcfile, char const* const& function, int line) {
        return log(code, Level::WARNING, srcfile, function, line);
    }
    virtual RC warning(RC code) {
        return log(code, Level::WARNING);
    }

    virtual RC info(RC code, char const* const& srcfile, char const* const& function, int line) {
        return log(code, Level::INFO, srcfile, function, line);
    }
    virtual RC info(RC code) {
        return log(code, Level::INFO);
    }

    virtual ~ILogger() = 0;

private:
    ILogger(ILogger const&);
    ILogger &operator=(ILogger const&);

protected:
    ILogger() = default;

};

