// =================================================================== //
// Copyright (C) 2019 Kimura Ryo                                       //
//                                                                     //
// This Source Code Form is subject to the terms of the Mozilla Public //
// License, v. 2.0. If a copy of the MPL was not distributed with this //
// file, You can obtain one at http://mozilla.org/MPL/2.0/.            //
// =================================================================== //

#ifndef LIBBSDF_LOG_H
#define LIBBSDF_LOG_H

#include <iostream>

#define lbTrace lb::Log(lb::Log::Level::TRACE_MSG)
#define lbDebug lb::Log(lb::Log::Level::DEBUG_MSG)
#define lbInfo  lb::Log(lb::Log::Level::INFO_MSG)
#define lbWarn  lb::Log(lb::Log::Level::WARN_MSG)
#define lbError lb::Log(lb::Log::Level::ERROR_MSG)

namespace lb {

/*!
 * \class   Log
 * \brief   The Log class provides a simple logger.
 */
class Log
{
public:
    /*! \brief The severity level of logging. OFF is reserved for suppressing messages. */
    enum class Level {
        TRACE_MSG,
        DEBUG_MSG,
        INFO_MSG,
        WARN_MSG,
        ERROR_MSG,
        OFF_MSG
    };

    explicit Log(Level level = Level::TRACE_MSG) : level_(level) {}
    ~Log();

    /*! Gets the severity of notification. */
    static Level getNotificationLevel();

    /*! Sets the severity of notification. */
    static void setNotificationLevel(Level level);

    template <typename T>
    Log& operator<<(const T& message);

private:
    Level level_;
    static Level notificationLevel_;
};

inline Log::~Log()
{
    if (level_ >= notificationLevel_) {
        std::cout << std::endl;
    }
}

inline Log::Level Log::getNotificationLevel()
{
    return notificationLevel_;
}

inline void Log::setNotificationLevel(Level level)
{
    notificationLevel_ = level;
}

template <typename T>
Log& Log::operator<<(const T& message)
{
    if (level_ >= notificationLevel_) {
        std::cout << message;
    }
    return *this;
}

} // namespace lb

#endif // LIBBSDF_LOG_H
