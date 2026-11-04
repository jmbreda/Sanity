#ifndef SANITY_LOGGING_HPP
#define SANITY_LOGGING_HPP

#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

// Prints a timestamped message to stderr.
inline void logging_debug(const std::string &msg)
{
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;

    std::tm buf;
    #ifdef _WIN32
        localtime_s(&buf, &in_time_t);
    #else
        localtime_r(&in_time_t, &buf);
    #endif
    std::cerr << "[" << std::put_time(&buf, "%Y-%m-%d %H:%M:%S");
    std::cerr << ',' << std::setw(3) << std::setfill('0') << ms.count() << "] ";
    std::cerr << msg << std::endl;
}

// Formats a duration in the largest unit that fits (days, hours, minutes or seconds), with two decimals.
inline std::string format_duration(double seconds)
{
    double value = seconds;
    const char *unit = "seconds";
    if (seconds > 86400)
    {
        value = seconds / (3600. * 24.);
        unit = "days";
    }
    else if (seconds > 3600)
    {
        value = seconds / 3600.;
        unit = "hours";
    }
    else if (seconds > 60)
    {
        value = seconds / 60.;
        unit = "minutes";
    }
    std::stringstream stream;
    stream << std::fixed << std::setprecision(2) << value << " " << unit;
    return stream.str();
}

#endif
