/*
 * Copyright (c) <2016>, <Robert Reinecke>
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef GLOBAL_FLOW_LOGGING_HPP
#define GLOBAL_FLOW_LOGGING_HPP

#include "Sinks.hpp"
#include "iostream"

// This nasty hack is only necessary because WaterGAP doesn't stick to Macro guidelines :(
#pragma push_macro("ng")
#undef ng

#include <boost/log/support/date_time.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>

#define LOG(level) BOOST_LOG_SEV(GlobalFlow::Logging::global_logger::get(), level)


namespace GlobalFlow {

#define NUM_SEVERITY_LEVELS 6

    /**
     * @enum The logging levels for the global logger
     */
    enum custom_severity_level {
        debug = 0,
        userinfo,
        stateinfo,
        numerics,
        error,
        critical
    };

    template< typename CharT, typename TraitsT >
    std::basic_ostream< CharT, TraitsT >&
    operator<< (std::basic_ostream< CharT, TraitsT >& strm, custom_severity_level lvl) {
        const char* severity_level_str[NUM_SEVERITY_LEVELS] = {
                "DEBUG",
                "USER INFO",
                "STATE INFO",
                "NUMERIC",
                "ERROR",
                "CRITICAL"
        };
        const char* str = severity_level_str[lvl];
        if (lvl < NUM_SEVERITY_LEVELS && lvl >= 0)
            strm << str;
        else
            strm << static_cast< int >(lvl);
        return strm;
    }

    namespace Logging {
        BOOST_LOG_ATTRIBUTE_KEYWORD(severity, "Severity", custom_severity_level);

        inline void initLogFile() {
            boost::gregorian::date d = boost::gregorian::day_clock::universal_day();
            std::stringstream ss;
            ss <<  "output_" << d.day() << d.month() << d.year() << ".log";
            std::string logfilename = ss.str();
            // create sink to logfile
            boost::shared_ptr<text_sink> sink = boost::make_shared<text_sink>();
            sink->locked_backend()->add_stream(boost::make_shared<std::ofstream>(logfilename.c_str()));

            // flush
            sink->locked_backend()->auto_flush(true);

            // format sink
            sink->set_formatter
                    (
                            expr::format("%1%: <%2%> %3%")
                            % expr::attr<unsigned int>("LineID")
                            % expr::attr<custom_severity_level>("Severity")
                            % expr::smessage
                    );

            //
            //sink->set_filter(severity == debug | logboost::trivial::severity == userinfo |
            //                 logboost::trivial::severity == numerics | logboost::trivial::severity == error |
            //                 logboost::trivial::severity == critical);

            // register sink
            logboost::core::get()->add_sink(sink);
        }

        inline void initCoutLog() {
            // create sink to stdout
            boost::shared_ptr<text_sink> sink = boost::make_shared<text_sink>();
            sink->locked_backend()->add_stream(boost::shared_ptr<std::ostream>(&std::clog, boost::null_deleter()));

            // flush
            sink->locked_backend()->auto_flush(true);

            // format sink
            // Currently looks like that:
            // YYYY-MM-DD HH:MI:SS: <level> A message
            sink->set_formatter(
                    expr::stream
                            << expr::format_date_time<boost::posix_time::ptime>("TimeStamp", "%Y-%m-%d %H:%M:%S")
                            << ": <" << expr::attr<custom_severity_level>("Severity")
                            << "> " << expr::smessage
            );

            // filter no logging of debug and numerics messages
            sink->set_filter(severity != custom_severity_level::debug and severity != custom_severity_level::numerics);

            // register sink
            logboost::core::get()->add_sink(sink);
        }

        // Initializing global boost::log logger
        typedef boost::log::sources::severity_channel_logger_mt<custom_severity_level, std::string> global_logger_type;
        BOOST_LOG_INLINE_GLOBAL_LOGGER_INIT(global_logger, global_logger_type) {
            boost::log::add_common_attributes();
            initLogFile();
            initCoutLog();
            return global_logger_type(boost::log::keywords::channel = "global_logger");
        }

} //ns Logging
} //ns GlobalFlow
#pragma pop_macro("ng")
#endif //LOGGING_HPP
