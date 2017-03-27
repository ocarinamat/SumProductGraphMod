/* 
 * File:   Log.h
 * Author: mdesana
 *
 * Created on May 15, 2014, 5:27 PM
 */

//this is REALLY a simple log

//we put this ifdef because we want to prevent usage of this class from the code when USE_DEBUG_LOG is off
//NOTE: the log must always be used between #ifdef USE_DEBUG_LOG ... #endif directives

#ifdef USE_DEBUG_LOG 

#pragma once

#include <fstream>
#include <sstream>
#include "Params.h"

// a minimal log represented as a singleton

class Log
{
    std::ofstream m_file;
    int m_indentationLevel;

    static Log* m_log;
    static Log* GetInstance();

    void WriteAux(const std::stringstream& ss, bool endline, int level)
    {
        if (level < ENABLED_LOG_LEVEL)
            return;

        std::string indent;
        for (int i = 0; i < m_indentationLevel; i++)
            indent += ".    ";
        
        indentedOutput(m_file, ss.str().c_str(),indent.c_str());
//        m_file << ss.str();

        if (endline)
        {
            m_file << std::endl;
        }
        m_file.flush(); //todo remove if too slow
    }

    static void indentedOutput(std::ostream &outStream, const char *message, const char *indentation)
    {
        bool newline;
        while (char cur = *message)
        {
            if (newline)
            {
                outStream << indentation;
                newline = false;
            }
            outStream << cur;
            if (cur == '\n')
            {
                newline = true;
            }
            ++message;
        }
    }

public:

    ~Log()
    {
        m_file.close();
    }

    static void Indent()
    {
        GetInstance()->m_indentationLevel++;
    }

    static void Dedent()
    {
        if (GetInstance()->m_indentationLevel > 0)
            GetInstance()->m_indentationLevel--;
        else
            std::cout << "\nWARNING trying to dedent log with no indentation\n";
    }

    static void Write(const std::stringstream& s, int level = LOG_LEVEL_LOW)
    {
        Log::GetInstance()->WriteAux(s, true, level);
    }

    static void Write(const std::string& s, int level = LOG_LEVEL_LOW)
    {
        std::stringstream sstr;
        sstr << s;
        Log::GetInstance()->WriteAux(sstr, true, level);
    }

    static void Endl(int level = LOG_LEVEL_LOW)
    {
        std::stringstream sstr;
        Log::GetInstance()->WriteAux(sstr, true, level);
    }

    static void Write(Real s, int level = LOG_LEVEL_LOW)
    {
        std::stringstream sstr;
        sstr << s;
        Log::GetInstance()->WriteAux(sstr, true, level);
    }

    static void Write(int s, int level = LOG_LEVEL_LOW)
    {
        std::stringstream sstr;
        sstr << s;
        Log::GetInstance()->WriteAux(sstr, true, level);
    }

    //        static void Write(const std::stringstream& s, int level = LOG_LEVEL_LOW)
    //    {
    //        Log::GetInstance()->WriteAux(s, false, level);
    //    }
    //
    //    static void Write(const std::string& s, int level = LOG_LEVEL_LOW)
    //    {
    //        std::stringstream sstr;
    //        sstr << s;
    //        Log::GetInstance()->WriteAux(sstr, false, level);
    //    }
    //
    //    static void Write(Real s, int level = LOG_LEVEL_LOW)
    //    {
    //        std::stringstream sstr;
    //        sstr << s;
    //        Log::GetInstance()->WriteAux(sstr, false, level);
    //    }
    //
    //    static void Write(int s, int level = LOG_LEVEL_LOW)
    //    {
    //        std::stringstream sstr;
    //        sstr << s;
    //        Log::GetInstance()->WriteAux(sstr, false, level);
    //    }
};


#endif /* USE_DEBUG_LOG */