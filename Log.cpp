# include "Common.h"

#ifdef USE_DEBUG_LOG

Log* Log::m_log = NULL;
Log* Log::GetInstance()
{
    if (m_log == NULL)
    {
        m_log = new Log();
        m_log->m_file.open("log/SPN_log.txt");
        m_log->m_indentationLevel=0;
        return m_log;
    }
    else
        return m_log;
}

#endif