/*******************************************************************************
 SINDER - Matt Rasmussen
 common.cpp
 11/19/05

*******************************************************************************/

#include <sys/times.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string>
#include <vector>
#include "common.h"


namespace Sinder {

using namespace std;


vector<string> g_logStages;
vector<clock_t> g_logTimes;


void Error(const char *fmt, ...)
{
   va_list ap;   
   va_start(ap, fmt);
   
   fprintf(stderr, "error: ");
   vfprintf(stderr, fmt, ap);
   fprintf(stderr, "\n");
   
   va_end(ap);
}

void Log(int level, const char *fmt, ...)
{
   va_list ap;   
   va_start(ap, fmt);
   
   for (unsigned int i=0; i<g_logStages.size(); i++)
      fprintf(stderr, "  ");
   
   vfprintf(stderr, fmt, ap);
   fprintf(stderr, "\n");
   
   va_end(ap);
}

void PushStage(const char *stage)
{
   for (unsigned int i=0; i<g_logStages.size(); i++)
      fprintf(stderr, "  ");
   
   fprintf(stderr, "BEGIN: %s\n", stage);
   g_logStages.push_back(string(stage));
   
   // timing
   struct tms buf;
   times(&buf);
   g_logTimes.push_back(buf.tms_utime);
}

float PopStage()
{
   for (unsigned int i=0; i<g_logStages.size()-1; i++)
      fprintf(stderr, "  ");
   
   struct tms buf;
   times(&buf);
   float elapse = float(buf.tms_utime - g_logTimes.back()) / sysconf(_SC_CLK_TCK);
      
   fprintf(stderr, "END:   %s [%.3f s]\n", 
           g_logStages.back().c_str(), elapse);
   g_logStages.pop_back();
   g_logTimes.pop_back();
   
   return elapse;
}


};
