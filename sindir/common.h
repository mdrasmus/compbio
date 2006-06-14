/*******************************************************************************
 SINDER - Matt Rasmussen
 common.h
 11/19/05

*******************************************************************************/

namespace Sinder {


#define max(a,b)  (((a)>(b)) ? (a) : (b))
#define min(a,b)  (((a)<(b)) ? (a) : (b))


void Error(const char *fmt, ...);
void Log(int level, const char *fmt, ...);
void PushStage(const char *stage);
float PopStage();



};
