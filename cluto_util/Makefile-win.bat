echo off

set WINTOOL="c:\Program Files\Microsoft Visual C++ Toolkit 2003"
set WINSDK="c:\Program Files\Microsoft SDK"

set CC=cl.exe
set LD=link.exe
set INCLUDES="-I%WINTOOL%\include -I%WINSDK%\include -Iinclude"
set CPPFLAGS="%INCLUDES%" /EHsc -Dsnprintf=_snprintf -D_WIN32 /Od -I.
set LIB=C:\Program Files\Microsoft Visual C++ Toolkit 2003\lib;C:\ProgramFiles\Microsoft SDK\lib;
set LIBS=libcluto.lib

set STATS_OBJS=cluto_stats.obj cluto_io.obj
set REORDER_OBJS=cluto_reorder.obj cluto_io.obj cluto_tree.obj
set GRAPH_OBJS=cluto_graphify.obj cluto_io.obj

echo compiling cpp files...
%CC% %CPPFLAGS% -c cluto_graphify.cpp
%CC% %CPPFLAGS% -c cluto_io.cpp
%CC% %CPPFLAGS% -c cluto_reorder.cpp
%CC% %CPPFLAGS% -c cluto_stats.cpp
%CC% %CPPFLAGS% -c cluto_tree.cpp


echo linking cluto_stats...
%LD% %STATS_OBJS% %LIBS% /OUT:cluto_stats.exe

echo linking cluto_reorder...
%LD% %REORDER_OBJS% %LIBS% /OUT:cluto_reorder.exe

echo linking cluto_graphify...
%LD% %GRAPH_OBJS% %LIBS% /OUT:cluto_graphify.exe

