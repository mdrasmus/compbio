#ifndef SPIDIR_CONFIG_PARAM_H
#define SPIDIR_CONFIG_PARAM_H

#include <libgen.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>


using namespace std;

enum {
    OPTION_ARG,
    OPTION_COMMENT
};


class ConfigParamBase
{
public:
    ConfigParamBase(string shortarg, string longarg, string argstr,
                    string help="") :
        kind(OPTION_ARG),
        shortarg(shortarg),
        longarg(longarg),
        argstr(argstr),
        help(help)
    {
    }
    
    virtual ~ConfigParamBase()
    {}
    
    virtual int parse(int argc, const char **argv)
    {
        return -1;
    }
    
    int kind;
    string shortarg;
    string longarg;
    string argstr;
    string help;
};


class ConfigParamComment : public ConfigParamBase
{
public:
    ConfigParamComment(string msg) :
        ConfigParamBase("", "", "", ""),
        msg(msg)
    {
        kind = OPTION_COMMENT;
    }
    
    virtual ~ConfigParamComment()
    {}
    
    virtual int parse(int argc, const char **argv)
    {
        return -1;
    }
    
    string msg;
};


template <class T>
class ConfigParam : public ConfigParamBase
{
public:
    ConfigParam(string shortarg, string longarg, string argstr,
                T *value, string help) :
        ConfigParamBase(shortarg, longarg, argstr, help),
        value(value),
        hasDefault(false)
    {
    }
    
    ConfigParam(string shortarg, string longarg, string argstr,
                T *value, T defaultValue, string help) :
        ConfigParamBase(shortarg, longarg, argstr, help),
        value(value),
        defaultValue(defaultValue),
        hasDefault(true)
    {
        *value = defaultValue;
    }    
    
    virtual int parse(int argc, const char **argv)
    {
        if (argc > 0) {
            *value = T(argv[0]);
            return 1;
        } else {
            fprintf(stderr, "error: argument value expected\n");
            return -1;
        }
    }
    
    T *value;
    T defaultValue;
    bool hasDefault;
};


class ConfigSwitch : public ConfigParamBase
{
public:
    ConfigSwitch(string shortarg, string longarg, 
                bool *value, string help) :
        ConfigParamBase(shortarg, longarg, "", help),
        value(value)
    {
        *value = false;
    }
    
    virtual int parse(int argc, const char **argv)
    {
        *value = true;
        return 0;
    }
    
    bool *value;
};


template <>
int ConfigParam<int>::parse(int argc, const char **argv)
{
    if (argc > 0) {
        if (sscanf(argv[0], "%d", value) != 1) {
            fprintf(stderr, "error: int expected '%s'\n", argv[0]);        
            return -1;
        }
        return 1;
    } else {
        fprintf(stderr, "error: argument value expected\n");    
        return -1;
    }
}


template <>
int ConfigParam<float>::parse(int argc, const char **argv)
{
    if (argc > 0) {
        if (sscanf(argv[0], "%f", value) != 1) {
            fprintf(stderr, "error: float expected '%s'\n", argv[0]);
            return -1;
        }
        return 1;
    } else {
        fprintf(stderr, "error: argument value expected\n");
        return -1;
    }
}


class ConfigParser
{
public:
    ConfigParser()
    {}
    
    ~ConfigParser()
    {
        for (unsigned int i=0; i<rules.size(); i++) {
            delete rules[i];
        }
    }
    
    
    bool parse(int argc, const char **argv)
    {
        assert(argc > 0);
        
        char *argdup = strdup(argv[0]);
        prog = basename(argdup); 
        free(argdup);
        
        if (argc < 2)
            return false;
        
        int i;
        for (i=1; i<argc; i++) {
            bool parsed = false;
            
            // detect stop parsing options
            if (!strcmp(argv[i], "--")) {
                i++;
                break;
            }
            
            for (unsigned int j=0; j<rules.size(); j++) {
                if (argv[i] == rules[j]->shortarg ||
                    argv[i] == rules[j]->longarg)
                {
                    int consume = rules[j]->parse(argc - i-1, &argv[i+1]);
                    if (consume == -1)
                        return false;
                    i += consume;
                    parsed = true;
                    break;
                }
            }
            
            // report error if no arguments are matched
            if (!parsed) {
                if (argv[i][0] == '-') {
                    fprintf(stderr, "unknown option '%s'\n", argv[i]);
                    return false;
                } else {
                    // no more options
                    break;
                }
            }
        }
        
        
        // remaining arguments are "rest" arguments
        for (; i<argc; i++)
            rest.push_back(argv[i]);
        
        return true;
    }
    
    void printHelp(FILE *stream=stderr)
    {
        fprintf(stream, "Usage: %s [OPTION]\n\n", prog.c_str());
    
        for (unsigned int i=0; i<rules.size(); i++) {
            if (rules[i]->kind == OPTION_ARG)
                fprintf(stream, "  %s,%s  %s\n    %s\n\n", 
                        rules[i]->shortarg.c_str(), rules[i]->longarg.c_str(), 
                        rules[i]->argstr.c_str(), rules[i]->help.c_str());
            else if (rules[i]->kind == OPTION_COMMENT)
                fprintf(stream, "%s\n", 
                        ((ConfigParamComment*) rules[i])->msg.c_str());
            else
                assert(0);
        }
    }
    
    void add(ConfigParamBase *rule)
    {
        rules.push_back(rule);
    }
    
    
    string prog;
    vector<ConfigParamBase*> rules;
    vector<string> rest;
};


#endif // SPIDIR_CONFIG_PARAM_H
