# sexp_lex.py. This file automatically created by PLY (version 3.2). Don't edit!
_tabversion   = '3.2'
_lextokens    = {'SYMBOL': 1, 'BOOL': 1, 'STRING': 1, 'NUMBER': 1}
_lexreflags   = 0
_lexliterals  = '()[]{}'
_lexstateinfo = {'INITIAL': 'inclusive'}
_lexstatere   = {'INITIAL': [('(?P<t_STRING>"(\\\\\\\\|\\\\"|[^"\\\\]+)*")|(?P<t_SYMBOL>[^()\\[\\]{}" \\#\\t\\r\\n0-9;][^()\\[\\]{}" \\t\\r\\n]*)|(?P<t_NUMBER>[+-]?(\\d+\\.?|\\.\\d)(\\d+([eE][+-]?\\d+)?)?)|(?P<t_BOOL>\\#t|\\#f)|(?P<t_ignore_COMMENT>;[^\\n]*\\n)', [None, ('t_STRING', 'STRING'), None, ('t_SYMBOL', 'SYMBOL'), ('t_NUMBER', 'NUMBER'), None, None, None, ('t_BOOL', 'BOOL'), (None, None)])]}
_lexstateignore = {'INITIAL': ' \t\r\n'}
_lexstateerrorf = {'INITIAL': 't_error'}
