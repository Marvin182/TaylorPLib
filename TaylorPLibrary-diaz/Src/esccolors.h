#ifndef ESCCOLORS_H
#define ESCCOLORS_H

#define ESCSEQUENCES

#ifdef ESCSEQUENCES
const char *const clearscreen = "\033[2J";
const char *const movetopleft = "\033[1;1H";
const char *const normal      = "\033[0m";
const char *const black       = "\033[0;30m";
const char *const red         = "\033[0;31m";
const char *const green       = "\033[0;32m";
const char *const yellow      = "\033[0;33m";
const char *const blue        = "\033[0;34m";
const char *const magenta     = "\033[0;35m";
const char *const cyan        = "\033[0;36m";
const char *const white       = "\033[0;37m";
#else
const char *const clearscreen = "";
const char *const movetopleft = "";
const char *const normal      = "";
const char *const black       = "";
const char *const red         = "";
const char *const green       = "";
const char *const yellow      = "";
const char *const blue        = "";
const char *const magenta     = "";
const char *const cyan        = "";
const char *const white       = "";
#endif

#endif
