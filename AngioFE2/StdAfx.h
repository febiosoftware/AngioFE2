// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"

//keep these in alphabetical order
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <list>
#include <random>
#include <regex>
#include <set>
#include <unordered_map>
#include <vector>

#ifdef WIN32
#define WIN32_LEAN_AND_MEAN             // Exclude rarely-used stuff from Windows headers
#define NOMINMAX //use the c++ std versions of min and max
// Windows Header Files:
#include <windows.h>
#endif


