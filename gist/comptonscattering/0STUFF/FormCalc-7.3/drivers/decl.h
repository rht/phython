* decl.h
* these declarations are included "everywhere"
* this file is part of FormCalc
* last modified 9 Aug 11 th


#ifndef DECL_H
#define DECL_H

* declarations for the whole file (e.g. preprocessor defs)

#include "types.h"
#include "user.h"

#else

* declarations for every subroutine

#include "const.h"
#include "user.h"
#include "util.h"
#include "looptools.h"
#include "renconst.h"

#endif

