/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * metis.h
 *
 * This file includes all necessary header files
 *
 * Started 8/27/94
 * George
 *
 * $Id: metis.h,v 1.2 2007-07-11 01:55:06 mbsteph Exp $
 */

#ifndef METIS_H
#define METIS_H

#include <stdio.h>
#ifdef __STDC__
#include <stdlib.h>
#else
#include <malloc.h>
#endif
#include <strings.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <time.h>

#ifdef DMALLOC
#include <dmalloc.h>
#endif

#include "defs.h"
#include "struct.h"
#include "macros.h"
#include "rename.h"
#include "proto.h"

#endif /* METIS_H */
