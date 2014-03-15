// -*- mode: C; c-indent-level: 2; c-basic-offset: 2; tab-width: 8 -*-
//
// Copyright (C) 2009-2014 Roberto Bertolusso and Marek Kimmel
//
// This file is part of sbioPN.
//
// sbioPN is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// sbioPN is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with sbioPN. If not, see <http://www.gnu.org/licenses/>.

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#include "quicksort.h"
#include "helper.h"

//#define rb_double long double
#define rb_double double

//#define RB_PRINT_INCR_INFO
//#define RB_SAVE_INCR_INFO

#define a21 1./5.
#define c2 1./5.

#define a31 3./40.
#define a32 9./40.
#define c3 3./10.

#define a41 44./45.
#define a42 -56./15.
#define a43 32./9.
#define c4 4./5.

#define a51 19372./6561.
#define a52 -25360./2187.
#define a53 64448./6561.
#define a54 -212./729.
#define c5 8./9.

#define a61 9017./3168.
#define a62 -355./33.
#define a63 46732./5247.
#define a64 49./176.
#define a65 -5103./18656.
#define c6 1.

#define a71 35./384.
#define a73 500./1113.
#define a74 125./192.
#define a75 -2187./6784.
#define a76 11./84.
#define c7 1.

#define b41 5179./57600.
#define b43 7571./16695.
#define b44 393./640.
#define b45 -92097./339200.
#define b46 187./2100.
#define b47 1./40.

#define b51 35./384.
#define b53 500./1113.
#define b54 125./192.
#define b55 -2187./6784.
#define b56 11./84.

#define err_pow 1./6.

#define INCR_TO_SAVE 200000

#define RB_MEMORY
#define RB_TIME
#define RB_SUBTIME

#ifdef RB_TIME
#include <time.h>
#endif

struct transition {
  // DO NOT CHANGE ORDER OF ELEMENTS!! (Optimized for total size of structure)
  char cPreNZ_len;
  char cSNZ_len;
  char cType_of_h;
  char cInGroup;
  int iPositionInGroup;
  // NOTE: 4 chars + 1 int use 8 bytes so no space is wasted in 64 bit machines

  int *piPreNZ_VoxelPlace;
  char *pcPreNZ_value;

  int *piSNZ_VoxelPlace;
  char *pcSNZ_value;

  union {
    double   dTrCnt;
    DL_FUNC *pDL_FUNC;
    SEXP     SEXP;
  } h;

  int *piHazardsToMod;

  int iHazardsToMod_count;
  int iLocalVoxelPlace_offset;

  char cSlow;
  char cSlowUpdatedByFastVoxelTransition;
  // Two bytes wasted.
};
typedef struct transition transition_el;


SEXP getListElement(SEXP list, const char *str);
/*
{
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i, len = length(list);
  
  for (i = 0; i < len; i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  return elmt;
}
*/

typedef enum {
  HZ_DOUBLE,
  HZ_CFUNCTION,
  HZ_RFUNCTION
} HZ_type;

#define MIN_INCR 1e-20

/*
#include <pthread.h>
//#include <stdio.h>
#define NUM_THREADS     5

void *PrintHello(void *threadid)
{
   long tid;
   tid = (long)threadid;
   printf("Hello World! It's me, thread #%ld!\n", tid);
   pthread_exit(NULL);
}
*/



SEXP sHaseltineRawlings(SEXP model, SEXP T, SEXP delta, SEXP runs,
			SEXP place, SEXP transition, SEXP ect, SEXP burnRnd, SEXP rho)
{
  SEXP sexpTmp;
  int *piTmp;
  char cTmp;

  /*
   pthread_t threads[NUM_THREADS];
   int rc;
   long t;
   for(t=0; t<NUM_THREADS; t++){
      printf("In main: creating thread %ld\n", t);
      rc = pthread_create(&threads[t], NULL, PrintHello, (void *)t);
      if (rc){
         printf("ERROR; return code from pthread_create() is %d\n", rc);
         exit(-1);
      }
   }
   //   pthread_exit(NULL);
   */
  int iVoxelPlaceIdx;
  int iVoxelTransitionIdx;

  int iBurnRnd = *INTEGER(burnRnd);

  if (iBurnRnd) {
    Rprintf("\nBurning %d random numbers...", iBurnRnd);
    GetRNGstate();
    int k;
    for (k = 0; k < iBurnRnd; k++) {
      unif_rand();
    }
    PutRNGstate();
    Rprintf(" done!\n", iBurnRnd);
  }

#ifdef RB_TIME
  clock_t c0, c1;
  c0 = clock();
  clock_t clkStep = 0;
  int iElSteps = 0, iTotalStepsOver50, iTotalStepsOver10;
#endif
#ifdef RB_MEMORY
  double dUsedMemAcum = 0, dThisMem;
#endif

  SEXP VoxelFamily = getListElement(model, "VoxelFamily");

  if (VoxelFamily == R_NilValue)
    error("VoxelFamily not provided!\n");

  // Process groups of Voxels
  int iVoxelFamily_len = length(VoxelFamily);
  int iVoxelFamily;
  
  int iVoxelSlowTransitions = 0,  iVoxelFastTransitions = 0;
  int iVoxelTransitions = 0;
  int iVoxelPlaces = 0;
  // First: Find out how many real Places and Transitions there are.
  for (iVoxelFamily = 0; iVoxelFamily < iVoxelFamily_len; iVoxelFamily++) {
    SEXP this_VoxelFamily = VECTOR_ELT(VoxelFamily, iVoxelFamily);

    if ((sexpTmp = getListElement(this_VoxelFamily, "pre")) == R_NilValue)
      error("pre not provided!\n");
    piTmp = INTEGER(getAttrib(sexpTmp, R_DimSymbol));
    int iOrigTransitions = piTmp[0], iOrigPlaces = piTmp[1];

    if ((sexpTmp = getListElement(this_VoxelFamily, "Voxels")) == R_NilValue)
      error("Voxels not provided!\n");
    PROTECT(sexpTmp = coerceVector(sexpTmp, INTSXP));
    int iVoxels = INTEGER(sexpTmp)[0];
    UNPROTECT_PTR(sexpTmp);

    iVoxelPlaces += iOrigPlaces * iVoxels;

    if ((sexpTmp = getListElement(this_VoxelFamily, "TransitionInVoxelSubset")) == R_NilValue)
      error("TransitionInVoxelSubset not provided!\n");
    SEXP TransitionInVoxelSubset;
    PROTECT(TransitionInVoxelSubset = coerceVector(sexpTmp, INTSXP));
    int *piTransitionInVoxelSubset = INTEGER(TransitionInVoxelSubset);

    if ((sexpTmp = getListElement(this_VoxelFamily, "slow")) == R_NilValue)
      error("slow not provided!\n");
    SEXP slow;
    PROTECT(slow = coerceVector(sexpTmp, INTSXP));
    int *piSlow = INTEGER(slow);

    SEXP VoxelSubset;
    if ((VoxelSubset = getListElement(this_VoxelFamily, "VoxelSubset")) == R_NilValue)
      error("VoxelSubset not provided!\n");

    int iOrigTransition;
    for (iOrigTransition = 0; iOrigTransition < iOrigTransitions; iOrigTransition++) {
      int iThisVoxelTransitions = length(VECTOR_ELT(VoxelSubset, piTransitionInVoxelSubset[iOrigTransition] - 1));
      iVoxelTransitions += iThisVoxelTransitions;
      if (piSlow[iOrigTransition]) {
	iVoxelSlowTransitions += iThisVoxelTransitions;
      } else {
	iVoxelFastTransitions += iThisVoxelTransitions;
      }
    }
    UNPROTECT_PTR(slow);
    UNPROTECT_PTR(TransitionInVoxelSubset);

    SEXP JumpingPlace;
    if ((JumpingPlace = getListElement(this_VoxelFamily, "JumpingPlace")) != R_NilValue) {

      SEXP JumpingPattern;
      if ((JumpingPattern = getListElement(this_VoxelFamily, "JumpingPattern")) == R_NilValue)
	error("JumpingPattern not provided!\n");
      
      int iJumpingPlace;
      for (iJumpingPlace = 0; iJumpingPlace < length(JumpingPlace); iJumpingPlace++) {    
	sexpTmp = VECTOR_ELT(VECTOR_ELT(JumpingPlace, iJumpingPlace), 2);
	PROTECT(sexpTmp = coerceVector(sexpTmp, INTSXP));
	int iJumpingPattern = INTEGER(sexpTmp)[0] - 1;
	UNPROTECT_PTR(sexpTmp);

	sexpTmp = VECTOR_ELT(VECTOR_ELT(JumpingPlace, iJumpingPlace), 3);
	PROTECT(sexpTmp = coerceVector(sexpTmp, INTSXP));
	int iSlow = INTEGER(sexpTmp)[0];
	UNPROTECT_PTR(sexpTmp);
	
	sexpTmp = VECTOR_ELT(JumpingPattern, iJumpingPattern);
	piTmp = INTEGER(getAttrib(sexpTmp, R_DimSymbol));
	int iJumps = piTmp[0];
	
	iVoxelTransitions += iJumps;      
	if (iSlow) {
	  iVoxelSlowTransitions += iJumps;
	} else {
	  iVoxelFastTransitions += iJumps;
	}
      }
    }
  }
  Rprintf("VoxelTransitions: %d\tVoxelPlaces: %d\n", iVoxelTransitions, iVoxelPlaces);
  Rprintf("VoxelSlowTransitions: %d\tVoxelFastTransitions: %d\n", iVoxelSlowTransitions, iVoxelFastTransitions);

  transition_el *pVoxelTransition = (transition_el *) R_alloc(iVoxelTransitions, sizeof(transition_el));
#ifdef RB_MEMORY
  dUsedMemAcum += (dThisMem = ((double) sizeof(transition_el) * iVoxelTransitions)/1e6);
  Rprintf("pVoxelTransition: %.2f MB (%d x %d bytes)\n", dThisMem, iVoxelTransitions, sizeof(transition_el));
#endif

  transition_el *pVoxelSlowTransition = pVoxelTransition;
  transition_el *pVoxelFastTransition = pVoxelTransition + iVoxelSlowTransitions;

  int iVoxelTransition;

  int iVoxelTransitionsPreNZ_len = 0;
  int iVoxelTransitionsSNZ_len = 0;

  int iVoxelPlacesPreNZ_len = 0;
  int *piVoxelPlacePreNZ_len = (int *) R_alloc(iVoxelPlaces, sizeof(int));
#ifdef RB_MEMORY
  dUsedMemAcum += (dThisMem = ((double) sizeof(int) * iVoxelPlaces)/1e6);
  Rprintf("piVoxelPlacePreNZ_len: %.2f MB (%d x %d bytes)\n", dThisMem, iVoxelPlaces, sizeof(int));
#endif
  int **ppiVoxelPlacePreNZ_VoxelTransition = (int **) R_alloc(iVoxelPlaces, sizeof(int *));
#ifdef RB_MEMORY
  dUsedMemAcum += (dThisMem = ((double) sizeof(int *) * iVoxelPlaces)/1e6);
  Rprintf("ppiVoxelPlacePreNZ_VoxelTransition: %.2f MB (%d x %d bytes)\n", dThisMem, iVoxelPlaces, sizeof(int *));
#endif

  int iVoxelPlacesSNZ_len = 0;
  int *piVoxelPlaceSNZ_len = (int *) R_alloc(iVoxelPlaces, sizeof(int));
#ifdef RB_MEMORY
  dUsedMemAcum += (dThisMem = ((double) sizeof(int) * iVoxelPlaces)/1e6);
  Rprintf("piVoxelPlaceSNZ_len: %.2f MB (%d x %d bytes)\n", dThisMem, iVoxelPlaces, sizeof(int));
#endif
  int **ppiVoxelPlaceSNZ_VoxelTransition = (int **) R_alloc(iVoxelPlaces, sizeof(int *));
#ifdef RB_MEMORY
  dUsedMemAcum += (dThisMem = ((double) sizeof(int *) * iVoxelPlaces)/1e6);
  Rprintf("ppiVoxelPlaceSNZ_VoxelTransition: %.2f MB (%d x %d bytes)\n", dThisMem, iVoxelPlaces, sizeof(int *));
#endif

  int iVoxelPlace;
  for(iVoxelPlace = 0; iVoxelPlace < iVoxelPlaces; iVoxelPlace++) {
    piVoxelPlacePreNZ_len[iVoxelPlace] = piVoxelPlaceSNZ_len[iVoxelPlace] = 0;
  }

  int *piPreNZ_VoxelPlace = NULL, *piSNZ_VoxelPlace = NULL;

  char *pcPreNZ_value = NULL, *pcSNZ_value = NULL;

  int iVoxelCurrentTransition = 0;
  int iVoxelPlace_offset = 0;
  // Second:
  for (iVoxelFamily = 0; iVoxelFamily < iVoxelFamily_len; iVoxelFamily++) {
    //Rprintf("Voxel group: %d\n", iVoxelFamily);
    SEXP this_VoxelFamily = VECTOR_ELT(VoxelFamily, iVoxelFamily);

    sexpTmp = getListElement(this_VoxelFamily, "pre");
    SEXP pre;
    PROTECT(pre = coerceVector(sexpTmp, INTSXP));
    int *piPre = INTEGER(pre);
    piTmp = INTEGER(getAttrib(sexpTmp, R_DimSymbol));
    int iOrigTransitions = piTmp[0], iOrigPlaces = piTmp[1];
    //Rprintf("pre: rows: %d\tcols: %d\n", iOrigTransitions, iOrigPlaces);
    
    sexpTmp = getListElement(this_VoxelFamily, "Voxels");
    PROTECT(sexpTmp = coerceVector(sexpTmp, INTSXP));
    int iVoxels = INTEGER(sexpTmp)[0];
    UNPROTECT_PTR(sexpTmp);

    sexpTmp = getListElement(this_VoxelFamily, "post");
    SEXP post;
    PROTECT(post = coerceVector(sexpTmp, INTSXP));
    int *piPost = INTEGER(post);
    
    piPreNZ_VoxelPlace = (int *) Realloc(piPreNZ_VoxelPlace, iOrigPlaces, int);
    pcPreNZ_value = (char *) Realloc(pcPreNZ_value, iOrigPlaces, char);
    piSNZ_VoxelPlace = (int *) Realloc(piSNZ_VoxelPlace, iOrigPlaces, int);
    pcSNZ_value = (char *) Realloc(pcSNZ_value, iOrigPlaces, char);

    sexpTmp = getListElement(this_VoxelFamily, "TransitionInVoxelSubset");
    SEXP TransitionInVoxelSubset;
    PROTECT(TransitionInVoxelSubset = coerceVector(sexpTmp, INTSXP));
    int *piTransitionInVoxelSubset = INTEGER(TransitionInVoxelSubset);

    SEXP VoxelSubset;
    VoxelSubset = getListElement(this_VoxelFamily, "VoxelSubset");

    SEXP h;
    h = getListElement(this_VoxelFamily, "h");

    int iOrigTransition, iOrigPlace;
    for (iOrigTransition = 0; iOrigTransition < iOrigTransitions; iOrigTransition++) {      
      char cPreNZ_len = 0;
      char cSNZ_len = 0;
      for (iOrigPlace = 0; iOrigPlace < iOrigPlaces; iOrigPlace++) {
	if ((cTmp = piPre[iOrigTransition + iOrigTransitions * iOrigPlace])) {
	  piPreNZ_VoxelPlace[(int) cPreNZ_len] = iOrigPlace;
	  pcPreNZ_value[(int) cPreNZ_len++] = cTmp;
	}
	if ((cTmp = piPost[iOrigTransition + iOrigTransitions * iOrigPlace] - 
	     piPre[iOrigTransition + iOrigTransitions * iOrigPlace])) {
	  piSNZ_VoxelPlace[(int) cSNZ_len] = iOrigPlace;
	  pcSNZ_value[(int) cSNZ_len++] = cTmp;
	}
      }

      // For the current transition, find the subset of Voxels in which it applies.
      sexpTmp = VECTOR_ELT(VoxelSubset, piTransitionInVoxelSubset[iOrigTransition] - 1);
      int iVoxelSubset_len = length(sexpTmp);
      SEXP VoxelSubset;
      PROTECT(VoxelSubset = coerceVector(sexpTmp, INTSXP));
      int *piVoxelSubset = INTEGER(VoxelSubset); // Vector with subset of Voxels

      int iVoxelSubset;
      // Create a corresponding transition on each of the Voxels of this subset.
      for (iVoxelSubset = 0; iVoxelSubset < iVoxelSubset_len; iVoxelSubset++) {
	int iVoxel = piVoxelSubset[iVoxelSubset] - 1;
	int iLocalVoxelPlace_offset = iVoxelPlace_offset + iOrigPlaces * iVoxel;
	//Rprintf("Voxel: %d\tTransition: %d\tVoxelTransition: %d\n", iVoxel, iOrigTransition, iVoxelCurrentTransition);

	iVoxelTransitionsPreNZ_len += cPreNZ_len;
	for(iVoxelPlaceIdx = 0; iVoxelPlaceIdx < cPreNZ_len; iVoxelPlaceIdx++) {
	  iVoxelPlace = iLocalVoxelPlace_offset + piPreNZ_VoxelPlace[iVoxelPlaceIdx];
	  piVoxelPlacePreNZ_len[iVoxelPlace]++;
	  iVoxelPlacesPreNZ_len++;
	}

	iVoxelTransitionsSNZ_len += cSNZ_len;
	for(iVoxelPlaceIdx = 0; iVoxelPlaceIdx < cSNZ_len; iVoxelPlaceIdx++) {
	  iVoxelPlace = iLocalVoxelPlace_offset + piSNZ_VoxelPlace[iVoxelPlaceIdx];
	  piVoxelPlaceSNZ_len[iVoxelPlace]++;
	  iVoxelPlacesSNZ_len++;
	}
	iVoxelCurrentTransition++;
      }
      UNPROTECT_PTR(VoxelSubset);
    }
    UNPROTECT_PTR(TransitionInVoxelSubset);
    UNPROTECT_PTR(post);
    UNPROTECT_PTR(pre);

    // Manage jumps.
    SEXP JumpingPlace;
    if ((JumpingPlace = getListElement(this_VoxelFamily, "JumpingPlace")) != R_NilValue) {
      SEXP JumpingPattern;
      JumpingPattern = getListElement(this_VoxelFamily, "JumpingPattern");
      
      int iJumpingPlace;
      for (iJumpingPlace = 0; iJumpingPlace < length(JumpingPlace); iJumpingPlace++) {
	SEXP this_JumpingPlace = VECTOR_ELT(JumpingPlace, iJumpingPlace);

	sexpTmp = VECTOR_ELT(this_JumpingPlace, 0);
	PROTECT(sexpTmp = coerceVector(sexpTmp, INTSXP));
	int iOrigPlace = INTEGER(sexpTmp)[0] - 1;
	UNPROTECT_PTR(sexpTmp);

	sexpTmp = VECTOR_ELT(this_JumpingPlace, 2);
	PROTECT(sexpTmp = coerceVector(sexpTmp, INTSXP));
	int iJumpingPattern = INTEGER(sexpTmp)[0] - 1;
	UNPROTECT_PTR(sexpTmp);
	
	sexpTmp = VECTOR_ELT(JumpingPattern, iJumpingPattern);
	piTmp = INTEGER(getAttrib(sexpTmp, R_DimSymbol));
	int iJumps = piTmp[0];
	SEXP jump;
	PROTECT(jump = coerceVector(sexpTmp, INTSXP));
	int *piJump = INTEGER(jump);
	int iJump;
	for (iJump = 0; iJump < iJumps; iJump++) {
	  int iFromVoxelPlace = iVoxelPlace_offset + iOrigPlaces * (piJump[iJump] - 1) + iOrigPlace;
	  int iToVoxelPlace = iVoxelPlace_offset + iOrigPlaces * (piJump[iJump + iJumps] - 1) + iOrigPlace;
	  char cToVoxelPlacePost_value = piJump[iJump + iJumps*2];

	  iVoxelTransitionsPreNZ_len++;

	  piVoxelPlacePreNZ_len[iFromVoxelPlace]++;
	  iVoxelPlacesPreNZ_len++;

	  iVoxelTransitionsSNZ_len++;

	  piVoxelPlaceSNZ_len[iFromVoxelPlace]++;
	  iVoxelPlacesSNZ_len++;
	  if (cToVoxelPlacePost_value) {
	    iVoxelTransitionsSNZ_len++;

	    piVoxelPlaceSNZ_len[iToVoxelPlace]++;
	    iVoxelPlacesSNZ_len++;
	  }

	  iVoxelCurrentTransition++;

	  //Rprintf("%d\t%d\n", iFromVoxelPlace, iToVoxelPlace);
	}
	UNPROTECT_PTR(jump);
      }
    }
    iVoxelPlace_offset += iOrigPlaces * iVoxels;
  }

  int *piThisPreNZ = (int *) R_alloc(iVoxelPlacesPreNZ_len, sizeof(int));
#ifdef RB_MEMORY
  dUsedMemAcum += (dThisMem = ((double) sizeof(int) * iVoxelPlacesPreNZ_len)/1e6);
  Rprintf("ppiVoxelPlacePreNZ_VoxelTransition elements: %.2f MB (%d x %d bytes)\n", dThisMem, iVoxelPlacesPreNZ_len, sizeof(int));
#endif

  int *piThisSNZ = (int *) R_alloc(iVoxelPlacesSNZ_len, sizeof(int));
#ifdef RB_MEMORY
  dUsedMemAcum += (dThisMem = ((double) sizeof(int) * iVoxelPlacesSNZ_len)/1e6);
  Rprintf("ppiVoxelPlaceSNZ_VoxelTransition elements: %.2f MB (%d x %d bytes)\n", dThisMem, iVoxelPlacesSNZ_len, sizeof(int));
#endif

  for(iVoxelPlace = 0; iVoxelPlace < iVoxelPlaces; iVoxelPlace++) {
    ppiVoxelPlacePreNZ_VoxelTransition[iVoxelPlace] = piThisPreNZ;
    ppiVoxelPlaceSNZ_VoxelTransition[iVoxelPlace] = piThisSNZ;

    piThisPreNZ += piVoxelPlacePreNZ_len[iVoxelPlace];
    piThisSNZ += piVoxelPlaceSNZ_len[iVoxelPlace];

    piVoxelPlacePreNZ_len[iVoxelPlace] = 0;
    piVoxelPlaceSNZ_len[iVoxelPlace] = 0;
  }

  piThisPreNZ = (int *) R_alloc(iVoxelTransitionsPreNZ_len, sizeof(int));
#ifdef RB_MEMORY
  dUsedMemAcum += (dThisMem = ((double) sizeof(int) * iVoxelTransitionsPreNZ_len)/1e6);
  Rprintf("piPreNZ_VoxelPlace elements: %.2f MB (%d x %d bytes)\n", dThisMem, iVoxelTransitionsPreNZ_len, sizeof(int));
#endif

  char *pcThisPreNZ_value = (char *) R_alloc(iVoxelTransitionsPreNZ_len, sizeof(char));
#ifdef RB_MEMORY
  dUsedMemAcum += (dThisMem = ((double) sizeof(char) * iVoxelTransitionsPreNZ_len)/1e6);
  Rprintf("pcPreNZ_value total: %.2f MB (%d x %d bytes)\n", dThisMem, iVoxelTransitionsPreNZ_len, sizeof(char));
#endif
  
  piThisSNZ = (int *) R_alloc(iVoxelTransitionsSNZ_len, sizeof(int));
#ifdef RB_MEMORY
  dUsedMemAcum += (dThisMem = ((double) sizeof(int) * iVoxelTransitionsSNZ_len)/1e6);
  Rprintf("piSNZ_VoxelPlace elements: %.2f MB (%d x %d bytes)\n", dThisMem, iVoxelTransitionsSNZ_len, sizeof(int));
#endif

  char *pcThisSNZ_value = (char *) R_alloc(iVoxelTransitionsSNZ_len, sizeof(char));
#ifdef RB_MEMORY
  dUsedMemAcum += (dThisMem = ((double) sizeof(char) * iVoxelTransitionsSNZ_len)/1e6);
  Rprintf("pcSNZ_value total: %.2f MB (%d x %d bytes)\n", dThisMem, iVoxelTransitionsSNZ_len, sizeof(char));
#endif

  iVoxelCurrentTransition = 0;
  int iVoxelCurrentSlowTransition = 0, iVoxelCurrentFastTransition = 0;
  iVoxelPlace_offset = 0;
  // Third:
  for (iVoxelFamily = 0; iVoxelFamily < iVoxelFamily_len; iVoxelFamily++) {
    //Rprintf("Voxel group: %d\n", iVoxelFamily);
    SEXP this_VoxelFamily = VECTOR_ELT(VoxelFamily, iVoxelFamily);

    sexpTmp = getListElement(this_VoxelFamily, "pre");
    SEXP pre;
    PROTECT(pre = coerceVector(sexpTmp, INTSXP));
    int *piPre = INTEGER(pre);
    piTmp = INTEGER(getAttrib(sexpTmp, R_DimSymbol));
    int iOrigTransitions = piTmp[0], iOrigPlaces = piTmp[1];
    //Rprintf("pre: rows: %d\tcols: %d\n", iOrigTransitions, iOrigPlaces);
    
    sexpTmp = getListElement(this_VoxelFamily, "Voxels");
    PROTECT(sexpTmp = coerceVector(sexpTmp, INTSXP));
    int iVoxels = INTEGER(sexpTmp)[0];
    UNPROTECT_PTR(sexpTmp);

    sexpTmp = getListElement(this_VoxelFamily, "post");
    SEXP post;
    PROTECT(post = coerceVector(sexpTmp, INTSXP));
    int *piPost = INTEGER(post);
    
    piPreNZ_VoxelPlace = (int *) Realloc(piPreNZ_VoxelPlace, iOrigPlaces, int);
    pcPreNZ_value = (char *) Realloc(pcPreNZ_value, iOrigPlaces, char);
    piSNZ_VoxelPlace = (int *) Realloc(piSNZ_VoxelPlace, iOrigPlaces, int);
    pcSNZ_value = (char *) Realloc(pcSNZ_value, iOrigPlaces, char);

    sexpTmp = getListElement(this_VoxelFamily, "TransitionInVoxelSubset");
    SEXP TransitionInVoxelSubset;
    PROTECT(TransitionInVoxelSubset = coerceVector(sexpTmp, INTSXP));
    int *piTransitionInVoxelSubset = INTEGER(TransitionInVoxelSubset);

    sexpTmp = getListElement(this_VoxelFamily, "slow");
    SEXP slow;
    PROTECT(slow = coerceVector(sexpTmp, INTSXP));
    int *piSlow = INTEGER(slow);

    SEXP VoxelSubset;
    VoxelSubset = getListElement(this_VoxelFamily, "VoxelSubset");

    SEXP h;
    h = getListElement(this_VoxelFamily, "h");

    int iOrigTransition, iOrigPlace;
    for (iOrigTransition = 0; iOrigTransition < iOrigTransitions; iOrigTransition++) {      
      char cPreNZ_len = 0;
      char cSNZ_len = 0;
      for (iOrigPlace = 0; iOrigPlace < iOrigPlaces; iOrigPlace++) {
	if ((cTmp = piPre[iOrigTransition + iOrigTransitions * iOrigPlace])) {
	  piPreNZ_VoxelPlace[(int) cPreNZ_len] = iOrigPlace;
	  pcPreNZ_value[(int) cPreNZ_len++] = cTmp;
	}
	if ((cTmp = piPost[iOrigTransition + iOrigTransitions * iOrigPlace] - 
	     piPre[iOrigTransition + iOrigTransitions * iOrigPlace])) {
	  piSNZ_VoxelPlace[(int) cSNZ_len] = iOrigPlace;
	  pcSNZ_value[(int) cSNZ_len++] = cTmp;
	}
      }

      // For the current transition, find the subset of Voxels in which it applies.
      sexpTmp = VECTOR_ELT(VoxelSubset, piTransitionInVoxelSubset[iOrigTransition] - 1);
      int iVoxelSubset_len = length(sexpTmp);
      SEXP VoxelSubset;
      PROTECT(VoxelSubset = coerceVector(sexpTmp, INTSXP));
      int *piVoxelSubset = INTEGER(VoxelSubset); // Vector with subset of Voxels

      int iVoxelSubset;
      // Create a corresponding transition on each of the Voxels of this subset.
      for (iVoxelSubset = 0; iVoxelSubset < iVoxelSubset_len; iVoxelSubset++) {
	int iVoxel = piVoxelSubset[iVoxelSubset] - 1;
	int iLocalVoxelPlace_offset = iVoxelPlace_offset + iOrigPlaces * iVoxel;
	//Rprintf("Voxel: %d\tTransition: %d\tVoxelTransition: %d\n", iVoxel, iOrigTransition, iVoxelCurrentTransition);

	iVoxelCurrentTransition = (piSlow[iOrigTransition] ? iVoxelCurrentSlowTransition++ : iVoxelSlowTransitions + iVoxelCurrentFastTransition++);

	transition_el *pThisVoxelTransition = &pVoxelTransition[iVoxelCurrentTransition];

	if ((pThisVoxelTransition->cSlow = (char) piSlow[iOrigTransition]))
	  pThisVoxelTransition->cSlowUpdatedByFastVoxelTransition = 0;

	if (inherits(sexpTmp = VECTOR_ELT(h, iOrigTransition), "NativeSymbol")) {
	  pThisVoxelTransition->h.pDL_FUNC = (void *) R_ExternalPtrAddr(sexpTmp);
	  pThisVoxelTransition->cType_of_h = HZ_CFUNCTION;    
	} else if (isNumeric(sexpTmp)){
	  pThisVoxelTransition->h.dTrCnt = REAL(sexpTmp)[0];
	  pThisVoxelTransition->cType_of_h = HZ_DOUBLE;
	} else  if (isFunction(sexpTmp)) {
	  //SET_VECTOR_ELT(sexpFunction, iOrigTransition, lang1(sexpTmp));
	  //piHzType[iOrigTransition] = HZ_RFUNCTION;
	  ;
	} else {
	  error("Unrecognized transition function type\n");
	}

	pThisVoxelTransition->iLocalVoxelPlace_offset = iLocalVoxelPlace_offset;
	pThisVoxelTransition->cPreNZ_len = cPreNZ_len;
	pThisVoxelTransition->piPreNZ_VoxelPlace = piThisPreNZ;
	piThisPreNZ += cPreNZ_len;
	pThisVoxelTransition->pcPreNZ_value = pcThisPreNZ_value;
	pcThisPreNZ_value += cPreNZ_len;
	for(iVoxelPlaceIdx = 0; iVoxelPlaceIdx < cPreNZ_len; iVoxelPlaceIdx++) {
	  iVoxelPlace = iLocalVoxelPlace_offset + piPreNZ_VoxelPlace[iVoxelPlaceIdx];
	  pThisVoxelTransition->piPreNZ_VoxelPlace[iVoxelPlaceIdx] = iVoxelPlace;
	  pThisVoxelTransition->pcPreNZ_value[iVoxelPlaceIdx] = pcPreNZ_value[iVoxelPlaceIdx];

	  ppiVoxelPlacePreNZ_VoxelTransition[iVoxelPlace][piVoxelPlacePreNZ_len[iVoxelPlace]++] = iVoxelCurrentTransition;
	}

	pThisVoxelTransition->cSNZ_len = cSNZ_len;
	pThisVoxelTransition->piSNZ_VoxelPlace = piThisSNZ;
	piThisSNZ += cSNZ_len;
	pThisVoxelTransition->pcSNZ_value = pcThisSNZ_value;
	pcThisSNZ_value += cSNZ_len;
	for(iVoxelPlaceIdx = 0; iVoxelPlaceIdx < cSNZ_len; iVoxelPlaceIdx++) {
	  iVoxelPlace = iLocalVoxelPlace_offset + piSNZ_VoxelPlace[iVoxelPlaceIdx];
	  pThisVoxelTransition->piSNZ_VoxelPlace[iVoxelPlaceIdx] = iVoxelPlace;
	  pThisVoxelTransition->pcSNZ_value[iVoxelPlaceIdx] = pcSNZ_value[iVoxelPlaceIdx];

	  ppiVoxelPlaceSNZ_VoxelTransition[iVoxelPlace][piVoxelPlaceSNZ_len[iVoxelPlace]++] = iVoxelCurrentTransition;
	}
      }
      UNPROTECT_PTR(VoxelSubset);
    }
    UNPROTECT_PTR(slow);
    UNPROTECT_PTR(TransitionInVoxelSubset);
    UNPROTECT_PTR(post);
    UNPROTECT_PTR(pre);

    // Manage jumps.
    SEXP JumpingPlace;
    if ((JumpingPlace = getListElement(this_VoxelFamily, "JumpingPlace")) != R_NilValue) {

      SEXP JumpingPattern;
      JumpingPattern = getListElement(this_VoxelFamily, "JumpingPattern");
      
      int iJumpingPlace;
      for (iJumpingPlace = 0; iJumpingPlace < length(JumpingPlace); iJumpingPlace++) {
	SEXP this_JumpingPlace = VECTOR_ELT(JumpingPlace, iJumpingPlace);

	sexpTmp = VECTOR_ELT(this_JumpingPlace, 0);
	PROTECT(sexpTmp = coerceVector(sexpTmp, INTSXP));
	int iOrigPlace = INTEGER(sexpTmp)[0] - 1;
	UNPROTECT_PTR(sexpTmp);

	sexpTmp = VECTOR_ELT(this_JumpingPlace, 1);
	PROTECT(sexpTmp = coerceVector(sexpTmp, REALSXP));
	double dH = REAL(sexpTmp)[0];
	UNPROTECT_PTR(sexpTmp);

	sexpTmp = VECTOR_ELT(this_JumpingPlace, 2);
	PROTECT(sexpTmp = coerceVector(sexpTmp, INTSXP));
	int iJumpingPattern = INTEGER(sexpTmp)[0] - 1;
	UNPROTECT_PTR(sexpTmp);

	sexpTmp = VECTOR_ELT(VECTOR_ELT(JumpingPlace, iJumpingPlace), 3);
	PROTECT(sexpTmp = coerceVector(sexpTmp, INTSXP));
	int iSlow = INTEGER(sexpTmp)[0];
	UNPROTECT_PTR(sexpTmp);
	
	sexpTmp = VECTOR_ELT(JumpingPattern, iJumpingPattern);
	piTmp = INTEGER(getAttrib(sexpTmp, R_DimSymbol));
	int iJumps = piTmp[0];
	SEXP jump;
	PROTECT(jump = coerceVector(sexpTmp, INTSXP));
	int *piJump = INTEGER(jump);
	int iJump;
	for (iJump = 0; iJump < iJumps; iJump++) {
	  int iFromVoxelPlace = iVoxelPlace_offset + iOrigPlaces * (piJump[iJump] - 1) + iOrigPlace;
	  int iToVoxelPlace = iVoxelPlace_offset + iOrigPlaces * (piJump[iJump + iJumps] - 1) + iOrigPlace;
	  char cToVoxelPlacePost_value = piJump[iJump + iJumps*2];

	iVoxelCurrentTransition = (iSlow ? iVoxelCurrentSlowTransition++ : iVoxelSlowTransitions + iVoxelCurrentFastTransition++);

	  transition_el *pThisVoxelTransition = &pVoxelTransition[iVoxelCurrentTransition];
	  
	  if ((pThisVoxelTransition->cSlow = (char) iSlow))
	    pThisVoxelTransition->cSlowUpdatedByFastVoxelTransition = 0;

	  pThisVoxelTransition->cPreNZ_len = 1;
	  pThisVoxelTransition->piPreNZ_VoxelPlace = piThisPreNZ;
	  piThisPreNZ++;
	  pThisVoxelTransition->pcPreNZ_value = pcThisPreNZ_value;
	  pcThisPreNZ_value++;

	  pThisVoxelTransition->h.dTrCnt = dH;
	  pThisVoxelTransition->cType_of_h = HZ_DOUBLE;
	  pThisVoxelTransition->iLocalVoxelPlace_offset = iVoxelPlace_offset + iOrigPlaces * (piJump[iJump] - 1);

	  pThisVoxelTransition->piPreNZ_VoxelPlace[0] = iFromVoxelPlace;
	  pThisVoxelTransition->pcPreNZ_value[0] = 1;

	  ppiVoxelPlacePreNZ_VoxelTransition[iFromVoxelPlace][piVoxelPlacePreNZ_len[iFromVoxelPlace]++] = iVoxelCurrentTransition;

	  pThisVoxelTransition->cSNZ_len = 1;
	  pThisVoxelTransition->piSNZ_VoxelPlace = piThisSNZ;
	  piThisSNZ++;
	  pThisVoxelTransition->pcSNZ_value = pcThisSNZ_value;
	  pcThisSNZ_value++;
	  pThisVoxelTransition->piSNZ_VoxelPlace[0] = iFromVoxelPlace;
	  pThisVoxelTransition->pcSNZ_value[0] = -1;
	  ppiVoxelPlaceSNZ_VoxelTransition[iFromVoxelPlace][piVoxelPlaceSNZ_len[iFromVoxelPlace]++] = iVoxelCurrentTransition;
	  if (cToVoxelPlacePost_value) {
	    pThisVoxelTransition->cSNZ_len++;
	    piThisSNZ++;
	    pcThisSNZ_value++;
	    pThisVoxelTransition->piSNZ_VoxelPlace[1] = iToVoxelPlace;
	    pThisVoxelTransition->pcSNZ_value[1] = cToVoxelPlacePost_value;
	    ppiVoxelPlaceSNZ_VoxelTransition[iToVoxelPlace][piVoxelPlaceSNZ_len[iToVoxelPlace]++] = iVoxelCurrentTransition;	  
	  }
	  //Rprintf("%d\t%d\n", iFromVoxelPlace, iToVoxelPlace);
	}
	UNPROTECT_PTR(jump);
      }
    }
    iVoxelPlace_offset += iOrigPlaces * iVoxels;
  }

  Free(piPreNZ_VoxelPlace);
  Free(pcPreNZ_value);
  Free(piSNZ_VoxelPlace);
  Free(pcSNZ_value);

  //  Rprintf("VoxelTransitions: %d\tVoxelCurrentTransition: %d\n", iVoxelTransitions, iVoxelCurrentTransition);

  // For the initial calculation of all hazards...
  int *piHazardsToMod_start = (int *) R_alloc(iVoxelTransitions, sizeof(int));
#ifdef RB_MEMORY
  dUsedMemAcum += (dThisMem = ((double) sizeof(int) * iVoxelTransitions)/1e6);
  Rprintf("piHazardsToMod_start: %.2f MB (%d x %d bytes)\n", dThisMem, iVoxelTransitions, sizeof(int));
#endif
  
  int iHazardsToMod_count = 0;
  // Identify slow transitions (only counts) that need the hazards recalculated
  // when slow transition fire or depending on fast reactions.
  for(iVoxelTransition = 0; iVoxelTransition < iVoxelTransitions; iVoxelTransition++) {
    transition_el *pThisVoxelTransition = &pVoxelTransition[iVoxelTransition];
    int iHazardToCompTot = 0;
    for(iVoxelPlaceIdx = 0; iVoxelPlaceIdx < pThisVoxelTransition->cSNZ_len; iVoxelPlaceIdx++) {
      int iVoxelPlace = pThisVoxelTransition->piSNZ_VoxelPlace[iVoxelPlaceIdx];
      
      int iVoxelTransitionIdx2, iVoxelTransition2;
      for(iVoxelTransitionIdx2 = 0; iVoxelTransitionIdx2 < piVoxelPlacePreNZ_len[iVoxelPlace]; iVoxelTransitionIdx2++) {
	iVoxelTransition2 = ppiVoxelPlacePreNZ_VoxelTransition[iVoxelPlace][iVoxelTransitionIdx2];
	if (iVoxelTransition2 < iVoxelSlowTransitions) {
	  if(iVoxelTransition < iVoxelSlowTransitions) {
	    // SlowVoxelTransition updated by a SlowVoxelTransition
	    int iAddThis = TRUE;
	    int iVoxelTransitionIdx3;
	    for (iVoxelTransitionIdx3 = 0; iVoxelTransitionIdx3 < iHazardToCompTot; iVoxelTransitionIdx3++) {
	      if(piHazardsToMod_start[iVoxelTransitionIdx3] == iVoxelTransition2) {
		iAddThis = FALSE;
		break;
	      }
	    }	    
	    if (iAddThis)
	      piHazardsToMod_start[iHazardToCompTot++] = iVoxelTransition2;
	  } else {
	    // SlowVoxelTransition updated by a FastVoxelTransition
	    // Mark it generally as all this kind of SlowVoxelTransitions
	    // are updated together once a deterministic step is performed.
	    pVoxelTransition[iVoxelTransition2].cSlowUpdatedByFastVoxelTransition = 1;
	  }
	}
      }
    }
    iHazardsToMod_count += iHazardToCompTot;
  }  
  int *piThisHazardsToMod = (int *) R_alloc(iHazardsToMod_count, sizeof(int));
#ifdef RB_MEMORY
  dUsedMemAcum += (dThisMem = ((double) sizeof(int) * iHazardsToMod_count)/1e6);
  Rprintf("piHazardsToMod elements: %.2f MB (%d x %d bytes)\n", dThisMem, iHazardsToMod_count, sizeof(int));
#endif
  int iSlowHazardsToModFromFast_len = 0;
  int *piSlowHazardsToModFromFast = (int *) R_alloc(iVoxelSlowTransitions, sizeof(int));
#ifdef RB_MEMORY
  dUsedMemAcum += (dThisMem = ((double) sizeof(int) * iVoxelSlowTransitions)/1e6);
  Rprintf("piSlowHazardsToModFromFast: %.2f MB (%d x %d bytes)\n", dThisMem, iVoxelSlowTransitions, sizeof(int));
#endif
  //int iHazardsToMod_count2 = 0;
  // Identify slow transitions that need the hazards recalculated
  for(iVoxelTransition = 0; iVoxelTransition < iVoxelSlowTransitions; iVoxelTransition++) {
    transition_el *pThisVoxelTransition = &pVoxelSlowTransition[iVoxelTransition];
    if(pThisVoxelTransition->cSlowUpdatedByFastVoxelTransition) {
      piSlowHazardsToModFromFast[iSlowHazardsToModFromFast_len++] = iVoxelTransition;
    }
    int iHazardToCompTot = 0;
    for(iVoxelPlaceIdx = 0; iVoxelPlaceIdx < pThisVoxelTransition->cSNZ_len; iVoxelPlaceIdx++) {
      int iVoxelPlace = pThisVoxelTransition->piSNZ_VoxelPlace[iVoxelPlaceIdx];

      int iVoxelTransitionIdx2, iVoxelTransition2;
      for(iVoxelTransitionIdx2 = 0; iVoxelTransitionIdx2 < piVoxelPlacePreNZ_len[iVoxelPlace]; iVoxelTransitionIdx2++) {
	iVoxelTransition2 = ppiVoxelPlacePreNZ_VoxelTransition[iVoxelPlace][iVoxelTransitionIdx2];
	if (iVoxelTransition2 < iVoxelSlowTransitions) {
	  // SlowVoxelTransition updated by a SlowVoxelTransition
	  int iAddThis = TRUE;
	  int iVoxelTransitionIdx3;
	  for (iVoxelTransitionIdx3 = 0; iVoxelTransitionIdx3 < iHazardToCompTot; iVoxelTransitionIdx3++) {
	    if(piHazardsToMod_start[iVoxelTransitionIdx3] == iVoxelTransition2) {
	      iAddThis = FALSE;
	      break;
	    }
	  }	    
	  if (iAddThis)
	    piHazardsToMod_start[iHazardToCompTot++] = iVoxelTransition2;
	}
      }
    }
    pThisVoxelTransition->iHazardsToMod_count = iHazardToCompTot;
    pThisVoxelTransition->piHazardsToMod = piThisHazardsToMod;
    piThisHazardsToMod += iHazardToCompTot;

    //Rprintf("%d\t%d\n", iVoxelTransition < iVoxelSlowTransitions, iHazardToCompTot);
    for (iVoxelTransitionIdx = 0; iVoxelTransitionIdx < iHazardToCompTot; iVoxelTransitionIdx++) {
      pThisVoxelTransition->piHazardsToMod[iVoxelTransitionIdx] = piHazardsToMod_start[iVoxelTransitionIdx];
      //      iHazardsToMod_count2++;
    }
  }

  //Rprintf("%d\t%d\n", iHazardsToMod_count, iHazardsToMod_count2);
  //  error("");
  int iVoxelFastPlaceIdx_len = 0;
  int *piVoxelFastPlaceIdx = (int *) R_alloc(iVoxelPlaces, sizeof(int));
#ifdef RB_MEMORY
  dUsedMemAcum += (dThisMem = ((double) sizeof(int) * iVoxelPlaces)/1e6);
  Rprintf("piVoxelFastPlaceIdx: %.2f MB (%d x %d bytes)\n", dThisMem, iVoxelPlaces, sizeof(int));
#endif

  for(iVoxelPlace = 0; iVoxelPlace < iVoxelPlaces; iVoxelPlace++) {
    for (iVoxelTransitionIdx = 0; iVoxelTransitionIdx < piVoxelPlaceSNZ_len[iVoxelPlace]; iVoxelTransitionIdx++) {
      if (ppiVoxelPlaceSNZ_VoxelTransition[iVoxelPlace][iVoxelTransitionIdx] >= iVoxelSlowTransitions) {
	piVoxelFastPlaceIdx[iVoxelFastPlaceIdx_len++] = iVoxelPlace;
	break;
      }
    }
  }

  // For the initial calculation of all slow hazards...
  for(iVoxelTransition = 0; iVoxelTransition < iVoxelSlowTransitions; iVoxelTransition++) {
    piHazardsToMod_start[iVoxelTransition] = iVoxelTransition;
  }

  rb_double *pdK1 = (rb_double *) R_alloc(7*iVoxelPlaces, sizeof(rb_double));
#ifdef RB_MEMORY
  dUsedMemAcum += (dThisMem = ((double) sizeof(rb_double) * 7*iVoxelPlaces)/1e6);
  Rprintf("pdK1 to pdK7: %.2f MB (%d x %d bytes)\n", dThisMem, 7*iVoxelPlaces, sizeof(rb_double));
#endif
  rb_double *pdK2 = pdK1 + iVoxelPlaces;
  rb_double *pdK3 = pdK2 + iVoxelPlaces;
  rb_double *pdK4 = pdK3 + iVoxelPlaces;
  rb_double *pdK5 = pdK4 + iVoxelPlaces;
  rb_double *pdK6 = pdK5 + iVoxelPlaces;
  rb_double *pdK7 = pdK6 + iVoxelPlaces;
  rb_double *pdYdot = 0;
  
  double dEct = *REAL(ect);
  
  SEXP sexpTmpCrntMarking;
  PROTECT(sexpTmpCrntMarking = allocVector(REALSXP, iVoxelPlaces));
  // double *pdTmpCrntMarking = REAL(sexpTmpCrntMarking);
  rb_double *pdTmpCrntMarking = (rb_double *) R_alloc(iVoxelPlaces, sizeof(rb_double));
  
  SEXP sexpCrntMarking;
  PROTECT(sexpCrntMarking = allocVector(REALSXP, iVoxelPlaces));
    double *pdCrntMarking = REAL(sexpCrntMarking);
    //rb_double *pdCrntMarking = (rb_double *) R_alloc(iVoxelPlaces, sizeof(rb_double));

  rb_double *pdBakCrntMarking = (rb_double *) R_alloc(iVoxelPlaces, sizeof(rb_double));

  rb_double *pdCrntDiffMarking = (rb_double *) R_alloc(iVoxelPlaces, sizeof(rb_double));

  double dDelta = *REAL(delta);
  int iTotalSteps, iSectionSteps;
  double dT = 0;
  void *pCManage_time = 0;
  SEXP sexpRManage_time = 0;
  if (inherits(T, "NativeSymbol")) {
    pCManage_time = (void *) R_ExternalPtrAddr(T);
    dT = ((double(*)(double, rb_double *)) pCManage_time)(-1, pdCrntMarking);
  } else if (isNumeric(T)){
    dT = *REAL(T);
  } else if (isFunction(T)) {
    PROTECT(sexpRManage_time = lang1(T));
    
    defineVar(install("y"), sexpCrntMarking, rho);
    PROTECT(sexpTmp = allocVector(REALSXP, 1));
    *REAL(sexpTmp) = -1;
    defineVar(install("StartTime"), sexpTmp, rho);
    UNPROTECT_PTR(sexpTmp);
    dT = *REAL(VECTOR_ELT(eval(sexpRManage_time, rho),0));
  } else {
    error("Unrecognized time function type\n");
  }

  iTotalSteps = iSectionSteps = (int)(dT / dDelta) + 1;
#ifdef RB_TIME
  iTotalStepsOver50 = iTotalSteps/50;
  iTotalStepsOver10 = iTotalSteps/10;
#endif
  // Hazard vector
  double *pdHazard = (double *) R_alloc(iVoxelSlowTransitions, sizeof(double));
  double *pdBakHazard = (double *) R_alloc(iVoxelSlowTransitions, sizeof(double));

  int iRun, iRuns = *INTEGER(runs);
  
  SEXP sexpRun;
  PROTECT(sexpRun = allocVector(VECSXP, iRuns));
  
  int iTotalUsedRandomNumbers = 0;
  
  // DiscTime Vector
  SEXP sexpD_time;
  PROTECT(sexpD_time = allocVector(REALSXP, iTotalSteps));
  double *pdDiscTime = REAL(sexpD_time);
  double dTmp = 0;
  int k;
  for (k = 0; k < iTotalSteps; k++) {
    pdDiscTime[k] = dTmp;
    dTmp += dDelta;
  }
  
  SEXP sexpMarkingRowNames;
  PROTECT(sexpMarkingRowNames = allocVector(INTSXP, iTotalSteps));
  piTmp = INTEGER(sexpMarkingRowNames);
  for (k = 0; k < iTotalSteps; k++)
    piTmp[k] = k+1;
  
  double **ppdMarking = (double **) R_alloc(iVoxelPlaces, sizeof(double *));

#ifdef RB_SAVE_INCR_INFO
  double *pdIncr = (double *) R_alloc(INCR_TO_SAVE, sizeof(double));
  double *pdIncrTime = (double *) R_alloc(INCR_TO_SAVE, sizeof(double));
  double *pdAcumHazard = (double *) R_alloc(INCR_TO_SAVE, sizeof(double));
  double *pdIntHazard = (double *) R_alloc(INCR_TO_SAVE, sizeof(double));
  double *pdIntHazardTime = (double *) R_alloc(INCR_TO_SAVE, sizeof(double));
#endif

  int *piOrderedTransition = (int *) R_alloc(iVoxelSlowTransitions, sizeof(int));

#ifdef RB_TIME
  c1 = clock();
  Rprintf ("Elapsed CPU time: %f\n", (float) (c1 - c0)/CLOCKS_PER_SEC);
  c0 = c1;
#endif

  GetRNGstate();
  for (iRun = 0; iRun < iRuns; iRun++) {
    
#ifdef RB_SAVE_INCR_INFO
    int iTotAccpIncr = 0, iTotRejIncr = 0, iTotGoBackIncr = 0, iTotIntHazardTime = 0, iTotIncrTime = 0;
    double dSumAccpIncr = 0;
    double dSumSqAccpIncr = 0;
#endif
    
    int iUsedRandomNumbers = 0;
    Rprintf("%d ", iRun+1);
    
    // Totals for kind of transition vector
    SEXP sexpTotXTransition;
    PROTECT(sexpTotXTransition = allocVector(INTSXP, iVoxelTransitions));
    int *piTotTransitions = INTEGER(sexpTotXTransition);
    
    for(iVoxelTransition = 0; iVoxelTransition < iVoxelTransitions; iVoxelTransition++) {
      piTotTransitions[iVoxelTransition] = 0;
    }
    for(iVoxelTransition = 0; iVoxelTransition < iVoxelSlowTransitions; iVoxelTransition++) {
      piOrderedTransition[iVoxelTransition] = iVoxelTransition;
    }
    int iTillResort = 1000, iTotResort = 0;

    // Totals for Infinite Norm error per place
    SEXP sexpInfNormErrorXVoxelPlace;
    PROTECT(sexpInfNormErrorXVoxelPlace = allocVector(INTSXP, iVoxelPlaces));
    int *piInfNormErrorXVoxelPlace = INTEGER(sexpInfNormErrorXVoxelPlace);
    
    for(iVoxelPlace = 0; iVoxelPlace < iVoxelPlaces; iVoxelPlace++) {
      piInfNormErrorXVoxelPlace[iVoxelPlace] = 0;
    }
    
    SEXP sexpMarking;
    PROTECT(sexpMarking = allocVector(VECSXP, iVoxelPlaces));
    // NOTE: Need to address next line (expansion of place names!)
    //setAttrib(sexpMarking, R_NamesSymbol, place);
    //setAttrib(sexpMarking, R_RowNamesSymbol, sexpMarkingRowNames);
    //setAttrib(sexpMarking, R_ClassSymbol, ScalarString(mkChar("data.frame")));

    // Setup initial state
    int iVoxelCurrentPlace = 0;
    for (iVoxelFamily = 0; iVoxelFamily < iVoxelFamily_len; iVoxelFamily++) {
      SEXP this_VoxelFamily = VECTOR_ELT(VoxelFamily, iVoxelFamily);
      
      sexpTmp = getListElement(this_VoxelFamily, "pre");
      piTmp = INTEGER(getAttrib(sexpTmp, R_DimSymbol));
      int iOrigPlaces = piTmp[1];
      
      sexpTmp = getListElement(this_VoxelFamily, "Voxels");
      PROTECT(sexpTmp = coerceVector(sexpTmp, INTSXP));
      int iVoxels = INTEGER(sexpTmp)[0];
      UNPROTECT_PTR(sexpTmp);
      
      if ((sexpTmp = getListElement(this_VoxelFamily, "M")) == R_NilValue)
	error("M not provided!\n");
      SEXP M;
      PROTECT(M = coerceVector(sexpTmp, REALSXP));
      double *pdM = REAL(M);
      
      //Rprintf("\n");
      int iVoxel, iOrigPlace;
      for (iVoxel = 0; iVoxel < iVoxels; iVoxel++) {
	for (iOrigPlace = 0; iOrigPlace < iOrigPlaces; iOrigPlace++) {
	  SET_VECTOR_ELT(sexpMarking, iVoxelCurrentPlace, sexpTmp = allocVector(REALSXP, iTotalSteps));
	  ppdMarking[iVoxelCurrentPlace] = REAL(sexpTmp);
	  
	  pdCrntMarking[iVoxelCurrentPlace] = pdM[iVoxel + iVoxels * iOrigPlace];
	  
	  iVoxelCurrentPlace++;
	}
      }
      UNPROTECT_PTR(M);
    }
    
    double dAcumHazard = 0;
    for(iVoxelTransition = 0; iVoxelTransition < iVoxelSlowTransitions; iVoxelTransition++) {
      pdHazard[iVoxelTransition] = 0;
    }

    double dTime, dTarget = 0;
    int iTotTransitions = 0;
    
    double dIncr = MIN_INCR, dStartIncr;
    double dAbsoluteMaxError = 0;

    int iStep = 0;
    int iAcceptedIncr = 0, iRejectedIncr = 0, iInterruptCnt = 100000;
    double dNewHazard = 0;
#ifdef RB_TIME
    clkStep = clock();
    iElSteps = 0;
#endif
    do {
      if (pCManage_time || sexpRManage_time) {
	double dEnd = 0;
	if (pCManage_time) {
	  dEnd = ((double(*)(double, rb_double *)) pCManage_time)(dTarget, pdCrntMarking);
	} else {
	  defineVar(install("y"), sexpCrntMarking, rho);
	  PROTECT(sexpTmp = allocVector(REALSXP, 1));
	  *REAL(sexpTmp) = dTarget;
	  defineVar(install("StartTime"), sexpTmp, rho);
	  UNPROTECT_PTR(sexpTmp);
	  
	  sexpTmp = eval(sexpRManage_time, rho);
	  dEnd = *REAL(VECTOR_ELT(sexpTmp,0));
	  for(iVoxelPlace = 0; iVoxelPlace < iVoxelPlaces; iVoxelPlace++) {
	    pdCrntMarking[iVoxelPlace] = REAL(VECTOR_ELT(sexpTmp,1))[iVoxelPlace];
	  }
	}
	iSectionSteps = (int)(dEnd / dDelta) + 1;
      }
      
      for(iVoxelPlace = 0; iVoxelPlace < iVoxelPlaces; iVoxelPlace++) {
	ppdMarking[iVoxelPlace][iStep] = pdBakCrntMarking[iVoxelPlace] = pdTmpCrntMarking[iVoxelPlace] = pdCrntMarking[iVoxelPlace];
      }
      
      dTime = dTarget;
      dTarget += dDelta;

      dStartIncr = dIncr;
      
      // For the calculation of all hazards...
      int *piHazardsToMod = piHazardsToMod_start;
      int iHazardsToMod_count = iVoxelSlowTransitions;

	/*
      dAcumHazard = 0;
      for(iVoxelTransition = 0; iVoxelTransition < iVoxelSlowTransitions; iVoxelTransition++) {
	pdHazard[iVoxelTransition] = 0;
      }
	*/
      
      do {    
	// Get hazards only for the transitions associated with
	// places whose quantities changed in the last step.
	for(iVoxelTransitionIdx = 0; iVoxelTransitionIdx < iHazardsToMod_count; iVoxelTransitionIdx++) {
	  iVoxelTransition = piHazardsToMod[iVoxelTransitionIdx];
	  transition_el *pThisVoxelTransition = &pVoxelTransition[iVoxelTransition];

	  switch(pThisVoxelTransition->cType_of_h) {
	  case HZ_DOUBLE:
	    dNewHazard = pThisVoxelTransition->h.dTrCnt;
	    for(iVoxelPlaceIdx = 0; iVoxelPlaceIdx < pThisVoxelTransition->cPreNZ_len; iVoxelPlaceIdx++) {
	      char cPre = pThisVoxelTransition->pcPreNZ_value[iVoxelPlaceIdx];
	      // NOTE!!! Decide below between int or double
	      int iMarking = pdCrntMarking[pThisVoxelTransition->piPreNZ_VoxelPlace[iVoxelPlaceIdx]];
	      if (iMarking < cPre) {
		dNewHazard = 0;
		goto update_hazard_1;
	      }
	      for (k = 0; k < cPre; k++)
		dNewHazard *= (iMarking - k) / (double)(k+1);
	    }	    
	    break;
	  case HZ_CFUNCTION:	
	    dNewHazard = ((double(*)(double, rb_double *)) pThisVoxelTransition->h.pDL_FUNC)(dTime, pdCrntMarking + pThisVoxelTransition->iLocalVoxelPlace_offset);
	    break;
	  case HZ_RFUNCTION:
	    defineVar(install("y"), sexpCrntMarking, rho);
	    //dNewHazard = REAL(eval(VECTOR_ELT(sexpFunction, iVoxelTransition), rho))[0];
	    break;
	  }
	update_hazard_1:
	  dAcumHazard += dNewHazard - pdHazard[iVoxelTransition];
	  pdHazard[iVoxelTransition] = dNewHazard;
	}
	
	double dLogRandom;
	if (iVoxelSlowTransitions) {
	  dLogRandom = log(unif_rand());
	  iUsedRandomNumbers++;
	} else
	  dLogRandom = -DBL_MAX;
	double dIntHazard = 0;
	
	dIncr = dStartIncr;
	//dIncr = 1;
	dStartIncr = -1;
	
	int iRKs;
	double dEvalTime = 0;
	int iContinue = TRUE;
	int iStartRK = 1;
	do {
	  if (dTime + dIncr > dTarget) {
	    dIncr = dTarget - dTime;
	  }
	  for(iRKs = iStartRK; iRKs < 8; iRKs++) {
	    switch (iRKs) {
	    case 1:
	      dEvalTime = dTime;
	      pdYdot = pdK1;
	      for(iVoxelPlaceIdx = 0; iVoxelPlaceIdx < iVoxelFastPlaceIdx_len; iVoxelPlaceIdx++) {
		iVoxelPlace = piVoxelFastPlaceIdx[iVoxelPlaceIdx];
		pdYdot[iVoxelPlace] = 0;
		pdTmpCrntMarking[iVoxelPlace] = pdCrntMarking[iVoxelPlace];
	      }
	      break;
	    case 2:
	      dEvalTime = dTime + c2 * dIncr;
	      pdYdot = pdK2;
	      for(iVoxelPlaceIdx = 0; iVoxelPlaceIdx < iVoxelFastPlaceIdx_len; iVoxelPlaceIdx++) {
		iVoxelPlace = piVoxelFastPlaceIdx[iVoxelPlaceIdx];
		pdYdot[iVoxelPlace] = 0;
		pdTmpCrntMarking[iVoxelPlace] = pdCrntMarking[iVoxelPlace] +
		  a21*pdK1[iVoxelPlace] * dIncr;
	      }
	      break;
	    case 3:
	      dEvalTime = dTime + c3 * dIncr;
	      pdYdot = pdK3;
	      for(iVoxelPlaceIdx = 0; iVoxelPlaceIdx < iVoxelFastPlaceIdx_len; iVoxelPlaceIdx++) {
		iVoxelPlace = piVoxelFastPlaceIdx[iVoxelPlaceIdx];
		pdYdot[iVoxelPlace] = 0;
		pdTmpCrntMarking[iVoxelPlace] = pdCrntMarking[iVoxelPlace] +
		  (a31*pdK1[iVoxelPlace] + a32*pdK2[iVoxelPlace]) * dIncr;
	      }
	      break;
	    case 4:
	      dEvalTime = dTime + c4 * dIncr;
	      pdYdot = pdK4;
	      for(iVoxelPlaceIdx = 0; iVoxelPlaceIdx < iVoxelFastPlaceIdx_len; iVoxelPlaceIdx++) {
		iVoxelPlace = piVoxelFastPlaceIdx[iVoxelPlaceIdx];
		pdYdot[iVoxelPlace] = 0;
		pdTmpCrntMarking[iVoxelPlace] = pdCrntMarking[iVoxelPlace] +
		  (a41*pdK1[iVoxelPlace] + a42*pdK2[iVoxelPlace] + a43*pdK3[iVoxelPlace]) * dIncr;
	      }
	      break;
	    case 5:
	      dEvalTime = dTime + c5 * dIncr;
	      pdYdot = pdK5;
	      for(iVoxelPlaceIdx = 0; iVoxelPlaceIdx < iVoxelFastPlaceIdx_len; iVoxelPlaceIdx++) {
		iVoxelPlace = piVoxelFastPlaceIdx[iVoxelPlaceIdx];
		pdYdot[iVoxelPlace] = 0;
		pdTmpCrntMarking[iVoxelPlace] = pdCrntMarking[iVoxelPlace] +
		  (a51*pdK1[iVoxelPlace] + a52*pdK2[iVoxelPlace] + a53*pdK3[iVoxelPlace] +
		   a54*pdK4[iVoxelPlace]) * dIncr;
	      }
	      break;
	    case 6:
	      dEvalTime = dTime + c6 * dIncr;
	      pdYdot = pdK6;
	      for(iVoxelPlaceIdx = 0; iVoxelPlaceIdx < iVoxelFastPlaceIdx_len; iVoxelPlaceIdx++) {
		iVoxelPlace = piVoxelFastPlaceIdx[iVoxelPlaceIdx];
		pdYdot[iVoxelPlace] = 0;
		pdTmpCrntMarking[iVoxelPlace] = pdCrntMarking[iVoxelPlace] +
		  (a61*pdK1[iVoxelPlace] + a62*pdK2[iVoxelPlace] + a63*pdK3[iVoxelPlace] +
		   a64*pdK4[iVoxelPlace] + a65*pdK5[iVoxelPlace]) * dIncr;
	      }
	      break;
	    case 7:
	      dEvalTime = dTime + c7 * dIncr;
	      pdYdot = pdK7;
	      for(iVoxelPlaceIdx = 0; iVoxelPlaceIdx < iVoxelFastPlaceIdx_len; iVoxelPlaceIdx++) {
		iVoxelPlace = piVoxelFastPlaceIdx[iVoxelPlaceIdx];
		pdYdot[iVoxelPlace] = 0;
		pdTmpCrntMarking[iVoxelPlace] = pdCrntMarking[iVoxelPlace] +
		  (a71*pdK1[iVoxelPlace] + a73*pdK3[iVoxelPlace] + a74*pdK4[iVoxelPlace] +
		   a75*pdK5[iVoxelPlace] + a76*pdK6[iVoxelPlace]) * dIncr;
	      }
	      break;
	    }

	    for(iVoxelTransition = 0; iVoxelTransition < iVoxelFastTransitions; iVoxelTransition++) {
	      transition_el *pThisVoxelTransition = &pVoxelFastTransition[iVoxelTransition];
	      switch(pThisVoxelTransition->cType_of_h) {
	      case HZ_DOUBLE:
		dNewHazard = pThisVoxelTransition->h.dTrCnt;
		for(iVoxelPlaceIdx = 0; iVoxelPlaceIdx < pThisVoxelTransition->cPreNZ_len; iVoxelPlaceIdx++) {
		  for (k = 0; k < pThisVoxelTransition->pcPreNZ_value[iVoxelPlaceIdx]; k++)
		    dNewHazard *= pdTmpCrntMarking[pThisVoxelTransition->piPreNZ_VoxelPlace[iVoxelPlaceIdx]];
		}
		/*
		if (dNewHazard < dEct)
		  continue;
		*/
		break;
	      case HZ_CFUNCTION:	
		dNewHazard = ((double(*)(double, rb_double *)) pThisVoxelTransition->h.pDL_FUNC)(dEvalTime, pdTmpCrntMarking + pThisVoxelTransition->iLocalVoxelPlace_offset);
		break;
	      case HZ_RFUNCTION:
		// defineVar(install("y"), sexpCrntMarking, rho);
		//dNewHazard = REAL(eval(VECTOR_ELT(sexpFunction, iVoxelTransition), rho))[0];
		break;
	      }

	      for(iVoxelPlaceIdx = 0; iVoxelPlaceIdx < pThisVoxelTransition->cSNZ_len; iVoxelPlaceIdx++) {
		pdYdot[pThisVoxelTransition->piSNZ_VoxelPlace[iVoxelPlaceIdx]] += pThisVoxelTransition->pcSNZ_value[iVoxelPlaceIdx] * dNewHazard;
	      }
	    }
	  }
	  rb_double dInfNormError = 0, dInfNormMarking = 0;
	  int iInfNormErrorVoxelPlace = -1;
	  for(iVoxelPlaceIdx = 0; iVoxelPlaceIdx < iVoxelFastPlaceIdx_len; iVoxelPlaceIdx++) {
	    iVoxelPlace = piVoxelFastPlaceIdx[iVoxelPlaceIdx];
	    rb_double dYj = 
	      (b41*pdK1[iVoxelPlace] + b43*pdK3[iVoxelPlace] + b44*pdK4[iVoxelPlace] + b45*pdK5[iVoxelPlace] + b46*pdK6[iVoxelPlace] + b47*pdK7[iVoxelPlace]) * dIncr;
	    rb_double dZj = 
	      (pdCrntDiffMarking[iVoxelPlace] = (b51*pdK1[iVoxelPlace] + b53*pdK3[iVoxelPlace] + b54*pdK4[iVoxelPlace] + b55*pdK5[iVoxelPlace] + b56*pdK6[iVoxelPlace]) * dIncr);
	    
	    rb_double dThisError;
	    if ((dThisError = dYj-dZj) < 0.)
	      dThisError = -dThisError;
	    if (dThisError > dInfNormError) {
	      dInfNormError = dThisError;
	      iInfNormErrorVoxelPlace = iVoxelPlace;
	    }
	    rb_double dThisMarking;	    
	    if ((dThisMarking = pdCrntMarking[iVoxelPlace]) < 0.)
	      dThisMarking = -dThisMarking;
	    if (dThisMarking > dInfNormMarking)
	      dInfNormMarking = dThisMarking;
	  }
	  if (dInfNormMarking > 1.)
	    dInfNormMarking = 1.;
	  rb_double dTau;
	  if ((dTau = dEct*dInfNormMarking) == 0.)
	    dTau = MIN_INCR;

	  if (dInfNormError == 0.) {
	    if (dTau == MIN_INCR)
	      dInfNormError = MIN_INCR*1e-4;
	    else
	      dInfNormError = MIN_INCR;
	  }

	  piInfNormErrorXVoxelPlace[iInfNormErrorVoxelPlace]++;
	  if (dInfNormError > dTau) {
	    // Current increment is rejected.
	    // Try a new one and retry integration step
	    dIncr = .8 * dIncr * R_pow(dTau/dInfNormError,err_pow);
	    iRejectedIncr++;
	    continue;
	  }

	  iAcceptedIncr++;
	  if (dInfNormError > dAbsoluteMaxError)
	    dAbsoluteMaxError = dInfNormError;

	  for(iVoxelPlaceIdx = 0; iVoxelPlaceIdx < iVoxelFastPlaceIdx_len; iVoxelPlaceIdx++) {
	    iVoxelPlace = piVoxelFastPlaceIdx[iVoxelPlaceIdx];
	    pdBakCrntMarking[iVoxelPlace] = pdCrntMarking[iVoxelPlace];
	    pdCrntMarking[iVoxelPlace] += pdCrntDiffMarking[iVoxelPlace];
	    pdK1[iVoxelPlace] = pdK7[iVoxelPlace];
	  }
	  if (iStartRK == 1)
	    iStartRK = 2;
	  
	  double dPrevAcumHazard = dAcumHazard;
	  for(iVoxelTransitionIdx = 0; iVoxelTransitionIdx < iSlowHazardsToModFromFast_len; iVoxelTransitionIdx++) {
	    iVoxelTransition = piSlowHazardsToModFromFast[iVoxelTransitionIdx];
	    transition_el *pThisVoxelTransition = &pVoxelTransition[iVoxelTransition];

	    switch(pThisVoxelTransition->cType_of_h) {
	    case HZ_DOUBLE:
	      dNewHazard = pThisVoxelTransition->h.dTrCnt;
	      for(iVoxelPlaceIdx = 0; iVoxelPlaceIdx < pThisVoxelTransition->cPreNZ_len; iVoxelPlaceIdx++) {
		char cPre = pThisVoxelTransition->pcPreNZ_value[iVoxelPlaceIdx];
		// NOTE!!! Decide below between int or double
		int iMarking = pdCrntMarking[pThisVoxelTransition->piPreNZ_VoxelPlace[iVoxelPlaceIdx]];
		if (iMarking < cPre) {
		  dNewHazard = 0;
		  goto update_hazard_2;
		}
		for (k = 0; k < cPre; k++)
		  dNewHazard *= (iMarking - k) / (double)(k+1);
	      }	    
	      break;
	    case HZ_CFUNCTION:	
	      dNewHazard = ((double(*)(double, rb_double *)) pThisVoxelTransition->h.pDL_FUNC)(dTime, pdCrntMarking + pThisVoxelTransition->iLocalVoxelPlace_offset);
	      break;
	    case HZ_RFUNCTION:
	      defineVar(install("y"), sexpCrntMarking, rho);
	      //dNewHazard = REAL(eval(VECTOR_ELT(sexpFunction, iVoxelTransition), rho))[0];
	      break;
	    }
	update_hazard_2:
	    dAcumHazard += dNewHazard - (pdBakHazard[iVoxelTransition] = pdHazard[iVoxelTransition]);
	    pdHazard[iVoxelTransition] = dNewHazard;
	  }
	  
	  double dIncrIntHazard = (dPrevAcumHazard + dAcumHazard) * dIncr / 2;
	  double dDiff = dIntHazard + dIncrIntHazard + dLogRandom;
	  if (fabs(dDiff) < dEct) {
	    // Next check is needed because once in a while
	    // you will end here on a first attempt
	    // and not after having entered the else below
	    // on the previous iteration, and you want to
	    // avoid having dStartIncr < 0 which will be
	    // assigned to dIncr for the next cycle...
	    if (dStartIncr < 0)
	      dStartIncr = dIncr;
	    iContinue = FALSE;
	    dIntHazard += dIncrIntHazard;
	  } else if (dDiff < 0) {
	    dIntHazard += dIncrIntHazard;
	  } else {
	    if (dStartIncr < 0)
	      dStartIncr = dIncr;
	    dIncr *= - (dIntHazard + dLogRandom) / dIncrIntHazard;
	    
	    // Invalidate results
	    for(iVoxelPlaceIdx = 0; iVoxelPlaceIdx < iVoxelFastPlaceIdx_len; iVoxelPlaceIdx++) {
	      iVoxelPlace = piVoxelFastPlaceIdx[iVoxelPlaceIdx];
	      pdCrntMarking[iVoxelPlace] = pdBakCrntMarking[iVoxelPlace];
	    }
	    dAcumHazard = dPrevAcumHazard;
	    for(iVoxelTransitionIdx = 0; iVoxelTransitionIdx < iSlowHazardsToModFromFast_len; iVoxelTransitionIdx++) {
	      iVoxelTransition = piSlowHazardsToModFromFast[iVoxelTransitionIdx];
	      pdHazard[iVoxelTransition] = pdBakHazard[iVoxelTransition];
	    }
	    iStartRK = 1;
	    continue;
	  }

	  // Update clock according to the last accepted increment
	  dTime += dIncr;

	  // Check if current state needs to be saved
	  if (dTime == dTarget) {
	    ++iStep;
	    // Update the state for the fixed incremented time.
	    for(iVoxelPlace = 0; iVoxelPlace < iVoxelPlaces; iVoxelPlace++) {
	      ppdMarking[iVoxelPlace][iStep] = pdCrntMarking[iVoxelPlace];
	    }
	    if (iStep == iSectionSteps - 1)
	      goto EXIT_LOOP;
#ifdef RB_TIME
	    iElSteps++;
	    if ((iStep % iTotalStepsOver50) == 0) {
	      clock_t clkTmp2 = clock();
	      double dSecondsToFinish = (double) (clkTmp2- clkStep)/CLOCKS_PER_SEC/iElSteps*(iTotalSteps-iStep);
	      Rprintf ("\t");
	      PrintfTime(dSecondsToFinish);
	      clkStep = clkTmp2;
	      iElSteps = 0;
	    }
	    if ((iStep % iTotalStepsOver10) == 0) {
	      Rprintf(".\n");
	    }
#endif
	    dTarget += dDelta;
	    
	    // Force check if user interrupted
	    iInterruptCnt = 1;
	  }
	  if (! --iInterruptCnt) {
	    // Allow user interruption
	    R_CheckUserInterrupt();
	    iInterruptCnt = 100000;
	  }
	  // Set new increment for next integration step
	  dIncr = .8 * dIncr * R_pow(dTau/dInfNormError,err_pow);
	} while (iContinue);

	while (!--iTillResort) {
	  quicksort(piOrderedTransition, piTotTransitions, 0, iVoxelSlowTransitions-1);
	  switch (iTotResort++) {
	  case 0:
	    iTillResort = 10000;
	    break;
	  case 1:
	    iTillResort = 100000;
	    break;
	  default:
	    iTillResort = 1000000;
	  }
	}
	double dPartialAcumHazard = 0;
	// Find out which transition happened
	double dRnd = runif(0, dAcumHazard);
	iUsedRandomNumbers++;
	for(iVoxelTransitionIdx = 0; iVoxelTransitionIdx < iVoxelSlowTransitions; iVoxelTransitionIdx++) {
	  iVoxelTransition = piOrderedTransition[iVoxelTransitionIdx];
	  if (dRnd < (dPartialAcumHazard += pdHazard[iVoxelTransition])) {
	    piTotTransitions[iVoxelTransition]++;
	    
	    transition_el *pThisVoxelTransition = &pVoxelTransition[iVoxelTransition];
	    piHazardsToMod = pThisVoxelTransition->piHazardsToMod;
	    iHazardsToMod_count = pThisVoxelTransition->iHazardsToMod_count;
	    for(iVoxelPlaceIdx = 0; iVoxelPlaceIdx < pThisVoxelTransition->cSNZ_len; iVoxelPlaceIdx++) {
	      iVoxelPlace = pThisVoxelTransition->piSNZ_VoxelPlace[iVoxelPlaceIdx];
	      // Update the state
	      if ((pdCrntMarking[iVoxelPlace] += pThisVoxelTransition->pcSNZ_value[iVoxelPlaceIdx]) < 0)
		pdCrntMarking[iVoxelPlace] = 0;
	      pdBakCrntMarking[iVoxelPlace] = pdTmpCrntMarking[iVoxelPlace] = pdCrntMarking[iVoxelPlace];
	    }
	    break;
	  }
	}
	++iTotTransitions;
      } while (TRUE);
    EXIT_LOOP:;
      Rprintf("|\n");
    } while (iSectionSteps < iTotalSteps);
    iTotalUsedRandomNumbers += iUsedRandomNumbers;
    Rprintf("\t%d\t%d\t%d\t%e\t%d\t%d", iTotTransitions, iUsedRandomNumbers, iTotalUsedRandomNumbers, dAbsoluteMaxError, iAcceptedIncr, iRejectedIncr);
    
#ifdef RB_SUBTIME
    c1 = clock();
    Rprintf ("\t To go: ");
    PrintfTime((double) (c1 - c0)/CLOCKS_PER_SEC/(iRun+1)*(iRuns-iRun-1));
#endif
    Rprintf("\n");

    SEXP sexpTotTransitions;
    PROTECT(sexpTotTransitions = allocVector(INTSXP, 1));
    INTEGER(sexpTotTransitions)[0] = iTotTransitions;
    
    SEXP sexpUsedRandomNumbers;
    PROTECT(sexpUsedRandomNumbers = allocVector(INTSXP, 1));
    INTEGER(sexpUsedRandomNumbers)[0] = iUsedRandomNumbers;
    
    SEXP sexpThisRun;
#ifdef RB_SAVE_INCR_INFO
    if (iRun >= 10)
      PROTECT(sexpThisRun = allocVector(VECSXP, 5));
    else
      PROTECT(sexpThisRun = allocVector(VECSXP, 10));
#else
    PROTECT(sexpThisRun = allocVector(VECSXP, 5));
#endif
    
    SET_VECTOR_ELT(sexpThisRun, 0, sexpMarking);
    UNPROTECT_PTR(sexpMarking);
    SET_VECTOR_ELT(sexpThisRun, 1, sexpTotXTransition);
    UNPROTECT_PTR(sexpTotXTransition);
    SET_VECTOR_ELT(sexpThisRun, 2, sexpTotTransitions);
    UNPROTECT_PTR(sexpTotTransitions);
    SET_VECTOR_ELT(sexpThisRun, 3, sexpUsedRandomNumbers);
    UNPROTECT_PTR(sexpUsedRandomNumbers);
    SET_VECTOR_ELT(sexpThisRun, 4, sexpInfNormErrorXVoxelPlace);
    UNPROTECT_PTR(sexpInfNormErrorXVoxelPlace);

    //    PrintValue(sexpThisRun);
    /*
    int *pttt = (int *) R_alloc(2000000, sizeof(int));
    PROTECT(sexpTmp = allocVector(REALSXP, 2000000));
      UNPROTECT_PTR(sexpTmp);
    */
#ifdef RB_SAVE_INCR_INFO
    if (iRun < 10) {
      SEXP sexpTmp;
      double *pdTmp;

      PROTECT(sexpTmp = allocVector(REALSXP, iTotIncrTime));
      pdTmp = REAL(sexpTmp);
      int i;
      for (i = 0; i < iTotIncrTime; i++)
	pdTmp[i] = pdIncr[i];
      SET_VECTOR_ELT(sexpThisRun, 5, sexpTmp);
      UNPROTECT_PTR(sexpTmp);
      
      PROTECT(sexpTmp = allocVector(REALSXP, iTotIncrTime));
      pdTmp = REAL(sexpTmp);
      for (i = 0; i < iTotIncrTime; i++)
	pdTmp[i] = pdIncrTime[i];
      SET_VECTOR_ELT(sexpThisRun, 6, sexpTmp);
      UNPROTECT_PTR(sexpTmp);
      
      PROTECT(sexpTmp = allocVector(REALSXP, iTotIntHazardTime));
      pdTmp = REAL(sexpTmp);
      for (i = 0; i < iTotIntHazardTime; i++)
	pdTmp[i] = pdAcumHazard[i];
      SET_VECTOR_ELT(sexpThisRun, 7, sexpTmp);
      UNPROTECT_PTR(sexpTmp);
      
      PROTECT(sexpTmp = allocVector(REALSXP, iTotIntHazardTime));
      pdTmp = REAL(sexpTmp);
      for (i = 0; i < iTotIntHazardTime; i++)
	pdTmp[i] = pdIntHazard[i];
      SET_VECTOR_ELT(sexpThisRun, 8, sexpTmp);
      UNPROTECT_PTR(sexpTmp);
      
      PROTECT(sexpTmp = allocVector(REALSXP, iTotIntHazardTime));
      pdTmp = REAL(sexpTmp);
      for (i = 0; i < iTotIntHazardTime; i++)
	pdTmp[i] = pdIntHazardTime[i];
      SET_VECTOR_ELT(sexpThisRun, 9, sexpTmp);
      UNPROTECT_PTR(sexpTmp);
    }
#endif
    //    PrintValue(sexpThisRun);
    
    SEXP sexpNames;
#ifdef RB_SAVE_INCR_INFO
    if (iRun >= 10)
      PROTECT(sexpNames = allocVector(VECSXP, 5));
    else
      PROTECT(sexpNames = allocVector(VECSXP, 10));
#else
    PROTECT(sexpNames = allocVector(VECSXP, 5));
#endif
    SET_VECTOR_ELT(sexpNames, 0, mkChar("M"));
    SET_VECTOR_ELT(sexpNames, 1, mkChar("transitions"));
    SET_VECTOR_ELT(sexpNames, 2, mkChar("tot.transitions"));
    SET_VECTOR_ELT(sexpNames, 3, mkChar("tot.rnd"));
    SET_VECTOR_ELT(sexpNames, 4, mkChar("max.inf.norm"));
#ifdef RB_SAVE_INCR_INFO
    if (iRun < 10) {
      SET_VECTOR_ELT(sexpNames, 5, mkChar("incr"));
      SET_VECTOR_ELT(sexpNames, 6, mkChar("incr.time"));
      SET_VECTOR_ELT(sexpNames, 7, mkChar("hazard"));
      SET_VECTOR_ELT(sexpNames, 8, mkChar("int.hazard"));
      SET_VECTOR_ELT(sexpNames, 9, mkChar("int.hazard.time"));
    }
#endif
    setAttrib(sexpThisRun, R_NamesSymbol, sexpNames);
    UNPROTECT_PTR(sexpNames);

    SET_VECTOR_ELT(sexpRun, iRun, sexpThisRun);
    UNPROTECT_PTR(sexpThisRun);
  }
  PutRNGstate();

#ifdef RB_MEMORY
  Rprintf("Total Memory used: %.2f MB\n", dUsedMemAcum);
#endif

  SEXP sexpAns;
  PROTECT(sexpAns = allocVector(VECSXP, 4));
  SET_VECTOR_ELT(sexpAns, 0, place);
  SET_VECTOR_ELT(sexpAns, 1, transition);
  SET_VECTOR_ELT(sexpAns, 2, sexpD_time);
  UNPROTECT_PTR(sexpD_time);
  SET_VECTOR_ELT(sexpAns, 3, sexpRun);
  UNPROTECT_PTR(sexpRun);

  SEXP sexpNames;
  PROTECT(sexpNames = allocVector(VECSXP, 4));
  SET_VECTOR_ELT(sexpNames, 0, mkChar("place"));
  SET_VECTOR_ELT(sexpNames, 1, mkChar("transition"));
  SET_VECTOR_ELT(sexpNames, 2, mkChar("dt"));
  SET_VECTOR_ELT(sexpNames, 3, mkChar("run"));
  setAttrib(sexpAns, R_NamesSymbol, sexpNames);
  UNPROTECT_PTR(sexpNames);

#ifdef RB_TIME
  c1 = clock();
  double dCpuTime = (double) (c1 - c0)/CLOCKS_PER_SEC;
  Rprintf ("Elapsed CPU time: ");
  PrintfTime(dCpuTime);
  Rprintf ("\t(%fs)\n", dCpuTime);
#endif

  if (sexpRManage_time)
    UNPROTECT_PTR(sexpRManage_time);
  //UNPROTECT_PTR(sexpFunction);
  UNPROTECT_PTR(sexpMarkingRowNames);
  UNPROTECT_PTR(sexpTmpCrntMarking);
  UNPROTECT_PTR(sexpCrntMarking);
  UNPROTECT_PTR(sexpAns);

  return(sexpAns);
}
