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

#include "helper.h"

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
  // Next int wastes 4 bytes per element.
  int iHazardsToMod_count;
  int iLocalVoxelPlace_offset;

};
typedef struct transition transition_el;

// This structure has 7 unused bytes per element.
struct tree_el {
  char iGroup;
  double dPartialAcumHazard;
  struct tree_el *parent, *left, *right;
};
typedef struct tree_el node;

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


SEXP sGillespieDirectCR(SEXP model, SEXP T, SEXP delta, SEXP runs,
			SEXP place, SEXP transition, SEXP rho)
{
  SEXP sexpTmp;
  int *piTmp;
  char cTmp;

  int iVoxelPlaceIdx;
  int iVoxelTransitionIdx;

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

    SEXP VoxelSubset;
    if ((VoxelSubset = getListElement(this_VoxelFamily, "VoxelSubset")) == R_NilValue)
      error("VoxelSubset not provided!\n");

    int iOrigTransition;
    for (iOrigTransition = 0; iOrigTransition < iOrigTransitions; iOrigTransition++) {      
      iVoxelTransitions += length(VECTOR_ELT(VoxelSubset, piTransitionInVoxelSubset[iOrigTransition] - 1));      
    }
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
	
	sexpTmp = VECTOR_ELT(JumpingPattern, iJumpingPattern);
	piTmp = INTEGER(getAttrib(sexpTmp, R_DimSymbol));
	int iJumps = piTmp[0];
	
	iVoxelTransitions += iJumps;      
      }
    }
  }
  Rprintf("VoxelTransitions: %d\tVoxelPlaces: %d\n", iVoxelTransitions, iVoxelPlaces);

  transition_el *pVoxelTransition = (transition_el *) R_alloc(iVoxelTransitions, sizeof(transition_el));
#ifdef RB_MEMORY
  dUsedMemAcum += (dThisMem = ((double) sizeof(transition_el) * iVoxelTransitions)/1e6);
  Rprintf("pVoxelTransition: %.2f MB (%d x %d bytes)\n", dThisMem, iVoxelTransitions, sizeof(transition_el));
#endif

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

	transition_el *pThisVoxelTransition = &pVoxelTransition[iVoxelCurrentTransition];

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

	sexpTmp = VECTOR_ELT(this_JumpingPlace, 1);
	PROTECT(sexpTmp = coerceVector(sexpTmp, REALSXP));
	double dH = REAL(sexpTmp)[0];
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

	  transition_el *pThisVoxelTransition = &pVoxelTransition[iVoxelCurrentTransition];
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
	  iVoxelCurrentTransition++;
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
  // Identify the transitions that need the hazards recalculated
  for(iVoxelTransition = 0; iVoxelTransition < iVoxelTransitions; iVoxelTransition++) {
    transition_el *pThisVoxelTransition = &pVoxelTransition[iVoxelTransition];
    int iHazardToCompTot = 0;
    for(iVoxelPlaceIdx = 0; iVoxelPlaceIdx < pThisVoxelTransition->cSNZ_len; iVoxelPlaceIdx++) {
      int iVoxelPlace = pThisVoxelTransition->piSNZ_VoxelPlace[iVoxelPlaceIdx];
      
      int iVoxelTransitionIdx2, iVoxelTransition2;
      for(iVoxelTransitionIdx2 = 0; iVoxelTransitionIdx2 < piVoxelPlacePreNZ_len[iVoxelPlace]; iVoxelTransitionIdx2++) {
	iVoxelTransition2 = ppiVoxelPlacePreNZ_VoxelTransition[iVoxelPlace][iVoxelTransitionIdx2];
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
    iHazardsToMod_count += iHazardToCompTot;
  }  
  int *piThisHazardsToMod = (int *) R_alloc(iHazardsToMod_count, sizeof(int));
#ifdef RB_MEMORY
  dUsedMemAcum += (dThisMem = ((double) sizeof(int) * iHazardsToMod_count)/1e6);
  Rprintf("piHazardsToMod elements: %.2f MB (%d x %d bytes)\n", dThisMem, iHazardsToMod_count, sizeof(int));
#endif
  // Identify the transitions that need the hazards recalculated
  for(iVoxelTransition = 0; iVoxelTransition < iVoxelTransitions; iVoxelTransition++) {
    transition_el *pThisVoxelTransition = &pVoxelTransition[iVoxelTransition];
    int iHazardToCompTot = 0;
    for(iVoxelPlaceIdx = 0; iVoxelPlaceIdx < pThisVoxelTransition->cSNZ_len; iVoxelPlaceIdx++) {
      int iVoxelPlace = pThisVoxelTransition->piSNZ_VoxelPlace[iVoxelPlaceIdx];

      int iVoxelTransitionIdx2, iVoxelTransition2;
      for(iVoxelTransitionIdx2 = 0; iVoxelTransitionIdx2 < piVoxelPlacePreNZ_len[iVoxelPlace]; iVoxelTransitionIdx2++) {
	iVoxelTransition2 = ppiVoxelPlacePreNZ_VoxelTransition[iVoxelPlace][iVoxelTransitionIdx2];
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
    pThisVoxelTransition->iHazardsToMod_count = iHazardToCompTot;
    pThisVoxelTransition->piHazardsToMod = piThisHazardsToMod;
    piThisHazardsToMod += iHazardToCompTot;

    for (iVoxelTransitionIdx = 0; iVoxelTransitionIdx < iHazardToCompTot; iVoxelTransitionIdx++) {
      pThisVoxelTransition->piHazardsToMod[iVoxelTransitionIdx] = piHazardsToMod_start[iVoxelTransitionIdx];
    }
  }

  // For the initial calculation of all hazards...
  for(iVoxelTransition = 0; iVoxelTransition < iVoxelTransitions; iVoxelTransition++) {
    piHazardsToMod_start[iVoxelTransition] = iVoxelTransition;
  }

  SEXP sexpCrntMarking;
  PROTECT(sexpCrntMarking = allocVector(REALSXP, iVoxelPlaces));
  double *pdCrntMarking = REAL(sexpCrntMarking);
  
  double dDelta = *REAL(delta);
  int iTotalSteps, iSectionSteps;
  double dT = 0;
  void *pCManage_time = 0;
  SEXP sexpRManage_time = 0;
  if (inherits(T, "NativeSymbol")) {
    pCManage_time = (void *) R_ExternalPtrAddr(T);
    dT = ((double(*)(double, double *)) pCManage_time)(-1, pdCrntMarking);
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
  
  double *pdHazard = (double *) R_alloc(iVoxelTransitions, sizeof(double));
  int iRun, iRuns = *INTEGER(runs);
  
  SEXP sexpRun;
  PROTECT(sexpRun = allocVector(VECSXP, iRuns));
  
  int iTotalUsedRandomNumbers = 0;
  
  // DiscTime Vector
  SEXP sexpD_time;
  PROTECT(sexpD_time = allocVector(REALSXP, iTotalSteps));
  double *pdDiscTime = REAL(sexpD_time);
  int k;
  double dTmp = 0;
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
  
  int iLevels = 7;
  int iGroups = pow(2, iLevels - 1);
  // Group holding the transitions that lie between boundaries
  int **ppiGroup = (int **) R_alloc(iGroups, sizeof(int *));
  // Number of transition each group has
  int *piGroupElm = (int *) R_alloc(iGroups, sizeof(int));
  // Total propensity hazard for each group
  int *piTotGroupTransitions = (int *) R_alloc(iGroups, sizeof(int));

  int iGroup, iGroupsWithElements = 0;
  for (iGroup = 0; iGroup < iGroups; iGroup++) {
    //ppiGroup[iGroup] = NULL;
    ppiGroup[iGroup] = (int *) R_alloc(iVoxelTransitions, sizeof(int));
  }

  node **ppnodeLevel = (node **) R_alloc(iLevels, sizeof(node *));
  int iLevel, iNode;
  int iNodesPerLevel = 1;
  for (iLevel = 0; iLevel < iLevels; iLevel++) {
    ppnodeLevel[iLevel] = (node *) R_alloc(iNodesPerLevel, sizeof(node));
    iNodesPerLevel *= 2;
  }
  node *pnodeRoot = &ppnodeLevel[0][0];
  pnodeRoot->parent = 0;
  node *pnodeGroup = ppnodeLevel[iLevels-1];

  iNodesPerLevel = 1;
  for (iLevel = 0; iLevel < iLevels; iLevel++) {
    for (iNode = 0; iNode < iNodesPerLevel; iNode++) {
      if (iLevel < iLevels-1) {
	ppnodeLevel[iLevel][iNode].iGroup = -1;
	ppnodeLevel[iLevel][iNode].left = &ppnodeLevel[iLevel+1][iNode*2];
	ppnodeLevel[iLevel][iNode].right = &ppnodeLevel[iLevel+1][iNode*2+1];
	ppnodeLevel[iLevel+1][iNode*2].parent = ppnodeLevel[iLevel+1][iNode*2+1].parent =
	  &ppnodeLevel[iLevel][iNode];
      } else {
	ppnodeLevel[iLevel][iNode].iGroup = iNode;
	ppnodeLevel[iLevel][iNode].left = ppnodeLevel[iLevel][iNode].right = 0;
      }
    }
    iNodesPerLevel *= 2;
  }

  double dNewHazard = 0;
  // Find minimum propensity
  double dMinHazard = DBL_MAX;
  for(iVoxelTransition = 0; iVoxelTransition < iVoxelTransitions; iVoxelTransition++) {
    transition_el *pThisVoxelTransition = &pVoxelTransition[iVoxelTransition];
    switch(pThisVoxelTransition->cType_of_h) {
    case HZ_DOUBLE:
      dNewHazard = pThisVoxelTransition->h.dTrCnt;
      int iVoxelPlaceIdx;
      for(iVoxelPlaceIdx = 0; 
	  iVoxelPlaceIdx < pThisVoxelTransition->cPreNZ_len; 
	  iVoxelPlaceIdx++) {
	for (k = 0; k < pThisVoxelTransition->pcPreNZ_value[iVoxelPlaceIdx]; k++)
	  dNewHazard *= (pThisVoxelTransition->pcPreNZ_value[iVoxelPlaceIdx] - k) / (double)(k+1);
      }
      if (dNewHazard > 0 && dNewHazard < dMinHazard)
	dMinHazard = dNewHazard;
      
      break;
    case HZ_CFUNCTION:	
      break;
    case HZ_RFUNCTION:
      break;
    }
  }
  
#ifdef RB_TIME
  c1 = clock();
  Rprintf ("Elapsed CPU time: %f\n", (float) (c1 - c0)/CLOCKS_PER_SEC);
  c0 = c1;
#endif

  GetRNGstate();
  for (iRun = 0; iRun < iRuns; iRun++) {
    
    int iUsedRandomNumbers = 0;
    Rprintf("%d ", iRun+1);
    
    // Totals for kind of transition vector
    SEXP sexpTotXTransition;
    PROTECT(sexpTotXTransition = allocVector(INTSXP, iVoxelTransitions));
    int *piTotTransitions = INTEGER(sexpTotXTransition);
    
    for(iVoxelTransition = 0; iVoxelTransition < iVoxelTransitions; iVoxelTransition++) {
      piTotTransitions[iVoxelTransition] = 0;
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
      
      Rprintf("\n");
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
    
    for(iVoxelTransition = 0; iVoxelTransition < iVoxelTransitions; iVoxelTransition++) {
      pdHazard[iVoxelTransition] = 0;
      pVoxelTransition[iVoxelTransition].cInGroup = -1;
    }

    for (iGroup = 0; iGroup < iGroups; iGroup++) {
      piGroupElm[iGroup] = 0;
      piTotGroupTransitions[iGroup] = 0;
    }
    
    iNodesPerLevel = 1;
    for (iLevel = 0; iLevel < iLevels; iLevel++) {
      for (iNode = 0; iNode < iNodesPerLevel; iNode++) {
	ppnodeLevel[iLevel][iNode].dPartialAcumHazard = 0;
      }
      iNodesPerLevel *= 2;
    }
    node *pnode;
    
    double dTime = 0, dTarget = 0;
    int iTotTransitions = 0;
    
    int iStep = 0;
    int iInterruptCnt = 10000000;
    double dNewHazard = 0;
#ifdef RB_TIME
    clkStep = clock();
    iElSteps = 0;
#endif
    do {
      if (pCManage_time || sexpRManage_time) {
	double dEnd = 0;
	if (pCManage_time) {
	  dEnd = ((double(*)(double, double *)) pCManage_time)(dTarget, pdCrntMarking);
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
	ppdMarking[iVoxelPlace][iStep] = pdCrntMarking[iVoxelPlace];
      }

      dTime = dTarget;
      dTarget += dDelta;
      
      // For the calculation of all hazards...
      int *piHazardsToMod = piHazardsToMod_start;
      int iHazardsToMod_count = iVoxelTransitions;
      
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
	      for (k = 0; k < pThisVoxelTransition->pcPreNZ_value[iVoxelPlaceIdx]; k++)
		dNewHazard *= (pdCrntMarking[pThisVoxelTransition->piPreNZ_VoxelPlace[iVoxelPlaceIdx]] - k) / (double)(k+1);
	    }	    
	    break;
	  case HZ_CFUNCTION:	
	    dNewHazard = ((double(*)(double, double *)) pThisVoxelTransition->h.pDL_FUNC)(dTime, pdCrntMarking + pThisVoxelTransition->iLocalVoxelPlace_offset);
	    break;
	  case HZ_RFUNCTION:
	    defineVar(install("y"), sexpCrntMarking, rho);
	    //dNewHazard = REAL(eval(VECTOR_ELT(sexpFunction, iVoxelTransition), rho))[0];
	    break;
	  }

	  double dDeltaHazard;
	  frexp(dNewHazard/dMinHazard, &iGroup);
	  if (iGroup-- > 0) {
	    // Transition belongs to a group
	    /*
	    // Check if the group has been initialized
	    if (!ppiGroup[iGroup]) {
	      ppiGroup[iGroup] = (int *) R_alloc(iVoxelTransitions, sizeof(int));
	      iGroupsWithElements++;
	    }
	    */
	    if (iGroup == pThisVoxelTransition->cInGroup) {
	      // Transitions will stay in same group as it was
	      dDeltaHazard = dNewHazard - pdHazard[iVoxelTransition];
	      pnode = &pnodeGroup[iGroup];
	      do {
		pnode->dPartialAcumHazard += dDeltaHazard;
	      } while ((pnode = pnode->parent));
	    } else if (pThisVoxelTransition->cInGroup != -1) {
	      // Transition was in another group and needs to be moved to the new one
	      int iOldGroup = pThisVoxelTransition->cInGroup;
	      int iOldPositionInGroup = pThisVoxelTransition->iPositionInGroup;
	      dDeltaHazard = -pdHazard[iVoxelTransition];
	      pnode = &pnodeGroup[iOldGroup];
	      do {
		pnode->dPartialAcumHazard += dDeltaHazard;
	      } while ((pnode = pnode->parent));
	      piGroupElm[iOldGroup]--; // Old group will have one less element
	      // Now, piGroupElm[iOldGroup] is the index to last transition in group
	      if (iOldPositionInGroup != piGroupElm[iOldGroup]) {
		// Transition is not the last in group,
		// put the last transition in place of the one to be removed
		ppiGroup[iOldGroup][iOldPositionInGroup] = 
		  ppiGroup[iOldGroup][piGroupElm[iOldGroup]];
		// Update position of previous last transition in group
		pVoxelTransition[ppiGroup[iOldGroup][iOldPositionInGroup]].iPositionInGroup =
		  iOldPositionInGroup;
	      }
	      dDeltaHazard = dNewHazard;
	      pnode = &pnodeGroup[iGroup];
	      do {
		pnode->dPartialAcumHazard += dDeltaHazard;
	      } while ((pnode = pnode->parent));
	      pThisVoxelTransition->cInGroup = iGroup;
	      pThisVoxelTransition->iPositionInGroup = piGroupElm[iGroup];
	      ppiGroup[iGroup][piGroupElm[iGroup]++] = iVoxelTransition;
	    } else if (pThisVoxelTransition->cInGroup == -1) { // Transition was in no group
	      dDeltaHazard = dNewHazard;
	      pnode = &pnodeGroup[iGroup];
	      do {
		pnode->dPartialAcumHazard += dDeltaHazard;
	      } while ((pnode = pnode->parent));
	      pThisVoxelTransition->cInGroup = iGroup;
	      pThisVoxelTransition->iPositionInGroup = piGroupElm[iGroup];
	      ppiGroup[iGroup][piGroupElm[iGroup]++] = iVoxelTransition;
	    } else {
	    error("ERROR: Option not considered 1\n");
	    }
	  } else if (pThisVoxelTransition->cInGroup != -1) {
	    // Transition will not belong to any group and needs to be removed from old
	    int iOldGroup = pThisVoxelTransition->cInGroup;
	    int iOldPositionInGroup = pThisVoxelTransition->iPositionInGroup;
	    dDeltaHazard = -pdHazard[iVoxelTransition];
	    pnode = &pnodeGroup[iOldGroup];
	    do {
	      pnode->dPartialAcumHazard += dDeltaHazard;
	    } while ((pnode = pnode->parent));
	    piGroupElm[iOldGroup]--; // Old group will have one less element
	    // Now, piGroupElm[iOldGroup] is the index to last transition in group
	    if (iOldPositionInGroup != piGroupElm[iOldGroup]) {
	      // Transition is not the last in group,
	      // put the last transition in place of the one to be removed
	      ppiGroup[iOldGroup][iOldPositionInGroup] = 
		ppiGroup[iOldGroup][piGroupElm[iOldGroup]];
	      // Update position of previous last transition in group
	      pVoxelTransition[ppiGroup[iOldGroup][iOldPositionInGroup]].iPositionInGroup =
		iOldPositionInGroup;
	    }
	    pThisVoxelTransition->cInGroup = -1;
	  }
	  pdHazard[iVoxelTransition] = dNewHazard;
	}
	
	// Get Time to transition
	if (pnodeRoot->dPartialAcumHazard > 0) {
	  dTime += exp_rand() / pnodeRoot->dPartialAcumHazard;
	  //dTime -= log(unif_rand()) / pnodeRoot->dPartialAcumHazard;
	  iUsedRandomNumbers++;
	} else
	  dTime = dT + 1;	

	while (dTime >= dTarget) {
	  ++iStep;
	  // Update the state for the fixed incremented time.
	  for(iVoxelPlace = 0; iVoxelPlace < iVoxelPlaces; iVoxelPlace++)
	    ppdMarking[iVoxelPlace][iStep] = pdCrntMarking[iVoxelPlace];
	  if (iStep == iSectionSteps - 1)
	    goto EXIT_LOOP;
#ifdef RB_TIME
	    iElSteps++;
	    if ((iStep % iTotalStepsOver50) == 0) {
	      double clkTmp2 = clock();
	      double dSecondsToFinish = (double) (clkTmp2 - clkStep)/CLOCKS_PER_SEC/iElSteps*(iTotalSteps-iStep);
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
	  iInterruptCnt = 10000000;
	}
	do {
	  // Find group containing firing transition
	  double dRnd = unif_rand() * pnodeRoot->dPartialAcumHazard;
	  iUsedRandomNumbers++;
	  pnode = pnodeRoot;
	  do {
	    if (dRnd < pnode->left->dPartialAcumHazard) {
	      pnode = pnode->left;
	    } else {
	      dRnd -= pnode->left->dPartialAcumHazard;
	      pnode = pnode->right;
	    }	      
	  } while (pnode->left);
	  // Next check is because
	  // once in a while it is generated a number that goes past
	  // the last group or selects a group with zero elements
	  // due to accumulated truncation errors.
	  // Discard this random number and try again.
	} while (piGroupElm[iGroup = pnode->iGroup] == 0);

	int iGroupElm = piGroupElm[iGroup];
	double dMaxInGroup = dMinHazard * pow(2, iGroup + 1);
	// Find transition in group
	while (1) {
	  if (! --iInterruptCnt) {
	    // Allow user interruption
	    R_CheckUserInterrupt();
	    iInterruptCnt = 10000000;
	  }
	  int iVoxelTransitionIdx = (int) (unif_rand() * iGroupElm);
	  iUsedRandomNumbers++;
	  iVoxelTransition = ppiGroup[iGroup][iVoxelTransitionIdx];
	  iUsedRandomNumbers++;
	  if (pdHazard[iVoxelTransition] > unif_rand() * dMaxInGroup) {
	    piTotTransitions[iVoxelTransition]++;
	    
	    transition_el *pThisVoxelTransition = &pVoxelTransition[iVoxelTransition];
	    piHazardsToMod = pThisVoxelTransition->piHazardsToMod;
	    iHazardsToMod_count = pThisVoxelTransition->iHazardsToMod_count;
	    for(iVoxelPlaceIdx = 0; iVoxelPlaceIdx < pThisVoxelTransition->cSNZ_len; iVoxelPlaceIdx++) {
	      
	      // Update the state
	      pdCrntMarking[pThisVoxelTransition->piSNZ_VoxelPlace[iVoxelPlaceIdx]] += pThisVoxelTransition->pcSNZ_value[iVoxelPlaceIdx];
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
    Rprintf("\t%d\t%d\t%d", iTotTransitions, iUsedRandomNumbers, iTotalUsedRandomNumbers);
#ifdef RB_SUBTIME
    c1 = clock();
    Rprintf ("\t To go: ");
    PrintfTime((double) (c1 - c0)/CLOCKS_PER_SEC/(iRun+1)*(iRuns-iRun-1));
#endif
    Rprintf ("\n");
    
    SEXP sexpTotTransitions;
    PROTECT(sexpTotTransitions = allocVector(INTSXP, 1));
    INTEGER(sexpTotTransitions)[0] = iTotTransitions;

    SEXP sexpThisRun;
    PROTECT(sexpThisRun = allocVector(VECSXP, 3));

    SET_VECTOR_ELT(sexpThisRun, 0, sexpMarking);
    UNPROTECT_PTR(sexpMarking);
    SET_VECTOR_ELT(sexpThisRun, 1, sexpTotXTransition);
    UNPROTECT_PTR(sexpTotXTransition);
    SET_VECTOR_ELT(sexpThisRun, 2, sexpTotTransitions);
    UNPROTECT_PTR(sexpTotTransitions);

    SEXP sexpNames;
    PROTECT(sexpNames = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(sexpNames, 0, mkChar("M"));
    SET_VECTOR_ELT(sexpNames, 1, mkChar("transitions"));
    SET_VECTOR_ELT(sexpNames, 2, mkChar("tot.transitions"));
    setAttrib(sexpThisRun, R_NamesSymbol, sexpNames);
    UNPROTECT_PTR(sexpNames);

    SET_VECTOR_ELT(sexpRun, iRun, sexpThisRun);
    UNPROTECT_PTR(sexpThisRun);
  }
  PutRNGstate();

#ifdef RB_MEMORY
  dUsedMemAcum += (dThisMem = ((double) sizeof(int) * iGroupsWithElements * iVoxelTransitions)/1e6);
  Rprintf("ppiGroup elements (%d/%d): %.2f MB (%d x %d bytes)\n", iGroupsWithElements, iGroups, dThisMem, iGroupsWithElements * iVoxelTransitions, sizeof(int));
#endif
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
  //  UNPROTECT_PTR(sexpFunction);
  UNPROTECT_PTR(sexpMarkingRowNames);
  UNPROTECT_PTR(sexpCrntMarking);
  UNPROTECT_PTR(sexpAns);
  return(sexpAns);
   
}
