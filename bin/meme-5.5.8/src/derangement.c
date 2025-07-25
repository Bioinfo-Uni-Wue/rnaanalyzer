/********************************************************************
 * FILE: tgene.c
 * AUTHOR: Timothy Bailey
 * CREATE DATE: 06/05/2024
 * PROJECT: MEME suite
 * COPYRIGHT: 2024, Timothy Bailey
 ********************************************************************/
#ifdef MAIN
#define DEFINE_GLOBALS
#endif
#include <stdlib.h>
#include "utils.h"

//
// Return then number of derangements of the N integers from 0 to N-1.
//
int get_num_derangements(int n) {
  int i;

  // Easy cases.
  if (n==0) return(0);
  if (n==1 || n==2) return(n-1);
  
  // Store previous two values for f(n).
  int p0 = 0;
  int p1 = 1;

  // Iteration:
  for (i=3; i<=n; i++) {
    int f_i = (i-1) * (p0 + p1);
    p0 = p1;
    p1 = f_i;
  }

  // Return result.
  return(p1);
} // get_num_derangements

//
// Get all derangements recursively.
//
void get_derangements_helper(
  int n, 		// original N
  int value,		// value to add to derangement
  int *cur_derangement,	// current derangement
  int **derangements,	// final derangements
  int *n_der		// current number of final derangements
) {
  int i, j;

  for (i=0; i<n; i++) {
    // See if position i is unoccupied.
    if (i == value || cur_derangement[i] != -1) continue;
    // place value at position i
    cur_derangement[i] = value;
    if (value == n-1) {
      // Save final derangement
      for (j=0; j<n; j++) {
        derangements[*n_der][j] = cur_derangement[j];
      }
      (*n_der)++;
    } else {
      // recurse
      get_derangements_helper(n, value+1, cur_derangement, derangements, n_der);
    }
    cur_derangement[i] = -1;
  }
} // get_derangements_helper

//
// Return all possible derangements of the N integers from 0 to N-1
// as a list of lists.
//

int **get_derangements(
  int n,
  int num_der
) {
  int i;

  // Create the array to hold the final derangements.
  int **derangements = mm_malloc(sizeof(int *) * num_der);
  for (i=0; i<num_der; i++) {
    derangements[i] = mm_malloc(sizeof(int) * n);
  }
  
  // Create the array for the current derangement and initialize it to -1's.
  int *cur_derangement = mm_malloc(sizeof(int *) * n);
  for (i=0; i<n; i++) cur_derangement[i] = -1;

  // Get derangements recursively.
  int n_der = 0;
  get_derangements_helper(n, 0, cur_derangement, derangements, &n_der);

  return(derangements);
} // get_derangements

#ifdef MAIN
/************************************************************************/
/*
	Compute the number of derangements of N integers from 0 to N-1.
*/
/************************************************************************/
#define BUFSIZE 100
int main(int argc, char **argv) {
  int n, i, j;

  if (argc != 2) {
    fprintf(stderr, "Usage: derangement <N>\n"
	"\tPrint all the derangements of N integers from 0 to N-1.\n"
    );
    return(1);
  }
  n = atoi(argv[1]);

  int num_der = get_num_derangements(n);
  printf("N: %d Number of Derangements: %d\n", n, num_der);

  int **derangements = get_derangements(n, num_der);

  // Print the derangements.
  for (i=0; i<num_der; i++) {
    fprintf(stderr, "derangement %d: ", i+1);
    for (j=0; j<n; j++) fprintf(stderr, " %d", derangements[i][j]);
    fprintf(stderr, "\n");
  }

  return(0);
} // main
#endif
