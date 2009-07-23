/*  File: fastqcheck.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) Genome Research Limited, 2006
 *-------------------------------------------------------------------
 * Description: check and return basic stats for a fastq file
 * Exported functions:
 * HISTORY:
 * Last edited: Oct 26 13:58 2006 (rd)
 * Created: Tue May  9 01:05:21 2006 (rd)
 *-------------------------------------------------------------------
 */

#include <stdio.h>
#include <math.h>
#include "readseq.h"
#define _A 0
#define _C 1
#define _G 2
#define _T 3
#define _N 4

#define MAX_LENGTH 10000

int main (int argc, char **argv)
{
  int i, j, length, lengthMax = 0, qMax = 0 ;
  int nseq = 0, total = 0 ;
  char *seq, *id ;
  unsigned char *qval ;
  FILE *fil ;
  static int sum[5], qsum[256] ;		/* 0 automatically */
  static int psum[MAX_LENGTH][5], pqsum[MAX_LENGTH][256], nlen[MAX_LENGTH] ;
  double erate;

  if (argc == 1)
    fil = stdin ;
  else if (!(fil = fopen (argv[1], "r")))
    { fprintf (stderr, "Failed to open fastq file %s\n", argv[1]) ;
      exit (-1) ;
    }

  while (readFastq (fil, dna2indexConv, &seq, &qval, &id, &length))
    { ++nseq ; ++nlen[length] ;
      total += length ;
      if (length > lengthMax) 
	{ lengthMax = length ;
	  if (length > MAX_LENGTH)
	    { fprintf (stderr, "read %s length = %d longer than MAX_LENGTH = %d; edit and recompile with larger MAX_LENGTH\n", id, length, MAX_LENGTH) ;
	      exit (-1) ;
	    }
	}
      for (i = 0 ; i < length ; ++i)
	{ ++sum[seq[i]] ; ++psum[i][seq[i]] ;
	  ++qsum[qval[i]] ; ++pqsum[i][qval[i]] ;
	  if (qval[i] > qMax) qMax = qval[i] ;
	}
      free (seq) ; free (qval) ; free (id) ;
    }

  printf ("%d sequences, %d total length", nseq, total) ;
  if (nseq)
    printf (", %.2f average, %d max", total/(float)nseq, lengthMax) ;
  printf ("\n") ;

  if (total)
    { printf ("Standard deviations at 0.25:  total %5.2f %%, per base %5.2f %%\n", 
	      100*(sqrt(0.25*(double)total)/total), 100*(sqrt(0.25*(double)nseq)/nseq)) ;
      printf ("            A    C    G    T    N ") ;
      for (i = 0 ; i <= qMax ; ++i) printf ("%4d",i) ;
      printf (" AQ\nTotal  ") ;
      printf ("  %4.1f %4.1f %4.1f %4.1f %4.1f ", 
	      100*((double)sum[_A]/total), 100*((double)sum[_C]/total),
	      100*((double)sum[_G]/total), 100*((double)sum[_T]/total),
	      100*((double)sum[_N]/total)) ;
      for (erate = j = 0 ; j <= qMax ; ++j) {
	printf ("%4d",(int)lrint(1000*((double)qsum[j]/total))) ;
	erate += pow(10, j/-10.0) * qsum[j];
      }
      printf(" %4.1f", -10*log(erate/total)/log(10));
      for (i = 0 ; i < lengthMax ; ++i)
	{ nseq -= nlen[i] ;
	  printf ("\nbase %2d", i+1) ;
	  printf ("  %4.1f %4.1f %4.1f %4.1f %4.1f ",
		  100*((double)psum[i][_A]/nseq), 100*((double)psum[i][_C]/nseq),
		  100*((double)psum[i][_G]/nseq), 100*((double)psum[i][_T]/nseq),
		  100*((double)psum[i][_N]/nseq)) ;
	  for (erate = j = 0 ; j <= qMax ; ++j) {
	      printf ("%4d",(int)lrint(1000*((double)pqsum[i][j]/nseq))) ;
	      erate += pow(10, j/-10.0) * pqsum[i][j];
	  }
	  printf(" %4.1f", -10*log(erate/nseq)/log(10));
	}
      printf ("\n") ;
    }
  exit(0);

}
