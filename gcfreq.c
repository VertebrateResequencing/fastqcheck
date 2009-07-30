/*  File: gcfreq.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) Genome Research Limited, 2007
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Mar 10 02:19 2008 (rd)
 * Created: Wed Jan 31 23:08:40 2007 (rd)
 *-------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "readseq.h"

static int KMER = 35 ;
static int isHuman = 0 ;
static int isFastq = 0 ;
static int *refHist = 0 ;

static int human35[] = { 50000, 40666, 75338, 151950, 319678, 643075, 1186889,
1981812, 3010924, 4185652, 5349327, 6336269, 6976056, 7164402, 6925589, 6376732,
5681932, 5071200, 4602808, 4222855, 3709856, 2931026, 2039978, 1281887, 755779,
426143, 234189, 130134, 76483, 47335, 27610, 14042, 6268, 2442, 873, 238 } ;

static int human70[] = {3174, 2481, 3543, 4155, 5025, 6154, 8414, 11280, 17169, 27196, 
44990, 75012, 121075, 190836, 288382, 417859, 577470, 768482, 974004, 1187865, 1399913, 1594720, 
1766654, 1908469, 2011838, 2073520, 2089823, 2065318, 2012702, 1931333, 1830482, 1717664, 
1596210, 1466033, 1331400, 1212511, 1114904, 1042194, 992070, 938138, 864086, 755260, 619107, 
484930, 368407, 278431, 208731, 154504, 114709, 82820, 60186, 43866, 33664, 25736, 19957, 16018, 
12930, 10889, 8621, 6130, 4284, 2630, 1548, 847, 474, 270, 116, 60, 34, 13, 15 } ;

void report (char *prefix, int *hist, int *refHist) ;

int main (int argc, char *argv[])
{ 
  char *seq, *s ;
  int length, i ;
  int hist[256] ;
  char *prefix ;
  FILE *fil = 0 ;

  --argc ; ++argv ;

  while (argc && **argv == '-')
    if (!strcmp (*argv, "-k") && argc >= 2)
      { KMER = atoi (argv[1]) ;
	if (KMER < 1 || KMER >= 256)
	  { fprintf (stderr, "ERROR: k must be between 1 and 255\n") ; exit (-1) ; }
	argc -= 2 ; argv += 2 ;
      }
    else if (!strcmp (*argv, "-H"))
      { isHuman = 1 ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-q"))
      { isFastq = 1 ; --argc ; ++argv ; }
    else
      { fprintf (stderr, "unknown option %s - usage: gcfreq [-q] [-k n] [-H] [filename]\n", 
		 *argv) ; 
	exit (-1) ; 
      }

  if (isHuman)
    { if (KMER == 35)
	refHist = human35 ;
      else if (KMER == 70)
	refHist = human70 ;
      else
	{ fprintf (stderr, "when -H, k = %d must be 35 or 70\n", KMER) ; exit (-1) ; }
    }

  if (refHist)
    { int refTot = 0 ;
      double f, sum = 0 ;
      for (i = 0 ; i <= KMER ; ++i) refTot += refHist[i] ;
      f = 1.0 / refTot ;
      printf ("human %10d", refTot) ; 
      for (i = 0 ; i <= KMER ; ++i) printf (" %8.3g", f*refHist[i]) ;
      printf ("\n") ;
      printf ("human %10d", refTot) ; 
      for (i = 0 ; i <= KMER ; ++i) printf (" %8.3g", sum += f*refHist[i]) ;
      printf ("\n") ;
      fflush (stdout) ;
    }

/* now read files */

  if (!argc)
      { fprintf (stderr, "missing filenames - usage: gcfreq [-q] [-k n] [-H] [filename]+\n") ; 
	exit (-1) ; 
      }
  else
    while (argc)
      { 
	if (!prefix || strncmp (prefix, *argv, strlen(prefix)))
	  { if (prefix)
	      report (prefix, hist, refHist) ;
	    { char *cp ;
	      prefix = strdup (*argv) ; 
	      for (cp = prefix ; *cp && *cp >= '0' && *cp <= '9'; ++cp) ;
	      *cp = 0 ;
	    }
	    for (i = 0 ; i <= KMER ; ++i) hist[i] = 0 ;
	}
	
	if (!(fil = fopen (*argv, "r")))
	  { fprintf (stderr, "failed to open file %s\n", *argv) ; exit (-1) ; }
	while (isFastq ? readFastq (fil, dna2textConv, &seq, 0, 0, &length)
	               : readSequence (fil, dna2textConv, &seq, 0, 0, &length))
	  { s = seq ;
	    while (length >= KMER)
	      { int nGC = 0 ;
		for (i = KMER ; i-- ; length--, s++)
		  switch (*s)
		    {
		    case 'C': case 'G': ++nGC ; break ;
		    case 'A': case 'T': break ;
		    default: --length ; ++s ; goto startAgain ;
		    }
		++hist[nGC] ;
	      startAgain: continue ;
	      }
	    free (seq) ;
	  }

	fprintf (stderr, "done file %s\n", *argv) ;
	fclose (fil) ;
	--argc ; ++argv ;
      }

/* and report */

  report (prefix, hist, refHist) ;
}

void report (char *prefix, int *hist, int *refHist)
{
  int i, tot = 0, refTot = 0 ; 
  double f ;

  printf ("%5s", prefix) ;

  for (i = 0 ; i <= KMER ; ++i) 
    { tot +=  hist[i] ; 
      if (refHist) refTot += refHist[i] ; 
    }
  printf (" %10d", tot) ;
  if (!tot)
    return ;

  if (refHist)
    { f = refTot / (double) tot ;
      for (i = 0 ; i <= KMER ; ++i)
	printf (" %8.3g", (f*hist[i])/refHist[i]) ;
      printf ("\n") ;
    }
  else
    { f = 1.0 / (double) tot ;
      for (i = 0 ; i <= KMER ; ++i)
	printf (" %8.3g", f*hist[i]) ;
      printf ("\n") ;
    }

  printf ("\n") ;
  fflush (stdout) ;
}

/**************** end of file *****************/
