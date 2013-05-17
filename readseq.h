/*  File: readseq.h
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last code change: Sep 24 21:23 2008 (rd)
 * Created: some time in 1993
 *-------------------------------------------------------------------
 * Copyright (c) 1993, 2008 Genome Research Limited.
 *
 * License:
 * This file is part of fastqcheck.
 *
 * fastqcheck is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see
 * L<http://www.gnu.org/licenses/>.
 *-------------------------------------------------------------------
 */

extern int readSequence (FILE *fil, int *conv,
			 char **seq, char **id, char **desc, int *length) ;
				/* read next sequence from file */
extern int readFastq (FILE *fil, int *conv, 
		      char **seq, unsigned char **qval, char **id, int *length) ;
				/* read next fastq entry from file */
extern int writeSequence (FILE *fil, int *conv, 
			  char *seq, char *id, char *desc, int len) ;
				/* write sequence to file, using convert */
extern int writeFastq (FILE *fil, int *conv,
		       char *seq, char *qval, char *id, int len) ;
				/* write fastq to file, using convert */
extern int seqConvert (char *seq, int *length, int *conv) ;
				/* convert in place - can shorten */
extern int readMatrix (char *name, int *conv, int** *mat) ;

extern int dna2textConv[] ;
extern int dna2textAmbig2NConv[] ;
extern int dna2indexConv[] ;
extern int dna2binaryConv[] ;
static const char index2char[] = "acgtn" ;
static const char binary2char[] = "-ACMGRSVTWYHKDBN" ;
extern int aa2textConv[] ;
extern int aa2indexConv[] ;
static const char index2aa[] = "ACDEFGHIKLMNPQRSTVWYX*" ;
extern int noConv[] ;

/***** end of file *****/
