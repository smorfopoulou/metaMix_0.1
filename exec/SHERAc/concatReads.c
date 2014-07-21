/*
 *  concatReads.c
 *  SHERAc
 *
 *  Created by Martin Pollard on Thursday, 17 March 2011.
 *  Copyright 2011 Martin Pollard. All rights reserved.
 *
 */


#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <limits.h>
#include "concatReads.h"
#include "svn_version.h"

int main(int argc, char** argv)
{
	// Parse command line inputs
	const struct option longopts[] = {
		{"debug", no_argument, NULL, 0},
		{"numseqs2process", required_argument, NULL, 0},
		{"numseqs2skip", required_argument, NULL, 0},
		{"outputFolder", required_argument, NULL, 0},
		{"baseFreqFile", required_argument, NULL, 0},
		{"binomialTableFile", required_argument, NULL, 0},
		{"adaptersFile", required_argument, NULL, 0},
		{"adaptorsFile", required_argument, NULL, 0},
		{"minOverlap", required_argument, NULL, 0},
		{"version", no_argument, NULL, 0},
		{"maskThreshold", required_argument, NULL, 0},
		{"matchThreshold", required_argument, NULL, 0},
		{"help", no_argument, NULL, 'h'},
		{NULL, 0, NULL, 0 }
	};
	char* forwardreadsfn = NULL;
	char* reversereadsfn = NULL;
	char* outputprefix = NULL;

	int option_index = 0;
	const char* optionstring = "h?";
	int option = -1;
	// extract --options
	while ((option = getopt_long(argc, argv, optionstring, longopts, &option_index)) != -1)
	{
		switch (option)
		{
		case '?':
		case 'h':
			printf(
					"SHERAc\r\n"
					"usage: concatReads [options] filea fileb resultset\r\n"
					"\r\n"
					"The following options are valid:\r\n"
					"\t--adaptersFile - sets which file contains the sequence adaptors\r\n"
					"\t--maskThreshold - quality level below which a base is masked\r\n"
					"\t--matchThreshold - score statistic threshold below which a match is rejected\r\n"
					"\t--minOverlap - the minimum number of bases which overlap to be considered a match\r\n"
					"\t--version - displays the version of SHERAc you're running\r\n"
					);
			return EXIT_SUCCESS;
		case 0:
			if (!strcmp(longopts[option_index].name,"adaptersFile") ||
					!strcmp(longopts[option_index].name,"adaptorsFile"))
			{
				adaptersFilename = strdup(optarg);
			}
			else if (!strcmp(longopts[option_index].name,"minOverlap"))
			{
				int overlap = atoi(optarg);
				if (0 < overlap &&  255 > overlap)
				{
					minOverlap = overlap;
				}
				else
				{
					const char* error = "Invalid minimum overlap\r\n";
					write( STDERR_FILENO, error, strlen(error));
					return EXIT_FAILURE;
				}
			}
			else if (!strcmp(longopts[option_index].name,"version"))
			{
				printf("SHERAc svn version:%s\r\n",svn_version());
				return EXIT_SUCCESS;
			}
			else if (!strcmp(longopts[option_index].name,"maskThreshold"))
			{
				int thresh = atoi(optarg);
				if (0 < thresh &&  255 > thresh)
				{
					QUAL_THRESHOLD = thresh;
				}
				else
				{
					const char* error = "Invalid mask quality threshold\r\n";
					write( STDERR_FILENO, error, strlen(error));
					return EXIT_FAILURE;
				}
			}
			else if (!strcmp(longopts[option_index].name,"matchThreshold"))
			{
				float thresh = atoi(optarg);
				if (0.0f < thresh && 1.0f > thresh)
				{
					THRESHOLD = thresh;
				}
				else
				{
					const char* error = "Invalid match quality threshold\r\n";
					write( STDERR_FILENO, error, strlen(error));
					return EXIT_FAILURE;
				}
			} 
			break;
		default:
			// Should never get here
			abort();
			break;
		}
	}
	
	// Set adapters filename to default if not set
	if (adaptersFilename == NULL)
	{
		adaptersFilename = strdup(defaultAdaptersFilename);
	}
	
	// Get file names these should be all that's left in argc and argv
	if ( optind + 2 <= argc )
	{
		forwardreadsfn = argv[optind];
		reversereadsfn = argv[optind+1];
		if ( optind + 3 <= argc )
		{
			outputprefix = argv[optind+2];
		}
		else {
			outputprefix = "output";
		}
	}
	else {
		const char* error = "Error in parsing command line arguments you must provide the forward and reverse file names\r\n";
		write( STDERR_FILENO, error, strlen(error));
		return EXIT_FAILURE;
	}
	
	char* summaryFilename = NULL;
	char* dataFilename = NULL;
	char* dataLeftFilename = NULL;
	char* dataRightFilename = NULL;

	asprintf(&summaryFilename, summaryFilenamePattern, outputprefix);
	asprintf(&dataFilename, dataFilenamePattern, outputprefix);
	asprintf(&dataLeftFilename, dataLeftFilenamePattern, outputprefix);
	asprintf(&dataRightFilename, dataRightFilenamePattern, outputprefix);
	
	// Open files
	FILE* summary = fopen(summaryFilename, "w");  // This file tells us what went in to our file and what came out
	fprintf(summary, "SHERAc svn version:%s\r\n\r\n",svn_version());
	fprintf(summary, "Forward Read File: %s\r\n",forwardreadsfn);
	fprintf(summary, "Reverse Read File: %s\r\n",reversereadsfn);

	// First fetch adaptors
	char* fwdadaptor;
	int fwdadaptors;
	char* revadaptor;
	int revadaptors;
	size_t fwdadaptors_raw, revadaptors_raw;
	FILE* adaptfile = fopen(adaptersFilename,"r");
	if (adaptfile == NULL)
	{
		const char* error = "Unable to open adapters file\r\n";
		write( STDERR_FILENO, error, strlen(error));
		return EXIT_FAILURE;
	}
	size_t dummy = 0;
	free(readLine(adaptfile, &dummy));
	fwdadaptor = readLine(adaptfile, &fwdadaptors_raw);
	free(readLine(adaptfile, &dummy));
	revadaptor = readLine(adaptfile, &revadaptors_raw);
	fclose(adaptfile);
	
	if ( fwdadaptors_raw < INT_MAX && revadaptors_raw < INT_MAX )
	{
		fwdadaptors = fwdadaptors_raw;
		revadaptors = revadaptors_raw;
	}
	else {
		const char* error = "Adaptor file too big\r\n";
		write( STDERR_FILENO, error, strlen(error));
		return EXIT_FAILURE;
	}

	// Remove \n's
	fwdadaptors--;
	revadaptors--;

	// Show what adaptors are
	fprintf(summary, "Adaptor File: %s\r\n", adaptersFilename );
	fprintf(summary, "Forward Adaptor: %.*s\r\n", fwdadaptors, fwdadaptor );
	fprintf(summary,  "Reverse Adaptor: %.*s\r\n", revadaptors, revadaptor );
	char* revfcadaptor = flipAndCompAdaptor( revadaptor, revadaptors );
	fprintf(summary,  "Rev FC Adaptor:  %.*s\r\n", revadaptors, revfcadaptor);
	
	// Files are quite large but accessed sequentially
	int fwdfd = open(forwardreadsfn, O_RDONLY);
	int revfd = open(reversereadsfn, O_RDONLY);
	if (fwdfd == -1)
	{
		perror("Cannot open forward reads file");
		return -1;
	}
	if (revfd == -1)
	{
		perror("Cannot open reverse reads file");
		return -1;
	}
	
	// mmap them
	struct stat fwdstat;
	struct stat revstat;
	char* fwddata = NULL;
	char* revdata = NULL;
	fstat(fwdfd, &fwdstat);
	fstat(revfd, &revstat);
	fwddata = mmap(0, fwdstat.st_size, PROT_READ, MAP_SHARED, fwdfd, 0);
	if (fwddata == MAP_FAILED)
	{
		perror("Forward Data mmap failed");
		return -1;
	}
	// Tell OS to readahead cause we're only really going forward
	madvise(fwddata, fwdstat.st_size,MADV_SEQUENTIAL);
	
	revdata = mmap(0, revstat.st_size, PROT_READ, MAP_SHARED, revfd, 0);
	if (revdata == MAP_FAILED)
	{
		perror("Reverse Data mmap failed");
		return -1;
	}
	// Tell OS to readahead cause we're only really going forward
	madvise(revdata, revstat.st_size,MADV_SEQUENTIAL);
	
	char* fwdtop = fwddata + (fwdstat.st_size - 1);
	char* revtop = revdata + (revstat.st_size - 1);
	
	// Now open output files
#ifdef FASTA
	FILE* outputdata = fopen("output.fa", "w");
	FILE* outputquala = fopen("output.quala", "w");
#else
	FILE* outputdata = fopen(dataFilename, "w");
#endif

	FILE* outputleft = fopen(dataLeftFilename, "w");
	FILE* outputright = fopen(dataRightFilename, "w");
	
	// Read data
	char* fiter = fwddata;
	char* riter = revdata;
	int sequencesprocessed = 0;
	int sequencesoverlapping = 0;
	int sequencesrejectednool = 0;
	int sequencesrejectedbadscore = 0;
	int sequencesrejectedtooshort = 0;
	int sequencesbad = 0;
	
	while (fiter < fwdtop)
	{
		// First read labels
		char* labelf = fwddata;
		char* labelr = revdata;
		int labelfs;
		int labelrs;
		fetchLine(&fiter, &labelf, &labelfs, fwdtop);
		fetchLine(&riter, &labelr, &labelrs, revtop);
		
		// Next read sequence
		char* fwdseq_raw;
		char* revseq_raw;
		int fwdseqs_raw;
		int revseqs_raw;
		fetchLine(&fiter, &fwdseq_raw, &fwdseqs_raw, fwdtop);
		fetchLine(&riter, &revseq_raw, &revseqs_raw, revtop);
		
		// Read qual label
		char* fwdql;
		char* revql;
		int fwdqls;
		int revqls;
		fetchLine(&fiter, &fwdql, &fwdqls, fwdtop);
		fetchLine(&riter, &revql, &revqls, revtop);
		
		// Now read quality
		char* fwdqual_raw;
		char* revqual_raw;
		int fwdquals;
		int revquals;
		fetchLine(&fiter, &fwdqual_raw, &fwdquals, fwdtop);
		fetchLine(&riter, &revqual_raw, &revquals, revtop);
		
		// Sanity check
		if ( fwdseqs_raw != fwdquals || revseqs_raw != revquals )
		{
			// Maybe log this
			sequencesbad++;
			// Echo the read to the fail bin
			fprintf(outputleft, "%.*s\n%.*s\n%.*s\n%.*s\n", labelfs, labelf, fwdseqs_raw, fwdseq_raw, fwdqls, fwdql, fwdquals, fwdqual_raw);
			fprintf(outputright, "%.*s\n%.*s\n%.*s\n%.*s\n", labelrs, labelr, revseqs_raw, revseq_raw, revqls, revql, revquals, revqual_raw);			
			continue;  // skip it
		}

		// Pre-process the sequence masking bad reads
		char* fwdseq;
		char* revseq;
		unsigned char* fwdqual;
		unsigned char* revqual;
		int fwdseqs = fwdseqs_raw;
		int revseqs = revseqs_raw;
		processMaskQual(fwdseq_raw, fwdqual_raw, &fwdseqs, &fwdseq, &fwdqual);
		processMaskQual(revseq_raw, revqual_raw, &revseqs, &revseq, &revqual);
		
		// flip and compliment the reverse sequence
		char* revseqfc;
		unsigned char* revqualfc;
		flipAndComp(revseq, revqual, revseqs, &revseqfc, &revqualfc);
		free(revseq);
		free(revqual);
		
		
		// Align reads
		if (fwdseqs < 10 || revseqs < 10)
		{
			sequencesrejectedtooshort++;
			// Echo the read to the fail bin
			fprintf(outputleft, "%.*s\n%.*s\n%.*s\n%.*s\n", labelfs, labelf, fwdseqs_raw, fwdseq_raw, fwdqls, fwdql, fwdquals, fwdqual_raw);
			fprintf(outputright, "%.*s\n%.*s\n%.*s\n%.*s\n", labelrs, labelr, revseqs_raw, revseq_raw, revqls, revql, revquals, revqual_raw);			
		}
		else
		{
			int overlapLen;
			int* scoreVector = scoreThis(fwdseq, revseqfc, fwdseqs, revseqs, fwdadaptor, revadaptor, fwdadaptors, revadaptors, &overlapLen);
			
			// Compute max score and element containing it
			int maxScore;
			int maxScoreIndex;

			calcMaxScore(scoreVector, overlapLen, &maxScore, &maxScoreIndex);
			
			// Begin output function
			float scorestat = calcStatistic(maxScore, maxScoreIndex, scoreVector, overlapLen);

			if ( scorestat >= THRESHOLD  && maxScoreIndex != -1)
			{
				sequencesoverlapping++;
#ifdef FASTA
				outputFastA(outputdata, outputquala, labelf, fwdseq, revseqfc, fwdqual, revqualfc, labelfs, fwdseqs, revseqs, scorestat, maxScore, maxScoreIndex);
#else
				outputFastX(outputdata, labelf, fwdseq, revseqfc, fwdqual, revqualfc, labelfs, fwdseqs, revseqs, scorestat, maxScore, maxScoreIndex);
#endif
			}
			else {
				if ( maxScoreIndex == -1 )
				{
					sequencesrejectednool++;
				}
				else {
					sequencesrejectedbadscore++;
				}
				// Echo the read to the fail bin
				fprintf(outputleft, "%.*s\n%.*s\n%.*s\n%.*s\n", labelfs, labelf, fwdseqs_raw, fwdseq_raw, fwdqls, fwdql, fwdquals, fwdqual_raw);
				fprintf(outputright, "%.*s\n%.*s\n%.*s\n%.*s\n", labelrs, labelr, revseqs_raw, revseq_raw, revqls, revql, revquals, revqual_raw);			
			}

			// End output function
			
			free(scoreVector);
		}
		free(fwdqual);
		free(revqualfc);
		free(revseqfc);
		free(fwdseq);

		sequencesprocessed++;
#ifdef SHERA_DEBUG
		if (sequencesprocessed == 100) break;
#endif
	}
	
	writeSummary(summary, sequencesprocessed,sequencesoverlapping, sequencesrejectednool, sequencesrejectedbadscore, sequencesbad, sequencesrejectedtooshort);
	
	// Cleanup
	free(summaryFilename);
	free(dataFilename);
	free(dataLeftFilename);
	free(dataRightFilename);	
	free(fwdadaptor);
	free(revadaptor);
	free(revfcadaptor);

	// Close files
	fclose(outputright);
	fclose(outputleft);
#ifdef FASTA
	fclose(outputquala);
#endif
	fclose(outputdata);
	fclose(summary);
	munmap(revdata, revstat.st_size);
	munmap(fwddata, fwdstat.st_size);
	close(revfd);
	close(fwdfd);
	
	return 0;
}

char* readLine(FILE* input, size_t* size)
{
#ifdef __linux__
	char* data = NULL;
	*size = getline(&data, size, input);

#else
	char* inbuf = fgetln(input, size);
	char* data = malloc(*size);
	memcpy(data, inbuf, *size);
#endif
	return data;
}

// Writes out summary of what we've done
void writeSummary(FILE* summary, const int sequencesprocessed,const int sequencesoverlapping, const int sequencesrejectednool,
				  const int sequencesrejectedbadscore, const int sequencesbad, const int sequencesrejectedtooshort)
{
	fprintf(summary, "Sequences Processed: %d\nSequences Not Processed Bad Size: %d\nSequences Overlapping: %d\nSequences Rejected No Overlap: %d\nSequences Rejected Bad Score: %d\nSequences Rejected too short: %d\n", sequencesprocessed, sequencesbad, sequencesoverlapping, sequencesrejectednool, sequencesrejectedbadscore, sequencesrejectedtooshort);
}

int fetchLine(char** iter, char** begin, int* size, const char* top)
{
	if (*iter >= top ) {return -1;}
	*begin = *iter;
	while ( *iter != top && **iter != '\n' ) { (*iter)++; } 
	*size = *iter-*begin;
	// Now chomp the end
	while (*size != 0 && ((*begin)[*size] == '\n' || (*begin)[*size] =='\r' || (*begin)[*size] ==' ')) { (*size)--; }
	(*iter)++;
	(*size)++;
	return 0;
}

void flipAndComp(const char* inputseq, const unsigned char* inputqual, const int seqsize, char** outputseq,
				  unsigned char** outputqual)
{
	char* seq = (char*)malloc(seqsize);
	unsigned char* qual = (unsigned char*)malloc(seqsize);

	for (int i = 0; i < seqsize; i++)
	{
		switch (inputseq[i]) {
			case 'A':
				seq[seqsize-i-1] = 'T';
				break;
			case 'T':
				seq[seqsize-i-1] = 'A';
				break;
			case 'G':
				seq[seqsize-i-1] = 'C';
				break;
			case 'C':
				seq[seqsize-i-1] = 'G';
				break;
			case 'N':
				seq[seqsize-i-1] = 'N';
				break;
			default:
				seq[seqsize-i-1] = 'N';
				break;
		}
		qual[seqsize-i-1] = inputqual[i];
		
	}
	*outputseq = seq;
	*outputqual = qual;
}

void processMaskQual(const char* inputseq, const char* inputqual, int* seqsize, char** outputseq,
					 unsigned char** outputqual)
{
	// Trim the end of 'N' and bases with quality below threshold
	int iter = *seqsize - 1;
	while ( iter > 0 && (inputseq[ iter ] == 'N' || (inputqual[iter] - QUAL_ADJUST) < QUAL_THRESHOLD)) { iter--; }
	*seqsize = iter + 1;
	
	char* seq = (char*)malloc(*seqsize);
	unsigned char* qual = (unsigned char*)malloc(*seqsize);

	for (int i = 0; i < *seqsize; i++)
	{
		qual[i] = inputqual[i] - QUAL_ADJUST;

		if (qual[i] >= QUAL_THRESHOLD)
		{
			seq[i] = inputseq[i];
		}
		else {
			seq[i] = 'N';
		}
	}
	*outputseq = seq;
	*outputqual = qual;
}

// This could probably be done slightly more cleanly
float calcStatistic(const int bestscore, const int bestscoreIdx, const int* scores, const int scoresSize)
{
	// extract worst score from entire score array
	int worstScore = bestscore;
	for (int i = 0; i < scoresSize; i++)
	{
		worstScore = scores[i] < worstScore ? scores[i] : worstScore;
	}
	
	// create localscore array of 7 elements
	int localscore[7];
	// foreach element in scores from bestOverlap-3 to bestOverlap+3 (7 elements)
	int iter = 0;
	for (int j = bestscoreIdx-3; j <= bestscoreIdx+3; j++)
	{
		if ( j > 0 && j < scoresSize)
		{
			// if there is an element there put in in array
			localscore[iter] = scores[j];
		}
		else {
			// if not put worstscore in there
			localscore[iter] = worstScore;
		}
		iter++;
	}
	// Already have best score cached
	// So move it out the way
	localscore[3] = localscore[6];
	
	// Extract second best score and local worst score
	int secondbest = worstScore;
	int localworst = bestscore;
	for (int k = 0; k < 6; k++)
	{
		secondbest = localscore[k] > secondbest ? localscore[k] : secondbest;
		localworst = localscore[k] < localworst ? localscore[k] : localworst;
	}

	// take the absolute value of the best score - the second best score / best score-the worst score (the +0.00001 is to prevent divided by zeros
	return fabs( (bestscore - secondbest)/(bestscore-localworst+0.00001f));
}

int scoreMatch(const char fwd, const char rev)
{
	if (fwd == 'N' || rev == 'N')
		return 0;
	if (fwd == rev)
	{
		return 1;
	}
	return -1;
}

// Computes the max score and element containing it
void calcMaxScore(int* scoreVector, const int overlapLen, int* maxScore, int* maxScoreIndex)
{
	*maxScore = 0;
	*maxScoreIndex = -1;
	for (int scorei = 0; scorei < overlapLen; scorei++)
	{
		if (scoreVector[scorei] > *maxScore)
		{
			*maxScore = scoreVector[scorei];
			*maxScoreIndex = scorei;
		}
	}	
}

char* flipAndCompAdaptor(const char* input, const int size)
{
	char* revseqfc = (char*)malloc(size);
	for (int i = 0; i < size; i++)
	{
		switch (input[i]) {
			case 'A':
				revseqfc[size-i-1] = 'T';
				break;
			case 'T':
				revseqfc[size-i-1] = 'A';
				break;
			case 'G':
				revseqfc[size-i-1] = 'C';
				break;
			case 'C':
				revseqfc[size-i-1] = 'G';
				break;
			case 'N':
				revseqfc[size-i-1] = 'N';
				break;
			default:
				revseqfc[size-i-1] = 'N';
				break;
		}
		
	}
	return revseqfc;
}

int* scoreThis(const char* fwdseq, const char* revseq, const int fwdseqs, const int revseqs, const char* fwdadaptor,
				const char* revadaptor, const int fwdadaptors, const int revadaptors, int* overlapLen)
{
	int maxOverlap = fwdseqs < revseqs? fwdseqs:revseqs;
	int maxAdaptorOverlap = (fwdseqs+fwdadaptors) < (revseqs + revadaptors) ? (fwdseqs+fwdadaptors) : (revseqs + revadaptors);
	maxAdaptorOverlap = maxOverlap * 2 < maxAdaptorOverlap ? maxOverlap * 2 : maxAdaptorOverlap;
	*overlapLen = maxAdaptorOverlap-minOverlap + 1;
	int* scoreVector = malloc((*overlapLen)*sizeof(int));
	
	// Check standard alignments
	int overlap = minOverlap;
	
	for (; overlap<=maxOverlap; overlap++)
	{
		int score = 0;
		
		for (int j=0; j<overlap; j++)
		{
			score += scoreMatch(fwdseq[fwdseqs-overlap+j], revseq[j]);
		}
		scoreVector[overlap-minOverlap] = score;
	}
	
	// Check adaptor alignments
	for (; overlap <= maxAdaptorOverlap; overlap++)
	{
		int score = 0;
		int lap = 0;
		// Calcuate where to start comparing the 
		// (overlap -fwdseqs) gives us how much of the adapter we need to compare
		int fwdadaptorBase = fwdadaptors- (overlap - fwdseqs);
		
		// fwdadaptor <-> revseq
		for (; (fwdadaptorBase+lap) < fwdadaptors && lap < revseqs; lap++)
		{
			score += scoreMatch(fwdadaptor[fwdadaptorBase+lap], revseq[lap]);
		}
		
		int fwdadj = lap;
		// fwdseq <-> revseq
		for (; lap< overlap && lap < revseqs; lap++)
		{
			score += scoreMatch(fwdseq[lap-fwdadj], revseq[lap]);
		}
		
		int revadj = lap;
		// fwdseq <-> revadaptor
		for (; lap< overlap; lap++)
		{
			score += scoreMatch(fwdseq[lap-fwdadj], revadaptor[lap-revadj]);
		}
		scoreVector[overlap-minOverlap] = score;
	}
	return scoreVector;
}

void outputFastX(FILE* outputdata, const char* fwdlabel,
				 const char* fwdseq, const char* revseqfc,
				 const unsigned char* fwdqual, const unsigned char* revqualfc,
				 const int fwdlabels, const int fwdseqs, const int revseqs,
				 const float scorestat, const int maxScore, const int maxScoreIndex)
{
	// write modified label
	// append overlap length 'bp', bestScore, statistic
	fprintf(outputdata, "@%.*s_%dbp_%d.0_%.2f\n", (int)(fwdlabels-3)>0?(int)(fwdlabels-3):0, (fwdlabel+1), minOverlap+maxScoreIndex, maxScore, scorestat);
	
	int maxOverlap = fwdseqs < revseqs? fwdseqs:revseqs;
	int trueoverlap = (maxScoreIndex+minOverlap);
	int fwdtowrite = fwdseqs-trueoverlap;
	int revtowrite = revseqs-trueoverlap;
	
	int fwdwriteadj = fwdtowrite >= 0 ? fwdtowrite : 0;
	
	int outputquala_datas = fwdtowrite+revtowrite+trueoverlap;
	char* outputquala_data = malloc(outputquala_datas);
	
	if (fwdtowrite > 0)
	{
		fprintf(outputdata, "%.*s", fwdtowrite, fwdseq);
		memcpy(outputquala_data, fwdqual, fwdtowrite);
		for (int i = 0; i < fwdtowrite; i++)
		{
			outputquala_data[i] += QUAL_ADJUST;
		}					
	}
	
	if (trueoverlap <= maxOverlap)
	{
		for (int i = 0; i < trueoverlap; i++)
		{
			if (fwdseq[fwdseqs-trueoverlap+i]==revseqfc[i])
			{
				fprintf(outputdata,"%c",revseqfc[i]);
				unsigned char qual = fwdqual[fwdseqs-trueoverlap+i]+revqualfc[i];
				qual = qual > QUAL_MAX ? QUAL_MAX : qual;
				outputquala_data[fwdwriteadj+i] = qual+QUAL_ADJUST;
			}
			else {
				fprintf(outputdata,"N");
				outputquala_data[fwdwriteadj+i] = QUAL_ADJUST;
			}
		}
	}
	else
	{
		int fwdadj = trueoverlap-revseqs;
		int revadj = trueoverlap - fwdseqs;
		
		if ( revadj > revseqs || fwdadj > fwdseqs)
		{
		}
		else
		{
			for (int j = 0; j < revseqs-revadj && j < fwdseqs; j++)
			{
				if (fwdseq[j]==revseqfc[j+revadj])
				{
					fprintf(outputdata,"%c",fwdseq[j]);
					int qual = fwdqual[j]+revqualfc[j+revadj];
					qual = qual > QUAL_MAX ? QUAL_MAX : qual;
					outputquala_data[fwdwriteadj+j] = qual+QUAL_ADJUST;
				}
				else {
					fprintf(outputdata,"N");
					outputquala_data[fwdwriteadj+j] = QUAL_ADJUST;
				}
			}
		}
	}
	
	
	if (revtowrite > 0)
	{
		fprintf(outputdata, "%.*s", revtowrite, revseqfc+trueoverlap);
		memcpy(outputquala_data+trueoverlap+fwdtowrite,revqualfc+trueoverlap, revtowrite);
		for (int i = 0; i < revtowrite; i++)
		{
			(outputquala_data+trueoverlap+fwdtowrite)[i] += QUAL_ADJUST;
		}
	}
	fprintf(outputdata, "\n+\n%.*s\n",outputquala_datas,outputquala_data);
	free(outputquala_data);
}

void outputFastA(FILE* outputdata, FILE* outputquala, const char* fwdlabel,
				 const char* fwdseq, const char* revseqfc,
				 const unsigned char* fwdqual, const unsigned char* revqualfc,
				 const int fwdlabels, const int fwdseqs, const int revseqs,
				 const float scorestat, const int maxScore, const int maxScoreIndex)
{
	// write modified label to both files
	// append overlap length 'bp', bestScore, statistic
	fprintf(outputdata, ">%.*s_%dbp_%d.0_%.2f\n", (int)(fwdlabels-3)>0?(int)(fwdlabels-3):0, (fwdlabel+1), minOverlap+maxScoreIndex, maxScore, scorestat);
	fprintf(outputquala, ">%.*s_%dbp_%d.0_%.2f\n", (int)(fwdlabels-3)>0?(int)(fwdlabels-3):0, (fwdlabel+1), minOverlap+maxScoreIndex, maxScore, scorestat);
	
	int maxOverlap = fwdseqs < revseqs? fwdseqs:revseqs;
	int trueoverlap = (maxScoreIndex+minOverlap);
	int fwdtowrite = fwdseqs-trueoverlap;
	int revtowrite = revseqs-trueoverlap;
	
	if (fwdtowrite > 0)
	{
		fprintf(outputdata, "%.*s", fwdtowrite, fwdseq);
		for (int fwdqualiter = 0; fwdqualiter < fwdtowrite; fwdqualiter++)
		{
			fprintf(outputquala,"%u ", fwdqual[fwdqualiter]);
		}
	}
	
	if (trueoverlap <= maxOverlap)
	{
		for (int i = 0; i < trueoverlap; i++)
		{
			if (fwdseq[fwdseqs-trueoverlap+i]==revseqfc[i])
			{
				fprintf(outputdata,"%c",revseqfc[i]);
				int qual = fwdqual[fwdseqs-trueoverlap+i]+revqualfc[i];
				qual = qual > QUAL_MAX ? QUAL_MAX : qual;
				fprintf(outputquala,"%u ", qual);
			}
			else {
				fprintf(outputdata,"N");
				fprintf(outputquala, "0 ");
			}
		}
	}
	else
	{
		int fwdadj = trueoverlap-revseqs;
		int revadj = trueoverlap - fwdseqs;
		
		if ( revadj > revseqs || fwdadj > fwdseqs)
		{
		}
		else
		{
			for (int j = 0; j < revseqs-revadj && j < fwdseqs; j++)
			{
				if (fwdseq[j]==revseqfc[j+revadj])
				{
					fprintf(outputdata,"%c",fwdseq[j]);
					int qual = fwdqual[j]+revqualfc[j+revadj];
					qual = qual > QUAL_MAX ? QUAL_MAX : qual;
					fprintf(outputquala,"%u ", qual);
					
				}
				else {
					fprintf(outputdata,"N");
					fprintf(outputquala, "0 ");
				}
			}
		}
	}
	
	
	if (revtowrite > 0)
	{
		fprintf(outputdata, "%.*s", revtowrite, revseqfc+trueoverlap);
		
		for (int revqualiter = 0; revqualiter < revtowrite; revqualiter++)
		{
			fprintf(outputquala,"%u ", revqualfc[trueoverlap+revqualiter]);
		}
		
	}
	fprintf(outputdata, "\n");
	fprintf(outputquala, "\n");	
}
