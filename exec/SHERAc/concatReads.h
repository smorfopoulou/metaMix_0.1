/*
 *  concatReads.h
 *  SHERAc
 *
 *  Created by Martin Pollard on 14/04/2011.
 *  Copyright 2011 Martin Pollard. All rights reserved.
 *
 */

// Forward declarations
char* readLine(FILE* input, size_t* size);
int fetchLine(char** iter, char** begin, int* size, const char* top);
void processMaskQual(const char* inputseq, const char* inputqual, int* seqsize, char** outputseq,
					 unsigned char** outputqual);
void flipAndComp(const char* inputseq, const unsigned char* inputqual, const int seqsize, char** outputseq,
				 unsigned char** outputqual);
inline int scoreMatch(const char fwd, const char rev);
float calcStatistic(const int bestscore, const int bestscoreIdx, const int* scores, const int scoresSize);
void calcMaxScore(int* scoreVector, const int overlapLen, int* maxScore, int* maxScoreIndex);
void writeSummary(FILE* summary, const int sequencesprocessed, const int sequencesoverlapping, const int sequencesrejectednool,
				  const int sequencesrejectedbadscore, const int sequencesbad, const int sequencesrejectedtooshort);
char* flipAndCompAdaptor(const char* input, const int size);
int* scoreThis(const char* fwdseq, const char* revseqfc, const int fwdseqs, const int revseqs, const char* fwdadaptor,
			   const char* revadaptor, const int fwdadaptors, const int revadaptors, int* overlapLen);
void outputFastX(FILE* outputdata, const char* fwdlabel,
				 const char* fwdseq, const char* revseqfc,
				 const unsigned char* fwdqual, const unsigned char* revqualfc,
				 const int fwdlabels, const int fwdseqs, const int revseqs,
				 const float scorestat, const int maxScore, const int maxScoreIndex);
void outputFastA(FILE* outputdata, FILE* outputquala, const char* fwdlabel,
				 const char* fwdseq, const char* revseqfc,
				 const unsigned char* fwdqual, const unsigned char* revqualfc,
				 const int fwdlabels, const int fwdseqs, const int revseqs,
				 const float scorestat, const int maxScore, const int maxScoreIndex);


// Constants
const char* summaryFilenamePattern = "%s.summary";
const char* dataFilenamePattern = "%s.fq";
const char* dataLeftFilenamePattern = "%s_single_1.fq";
const char* dataRightFilenamePattern = "%s_single_2.fq";
const int QUAL_ADJUST = 33;

// Parameters
const char* defaultAdaptersFilename = "support/adapters.fa";
char* adaptersFilename = NULL;
float THRESHOLD = 0.50f;
int QUAL_THRESHOLD = 10;
int QUAL_MAX = 40;
int minOverlap = 10;
