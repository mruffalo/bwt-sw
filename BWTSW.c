/*

   BWTSW.c		BWT-SW - Local Alignment by BWT-index (100% sensitivity)

   Copyright (C) 2006, Wong Chi Kwong.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "iniparser.h"
#include "MiscUtilities.h"
#include "MemManager.h"
#include "Timing.h"
#include "TextConverter.h"
#include "BWT.h"
#include "HSP.h"
#include "HSPstatistic.h"
#include "blast_dust.h"
#include "Socket.h"

// Query input for communication with server process
typedef struct QueryInput {
	char DatabaseName[MAX_FILENAME_LEN+1];
	char QueryFileName[MAX_FILENAME_LEN+1];
	char OutputFileName[MAX_FILENAME_LEN+1];
	char AlignFileName[MAX_FILENAME_LEN+1];
	char TimingFileName[MAX_FILENAME_LEN+1];
	int OutputFormat;
	int QueryStrand;
	int MaskLowerCase;
	int DustMasking;
	int DustLevel;
	int DustWindow;
	int match;
	int mismatch;
	int gapOpen;
	int gapExt;
	double expectationValue;
} QueryInput;

#define SERVER_SOCKET_NAME	"BWTSW_SOCKET"


// Database and ini
dictionary *ParseInput(int argc, char** argv);
void ParseIniFile(char *iniFileName);
void ProcessIni();
void ValidateIni();
void PrintIni();
void PrintShortDesc();
void PrintHelp();
void ProcessFileName(char *outputFileName, const char *inputFileName, const char *databaseLocation, const char *databaseName);

// Query parameters
void ParseQueryParameterFile(char *queryParameterFileName);
void ValidateQueryParameters();
void PrintQueryParameters();

// QSort
int FinalHitAllContextDbSeqIndexScorePosTextOrder(const void *gappedHitList, const int index1, const int index2);
int HitListPosTextDescOrder(const void *hitList, const int index1, const int index2);

#define BWTSW_MAX_QUERY_LENGTH		268435456	// 256M
#define SMALL_QUERY_ALLOCATION_UNIT	65536
#define LARGE_QUERY_ALLOCATION_UNIT	1048576

#define NON_INPUT_VALUE	99999999

	// Parameters
	int Confirmation = FALSE;
	int LoadIndex = FALSE;
	int UnloadIndex = FALSE;
    	
	// Memory parameters
	int PoolSize = 2097152;				// 2M  - fixed; not configurable through ini
	int WorkingMemorySize = 67108864;	// 64M - good for 8M hit; configurable through ini
	int AlignmentMemorySize = 4194304;	// 4M

	// Performance Statistics parameters
	int PrintProgress = FALSE;
	int PrintPerformanceStatistics = FALSE;
	int PrintHistogram = FALSE;
	double HistogramMinEvalue = 1e-24;
	int PrintDPStatistics = FALSE;
	int PrintProgressDepth = 0;

	// Database parameters
	char DatabaseName[MAX_FILENAME_LEN+1] = "";
	char DatabaseLocation[MAX_FILENAME_LEN+1] = "./";

	// Query Parameter file
	int AlternateQueryParameter = FALSE;
	char QueryParameterFileName[MAX_FILENAME_LEN+1] = "";

	// Input and output file
	char QueryFileName[MAX_FILENAME_LEN+1] = "";
	char OutputFileName[MAX_FILENAME_LEN+1] = "";
	char AlignFileName[MAX_FILENAME_LEN+1] = "";
	char TimingFileName[MAX_FILENAME_LEN+1] = "";

	// Output parameters
	int OutputFormat = OUTPUT_PAIRWISE;

	// DatabaseFiles parameters
	char AnnotationFileName[MAX_FILENAME_LEN+1] = "*.ann";
	char AmbiguityFileName[MAX_FILENAME_LEN+1] = "*.amb";
	char PackedDNAFileName[MAX_FILENAME_LEN+1] = "*.pac";
	char BWTCodeFileName[MAX_FILENAME_LEN+1] = "*.bwt";
	char BWTOccValueFileName[MAX_FILENAME_LEN+1] = "*.fmv";
	char SaValueFileName[MAX_FILENAME_LEN+1] = "*.sa";

	// Strand parameters
	int QueryStrand = QUERY_BOTH_STRAND;

	// Filtering parameters
	int MaskLowerCase = FALSE;
	int DustMasking = TRUE;
	int DustLevel = 20;
	int DustWindow = 64;
	char *MaskOption;

	// Extension Statistic parameters 
	extern STAT_scoreBlock stat_scoreBlock;
	extern double stat_expectationValue;
	double InputEValue = (double)NON_INPUT_VALUE;
	int InputMatch = NON_INPUT_VALUE;
	int InputMismatch = NON_INPUT_VALUE;
	int InputGapOpen = NON_INPUT_VALUE;
	int InputGapExtend = NON_INPUT_VALUE;

int main(int argc, char** argv) {

	// Program input
	dictionary *programInput;

	// Context
	Context context[2];
	ContextInfo contextInfo[2];
	int numOfContext;
	int currentContext;

	// Query patterns
	FILE *queryFile;
	int numOfQueryPattern, queryPatternNameLength, queryPatternLength, numUnmaskedChar;
	char queryPatternName[MAX_SEQ_NAME_LENGTH + 1];
	unsigned char *queryPattern;
	int queryPatternAllocated = 0;
	char filename[MAX_FILENAME_LEN];

	// Character frequency and probabilities
	unsigned int queryCharFreq[16];
	double queryCharProb[16];
	unsigned int dbCharFreq[16];
	double dbCharProb[16];
	int scoringMatrix[16][16];

	// Output file
	FILE *outputFile;
	FILE *alignFile;
	FILE *timingFile;

	// Character map
	unsigned char charMap[256];
	unsigned char complementMap[256];

	// Index
	BWT *bwt;
	HSP *hsp;

	// Communicating to server process
	Socket *bwtServerSocket;
	QueryInput QueryInput;

	// Working variables
	char queryChar, c;
	int i, j, k, n;
	MMPool *mmPool;
	MMPool *workingMemoryPool;
	MMPool *alignmentPool;

	char* workingMemory;
	char* workingMemoryPointer;
	unsigned int workingMemoryUsed;
	HitList *hitList;
	GappedHitListWithAlignment* __restrict finalHitListAllContext;
	int *dbScore;
	int *dbOrder;
	int *dbRank;

	int depthBit;
	unsigned int depthMask;

	BWTSaRetrievalStatistics bwtSaRetrievalStatistics;
	BWTDPStatistics bwtDPStatistics;
	int numOfSaIndexGroup;
	SaIndexGroupNew *saIndexGroup;

	int numOfGappedHitWithTraceback;
	int numOfHitForQuery;

	int numOfQueryPos;
	unsigned int *queryPosList;
	int *queryPosListIndex;
	int numOfDPPoint;

	// Performance statistics
	double startTime;
	double elapsedTime = 0, totalQueryTime = 0;
	double databaseLoadTime, totalSearchTime, totalTextDecodeTime, totalUngappedExtensionTime, totalGappedExtensionTime1, totalGappedExtensionTime2, totalPrintTime;
	double queryStartTime;
	long long totalSaIndexRange, totalHitGenerated, totalUngappedHit, totalUniqueUngappedHit;
	long long totalGappedHit, totalNonCrossBoundaryGappedHit, totalUniqueGappedHit;
	long long totalGappedHitWithTraceback, totalFinalHit;
	long long totalNumOfUnmaskedChar, totalNumOfChar;
	long long totalNumOfDPPoint;

	// Count Histogram statistic by score
	Histogram *histogram;


	// Program input
	programInput = ParseInput(argc, argv);
	PrintShortDesc();

	if (UnloadIndex) {
		if (Confirmation) {
			printf("Unloading BWT-SW server process. Press Y to go or N to cancel. ");
			queryChar = (char)getchar();
			while (queryChar != 'y' && queryChar != 'Y' && queryChar != 'n' && queryChar!= 'N') {
				queryChar = (char)getchar();
			}
			if (queryChar == 'n' || queryChar == 'N') {
				exit(0);
			}
		}
		QueryInput.DatabaseName[0] = '\0';	// An exit message
		bwtServerSocket = SocketInitiateConnection(LOCAL_SOCKET, SERVER_SOCKET_NAME);
		if (bwtServerSocket == NULL) {
			fprintf(stderr, "Server process does not exist!\n");
			exit(1);
		}
		if (SocketSend(bwtServerSocket, &QueryInput, sizeof(QueryInput)) < 0) {
			fprintf(stderr, "BWTSW(): Socket communication error!\n");
			exit(1);
		}
		SocketEndConnection(bwtServerSocket);
		exit(0);
	}

	if (LoadIndex) {

		// Ini
		if (strcmp(argv[0] + strlen(argv[0]) - 4, ".exe") == 0) {
			*(argv[0] + strlen(argv[0]) - 4) = '\0';
		}
		sprintf(filename, "%s.ini", argv[0]);			// First default ini is <program name>.ini
		ParseIniFile(filename);
		sprintf(filename, "%s.ini", DatabaseName);		// Second default ini is <database name>.ini
		ParseIniFile(filename);
		printf("\n");
		ProcessIni();
		ValidateIni();
		PrintIni();

		if (Confirmation) {
			printf("Loading BWT-SW server process. Press Y to go or N to cancel. ");
			c = (char)getchar();
			while (c != 'y' && c != 'Y' && c != 'n' && c!= 'N') {
				c = (char)getchar();
			}
			if (c == 'n' || c == 'N') {
				exit(0);
			}
		}

		// Trying to communicate with server process
		bwtServerSocket = SocketInitiateConnection(LOCAL_SOCKET, SERVER_SOCKET_NAME);
		if (bwtServerSocket != NULL) {
			QueryInput.DatabaseName[0] = '\0';	// An exit message
			if (SocketSend(bwtServerSocket, &QueryInput, sizeof(QueryInput)) < 0) {
				fprintf(stderr, "BWTSW(): Socket communication error!\n");
				exit(1);
			}
			SocketEndConnection(bwtServerSocket);
			exit(0);

		}

	} else {

		// Query

		// Ini
		if (strcmp(argv[0] + strlen(argv[0]) - 4, ".exe") == 0) {
			*(argv[0] + strlen(argv[0]) - 4) = '\0';
		}
		sprintf(filename, "%s.ini", argv[0]);			// First default ini is <program name>.ini
		ParseIniFile(filename);
		sprintf(filename, "%s.ini", DatabaseName);		// Second default ini is <database name>.ini
		ParseIniFile(filename);
		printf("\n");
		ProcessIni();
		ValidateIni();
		PrintIni();

		// Query parameters
		sprintf(filename, "%s.parm", argv[0]);			// First default query parameter is <program name>.parm
		ParseQueryParameterFile(filename);
		sprintf(filename, "%s.parm", DatabaseName);		// Second default query parameter is <database name>.parm
		ParseQueryParameterFile(filename);
		if (AlternateQueryParameter) {
			ParseQueryParameterFile(QueryParameterFileName);	// User supplied alternate query parameter
		}
		printf("\n");
		ValidateQueryParameters();
		PrintQueryParameters();

		if (Confirmation == TRUE) {
			printf("Press Y to go or N to cancel. ");
			c = (char)getchar();
			while (c != 'y' && c != 'Y' && c != 'n' && c!= 'N') {
				c = (char)getchar();
			}
			if (c == 'n' || c == 'N') {
				exit(0);
			}
		}

	}

	if (!LoadIndex) {
		// Trying to communicate with server process
		bwtServerSocket = SocketInitiateConnection(LOCAL_SOCKET, SERVER_SOCKET_NAME);
		if (bwtServerSocket != NULL) {
			// Send request to server process
			memcpy(QueryInput.DatabaseName, DatabaseName, MAX_FILENAME_LEN+1);
			memcpy(QueryInput.QueryFileName, QueryFileName, MAX_FILENAME_LEN+1);
			memcpy(QueryInput.OutputFileName, OutputFileName, MAX_FILENAME_LEN+1);
			memcpy(QueryInput.AlignFileName, AlignFileName, MAX_FILENAME_LEN+1);
			memcpy(QueryInput.TimingFileName, TimingFileName, MAX_FILENAME_LEN+1);
			QueryInput.OutputFormat = OutputFormat;
			QueryInput.QueryStrand = QueryStrand;
			QueryInput.MaskLowerCase = MaskLowerCase;
			QueryInput.DustMasking = DustMasking;
			QueryInput.DustLevel = DustLevel;
			QueryInput.DustWindow = DustWindow;
			QueryInput.match = stat_scoreBlock.match;
			QueryInput.mismatch = stat_scoreBlock.mismatch;
			QueryInput.gapOpen = stat_scoreBlock.gapOpen;
			QueryInput.gapExt = stat_scoreBlock.gapExt;
			QueryInput.expectationValue = stat_expectationValue;

			if (SocketSend(bwtServerSocket, &QueryInput, sizeof(QueryInput)) < 0) {
				fprintf(stderr, "BWTSW(): Socket communication error!\n");
				exit(1);
			}

			// Redirect output from server
			SocketSetRedirect(bwtServerSocket);
			while ((n = SocketRedirect()) > 0) {
			}

			// End connection and exit
			SocketEndConnection(bwtServerSocket);

			if (n == 0) {
				exit(0);
			} else {
				exit(1);
			}

		}
	}

	// Measure searching performance
	startTime = setStartTime();

	// Initialize memory manager
	MMMasterInitialize(3, 0, FALSE, NULL);
	mmPool = MMPoolCreate(PoolSize);

	// Initialize character map
	HSPFillCharMap(charMap);
	HSPFillComplementMap(complementMap);
	HSPFillScoringMatrix(scoringMatrix, stat_scoreBlock.match, stat_scoreBlock.mismatch, 0);

	if (!LoadIndex) {

		// Check and open query file
		queryFile = (FILE*)fopen64(QueryFileName, "r");
		if (queryFile == NULL) {
			fprintf(stderr, "Cannot open query file %s!\n", QueryFileName);
			exit(1);
		}

		// Open output file
		outputFile = (FILE*)fopen64(OutputFileName, "w");
		if (outputFile == NULL) {
			fprintf(stderr, "Cannot open output file %s!\n", OutputFileName);
			exit(1);
		}

		// Open extra alignment file
		if (AlignFileName[0] != ' ' && AlignFileName[0] != '-' && AlignFileName[0] != '\0') {
			alignFile = fopen(AlignFileName, "w");
			if (alignFile == NULL) {
				fprintf(stderr, "Cannot open alignment file %s!\n", AlignFileName);
				exit(1);
			}
		} else {
			alignFile = NULL;
		}

		// Open timing file
		if (TimingFileName[0] != ' ' && TimingFileName[0] != '-' && TimingFileName[0] != '\0') {
			timingFile = fopen(TimingFileName, "a");
			if (timingFile == NULL) {
				fprintf(stderr, "Cannot open timing file %s!\n", TimingFileName);
				exit(1);
			}
		} else {
			timingFile = NULL;
		}

	}

	// Load Database
	printf("Loading Database..");
	fflush(stdout);

	bwt = BWTLoad(mmPool, BWTCodeFileName, BWTOccValueFileName, SaValueFileName, NULL, NULL, NULL);
	hsp = HSPLoad(mmPool, PackedDNAFileName, AnnotationFileName, AmbiguityFileName);
	if (bwt->textLength != hsp->dnaLength) {
		fprintf(stderr, "BWT-SW: Database length inconsistent!\n");
		exit(1);
	}

	// Initialize character frequencies and probabilities
	for (i=0; i<16; i++) {
		dbCharFreq[i] = 0;
		dbCharProb[i] = 0.0;
	}
	for (i=0; i<4; i++) {
		dbCharProb[i] = 0.25;
	}

	histogram = HSPAllocateHistogram(HistogramMinEvalue, stat_expectationValue);

	contextInfo[0].complemented = FALSE;
	contextInfo[0].reversed = FALSE;
	strcpy(contextInfo[0].name, "Plus");
	contextInfo[1].complemented = TRUE;
	contextInfo[1].reversed = TRUE;
	strcpy(contextInfo[1].name, "Minus");

	// Allocate DB sequence score array
	dbScore = MMUnitAllocate(hsp->numOfSeq * sizeof(int));
	dbOrder = MMUnitAllocate(hsp->numOfSeq * sizeof(int));
	dbRank = MMUnitAllocate(hsp->numOfSeq * sizeof(int));

	// Allocate working memory for hit generation
	workingMemoryPool = MMPoolCreate(WorkingMemorySize);
	workingMemory = (char*)(workingMemoryPool + 1);

	// Allocate memory for storing alignments
	alignmentPool = MMPoolCreate(AlignmentMemorySize);

	// Allocate memory for BWT statistics
	BWTAllocateDPStatistics(&bwtDPStatistics);

	// Allocate memory for query
	queryPattern = MMUnitAllocate(SMALL_QUERY_ALLOCATION_UNIT);
	queryPatternAllocated = SMALL_QUERY_ALLOCATION_UNIT;

	// Finished loading
	databaseLoadTime = getElapsedTime(startTime);

	printf("Finished.\n");
	printf("There are %d sequences in the database. Database size is %u.\n\n", hsp->numOfSeq, hsp->dnaLength);
	printf("Elapsed time = %9.4f s", databaseLoadTime);
	printf("\n");
	fflush(stdout);

	for (;;) {	// Loop until unload request is received if this is a server process
				// or exit if this is a client process

		if (LoadIndex) {

			databaseLoadTime = 0;
			
			// Accept connection
			bwtServerSocket = SocketCreate(LOCAL_SOCKET, SERVER_SOCKET_NAME);
			if (bwtServerSocket == NULL) {
				fprintf(stderr, "BWTSW(): Cannot create socket!\n");
				fprintf(stderr, "Server process is unloaded.\n");
				exit(1);
			}
			if (SocketAcceptConnection(bwtServerSocket) < 0) {
				fprintf(stderr, "BWTSW(): Error accepting connection!\n");
				fprintf(stderr, "Server process is unloaded.\n");
				exit(1);
			}
			if (SocketReceive(bwtServerSocket, &QueryInput, sizeof(QueryInput)) < 0) {
				fprintf(stderr, "BWTSW(): Error receiving request!\n");
				fprintf(stderr, "Server process is unloaded.\n");
				exit(1);
			}

			// Prepare to redirect output to socket
			SocketSetRedirect(bwtServerSocket);

			if (QueryInput.DatabaseName[0] == '\0') {
				// Unload message received
				printf("BWT-SW server process is unloaded.\n");
				SocketFree(bwtServerSocket);
				exit(0);
			}

			if (strncmp(QueryInput.DatabaseName, DatabaseName, MAX_FILENAME_LEN) != 0) {
				Socketfprintf(stderr, "Server is loaded with another database.\n");
				SocketEndConnection(bwtServerSocket);
				continue;
			}
			memcpy(QueryFileName, QueryInput.QueryFileName, MAX_FILENAME_LEN+1);
			memcpy(OutputFileName, QueryInput.OutputFileName, MAX_FILENAME_LEN+1);
			memcpy(AlignFileName, QueryInput.AlignFileName, MAX_FILENAME_LEN+1);
			memcpy(TimingFileName, QueryInput.TimingFileName, MAX_FILENAME_LEN+1);
			OutputFormat = QueryInput.OutputFormat;
			QueryStrand = QueryInput.QueryStrand;
			MaskLowerCase = QueryInput.MaskLowerCase;
			DustMasking = QueryInput.DustMasking;
			DustLevel = QueryInput.DustLevel;
			DustWindow = QueryInput.DustWindow;
			stat_scoreBlock.match = QueryInput.match;
			stat_scoreBlock.mismatch = QueryInput.mismatch;
			stat_scoreBlock.gapOpen = QueryInput.gapOpen;
			stat_scoreBlock.gapExt = QueryInput.gapExt;
			stat_expectationValue = QueryInput.expectationValue;

			// Check and open query file
			queryFile = (FILE*)fopen64(QueryFileName, "r");
			if (queryFile == NULL) {
				Socketfprintf(stderr, "Cannot open query file %s!\n", QueryFileName);
				SocketEndConnection(bwtServerSocket);
				continue;
			}

			// Open output file
			outputFile = (FILE*)fopen64(OutputFileName, "w");
			if (outputFile == NULL) {
				Socketfprintf(stderr, "Cannot open output file %s!\n", OutputFileName);
				SocketEndConnection(bwtServerSocket);
				continue;
			}

			// Open extra alignment file
			if (AlignFileName[0] != ' ' && AlignFileName[0] != '-' && AlignFileName[0] != '\0') {
				alignFile = fopen(AlignFileName, "w");
				if (alignFile == NULL) {
					Socketfprintf(stderr, "Cannot open alignment file %s!\n", AlignFileName);
					SocketEndConnection(bwtServerSocket);
					continue;
				}
			} else {
				alignFile = NULL;
			}

			// Open timing file
			if (TimingFileName[0] != ' ' && TimingFileName[0] != '-' && TimingFileName[0] != '\0') {
				timingFile = fopen(TimingFileName, "a");
				if (timingFile == NULL) {
					Socketfprintf(stderr, "Cannot open timing file %s!\n", TimingFileName);
					SocketEndConnection(bwtServerSocket);
					continue;
				}
			} else {
				timingFile = NULL;
			}

		}

		// Process query

		startTime = setStartTime();
		totalQueryTime = 0;

		HSPPrintHeader(outputFile, OutputFormat, DatabaseName, hsp->numOfSeq, hsp->dnaLength);
		if (alignFile != NULL) {
			HSPPrintHeader(alignFile, OUTPUT_PAIRWISE, DatabaseName, hsp->numOfSeq, hsp->dnaLength);
		}

		// Initialize statistics variables
		totalSearchTime = 0;
		totalTextDecodeTime = 0;
		totalUngappedExtensionTime = 0;
		totalGappedExtensionTime1 = 0;
		totalGappedExtensionTime2 = 0;
		totalPrintTime = 0;

		totalSaIndexRange = 0;
		totalHitGenerated = 0;
		totalUngappedHit = 0;
		totalUniqueUngappedHit = 0;
		totalGappedHit = 0;
		totalNonCrossBoundaryGappedHit = 0;
		totalUniqueGappedHit = 0;
		totalGappedHitWithTraceback = 0;
		totalFinalHit = 0;
		totalNumOfUnmaskedChar = 0;
		totalNumOfChar = 0;

		totalNumOfDPPoint = 0;

		HSPInitializeHistogram(histogram);


		numOfQueryPattern = 0;

		BWTInitializeSaRetrievalStatistics(&bwtSaRetrievalStatistics);
		BWTInitializeDPStatistics(&bwtDPStatistics);

		queryChar = (char)getc(queryFile);
		while (!feof(queryFile)) {

			// Parse input file in FASTA format

			while (!feof(queryFile) && queryChar != '>') {
				queryChar = (char)getc(queryFile);
			}
			if (feof(queryFile)) {
				break;
			}
			queryChar = (char)getc(queryFile);

			queryPatternNameLength = 0;
			while (!feof(queryFile) && queryPatternNameLength < MAX_SEQ_NAME_LENGTH && queryChar != '\n') {
				// Get query pattern name
				queryPatternName[queryPatternNameLength] = queryChar;
				queryPatternNameLength++;
				queryChar = (char)getc(queryFile);
			}
			queryPatternName[queryPatternNameLength] = '\0';
			while (!feof(queryFile) && queryChar != '\n') {
				queryChar = (char)getc(queryFile);
			}
			if (feof(queryFile)) {
				break;
			}
			queryChar = (char)getc(queryFile);
			queryPatternLength = 0;
			while (!feof(queryFile) && queryPatternLength < BWTSW_MAX_QUERY_LENGTH && queryChar != '>') {
				// Get query pattern
				if (queryChar != '\n') {
					if (queryPatternLength >= queryPatternAllocated) {
						if (queryPatternLength < LARGE_QUERY_ALLOCATION_UNIT) {
							queryPattern = MMUnitReallocate(queryPattern, queryPatternAllocated + SMALL_QUERY_ALLOCATION_UNIT, queryPatternAllocated);
							queryPatternAllocated += SMALL_QUERY_ALLOCATION_UNIT;
						} else {
							queryPattern = MMUnitReallocate(queryPattern, queryPatternAllocated + LARGE_QUERY_ALLOCATION_UNIT, queryPatternAllocated);
							queryPatternAllocated += LARGE_QUERY_ALLOCATION_UNIT;
						}
					}
					queryPattern[queryPatternLength] = queryChar;
					queryPatternLength++;
				}
				queryChar = (char)getc(queryFile);
			}

			queryStartTime = totalQueryTime;

			// Mask query pattern; BWTSW only support hard masking

			if (!MaskLowerCase) {
				// Change lower case in query input to upper case
				for (i=0; i<queryPatternLength; i++) {
					if (queryPattern[i] >= 'a' && queryPattern[i] <= 'z') {
						queryPattern[i] += (unsigned char)('A' - 'a');
					}
				}
			}

			// Apply dust masking; masked characters are set to lower case
			if (DustMasking) {
				blast_dust (queryPattern, queryPatternLength, DustLevel, DustWindow, 1);
			}

			// Initialize frequecies
			for (i=0; i<16; i++) {
				queryCharFreq[i] = 0;
			}

			// Convert text to code and mask the lowercase
			numUnmaskedChar = 0;
			for (i=0; i<queryPatternLength; i++) {
				if (queryPattern[i] >= 'a' && queryPattern[i] <= 'z') {
					queryPattern[i] = lowercaseDnaCharIndex;
				} else {
					queryPattern[i] = charMap[queryPattern[i]];
					numUnmaskedChar++;
				}
				queryCharFreq[queryPattern[i]]++;
			}

			// Calculate probabilities
			for (i=0; i<16; i++) {
				queryCharProb[i] = (double)queryCharFreq[i] / (queryPatternLength - queryCharFreq[lowercaseDnaCharIndex]);
			}
			queryCharProb[lowercaseDnaCharIndex] = 0;

			numOfContext = 0;

			if (QueryStrand == QUERY_BOTH_STRAND || QueryStrand == QUERY_POS_STRAND) {
				context[numOfContext].contextNum = 0;
				context[numOfContext].numOfHit = 0;
				numOfContext++;
			}

			if (QueryStrand == QUERY_BOTH_STRAND || QueryStrand == QUERY_NEG_STRAND) {
				context[numOfContext].contextNum = 1;
				context[numOfContext].numOfHit = 0;
				numOfContext++;
			}

			numOfQueryPattern++;

			// Initialize statistics
			initializeHSPstatistic(hsp->dnaLength, hsp->numOfSeq, hsp->minSeqLength, dbCharProb, 
								   queryPatternLength, numOfContext, queryCharProb, scoringMatrix);

			if (PrintProgress) {
				Socketprintf("%s : Cutoff score = %d, no. of char = %d/%d\n", queryPatternName, (int)calcGapCutoffScore(), numUnmaskedChar, queryPatternLength);
				fflush(stdout);
			}

			numOfHitForQuery = 0;

			MMPoolReset(alignmentPool);

			// Process all contexts
			for (currentContext=0; currentContext<numOfContext; currentContext++) {

				// Handle reverse complement in place to save memory
				if (contextInfo[context[currentContext].contextNum].reversed && contextInfo[context[currentContext].contextNum].complemented) {
					for (i=0; i<queryPatternLength/2; i++) {
						c = charMap[complementMap[dnaChar[queryPattern[i]]]];
						queryPattern[i] = charMap[complementMap[dnaChar[queryPattern[queryPatternLength - 1 - i]]]];
						queryPattern[queryPatternLength - 1 - i] = c;
					}
					if (queryPatternLength % 2 == 1) {
						queryPattern[queryPatternLength/2] = charMap[complementMap[dnaChar[queryPattern[queryPatternLength/2]]]];
					}
				}

				workingMemoryUsed = 0;
				MMPoolReset(workingMemoryPool);


				saIndexGroup = (SaIndexGroupNew*)(workingMemory + workingMemoryUsed);

				numOfSaIndexGroup = BWTGappedDPDBvsQuery(bwt, queryPattern, queryPatternLength,
												 workingMemory + workingMemoryUsed, WorkingMemorySize, &numOfQueryPos,
												stat_scoreBlock.match, stat_scoreBlock.mismatch,
												-stat_scoreBlock.gapOpen, -stat_scoreBlock.gapExt,
												(int)calcGapCutoffScore(), 
												&bwtDPStatistics,
												PrintProgressDepth);
											
				if (numOfSaIndexGroup > 0) {

					// Reconstruct the result from DP on index
					queryPosList = MMUnitAllocate(numOfQueryPos * sizeof(int));
					queryPosListIndex = MMUnitAllocate((numOfSaIndexGroup+1) * sizeof(int));
					j = 0;
					numOfDPPoint = 0;
					workingMemoryPointer = workingMemory + workingMemoryUsed;
					for (i=0; i<numOfSaIndexGroup; i++) {
						saIndexGroup[i].startSaIndex = *(unsigned int*)workingMemoryPointer;
						workingMemoryPointer += sizeof(unsigned int);
						saIndexGroup[i].numOfMatch = *(unsigned int*)workingMemoryPointer;
						workingMemoryPointer += sizeof(unsigned int);
						saIndexGroup[i].posQuery = *(unsigned int*)workingMemoryPointer;
						workingMemoryPointer += sizeof(unsigned int);
						saIndexGroup[i].info = *(unsigned int*)workingMemoryPointer;
						workingMemoryPointer += sizeof(unsigned int);
						queryPosListIndex[i] = j;
						n = *(int*)workingMemoryPointer;
						numOfDPPoint += saIndexGroup[i].numOfMatch * n;
						workingMemoryPointer += sizeof(int);
						for (k=0; k<n; k++) {
							queryPosList[j] = *(unsigned int*)workingMemoryPointer;
							workingMemoryPointer += sizeof(unsigned int);
							j++;
						}
					}
					queryPosListIndex[i] = j;
					totalNumOfDPPoint += numOfDPPoint;

					workingMemoryUsed += numOfSaIndexGroup * sizeof(SaIndexGroupNew);
					totalSaIndexRange += numOfSaIndexGroup;

				}

				if (PrintPerformanceStatistics) {
					elapsedTime = getElapsedTime(startTime) - totalQueryTime;
					totalQueryTime += elapsedTime;
					totalSearchTime += elapsedTime;
				}

				if (numOfSaIndexGroup > 0) {

					hitList = (HitList*)(workingMemory + workingMemoryUsed);

					context[currentContext].numOfHit = BWTDPHit(bwt, saIndexGroup, numOfSaIndexGroup,
																0, &i,
																FALSE, workingMemory + workingMemoryUsed, WorkingMemorySize - workingMemoryUsed, 
																&bwtSaRetrievalStatistics);
					context[currentContext].numOfDPPoint = numOfDPPoint;

					if (i < numOfSaIndexGroup) {
						fprintf(stderr, "Not enough memory allocated!\n");
						exit(1);
					}

					depthBit = ceilLog2(BWTDP_MAX_SUBSTRING_LENGTH);
					depthMask = ALL_ONE_MASK >> (BITS_IN_WORD - depthBit);
					for (i=0; i<context[currentContext].numOfHit; i++) {
						hitList[i].posText += hitList[i].info & depthMask;	// Adjust text position with length of substring
						hitList[i].info = (hitList[i].info >> depthBit);	// Recover SA index group index
					}
					totalHitGenerated += context[currentContext].numOfHit;

				} else {

					context[currentContext].numOfHit = 0;

				}

				if (PrintPerformanceStatistics) {
					elapsedTime = getElapsedTime(startTime) - totalQueryTime;
					totalQueryTime += elapsedTime;
					totalTextDecodeTime += elapsedTime;
				}

				// Allocate memory for gapped hit
				// hitList is copied to the end of allocated memory
				// and working memory is reused
				if (context[currentContext].numOfHit > 0) {

					QSort(hitList, context[currentContext].numOfHit, sizeof(HitList), HitListPosTextDescOrder);

					context[currentContext].gappedHitList = MMUnitAllocate(context[currentContext].numOfDPPoint * sizeof(GappedHitList));
					memcpy((char*)context[currentContext].gappedHitList
								  + context[currentContext].numOfDPPoint * sizeof(GappedHitList) 
								  - context[currentContext].numOfHit * sizeof(HitList), 
								  hitList, context[currentContext].numOfHit * sizeof(HitList));
					hitList = (HitList*)((char*)context[currentContext].gappedHitList 
															+ context[currentContext].numOfDPPoint * sizeof(GappedHitList) 
															- context[currentContext].numOfHit * sizeof(HitList));

					numOfGappedHitWithTraceback = HSPDPDBvsQuery(hsp->packedDNA, hitList, context[currentContext].numOfHit,
																	  hsp->seqOffset, hsp->ambiguity, hsp->numOfSeq,
																	  queryPattern, queryPatternLength,
																	  queryPosListIndex, queryPosList,
																	  context[currentContext].gappedHitList, context[currentContext].numOfDPPoint,
																	  workingMemoryPool, alignmentPool,
																	  stat_scoreBlock.match, stat_scoreBlock.mismatch,
																	  -stat_scoreBlock.gapOpen, -stat_scoreBlock.gapExt,
																	  calcGapCutoffScore(), stat_expectationValue);

					totalGappedHitWithTraceback += numOfGappedHitWithTraceback;

				} else {

					numOfGappedHitWithTraceback = 0;

				}

				if (numOfSaIndexGroup > 0) {
					// Free memory for DP points on query
					MMUnitFree(queryPosList, numOfQueryPos * sizeof(int));
					MMUnitFree(queryPosListIndex, (numOfSaIndexGroup+1) * sizeof(int));
				}

				// Final filter
				context[currentContext].numOfFinalHit = HSPFinalFilter(context[currentContext].gappedHitList, numOfGappedHitWithTraceback);
				totalFinalHit += context[currentContext].numOfFinalHit;

				if (PrintPerformanceStatistics) {
					elapsedTime = getElapsedTime(startTime) - totalQueryTime;
					totalQueryTime += elapsedTime;
					totalGappedExtensionTime2 += elapsedTime;
				}

				numOfHitForQuery += context[currentContext].numOfFinalHit;

				// Revert reverse complement in place
				if (contextInfo[context[currentContext].contextNum].reversed && contextInfo[context[currentContext].contextNum].complemented) {
					for (i=0; i<queryPatternLength/2; i++) {
						c = charMap[complementMap[dnaChar[queryPattern[i]]]];
						queryPattern[i] = charMap[complementMap[dnaChar[queryPattern[queryPatternLength - 1 - i]]]];
						queryPattern[queryPatternLength - 1 - i] = c;
					}
					if (queryPatternLength % 2 == 1) {
						queryPattern[queryPatternLength/2] = charMap[complementMap[dnaChar[queryPattern[queryPatternLength/2]]]];
					}
				}

			}

			// Combine final hits for all contexts from a query
			for (i=0; i<hsp->numOfSeq; i++) {
				dbScore[i] = 0;
				dbOrder[i] = i;
			}

			// Get the best score for each DB sequence
			for (currentContext=0; currentContext<numOfContext; currentContext++) {
				for (j=0; j<context[currentContext].numOfFinalHit; j++) {
					if (context[currentContext].gappedHitList[j].score > dbScore[context[currentContext].gappedHitList[j].dbSeqIndex]) {
						dbScore[context[currentContext].gappedHitList[j].dbSeqIndex] = context[currentContext].gappedHitList[j].score;
					}
				}
			}

			// Sort DB sequence according to dbScore
			for (i=1; i<hsp->numOfSeq; i++) {
				for (j=i; j>0 && dbScore[dbOrder[j-1]] < dbScore[dbOrder[j]]; j--) {
					swap(dbOrder[j-1], dbOrder[j], n);
				}
			}
			for (i=0; i<hsp->numOfSeq; i++) {
				dbRank[dbOrder[i]] = i;
			}

			for (currentContext=0; currentContext<numOfContext; currentContext++) {
				// Replace dbSeqIndex with dbRank + context
				for (j=0; j<context[currentContext].numOfFinalHit; j++) {
					context[currentContext].gappedHitList[j].dbSeqIndex = dbRank[context[currentContext].gappedHitList[j].dbSeqIndex] | (context[currentContext].contextNum << (BITS_IN_WORD - CONTEXT_BIT));
				}
			}

			HSPPrintQueryHeader(outputFile, OutputFormat, queryPatternName, queryPatternLength,
								dbOrder, dbScore, hsp->annotation, hsp->numOfSeq);
			if (alignFile != NULL) {
				HSPPrintQueryHeader(alignFile, OUTPUT_PAIRWISE, queryPatternName, queryPatternLength,
									dbOrder, dbScore, hsp->annotation, hsp->numOfSeq);
			}

			MMPoolReset(workingMemoryPool);

			if (numOfHitForQuery > 0) {

				finalHitListAllContext = MMUnitAllocate(numOfHitForQuery * sizeof(GappedHitListWithAlignment));
				n = 0;
				for (currentContext=0; currentContext<numOfContext; currentContext++) {
					if (context[currentContext].numOfFinalHit > 0) {
						memcpy(finalHitListAllContext + n, context[currentContext].gappedHitList, context[currentContext].numOfFinalHit * sizeof(GappedHitListWithAlignment));
						n += context[currentContext].numOfFinalHit;
					}
					if (context[currentContext].numOfHit > 0) {
						MMUnitFree(context[currentContext].gappedHitList, context[currentContext].numOfDPPoint * sizeof(GappedHitListWithAlignment));
					}
				}

				// Sort final hit by dbSeqIndex + score
				QSort(finalHitListAllContext, numOfHitForQuery, sizeof(GappedHitList), FinalHitAllContextDbSeqIndexScorePosTextOrder);

				// Store evalue histogram
				HSPCountEvalueToHistogram(histogram, finalHitListAllContext, numOfHitForQuery, TRUE);

				// Output results
				HSPPrintAlignment(workingMemoryPool, outputFile, finalHitListAllContext, numOfHitForQuery,
								  OutputFormat, contextInfo, 
								  charMap, complementMap,
								  queryPatternName, queryPattern, queryPatternLength,
								  dbOrder, hsp->seqOffset, hsp->annotation);
				if (alignFile != NULL) {
					HSPPrintAlignment(workingMemoryPool, alignFile, finalHitListAllContext, numOfHitForQuery,
									  OUTPUT_PAIRWISE, contextInfo, 
									  charMap, complementMap,
									  queryPatternName, queryPattern, queryPatternLength,
									  dbOrder, hsp->seqOffset, hsp->annotation);
				}

				// Alignment and auxiliary text are freed through alignmentPool
				MMUnitFree(finalHitListAllContext, numOfHitForQuery * sizeof(GappedHitListWithAlignment));

			} else {

				HSPPrintNoAlignment(outputFile, OutputFormat);
				if (alignFile != NULL) {
					HSPPrintNoAlignment(alignFile, OUTPUT_PAIRWISE);
				}

				if (PrintProgress) {
					Socketprintf("No alignment is found!\n");
					fflush(stdout);
				}

			}

			if (PrintPerformanceStatistics) {
				elapsedTime = getElapsedTime(startTime) - totalQueryTime;
				totalQueryTime += elapsedTime;
				totalPrintTime += elapsedTime;
			}

			if (timingFile != NULL) {
				fprintf(timingFile, "%s, %s, %9.4f, %d\n", "BWT-SW", queryPatternName, totalQueryTime - queryStartTime, numUnmaskedChar);
			}

			totalNumOfUnmaskedChar += numUnmaskedChar;
			totalNumOfChar += queryPatternLength;

		}
		
		HSPPrintTrailer(outputFile, OutputFormat, DatabaseName, hsp->numOfSeq, hsp->dnaLength);
		if (alignFile != NULL) {
			HSPPrintTrailer(alignFile, OUTPUT_PAIRWISE, DatabaseName, hsp->numOfSeq, hsp->dnaLength);
		}

		fclose(queryFile);
		fclose(outputFile);
		if (alignFile != NULL) {
			fclose(alignFile);
		}
		if (timingFile != NULL) {
			fclose(timingFile);
		}

		// printing the histogram on the gap hit count, based on the evalue
		if (PrintHistogram) {
			HSPPrintHistogram(stdout, histogram);
		}

		if (PrintDPStatistics) {
			BWTPrintDPStatistics(stdout, &bwtDPStatistics);
		}

		Socketprintf("\nFinished query of %d sequences. Total no. of characters = %lld/%lld.\n\n", numOfQueryPattern, totalNumOfUnmaskedChar, totalNumOfChar);
		if (PrintPerformanceStatistics) {
			Socketprintf("Search time             = %9.4f s", totalSearchTime);
			Socketprintf("   SA index range = %9lld\n", totalSaIndexRange);
			Socketprintf("Text decode time        = %9.4f s", totalTextDecodeTime);
			Socketprintf("   Hit generated  = %9lld (%lld)\n", totalHitGenerated, totalNumOfDPPoint);
			Socketprintf("Final DP with traceback = %9.4f s", totalGappedExtensionTime2);
			Socketprintf("   Final hit      = %9lld (%lld)\n", totalFinalHit, totalGappedHitWithTraceback);
			Socketprintf("Print time              = %9.4f s", totalPrintTime);
			Socketprintf("\n");
		} else {
			elapsedTime = getElapsedTime(startTime) - totalQueryTime;
			totalQueryTime += elapsedTime;
		}
		Socketprintf("\nDatabase load time      = %9.4f s", databaseLoadTime);
		Socketprintf("\n");
		Socketprintf("Total query time        = %9.4f s", totalQueryTime);
		Socketprintf("\n");
		Socketprintf("Total elapsed time      = %9.4f s", totalQueryTime + databaseLoadTime);
		Socketprintf("\n\n");

		if (LoadIndex) {
			// End connection
			Socketprintf("");	// Special message for client to exit with 0 return code
			SocketEndConnection(bwtServerSocket);
		} else {
			// Exit loop if this is a client
			break;
		}

	}

	MMUnitFree(dbScore, hsp->numOfSeq * sizeof(int));
	MMUnitFree(dbOrder, hsp->numOfSeq * sizeof(int));
	MMUnitFree(dbRank, hsp->numOfSeq * sizeof(int));
	MMPoolFree(workingMemoryPool);
	MMPoolFree(alignmentPool);

	BWTFreeDPStatistics(&bwtDPStatistics);

	MMUnitFree(queryPattern, queryPatternAllocated);

	HSPFree(mmPool, hsp);
	BWTFree(mmPool, bwt);

	MMPoolFree(mmPool);

//	MMMasterPrintReport(stdout, TRUE, TRUE, FALSE);

	if (LoadIndex) {
		SocketFree(bwtServerSocket);
	}

	iniparser_freedict(programInput);

	return 0;

}

dictionary *ParseInput(int argc, char** argv) {

	dictionary *programInput;
	char t1[3] = "-c";	// specify that this is a boolean type parameter; no following argument
	char t2[3] = "-U";	// specify that this is a boolean type parameter; no following argument
	char t3[3] = "-L";	// specify that this is a boolean type parameter; no following argument
	char t4[3] = "-X";	// specify that this is a boolean type parameter; no following argument
	char *d[4];

	char *tempString;
	int len;
	int i;

	d[0] = t1;
	d[1] = t2;
	d[2] = t3;
	d[3] = t4;

	programInput = paraparser_load(argc, argv, 4, d);	// 4 parameters are boolean type

	// Whether confirmation is needed
	Confirmation = iniparser_find_entry(programInput, "parameter:-c");
	LoadIndex = iniparser_find_entry(programInput, "parameter:-L");
	UnloadIndex = iniparser_find_entry(programInput, "parameter:-X");

	if (LoadIndex && UnloadIndex) {
		PrintHelp();
		exit(1);
	}

	if (UnloadIndex) {
		return programInput;
	}

	// Get database, query name and output file name
	if (!iniparser_find_entry(programInput, "parameter:-d")) {
		// Database name may be entered through argument
		if (!iniparser_find_entry(programInput, "argument:1")) {
			PrintHelp();
			exit(1);
		}
		iniparser_copystring(programInput, "argument:1", DatabaseName, DatabaseName, MAX_FILENAME_LEN);
		if (LoadIndex) {
			return programInput;
		}
		if (!iniparser_find_entry(programInput, "parameter:-i")) {
			// Query name may be entered through argument
			if (!iniparser_find_entry(programInput, "argument:2")) {
				PrintHelp();
				exit(1);
			}
			iniparser_copystring(programInput, "argument:2", QueryFileName, QueryFileName, MAX_FILENAME_LEN);
			if (!iniparser_find_entry(programInput, "parameter:-o")) {
				// Output file name may be entered through argument
				if (!iniparser_find_entry(programInput, "argument:3")) {
					PrintHelp();
					exit(1);
				}
				iniparser_copystring(programInput, "argument:3", OutputFileName, OutputFileName, MAX_FILENAME_LEN);
			} else {
				iniparser_copystring(programInput, "parameter:-o", OutputFileName, OutputFileName, MAX_FILENAME_LEN);
			}
		} else {
			iniparser_copystring(programInput, "parameter:-i", QueryFileName, QueryFileName, MAX_FILENAME_LEN);
			if (!iniparser_find_entry(programInput, "parameter:-o")) {
				// Output file name may be entered through argument
				if (!iniparser_find_entry(programInput, "argument:2")) {
					PrintHelp();
					exit(1);
				}
				iniparser_copystring(programInput, "argument:2", OutputFileName, OutputFileName, MAX_FILENAME_LEN);
			} else {
				iniparser_copystring(programInput, "parameter:-o", OutputFileName, OutputFileName, MAX_FILENAME_LEN);
			}
		}
	} else {
		iniparser_copystring(programInput, "parameter:-d", DatabaseName, DatabaseName, MAX_FILENAME_LEN);
		if (LoadIndex) {
			return programInput;
		}
		if (!iniparser_find_entry(programInput, "parameter:-i")) {
			// Query name may be entered through argument
			if (!iniparser_find_entry(programInput, "argument:1")) {
				PrintHelp();
				exit(1);
			}
			iniparser_copystring(programInput, "argument:1", QueryFileName, QueryFileName, MAX_FILENAME_LEN);
			if (!iniparser_find_entry(programInput, "parameter:-o")) {
				// Output file name may be entered through argument
				if (!iniparser_find_entry(programInput, "argument:2")) {
					PrintHelp();
					exit(1);
				}
				iniparser_copystring(programInput, "argument:2", OutputFileName, OutputFileName, MAX_FILENAME_LEN);
			} else {
				iniparser_copystring(programInput, "parameter:-o", OutputFileName, OutputFileName, MAX_FILENAME_LEN);
			}
		} else {
			iniparser_copystring(programInput, "parameter:-i", QueryFileName, QueryFileName, MAX_FILENAME_LEN);
			if (!iniparser_find_entry(programInput, "parameter:-o")) {
				// Output file name may be entered through argument
				if (!iniparser_find_entry(programInput, "argument:1")) {
					PrintHelp();
					exit(1);
				}
				iniparser_copystring(programInput, "argument:1", OutputFileName, OutputFileName, MAX_FILENAME_LEN);
			} else {
				iniparser_copystring(programInput, "parameter:-o", OutputFileName, OutputFileName, MAX_FILENAME_LEN);
			}
		}
	}

	iniparser_copystring(programInput, "parameter:-time", TimingFileName, TimingFileName, MAX_FILENAME_LEN);
	iniparser_copystring(programInput, "parameter:-align", AlignFileName, AlignFileName, MAX_FILENAME_LEN);

	// Whether search parameter file is specified
	AlternateQueryParameter = iniparser_find_entry(programInput, "parameter:-p");
	if (AlternateQueryParameter) {
		tempString = iniparser_getstr(programInput, "parameter:-p");
		len = (int)strlen(tempString);
		if (len > 5 && strncmp(tempString + len - 5, ".parm", 5) == 0) {
			if (len > MAX_FILENAME_LEN) {
				fprintf(stderr, "Search parameter file name is too long!\n");
				exit(1);
			}
			sprintf(QueryParameterFileName, "%s", tempString);
		} else {
			if (strlen(tempString) + 5 > MAX_FILENAME_LEN) {
				fprintf(stderr, "Search parameter file name is too long!\n");
				exit(1);
			}
			sprintf(QueryParameterFileName, "%s.parm", tempString);
		}
	}

	QueryStrand = iniparser_getint(programInput, "parameter:-S", QueryStrand);

	MaskLowerCase = iniparser_find_entry(programInput, "parameter:-U");
	MaskOption = iniparser_getstr(programInput, "parameter:-F");
	if (MaskOption != NULL) {
		for (i=0; i<(int)strlen(MaskOption); i++) {
			if (MaskOption[i] == 'T' || MaskOption[i] == 't' || MaskOption[i] == 'Y' || MaskOption[i] == 'y' || MaskOption[i] == '1') {
				DustMasking = TRUE;
			}
			if (MaskOption[i] == 'F' || MaskOption[i] == 'f' || MaskOption[i] == 'N' || MaskOption[i] == 'n' || MaskOption[i] == '0') {
				DustMasking = FALSE;
			}
		}
	}

	InputEValue = iniparser_getdouble(programInput,"parameter:-e", InputEValue);
	InputMatch = iniparser_getint(programInput,"parameter:-r", InputMatch);
	InputMismatch = iniparser_getint(programInput,"parameter:-q", InputMismatch);
	InputGapOpen = iniparser_getint(programInput,"parameter:-G", InputGapOpen);
	InputGapExtend = iniparser_getint(programInput,"parameter:-E", InputGapExtend);

	OutputFormat = iniparser_getint(programInput, "parameter:-m", OutputFormat);

	return programInput;

}

void ParseIniFile(char *iniFileName) {

	dictionary *ini;

	printf("Loading %s ..", iniFileName);
	ini = iniparser_load(iniFileName, FALSE);
	if (ini == NULL) {
		printf("not found.\n");
		return;
	}
	printf("done.\n");

	// Memory parameters
	WorkingMemorySize = iniparser_getuint(ini, "Memory:WorkingMemorySize", WorkingMemorySize);

	// Performance statistics parameters
	PrintProgress = iniparser_getboolean(ini, "PerformanceStatistics:PrintProgress", PrintProgress);
	PrintPerformanceStatistics = iniparser_getboolean(ini, "PerformanceStatistics:PrintPerformanceStatistics", PrintPerformanceStatistics);
	PrintHistogram = iniparser_getboolean(ini, "PerformanceStatistics:PrintHistogram", PrintHistogram);
	HistogramMinEvalue = iniparser_getdouble(ini, "PerformanceStatistics:HistogramMinEvalue", HistogramMinEvalue);
	PrintDPStatistics = iniparser_getboolean(ini, "PerformanceStatistics:PrintDPStatistics", PrintDPStatistics);
	PrintProgressDepth = iniparser_getuint(ini, "PerformanceStatistics:PrintProgressDepth", PrintProgressDepth);

	// Database parameters
	iniparser_copystring(ini, "Database:Location", DatabaseLocation, DatabaseLocation, MAX_FILENAME_LEN);
	iniparser_copystring(ini, "Database:AnnotationFileName", AnnotationFileName, AnnotationFileName, MAX_FILENAME_LEN);
	iniparser_copystring(ini, "Database:AmbiguityFileName", AmbiguityFileName, AmbiguityFileName, MAX_FILENAME_LEN);
	iniparser_copystring(ini, "Database:PackedDNAFileName", PackedDNAFileName, PackedDNAFileName, MAX_FILENAME_LEN);
	iniparser_copystring(ini, "Database:BWTCodeFileName", BWTCodeFileName, BWTCodeFileName, MAX_FILENAME_LEN);
	iniparser_copystring(ini, "Database:BWTOccValueFileName", BWTOccValueFileName, BWTOccValueFileName, MAX_FILENAME_LEN);
	iniparser_copystring(ini, "Database:SaValueFileName", SaValueFileName, SaValueFileName, MAX_FILENAME_LEN);

	iniparser_freedict(ini);

}

void ProcessIni() {

	ProcessFileName(AnnotationFileName, AnnotationFileName, DatabaseLocation, DatabaseName);
	ProcessFileName(AmbiguityFileName, AmbiguityFileName, DatabaseLocation, DatabaseName);
	ProcessFileName(PackedDNAFileName, PackedDNAFileName, DatabaseLocation, DatabaseName);
	ProcessFileName(BWTCodeFileName, BWTCodeFileName, DatabaseLocation, DatabaseName);
	ProcessFileName(BWTOccValueFileName, BWTOccValueFileName, DatabaseLocation, DatabaseName);
	ProcessFileName(SaValueFileName, SaValueFileName, DatabaseLocation, DatabaseName);

}

void ValidateIni() {

	if (AnnotationFileName[0] == '\0') {
		fprintf(stderr, "Annotation file is not specified!\n");
		exit(1);
	}
	if (AmbiguityFileName[0] == '\0') {
		fprintf(stderr, "Ambiguity file is not specified!\n");
		exit(1);
	}	
	if (PackedDNAFileName[0] == '\0') {
		fprintf(stderr, "Packed DNA file is not specified!\n");
		exit(1);
	}
	if (BWTCodeFileName[0] == '\0') {
		fprintf(stderr, "BWT Code file is not specified!\n");
		exit(1);
	}
	if (BWTOccValueFileName[0] == '\0') {
		fprintf(stderr, "BWT Occ value file is not specified!\n");
		exit(1);
	}
	if (SaValueFileName[0] == '\0') {
		fprintf(stderr, "SA value file is not specified!\n");
		exit(1);
	}

}


void PrintIni() {

	char boolean[2];

	boolean[0] = 'N';
	boolean[1] = 'Y';

	printf("Annotation file          : %s\n", AnnotationFileName);
	printf("Ambiguity file           : %s\n", AmbiguityFileName);
	printf("Packed DNA file          : %s\n", PackedDNAFileName);
	printf("BWT Code file            : %s\n", BWTCodeFileName);
	printf("BWT Occ value file       : %s\n", BWTOccValueFileName);
	printf("SA value file            : %s\n", SaValueFileName);
	printf("\n");

}

void PrintShortDesc() {

	printf("BWT-SW v1.0, Copyright (C) 2006, Wong Chi Kwong.\n");
	printf("BWT-SW comes with ABSOLUTELY NO WARRENTY.\n");
	printf("BWT-SW is free software, and you are welcome to\n");
	printf("redistribute it under certain conditions.\n");
	printf("For details type bwtsw.\n");
	printf("\n");

}

void PrintHelp() {

	printf("BWT-SW v1.0, Copyright (C) 2006, Wong Chi Kwong.\n");
	printf("\n");

	printf("This program is free software; you can redistribute it and/or\n");
	printf("modify it under the terms of the GNU General Public License\n");
	printf("as published by the Free Software Foundation; either version 2\n");
	printf("of the License, or (at your option) any later version.\n");
	printf("\n");

	printf("This program is distributed in the hope that it will be useful,\n");
	printf("but WITHOUT ANY WARRANTY; without even the implied warranty of\n");
	printf("MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n");
	printf("GNU General Public License for more details.\n");
	printf("\n");

	printf("You should have received a copy of the GNU General Public License\n");
	printf("along with this program; if not, write to the Free Software\n");
	printf("Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.\n");
	printf("\n");

	printf("Syntax: bwtsw <database> <query> <output>\n");
	printf("              [-c Confirm]\n");
	printf("              [-L Load server process]\n");
	printf("              [-X Unload server process]\n");
	printf("              [-U Mask lower case]\n");
	printf("               -F T/F Mask query by Dust\n");
	printf("               -S <strand> : 1 - Upper; 2 - Lower; 3 - Both (default)\n");
	printf("               -e <min E-Value>\n");
	printf("               -r <reward for a match>\n");
	printf("               -q <penalty for a mismatch>\n");
	printf("               -G <cost to open a gap>\n");
	printf("               -E <cost to extend a gap>\n");
	printf("               -m <output format> : 0 - Pairwise (default); 8 - Tabular;\n");
	printf("                                    9 - Tabular with comment\n");
	printf("               -p <query parameter file>\n");

}

void ParseQueryParameterFile(char *queryParameterFileName) {

	dictionary *queryParameter;

	printf("Loading %s ..", queryParameterFileName);
	queryParameter = iniparser_load(queryParameterFileName, FALSE);
	if (queryParameter == NULL) {
		printf("not found.\n");
		return;
	}
	printf("done.\n");

	// Hit Scoring parameters
	stat_scoreBlock.match = iniparser_getint(queryParameter, "HitScoring:Match", stat_scoreBlock.match);
	stat_scoreBlock.mismatch = iniparser_getint(queryParameter, "HitScoring:Mismatch", stat_scoreBlock.mismatch);
	stat_scoreBlock.gapOpen = iniparser_getint(queryParameter, "HitScoring:GapOpen",stat_scoreBlock.gapOpen);
	stat_scoreBlock.gapExt = iniparser_getint(queryParameter, "HitScoring:GapExtension", stat_scoreBlock.gapExt);

	if (InputMatch != NON_INPUT_VALUE) {
		stat_scoreBlock.match = InputMatch;
	}
	if (InputMismatch != NON_INPUT_VALUE) {
		stat_scoreBlock.mismatch = InputMismatch;
	}
	if (InputGapOpen != NON_INPUT_VALUE) {
		stat_scoreBlock.gapOpen = InputGapOpen;
	}
	if (InputGapExtend != NON_INPUT_VALUE) {
		stat_scoreBlock.gapExt = InputGapExtend;
	}

	// Dust parameters
	DustLevel = iniparser_getint(queryParameter, "Dust:DustLevel", DustLevel);
	DustWindow = iniparser_getint(queryParameter, "Dust:DustWindow", DustWindow);

	// Expectation parameters
	stat_expectationValue = iniparser_getdouble(queryParameter,"ExpectationValue:ExpectationValue",stat_expectationValue);
	if (InputEValue < NON_INPUT_VALUE / 10) {
		stat_expectationValue = InputEValue;
	}

	iniparser_freedict(queryParameter);

}

void ValidateQueryParameters() {

	if (QueryStrand != QUERY_BOTH_STRAND && QueryStrand != QUERY_POS_STRAND && QueryStrand != QUERY_NEG_STRAND) {
		fprintf(stderr, "Query strand must be 1=Pos strand only; 2=Neg strand only; 3=Both strands!\n");
		exit(1);
	}
	if (OutputFormat != OUTPUT_PAIRWISE && OutputFormat != OUTPUT_TABULAR && OutputFormat != OUTPUT_TABULAR_COMMENT) {
		fprintf(stderr, "Only -m 0, -m 8 and -m 9 output formats are supported!\n");
		exit(1);
	}

	if (stat_scoreBlock.match <= 0) {
		fprintf(stderr, "Reward for a match must be positive!\n");
		exit(1);
	}
	if (stat_scoreBlock.mismatch >= 0) {
		fprintf(stderr, "Penalty for a mismatch must be negative!\n");
		exit(1);
	}
	if ((-stat_scoreBlock.mismatch) / stat_scoreBlock.match < 3) {
		fprintf(stderr, "Penalty for a mismatch must be at least 3 times the reward for a match!\n");
		exit(1);
	}
	if (stat_scoreBlock.gapOpen <= 0) {
		fprintf(stderr, "Cost for opening a gap must be positive!\n");
		exit(1);
	}
	if (stat_scoreBlock.gapExt <= 0) {
		fprintf(stderr, "Cost for extending a gap must be positive!\n");
		exit(1);
	}
	if (stat_scoreBlock.gapExt * 2 < -1 * stat_scoreBlock.mismatch) {
		fprintf(stderr, "Penaly of a mismatch must be at most the cost for extending a gap x 2!\n");
		exit(1);
	}

	if (stat_expectationValue <= 0.00000000000000000001)  {
		fprintf(stderr,"Expectation value must be a positive value > 0\n");
		exit(1);
	}

}

void PrintQueryParameters() {

	char boolean[2];

	boolean[0] = 'N';
	boolean[1] = 'Y';

	if (QueryStrand == QUERY_BOTH_STRAND) {
		printf("Query              :   Both strands\n");
	}
	if (QueryStrand == QUERY_POS_STRAND) {
		printf("Query              :   Positive strand only\n");
	}
	if (QueryStrand == QUERY_NEG_STRAND) {
		printf("Query              :   Negative strand only\n");
	}
	printf("\n");

	if (DustMasking) {
		printf("Dust masking       :   %c (Level = %d, Window = %d, Word = 3)\n", boolean[DustMasking], DustLevel, DustWindow);
	} else {
		printf("Dust masking       :   %c\n", boolean[DustMasking]);
	}
	printf("Mask lower case    :   %c\n", boolean[MaskLowerCase]);
	printf("\n");

	printf("Match Score        :   %d\n", stat_scoreBlock.match);
	printf("Mismatch Score     :   %d\n", stat_scoreBlock.mismatch);
	printf("Gap Open           :   %d\n", stat_scoreBlock.gapOpen);
	printf("Gap Extension      :   %d\n", stat_scoreBlock.gapExt);

	printf("\n");

	printf("Expectation Value  :   %5.2f\n", stat_expectationValue);

	printf("\n");

}


void ProcessFileName(char *outputFileName, const char *inputFileName, const char *databaseLocation,
					 const char *databaseName) {

	char tempChar[MAX_FILENAME_LEN];
	unsigned int i;

	if (inputFileName == NULL || inputFileName[0] == '\0' || inputFileName[0] == ' ' || inputFileName[0] == '-') {
		if (outputFileName != inputFileName) {
			outputFileName[0] = '-';
		}
		return;
	}

	if (strlen(databaseLocation) + strlen(databaseName) + strlen(inputFileName) > MAX_FILENAME_LEN) {
		fprintf(stderr, "File length is too long!\n");
		exit(1);
	}

	strncpy(tempChar, inputFileName, MAX_FILENAME_LEN);

	// locate the *
	for (i=0; i<MAX_FILENAME_LEN; i++) {
		if (tempChar[i] == '*') {
			break;
		}
	}
	if (i<MAX_FILENAME_LEN) {
		tempChar[i] = '\0';
		sprintf(outputFileName, "%s%s%s%s", databaseLocation, tempChar, databaseName, tempChar + i + 1);
	} else {
		sprintf(outputFileName, "%s%s", databaseLocation, tempChar);
	}

}

int HitListPosTextDescOrder(const void *hitList, const int index1, const int index2) {

	if (((HitList*)hitList + index1)->posText < ((HitList*)hitList + index2)->posText) {
		return 1;
	} else {
		return -1;
	}

}

int FinalHitAllContextDbSeqIndexScorePosTextOrder(const void *gappedHitList, const int index1, const int index2) {

	unsigned int dbSeqIndex1, dbSeqIndex2;

	dbSeqIndex1 = ((GappedHitList*)gappedHitList + index1)->dbSeqIndex & CONTEXT_MASK;
	dbSeqIndex2 = ((GappedHitList*)gappedHitList + index2)->dbSeqIndex & CONTEXT_MASK;

	if (dbSeqIndex1 != dbSeqIndex2) {
		return dbSeqIndex1 - dbSeqIndex2;
	} else {
		if (((GappedHitList*)gappedHitList + index1)->score != ((GappedHitList*)gappedHitList + index2)->score) {
			return ((GappedHitList*)gappedHitList + index2)->score - ((GappedHitList*)gappedHitList + index1)->score;
		} else {
			if (((GappedHitList*)gappedHitList + index1)->posText != ((GappedHitList*)gappedHitList + index2)->posText) {
				if (((GappedHitList*)gappedHitList + index1)->posText > ((GappedHitList*)gappedHitList + index2)->posText) {
					return 1;
				} else {
					return -1;
				}
			} else {
				return 0;
			}
		}
	}

}
