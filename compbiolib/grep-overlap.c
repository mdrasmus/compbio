#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/********* TAKES TWO LISTS OF REGIONS AND OUTPUTS ALL OVERLAPS ********/

/* INPUT FORMAT: NAME CHR START END [extra text] */

/********* BOTH LISTS MUST BE SORTED BY CHR (in STRCMP ORDER) THEN START!! ***********/

/* BEST PERFORMANCE IF SECOND LIST HAS BIGGER REGIONS THAN SECOND *****/

#define MAXLINELEN			5000
#define MAXLINELENSTR		"5000"
#define MAXCHRLEN			200

#define MAX(a,b)	(((a) > (b)) ? (a) : (b))
#define MIN(a,b)	(((a) < (b)) ? (a) : (b))

/* NOTE: this does NOT check to make sure the chromosome is the same !! */
#define OVERLAP(start1,end1,start2,end2) (MIN(end1,end2) - MAX(start1,start2) + 1)

/* given two chromosomes and a position checks to see if the first is not after the second */
#define ISINORDER(chr1, pos1, chr2, pos2)	(strcmp(chr1, chr2) < 0 || (strcmp(chr1, chr2) == 0 && (pos1) <= (pos2)))

#define asserte1(b,a1)			do {if (!(b)) {fprintf(stderr, (a1)); exit(1);}} while(0)
#define asserte2(b,a1,a2)		do {if (!(b)) {fprintf(stderr, (a1), (a2)); exit(1);}} while(0)
#define asserte3(b,a1,a2,a3)	do {if (!(b)) {fprintf(stderr, (a1), (a2), (a3)); exit(1);}} while(0)

/* when not doing strand, the strand may be left unfilled.. we don't care about that */
#define FILL_FMT(f,n)	\
	(f) = (char*) malloc(sizeof(char)*(4*(n)+14));  \
	for (i=0; i<(4*(n)); i+=4) \
		strcpy((f)+i, "%*s "); \
	strcpy((f)+4*(n), "%s %Ld %Ld %c")

#define FILL_REGION(reg, fmt) if (sscanf((reg)->line, fmt, (reg)->chr, &(reg)->start, &(reg)->end, &(reg)->str) == 3 || ((reg)->str != '+' && (reg)->str != '-')) (reg)->str = 0

#define MATCH_MODE_NORMAL				0
#define MATCH_MODE_ONLYMATCHES			1
#define MATCH_MODE_ONLYNONMATCHES		2
#define MATCH_MODE_PREFIX				3

struct region
{
	char line[MAXLINELEN+1];
	char chr[MAXCHRLEN+1];
	long long int start;
	long long int end;
	char str;
	struct region* next;
};


int main (int argc, char** argv)
{
	FILE *fp1, *fp2;

	/* this could have probably been implemented somewhat more quickly with some sort of
	 * sliding array to reduce the number of required mallocs... there are some subtleties 
	 * that make this non-trivial, and it looks like this is already very fast, but 
	 * something to consider */

	/* the first one on this list is used as the new item we read in */
	struct region* rlist = (struct region*) malloc (sizeof(struct region));

	/* we use and keep track of this rlist_end ptr so that our output is guaranteed to be
	 * first sorted by fp1 then fp2 */
	struct region* rlist_end = rlist;

	/* indicates with fp2 is empty */
	unsigned char fp2_empty = 0;

	/* use these to ensure that our input files are both sorted */
	char fp1_last_chr[MAXCHRLEN+1] = "";
	long long int fp1_last_start = 0;
	char fp2_last_chr[MAXCHRLEN+1] = "";
	long long int fp2_last_start = 0;
	int fp1_skip_col = 1;
	int fp2_skip_col = 1;

	char *fp1_fmt, *fp2_fmt;

	int i, padding = 0;

	int match_mode = MATCH_MODE_NORMAL;
	char sep_char = '|';
	char out_fmt[100];

	unsigned char match_r1_olen = 0;
	unsigned char match_r2_olen = 0;
	unsigned char require_str = 0; /* 2 for opposite strand */

	rlist->next = NULL;

	if (argc<3)
	{
		printf("USAGE: %s [options] <F1> <F2> > <OVERLAPS>\n", argv[0]);
		printf("Input format (F1 and F2): \n");
		printf("   [Skip Columns]Chr Start End [strand+-] [addl text] (whitespace separated)\n");
		printf("   Lines must be sorted first by Chr, then by Start\n");
		printf("   - may be used to indicate usage of stdin\n");
		printf("Output format: \n");
		printf("   L1|L2|len(R1)|len(R2)|OLen|Chr|OStart|OEnd\n");
		printf("   L1, L2 are the lines from F1 and F2, respectively\n");
		printf("   len(R1), len(R2) are the corresponding lengths\n");
		printf("   OLen, Chr, OStart, OEnd indicate the length and coordinates of the overlap (including the distance padding to S2)\n");
		printf("Options: \n");
		printf("-s char: separating character [default: |]\n");
		printf("-o: just show matching lines in F1\n");
		printf("-v: just show non-matching lines in F1\n");
		printf("-t/-T: require same/opposite strand match (note: if no strand is indicated, neither will match)\n");
		printf("-p: print all lines in F1 with 0/1 prefix indicating not match/match\n");
		printf("-w num: consider regions that are w distance apart as overlapping\n");
		printf("-c num, -c1 num, -c2 num: number of columns to skip for before the Chr for all input, F1 and F2, respectively [default: 1]\n");
		printf("-r1, -r2: require that OLen equal len(R1) or len(R2), respectively\n");

		return 1;
	}

	if (strcmp(argv[argc-2],"-") != 0)
		asserte2((fp1 = fopen(argv[argc-2], "r"))!=NULL, "Cannot open: %s\n", argv[argc-2]);
	else
		fp1 = stdin;

	if (strcmp(argv[argc-1],"-") != 0)
		asserte2((fp2 = fopen(argv[argc-1], "r"))!=NULL, "Cannot open: %s\n", argv[argc-1]);
	else
		fp2 = stdin;

	asserte1(fp1 != fp2, "Both input file streams cannot be stdin.\n");

	for (i=1; i<argc-2; i++)
		if (strcmp(argv[i],"-o") == 0)
			match_mode = MATCH_MODE_ONLYMATCHES;
		else if (strcmp(argv[i],"-v") == 0)
			match_mode = MATCH_MODE_ONLYNONMATCHES;
		else if (strcmp(argv[i],"-p") == 0)
			match_mode = MATCH_MODE_PREFIX;
		else if (strcmp(argv[i],"-w") == 0)
			padding = atoi(argv[++i]);
		else if (strcmp(argv[i],"-s") == 0)
			sep_char = argv[++i][0];
		else if (strcmp(argv[i],"-c") == 0)
			fp1_skip_col = fp2_skip_col = atoi(argv[++i]);
		else if (strcmp(argv[i],"-c1") == 0)
			fp1_skip_col = atoi(argv[++i]);
		else if (strcmp(argv[i],"-c2") == 0)
			fp2_skip_col = atoi(argv[++i]);
		else if (strcmp(argv[i],"-r1") == 0)
			match_r1_olen = 1; 
		else if (strcmp(argv[i],"-r2") == 0)
			match_r2_olen = 1; 
		else if (strcmp(argv[i],"-t") == 0)
			require_str = 1; 
		else if (strcmp(argv[i],"-T") == 0)
			require_str = 2; 
		else
			asserte1(0, "Invalid command line argument.\n");

	asserte1(padding >= 0, "-w argument must be >= 0.\n");

	/* fill out output format.. notice we allocate plenty of space for out_fmt above */
	if (match_mode == MATCH_MODE_NORMAL)
	{
		strcpy(out_fmt, "%s %s %Ld %Ld %Ld %s %Ld %Ld\n");

		for (i=0; out_fmt[i] != '\n'; i++)
			if (out_fmt[i] == ' ')
				out_fmt[i] = sep_char;
	}
	else if (match_mode == MATCH_MODE_PREFIX)
	{
		strcpy(out_fmt, "%d %s\n");
		out_fmt[2] = sep_char;
	}
	else
		strcpy(out_fmt, "%s\n");

	/* fill in the formats for input*/
	FILL_FMT(fp1_fmt,fp1_skip_col);
	FILL_FMT(fp2_fmt,fp2_skip_col);

	while (fscanf(fp1, "%" MAXLINELENSTR "[^\n]\n", rlist->line) != EOF)
	{
		/* used for -v and -o modes of operation */
		int has_match = 0;

		/* used for iterating through rlist */
		struct region* rlist_ptr = NULL;

		FILL_REGION(rlist, fp1_fmt);

		/* make sure that the data point makes sense */
		asserte3(rlist->start <= rlist->end, "%s:%s has start follow end.\n", argv[argc-2], rlist->line);

		/* check to make sure fp1 is sorted */
		asserte3(ISINORDER(fp1_last_chr, fp1_last_start, rlist->chr, rlist->start), "%s:%s is out of order.\n", argv[argc-2], rlist->line);
		strcpy(fp1_last_chr, rlist->chr);
		fp1_last_start = rlist->start;

		/* scan through the list, eliminate unneeded entries */
		for (rlist_ptr=rlist; rlist_ptr->next != NULL; )
			if (ISINORDER(rlist_ptr->next->chr, rlist_ptr->next->end+1, rlist->chr, rlist->start))
			{
				/* this item will never be overlapped with... delete it */
				struct region* temp = rlist_ptr->next->next;

				if (rlist_ptr->next == rlist_end)
					rlist_end = rlist_ptr;

				free(rlist_ptr->next);
				rlist_ptr->next = temp;
			}
			else
				rlist_ptr=rlist_ptr->next;

		/* read in new regions until a region that starts after the end of this region *
		 * is read in */
		while (
				/* if this file is empty, we cannot scan for any more new items */
				!fp2_empty && (
				
				/* if the list is empty... read in more */
				rlist->next == NULL || 
				
				/* or if what we read in from F1 ends after our stack (i.e. the
				 * last stack item starts before one after the item we just read in)
				 * we want them out of order so we can ensure we have nothing left
				 * in fp2 that will overlap with rlist */
				ISINORDER(rlist_end->chr, rlist_end->start, rlist->chr, rlist->end)
				)
		)
		{
			rlist_end->next = (struct region*) malloc (sizeof(struct region));

			rlist_end->next->next = NULL;
			
			if (fscanf(fp2, "%[^\n]\n", rlist_end->next->line) == EOF)
			{
				fp2_empty = 1;
				free(rlist_end->next);
				rlist_end->next = NULL;
			}
			else
			{
				FILL_REGION(rlist_end->next, fp2_fmt);

				/* make sure that the data point makes sense */
				asserte3(rlist_end->next->start <= rlist_end->next->end, "%s:%s has start follow end.\n", argv[argc-1], rlist_end->next->line);

				/* check to make sure fp2 is sorted */
				asserte3(ISINORDER(fp2_last_chr, fp2_last_start, rlist_end->next->chr, rlist_end->next->start), "%s:%s is out of order.\n", argv[argc-1], rlist_end->next->line);
				strcpy(fp2_last_chr, rlist_end->next->chr);
				fp2_last_start = rlist_end->next->start;

				/* add padding for the padding */
				rlist_end->next->start -= padding;
				rlist_end->next->end += padding;

				/* throw out if this entry is already behind us */
				if (ISINORDER(rlist_end->next->chr, rlist_end->next->end+1, rlist->chr, rlist->start))
				{
					free(rlist_end->next);
					rlist_end->next = NULL;
				}
				else
					rlist_end = rlist_end->next;
			}
		}

		/* scan through the list and print overlaps */
		for (rlist_ptr=rlist; rlist_ptr->next != NULL; rlist_ptr=rlist_ptr->next)
			/* this check is important but very rarely encountered! */
			if (strcmp(rlist_ptr->next->chr, rlist->chr) == 0)
			{
				long long int overlap_size = OVERLAP(rlist_ptr->next->start, rlist_ptr->next->end, rlist->start, rlist->end);
				if (overlap_size > 0
						&& (!match_r1_olen || overlap_size == (rlist->end - rlist->start + 1))
						&& (!match_r2_olen || overlap_size == (rlist_ptr->next->end - rlist_ptr->next->start + 1 - 2 * padding))
						&& (!require_str || (rlist->str && rlist_ptr->next->str && ((require_str == 1) == (rlist->str == rlist_ptr->next->str))))
				)
				{
					
					has_match = 1;

					if (match_mode == MATCH_MODE_NORMAL)
						printf(out_fmt, rlist->line, rlist_ptr->next->line, rlist->end - rlist->start + 1, rlist_ptr->next->end - rlist_ptr->next->start + 1 - 2 * padding, overlap_size, rlist_ptr->next->chr, MAX(rlist_ptr->next->start,rlist->start), MIN(rlist_ptr->next->end,rlist->end));
					else
						/* we can short circuit in match_mode 1, 2 and 3 */
						break;
				}
			}

		if ((match_mode == MATCH_MODE_ONLYMATCHES && has_match) || (match_mode == MATCH_MODE_ONLYNONMATCHES && !has_match))
			printf(out_fmt, rlist->line);
		else if (match_mode == MATCH_MODE_PREFIX)
			printf(out_fmt, has_match ? 1 : 0, rlist->line);
	}

	free(fp1_fmt);
	free(fp2_fmt);

	fclose(fp1);
	fclose(fp2);
	while (rlist != NULL)
	{
		rlist_end = rlist;
		rlist = rlist->next;
		free(rlist_end);
	}

	return 0;
}

