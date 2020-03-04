#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <ctype.h>
#include <sys/types.h>

#define VERSION "0.3.0"
#define EXENAME "cgmlst-dists"
#define GITHUB_URL "https://github.com/tseemann/cgmlst-dists"
//#define DEBUG

const int MAX_LINE = 1E5;
const int MAX_ASM  = 1E5;
const char* DELIMS = "\n\r\t ";
const int IGNORE_ALLELE = 0;


//------------------------------------------------------------------------
void show_help(int retcode)
{
  FILE* out = (retcode == EXIT_SUCCESS ? stdout : stderr);

  static const char str[] = {
      "SYNOPSIS\n  Pairwise CG-MLST distance matrix from allele call tables\n"
      "USAGE\n  %s [options] chewbbaca.tab > distances.tsv\n"
      "OPTIONS\n"
      "  -h\tShow this help\n"
      "  -v\tPrint version and exit\n"
      "  -q\tQuiet mode; do not print progress information\n"
      "  -c\tUse comma instead of tab in output\n"
      "  -m N\tOutput: 1=lower-tri 2=upper-tri 3=full [3]\n"
//      "  -t N\tNumber of threads to use [1]\n"
      "URL\n  %s\n"};
  fprintf(out, str, EXENAME, GITHUB_URL);
  exit(retcode);
}

//------------------------------------------------------------------------
int distance(const int* restrict a, const int* restrict b, size_t len)
{
  int diff=0;
  for (size_t i=0; i < len; i++) {
    if (a[i] != b[i] && a[i] != IGNORE_ALLELE && b[i] != IGNORE_ALLELE) {
      diff++;
    }
  }
  return diff;
}

//------------------------------------------------------------------------
void* calloc_safe(size_t nmemb, size_t size) 
{
  void* ptr = calloc(nmemb, size);
  if (ptr == NULL) {
    fprintf(stderr, "ERROR: could not allocate %ld kb RAM\n", (nmemb*size)>>10);
    exit(EXIT_FAILURE);
  }
  return ptr;
}

//------------------------------------------------------------------------

void cleanup_line(char* str)
{
  char* s = str;
#ifdef DEBUG
  fprintf(stderr, "BEFORE: %s", str);
#endif
  // skip over first column (ID)
  while (*s != 0 && *s != '\t') s++;

  // replace alpha with space so atoi() works
  while (*s++) {
    if (isalpha(*s)) *s = ' ';
  }
#ifdef DEBUG
  fprintf(stderr, "AFTER : %s", str);
#endif
}

//------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  // parse command line parameters
  int opt, quiet = 0, csv = 0, threads = 1, mode = 3;
  while ((opt = getopt(argc, argv, "hqcvm:t:")) != -1) {
    switch (opt) {
      case 'h': show_help(EXIT_SUCCESS); break;
      case 'q': quiet = 1; break;
      case 'c': csv = 1; break;
      case 't': threads = atoi(optarg); break;
      case 'm': mode = atoi(optarg); break;
      case 'v': printf("%s %s\n", EXENAME, VERSION); exit(EXIT_SUCCESS);
      default: show_help(EXIT_FAILURE);
    }
  }

  // require a filename argument
  if (optind >= argc) {
    show_help(EXIT_FAILURE);
    return 0;
  }
  const char* infile = argv[optind];

  // parameter sanity check
  if (threads < 1) threads = 1;

  // say hello
  if (!quiet) {
    fprintf(stderr, "This is %s %s\n", EXENAME, VERSION);
    // fprintf(stderr, "Using %d threads (threads>1 currently unsupported).\n", threads);
  }

  // read file one line at a time
  FILE* in = (FILE*) fopen(infile, "r");
  if (! in) {
    fprintf(stderr, "ERROR: can not open file '%s'\n", infile);
    exit(EXIT_FAILURE);
  }

  // allocate RAM
  char* buf = (char*) calloc_safe( MAX_LINE, sizeof(char) );

  char** id  = (char**) calloc_safe( MAX_ASM, sizeof(char*) );
  int** call = (int**) calloc_safe( MAX_ASM, sizeof(int*) );

  int row = -1;
  int ncol = 0;
   
  while (fgets(buf, MAX_LINE, in))
  {
    // cleanup non-numerics in NON-HEADER lines
    if (row >=0) cleanup_line(buf);

    // scan for tab separated values
    char* save;
    char* s = strtok_r(buf, DELIMS, &save);
    int col = -1;
    while (s) {
      //fprintf(stderr, "DEBUG: row=%d col=%d s='%s'\n", row, col, s);
      if (row >= 0) {
        if (col < 0) {
          id[row] = strdup(s);
          call[row] = (int*) calloc_safe(ncol, sizeof(int*));
        }
        else {
          call[row][col] = atoi(s);
        }
      }
      else {
        // just parsing columns on first row
      }
      col++;
      s = strtok_r(NULL, DELIMS, &save);
    }
    row++;
//    if (!quiet) fprintf(stderr, "row %d has %d cols\n", row, col) ;
    if (!quiet) fprintf(stderr, "\rLoading row %d", row);
    if (row==0) ncol = col;
  }
  int nrow = row;
  fclose(in);
  
  // what we collected
  if (!quiet) fprintf(stderr, "\rLoaded %d samples x %d allele calls\n", nrow, ncol);

  // build an output matrix (one dimensional j*nrow+i access)
  int* dist = calloc_safe(nrow*nrow, sizeof(int));
  
  for (int j=0; j < nrow; j++) {
    if (!quiet) fprintf(stderr, "\rCalculating distances: %.2f%%", (j+1)*100.0/nrow);
    for (int i=0; i < j; i++) {
      int d = distance(call[j], call[i], ncol);
      dist[j*nrow+i] = dist[i*nrow+j] = d;  // matrix is diagonal symetric
    }
  }
  if (!quiet) fprintf(stderr, "\nWriting distance matrix to stdout...\n");

  // separator choice
 
  // Print header row
  char sep = csv ? ',' : '\t';
  for (int j=0; j < nrow; j++) {
    if (j==0) printf(EXENAME);
    printf("%c%s", sep, id[j]);
  }
  printf("\n");

  // Print matrix
  for (int j=0; j < nrow; j++) {
    printf("%s", id[j]);
    int start = (mode & 1) ?    0 : j   ;  // upper?
    int end   = (mode & 2) ? nrow : j+1 ;  // lower?
    for (int i=start; i < end; i++) {
      printf("%c%d", sep, dist[j*nrow + i]);
    }
    printf("\n");
  }

  // free RAM
  for (int i=0; i < nrow; i++) {
    free( id[i] );
    free( call[i] );
  }
  free(id);
  free(call);
  free(dist);
  free(buf);

  if (!quiet) fprintf(stderr, "\nDone.\n");

  return 0;
}

//------------------------------------------------------------------------

