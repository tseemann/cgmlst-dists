#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <ctype.h>
#include <sys/types.h>
#include <inttypes.h>

#define VERSION "0.4.0"
#define EXENAME "cgmlst-dists"
#define GITHUB_URL "https://github.com/genpat-it/cgmlst-dists"
//#define DEBUG

const int32_t MAX_LINE = 1E5;
const int32_t MAX_ASM  = 1E5;
const char* DELIMS = "\n\r\t";
const uint8_t IGNORE_ALLELE = 0;
const char REPLACE_CHAR = ' ';

//------------------------------------------------------------------------
void show_help(uint32_t retcode)
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
      "  -x N\tStop calculating beyond this distance [9999]\n"
//      "  -t N\tNumber of threads to use [1]\n"
      "URL\n  %s\n"};
  fprintf(out, str, EXENAME, GITHUB_URL);
  exit(retcode);
}

//------------------------------------------------------------------------
uint32_t distance(const uint32_t* restrict a, const uint32_t* restrict b, size_t len, uint32_t maxdiff)
{
  uint32_t diff=0;
  for (size_t i=0; i < len; i++) {
    if (a[i] != b[i] && a[i] != IGNORE_ALLELE && b[i] != IGNORE_ALLELE) {
      diff++;
      if (diff >= maxdiff) return maxdiff;
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

uint32_t str_replace(char* str, char* old, char* new)
{
  size_t sl = strlen(str);
  size_t ol = strlen(old);
  size_t nl = strlen(new);
  if (ol < 1 || nl < 1 || ol != nl || sl < ol) {
    fprintf(stderr, "ERROR: str_replace(%lu,%lu,%lu)\n", sl, ol, nl);
    exit(EXIT_FAILURE);
  }

  // char *strstr(const char *haystack, const char *needle);
  char* p = NULL;
  while ( (p = strstr(str, old)) != NULL ) {
    // char *strncpy(char *dest, const char *src, size_t n);
    strncpy(p, new, nl);
  }
  
  return sl;
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

  // CHewbacca codes
  // LINF NIPH INF-nnn PLOT3 PLOT5 ASM
  // these two are special as they end in numbers
  // and don't want to confuse with INF-xxx
  str_replace(s, "PLOT3", "    0");
  str_replace(s, "PLOT5", "    0");

  // replace alpha with space so atoi() works
  while (*s++) {
    if (isalpha(*s)) {
      *s = REPLACE_CHAR;
    }
  }
#ifdef DEBUG
  fprintf(stderr, "AFTER : %s", str);
#endif
}

//------------------------------------------------------------------------
int32_t main(int argc, char* argv[])
{
  // parse command line parameters
  int32_t opt, quiet = 0, csv = 0, threads = 1, mode = 3, maxdiff = 9999;
  while ((opt = getopt(argc, argv, "hqcvm:t:x:")) != -1) {
    switch (opt) {
      case 'h': show_help(EXIT_SUCCESS); break;
      case 'q': quiet = 1; break;
      case 'c': csv = 1; break;
      case 't': threads = atoi(optarg); break;
      case 'm': mode = atoi(optarg); break;
      case 'x': maxdiff = atoi(optarg); break;
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
  uint32_t** call = (uint32_t**) calloc_safe( MAX_ASM, sizeof(uint32_t*) );

  int32_t row = -1;
  uint32_t ncol = 0;
   
  while (fgets(buf, MAX_LINE, in))
  {
    // cleanup non-numerics in NON-HEADER lines
    if (row >=0) cleanup_line(buf);

    // scan for tab separated values
    char* save;
    char* s = strtok_r(buf, DELIMS, &save);
    int32_t col = -1;
    while (s) {
      //fprintf(stderr, "DEBUG: row=%d col=%d s='%s'\n", row, col, s);
      if (row >= 0) {
        if (col < 0) {
          if (strlen(s)==0) {
            fprintf(stderr, "row %d has an empty ID in first column\n", row+1);
            exit(EXIT_FAILURE);
          }
          id[row] = strdup(s);
          call[row] = (uint32_t*) calloc_safe(ncol, sizeof(uint32_t*));
        }
        else {
          // INF-xxxx are returned as -ve numbers
          call[row][col] = abs( atoi(s) );
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
    if (!quiet) fprintf(stderr, "\rLoaded row %d", row);
    if (row==0) ncol = col;
    if (row >= MAX_ASM) {
      fprintf(stderr, "Too many rows, can only handle %d\n", MAX_ASM);
      exit(EXIT_FAILURE);
    }
    
    if (col != ncol) {
      fprintf(stderr, "\nERROR: row %d had %d cols, expected %d\n", row+1, col+1, ncol);
      exit(-1);
    }
  }
  int64_t nrow = row;
  fclose(in);
  
  // what we collected
  if (!quiet) fprintf(stderr, "\rLoaded %ld samples x %d allele calls\n", nrow, ncol);

  // build an output matrix (one dimensional j*nrow+i access)
  uint32_t* dist = calloc_safe(nrow*nrow, sizeof(uint32_t));
  
  for (int64_t j=0; j < nrow; j++) {
    if (!quiet) fprintf(stderr, "\rCalculating distances: %.2f%%", (j+1)*100.0/nrow);
    for (int64_t i=0; i < j; i++) {
      uint32_t d = distance(call[j], call[i], ncol, maxdiff);
      dist[j*nrow+i] = dist[i*nrow+j] = d;  // matrix is diagonal symetric
    }
  }
  if (!quiet) fprintf(stderr, "\nWriting distance matrix to stdout...\n");

  // separator choice
 
  // Print header row
  char sep = csv ? ',' : '\t';
  for (uint32_t j=0; j < nrow; j++) {
    if (j==0) printf(EXENAME);
    printf("%c%s", sep, id[j]);
  }
  printf("\n");

  // Print matrix
  for (int64_t j=0; j < nrow; j++) {
    printf("%s", id[j]);
    uint32_t start = (mode & 1) ?    0 : j   ;  // upper?
    uint32_t end   = (mode & 2) ? nrow : j+1 ;  // lower?
    for (int64_t i=start; i < end; i++) {
      printf("%c%d", sep, dist[j*nrow + i]);
    }
    printf("\n");
  }

  // free RAM
  for (uint32_t i=0; i < nrow; i++) {
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