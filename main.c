#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <ctype.h>
#include <sys/types.h>

#define VERSION "0.4.0"
#define EXENAME "cgmlst-dists"
#define GITHUB_URL "https://github.com/tseemann/cgmlst-dists"
#define DEBUG

const int MAX_LINE = 1E5;
const int MAX_ASM  = 1E5;
const char* DELIMS = "\n\r\t";
const int IGNORE_ALLELE = 0;
const char REPLACE_CHAR = ' ';

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
      "  -x N\tStop calculating beyond this distance [9999]\n"
      "  -H\tAdd this flag if chewbbaca was run with the --hash-profiles parameter\n"
//      "  -t N\tNumber of threads to use [1]\n"
      "URL\n  %s\n"};
  fprintf(out, str, EXENAME, GITHUB_URL);
  exit(retcode);
}

//------------------------------------------------------------------------
int distance(const long long* restrict a, const long long* restrict b, size_t len, int maxdiff)
{
  int diff=0;
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

int str_replace(char* str, char* old, char* new)
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

#include <stdio.h>
#include <string.h>

long long sha1_to_int(char *sha1_str) {
    // This overflows which increases the chance of a collision.
    // This is highly unlikely this context though.
    int i = 0;
    long long value = 0;
    char *p;

    // convert each hex digit to its integer value and accumulate
    for (i = 0, p = sha1_str; i < 40; i++, p++) {
        int digit;
        if (*p >= '0' && *p <= '9') {
          digit = *p - '0';
        } else if (*p >= 'a' && *p <= 'f') {
          digit = *p - 'a' + 10;
        } else if (*p >= 'A' && *p <= 'F') {
          digit = *p - 'A' + 10;
        } else {
          // If anything other than a hex digit is found the end of the hash is reached
          return value;
        }
        value = (value << 4) | digit;
    }
    return value;
}


//------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  // parse command line parameters
  int opt, quiet = 0, csv = 0, threads = 1, mode = 3, maxdiff = 9999, use_hashed_profiles = 0;
  while ((opt = getopt(argc, argv, "hqcvHm:t:x:")) != -1) {
    switch (opt) {
      case 'h': show_help(EXIT_SUCCESS); break;
      case 'q': quiet = 1; break;
      case 'c': csv = 1; break;
      case 'H': use_hashed_profiles = 1; break;
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
  long long** call = (long long**) calloc_safe( MAX_ASM, sizeof(long long*) );

  int row = -1;
  int ncol = 0;
   
  while (fgets(buf, MAX_LINE, in))
  {
    // cleanup non-numerics in NON-HEADER lines
    // No need to do this when hashed chewbbaca output is used.
    // '-' and 'NA' get evaluated to 0 by sha1_to_int.
    if (!use_hashed_profiles && row >=0){
      cleanup_line(buf);
    }
    
    // scan for tab separated values
    char* save;
    char* s = strtok_r(buf, DELIMS, &save);
    int col = -1;
    while (s) {
      //fprintf(stderr, "DEBUG: row=%d col=%d s='%s'\n", row, col, s);
      if (row >= 0) {
        if (col < 0) {
          if (strlen(s)==0) {
            fprintf(stderr, "row %d has an empty ID in first column\n", row+1);
            exit(EXIT_FAILURE);
          }
          id[row] = strdup(s);
          call[row] = (long long*) calloc_safe(ncol, sizeof(long long*));
        }
        else {
          // INF-xxxx are returned as -ve numbers
          if (use_hashed_profiles) {
            call[row][col] = sha1_to_int(s) ;
          }
          else{
            call[row][col] = abs( atoi(s) );
          }
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
  int nrow = row;
  fclose(in);
  
  // what we collected
  if (!quiet) fprintf(stderr, "\rLoaded %d samples x %d allele calls\n", nrow, ncol);

  // build an output matrix (one dimensional j*nrow+i access)
  int* dist = calloc_safe(nrow*nrow, sizeof(int));
  
  for (int j=0; j < nrow; j++) {
    if (!quiet) fprintf(stderr, "\rCalculating distances: %.2f%%", (j+1)*100.0/nrow);
    for (int i=0; i < j; i++) {
      int d = distance(call[j], call[i], ncol, maxdiff);
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
