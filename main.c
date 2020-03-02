#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <sys/stat.h>
#include <sys/types.h>

#define VERSION "0.1.0"
#define EXENAME "cgmlst-dists"
#define GITHUB_URL "https://github.com/tseemann/cgmlst-dists"

const int MAX_LINE = 1E5;
const int MAX_ASM  = 1E5;
const char* DELIMS = "\n\r\t ";
const int IGNORE_ALLELE = 0;

//------------------------------------------------------------------------
int distance(const int* restrict a, const int* restrict b, size_t L)
{
  int diff=0;
  for (size_t i=0; i < L; i++) {
    if (a[i] != b[i] && a[i] != IGNORE_ALLELE && b[i] != IGNORE_ALLELE) {
      diff++;
    }
  }
  return diff;
}

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
      "URL\n  %s\n"};
  fprintf(out, str, EXENAME, GITHUB_URL);
  exit(retcode);
}

//------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  // parse command line parameters
  int opt, quiet = 0, csv = 0;
  while ((opt = getopt(argc, argv, "hqcv")) != -1) {
    switch (opt) {
      case 'h': show_help(EXIT_SUCCESS); break;
      case 'q': quiet = 1; break;
      case 'c': csv = 1; break;
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

  // say hello
  if (!quiet)
    fprintf(stderr, "This is %s %s\n", EXENAME, VERSION);


  // read file one line at a time
  FILE* in = (FILE*) fopen(infile, "r");
  if (! in) {
    fprintf(stderr, "ERROR: can not open file '%s'\n", infile);
    exit(EXIT_FAILURE);
  }

  // allocate RAM
  char* buf = (char*) calloc( MAX_LINE, sizeof(char) );

  char** id  = (char**) calloc( MAX_ASM, sizeof(char*) );
  int** call = (int**) calloc( MAX_ASM, sizeof(int*) );

  int row = -1;
  int ncol = 0;
   
  while (fgets(buf, MAX_LINE, in)) {
    char* save;
    char* s = strtok_r(buf, DELIMS, &save);
    int col = -1;
    while (s) {
      //fprintf(stderr, "DEBUG: row=%d col=%d s='%s'\n", row, col, s);
      if (row >= 0) {
        if (col < 0) {
          id[row] = strdup(s);
          call[row] = (int*) calloc(ncol, sizeof(int*));
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
  
  // Print header row
  char sep = csv ? ',' : '\t';
  for (int j=0; j < nrow; j++) {
    if (j==0) printf(EXENAME);
    printf("%c%s", sep, id[j]);
  }
  printf("\n");
    
  // Calculate and print distances as we go, one row at a time
  // FIXME: fill a matrix and halve the work, and parallelize this
  for (int j=0; j < nrow; j++) {
    if (!quiet) fprintf(stderr, "\rCalulating distances... %.2f%%", j*100.0/(nrow-1));
    printf("%s", id[j]);
    for (int i=0; i < nrow; i++) {
      int d = distance(call[j], call[i], ncol);
      printf("%c%d", sep, d);
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
  free(buf);

  if (!quiet) fprintf(stderr, "\nDone.\n");

  return 0;
}

//------------------------------------------------------------------------

