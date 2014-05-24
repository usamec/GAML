#include "input_output.h"

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

void OutputPathsToConsole(const vector<vector<int>>& paths, Graph& gr, int threshold) {
  for (int i = 0; i < paths.size(); i++) {
    printf("(");
    for (int j = 0; j < paths[i].size(); j++) {
      if (paths[i][j] >= 0 && gr.nodes[paths[i][j]]->s.length() > threshold) printf(ANSI_COLOR_GREEN);
      printf("%d" ANSI_COLOR_RESET "%c", paths[i][j], j + 1 == paths[i].size() ? ')' : ',');
    }
    printf(" ");
  }
}

void OutputPathsToFile(const vector<vector<int>>& paths, Graph& gr, int kmer, int threshold, const string& filename) {
  FILE *foc = fopen((filename+".walks").c_str(), "w");
  fclose(foc);
  FILE *fo = fopen((filename+".fasta").c_str(), "w");
  fclose(fo);
  for (int i = 0; i < paths.size(); i++) {
    printf("(");
    for (int j = 0; j < paths[i].size(); j++) {
      printf("%d%c", paths[i][j], j + 1 == paths[i].size() ? ')' : ',');
    }
    printf(" ");
    gr.OutputPathA(paths[i], kmer, filename.c_str(), i);
    gr.OutputPathC(paths[i], kmer, ("c"+filename).c_str(), i);
  }
  fo = fopen((filename+".onlylarge.fasta").c_str(), "w");
  fclose(fo);
  for (int i = 0; i < paths.size(); i++) {
    gr.OutputPathAT(paths[i], kmer, ("s"+filename).c_str(), i, threshold);
  }
}
