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
  string walksfilename = filename+".walks";
  string outputfilename = filename+".fasta";
  string largefilename = filename+".onlylarge.fasta";
  FILE *foc = fopen(walksfilename.c_str(), "w");
  fclose(foc);
  FILE *fo = fopen(outputfilename.c_str(), "w");
  fclose(fo);
  for (int i = 0; i < paths.size(); i++) {
    printf("(");
    for (int j = 0; j < paths[i].size(); j++) {
      printf("%d(%d)%c", paths[i][j], paths[i][j] >= 0 ? gr.nodes[paths[i][j]]->s.size() : 0, j + 1 == paths[i].size() ? ')' : ',');
    }
    printf(" ");
    gr.OutputPathA(paths[i], kmer, outputfilename.c_str(), i);
    gr.OutputPathC(paths[i], kmer, walksfilename.c_str(), i);
  }
  printf("onlylarge\n");
  fo = fopen((filename+".onlylarge.fasta").c_str(), "w");
  fclose(fo);
  for (int i = 0; i < paths.size(); i++) {
    gr.OutputPathAT(paths[i], kmer, largefilename.c_str(), i, threshold);
  }
}
