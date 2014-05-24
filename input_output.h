#ifndef INPUT_OUTPUT_H__
#define INPUT_OUTPUT_H__

#include "graph.h"

void OutputPathsToConsole(const vector<vector<int>>& paths, Graph& gr, int threshold);
void OutputPathsToFile(const vector<vector<int>>& paths, Graph& gr, int kmer, int threshold, const string& filename);

#endif 
