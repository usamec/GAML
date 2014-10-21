#include "graph.h"
#include <cstdio>
#include <string>
#include <fstream>
#include <cassert>
#include <boost/algorithm/string.hpp>
#include <cstdlib>
#include <queue>
#include <deque>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <map>
#include <set>
#include <ctime>
#include <chrono>
#include <sys/timeb.h>
#include "unordered_map.hpp"
#include "utility.h"

using namespace std;
using namespace boost;

const double kThresholdProb = 1e-35;
const double kThresholdProb2 = 1e-15;
const double kSmooth = 1;
const int kMinSubpathLength = 300;
const char kThreads[] = "-p 4";
const char kThreadsBlasr[] = "-nproc 16";
const char kContigSeparator = '\n';
const int kMinAnchorLen = 80;
const int kBorderLen = 60;
const int kIndexKmer = 15;

extern string gBowtiePath;
extern string gBlasrPath;

unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
//default_random_engine generator(seed1);
default_random_engine generator(47);

int getMilliCount(){
  static int last = 0;
  timeb tb;
  ftime(&tb);
  int nCount = tb.millitm + (tb.time & 0xfffff) * 1000;
  int ret = nCount - last;
  last = nCount;
  return ret;
}

bool LoadGraph(const string& filename, Graph& gr) {
  ifstream f(filename.c_str());
  if (!f.is_open()) {
    return false;
  }

  string l;
  getline(f, l);

  vector<string> fl;
  split(fl, l, is_any_of("\t"));

  int n = StringToInt(fl[0]);
  gr.nodes.resize(2*n);

  for (int i = 0; i < n; i++) {
    // TODO: check this line
    getline(f, l);
    Node *n1 = new Node;
    Node *n2 = new Node;
    n1->id = 2*i;
    n2->id = 2*i+1;
    n1->inv = n2;
    n2->inv = n1;
    getline(f, n1->s);
    getline(f, n2->s);
    gr.nodes[2*i] = n1;
    gr.nodes[2*i+1] = n2;
  }
  int narcs = 0;
  while (getline(f, l)) {
    if (l.substr(0, 3) == "ARC") {
      vector<string> al;
      split(al, l, is_any_of("\t"));
      int source = ConvertNodeId(StringToInt(al[1]));
      int dest = ConvertNodeId(StringToInt(al[2]));
      gr[source]->next.push_back(gr[dest]);
      gr[source]->next_prob.push_back(kSmooth);
      gr[InvertNode(dest)]->next.push_back(gr[InvertNode(source)]);
      gr[InvertNode(dest)]->next_prob.push_back(kSmooth);
      narcs++;
    }
  }

  f.close();

  gr.CalcProbSums();
  gr.CalcNormalizeMap();
  printf("Loaded %d nodes %d arcs\n", n, narcs);
  return true;
}

void Graph::CalcReachabilityLimit(int max_dist) {
  printf("reach limit start\n");
  reach_limit_.resize(nodes.size());
  for (int i = 0; i < nodes.size(); i++) {
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> fr;
    fr.push(make_pair(0, i));
    vector<int> final_dist(nodes.size(), -1);
    vector<int> tmp_dist(nodes.size(), 2*max_dist);
    vector<int> prev(nodes.size(), -1);
    tmp_dist[i] = 0;
    prev[i] = -2;
    while (!fr.empty()) {
      int d = fr.top().first;
      int x = fr.top().second;
      fr.pop();
      if (final_dist[x] != -1) {
        continue;
      }
      final_dist[x] = d;
      int nd = d;
      if (x != i) {
        vector<int> pp;
        int cur = prev[x];
        while (cur != i) {
          pp.push_back(cur);
          cur = prev[cur];
        }
        reverse(pp.begin(), pp.end());
        reach_limit_[i][x] = pp;

        nd += nodes[x]->s.length();
      }
      for (int j = 0; j < nodes[x]->next.size(); j++) {
        int nx = nodes[x]->next[j]->id;
        if (tmp_dist[nx] > nd && nd <= max_dist) {
          tmp_dist[nx] = nd;
          prev[nx] = x;
          fr.push(make_pair(nd, nx));
        }
      }
    }
  }
  printf("reach limit end\n");
}

void Graph::CalcReachabilityBig(int threshold) {
  printf("reach start\n");
  reach_sets_.resize(nodes.size());
  reach_big_.resize(nodes.size());
  for (int i = 0; i < nodes.size(); i++) {
    queue<int> fr;
    vector<bool> visited(nodes.size());
    vector<int> prev(nodes.size(), -1);
    visited[i] = true;
    fr.push(i);
    while (!fr.empty()) {
      int x = fr.front();
      fr.pop();
      reach_sets_[x].insert(i);
      if (nodes[x]->s.length() > threshold && nodes[i]->s.length() > threshold && x != i) {
        vector<int> pp;
        int cur = prev[x];
        while (cur != i) {
          pp.push_back(cur);
          cur = prev[cur];
        }
        reverse(pp.begin(), pp.end());
        reach_big_[i][x] = pp;
      }
      if (nodes[x]->s.length() > threshold && x != i) continue;
      for (int j = 0; j < nodes[x]->next.size(); j++) {
        int ni = nodes[x]->next[j]->id;
        if (visited[ni]) continue;
        visited[ni] = true;
        prev[ni] = x;
        fr.push(ni);
      }
    }
  }
  printf("reach end\n");
}

void Graph::CalcReachability() {
  printf("reach start\n");
  reachability_.resize(nodes.size());
  for (int i = 0; i < nodes.size(); i++) {
    queue<int> fr;
    vector<bool> visited(nodes.size());
    visited[i] = true;
    fr.push(i);
    while (!fr.empty()) {
      int x = fr.front();
      fr.pop();
      reachability_[i]+=nodes[x]->s.length();
      for (int j = 0; j < nodes[x]->next.size(); j++) {
        int ni = nodes[x]->next[j]->id;
        if (visited[ni]) continue;
        visited[ni] = true;
        fr.push(ni);
      }
    }
  }
  for (int i = 0; i < reachability_.size(); i++)
    reachability_[i] = pow(reachability_[i], 2);
  reach_sum_ = accumulate(reachability_.begin(), reachability_.end(), 0.0);
  printf("reach end %lf\n", reach_sum_);

  printf("reach self start\n");
  int num_self_find = 0;
  reach_self_.resize(nodes.size());
  for (int i = 0; i < nodes.size(); i++) {
    vector<vector<int> > cands({vector<int>({i})});
    for (int s = 0; s < 4; s++) {
      vector<vector<int> > cands2;
      for (auto &c: cands) {
        for (int j = 0; j < nodes[c.back()]->next.size(); j++) {
          if (nodes[c.back()]->next[j]->id == i) {
            printf("add %d ", i);
            for (auto &cc: c)
              printf("%d ", cc);
            printf("\n");
            reach_self_[i].push_back(c);
            num_self_find++;
          } else {
            vector<int> p(c);
            p.push_back(nodes[c.back()]->next[j]->id);
            cands2.push_back(p);
          }
        }
      }
      cands = cands2;
    }
  }
  printf("reach self done %d\n", num_self_find);
}

int Graph::SampleVertexByReach() const {
  uniform_real_distribution<double> dist(0.0, reach_sum_);
  double r = dist(generator);
  double ss = 0;
  for (int i = 0; i < reachability_.size(); i++) {
    ss += reachability_[i];
    if (ss > r) {
      return i;
    }
  }
  return reachability_.size() - 1;
}

void Graph::OutputPath(const vector<int>& path, int kmer) {
  FILE *f = fopen("outputx.fas", "w");
  int id = 0;
  fprintf(f, ">tmp\n");
  fprintf(f, "%s", nodes[path[0]]->GetNodeSeq(kmer).c_str());
  for (int i = 1; i < path.size(); i++) {
    fprintf(f, "%s", nodes[path[i]]->s.c_str());
  }
  fprintf(f, "\n");
  fclose(f);
}

void Graph::OutputPath(const vector<int>& path, int kmer, string filename) {
  FILE *f = fopen(filename.c_str(), "w");
  int id = 0;
  fprintf(f, ">tmp\n");
  fprintf(f, "%s", nodes[path[0]]->GetNodeSeq(kmer).c_str());
  for (int i = 1; i < path.size(); i++) {
    fprintf(f, "%s", nodes[path[i]]->s.c_str());
  }
  fprintf(f, "\n");
  fclose(f);
}

void Graph::OutputPathAT(const vector<int>& path, int kmer, string filename, int cid, int threshold) {
  FILE *f = fopen(filename.c_str(), "a");
  int id = 0;
  fprintf(f, ">tmp%d\n", cid);
  fprintf(f, "%s", nodes[path[0]]->GetNodeSeq(kmer).c_str());
  for (int i = 1; i < path.size(); i++) {
    if (path[i] < 0) {
      for (int j = 0; j < -path[i]; j++) {
        fprintf(f, "N");
      }
    }
    else if (nodes[path[i]]->s.length() > threshold) {
      fprintf(f, "%s", nodes[path[i]]->s.c_str());
    } else {
      for (int j = 0; j < nodes[path[i]]->s.length(); j++) {
        fprintf(f, "N");
      }
    }
  }
  fprintf(f, "\n");
  fclose(f);
}

void Graph::OutputPathC(const vector<int>& path, int kmer, string filename, int cid) {
  FILE *f = fopen(filename.c_str(), "a");
  int id = 0;
  fprintf(f, ">tmp%d-", cid);
  int pos = 0;
  for (int i = 0; i < path.size(); i++) {
    fprintf(f, "%d(%d)%c", path[i], pos, i + 1 == path.size() ? '\n' : '-');
    if (path[i] >= 0)
      pos += nodes[path[i]]->s.length();
    else
      pos += -path[i];
  }
  fclose(f);
}

void Graph::OutputPathA(const vector<int>& path, int kmer, string filename, int cid) {
  FILE *f = fopen(filename.c_str(), "a");
  int id = 0;
  fprintf(f, ">tmp%d\n", cid);
/*  fprintf(f, ">c%d-", cid);
  int pos = 0;
  for (int i = 0; i < path.size(); i++) {
    fprintf(f, "%d(%d)%c", path[i], pos, i + 1 == path.size() ? '\n' : '-');
    pos += nodes[path[i]]->s.length();
  }*/
  fprintf(f, "%s", nodes[path[0]]->GetNodeSeq(kmer).c_str());
  for (int i = 1; i < path.size(); i++) {
    if (path[i] < 0) {
      for (int j = 0; j < -path[i]; j++) {
        fprintf(f, "N");
      }      
    } else {
      fprintf(f, "%s", nodes[path[i]]->s.c_str());
    }
  }
  fprintf(f, "\n");
  fclose(f);
}

void ReadSet::ClearPositions() {
  positions_.resize(reads_num_);
  for (int i = 0; i < positions_.size(); i++) {
    positions_[i].clear();
  }
}

vector<vector<pair<int, pair<int, int> > > >& ReadSet::GetPositionsSlow(
    const Graph& gr, const vector<int>& path, int& total_len) {
  char tmpname1[L_tmpnam+6], tmpname2[L_tmpnam], tmpname3[L_tmpnam],
       tmpname4[L_tmpnam+8];
  tmpnam(tmpname1);
  strcat(tmpname1, ".fas");
  tmpnam(tmpname2);
  tmpnam(tmpname3);
  tmpnam(tmpname4);
  strcat(tmpname4, ".fastq");
  FILE *f = fopen(tmpname1, "w");
  int id = 0;
  fprintf(f, ">tmp\n");
  string seq;
  total_len = gr.nodes[path[0]]->s.length();
  seq = gr.nodes[path[0]]->s;
  fprintf(f, "%s", gr.nodes[path[0]]->s.c_str());
  for (int i = 1; i < path.size(); i++) {
    fprintf(f, "%s", gr.nodes[path[i]]->s.c_str());
    total_len += gr.nodes[path[i]]->s.length();
    seq += gr.nodes[path[i]]->s;
  }
  fprintf(f, "\n");
  fclose(f);
  for (int i = 0; i < positions_.size(); i++) {
    positions_[i].clear();
  }
  
  string reads_filename = filename_;

  printf("slow files %s %s %s %s\n", tmpname1, tmpname2, tmpname3, tmpname4);
  unordered_set<int> read_cands;
  read_index_.GetReadCands(seq, read_cands);

  ofstream of(tmpname4, ios_base::out | ios_base::trunc);
  for (auto &e: read_cands) {
    of << "@" << read_map_inv_[e] << endl;
    of << read_seqs_[e] << endl;
    of << "+" << endl;
    of << read_seqs_[e] << endl;
  }
  of.close();
  reads_filename = tmpname4;

  string cmd = gBowtiePath + "/bowtie2-build ";
  cmd += tmpname1;
  cmd += " ";
  cmd += tmpname2;
  cmd += ">/dev/null 2>/dev/null";
  system(cmd.c_str());
  
  cmd = gBowtiePath + "/bowtie2 -a -x ";
  cmd += tmpname2;
  cmd += " -q -U ";
  cmd += reads_filename;
  cmd += " --very-sensitive --sam-no-hd --no-unal -k 10000 ";
  cmd += " --ignore-quals ";
  cmd += kThreads;
  cmd += " -S ";
  cmd += tmpname3;
  cmd += ">/dev/null 2>/dev/null";
  system(cmd.c_str());

  ifstream fi(tmpname3);
  string l;
  while (getline(fi, l)) {
    vector<string> parts;
    split(parts, l, is_any_of("\t"));
    int read_id = GetReadId(parts[0]);
    if (StringToInt(parts[1]) & 0x4) {
      continue;
    }
    int orientation = 0;
    if (StringToInt(parts[1]) & 0x10) {
      orientation = 1;
    }
    int edit_dist = -1;
    for (int i = 11; i < parts.size(); i++) {
      if (parts[i].substr(0, 5) == "NM:i:") {
        edit_dist = StringToInt(parts[i].substr(5));
        break;
      }
    }
    if (edit_dist != -1) {
      int pos = StringToInt(parts[3]);
      if (read_id >= positions_.size()) {
        positions_.resize(read_id+1);
      }
      positions_[read_id].push_back(make_pair(pos, make_pair(edit_dist, orientation)));
    }
  }

  remove(tmpname1);
  remove(tmpname2);
  remove(tmpname3);

  return positions_;
}

vector<vector<pair<int, pair<int, int> > > >& ReadSet::GetPositions() {
  return positions_;
}

void ReadSet::PrecomputeAlignmentForPaths(const vector<vector<int>>& paths, const Graph& gr) {
  unordered_set<vector<int>> subpaths_precomp;
  for (auto &path: paths) {
    for (int i = 0; i < path.size(); i++) {
      if (path[i] < 0) continue;
      int cur_len = 0;
      cur_len = gr.nodes[path[i]]->s.length();
      int cur_seq_len = 0;
      vector<int> cur_seq, cur_seq2({-1});
      cur_seq.push_back(path[i]);
      cur_seq2.push_back(path[i]);
      for (int j = i+1; j < path.size(); j++) {
        if (path[j] < 0) break;
        cur_seq_len += gr.nodes[path[j]]->s.length();
        cur_seq.push_back(path[j]);
        cur_seq2.push_back(path[j]);
        if (cur_seq_len > kMinSubpathLength) {
          break;
        }
      }
      if (aligment_cache_.count(cur_seq) == 0) {
        subpaths_precomp.insert(cur_seq);
        subpaths_precomp.insert(InvertPath(cur_seq));
      }
  //    if (aligment_cache_.count(cur_seq2) == 0) {
  //      subpaths_precomp.insert(cur_seq2);
  //    }
    }
  }
  if (!subpaths_precomp.empty()) {
    printf("mass precomp start\n");
    PrecomputeAligmentForSubpaths(gr, USetToVector(subpaths_precomp));
  }
}

void ReadSet::GetSubpathsFromPath(
    const vector<int>& path, const Graph& gr, unordered_set<vector<int>>& subpaths_precomp) {
  for (int i = 0; i < path.size(); i++) {
    int cur_len = 0;
    cur_len = gr.nodes[path[i]]->s.length();
    int cur_seq_len = 0;
    vector<int> cur_seq, cur_seq2({-1});
    cur_seq.push_back(path[i]);
    cur_seq2.push_back(path[i]);
    for (int j = i+1; j < path.size(); j++) {
      cur_seq_len += gr.nodes[path[j]]->s.length();
      cur_seq.push_back(path[j]);
      cur_seq2.push_back(path[j]);
      if (cur_seq_len > kMinSubpathLength) {
        break;
      }
    }
    if (aligment_cache_.count(cur_seq) == 0) {
      subpaths_precomp.insert(cur_seq);
    }
//    if (aligment_cache_.count(cur_seq2) == 0) {
//      subpaths_precomp.insert(cur_seq2);
//    }
  }
}

vector<vector<pair<int, pair<int, int> > > >& ReadSet::AddPositions(
    const Graph& gr, const vector<int>& path,
    int& total_len, int st) {
//  printf("calc score\n");
  // Precomputation at once
  unordered_set<vector<int> > subpaths_precomp;
  GetSubpathsFromPath(path, gr, subpaths_precomp);
  if (!subpaths_precomp.empty()) {
    PrecomputeAligmentForSubpaths(gr, USetToVector(subpaths_precomp));
  }
 
  int cur_pos = st;
//  total_len = 0;
  for (int i = 0; i < path.size(); i++) {
    int cur_len = 0;
    cur_len = gr.nodes[path[i]]->s.length();
    total_len += cur_len;

    int cur_seq_len = 0;
    vector<int> cur_seq;
    cur_seq.push_back(path[i]);
    for (int j = i+1; j < path.size(); j++) {
      cur_seq_len += gr.nodes[path[j]]->s.length();
      cur_seq.push_back(path[j]);
      if (cur_seq_len > kMinSubpathLength) {
        break;
      }
    }

    const vector<Aligment>& align = GetAligmentForSubpath(
        gr, cur_seq);

    for (auto& al: align) {
      bool found = false;
      // TODO: optimize this (maybe)
      for (int j = 0; j < positions_[al.read_id].size(); j++) {
        if (positions_[al.read_id][j].first == al.position + cur_pos) {
          positions_[al.read_id][j].second = make_pair(al.edit_dist, al.orientation);
          found = true;
          break;
        }
      }
      if (found) continue;
      positions_[al.read_id].push_back(
          make_pair(al.position + cur_pos, make_pair(al.edit_dist, al.orientation)));
    }
    cur_pos += gr.nodes[path[i]]->s.length();
  }
  return positions_;
}

vector<vector<pair<int, pair<int, int> > > >& ReadSet::GetPositions(
    const Graph& gr, const vector<int>& path,
    int& total_len) {
  positions_.resize(reads_num_);
  for (int i = 0; i < reads_num_; i++)
    positions_[i].clear();
//  printf("calc score\n");
  // Precomputation at once
  unordered_set<vector<int> > subpaths_precomp;
  GetSubpathsFromPath(path, gr, subpaths_precomp);
  if (!subpaths_precomp.empty()) {
    PrecomputeAligmentForSubpaths(gr, USetToVector(subpaths_precomp));
  }
 
  int cur_pos = 0;
  total_len = 0;
  for (int i = 0; i < path.size(); i++) {
    if (path[i] < 0) {
      cur_pos += -path[i];
      continue;
    }
    int cur_len = 0;
    cur_len = gr.nodes[path[i]]->s.length();
    total_len += cur_len;

    int cur_seq_len = 0;
    vector<int> cur_seq;
    cur_seq.push_back(path[i]);
    for (int j = i+1; j < path.size(); j++) {
      if (path[j] < 0)
        break;
      cur_seq_len += gr.nodes[path[j]]->s.length();
      cur_seq.push_back(path[j]);
      if (cur_seq_len > kMinSubpathLength) {
        break;
      }
    }

    const vector<Aligment>& align = GetAligmentForSubpath(
        gr, cur_seq);

    for (auto& al: align) {
      bool found = false;
      // TODO: optimize this (maybe)
      for (int j = 0; j < positions_[al.read_id].size(); j++) {
        if (positions_[al.read_id][j].first == al.position + cur_pos) {
          positions_[al.read_id][j].second = make_pair(al.edit_dist, al.orientation);
          found = true;
          break;
        }
      }
      if (found) continue;
      positions_[al.read_id].push_back(
          make_pair(al.position + cur_pos, make_pair(al.edit_dist, al.orientation)));
    }
    cur_pos += gr.nodes[path[i]]->s.length();
  }
  return positions_;
}

inline void PushIfNotVisited(
    int dist, int cur_genome_pos, int cur_read_pos,
    int read_pos, int genome_pos, int iteration,
    deque<pair<int, pair<int, int>>>&fr, vector<vector<int>>& visited) {
  int gp = cur_genome_pos - genome_pos + read_pos + 20;
  if (visited[cur_read_pos + 1][gp] != iteration) {
    fr.push_back(make_pair(dist, make_pair(cur_genome_pos, cur_read_pos)));
    visited[cur_read_pos + 1][gp] = iteration;
  }
}

inline void PushFrontIfNotVisited(
    int dist, int cur_genome_pos, int cur_read_pos,
    int read_pos, int genome_pos, int iteration,
    deque<pair<int, pair<int, int>>>&fr, vector<vector<int>>& visited) {
  int gp = cur_genome_pos - genome_pos + read_pos + 20;
  if (visited[cur_read_pos + 1][gp] != iteration) {
    fr.push_front(make_pair(dist, make_pair(cur_genome_pos, cur_read_pos)));
    visited[cur_read_pos + 1][gp] = iteration;
  }
}

// Errors, genome begin, genome end
pair<int, pair<int, int>> ProcessHit(int genome_pos, int read_pos, const string& read, const string& genome) {
  static deque<pair<int, pair<int, int>>> fr;
  static int iteration = 0;
  iteration++;
  static vector<vector<int>> visited(read.size() + 47, vector<int>(read.size() + 47));
  assert(read.substr(read_pos, kIndexKmer) == genome.substr(genome_pos, kIndexKmer));
  int error_limit = 6;
  // Forward
  int forward_errs = -1;
  fr.push_back(make_pair(0, make_pair(genome_pos + kIndexKmer, read_pos + kIndexKmer)));
  int end_pos = -1;
  while (!fr.empty()) {
    pair<int, pair<int, int>> x = fr.front();
    fr.pop_front();
    if (x.first > error_limit) { 
      fr.clear();
      break;
    }
    if (x.second.second == read.size()) {
      forward_errs = x.first;
      fr.clear();
      end_pos = x.second.first - 1;
      break;
    }
    if (genome[x.second.first] == read[x.second.second]) {
      if (x.second.first + 1 < genome.size() || x.second.second + 1 == read.size()) {
        PushFrontIfNotVisited(x.first, x.second.first + 1, x.second.second + 1,
                              read_pos, genome_pos, iteration, fr, visited);
      }
    } else {
      if (x.second.first + 1 < genome.size()) {
        PushIfNotVisited(x.first + 1, x.second.first + 1, x.second.second + 1,
                         read_pos, genome_pos, iteration, fr, visited);
        PushIfNotVisited(x.first + 1, x.second.first + 1, x.second.second,
                         read_pos, genome_pos, iteration, fr, visited);
      }
      PushIfNotVisited(x.first + 1, x.second.first, x.second.second + 1,
                       read_pos, genome_pos, iteration, fr, visited);
    }
  }
  if (forward_errs == -1) return make_pair(-1, make_pair(-1, -1));
  // Backward
  int backward_errs = -1;
  int begin_pos = -1;
  fr.push_back(make_pair(0, make_pair(genome_pos - 1, read_pos - 1)));
  while (!fr.empty()) {
    pair<int, pair<int, int>> x = fr.front();
    fr.pop_front();
    if (x.first > error_limit) {
      fr.clear();
      break;
    }
    if (x.second.second == -1) {
      backward_errs = x.first;
      begin_pos = x.second.first + 1;
      fr.clear();
      break;
    }
    if (genome[x.second.first] == read[x.second.second]) {
      if (x.second.first - 1 >= 0 || x.second.second - 1 == -1) {
        PushFrontIfNotVisited(x.first, x.second.first - 1, x.second.second - 1,
                              read_pos, genome_pos, iteration, fr, visited);
      }
    } else {
      if (x.second.first - 1 >= 0) {
        PushIfNotVisited(x.first + 1, x.second.first - 1, x.second.second - 1,
                         read_pos, genome_pos, iteration, fr, visited);
        PushIfNotVisited(x.first + 1, x.second.first - 1, x.second.second,
                         read_pos, genome_pos, iteration, fr, visited);
      }
      PushIfNotVisited(x.first + 1, x.second.first, x.second.second - 1,
                       read_pos, genome_pos, iteration, fr, visited);
    }
  }
  if (backward_errs == -1) return make_pair(-1, make_pair(-1, -1));
  return make_pair(backward_errs + forward_errs, make_pair(begin_pos, end_pos));
}


void ReadSet::AlignSubpathsInternal(
    const Graph& gr, const vector<vector<int>>& subpaths) {
  for (auto &path: subpaths) {
    set<Aligment> current;
/*    for (auto &e: golden) {
      printf("g %d %d %d %d\n", e.position, e.read_id, e.orientation, e.edit_dist);
    }*/
    unordered_map<int, vector<int>> read_cands;
    string seq;
    for (auto &p: path) {
      seq += gr.nodes[p]->s;
    }
    read_index_.GetReadCandsWithPoses(seq, read_cands);
    for (auto &e: read_cands) {
      for (auto &e2: e.second) {
        int genome_pos;
        string read_seq;
//        printf("e2 %d\n", e2);
        if (e2 > 0) {
          genome_pos = e2 - kIndexKmer + 1;
          read_seq = read_seqs_[e.first];
        } else {
          genome_pos = seq.size() - (-e2 + 1);
          read_seq = ReverseSeq(read_seqs_[e.first]);
        }
        int read_pos = -1;
        for (int i = 0; i + kIndexKmer - 1 < read_seq.length(); i++) {
          if (read_seq.substr(i, kIndexKmer) == seq.substr(genome_pos, kIndexKmer)) {
            read_pos = i;
            break;
          }
        }
        if (read_pos == -1) {
          printf("%d\n%s\n%s\n", genome_pos, seq.substr(max(0,genome_pos-20), kIndexKmer+40).c_str(),
              read_seq.c_str());
        }
        assert(read_pos != -1);
        pair<int, pair<int, int>> align_res = ProcessHit(genome_pos, read_pos, read_seq, seq);
//        printf("%d %d %d\n", align_res.first, align_res.second.first, align_res.second.second);

        if (align_res.first != -1) {
          Aligment al(align_res.second.first + 1, align_res.first, e.first, e2 > 0 ? 0 : 1);
          current.insert(al);
        }
      }
    }
    for (auto &e: current) {
      aligment_cache_[path].push_back(e);
    }
  }
}

void ReadSet::PrecomputeAligmentForSubpaths(
    const Graph& gr, const vector<vector<int> >& subpaths) {
  if (subpaths.empty()) return;
  for (auto &subpath: subpaths) {
    aligment_cache_[subpath] = vector<Aligment>();
  }

  if (!external_aligner_) {
    AlignSubpathsInternal(gr, subpaths);
    return;
  }

  char tmpname1[L_tmpnam+6], tmpname2[L_tmpnam], tmpname3[L_tmpnam],
       tmpname4[L_tmpnam+8];
  tmpnam(tmpname1);
  strcat(tmpname1, ".fas");
  tmpnam(tmpname2);
  tmpnam(tmpname3);
  tmpnam(tmpname4);
  strcat(tmpname4, ".fastq");
  printf("precomp %s %s %s %s; %d subpaths\n", tmpname1, tmpname2, tmpname3, tmpname4, 
         (int)subpaths.size());
  for (int i = 0; i < subpaths.size() && i < 2; i++) {
    for (int j = 0; j < subpaths[i].size(); j++)
      printf("%d ", subpaths[i][j]);
    printf("\n");
  }
  unordered_set<int> read_cands;
  FILE *f = fopen(tmpname1, "w");
  for (int i = 0; i < subpaths.size(); i++) {
    string seq;
    fprintf(f, ">");
    for (int j = 0; j < subpaths[i].size(); j++) {
      fprintf(f, "%d%c", subpaths[i][j], j + 1 == subpaths[i].size() ? '\n' : ';');
    }
    int start = 1;
    fprintf(f, "%s", gr.nodes[subpaths[i][0]]->s.c_str());
    seq += gr.nodes[subpaths[i][0]]->s;
    for (int j = start; j < subpaths[i].size(); j++) {
      fprintf(f, "%s", gr.nodes[subpaths[i][j]]->s.c_str());
      seq += gr.nodes[subpaths[i][j]]->s;
    }
    fprintf(f, "\n");
    read_index_.GetReadCands(seq, read_cands);
  }
  printf("rc size %d\n", read_cands.size());
  fclose(f);
  ofstream of(tmpname4, ios_base::out | ios_base::trunc);
  for (auto &e: read_cands) {
    of << "@" << read_map_inv_[e] << endl;
    of << read_seqs_[e] << endl;
    of << "+" << endl;
    of << read_seqs_[e] << endl;
  }
  of.close();

  string cmd = gBowtiePath + "/bowtie2-build ";
  cmd += tmpname1;
  cmd += " ";
  cmd += tmpname2;
  cmd += ">/dev/null 2>/dev/null";
  system(cmd.c_str());
  
  cmd = gBowtiePath + "/bowtie2 -a -x ";
  cmd += tmpname2;
  cmd += " -q -U ";
  cmd += tmpname4;
  cmd += " --very-sensitive --sam-no-hd -k 10000 --reorder ";
  cmd += " --no-unal --ignore-quals ";
  if (read_cands.size() > 500) 
    cmd += kThreads;
  cmd += " -S ";
  cmd += tmpname3;
  cmd += ">/dev/null 2>/dev/null";
  system(cmd.c_str());
 
  ifstream fi(tmpname3);
  assert(fi);
  string l;
  int n_als = 0;
  while (getline(fi, l)) {
    vector<string> parts;
    split(parts, l, is_any_of("\t"));
    int read_id = GetReadId(parts[0]);
    if (StringToInt(parts[1]) & 0x4) {
      continue;
    }
    int orientation = 0;
    if (StringToInt(parts[1]) & 0x10) {
      orientation = 1;
    }
    int edit_dist = -1;
    for (int i = 11; i < parts.size(); i++) {
      if (parts[i].substr(0, 5) == "NM:i:") {
        edit_dist = StringToInt(parts[i].substr(5));
        break;
      }
    }
    string &subpath_str = parts[2];
    vector<string> subpath_parts;
    split(subpath_parts, subpath_str, is_any_of(";"));
    vector<int> subpath(subpath_parts.size());
    transform(subpath_parts.begin(), subpath_parts.end(), subpath.begin(), StringToInt);
    if (edit_dist != -1) {
      int pos = StringToInt(parts[3]);
      aligment_cache_[subpath].push_back(Aligment(pos, edit_dist, read_id, orientation));
    }
    n_als++;
  }
  for (auto& s: subpaths) {
    sort(aligment_cache_[s].begin(), aligment_cache_[s].end());
  }

  printf("precomp done %d\n", n_als);
  SaveAligments();
  remove(tmpname1);
  remove(tmpname2);
  remove(tmpname3);
}

void ReadSet::SaveAligments(bool force) {
  return;
  save_changes_++;
  printf("save ch %d\n", save_changes_);
  if (save_changes_ == 50 || force) {
    ofstream ofs(name_);
    boost::archive::binary_oarchive oa(ofs);
    oa << aligment_cache_;
    oa << read_lens_;
    oa << reads_num_;
    oa << read_map_;
    oa << read_map_inv_;
    printf("saved %d cached members\n", (int)aligment_cache_.size());
    save_changes_ = 0;
  }
}

void ReadSet::LoadAligments() {
  printf("loading aligments from %s\n", name_.c_str());
  ifstream ifs(name_);
  if (ifs.is_open()) {
    boost::archive::binary_iarchive ia(ifs);
    ia >> aligment_cache_;
    ia >> read_lens_;
    ia >> reads_num_;
    ia >> read_map_;
    ia >> read_map_inv_;
    printf("loaded %d cached members\n", (int)aligment_cache_.size());
    CalcMaxReadLen();
    load_success_ = true;
  }
}

void PacbioReadSet::SaveAligments() {
  save_changes_++;
  printf("save ch %d\n", save_changes_);
  if (save_changes_ == 50) {
    ofstream ofs(name_);
    boost::archive::binary_oarchive oa(ofs);
    oa << aligment_cache_;
    oa << read_lens_;
    oa << read_seq_;
    oa << reads_num_;
    oa << read_map_;
    oa << read_map_inv_;
    printf("saved %d cached members\n", (int)aligment_cache_.size());
    save_changes_ = 0;
  }
}

void PacbioReadSet::LoadAligments() {
  printf("loading aligments from %s\n", name_.c_str());
  ifstream ifs(name_);
  if (ifs.is_open()) {
    boost::archive::binary_iarchive ia(ifs);
    ia >> aligment_cache_;
    ia >> read_lens_;
    ia >> read_seq_;
    ia >> reads_num_;
    ia >> read_map_;
    ia >> read_map_inv_;
    printf("loaded %d cached members\n", (int)aligment_cache_.size());
    CalcMaxReadLen();
    load_success_ = true;
  }
}

void PacbioReadSet::NormalizeCache(const Graph& gr) {
  printf("normalize start\n");
  unordered_set<vector<int> > keys;
  for (auto &e: aligment_cache_)
    keys.insert(e.first);
  for (auto e: keys) {
    vector<int> path = e;
    gr.NormalizePath(path);
    aligment_cache_[path] = aligment_cache_[e];
  }
  printf("normalize done\n");
}

void ReadIndexTrivial::AddRead(const string& seq, int read_id) {
  unsigned long long curhash = 0; 
  for (int i = 0; i < kIndexKmer; i++) {
    curhash <<= 2;
    curhash += trans[seq[i]];
  }
  read_index_[curhash].push_back(read_id);
  for (int i = kIndexKmer; i < seq.length(); i++) {
    curhash <<= 2;
    curhash &= (1ll << (2*kIndexKmer)) - 1;
    curhash += trans[seq[i]];
    read_index_[curhash].push_back(read_id);
  }    
}

void ReadIndexTrivial::GetReadCands(const string& seq, unordered_set<int>& read_cands) {
  unsigned long long curhash = 0; 
  for (int i = 0; i < kIndexKmer; i++) {
    curhash <<= 2;
    curhash += trans[seq[i]];
  }
  if (read_index_.count(curhash)) {
    for (int j = 0; j < read_index_[curhash].size(); j++) {
      read_cands.insert(read_index_[curhash][j]);  
    }
  }
  for (int i = kIndexKmer; i < seq.size(); i++) {
    curhash <<= 2;
    curhash &= (1ll << (2*kIndexKmer)) - 1;
    curhash += trans[seq[i]];
    if (read_index_.count(curhash)) {
      for (int j = 0; j < read_index_[curhash].size(); j++) {
        read_cands.insert(read_index_[curhash][j]);  
      }
    }
  }
  string seqr = ReverseSeq(seq);
  curhash = 0; 
  for (int i = 0; i < kIndexKmer; i++) {
    curhash <<= 2;
    curhash += trans[seqr[i]];
  }
  if (read_index_.count(curhash)) {
    for (int j = 0; j < read_index_[curhash].size(); j++) {
      read_cands.insert(read_index_[curhash][j]);  
    }
  }
  for (int i = kIndexKmer; i < seq.size(); i++) {
    curhash <<= 2;
    curhash &= (1ll << (2*kIndexKmer)) - 1;
    curhash += trans[seqr[i]];
    if (read_index_.count(curhash)) {
      for (int j = 0; j < read_index_[curhash].size(); j++) {
        read_cands.insert(read_index_[curhash][j]);  
      }
    }
  }
}

void ReadIndexTrivial::PrintSizeInfo() {
  long long ss = 0;
  for (auto &e: read_index_) {
    ss += 1 + e.second.size();
  }
  printf("read index done, size %lld, %lld\n", read_index_.size(), ss);
}

void ReadIndexMinHash::PrintSizeInfo() {
  long long ss = 0;
  for (auto &e: read_index_) {
    ss += 1 + e.second.size();
  }
  printf("read index done, size %lld, %lld\n", read_index_.size(), ss);
}

unsigned long long ReadIndexMinHash::Hash(unsigned long long x) {
//  return (x << 7) ^ (x >> 5) ^ (x & 0xffaaffaaffaaffaaULL) ^ (x >> 23) ^ (x << 23);
  x = x ^ 0x2204abcd; 
/*  x ^= x >> 16;
  x *= 0x85ebca6b;
  x ^= x >> 13;
  x *= 0xc2b2ae35;
  x ^= x >> 16;*/
  return x;
}

unsigned long long ReadIndexMinHash::GetMinHashForSeq(const string& seq) {
  unsigned long long curhash = 0;
  unsigned long long minhash = 0;
  for (int i = 0; i < kIndexKmer; i++) {
    curhash <<= 2;
    curhash += trans[seq[i]];
  }
  minhash = max(minhash, Hash(curhash));
  for (int i = kIndexKmer; i < seq.length(); i++) {
    curhash <<= 2;
    curhash &= (1ll << (2*kIndexKmer)) - 1;
    curhash += trans[seq[i]];
    minhash = max(minhash, Hash(curhash));
  }  
  return minhash;
}

void ReadIndexMinHash::AddRead(const string& seq, int read_id) {
  unsigned long long minhash = GetMinHashForSeq(seq);
  read_index_[minhash].push_back(read_id);
  read_len = seq.length();
}

void ReadIndexMinHash::GetMinHashWithPoses(
    const string& seq, vector<pair<unsigned long long, int>>& mhs) {
  deque<pair<unsigned long long, int>> d;
  unsigned long long curhash = 0;
  for (int i = 0; i < kIndexKmer; i++) {
    curhash <<= 2;
    curhash += trans[seq[i]];
  }
  unsigned long long mh = Hash(curhash);
  d.push_back(make_pair(mh, kIndexKmer-1));
  unsigned long long last_mh = 0;
  for (int i = kIndexKmer; i < seq.length(); i++) {
    while (!d.empty() && d.front().second < i - read_len + kIndexKmer) {
      d.pop_front();
    }
    curhash <<= 2;
    curhash &= (1ll << (2*kIndexKmer)) - 1;
    curhash += trans[seq[i]];
    unsigned long long mh = Hash(curhash);
    while (!d.empty() && d.back().first < mh) {
      d.pop_back();
    }
    d.push_back(make_pair(mh,i));
    if (i >= read_len - 1) {
      unsigned long long mhx = d.front().first;
      if (i == read_len - 1 || mhx != last_mh) {
        mhs.push_back(make_pair(mhx, d.front().second));
        last_mh = mhx;
      }
    }
  }
}

void ReadIndexMinHash::GetReadCandsWithPoses(
    const string& seq, unordered_map<int, vector<int>>& read_cands) {
  vector<pair<unsigned long long, int>> mhsf;
  GetMinHashWithPoses(seq, mhsf);
  for (auto &e: mhsf) {
    if (read_index_.count(e.first)) {
      for (auto &e2: read_index_[e.first]) {
        assert(GetMinHashForSeq(seq.substr(e.second - kIndexKmer + 1, kIndexKmer))
               == e.first);
        read_cands[e2].push_back(e.second);
      }
    }
  }
  string seqr = ReverseSeq(seq);
  vector<pair<unsigned long long, int>> mhsr;
  GetMinHashWithPoses(seqr, mhsr);
  for (auto &e: mhsr) {
    if (read_index_.count(e.first)) {
      for (auto &e2: read_index_[e.first]) {
        read_cands[e2].push_back(-e.second);
      }
    }
  }
}

void ReadIndexMinHash::GetReadCands(const string& seq, unordered_set<int>& read_cands) {
  // Very stupid version
  for (int i = 0; i + read_len - 1 < seq.length(); i++) {
    string s = seq.substr(i, read_len);
    string sr = ReverseSeq(s);
    unsigned long long mh1 = GetMinHashForSeq(s);
    unsigned long long mh2 = GetMinHashForSeq(sr);
    if (read_index_.count(mh1)) {
      read_cands.insert(read_index_[mh1].begin(), read_index_[mh1].end());
    }
    if (read_index_.count(mh2)) {
      read_cands.insert(read_index_[mh2].begin(), read_index_[mh2].end());
    }
  }
}

void ReadSet::PrepareReadIndex() {
  printf("read index prepare\n");
  ifstream ifs(filename_);
  string l;
  while (getline(ifs, l)) {
    string nameline = l.substr(1);
    vector<string> nameparts;
    split(nameparts, nameline, is_any_of(" \t"));
    string name = nameparts[0];
    int read_id = GetReadId(name);
    string seq;
    getline(ifs, seq);
    read_seqs_[read_id] = seq;
    read_index_.AddRead(seq, read_id);
    getline(ifs, l);
    getline(ifs, l);
  }
  read_index_.PrintSizeInfo();
}

void ReadSet::PreprocessReads() {
  if (load_success_)
    return;

  printf("preprocessing reads\n");
  ifstream ifs(filename_);
  string l;
  unordered_set<int> read_lens;
  while (getline(ifs, l)) {
    string nameline = l.substr(1);
    vector<string> nameparts;
    split(nameparts, nameline, is_any_of(" \t"));
    string name = nameparts[0];
    int read_id = GetReadId(name);
    string seq;
    getline(ifs, seq);
    int read_len = seq.length();
    read_lens.insert(read_len);
    read_lens_[read_id] = read_len;
    getline(ifs, l);
    getline(ifs, l);
  }
  CalcMaxReadLen();
  load_success_ = true;
  printf("read lens: ");
  for (auto &e: read_lens) {
    printf("%d ", e);
  }
  printf("\n");
}

void PacbioReadSet::PreprocessReads() {
  if (load_success_)
      return;
  printf("preprocessing reads pacbio %s\n", filename_.c_str());
  ifstream ifs(filename_);
  assert(ifs);
  string l;
  while (getline(ifs, l)) {
    string nameline = l.substr(1);
    vector<string> nameparts;
    split(nameparts, nameline, is_any_of(" \t"));
    string name = nameparts[0];
    int read_id = GetReadId(name);
    string seq;
    getline(ifs, seq);
    int read_len = seq.length();
    read_seq_[read_id] = seq;
    read_lens_[read_id] = read_len;
    getline(ifs, l);
    getline(ifs, l);
  }
  CalcMaxReadLen();
  printf("preprocess done %d %d\n", (int)read_lens_.size(), max_read_len_);
  load_success_ = true;
}

void ReadSet::CalcMaxReadLen() {
  max_read_len_ = 0;
  for (int i = 0; i < read_lens_.size(); i++) {
    max_read_len_ = max(max_read_len_, read_lens_[i]);
  }
  match_probs_.resize(max_read_len_+7);
  mismatch_probs_.resize(max_read_len_+7);
  for (int i = 0; i < match_probs_.size(); i++) {
    match_probs_[i] = pow(match_prob_, i); 
    mismatch_probs_[i] = pow(mismatch_prob_, i); 
  }
}

void PacbioReadSet::CalcMaxReadLen() {
  max_read_len_ = 0;
  for (int i = 0; i < read_lens_.size(); i++) {
    max_read_len_ = max(max_read_len_, read_lens_[i]);
  }
}

const vector<Aligment>& ReadSet::GetAligmentForSubpath(
    const Graph& gr, const vector<int>& subpath) {
  if (aligment_cache_.count(subpath)) {
    return aligment_cache_[subpath];
  }
  assert(false);
}

void PositionsToReadProbs(
    int num_reads, const vector<vector<pair<int, pair<int, int > > > >& positions,
    const ReadSet& read_set, vector<double>& read_probs) {
  read_probs.clear();
  read_probs.resize(num_reads);
  for(int i = 0; i < positions.size(); i++) {
    for (auto& y: positions[i]) {
      read_probs[i] += read_set.mismatch_probs_[y.second.first] *
                       read_set.match_probs_[read_set.GetReadLen(i) - y.second.first];
    }
  }
}

double GetTotalProb(const vector<double>& read_probs, int total_len, int& zero_reads,
                    double min_prob_per_base, double min_prob_start, ReadSet& rs1,
                    ReadSet& rs2) {
  int total_c = 0;
  double total_prob = 0;
  if (total_len == 0) {
    total_len = 1;
  }
  zero_reads = 0;
  for (int i = 0; i < read_probs.size(); i++) {
    double prob = read_probs[i] / (2*total_len);
    double threshold = exp(min_prob_start + 
                           min_prob_per_base*(rs1.GetReadLen(i)+rs2.GetReadLen(i)));
    if (prob < threshold) {
      zero_reads++;
      prob = threshold;
    }
    total_prob += log(prob);
    total_c++;
  }
  return total_prob / total_c;
}

double GetTotalProb(const vector<double>& read_probs, int total_len, int& zero_reads,
                    double min_prob_per_base, double min_prob_start, ReadSet& rs) {
  int total_c = 0;
  double total_prob = 0;
  if (total_len == 0) {
    total_len = 1;
  }
  zero_reads = 0;
  for (int i = 0; i < read_probs.size(); i++) {
    double prob = read_probs[i] / (2*total_len);
    double threshold = exp(min_prob_start + min_prob_per_base*(rs.GetReadLen(i)));
    if (prob < threshold) {
      zero_reads++;
      prob = threshold;
    }
    total_prob += log(prob);
    total_c++;
  }
  return total_prob / total_c;
}
double GetTotalProb(const vector<double>& read_probs, int total_len, int& zero_reads,
                    double threshold) {
  int total_c = 0;
  double total_prob = 0;
  if (total_len == 0) {
    total_len = 1;
  }
  zero_reads = 0;
  for (auto &x: read_probs) {
    double prob = x / (2*total_len);
    if (prob < threshold) {
      zero_reads++;
      prob = threshold;
    }
    total_prob += log(prob);
    total_c++;
  }
  return total_prob / total_c;
}

double GetTotalProb(const vector<double>& read_probs, int total_len, int& zero_reads) {
  int total_c = 0;
  double total_prob = 0;
  if (total_len == 0) {
    total_len = 1;
  }
  zero_reads = 0;
  for (auto &x: read_probs) {
    double prob = x / (2*total_len);
    if (prob < kThresholdProb) {
      zero_reads++;
      prob = kThresholdProb;
    }
    total_prob += log10(prob);
    total_c++;
  }
  return total_prob / total_c;
}

double CalcScoreForPath(const Graph& gr, const vector<int>& path, int kmer,
                        ReadSet& read_set, bool use_caching) {
  int total_len;
  vector<vector<pair<int, pair<int, int> > > >& positions = 
      use_caching ? read_set.GetPositions(gr, path, total_len) :
      read_set.GetPositionsSlow(gr, path, total_len);
  vector<double> read_probs;
  PositionsToReadProbs(read_set.GetNumberOfReads(), positions, read_set, read_probs);

  int zero_reads;
  double total_prob = GetTotalProb(read_probs, total_len, zero_reads);
  printf("zero reads %d\n", zero_reads);
  return total_prob;
}

double GetInsertProbability(double insert_len, double insert_mean, double insert_std) {
  double z = (insert_len - insert_mean) / insert_std;
  double e = exp(-z*z/2.0);
  double C = sqrt(2*M_PI)*insert_std;
  return e/C;
}

double CalcScoreForPath(const Graph& gr, const vector<int>& path, int kmer, 
                        ReadSet& read_set1, ReadSet& read_set2, 
                        double insert_mean, double insert_std,
                        bool use_caching) {
  assert(read_set1.GetNumberOfReads() == read_set2.GetNumberOfReads());
  int total_len1, total_len2;
  vector<vector<pair<int, pair<int, int> > > >& positions1 = 
      use_caching ? read_set1.GetPositions(gr, path, total_len1) :
      read_set1.GetPositionsSlow(gr, path, total_len1);
  vector<vector<pair<int, pair<int, int> > > >& positions2 = 
      use_caching ? read_set2.GetPositions(gr, path, total_len2) :
      read_set2.GetPositionsSlow(gr, path, total_len2);

  assert(total_len1 == total_len2);

  vector<double> insert_probs((int)(insert_mean + 5*insert_std));
  for (int i = 0; i < insert_probs.size(); i++) {
    insert_probs[i] = GetInsertProbability(i, insert_mean, insert_std);
  }
  vector<double> read_probs(read_set1.GetNumberOfReads());
  for (int i = 0; i < read_set1.GetNumberOfReads(); i++) {
    for (auto &x: positions1[i]) {
      double p1 = read_set1.mismatch_probs_[x.second.first] *
                  read_set1.match_probs_[read_set1.GetReadLen(i) - x.second.first];
      for (auto &y: positions2[i]) {
        double p2 = read_set2.mismatch_probs_[y.second.first] *
                    read_set2.match_probs_[read_set2.GetReadLen(i) - y.second.first];
        if (x.second.second == y.second.second) continue;
        int dist;
        if (x.first < y.first) {
          dist = y.first - x.first - read_set1.GetReadLen(i);
        } else {
          dist = x.first - y.first - read_set2.GetReadLen(i);
        }
        double insprob;
        if (dist < insert_probs.size()) {
          insprob = insert_probs[dist];
        } else {
          insprob = GetInsertProbability(dist, insert_mean, insert_std);
        }
        read_probs[i] += p1*p2*insprob;
      }
    }
  }
  int zero_reads;
  double total_prob = GetTotalProb(read_probs, total_len1, zero_reads);
  printf("zero reads %d/%d\n", zero_reads, read_set1.GetNumberOfReads());
  return total_prob;
}

double CalcScoreForPaths(const Graph& gr, const vector<vector<int>>& paths, 
                         ReadSet& read_set1,
                         int &zero_reads, int &total_len,
                         bool use_caching, double no_cov_penalty,
                         double exp_cov_move,
                         double min_prob_per_base, double min_prob_start) {
//  printf("calc score\n");
  int total_len1 = 0;
  vector<double> read_probs(read_set1.GetNumberOfReads());
  read_set1.ClearPositions();
  int st = 0;
  // (position, type)
  // type: 3 - begin, 4 - end, 1 - start path, 2 - end path
  vector<pair<int, int> > events;
  
  for (auto &path: paths) {
    vector<vector<int>> ctgs;
    vector<int> gaps;
    int last = 0;
    for (int i = 0; i < path.size(); i++) {
      if (path[i] < 0) {
        gaps.push_back(-path[i]);
        ctgs.push_back(vector<int>(path.begin()+last, path.begin()+i));
        last = i+1;
      }
    }
    ctgs.push_back(vector<int>(path.begin()+last, path.end()));
    events.push_back(make_pair(st + total_len1, 1));
    for (int i = 0; i < ctgs.size(); i++) {
      if (i > 0) {
        total_len1 += gaps[i-1];
        events.push_back(make_pair(st + total_len1, 1));
      }
      read_set1.AddPositions(gr, ctgs[i], total_len1, st + total_len1);
    }
    st += 1000000;
  }
  vector<vector<pair<int, pair<int, int> > > >& positions1 = 
      read_set1.GetPositions();

  for (int i = 0; i < read_set1.GetNumberOfReads(); i++) {
    for (auto &x: positions1[i]) {
      double p1 = read_set1.mismatch_probs_[x.second.first] *
                  read_set1.match_probs_[read_set1.GetReadLen(i) - x.second.first];
      if (p1 > kThresholdProb2) {
        events.push_back(make_pair(x.first, read_set1.GetReadLen(i)));
      }
      read_probs[i] += p1;
    }
  }
  sort(events.begin(), events.end());
  int last_event_pos = 0;
  int last_fin = -1;
  int last_event_type = -1;
  int last_begin = 0;
  int bad_bases = 0;
  int bad_gaps = 0;
  int last_gap = 0;
  set<int> bc;
  int bn = 0;
  for (int i = 0; i < events.size(); i++) {
    if (events[i].second >= 3) {
      if (events[i].first > last_fin && (last_event_type >= 3)) {
//        printf("bad gap %d %d %d %d %d %d\n", events[i].first - last_begin, last_event_pos - last_begin, events[i].first -
//               last_event_pos, last_gap - last_begin, last_event_type, events[i].first - last_gap);
        bad_bases += events[i].first - last_fin;
//        printf("gap %d %d\n", events[i].first - last_event_pos, events[i].first - last_begin);
        bad_gaps++;
        bc.insert(bn);
      }
      last_fin = max(last_fin, (int)(events[i].first + events[i].second*exp_cov_move));
    }
    if (events[i].second == 1) {
      last_begin = events[i].first;
      last_event_pos = events[i].first;
      last_event_type = events[i].second;
      bn++;
    }
    if (events[i].second < -1) {
      last_gap = events[i].first;
      last_event_pos = events[i].first;
      last_event_type = events[i].second;
    }
  }
  if (no_cov_penalty > 0) {
    printf("bad (single) %d %d %d\n", bad_bases, bad_gaps, (int)bc.size());
  }
//  int zero_reads;
  double total_prob = GetTotalProb(read_probs, total_len1, zero_reads, 
                                   min_prob_per_base, min_prob_start, read_set1);
  total_len = total_len1;
//  printf("zero reads %d/%d %d\n", zero_reads, read_set1.GetNumberOfReads());
  return total_prob - bad_bases*no_cov_penalty;
}

double CalcScoreForPaths(const Graph& gr, const vector<vector<int>>& paths, 
                         ReadSet& read_set1, ReadSet& read_set2, 
                         double insert_mean, double insert_std,
                         int &zero_reads, int &total_len,
                         bool use_caching, double no_cov_penalty,
                         double exp_cov_move, bool use_all_to_cov,
                         double min_prob_per_base, double min_prob_start) {
//  printf("calc score\n");  
  printf("ins m %lf %lf %lf\n", insert_mean, no_cov_penalty, exp_cov_move);
  assert(read_set1.GetNumberOfReads() == read_set2.GetNumberOfReads());
  int total_len1 = 0, total_len2 = 0;
  vector<double> read_probs(read_set1.GetNumberOfReads());
  read_set1.ClearPositions();
  read_set2.ClearPositions();
  read_set1.PrecomputeAlignmentForPaths(paths, gr);
  read_set2.PrecomputeAlignmentForPaths(paths, gr);
  int st = 0;
  // (position, type)
  // type: 3 - begin, 4 - end, 1 - start path, 2 - end path
  vector<pair<int, int> > events;
  
  int overins = 0;
  for (auto &path: paths) {
/*    int tl1, tl2;
    vector<vector<pair<int, pair<int, int> > > >& positions1 = 
        use_caching ? read_set1.GetPositions(gr, path, kmer, tl1) :
        read_set1.GetPositionsSlow(gr, path, kmer, tl1);
    vector<vector<pair<int, pair<int, int> > > >& positions2 = 
        use_caching ? read_set2.GetPositions(gr, path, kmer, tl2) :
        read_set2.GetPositionsSlow(gr, path, kmer, tl2);
    total_len1 += tl1;
    total_len2 += tl2;*/
    vector<vector<int>> ctgs;
    vector<int> gaps;
    int last = 0;
    int scfl = 0;
    for (int i = 0; i < path.size(); i++) {
      if (path[i] < 0) {
        gaps.push_back(-path[i]);
        scfl += -path[i];
        ctgs.push_back(vector<int>(path.begin()+last, path.begin()+i));
        last = i+1;
      } else {
        scfl += gr.nodes[path[i]]->s.length();
      }
    }
    if (scfl > insert_mean) overins++;
    ctgs.push_back(vector<int>(path.begin()+last, path.end()));
    events.push_back(make_pair(st + total_len1, 1));
    for (int i = 0; i < ctgs.size(); i++) {
      if (i > 0) {
        int b = st + total_len1 + insert_mean - insert_std;
        total_len1 += gaps[i-1];
        total_len2 += gaps[i-1];
        int e = st + total_len1 + insert_mean + insert_std;
        events.push_back(make_pair(st + total_len1, 1));
      }
      read_set1.AddPositions(gr, ctgs[i], total_len1, st + total_len1);
      read_set2.AddPositions(gr, ctgs[i], total_len2, st + total_len2);
    }
    assert(total_len1 == total_len2);
    st += 1000000;
  }
  printf("overins %lf %d\n", insert_mean, overins);
  vector<vector<pair<int, pair<int, int> > > >& positions1 = 
      read_set1.GetPositions();
  vector<vector<pair<int, pair<int, int> > > >& positions2 = 
      read_set2.GetPositions();

  vector<double> insert_probs((int)(insert_mean + 5*insert_std));
  for (int i = 0; i < insert_probs.size(); i++) {
    insert_probs[i] = GetInsertProbability(i, insert_mean, insert_std);
  }
//  vector<int> dists;
  for (int i = 0; i < read_set1.GetNumberOfReads(); i++) {
    double threshold = exp(min_prob_start + 
                           min_prob_per_base*(read_set1.GetReadLen(i)+read_set2.GetReadLen(i)));
    for (auto &x: positions1[i]) {
      double p1 = read_set1.mismatch_probs_[x.second.first] *
                  read_set1.match_probs_[read_set1.GetReadLen(i) - x.second.first];
      for (auto &y: positions2[i]) {
        double p2 = read_set2.mismatch_probs_[y.second.first] *
                    read_set2.match_probs_[read_set2.GetReadLen(i) - y.second.first];
        if (x.second.second == y.second.second) continue;
        int dist;
        if (x.first < y.first) {
          if (x.second.second != 0 || y.second.second != 1) {
            continue;
          }
          dist = y.first - x.first + read_set2.GetReadLen(i);
        } else {
          if (x.second.second != 1 || y.second.second != 0) {
            continue;
          }
          dist = x.first - y.first + read_set1.GetReadLen(i);
        }
        double insprob;
        if (dist < insert_probs.size()) {
          insprob = insert_probs[dist];
        } else {
          insprob = GetInsertProbability(dist, insert_mean, insert_std);
        }
//        dists.push_back(dist);
        if (p1*p2*insprob > threshold) {
          events.push_back(make_pair(max(x.first, y.first), 3));
          if (use_all_to_cov) {
            events.push_back(make_pair(min(x.first, y.first), 3));
          }
        }
        read_probs[i] += p1*p2*insprob;
      }
    }
  }
/*  sort(dists.begin(), dists.end());
  printf("%d: %d %d %d %d %d\n", (int)dists.size(),
         dists[dists.size()/10],
         dists[dists.size()/4],
         dists[dists.size()/2],
         dists[3*dists.size()/4],
         dists[9*dists.size()/10]);*/
  sort(events.begin(), events.end());
  int last_event_pos = 0;
  int last_event_type = -1;
  int last_begin = 0;
  int bad_bases = 0;
  int bad_gaps = 0;
  int last_gap = 0;
  set<int> bc;
  int bn = 0;
  for (int i = 0; i < events.size(); i++) {
    if (events[i].second == 3) {
      if (events[i].first - last_event_pos > exp_cov_move && 
          (last_event_type == 3 || last_event_type < 0) && events[i].first - last_begin > insert_mean + 5*insert_std) {
//        printf("bad gap %d %d %d %d %d %d\n", events[i].first - last_begin, last_event_pos - last_begin, events[i].first -
//               last_event_pos, last_gap - last_begin, last_event_type, events[i].first - last_gap);
        bad_bases += events[i].first - last_event_pos;
//        printf("gap %d %d\n", events[i].first - last_event_pos, events[i].first - last_begin);
        bad_gaps++;
        bc.insert(bn);
      }
    }
    if (events[i].second == 1) {
      last_begin = events[i].first;
      bn++;
    }
    if (events[i].second < -1) {
      last_gap = events[i].first;
    }
    last_event_pos = events[i].first;
    last_event_type = events[i].second;
  }
  if (no_cov_penalty > 0) {
    printf("bad (%lf) %d %d %d\n", insert_mean, bad_bases, bad_gaps, (int)bc.size());
  }
//  int zero_reads;
  double total_prob = GetTotalProb(read_probs, total_len1, zero_reads, min_prob_per_base,
                                   min_prob_start, read_set1, read_set2);
  total_len = total_len1;
//  printf("zero reads %d/%d %d\n", zero_reads, read_set1.GetNumberOfReads());
  return total_prob - bad_bases*no_cov_penalty;
}

string ExpandCigar(const vector<pair<int, char>>& cigar) {
  string ret;
  for (int i = 0; i < cigar.size(); i++) {
    for (int j = 0; j < cigar[i].first; j++)
      ret += cigar[i].second;
  }
  return ret;
}

void GetCigarEnds(const string& cigar, int& bl, int& el) {
  for (int i = 0; i < cigar.length(); i++) {
    if (cigar[i] != 'I') {
      bl = i;
      break;
    }
  }
  for (int i = cigar.length()-1; i >= 0; i--) {
    if (cigar[i] != 'I') {
      el = cigar.length() - i;
      break;
    }
  }
}

void Uniquify(vector<pair<int, int> > &x) {
  if (x.empty()) return;
  int mi = x[0].first;
  int ma = x[0].first;
  for (auto &e: x) {
    mi = min(e.first, mi);
    ma = max(e.first, ma);
  }

  vector<pair<int, int> > ss(ma-mi+1, make_pair(1000000, -1000000));
  for (auto &e: x) {
    ss[e.first - mi].first = min(e.second, ss[e.first - mi].first);
    ss[e.first - mi].second = max(e.second, ss[e.first - mi].second);
  }
  x.clear();
  for (int i = mi; i <= ma; i++) {
    for (int j = ss[i-mi].first; j <= ss[i-mi].second; j++) {
      x.push_back(make_pair(i, j));
    }
  }
}

logdouble PacbioReadSet::AligmentProbability(
    const std::string &s1, const std::string &s2,
    const PacbioAligmentData& align_data, int band) const {
  string cigar = ExpandCigar(align_data.cigar);
  int bl, el;
  GetCigarEnds(cigar, bl, el);
  bl = min(bl, 200);
  el = min(el, 200);
  vector<pair<int, int> > positions;
  int currow = 0;
  int curcol = 0;
  positions.push_back(make_pair(0, 0));
  for (int i = -bl; i < 3; i++) {
    for (int j = 0; j < bl; j++) {
      positions.push_back(make_pair(i, j));
    }
  }
  for (int i = 0; i < cigar.length(); i++) {
    if (cigar[i] == 'M') {
      currow++;
      curcol++;
    } else if (cigar[i] == 'I') {
      curcol++;
    } else if (cigar[i] == 'D') {
      currow++;
    }
    positions.push_back(make_pair(currow, curcol));
  }
  for (int i = currow; i < currow + el; i++) {
    for (int j = curcol - el; j <= curcol; j++) {
      positions.push_back(make_pair(i, j));
    }
  }
  Uniquify(positions);
  vector<pair<int, int> > add_positions;
  for (auto &e: positions) {
    int bband = band;
    for (int i = -bband; i <= bband; i++) {
      for (int j = -bband; j <= bband; j++) {
        add_positions.push_back(make_pair(e.first+i, e.second+j));
      }
    }
  }
  for (auto &e: add_positions) {
    positions.push_back(e);
  }
  Uniquify(positions);

  int offset = positions[0].first;
  int sss = positions.back().first - offset + 1;
  vector<int> row_offsets(sss, positions.back().second + 1000000);
  vector<vector<logdouble> > results(sss);

  for (auto &e: positions) {
    row_offsets[e.first - offset] = min(row_offsets[e.first - offset], e.second);
  }
  for (auto &e: positions) {
    if (results[e.first - offset].size() <= e.second - row_offsets[e.first - offset]) {
      results[e.first - offset].resize(e.second - row_offsets[e.first - offset] + 1);
    }
  }

  logdouble ret = 0;
  char gap = '-';
  for (auto &e: positions) {
    if (e.second == 0) {
      results[e.first - offset][e.second - row_offsets[e.first - offset]] = 1;
    }
  }

  for (auto &e: positions) {
    if (e.second == 0) continue;
    if (e.second - 1 < 0 || e.second - 1 >= s2.length()) continue;
    if (e.first + align_data.posstart - 1 < 0 || e.first + align_data.posstart - 1 >= s1.length())
      continue;
    pair<int, int> e2 = make_pair(e.first-1, e.second-1);
    if (e2.first - offset >= 0 &&
        e2.second - row_offsets[e2.first - offset] >= 0 &&
        e2.second - row_offsets[e2.first - offset] < results[e2.first - offset].size()) {
      logdouble p = MatchProbability(s1[e.first + align_data.posstart-1], s2[e.second-1]);
      results[e.first - offset][e.second - row_offsets[e.first - offset]] +=
          results[e2.first - offset][e2.second - row_offsets[e2.first - offset]] * p;
    }
    e2 = make_pair(e.first-1, e.second);
    if (e2.first - offset >= 0 &&
        e2.second - row_offsets[e2.first - offset] >= 0 &&
        e2.second - row_offsets[e2.first - offset] < results[e2.first - offset].size()) {
      results[e.first - offset][e.second - row_offsets[e.first - offset]] +=
          results[e2.first - offset][e2.second - row_offsets[e2.first - offset]] *
          MatchProbability(s1[e.first + align_data.posstart-1], gap); 
    }
    e2 = make_pair(e.first, e.second-1);
    if (e2.first - offset >= 0 &&
        e2.second - row_offsets[e2.first - offset] >= 0 &&
        e2.second - row_offsets[e2.first - offset] < results[e2.first - offset].size()) {
      results[e.first - offset][e.second - row_offsets[e.first - offset]] +=
          results[e2.first - offset][e2.second - row_offsets[e2.first - offset]] *
          MatchProbability(gap, s2[e.second-1]);  
    }
    if (results[e.first - offset][e.second - row_offsets[e.first - offset]].logval != 
        results[e.first - offset][e.second - row_offsets[e.first - offset]].logval) {
      printf("wtf %d %d\n", e.first, e.second);
      assert(false);
    }
    if (e.second == s2.length()) {
      ret += results[e.first - offset][e.second - row_offsets[e.first - offset]];
    }
  }

/*  FILE *f = fopen("ap.dat", "a");
  int ll = align_data.slen;
  int bad = align_data.edit_dist + (align_data.slen - (align_data.send - align_data.sstart));
  int good = ll - bad;
  double rp = -0.13397892 * good -2.871402 * bad; 
  double rp2 = match_prob_.logval * good + mismatch_prob_.logval * bad; 
  fprintf(f, "%d %d %d %lf %lf %lf\n", align_data.edit_dist,
          align_data.edit_dist + (align_data.slen - (align_data.send - align_data.sstart)),
          align_data.slen, rp, rp2,
          ret.logval);
  fclose(f);*/

  return ret;
}

vector<vector<pair<int, logdouble> > >& PacbioReadSet::GetExactReadProbabilities(
    const Graph& gr, const vector<int>& path, int ps, int& total_len,
    int& total_len2) {
  string seq = gr.nodes[path[0]]->s;
  vector<int> pathnodesposes, pathnodesposesb;
  pathnodesposes.push_back(seq.length());
  pathnodesposesb.push_back(0);
  total_len2 = seq.length();
  int back_length = 0;
  for (int i = 1; i < path.size(); i++) {
    pathnodesposesb.push_back(seq.length());
    if (i < ps) {
      total_len2 += gr.nodes[path[i]]->s.length();
    } else {
      back_length += gr.nodes[path[i]]->s.length();
    }
    seq += gr.nodes[path[i]]->s;
    pathnodesposes.push_back(seq.length());
  }
  total_len = seq.length();
  total_len2 += min(max_read_len_ / 3, back_length);
//  printf("total len %d\n", total_len);

  vector<vector<int> > subpaths;
  vector<vector<int> > subpaths_inds;
//  bool missing = false;
  vector<pair<int, int> > missing;
  for (int i = 0; i < path.size(); i++) {
    vector<int> subpath, subpathind;
    for (int j = i; j < path.size(); j++) {
      subpath.push_back(path[j]);
      subpathind.push_back(j);
      int subpath_length = pathnodesposes[j] - pathnodesposesb[i]; 
      int first_length = pathnodesposes[i] - pathnodesposesb[i];
      if (aligment_cache_.count(subpath) == 0) {
        missing.push_back(make_pair(subpathind[0], subpathind.back()));
      }
      subpaths.push_back(subpath);
      subpaths_inds.push_back(subpathind);
      if (subpath_length - first_length > max_read_len_) {
        break;
      }
    }
  }
  if (!missing.empty()) {
    int lastmissend = -47;
    int lastmissbegin = -47; 
    sort(missing.begin(), missing.end());
    for (int i = 0; i < missing.size(); i++) {
      if (missing[i].first > lastmissend) {
        if (lastmissend != -47) {
          int tl;
          GetReadProbabilitiesSlow(
              gr, vector<int>(path.begin()+lastmissbegin, path.begin()+lastmissend+1),
              tl);
        }
        lastmissbegin = missing[i].first;
        lastmissend = missing[i].second;
      }
      lastmissend = max(lastmissend, missing[i].second);
    }
    if (lastmissend != -47) {
      int tl;
      GetReadProbabilitiesSlow(
          gr, vector<int>(path.begin()+lastmissbegin, path.begin()+lastmissend+1),
          tl);
    }
  }

//  printf("subpaths size %d\n", subpaths.size());
  int aa = 0;
  set<int> rr;
  for (int i = 0; i < positions_.size(); i++) {
    positions_[i].clear();
  }
  positions_.resize(reads_num_);

  subpaths.clear();
  for (int i = 0; i < path.size() && i < ps; i++) {
    vector<int> subpath, subpathind;
    for (int j = i; j < path.size(); j++) {
      subpath.push_back(path[j]);
      subpathind.push_back(j);
      int subpath_length = pathnodesposes[j] - pathnodesposesb[i]; 
      int first_length = pathnodesposes[i] - pathnodesposesb[i];
      if (aligment_cache_.count(subpath) == 0) {
        missing.push_back(make_pair(subpathind[0], subpathind.back()));
      }
      subpaths.push_back(subpath);
      subpaths_inds.push_back(subpathind);
      if (subpath_length - first_length > max_read_len_) {
        break;
      }
    }
  }
  for (int i = 0; i < subpaths.size(); i++) {
/*    assert(subpaths[i].size() > 0);
    assert(subpaths[i][0] >= 0);
    assert(subpaths[i][0] < pathnodesposesb.size());*/
//    int pos_begin = pathnodesposesb[subpaths[i][0]];
    assert(aligment_cache_.count(subpaths[i]));
    for (int j = 0; j < aligment_cache_[subpaths[i]].size(); j++) {
      aa++;
      PacbioAligment &al = aligment_cache_[subpaths[i]][j];
      rr.insert(al.read_id);
      positions_[al.read_id].push_back(make_pair(al.position, al.prob));
    }
  }
  return positions_;
}

vector<vector<pair<pair<int, int>, logdouble> > >& PacbioReadSet::GetReadProbabilities(
    const Graph& gr, const vector<int>& path, int& total_len) {
  string seq = gr.nodes[path[0]]->s;
  vector<int> pathnodesposes, pathnodesposesb;
  pathnodesposes.push_back(seq.length());
  pathnodesposesb.push_back(0);
/*  printf("path: ");
  for (int i = 0; i < path.size(); i++) {
    printf("%d ", path[i]);
  }
  printf("\n");*/
  for (int i = 1; i < path.size(); i++) {
    pathnodesposesb.push_back(seq.length());
    if (path[i] < 0) {
      for (int j = 0; j < -path[i]; j++)
        seq += "N";
    } else {
      seq += gr.nodes[path[i]]->s;
    }
    pathnodesposes.push_back(seq.length());
  }
  total_len = seq.length();
//  printf("total len %d\n", total_len);

  vector<vector<int> > subpaths;
  vector<vector<int> > subpaths_inds;
//  bool missing = false;
  vector<pair<int, int> > missing;
  for (int i = 0; i < path.size(); i++) {
    vector<int> subpath, subpathind;
    for (int j = i; j < path.size(); j++) {
      subpath.push_back(path[j]);
      subpathind.push_back(j);
      int subpath_length = pathnodesposes[j] - pathnodesposesb[i]; 
      int first_length = pathnodesposes[i] - pathnodesposesb[i];
      if (aligment_cache_.count(subpath) == 0) {
        missing.push_back(make_pair(subpathind[0], subpathind.back()));
      }
      subpaths.push_back(subpath);
      subpaths_inds.push_back(subpathind);
      if (subpath_length - first_length > max_read_len_) {
        break;
      }
    }
  }
  if (!missing.empty()) {
    int lastmissend = -47;
    int lastmissbegin = -47; 
    sort(missing.begin(), missing.end());
    for (int i = 0; i < missing.size(); i++) {
      if (missing[i].first > lastmissend) {
        if (lastmissend != -47) {
          int tl;
          GetReadProbabilitiesSlow(
              gr, vector<int>(path.begin()+lastmissbegin, path.begin()+lastmissend+1),
              tl);
        }
        lastmissbegin = missing[i].first;
        lastmissend = missing[i].second;
      }
      lastmissend = max(lastmissend, missing[i].second);
    }
    if (lastmissend != -47) {
      int tl;
      GetReadProbabilitiesSlow(
          gr, vector<int>(path.begin()+lastmissbegin, path.begin()+lastmissend+1),
          tl);
    }
  }

//  printf("subpaths size %d\n", subpaths.size());
  int aa = 0;
  set<int> rr;
  for (int i = 0; i < positions2_.size(); i++) {
    positions2_[i].clear();
  }
  positions2_.resize(reads_num_);
  for (int i = 0; i < subpaths.size(); i++) {
/*    assert(subpaths[i].size() > 0);
    assert(subpaths[i][0] >= 0);
    assert(subpaths[i][0] < pathnodesposesb.size());*/
    int pos_begin = pathnodesposesb[subpaths_inds[i][0]];
    assert(aligment_cache_.count(subpaths[i]));
    for (int j = 0; j < aligment_cache_[subpaths[i]].size(); j++) {
      aa++;
      PacbioAligment &al = aligment_cache_[subpaths[i]][j];
      rr.insert(al.read_id);
      positions2_[al.read_id].push_back(make_pair(make_pair(pos_begin + al.position,
                  pos_begin + al.position_end), al.prob));
    }
  }

  return positions2_;
}

void PacbioReadSet::ComputeAnchors(const Graph& gr) {
  string anchorsname = name_ + ".anchors";
  printf("loading anchors from %s\n", anchorsname.c_str());
  ifstream ifs(anchorsname);
  if (ifs.is_open()) {
    boost::archive::binary_iarchive ia(ifs);
    ia >> anchors_cache_;
    ia >> anchors_begin_;
    ia >> anchors_end_;
    printf("loaded %d anchors\n", (int)anchors_cache_.size());
    ifs.close();
  } else {
    char tmpname1[L_tmpnam+6], tmpname2[L_tmpnam];
    tmpnam(tmpname1);
    strcat(tmpname1, ".fas");
    tmpnam(tmpname2);
    printf("anchor files %s %s\n", tmpname1, tmpname2);
    FILE *f = fopen(tmpname1, "w");
    for (int i = 0; i < gr.nodes.size(); i++) {
      if (gr.nodes[i]->s.length() < kMinAnchorLen) continue;

      fprintf(f, ">%d\n", i);
      fprintf(f, "%s\n", gr.nodes[i]->s.c_str());
    }
    fclose(f);
    string cmd = gBlasrPath + "/blasr ";
    cmd += filename_;
    cmd += " ";
    cmd += tmpname1;
    cmd += " -sdpTupleSize 8 -guidedAlignBandSize 100 -nCandidates 50 ";
    cmd += "-minMatch 11 ";
    cmd += kThreadsBlasr;
    cmd += " >";
    cmd += tmpname2;
    system(cmd.c_str());
    ifstream fi(tmpname2);
    string l;
    while (getline(fi, l)) {
      vector<string> parts;
      split(parts, l, is_any_of(" "));
      int node_id = atoi(parts[1].c_str());
      int lastsep = 0;
      for (int i = 0; i < parts[0].length(); i++) {
        if (parts[0][i] == '/') {
          lastsep = i;
        }
      }
      string name = parts[0].substr(0, lastsep);
      int start = atoi(parts[6].c_str());
      int end = atoi(parts[7].c_str());
      anchors_cache_[node_id].insert(GetReadId(name));
      if (start <= 10) {
        anchors_begin_[node_id].insert(GetReadId(name));
      }
      if (end >= gr.nodes[node_id]->s.length() - 10) {
        anchors_end_[node_id].insert(GetReadId(name));
      }
    }

    ofstream ofs(anchorsname);
    boost::archive::binary_oarchive oa(ofs);
    oa << anchors_cache_;
    oa << anchors_begin_;
    oa << anchors_end_;
    printf("saved %d anchors\n", (int)anchors_cache_.size());
  }
  for (auto &e: anchors_begin_) {
    for (auto &x: e.second) {
      anchors_reverse_[x].insert(e.first);
    }
  }
}

int PacbioReadSet::GetGap(const Graph& gr, int first, int second, int read_id) {
  char tmpname1[L_tmpnam+6], tmpname2[L_tmpnam+6], tmpname3[L_tmpnam];
  tmpnam(tmpname1);
  strcat(tmpname1, ".fas");
  tmpnam(tmpname2);
  strcat(tmpname2, ".fq");
  tmpnam(tmpname3);
  printf("gap files %s %s %s\n", tmpname1, tmpname2, tmpname3);
  FILE *fn = fopen(tmpname1, "w");
  fprintf(fn, ">f\n");
  fprintf(fn, "%s\n", gr.nodes[first]->s.c_str());
  fprintf(fn, ">s\n");
  fprintf(fn, "%s\n", gr.nodes[second]->s.c_str());
  fclose(fn);
  unordered_set<int> rs;
  rs.insert(read_id);
  FilterReads(tmpname2, rs);
  string reads_filename = tmpname2;
  string cmd = gBlasrPath + "/blasr ";
  cmd += reads_filename;
  cmd += " ";
  cmd += tmpname1;
  cmd += " -sam -sdpTupleSize 8 -guidedAlignBandSize 100 -nCandidates 50 ";
  cmd += "-minMatch 11 ";
  cmd += kThreadsBlasr;
  cmd += " >";
  cmd += tmpname3;
  system(cmd.c_str());

  ifstream fi(tmpname3);
  PacbioAligmentData first_align, second_align;
  first_align.posend = -2000;
  second_align.posstart = 2000000;
  string l;
  int flen = gr.nodes[first]->s.length();
  int slen = gr.nodes[second]->s.length();
  while (getline(fi, l)) {
    if (l[0] == '@') {
      continue;
    }
    vector<string> parts;
    split(parts, l, is_any_of("\t"));
    int x = 0;
    if (parts[2][0] == 's') x = 1;
    int len = (x == 0 ? flen : slen);
    PacbioAligmentData align = ParseAligment(l, 2*len, false);
    if (x == 0 && align.posend > first_align.posend) 
      first_align = align;
    if (x == 1 && align.posstart < second_align.posstart) 
      second_align = align;
  }
  if (first_align.posend == -2000 || second_align.posstart == 2000000) {
    return -1;
  }
  if ((first_align.flags & 16) != (second_align.flags & 16)) {
    return -2;
  }
  if (second_align.posstart > 10) {
    printf("-3 fail %d %d\n", second_align.posstart, second_align.posend);
    return -3;
  }
  if (first_align.posend < flen - 10) {
    printf("-4 fail %d %d %d\n", first_align.posstart, first_align.posend, flen);
    return -4;
  }
  if (first_align.send > second_align.sstart) {
    return -5;
  }
  return flen - first_align.posend + second_align.posstart + 
         second_align.sstart - first_align.send;
}

vector<vector<pair<int, logdouble> > >& PacbioReadSet::GetReadProbabilitiesSlow(
    const Graph& gr, const vector<int>& path, int& total_len, bool save_to_cache) {
  char tmpname1[L_tmpnam+6], tmpname2[L_tmpnam+6], tmpname3[L_tmpnam];
  tmpnam(tmpname1);
  strcat(tmpname1, ".fas");
  tmpnam(tmpname2);
  strcat(tmpname2, ".fq");
  tmpnam(tmpname3);
  printf("pb slow files %s %s %s %d\n", tmpname1, tmpname2, tmpname3, path.size());
  FILE *f = fopen(tmpname1, "w");
  int id = 0;
  fprintf(f, ">tmp\n");
  string seq;
  if (path[0] >= 0) 
    seq = gr.nodes[path[0]]->s;
  else
    for (int j = 0; j < -path[0]; j++)
      seq += "N";
  vector<int> pathnodesposes, pathnodesposesb;
  pathnodesposes.push_back(seq.length());
  pathnodesposesb.push_back(0);
  for (int i = 1; i < path.size(); i++) {
    printf("path %d %d\n", path[i], gr.nodes.size());
    pathnodesposesb.push_back(seq.length());
    if (path[i] < 0) {
      for (int j = 0; j < -path[i]; j++)
        seq += "N";
    } else {
      seq += gr.nodes[path[i]]->s;
    }
    pathnodesposes.push_back(seq.length());
  }
  total_len = seq.length();
  printf("slow len %d\n", total_len);
  fprintf(f, "%s\n", seq.c_str());
  fclose(f);

  string seqrev = ReverseSeq(seq);
  string seqall = seq + kContigSeparator + seqrev;
  string reads_filename = filename_;

  unordered_set<int> read_filter;
  for (int i = 0; i < path.size(); i++) {
    if (path[i] >= 0) {
      for (auto &e: anchors_cache_[path[i]]) {
        read_filter.insert(e);
      }
    }
  }
  printf("read filter %d/%d\n", (int)read_filter.size(), GetNumberOfReads()); 
  if (!read_filter.empty()) {
    FilterReads(tmpname2, read_filter);
    reads_filename = tmpname2;
  }
  
  string cmd = gBlasrPath + "/blasr ";
  cmd += reads_filename;
  cmd += " ";
  cmd += tmpname1;
  cmd += " -sam -sdpTupleSize 8 -guidedAlignBandSize 100 -nCandidates 50 ";
  cmd += "-minMatch 11 ";
  cmd += kThreadsBlasr;
  cmd += " >";
  cmd += tmpname3;
  system(cmd.c_str());
  
  ifstream fi(tmpname3);
  string l;
  for (int i = 0; i < positions_.size(); i++) {
    positions_[i].clear();
  }
  positions_.resize(reads_num_);

  unordered_map<vector<int>, int> subpath_starts;
  unordered_set<vector<int>> dont_save;
  if (save_to_cache) {
    for (int i = 0; i < path.size(); i++) {
      vector<int> subpath;
      for (int j = i; j < path.size(); j++) {
        subpath.push_back(path[j]);
        int subpath_length = pathnodesposes[j] - pathnodesposesb[i]; 
        int first_length = pathnodesposes[i] - pathnodesposesb[i];
        if (aligment_cache_.count(subpath))
          dont_save.insert(subpath);
        else
          aligment_cache_[subpath].clear();
        subpath_starts[subpath] = i;
        if (subpath_length - first_length > max_read_len_) {
          break;
        }
      }
    }
  }
  set<int> rr;
  int good_out = 0, bad_out = 0, in_ok = 0, in_bad = 0;
  while (getline(fi, l)) {
    if (l[0] == '@') {
      continue;
    }
    PacbioAligmentData align = ParseAligment(l, seqall.length());
    assert(read_map_.count(align.name) > 0);
    int read_id = read_map_[align.name];
    int aligned_length = align.send - align.sstart;
    logdouble prob;
    for (int i = 2; i < 3; i++) {
      prob = AligmentProbability(seqall, read_seq_[read_id], align, i);
    }

    if (prob > GetMinReadProb(read_id) || true) {
      positions_[read_id].push_back(make_pair(align.tstart, prob));
      if (save_to_cache) {
        int it_begin = lower_bound(pathnodesposes.begin(), pathnodesposes.end(), 
                                   max(0, align.tstart - 5)) -
                       pathnodesposes.begin();
        int it_end = lower_bound(pathnodesposes.begin(), pathnodesposes.end(), 
                                 min(align.tstart + align.len + 5, (int)seq.length())) -
                     pathnodesposes.begin();
        assert(it_begin < path.size());
        assert(it_begin >= 0);
        assert(it_end < path.size());
        assert(it_end >= 0);
        vector<int> subpath(path.begin()+it_begin, path.begin()+it_end+1);
        int pos_begin = 0;
        if (it_begin > 0) {
          pos_begin = pathnodesposes[it_begin-1];
        }
        if (subpath_starts[subpath] == it_begin && dont_save.count(subpath) == 0) {
          assert(aligment_cache_.count(subpath));
          aligment_cache_[subpath].push_back(PacbioAligment(align.tstart - pos_begin, 
                                                            align.tend - pos_begin,
                                                            read_id,
                                                            prob));
        }
        rr.insert(read_id);
      }
    }
  }
  printf("rr %d\n", (int)rr.size());

  remove(tmpname1);
  remove(tmpname3);
  if (save_to_cache) {
    SaveAligments();
  }
  return positions_;
}

void PacbioReadSet::FilterReads(string out_filename, const unordered_set<int>& filter) {
  ifstream ifs(filename_);
  ofstream ofs(out_filename);

  string l1, l2, l3, l4;
  while (getline(ifs, l1)) {
    string nameline = l1.substr(1);
    vector<string> nameparts;
    split(nameparts, nameline, is_any_of(" \t"));
    string name = nameparts[0];
    int read_id = GetReadId(name);
    getline(ifs, l2);
    getline(ifs, l3);
    getline(ifs, l4);
    if (filter.count(read_id)) {
      ofs << l1 << endl;
      ofs << l2 << endl;
      ofs << l3 << endl;
      ofs << l4 << endl;
    }
  }
}

vector<vector<pair<int, logdouble> > >& PacbioReadSet::GetReadProbabilitiesAnchor(
    const Graph& gr, const vector<int>& path, int& total_len, int anchor) {
  assert(false);
  bool save_to_cache = true;
  printf("anch %d\n", save_to_cache);
  char tmpname1[L_tmpnam+6], tmpname2[L_tmpnam+6], tmpname3[L_tmpnam];
  tmpnam(tmpname1);
  strcat(tmpname1, ".fas");
  tmpnam(tmpname2);
  strcat(tmpname2, ".fq");
  tmpnam(tmpname3);
  printf("slow files %s %s %s\n", tmpname1, tmpname2, tmpname3);
  FILE *f = fopen(tmpname1, "w");
  int id = 0;
  fprintf(f, ">tmp\n");
  string seq = gr.nodes[path[0]]->s;
  vector<int> pathnodesposes, pathnodesposesb;
  pathnodesposes.push_back(seq.length());
  pathnodesposesb.push_back(0);
  for (int i = 1; i < path.size(); i++) {
    pathnodesposesb.push_back(seq.length());
    seq += gr.nodes[path[i]]->s;
    pathnodesposes.push_back(seq.length());
  }
  total_len = seq.length();
  printf("slow len %d\n", total_len);
  fprintf(f, "%s\n", seq.c_str());
  fclose(f);

  string seqrev = ReverseSeq(seq);
  string seqall = seq + kContigSeparator + seqrev;
  string reads_filename = filename_;
  if (anchors_cache_.count(path[anchor])) {
    FilterReads(tmpname2, anchors_cache_[path[anchor]]);
    reads_filename = tmpname2;
  }

  string cmd = gBlasrPath + "/blasr ";
  cmd += reads_filename;
  cmd += " ";
  cmd += tmpname1;
  cmd += " -sam -sdpTupleSize 8 -guidedAlignBandSize 100 -nCandidates 50 ";
  cmd += "-minMatch 11 ";
  cmd += kThreadsBlasr;
  cmd += " >";
  cmd += tmpname3;
  system(cmd.c_str());
  
  ifstream fi(tmpname3);
  string l;
  for (int i = 0; i < positions_.size(); i++) {
    positions_[i].clear();
  }
  positions_.resize(reads_num_);

  unordered_map<vector<int>, int> subpath_starts;
  if (save_to_cache) {
    for (int i = 0; i < path.size() && i <= anchor; i++) {
      vector<int> subpath;
      for (int j = i; j < path.size(); j++) {
        subpath.push_back(path[j]);
        int subpath_length = pathnodesposes[j] - pathnodesposesb[i]; 
        int first_length = pathnodesposes[i] - pathnodesposesb[i];
        aligment_cache_[subpath].clear();
        subpath_starts[subpath] = i;
        if (subpath_length - first_length > max_read_len_) {
          break;
        }
      }
    }
  }
  set<int> rr;
  while (getline(fi, l)) {
    if (l[0] == '@') {
      continue;
    }
    PacbioAligmentData align = ParseAligment(l, seqall.length());
    assert(read_map_.count(align.name) > 0);
    int read_id = read_map_[align.name];

    logdouble prob;
    for (int i = 2; i < 3; i++) {
      prob = AligmentProbability(seqall, read_seq_[read_id], align, i);
    }
    if (prob > GetMinReadProb(read_id) || true) {
      positions_[read_id].push_back(make_pair(align.tstart, prob));
      if (save_to_cache) {
        int it_begin = lower_bound(pathnodesposes.begin(), pathnodesposes.end(), 
                                   max(0, align.tstart - 5)) -
                       pathnodesposes.begin();
        int it_end = lower_bound(pathnodesposes.begin(), pathnodesposes.end(), 
                                 min(align.tstart + align.len + 5, (int)seq.length())) -
                     pathnodesposes.begin();
        assert(it_begin < path.size());
        assert(it_begin >= 0);
        assert(it_end < path.size());
        assert(it_end >= 0);
        if (it_begin <= anchor) {
          vector<int> subpath(path.begin()+it_begin, path.begin()+it_end+1);
          int pos_begin = 0;
          if (it_begin > 0) {
            pos_begin = pathnodesposes[it_begin-1];
          }
          assert(aligment_cache_.count(subpath));
          if (subpath_starts[subpath] == it_begin) {
            aligment_cache_[subpath].push_back(PacbioAligment(align.tstart - pos_begin,
                                                              align.tend - pos_begin,
                                                              read_id,
                                                              prob));
          }
          rr.insert(read_id);
        }
      }
    }
  }
  printf("rr %d\n", (int)rr.size());

  remove(tmpname1);
//  remove(tmpname3);
  if (save_to_cache) {
    SaveAligments();
  }
  return positions_;
}

PacbioReadSet::PacbioAligmentData PacbioReadSet::ParseAligment(
    const string& buf, int total_len, bool do_reverse) const {
  PacbioAligmentData ret;
  vector<string> parts;
  split(parts, buf, is_any_of("\t"));

  int lastsep = 0;
  for (int i = 0; i < parts[0].length(); i++) {
    if (parts[0][i] == '/') {
      lastsep = i;
    }
  }
  string name = parts[0].substr(0, lastsep);
  int posstart = atoi(parts[3].c_str());
  int flags = atoi(parts[1].c_str());
  int len = atoi(parts[8].c_str());
  int posend = posstart + len;
  int sstart = 0;
  int send = parts[9].length();
  int slen = parts[9].length();
  int edit_dist = 100000;

  for (int i = 11; i < parts.size(); i++) {
    if (parts[i][0] == 'X' && parts[i][1] == 'S') {
      sstart = atoi(parts[i].substr(5).c_str())-1;
    }
    if (parts[i][0] == 'X' && parts[i][1] == 'E') {
      send = atoi(parts[i].substr(5).c_str())-1;
    }
    if (parts[i][0] == 'X' && parts[i][1] == 'Q') {
      slen = atoi(parts[i].substr(5).c_str());
    }
    if (parts[i][0] == 'N' && parts[i][1] == 'M') {
      edit_dist = atoi(parts[i].substr(5).c_str());
    }
  }

  ret.tstart = posstart;
  ret.tend = posend;
  vector<pair<int, char> > cigar = ParseCigar(parts[5]);
  if ((flags & 16) && do_reverse) {
    int l = posend - posstart;
    posstart = total_len - posend;
    posend = posstart + l;
    reverse(cigar.begin(), cigar.end());
  } 
  if (send != slen) {
//    printf("add end %d\n", slen - send);
    cigar.push_back(make_pair(slen - send, 'I'));
/*    len += slen - send;
    posend += slen - send;
    assert(posend <= total_len);*/
  }
  if (sstart != 0) {
//    printf("add start %d\n", sstart);
    int match = min(sstart, posstart);
    int left = sstart - match;
    cigar.insert(cigar.begin(), make_pair(match, 'I'));
    if (left) {
      cigar.insert(cigar.begin(), make_pair(left, 'I'));
    }
/*    posstart -= match;
    assert(posstart >= 0);*/
  }
    
  ret.name = name;
  ret.flags = flags;
  ret.len = len;
  ret.posstart = posstart;
  ret.posend = posend;
  ret.sstart = sstart;
  ret.send = send;
  ret.slen = slen;
  ret.cigar = cigar;
  ret.edit_dist = edit_dist;
  return ret;
}

vector<pair<int, char> > PacbioReadSet::ParseCigar(const string& cigar) const {
  int start = 0;
  vector<pair<int, char> > ret;
  for (int i = 0; i < cigar.length(); i++) {
    if (cigar[i] < '0' || cigar[i] > '9') {
      if (cigar[i] == 'M' || cigar[i] == 'I' || cigar[i] == 'D') {
        int count = atoi(cigar.substr(start, i-start).c_str());
        ret.push_back(make_pair(count, cigar[i]));
        start = i+1;
      } else {
        printf("wtf %c\n", cigar[i]);
      }
    }
  }
  return ret;
}

void PositionsToReadProbsPacbio(
    int num_reads, const vector<vector<pair<int, logdouble> > >& positions,
    const PacbioReadSet& read_set, vector<logdouble>& read_probs) {
  read_probs.clear();
  read_probs.resize(num_reads);
  for(int i = 0; i < positions.size(); i++) {
    for (auto& y: positions[i]) {
      read_probs[i] += y.second;
    }
  }
}

void AddPositionsToReadProbsPacbio(
    const vector<vector<pair<pair<int, int>, logdouble> > >& positions,
    vector<logdouble>& read_probs) {
  for(int i = 0; i < positions.size(); i++) {
    for (auto& y: positions[i]) {
      read_probs[i] += y.second;
    }
  }
}

double GetTotalProbPacbio(const vector<logdouble>& read_probs, int total_len,
                          const PacbioReadSet& read_set, int& zero_reads,
                          double min_prob_per_base, double min_prob_start) {
  int total_c = 0;
  logdouble total_prob = 1;
  if (total_len == 0) {
    total_len = 1;
  }
  zero_reads = 0;
  FILE *f = fopen("rp.dat", "w");
  for (int i = 0; i < read_probs.size(); i++) {
    logdouble prob = read_probs[i];
    fprintf(f, "%s %lf\n", read_set.GetReadName(i).c_str(), prob.logval);
    logdouble mrp = logdouble(exp(min_prob_start)) *
                    (logdouble(exp(min_prob_per_base)) ^ read_set.GetReadLen(i));
    if (prob < mrp) {
      zero_reads++;
      prob = mrp;
    } else {
      double rq = prob.logval / read_set.GetReadLen(i);
    }
    total_prob *= prob;
    total_c++;
  }
  fclose(f);
  return total_prob.logval / total_c - log(2*total_len);
}

double GetTotalProbPacbio(const vector<logdouble>& read_probs, int total_len,
                          const PacbioReadSet& read_set, int& zero_reads,
                          double& non_zero_len) {
  int total_c = 0;
  logdouble total_prob = 1;
  if (total_len == 0) {
    total_len = 1;
  }
  zero_reads = 0;
  FILE *f = fopen("rp.dat", "w");
  non_zero_len = 0;
  for (int i = 0; i < read_probs.size(); i++) {
    logdouble prob = read_probs[i];
    fprintf(f, "%s %lf\n", read_set.GetReadName(i).c_str(), prob.logval);
    logdouble mrp = read_set.GetMinReadProb(i);
    if (prob < mrp) {
      zero_reads++;
      prob = mrp;
    } else {
      double rq = prob.logval / read_set.GetReadLen(i);
      non_zero_len += read_set.GetReadLen(i)*exp(rq);
    }
    total_prob *= prob;
    total_c++;
  }
  fclose(f);
  return total_prob.logval / total_c - log(2*total_len);
}

/*double CalcExactScoreForPacbio(const Graph& gr, vector<int> path, int kmer,
                          PacbioReadSet& read_set, int& zero_reads, int& total_len,
                          int ps, bool use_caching) {
  gr.NormalizePath(path);
  int tl2;
  vector<vector<pair<int, logdouble> > >& positions =
      read_set.GetExactReadProbabilities(gr, path, kmer, ps, total_len, tl2);
  vector<logdouble> read_probs;
  PositionsToReadProbsPacbio(read_set.GetNumberOfReads(), positions, read_set, read_probs);

  double nzl;
  double total_prob = GetTotalProbPacbio(read_probs, tl2, read_set, zero_reads, nzl);
  return total_prob;
}*/

/*double CalcScoreForPacbio(const Graph& gr, vector<int> path, int kmer,
                          PacbioReadSet& read_set, int& zero_reads, int& total_len, 
                          bool use_caching) {
  gr.NormalizePath(path);
  vector<vector<pair<int, logdouble> > >& positions =
      read_set.GetReadProbabilities(gr, path, kmer, total_len);
  vector<logdouble> read_probs;
  PositionsToReadProbsPacbio(read_set.GetNumberOfReads(), positions, read_set, read_probs);

  double nzl;
  double total_prob = GetTotalProbPacbio(read_probs, total_len, read_set, zero_reads, nzl);
  return total_prob;
}*/

/*double CalcScoreForPacbio2(const Graph& gr, vector<int> path, int kmer,
                           PacbioReadSet& read_set, int& zero_reads, int& total_len, 
                           bool use_caching) {
  gr.NormalizePath(path);
  vector<vector<pair<int, logdouble> > >& positions =
      read_set.GetReadProbabilities(gr, path, kmer, total_len);
  vector<logdouble> read_probs;
  PositionsToReadProbsPacbio(read_set.GetNumberOfReads(), positions, read_set, read_probs);

  double non_zero_len;
  double total_prob = GetTotalProbPacbio(read_probs, total_len, read_set, zero_reads,
                                         non_zero_len);
//  printf("zero reads %d/%d %lf\n", zero_reads, read_set.GetNumberOfReads(), 
//         1.*zero_reads/read_set.GetNumberOfReads());
  printf("total prob %lf\n", total_prob);
  return non_zero_len / total_len;
}*/

void InitReadProbs(int num_reads, vector<logdouble>& read_probs) {
  read_probs.clear();
  read_probs.resize(num_reads);
}

double CalcScoreForPacbio(const Graph& gr, vector<vector<int> > paths,
                          PacbioReadSet& read_set, int& zero_reads, int& total_len, 
                          bool use_caching, double no_cov_penalty,
                          double exp_cov_move,
                          double min_prob_per_base, double min_prob_start) {
  vector<logdouble> read_probs;
  InitReadProbs(read_set.GetNumberOfReads(), read_probs);
  total_len = 0;
  int st = 0;
  int bad_bases = 0;
  int bad_gaps = 0;
  int pn = 0;
  for (auto& path: paths) {
    gr.NormalizePath(path);
    vector<vector<int>> ctgs;
    vector<int> gaps;
    int last = 0;
    for (int i = 0; i < path.size(); i++) {
/*      if (path[i] < 0) {
        gaps.push_back(-path[i]);
        ctgs.push_back(vector<int>(path.begin()+last, path.begin()+i));
        last = i+1;
      }*/
    }
    ctgs.push_back(vector<int>(path.begin()+last, path.end()));
    for (int i = 0; i < ctgs.size(); i++) {
      vector<pair<int, int> > events;
      events.push_back(make_pair(-1000, 1));
      events.push_back(make_pair(2000, -3000));
      int tl;
      int pp = 0;
      for (int j = 0; j < ctgs[i].size(); j++) {
        if (ctgs[i][j] >= 0) {
          events.push_back(make_pair(pp, 1));
          int cl = gr.nodes[ctgs[i][j]]->s.length();
          events.push_back(make_pair(pp+cl, -cl));
          pp += cl;
        } else {
          pp += -ctgs[i][j];
        }
      }
      vector<vector<pair<pair<int, int>, logdouble> > >& positions =
          read_set.GetReadProbabilities(gr, ctgs[i], tl);
      for (int i = 0; i < positions.size(); i++) {
        for (auto &p: positions[i]) {
          if (p.second < read_set.GetMinReadProb(i)) continue;
          events.push_back(make_pair(p.first.first,
                                     1));
          events.push_back(make_pair(p.first.second,
                                     p.first.first - p.first.second));
        }
      }
      AddPositionsToReadProbsPacbio(positions, read_probs);
      total_len += tl;

      sort(events.begin(), events.end());
      multiset<int> inters;
      for (int j = 0; j < events.size(); j++) {
        if (events[j].second == 1) {
          inters.insert(events[j].first);
        }
        if (events[j].second != 1) {
          auto it = inters.find(events[j].first + events[j].second);
          inters.erase(it);
        }
        int good_start = tl-250;
        if (!inters.empty()) {
          int mm = *inters.begin();
          good_start = mm + exp_cov_move;
        }
        if (j + 1 < events.size()) {
          good_start = min(events[j+1].first, good_start);
        }
        good_start = min(good_start, tl - 250);
        if (good_start > max(2500, events[j].first)) {
          printf("ctg %d error %d-%d\n", pn, events[j].first, good_start);
          bad_bases += good_start - max(2500, events[j].first);
        }
         
      }
    }
    st += 1000000; pn++;
  }
  if (no_cov_penalty > 0) {
    printf("badp %d %d\n", bad_bases, bad_gaps);
  }

  double total_prob = GetTotalProbPacbio(read_probs, total_len, read_set, zero_reads,
                                         min_prob_per_base, min_prob_start);
  return total_prob - bad_bases*no_cov_penalty;
}
