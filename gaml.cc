#include <cstdio>
#include <string>
#include <fstream>
#include <cassert>
#include <boost/algorithm/string.hpp>
#include <map>
#include <cstdlib>
#include <ctime>
#include <queue>
#include <deque>
#include <set>
#include "graph.h"
#include "utility.h"
#include "input_output.h"
#include "moves.h"
#include "prob_calculator.h"
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

using namespace std;
using namespace boost;

string gBowtiePath = "../programs/bowtie2";
string gBlasrPath = "../programs/blasr/alignment/bin";

double ExtractDouble(const string& key, unordered_map<string, string>& cfg, double def) {
  if (cfg.count(key)) {
    return StringToDouble(cfg[key]);
  }
  return def;
}

int ExtractInt(const string& key, unordered_map<string, string>& cfg, int def) {
  if (cfg.count(key)) {
    return StringToInt(cfg[key]);
  }
  return def;
}

string ExtractString(const string &key, unordered_map<string, string>& cfg, string def) {
  if (cfg.count(key)) {
    return cfg["key"];
  }
  return def;
}

struct AssemblySettings {
 public:
  int threshold;
  string output_prefix;
  int max_iterations;
  bool do_postprocess;
  int extendadvp;
  int extendp;
  int breakp;
  int fixp;
  int localp;
  int fixlenp;
  double t0;
  AssemblySettings() {}
  AssemblySettings(unordered_map<string, string>& configs) {
    threshold = ExtractInt("long_contig_threshold", configs, 500);
    output_prefix = ExtractString("output_prefix", configs, "output");
    max_iterations = ExtractInt("max_iterations", configs, 50000);
    if (configs.count("do_proprocess")) {
      do_postprocess = true;
      max_iterations = 1;
    } else {
      do_postprocess = false;
    }
    extendadvp = ExtractInt("join_by_advice_p", configs, 20);
    extendp = ExtractInt("extend_p", configs, 10);
    breakp = ExtractInt("disconnect_p", configs, 10);
    fixp = ExtractInt("interchange_p", configs, 1);
    localp = ExtractInt("local_p", configs, 20);
    fixlenp = ExtractInt("fixlen_p", configs, 1);
    t0 = ExtractDouble("t0", configs, 0.008);
  }
};

void Optimize(Graph& gr, ProbCalculator& prob_calc, vector<vector<int>> paths,
    vector<pair<ReadSet*, ReadSet*>>& advice_paired,
    vector<PacbioReadSet*>& advice_pacbio,
    int longest_read, AssemblySettings& settings) {
  int threshold = settings.threshold;
  gr.CalcReachability();
  gr.CalcReachabilityBig(threshold);
  gr.CalcReachabilityLimit(2*longest_read);

  int total_len;
  int kmer = 47;

  vector<pair<int, int>> zeros;
  double cur_prob = prob_calc.CalcProb(paths, zeros, total_len);
  printf("start prob %lf len %d low prob reads", cur_prob, total_len);
  for (auto &e: zeros) {
    printf("%d/%d ", e.first, e.second);
  }
  printf("\n");
  OutputPathsToFile(paths, gr, kmer, threshold, settings.output_prefix);
  printf("\n");

  int itnum = 0;
  int start_len = total_len;
  double T = 5;
  double best_prob = cur_prob;
  vector<vector<int> > best_paths = paths;
  int num_start_paths = paths.size();
  while (true) {
    int clean = -1;
    unordered_map<int, vector<int> > locs;
    for (int i = 0; i < paths.size(); i++) {
      for (int j = 0; j < paths[i].size(); j++) {
        locs[paths[i][j]].push_back(i);
        locs[paths[i][j]^1].push_back(i);
      }
    }
    for (int i = 0; i < paths.size(); i++) {
      if (paths[i].size() > 1) {
        continue;
      }
      for (int j = 0; j < locs[paths[i][0]].size(); j++) {
        if (locs[paths[i][0]][j] != i) {
          clean = i;
        }
      }
    }
    if (clean == -1) {
      break;
    }
    printf("clean %d\n", clean);
    paths.erase(paths.begin() + clean);
  }

  map<int, int> last_rep;

  while (itnum <= settings.max_iterations) {
    // Start with checking repeated nodes and fix them
    
    vector<vector<int> > new_paths = paths;
    int extendadvp = settings.extendadvp;
    int extendp = settings.extendp;
    int breakp = settings.breakp;
    int fixp = settings.fixp;
    int localp = settings.localp;
    int fixlenp = settings.fixlenp;
    if (advice_pacbio.size() + advice_paired.size() == 0) {
      extendadvp = 0;
    }
    int r = rand() % (extendp + breakp + fixp + localp + extendadvp + fixlenp);
    bool was_local = false;
    bool was_break = false;
    int local_p, local_s, local_t;
    bool accept = false;
    bool force_best = false;
    
    if (settings.do_postprocess) {
      FixBigReps(new_paths, gr, threshold, true, prob_calc);
    } else {

      if (r < extendp) {
        if (!ExtendPaths(new_paths, gr, threshold, prob_calc)) {
          continue;
        }
      } else if (r < extendp + fixp) {
        if (!FixSomeBigReps(new_paths, gr, threshold, false, prob_calc)) {
          continue;
        }
      } else if (r < extendp + fixp + localp) {
        if (!LocalChange(new_paths, gr, threshold, local_p, local_s, local_t, prob_calc)) {
          continue;
        }
        if (local_p != -1) {
          was_local = true;
          printf("loc %d %d %d %d %d\n", new_paths[local_p][local_s], new_paths[local_p][local_t],
                 local_p, local_s, local_t);
        }
      } else if (r < extendp + fixp + localp + extendadvp) {
        printf("adv %d %d\n", advice_pacbio.size(), advice_paired.size());
        int r2 = rand() % (advice_pacbio.size() + advice_paired.size());
        if (r2 < advice_pacbio.size()) {
          PacbioReadSet* advice_set = advice_pacbio[rand()%advice_pacbio.size()];
          if (!ExtendPathsAdv(new_paths, gr, threshold, *advice_set, kmer, prob_calc)) {
            continue;
          }
        } else {
          pair<ReadSet*, ReadSet*> advice_set = advice_paired[rand()%advice_paired.size()];
          if (!ExtendPathsAdv(new_paths, gr, threshold, *advice_set.first, 
                              *advice_set.second, kmer, prob_calc)) {
            continue;
          }          
        }
        continue;
      } else if (r < extendp + fixp + localp + extendadvp + fixlenp) {
        if (!FixGapLength(new_paths, prob_calc)) {
          continue;
        }
      } else {
        if (!BreakPath(new_paths, gr, threshold)) {
          continue;
        }
        was_break = true;
      }
    }
    // Rep stats
    {
      unordered_map<int, int> counts;
      for (int i = 0; i < gr.nodes.size(); i+=2) {
        if (gr.nodes[i]->s.length() > threshold) {
          counts[i] = 0;
        }
      }
      for (int i = 0; i < paths.size(); i++) {
        for (int j = 0; j < paths[i].size(); j++) {
          if (paths[i][j] >= 0 && gr.nodes[paths[i][j]]->s.length() > threshold) {
            counts[(paths[i][j]/2)*2]++;
          }
        }
      }
      bool rep = false;
      for (auto &e: counts) {
        if (e.second > 1) {
          rep = true;
          printf("(%d: %dx %d) ", e.first, e.second, gr.nodes[e.first]->s.length());
        }
        if (e.second == 0) {
          new_paths.push_back(vector<int>({e.first}));
        }
      }
      if (rep) printf("\n");
    }

    while (true) {
      int clean = -1;
      unordered_map<int, vector<int> > locs;
      for (int i = 0; i < new_paths.size(); i++) {
        for (int j = 0; j < new_paths[i].size(); j++) {
          locs[new_paths[i][j]].push_back(i);
          locs[new_paths[i][j]^1].push_back(i);
        }
      }
      for (int i = 0; i < new_paths.size(); i++) {
        if (new_paths[i].size() > 1) {
          continue;
        }
        for (int j = 0; j < locs[new_paths[i][0]].size(); j++) {
          if (locs[new_paths[i][0]][j] != i) {
            clean = i;
          }
        }
      }
      if (clean == -1) {
        break;
      }
      if (was_local && clean < local_p) {
        local_p--;
      }
      printf("clean %d\n", clean);
      new_paths.erase(new_paths.begin() + clean);
    }
    itnum++;
    T = settings.t0 / log(itnum + 1);
    if (itnum % 100 == 0) {
      printf("cur best %lf: ", best_prob);
      OutputPathsToFile(best_paths, gr, kmer, threshold, settings.output_prefix);
      printf("\n");
    }
    OutputPathsToConsole(new_paths, gr, threshold);
    printf("\n");

    double new_prob = prob_calc.CalcProb(new_paths, zeros, total_len);

    if (new_prob > cur_prob || settings.do_postprocess) {
      if (was_local) {
        printf("local save\n");
        vector<int> pp;
        for (int i = local_s+1; i < local_t; i++) {
          pp.push_back(new_paths[local_p][i]);
        }
        int s = new_paths[local_p][local_s];
        int t = new_paths[local_p][local_t];
        printf("s t %d %d\n", s, t);
        if (gr.reach_big_[s].count(t)) {
          gr.reach_big_[s][t] = pp;
        }
        if (gr.reach_limit_[s].count(t)) {
          gr.reach_limit_[s][t] = pp;
        }
      }
      accept = true;
    } else if (was_break) {
      double prob = exp((new_prob - cur_prob) / T);
      uniform_real_distribution<double> dist(0.0, 1.0);
      double samp = dist(generator);
      if (samp < prob) {
        accept = true;
      }
    }
    if (accept) {
      printf("accept\n");
      cur_prob = new_prob;
      paths = new_paths;
    }
    if (new_prob > best_prob || force_best) {
      best_prob = new_prob;
      best_paths = new_paths;
    }
    time_t rawtime;
    struct tm * timeinfo;
    char buffer [80];

    time (&rawtime);
    timeinfo = localtime (&rawtime);

    strftime (buffer,80,"%H:%M:%S",timeinfo);
    printf("itnum %d temp %lf time %s new prob %lf %lf %lf len %d paths %d low prob reads ",
           itnum, T,
           buffer, new_prob,
           cur_prob, best_prob,
           total_len, new_paths.size());
    for (auto &e: zeros) {
      printf("%d/%d ", e.first, e.second);
    }
    printf("\n");
  }
  printf("cur best %lf: ", best_prob);
  OutputPathsToFile(best_paths, gr, kmer, threshold, settings.output_prefix);
  printf("\n");
}

struct Pos {
  int contig_pos;
  int node_pos;
  int node;
  int dist;
  vector<int> path;
  Pos() {}
  Pos(int cp, int np, int n, int d=0) : contig_pos(cp), node_pos(np), node(n), dist(d) {}
};

bool BaseEq(char a, char b) {
  if (a == b) return true;
  if (b == 'R' && (a == 'A' || a == 'G')) return true;
  if (b == 'Y' && (a == 'C' || a == 'T')) return true;
  if (b == 'K' && (a == 'G' || a == 'T')) return true;
  if (b == 'M' && (a == 'A' || a == 'C')) return true;
  if (b == 'S' && (a == 'C' || a == 'G')) return true;
  if (b == 'W' && (a == 'A' || a == 'T')) return true;
  return false;
}

bool AlignContig(const Graph& gr, int start, int target, 
                 const string& contig, vector<int>& path) {
  deque<Pos> fr;
  fr.push_back(Pos(0, gr.nodes[start]->s.length(), start));
  printf("0-1 BFS begin\n");
  int max_dist = 5;
  int mcp = 0;
  while (!fr.empty()) {
    Pos x = fr.front(); fr.pop_front();
    mcp = max(x.contig_pos, mcp);
    if (x.contig_pos > contig.length()) continue;
    if (x.dist < max_dist) {
      Pos nx(x.contig_pos+1, x.node_pos, x.node, x.dist+1);
      nx.path = x.path;
      fr.push_back(nx);
    }
    if (target == -1 && x.contig_pos == contig.length()) {
      path = x.path;
      return true;
    }
    if (x.node_pos == gr.nodes[x.node]->s.length()) {
      for (int j = 0; j < gr.nodes[x.node]->next.size(); j++) {
        int nnode = gr.nodes[x.node]->next[j]->id;
        if (nnode == target && x.contig_pos == contig.length()) {
          path = x.path;
          return true;
        }
        if (x.contig_pos >= contig.length()) continue;
        if (BaseEq(gr.nodes[nnode]->s[0], contig[x.contig_pos])) {
          Pos nx(x.contig_pos+1, 1, nnode, x.dist);
          nx.path = x.path;
          nx.path.push_back(nnode);
          fr.push_front(nx);
        } else if (x.dist < max_dist) {
          Pos nxs(x.contig_pos+1, 1, nnode, x.dist+1);
          nxs.path = x.path;
          nxs.path.push_back(nnode);
          fr.push_back(nxs);              
          Pos nxi(x.contig_pos, 1, nnode, x.dist+1);
          nxi.path = x.path;
          nxi.path.push_back(nnode);
          fr.push_back(nxi);              
        }
      }
    } else {
      if (x.contig_pos >= contig.length()) continue;
      if (BaseEq(gr.nodes[x.node]->s[x.node_pos], contig[x.contig_pos])) {
        Pos nx(x.contig_pos+1, x.node_pos+1, x.node, x.dist);
        nx.path = x.path;
        fr.push_front(nx);
      } else if (x.dist < max_dist) {
        Pos nxs(x.contig_pos+1, x.node_pos+1, x.node, x.dist+1);
        nxs.path = x.path;
        fr.push_back(nxs);
        Pos nxi(x.contig_pos, x.node_pos+1, x.node, x.dist+1);
        nxi.path = x.path;
        fr.push_back(nxi);
      }
    }
  }
  return false;
}

void AligmentToPath(
    const Graph&gr, const vector<pair<int, int>>& als, vector<vector<int>>& paths,
    const string& contig) {
  vector<int> cur_path;
  cur_path.push_back(als[0].second);
  int last = als[0].first + gr.nodes[als[0].second]->s.length();
  for (int i = 1; i < als.size(); i++) {
    int cur = als[i].first;
    printf("  last %d cur %d\n", last, cur);
    if (cur < last) {
      printf("PROBLEM\n");
      continue;
    }
    if (last < cur) {
      vector<pair<int, int>> runs;
      int current = 0;
      int beg = 0;
      for (int i = last; i < cur; i++) {
        if (contig[i] == 'N') {
          if (current == 0) beg = i;
          current += 1;
        } else {
          if (current > 4) {
            runs.push_back(make_pair(beg, i));
          }
          current = 0;
        }
      }
      if (current > 4) runs.push_back(make_pair(beg, cur));
      printf("runs %d\n", (int)runs.size());
      if (runs.size() > 1) {
        printf("wat %s\n", contig.substr(last, cur-last).c_str());
      }
      if (runs.empty()) {
        vector<int> found_path;
        bool found = AlignContig(gr, cur_path.back(), als[i].second, 
                                 contig.substr(last-1, cur-last), found_path);
        if (!found) {
          printf("not found %d %d %d\n", cur, last, cur-last);
          paths.push_back(cur_path);
          cur_path.clear();
        } else {
          printf("good found\n");
          cur_path.insert(cur_path.end(), found_path.begin(), found_path.end());
        }
      } else {
        cur_path.push_back(-(cur-last));
      }
    } 
    last = als[i].first + gr.nodes[als[i].second]->s.length();
    cur_path.push_back(als[i].second);
  }
  paths.push_back(cur_path);
}


void GetPaths(const Graph&gr, const string& contigs, vector<vector<int>>& paths) {
  unordered_map<string, string> ctgs;
  string buf = "";
  ifstream ifc(contigs);
  string lc;
  string last_name = "";
  while (getline(ifc, lc)) {
    if (lc[0] == '>') {
      if (!buf.empty()) {
        printf("add %s\n", last_name.c_str());
        ctgs[last_name] = buf;
      }
      buf = "";
      vector<string> np;
      string nn = lc.substr(1);
      split(np, nn, is_any_of(" "));
      last_name = np[0];
    } else {
      buf += lc;
    }
  }
  if (!buf.empty()) {
    printf("add %s\n", last_name.c_str());
    ctgs[last_name] = buf;
  }

  char tmpname1[L_tmpnam+6], tmpname2[L_tmpnam], tmpname3[L_tmpnam];
  tmpnam(tmpname1);
  strcat(tmpname1, ".fas");
  tmpnam(tmpname2);
  FILE *f = fopen(tmpname1, "w");
  for (int i = 0; i < gr.nodes.size(); i++) {
    if (gr.nodes[i]->s.length() >= 50) {
      fprintf(f, ">%d\n", i);
      fprintf(f, "%s\n", gr.nodes[i]->s.c_str());
    }
  }
  fclose(f);
  string cmd1 = "mummer/nucmer -maxmatch -p ";
  cmd1 += tmpname2;
  cmd1 += " "+ contigs + " " + tmpname1;
  printf("%s\n", cmd1.c_str());
  system(cmd1.c_str());
  string cmd2 = "mummer/show-coords ";
  cmd2 += tmpname2;
  cmd2 += ".delta >";
  cmd2 += tmpname2;
  cmd2 += ".coords";
  system(cmd2.c_str());

  ifstream fi(string(tmpname2)+".coords");
  string l;
  bool go = false;
  unordered_map<string, vector<pair<int, int>>> als;
  while (getline(fi, l)) {
    if (l[0] == '=') {
      go = true; continue;
    }
    if (!go) continue;
    vector<string> p;
    split(p, l, is_any_of(" \t"), token_compress_on);
    string contig = p[p.size()-2];
    int our_node = atoi(p[p.size()-1].c_str());
    double id = atof(p[p.size()-4].c_str());
    if (id < 99) continue;
    if (p[4] != "1" || atoi(p[5].c_str()) < gr.nodes[our_node]->s.length()-1) {
      printf("out %s: %s %s %d %d\n", contig.c_str(), p[4].c_str(), p[5].c_str(), our_node,
             (int)gr.nodes[our_node]->s.length());
      continue;
    }
    int place = atoi(p[1].c_str());
    als[contig].push_back(make_pair(place, our_node));
  }

  for (auto &e: als) {
    printf("%s:\n", e.first.c_str());
    sort(e.second.begin(), e.second.end());
    AligmentToPath(gr, e.second, paths, ctgs[e.first]);
    continue;
    vector<int> path;
    path.push_back(e.second[0].second);
    int last = e.second[0].first + gr.nodes[e.second[0].second]->s.length();
    for (int i = 1; i < e.second.size(); i++) {
      int cur = e.second[i].first;
      printf("  last %d cur %d\n", last, cur);
      if (cur < last) {
        printf("PROBLEM\n");
        continue;
      }
      while (last < cur) {
        int fit = -1;
        for (int j = 0; j < gr.nodes[path.back()]->next.size(); j++) {
          int nx = gr.nodes[path.back()]->next[j]->id;
          if (gr.nodes[nx]->s.length() > 500) continue;
          printf("  test:\n");
          printf("  %s\n", gr.nodes[nx]->s.c_str());
          printf("  %s\n", ctgs[e.first].substr(last-1, gr.nodes[nx]->s.length()).c_str());
          if (gr.nodes[nx]->s == ctgs[e.first].substr(last-1, gr.nodes[nx]->s.length())) {
            printf("good fit\n");
            fit = nx;
          }
        }
        if (fit == -1) {
          printf("no fit %s\n", ctgs[e.first].substr(last-1, cur-last).c_str());
          paths.push_back(path);
          path.clear();
          break;
        } else {
          path.push_back(fit);
          last += gr.nodes[fit]->s.length();
        }
      }
      last = e.second[i].first + gr.nodes[e.second[i].second]->s.length();
      path.push_back(e.second[i].second);
    }
    paths.push_back(path);
  }

  int alc = 0;
  for (auto &e: ctgs) {
    if (als.count(e.first) == 0) {
      printf("no al %s: %s\n", e.first.c_str(), e.second.c_str());
    } else {
      alc++;
    }
  }
  printf("paths size %d ctgs size %d %d\n", (int)paths.size(), (int)ctgs.size(), alc);
//  exit(1);
}

void ClipPaths(vector<vector<int>>& paths, const Graph& gr) {
  vector<vector<int>>out;
  for (auto &p: paths) {
    int b = -1, e = -1;
    for (int i = 0; i < p.size(); i++) {
      if (p[i] < 0) continue;
      if (gr.nodes[p[i]]->s.length() > 500) {
        e = i;
        if (b == -1) b = i;
      }
    }
    if (b == -1) continue;
    out.push_back(vector<int>(p.begin()+b, p.begin()+e+1));
  }
  paths=out;
}

bool ParseConfigLine(const string& line, string& key, string& value) {
  for (int i = 0; i < line.size(); i++) {
    if (line[i] == '=') {
      key = line.substr(0, i);
      value = line.substr(i+1);
      return true;
    }
  }
  return false;
}

bool LoadConfig(const string& config_file,
    unordered_map<string, string>& configs, 
    unordered_map<string, unordered_map<string, string>>& read_set_configs) {
  ifstream fi(config_file);
  if (fi.fail()) {
    printf("Failed to open config file\n");
    return false;
  }


  string current_read_set = "";

  string l;
  while (getline(fi, l)) {
    if (l.length() == 0) continue;
    if (l[0] == '[') {
      current_read_set = l.substr(1, l.length() - 2);
    } else if (l[0] >= 'a' && l[0] <= 'z') {
      string key, value;
      if (!ParseConfigLine(l, key, value)) {
        printf("Bad line in config file:\n%s\n", l.c_str());
        return false;
      }
      if (current_read_set.empty()) {
        configs[key] = value;
      } else {
        read_set_configs[current_read_set][key] = value;
      }
    }
  }

  return true;
}


void PrepareReadSetFromConfig(
    unordered_map<string, unordered_map<string, string>>& read_set_configs,
    vector<pair<SingleReadConfig, ReadSet*>>& single_reads,
    vector<pair<PairedReadConfig, pair<ReadSet*, ReadSet*>>>& paired_reads,
    vector<pair<SingleReadConfig, PacbioReadSet*>>& pacbio_reads) {
  for (auto &e: read_set_configs) {
    string cache_prefix = e.first;
    if (e.second.count("cache_prefix")) {
      cache_prefix = e.second["cache_prefix"];
    }

    if (e.second.count("type") == 0) {
      fprintf(stderr, "No type for read set %s, ignoring...\n", e.first.c_str());
      continue;
    }

    double weight = ExtractDouble("weight", e.second, 1);
    bool advice = false;
    if (e.second.count("advice")) {
      advice = true;
    }

    if (e.second["type"] == "single" || e.second["type"] == "pacbio") {
      if (e.second.count("filename") == 0) {
        fprintf(stderr, "Missing filename for read set %s, ignoring...\n",
                e.first.c_str());
        continue;
      }
      string filename = e.second["filename"];
      double mismatch_prob = ExtractDouble("mismatch_prob", e.second, 0.01);
      double match_prob = 1.0 - 4*mismatch_prob;
      double min_prob = ExtractDouble("min_prob_per_base", e.second, -0.7);
      double min_prob_start = ExtractDouble("min_prob_start", e.second, -10);

      double penalty_constant = ExtractDouble("penalty_constant", e.second, 0);
      double step = ExtractDouble("penalty_step", e.second, 50);
      SingleReadConfig cfg(penalty_constant, step, min_prob, min_prob_start, weight, advice);
      if (e.second["type"] == "single") {
        ReadSet* rs = new ReadSet(cache_prefix, filename, match_prob, mismatch_prob);
        single_reads.push_back(make_pair(cfg, rs)); 
      } else {
        PacbioReadSet* rs = new PacbioReadSet(cache_prefix, filename, match_prob,
                                              mismatch_prob);
        pacbio_reads.push_back(make_pair(cfg, rs));
      }
    } else if (e.second["type"] == "paired") {
      if (e.second.count("filename1") == 0) {
        fprintf(stderr, "Missing filename1 for read set %s, ignoring...\n",
                e.first.c_str());
        continue;
      }
      string filename1 = e.second["filename1"];
      if (e.second.count("filename2") == 0) {
        fprintf(stderr, "Missing filename2 for read set %s, ignoring...\n",
                e.first.c_str());
        continue;
      }
      string filename2 = e.second["filename2"];
      if (e.second.count("insert_mean") == 0) {
        fprintf(stderr, "Missing insert_mean for read set %s, ignoring...\n",
                e.first.c_str());
        continue;
      }
      if (e.second.count("insert_std") == 0) {
        fprintf(stderr, "Missing insert_std for read set %s, ignoring...\n",
                e.first.c_str());
        continue;
      }
      double insert_mean = StringToDouble(e.second["insert_mean"]);
      double insert_std = StringToDouble(e.second["insert_std"]);
      double mismatch_prob = ExtractDouble("mismatch_prob", e.second, 0.01);
      double match_prob = 1.0 - 4*mismatch_prob;
      double min_prob = ExtractDouble("min_prob_pre_base", e.second, -0.7);
      double min_prob_start = ExtractDouble("min_prob_start", e.second, -10);

      double penalty_constant = ExtractDouble("penalty_constant", e.second, 0);
      // TODO: fix this
      double step = insert_mean - ExtractDouble("penalty_step", e.second, 50);
      PairedReadConfig cfg(penalty_constant, step, insert_mean, insert_std, 
                           min_prob, min_prob_start, weight, advice);
      ReadSet* rs1 = new ReadSet(cache_prefix+"1", filename1, match_prob, mismatch_prob); 
      ReadSet* rs2 = new ReadSet(cache_prefix+"2", filename2, match_prob, mismatch_prob); 
      paired_reads.push_back(make_pair(cfg, make_pair(rs1, rs2)));
    } else {
      fprintf(stderr, "Unknown type %s for read set %s, ignoring...\n",
              e.second["type"].c_str(), e.first.c_str());
      continue;
    }
  }
}

template<class T, class TT>
void GetAdvice(const vector<pair<T, TT>>& input, vector<TT> output) {
  for (auto &e: input) {
    if (e.first.advice) {
      output.push_back(e.second);
    }
  }
}

void PrepareReads(
    vector<pair<SingleReadConfig, ReadSet*>>& single_reads,
    vector<pair<PairedReadConfig, pair<ReadSet*, ReadSet*>>>& paired_reads,
    vector<pair<SingleReadConfig, PacbioReadSet*>>& pacbio_reads,
    Graph& gr) {
  for (auto &e: pacbio_reads) {
    e.second->LoadAligments();
    e.second->PreprocessReads();
    e.second->NormalizeCache(gr);
    e.second->ComputeAnchors(gr);
  }

  for (auto &e: paired_reads) {
    e.second.first->LoadAligments();
    e.second.first->PreprocessReads();
    e.second.first->PrepareReadIndex();
    e.second.second->LoadAligments();
    e.second.second->PreprocessReads();
    e.second.second->PrepareReadIndex();
  }

  for (auto &e: single_reads) {
    e.second->LoadAligments();
    e.second->PreprocessReads();
    e.second->PrepareReadIndex();
  }
}

int GetLongestRead(
    vector<pair<SingleReadConfig, ReadSet*>>& single_reads,
    vector<pair<PairedReadConfig, pair<ReadSet*, ReadSet*>>>& paired_reads,
    vector<pair<SingleReadConfig, PacbioReadSet*>>& pacbio_reads) {
  int longest_read = 0;

  for (auto &e: single_reads) {
    for (int i = 0; i < e.second->GetNumberOfReads(); i++) {
      longest_read = max(longest_read, e.second->GetReadLen(i));
    }
  }
  
  for (auto &e: pacbio_reads) {
    for (int i = 0; i < e.second->GetNumberOfReads(); i++) {
      longest_read = max(longest_read, e.second->GetReadLen(i));
    }
  }

  for (auto &e: paired_reads) {
    longest_read = max(longest_read, (int)e.first.insert_mean);
  }
  return longest_read;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    printf("Missing config file!\nSyntax:\n./lama <config file>\n");
    return 1;
  }

  unordered_map<string, string> configs;
  unordered_map<string, unordered_map<string, string>> read_set_configs;

  if (!LoadConfig(argv[1], configs, read_set_configs)) {
    printf("Load config failed\n");
    return 1;
  }
  if (configs.count("graph") == 0) {
    fprintf(stderr, "Missing graph in config\n");
    return 1;
  }

  printf("%d %d\n", configs.size(), read_set_configs.size());
  vector<pair<SingleReadConfig, ReadSet*>> single_reads;
  vector<pair<PairedReadConfig, pair<ReadSet*, ReadSet*>>> paired_reads;
  vector<pair<SingleReadConfig, PacbioReadSet*>> pacbio_reads;
  PrepareReadSetFromConfig(read_set_configs, single_reads,
                           paired_reads, pacbio_reads);

  Graph gr;
  if (!LoadGraph(configs["graph"], gr)) {
    printf("Load graph failed\n");
    return 1;
  }

  ProbCalculator pc(single_reads, paired_reads, pacbio_reads, gr);

  vector<pair<ReadSet*, ReadSet*>> advice_paired;
  vector<PacbioReadSet*> advice_pacbio;
  GetAdvice(paired_reads, advice_paired);
  GetAdvice(pacbio_reads, advice_pacbio);
 
  PrepareReads(single_reads, paired_reads, pacbio_reads, gr);
  int longest_read = GetLongestRead(single_reads, paired_reads, pacbio_reads);

  //TODO: configure optimazation 

  vector<vector<int>> starting_paths;
  AssemblySettings settings(configs);

  if (configs.count("starting_assembly")) {
    GetPaths(gr, configs["starting_assembly"], starting_paths);
    ClipPaths(starting_paths, gr);
  } else {
    for (int i = 0; i < gr.nodes.size(); i+=2) {
      if (gr.nodes[i]->s.length() > settings.threshold)
        starting_paths.push_back(vector<int>({i}));
    }
  }

  Optimize(gr, pc, starting_paths, advice_paired, advice_pacbio, longest_read, settings); 
}


