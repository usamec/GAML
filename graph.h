#ifndef GRAPH_H__
#define GRAPH_H__

#include <string>
#include <boost/serialization/vector.hpp>
#include "unordered_map.hpp"
#include "unordered_set.hpp"
#include "logdouble.hpp"
#include <algorithm>
#include <random>
#include <cassert>

using namespace std;

extern const double kSmooth;
extern const char kContigSeparator;
extern default_random_engine generator;

typedef string Seq;

template <class T>
inline void hash_combine(std::size_t & seed, const T & v) {
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std {
  template<typename S, typename T> struct hash<pair<S, T>> {
    inline size_t operator()(const pair<S, T> & v) const {
      size_t seed = 0;
      ::hash_combine(seed, v.first);
      ::hash_combine(seed, v.second);
      return seed;
    }
  };

  template<typename T> struct hash<vector<T>> {
    inline size_t operator()(const vector<T>& v) const {
      size_t seed = 0;
      for (int i = 0; i < v.size(); i++) {
        ::hash_combine(seed, v[i]);
      }
      return seed;
    }
  };
}

inline int ConvertNodeId(int x) {
  if (x > 0)
    return 2*(x-1);
  else
    return 2*(-x-1)+1;
}

inline int InvertNode(int x) {
  return x ^ 1;
}
inline char ReverseBase(char a) {
  if (a == 'A') return 'T';
  if (a == 'C') return 'G';
  if (a == 'G') return 'C';
  if (a == 'T') return 'A';
}

inline Seq ReverseSeq(const Seq& x) {
  Seq ret;
  for (int i = x.length()-1; i >= 0; i--) {
    ret += ReverseBase(x[i]);
  }
  return ret;
}

class Node {
 public:
  int id;
  Seq s;
  vector<Node*> next;
  vector<double> next_prob;
  double next_sum;
  Node* inv;

  Seq GetNodeSeq(int kmer) const {
    if (inv->s.length() >= kmer-1) {
      Seq ss = inv->s.substr(inv->s.length() - kmer + 1);
      Seq retval = ReverseSeq(inv->s.substr(inv->s.length()-kmer+1)) + s;
      return retval;
    } else {
      Seq pred = ReverseSeq(inv->s);
      Node* cur = next[0];
      while (pred.length() < kmer-1) {
        pred += ReverseSeq(cur->inv->s);
        cur = cur->next[0];
      }
      pred = pred.substr(0, kmer-1);
      return pred + s;
    } 
  }

  int GetNodeLen(int kmer) const {
    return s.length() + kmer - 1;
  }

  void CalcProbSums() {
    next_sum = accumulate(next_prob.begin(), next_prob.end(), 0);
  }

  Node* SampleNext() const {
    if (next_prob.size() == 0) return NULL;

    uniform_real_distribution<double> dist(0.0, next_sum);
    double samp = dist(generator);
    double ss = 0;
    for (int i = 0; i < next_prob.size(); i++) {
      ss += next_prob[i];
      if (ss > samp || i == next_prob.size() - 1) {
        return next[i];
      }
    }
  }

  pair<Node*,double> SampleNextWithProb() const {
    if (next_prob.size() == 0) return make_pair((Node*)NULL, 0);

    uniform_real_distribution<double> dist(0.0, next_sum);
    double samp = dist(generator);
    double ss = 0;
    for (int i = 0; i < next_prob.size(); i++) {
      ss += next_prob[i];
      if (ss > samp || i == next_prob.size() - 1) {
        double prob = next_prob[i] / next_sum;
        return make_pair(next[i], prob);
      }
    }
  }

  // Precond - next.size() >= 2
  pair<Node*, double> SampleNextWithProbAndBan(int ban) const {
    double next_sum_ban = 0;
    for (int i = 0; i < next.size(); i++) {
      if (next[i]->id == ban) continue;
      next_sum_ban += next_prob[i];
    }
    uniform_real_distribution<double> dist(0.0, next_sum_ban);
    double samp = dist(generator);
    double ss = 0;
    for (int i = 0; i < next_prob.size(); i++) {
      if (next[i]->id == ban) continue;
      ss += next_prob[i];
      if (ss > samp || i == next_prob.size() - 1) {
        double prob = next_prob[i] / next_sum_ban;
        return make_pair(next[i], prob);
      }
    }
  }

  double GetNextProb(int next_id) const {
    for (int i = 0; i < next.size(); i++) {
      if (next[i]->id == next_id) {
        return next_prob[i] / next_sum; 
      }
    }
    assert(false);
    return 0;
  }

  double GetNextProbBan(int next_id, int ban) const {
    double next_sum_ban = 0;
    for (int i = 0; i < next.size(); i++) {
      if (next[i]->id == ban) continue;
      next_sum_ban += next_prob[i];
    }
    for (int i = 0; i < next.size(); i++) {
      if (next[i]->id == ban) continue;
      if (next[i]->id == next_id) {
        return next_prob[i] / next_sum_ban; 
      }
    }
    assert(false);
    return 0;
  }

  void InitProbs() {
    next_prob.clear();
    for (auto &x: next) {
      next_prob.push_back(kSmooth);
    }
  }

  // TODO: refactor!!!!
  void AddJump(int jump) {
    for (int i = 0; i < next.size(); i++) {
      if (next[i]->id == jump) {
        next_prob[i]++;
        return;
      }
    }
    assert(false);
  }

  bool HasNext(int next_id) {
    for (int i = 0; i < next.size(); i++) {
      if (next[i]->id == next_id) {
        return true;
      }
    }
    return false;
  }
};

struct Aligment {
  int position;
  int edit_dist;
  int read_id;
  int orientation;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
    ar & position;
    ar & edit_dist;
    ar & read_id;
    ar & orientation;
  }
  Aligment() {}
  Aligment(int pos, int ed, int read_id, int orientation) : 
        position(pos), edit_dist(ed), read_id(read_id), orientation(orientation) {}

  bool operator<(const Aligment& b) const {
    if (position == b.position) return read_id < b.read_id;
    return position < b.position;
  }
};

class Graph {
 public:
  vector<Node*> nodes;

  Node* operator[](int i) {
    return nodes[i];
  }

  vector<double> reachability_;
  double reach_sum_;
  vector<unordered_set<int> > reach_sets_;
  // from->to->path between
  vector<unordered_map<int, vector<int> > > reach_big_;
  vector<unordered_map<int, vector<int> > > reach_limit_;
  vector<vector<vector<int>>> reach_self_;

  vector<int> normalize_map;

  void CalcNormalizeMap() {
    normalize_map.resize(nodes.size());
    for (int i = 0; i < nodes.size(); i++) {
      normalize_map[i] = i;
    }
    for (int i = 0; i < nodes.size(); i++) {
      for (int j = 0; j < i; j++) {
        if (nodes[i]->s == nodes[j]->s) {
          normalize_map[i] = j;
          break;
        }
      }
    }
  }

  void NormalizePath(vector<int>& path) const {
    for (int i = 0; i < path.size(); i++) {
      if (path[i] >= 0)
        path[i] = normalize_map[path[i]];
    }
  }

  void CalcReachability();
  void CalcReachabilityBig(int threshold);
  void CalcReachabilityLimit(int max_dist);
  int SampleVertexByReach() const;

  void CalcProbSums() {
    for (auto &x: nodes) {
      x->CalcProbSums();
    }
  }

  void RecalculateProbsByPath(const vector<int>& path) {
    for (auto &x: nodes) {
      x->InitProbs();
    }

    for (int i = 1; i < path.size(); i++) {
      nodes[path[i-1]]->AddJump(path[i]);
      nodes[InvertNode(path[i])]->AddJump(InvertNode(path[i-1]));
    }

    CalcProbSums();
  }

  void OutputPath(const vector<int>& path, int kmer);
  void OutputPath(const vector<int>& path, int kmer, string filename);
  void OutputPathA(const vector<int>& path, int kmer, string filename, int cid);
  void OutputPathC(const vector<int>& path, int kmer, string filename, int cid);
  void OutputPathAT(const vector<int>& path, int kmer, string filename, int cid, int threshold);
 private:
};

bool LoadGraph(const string& filename, Graph& gr);

class ReadIndexTrivial {
 public:
  ReadIndexTrivial() {
    trans['A'] = 1;
    trans['T'] = 2;
    trans['C'] = 3;
    trans['G'] = 0;
  }
  void AddRead(const string& seq, int read_id);
  void GetReadCands(const string& seq, unordered_set<int>& read_cands);
  void PrintSizeInfo();
  unordered_map<unsigned long long, vector<int> > read_index_;
  char trans[256];
};

class ReadIndexMinHash {
 public:
  ReadIndexMinHash() {
    trans['A'] = 1;
    trans['T'] = 2;
    trans['C'] = 3;
    trans['G'] = 0;
  }
  void AddRead(const string& seq, int read_id);
  void GetMinHashWithPoses(const string& seq, vector<pair<unsigned long long, int>>& mhs);
  void GetReadCands(const string& seq, unordered_set<int>& read_cands);
  void GetReadCandsWithPoses(const string& seq, unordered_map<int, vector<int>>& read_cands);
  void PrintSizeInfo();
  unsigned long long Hash(unsigned long long x);
  unsigned long long GetMinHashForSeq(const string& seq);
  unordered_map<unsigned long long, vector<int> > read_index_;
  char trans[256];
  int read_len;
};

class ReadSet {
 public:
  // TODO: Calculate readlens from reads_file not from aligments
  ReadSet(const string& name, const string& filename, double match_prob, double mismatch_prob) : 
      save_changes_(0),
      reads_num_(0), name_(name), filename_(filename), match_prob_(match_prob),
      mismatch_prob_(mismatch_prob), load_success_(false), external_aligner_(false) {}

  void PreprocessReads();
  void PrepareReadIndex();

  // positions: read_id -> (position -> edit_dist, orientation)
  vector<vector<pair<int, pair<int, int> > > >& GetPositionsSlow(
      const Graph& gr, const vector<int>& path, int& total_len);
  // positions: read_id -> (position, (edit_dist, orientation))
  vector<vector<pair<int, pair<int, int> > > >& GetPositions(
      const Graph& gr, const vector<int>& path, int& total_len);
  vector<vector<pair<int, pair<int, int> > > >& AddPositions(
      const Graph& gr, const vector<int>& path, int& total_len, int st);
  vector<vector<pair<int, pair<int, int> > > >& GetPositions();

  void PrecomputeAlignmentForPaths(const vector<vector<int>>& paths, const Graph& gr);


  void ClearPositions();

  int GetNumberOfReads() const {
    return reads_num_;
  }

  int GetReadLen(int read_id) const {
    return read_lens_[read_id];
  
  }

  int save_changes_;

  void LoadAligments();
  void SaveAligments(bool force=false);
  
  double match_prob_;
  double mismatch_prob_;
  vector<double> match_probs_;
  vector<double> mismatch_probs_;
 private:
  const vector<Aligment>& GetAligmentForSubpath(
      const Graph& gr, const vector<int>& subpath);

  void PrecomputeAligmentForSubpaths(
      const Graph& gr, const vector<vector<int> >& subpaths);

  void AlignSubpathsInternal(
      const Graph& gr, const vector<vector<int> >& subpaths);

  int GetReadId(const string& read_name) {
    if (read_map_.count(read_name) == 0) {
      assert(load_success_ == false);
      int id = reads_num_;
      read_map_[read_name] = id;
      read_map_inv_[id] = read_name;
      reads_num_++;
      read_lens_.resize(reads_num_);
    }
    return read_map_[read_name];
  }

  void CalcMaxReadLen();

  void GetSubpathsFromPath(const vector<int>& path, const Graph& gr, unordered_set<vector<int>>& subpaths_precomp);

  int reads_num_;
  unordered_map<vector<int>, vector<Aligment> > aligment_cache_;
  unordered_map<string, int> read_map_;
  unordered_map<int, string> read_map_inv_;
  unordered_map<int, string> read_seqs_;
  vector<int> read_lens_;
  int max_read_len_;
  string name_;
  string filename_;
  bool load_success_;
  vector<vector<pair<int, pair<int, int> > > > positions_;
  ReadIndexMinHash read_index_;
  //ReadIndexTrivial read_index_;
  bool external_aligner_;
};

class PacbioReadSet {
 public:
  PacbioReadSet(const string& name, const string& filename, double match_prob, double mismatch_prob) : 
      save_changes_(0),
      reads_num_(0), name_(name), filename_(filename), match_prob_(match_prob),
      mismatch_prob_(mismatch_prob), min_match_prob_(1-2*(1-match_prob)), load_success_(false) {}

  int GetNumberOfReads() const {
    return reads_num_;
  }

  int GetReadLen(int read_id) const {
    return read_lens_[read_id];
  }

  void PreprocessReads();
  void ComputeAnchors(const Graph& gr);

  // read_id -> (position -> logprob)
  vector<vector<pair<int, logdouble> > >& GetReadProbabilitiesSlow(
      const Graph& gr, const vector<int>& path, int& total_len,
      bool save_to_cache=true);

  vector<vector<pair<int, logdouble> > >& GetReadProbabilitiesAnchor(
      const Graph& gr, const vector<int>& path, int& total_len,
      int anchor);

  vector<vector<pair<pair<int, int>, logdouble> > >& GetReadProbabilities(
      const Graph& gr, const vector<int>& path, int& total_len);

  vector<vector<pair<int, logdouble> > >& GetExactReadProbabilities(
      const Graph& gr, const vector<int>& path, int ps, int& total_len,
      int& total_len2);

  logdouble GetMinReadProb(int read_id) const {
    return (mismatch_prob_ ^ (read_lens_[read_id]*0.25)) *
           (match_prob_ ^ (read_lens_[read_id]*0.75));
  }
 
  string GetReadName(int read_id) const {
    auto it = read_map_inv_.find(read_id);
    return it->second;
  }

  void LoadAligments();
  void SaveAligments();
  void NormalizeCache(const Graph& gr);

  int GetMaxReadLen() const {
    return max_read_len_;
  }

  int GetGap(const Graph& gr, int first, int second, int read_id);
 private:
  int save_changes_;
  struct PacbioAligmentData {
    string name;
    int flags;
    int len;
    int posstart;
    int posend;
    int sstart;
    int send;
    int slen;
    int tstart;
    int tend;
    int edit_dist;
    vector<pair<int, char> > cigar;

    PacbioAligmentData() {}
  };

  struct PacbioAligment {
    int position;
    int position_end;
    int read_id;
    logdouble prob;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
      ar & position;
      ar & position_end;
      ar & prob.logval;
      ar & read_id;
    }
    PacbioAligment() {}
    PacbioAligment(int pos, int pos_end, int read_id, logdouble prob) : 
          position(pos), position_end(pos_end), read_id(read_id), prob(prob) {}

    bool operator<(const PacbioAligment& b) const {
      return position < b.position;
    }
  };

  void FilterReads(string out_filename, const unordered_set<int>& filter);

  int GetReadId(const string& read_name) {
    if (read_map_.count(read_name) == 0) {
      assert(load_success_ == false);
      int id = reads_num_;
      read_map_[read_name] = id;
      read_map_inv_[id] = read_name;
      reads_num_++;
      read_lens_.resize(reads_num_);
      read_seq_.resize(reads_num_);
    }
    return read_map_[read_name];
  }

  logdouble MatchProbability(const char c1, const char c2) const {
    if (c1 == kContigSeparator || c2 == kContigSeparator)
      return 0;

    if (c1 != c2) {
      return mismatch_prob_;
    } else {
      return match_prob_;
    }
  }

  void CalcMaxReadLen();

  PacbioAligmentData ParseAligment(const string& buf, int total_len, bool do_reverse=true) const;
  vector<pair<int, char> > ParseCigar(const string& cigar) const;
  logdouble AligmentProbability(
    const std::string &s1, const std::string &s2,
    const PacbioAligmentData& align_data, int band=2) const;
  int reads_num_;
  string name_;
  string filename_;
  logdouble match_prob_;
  logdouble mismatch_prob_;
  double min_match_prob_;
  bool load_success_;
  vector<int> read_lens_;
  int max_read_len_;
  unordered_map<string, int> read_map_;
  unordered_map<int, string> read_map_inv_;
  vector<vector<pair<int, logdouble> > > positions_;
  vector<vector<pair<pair<int, int>, logdouble> > > positions2_;
  vector<string> read_seq_;
  unordered_map<vector<int>, vector<PacbioAligment> > aligment_cache_;
 public:
  unordered_map<int, unordered_set<int> > anchors_cache_;
  unordered_map<int, unordered_set<int> > anchors_begin_;
  unordered_map<int, unordered_set<int> > anchors_end_;
  unordered_map<int, unordered_set<int> > anchors_reverse_;
};

double CalcScoreForPath(const Graph& gr, const vector<int>& path, int kmer,
                        ReadSet& read_set, bool use_caching = true);

double CalcScoreForPath(const Graph& gr, const vector<int>& path, int kmer,
                        ReadSet& read_set1, ReadSet& read_set2, 
                        double insert_mean, double insert_std,
                        bool use_caching = true);

double CalcScoreForPaths(const Graph& gr, const vector<vector<int>>& paths,
                         ReadSet& read_set1, ReadSet& read_set2, 
                         double insert_mean, double insert_std,
                         int& zero_reads, int& total_len,
                         bool use_caching = true,
                         double no_cov_penalty=0.0, double exp_cov_move=0.75,
                         bool use_all_to_cov=false,
                         double min_prob_per_base=-0.7, double min_prob_start=-10);


double CalcScoreForPaths(const Graph& gr, const vector<vector<int>>& paths,
                         ReadSet& read_set1, 
                         int& zero_reads, int& total_len,
                         bool use_caching = true,
                         double no_cov_penalty=0.0, double exp_cov_move=0.75,
                         double min_prob_per_base=-0.7, double min_prob_start=-10);


double CalcScoreForPacbio(const Graph& gr, vector<int> path,
                          PacbioReadSet& read_set, int& zero_reads,
                          int& total_len, bool use_caching = true);

double CalcExactScoreForPacbio(const Graph& gr, vector<int> path,
                               PacbioReadSet& read_set, int& zero_reads,
                               int& total_len, int ps, bool use_caching = true);

double CalcScoreForPacbio2(const Graph& gr, vector<int> path, int kmer,
                           PacbioReadSet& read_set, int& zero_reads,
                           int& total_len, 
                           bool use_caching = true);

double CalcScoreForPacbio(const Graph& gr, vector<vector<int> > paths,
                          PacbioReadSet& read_set, int& zero_reads,
                          int& total_len, bool use_caching = true,
                          double no_cov_penalty=0.0, double exp_cov_move=0.75,
                          double min_prob_per_base=-0.7, double min_prob_start=-10);


#endif
