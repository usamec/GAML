#include "moves.h"
#include <set>

bool BreakPath(vector<vector<int> >&new_paths, Graph& gr, int threshold) {
  printf("break\n");
  vector<pair<int, pair<int, int> > > options;
  for (int i = 0; i < new_paths.size(); i++) {
    if (new_paths[i].size() <= 1) continue;
    int last = -1;
    for (int j = 0; j < new_paths[i].size(); j++) {
      if (new_paths[i][j] >= 0 && gr.nodes[new_paths[i][j]]->s.length() > threshold) {
        if (last != -1) {
          options.push_back(make_pair(i, make_pair(last, j)));
        }
        last = j;
      }
    }
  }
  if (options.empty()) {
    return false;
  }

  int opt = rand()%options.size();
  int path_id = options[opt].first;
  vector<int> path = new_paths[path_id];
  bool has_minus = false;
  for (int i = options[opt].second.first; i <= options[opt].second.second; i++) {
    if (path[i] < 0)
      has_minus = true;
  }
  if (has_minus) {
    printf("break scaf\n");
  }
  new_paths.erase(new_paths.begin() + path_id);
  vector<int> path1(path.begin(), path.begin()+options[opt].second.first+1); 
  vector<int> path2(path.begin()+options[opt].second.second, path.end());
  new_paths.push_back(path1);
  new_paths.push_back(path2);

  return true;
}

bool LocalChange2(vector<vector<int> >& new_paths, Graph& gr, int threshold,
                  int path_id, int ps, int pt, ProbCalculator& prob_calc) {
  vector<int> path = new_paths[path_id];
  assert(gr.nodes[path[ps]]->s.length() > threshold);
  assert(gr.nodes[path[pt]]->s.length() > threshold);
  int elength = threshold;
  bool gap = false;
  for (int i = ps+1; i < pt; i++) {
    if (path[i] < 0) {
      elength += -path[i];
      gap = true;
    } else {
      elength += gr.nodes[path[i]]->s.length();
    }
  }
  printf("local 2 %d %d %d %d\n", ps, pt, path[ps], path[pt]);
  new_paths.erase(new_paths.begin() + path_id);
  new_paths.push_back(vector<int>(path.begin() + pt, path.end()));
  new_paths.push_back(vector<int>(path.begin(), path.begin() + ps + 1));

  int expect = path[pt];
  int max_extend = (pt - ps)*2;
  int total_added = 0;
  vector<int> last_path = new_paths.back();
  int start_size = last_path.size();
  while (last_path.back() != expect) {
    printf("lp %d ss %d me %d ta %d el %d exp %d\n", last_path.size(), start_size, max_extend, total_added,
        elength, expect);
    if ((last_path.size() > start_size + max_extend && gap == false) || total_added > 3*elength) {
      return false;
    }
    vector<vector<int> > cand_ends;
    vector<int> cand_add;
    for (int i = 0; i < 2; i++) {
      vector<int> cp = last_path;
      Node* next = NULL;
      int added_l = 0;
      while (true) {
        int fails = 0;
        while (true) {
          if (fails >= 20)
            return false;
          next = gr.nodes[cp.back()]->SampleNext();
          if (next == NULL) return false;
          fails++;
          if (next->s.length() > 2*elength && next->id != expect) {
            continue;
          }
          if (gr.reach_limit_[next->id].count(expect) || next->id == expect) {
            break;
          }
        }
        cp.push_back(next->id);
        if (next->id == expect) {
          break;
        }
        added_l += next->s.length();
        if (added_l > 200) {
          break;
        }
      }
      cand_ends.push_back(cp);
      cand_add.push_back(added_l);
    }
    vector<double> scores;
    for(int i = 0; i < cand_ends.size(); i++) {
      new_paths[new_paths.size()-1] = cand_ends[i];
      double score = prob_calc.CalcProb(new_paths);
      printf("ev %d: %lf\n", i, score);
      scores.push_back(score);
    }
    int best = 0;
    for (int i = 0; i < scores.size(); i++) {
      if (scores[i] > scores[best]) {
        best = i;
      }
    }
    last_path = cand_ends[best];
    total_added += cand_add[best];
    new_paths[new_paths.size()-1] = last_path;
  }
  assert(new_paths[new_paths.size()-1].back() == new_paths[new_paths.size()-2][0]);
  vector<int> op = new_paths.back();
  for (int i = 1; i < new_paths[new_paths.size()-2].size(); i++) {
    op.push_back(new_paths[new_paths.size()-2][i]);
  }
  new_paths[new_paths.size()-2] = op;
  new_paths.pop_back();
  return true;
}

bool FixMultiLocal(vector<vector<int> >& new_paths, Graph& gr, int threshold) {
  printf("fix multi local\n");
  int path_id = rand()%new_paths.size();
  vector<int> path = new_paths[path_id];
  unordered_map<int, vector<int> > poses;
  for (int i = 0; i < path.size(); i++) {
    if (path[i] < 0) continue;
    poses[path[i]].push_back(i);
  }
  vector<pair<int, pair<int, int>>> opts;
  for (auto &e: poses) {
    if (e.second.size() < 3) continue;
    for (int i = 2; i < e.second.size(); i++) {
      opts.push_back(make_pair(e.second[i-2], make_pair(e.second[i-1], e.second[i])));
    }
  }
  if (opts.size() == 0) return false;
  pair<int, pair<int, int>> opt = opts[rand()%opts.size()];
  vector<int> npath(path);
  int pp = opt.first;
  for (int i = opt.second.first; i < opt.second.second; i++, pp++) {
    npath[pp] = path[i];
  }
  for (int i = opt.first; i < opt.second.first; i++, pp++) {
    npath[pp] = path[i];
  }
  assert(pp == opt.second.second);
  new_paths[path_id] = npath;
  return true;
}

bool FixRep(vector<vector<int> >& new_paths, Graph& gr, int threshold) {
  printf("fix rep\n");
  int path_id = rand()%new_paths.size();
  vector<int> path = new_paths[path_id];
  unordered_map<int, vector<int> > poses;
  for (int i = 0; i < path.size(); i++) {
    if (path[i] < 0) continue;
    poses[path[i]].push_back(i);
  }
  vector<pair<int, int> > opts;
  for (auto &e: poses) {
    if (e.second.size() < 2) continue;
    for (int i = 1; i < e.second.size(); i++) {
      opts.push_back(make_pair(e.second[i-1], e.second[i]));
    }
  }
  if (opts.size() == 0) return false;
  pair<int, int> opt = opts[rand()%opts.size()];
  if (rand() % 4 == 0) { // double
    vector<int> path2(path.begin(), path.begin() + opt.second);
    path2.insert(path2.end(), path.begin() + opt.first, path.begin() + opt.second);
    path2.insert(path2.end(), path.begin() + opt.second, path.end());        
    new_paths[path_id] = path2;
    return true;
  } else {  // remove
    vector<int> path2(path.begin(), path.begin() + opt.first);
    path2.insert(path2.end(), path.begin() + opt.second, path.end());
    new_paths[path_id] = path2;
    return true;
  }
}

bool LocalChange(vector<vector<int> >& new_paths, Graph& gr, int threshold, int &path_id,
    int &xx, int &yy, ProbCalculator& prob_calc) {
/*  int r = rand() % 6;
  if (r <= 1) {
    path_id = -1;
    return FixRep(new_paths, gr, threshold);
  }
  if (r <= 3) {
    path_id = -1;
    return FixSelfLoops(new_paths, gr, threshold);
  }
  if (r <= 4) {
    path_id = -1;
    return FixMultiLocal(new_paths, gr, threshold);
  }*/
  vector<pair<int, pair<int, int> > > options;
  for (int i = 0; i < new_paths.size(); i++) {
    if (new_paths[i].size() <= 1) continue;
    vector<pair<int, int>> lp;
    int pos = 0;
    for (int j = 0; j < new_paths[i].size(); j++) {
      if (new_paths[i][j] >= 0 && gr.nodes[new_paths[i][j]]->s.length() > threshold) {
        lp.push_back(make_pair(pos, j));
      }
      if (new_paths[i][j] < 0) pos += -new_paths[i][j];
      else pos += gr.nodes[new_paths[i][j]]->s.length();
    }
/*    int last = -1;
    for (int j = 0; j < new_paths[i].size(); j++) {
      if (new_paths[i][j] >= 0 && gr.nodes[new_paths[i][j]]->s.length() > threshold) {
        if (last != -1 && j - last > 1) {
          options.push_back(make_pair(i, make_pair(last, j)));
        }
        last = j;
      }
    }*/
    for (int j = 1; j < lp.size(); j++) {
      options.push_back(make_pair(i, make_pair(lp[j-1].second, lp[j].second)));
      for (int k = 2; j-k >= 0; k++) {
        if (lp[j].first - lp[j-k].first < 5000) {
          options.push_back(make_pair(i, make_pair(lp[j-k].second, lp[j].second)));
        } else {
          break;
        }
      }
    }
  }
  if (options.empty()) {
    return false;
  }
  bool has_gap = false;
  int opt = rand() % options.size();
  printf("aaa %d %d %d %d\n", options[opt].first, options[opt].second.first,
         options[opt].second.second, new_paths.size());
  printf("%d\n", new_paths[options[opt].first].size());
  path_id = options[opt].first;
  printf("pick ");
  for (int i = options[opt].second.first; i <= options[opt].second.second; i++) {
    if (i > options[opt].second.first && new_paths[path_id][i-1] >= 0 && new_paths[path_id][i] >=0 ) {
      if (!gr.nodes[new_paths[path_id][i-1]]->HasNext(new_paths[path_id][i])) {
        printf("wtf %d %d: %d %d\n", path_id, i, new_paths[path_id][i-1], new_paths[path_id][i]);
      }
      assert(gr.nodes[new_paths[path_id][i-1]]->HasNext(new_paths[path_id][i]));
    }
    printf("%d", new_paths[path_id][i]);
    if (new_paths[path_id][i] >= 0)
      printf("(%d)", gr.nodes[new_paths[path_id][i]]->s.length());
    else
      has_gap = true;
    printf(" ");
  }
  printf("\n");
  if ((options[opt].second.second - options[opt].second.first > 7 || has_gap) && rand() % 2 <= 1) {
    path_id = -1;
    return LocalChange2(
        new_paths, gr, threshold, options[opt].first,
        options[opt].second.first, options[opt].second.second,
        prob_calc);
  }
  printf("local\n");
  vector<int> path = new_paths[path_id];
  int t = path[options[opt].second.second];
  int s = path[options[opt].second.first];
  xx = options[opt].second.first;
  vector<int> p2(path.begin(), path.begin()+options[opt].second.first+1);
  bool found = false;
  for (int extend = 0;
      extend < 2*(options[opt].second.second - options[opt].second.first + 1);
      extend++) {
    Node *next = NULL;
    int tries = 0;
    while (true) {
      tries++;
      if (tries > 100) return false;
      next = gr.nodes[p2.back()]->SampleNext();
      if (next == NULL) return false;
      if (gr.reach_limit_[next->id].count(t) || next->id == t) {
        break;
      }
    }
    if (next->id == t) {
      found = true;
      break;
    }
    p2.push_back(next->id);
  }
  if (!found) {
    return false;
  }
  yy = p2.size();
  for (int i = options[opt].second.second; i < path.size(); i++) {
    p2.push_back(path[i]);
  }
  new_paths[path_id] = p2;
  printf("done ");
  for (int i = xx; i <= yy; i++) {
    if (i > xx)
      assert(gr.nodes[new_paths[path_id][i-1]]->HasNext(new_paths[path_id][i]));
    printf("%d", new_paths[path_id][i]);
    if (new_paths[path_id][i] >= 0)
      printf("(%d)", gr.nodes[new_paths[path_id][i]]->s.length());
    printf(" ");
  }
  printf("\n");
  assert(new_paths[path_id][xx] == s);
  assert(new_paths[path_id][yy] == t);
  return true;
}

bool FixSelfLoops(vector<vector<int> >& new_paths, Graph& gr, int threshold) {
  printf("fix self\n");
  int path_id = rand()%new_paths.size();
  vector<int> path = new_paths[path_id];
/*  printf("path %d: ", path_id);
  for (auto &p: path)
    printf("%d ", p);
  printf("\n");*/
  vector<int> opts;
  for (int i = 0; i < path.size(); i++) {
    if (path[i] < 0) continue;
    if (gr.reach_self_[path[i]].size() > 0) {
//      printf("push %d ", path[i]);
      opts.push_back(i);
    }
  }
//  printf("\n");
  if (opts.size() == 0) return false;
  int opt = opts[rand()%opts.size()];
//  printf("try %d\n", path[opt]);
  vector<int> path2(path.begin(), path.begin()+opt);
  vector<int> ip = gr.reach_self_[path[opt]][rand()% gr.reach_self_[path[opt]].size()];
  path2.insert(path2.end(), ip.begin(), ip.end());
  path2.insert(path2.end(), path.begin()+opt, path.end());
/*  printf("new path %d: ", path_id);
  for (auto &p: path2)
    printf("%d ", p);
  printf("\n");*/
  new_paths[path_id] = path2;
  return true;
}

void ReversePath(vector<int>& path) {
  reverse(path.begin(), path.end());
  for (int i = 0; i < path.size(); i++) {
    if (path[i] >= 0) {
      path[i] ^= 1;
    }
  }
}

bool ExtendPathsAlt(vector<vector<int> >& paths, Graph& gr, int threshold) {
  for (int i = 0; i < paths.size(); i++) {
    if (rand() % 2 == 0) {
      reverse(paths[i].begin(), paths[i].end());
      for (int j = 0; j < paths[i].size(); j++) {
        if (paths[i][j] >= 0) {
          paths[i][j] ^= 1;
        }
      }
    }
  }

  int rp = rand()%paths.size();
  int rev = rand()%2;
  vector<int> path = paths[rp];
  paths.erase(paths.begin()+rp);
  if (rev) {
    reverse(path.begin(), path.end());
    for (int i = 0; i < path.size(); i++) {
      if (path[i] >= 0)
        path[i] ^= 1;
    }
  }
  unordered_map<int, vector<int> > path_ends;
  unordered_map<int, vector<pair<int, int> > > path_poses;
  for (int i = 0; i < paths.size(); i++) {
    path_ends[paths[i][0]].push_back(i+1);
    path_ends[paths[i].back()^1].push_back(-(i+1));
    for (int j = 1; j < paths[i].size() - 1; j++) {
      if (paths[i][j] >= 0 && gr.nodes[paths[i][j]]->s.length() > threshold) {
        path_poses[paths[i][j]].push_back(make_pair(i, j));
        path_poses[paths[i][j]^1].push_back(make_pair(i, j));
      }
    }
  }
  bool found = false;
  int join = 0;
  if (path_ends.count(path.back()) && path.size() > 1) {
    join = path_ends[path.back()][rand()%path_ends[path.back()].size()];
    found = true;
  }
  if (!found) {
    int add_length = 0;
    while (true) {
      vector<int> next_cand;                                                                    
      for (auto &e: gr.reach_big_[path.back()]) {                                               
        next_cand.push_back(e.first);                                                           
      }                                                                                         
      if (next_cand.empty() && add_length == 0) { 
        return false;
      }
      if (next_cand.empty()) {
        break;
      }
      int next = next_cand[rand()%next_cand.size()];
      int s = path.back();                                                                      
      for (int i = 0; i < gr.reach_big_[s][next].size(); i++) { 
        path.push_back(gr.reach_big_[s][next][i]);
        add_length += gr.nodes[path.back()]->s.length();
      }                                                                                         
      path.push_back(next);  
      add_length += gr.nodes[path.back()]->s.length();
      double p = exp(-add_length / 1000.0);
      uniform_real_distribution<double> dist(0.0, 1.0);
      double samp = dist(generator);
      if (samp > p) {
        break;
      }
    }
  }
  if (path_ends.count(path.back())) {
    join = path_ends[path.back()][rand()%path_ends[path.back()].size()];
    vector<int> join_path;
    int join_num;
    if (join < 0) {
      join_num = (-join)-1;
      join_path = paths[(-join)-1];
      reverse(join_path.begin(), join_path.end());
      for (int i = 0; i < join_path.size(); i++) {
        if (join_path[i] >= 0) 
          join_path[i] ^= 1;
      }
    } else {
      join_num = join-1;
      join_path = paths[join-1];
    }
    assert(path.back() == join_path[0]);
    for (int i = 1; i < join_path.size(); i++) {
      path.push_back(join_path[i]);
    }

    paths.erase(paths.begin() + join_num);
    paths.push_back(path);
    printf("alt extend early\n");
    return true;
  }
  if (path_poses[path.back()].empty()) {
    return false;
  }
  pair<int, int> pp = path_poses[path.back()][rand()%path_poses[path.back()].size()];
  if (paths[pp.first][pp.second] == path.back()) {
    vector<int> path2 = paths[pp.first];
    paths.erase(paths.begin() + pp.first);
    for (int i = pp.second + 1; i < path2.size(); i++) {
      path.push_back(path2[i]);
    }
    unordered_map<int, vector<int> > path_ends;
    unordered_map<int, vector<pair<int, int> > > path_poses;
    for (int i = 0; i < paths.size(); i++) {
      path_ends[paths[i][0]].push_back(i+1);
      path_ends[paths[i].back()^1].push_back(-(i+1));
      for (int j = 1; j < paths[i].size() - 1; j++) {
        if (paths[i][j] >= 0 && gr.nodes[paths[i][j]]->s.length() > threshold) {
          path_poses[paths[i][j]].push_back(make_pair(i, j));
          path_poses[paths[i][j]^1].push_back(make_pair(i, j));
        }
      }
    }
    
    path2.resize(pp.second+1);
    swap(path2, path);
    found = false;
    int join = 0;
    if (path_ends.count(path.back()) && path.size() > 1) {
      join = path_ends[path.back()][rand()%path_ends[path.back()].size()];
      found = true;
    }
    if (!found) {
      int add_length = 0;
      vector<int> path_zal = path;
      for (int tries = 0; tries < 5; tries++) {
        path = path_zal;
        while (true) {
          vector<int> next_cand;                                                                    
          for (auto &e: gr.reach_big_[path.back()]) {                                               
            next_cand.push_back(e.first);                                                           
          }                                                                                         
          if (next_cand.empty() && add_length == 0) { 
            return false;
          }
          if (next_cand.empty()) {
            break;
          }
          int next = next_cand[rand()%next_cand.size()];
          int s = path.back();                                                                      
          for (int i = 0; i < gr.reach_big_[s][next].size(); i++) { 
            path.push_back(gr.reach_big_[s][next][i]);
            add_length += gr.nodes[path.back()]->s.length();
          }                                                                                         
          path.push_back(next);  
          add_length += gr.nodes[path.back()]->s.length();
          double p = exp(-add_length / 1000.0);
          uniform_real_distribution<double> dist(0.0, 1.0);
          double samp = dist(generator);
          if (samp > p) {
            break;
          }
        }
      }
      if (path_ends.count(path.back())) {
        join = path_ends[path.back()][rand()%path_ends[path.back()].size()];
        vector<int> join_path;
        int join_num;
        if (join < 0) {
          join_num = (-join)-1;
          join_path = paths[(-join)-1];
          reverse(join_path.begin(), join_path.end());
          for (int i = 0; i < join_path.size(); i++) {
            if(join_path[i] >= 0)
              join_path[i] ^= 1;
          }
        } else {
          join_num = join-1;
          join_path = paths[join-1];
        }
        assert(path.back() == join_path[0]);
        for (int i = 1; i < join_path.size(); i++) {
          path.push_back(join_path[i]);
        }

        paths.erase(paths.begin() + join_num);
        paths.push_back(path);
        paths.push_back(path2);
        printf("2opt extend\n");
        return true;
      }
    }
    return false;
  } else {
    return false;
  }
  return false;
}

bool ExtendPaths(vector<vector<int> >& new_paths, Graph& gr, int threshold,
    ProbCalculator& prob_calc) {
  if (rand() % 7 == 0) {
    for (int i = 0; i < 5; i++) {
      vector<vector<int> > pp = new_paths;
      if (ExtendPathsAlt(pp, gr, threshold)) {
        new_paths = pp;
        return true;
      }
    }
    false;
  }
  printf("extend normal\n");
  bool found = false;
  int rp = rand()%new_paths.size();
  int rev = rand()%2;
  vector<int> path = new_paths[rp];
  int ps = path.size() - 1;
  if (rev == 1) {
    for (int i = 0; i < path.size(); i++) {
      if (path[i] >= 0)
        path[i] ^= 1;
    }
    reverse(path.begin(), path.end());
  }

  unordered_map<int, vector<int> > path_ends;
  unordered_map<int, vector<pair<int, int> > > path_poses;
  for (int i = 0; i < new_paths.size(); i++) {
    path_ends[new_paths[i][0]].push_back(i+1);
    path_ends[new_paths[i].back()^1].push_back(-(i+1));
    for (int j = 1; j < new_paths[i].size() - 1; j++) {
      if (new_paths[i][j] >= 0 && gr.nodes[new_paths[i][j]]->s.length() > threshold) {
        path_poses[new_paths[i][j]].push_back(make_pair(i, j));
        path_poses[new_paths[i][j]^1].push_back(make_pair(i, j));
      }
    }
  }

  int join = 0;
  pair<int, int> inner_join(-1, -1);
  if (path_ends.count(path.back()) && new_paths[rp].size() > 1) {
    join = path_ends[path.back()][rand()%path_ends[path.back()].size()];
    found = true;
  }
  if (!found) {
    int add_length = 0;
    while (true) {
      vector<int> next_cand;                                                                    
      for (auto &e: gr.reach_big_[path.back()]) {                                               
        next_cand.push_back(e.first);                                                           
      }                                                                                         
      if (next_cand.empty() && add_length == 0) { 
        return false;
      }
      if (next_cand.empty()) {
        break;
      }
      int next = next_cand[rand()%next_cand.size()];
      int s = path.back();                                                                      
      for (int i = 0; i < gr.reach_big_[s][next].size(); i++) { 
        path.push_back(gr.reach_big_[s][next][i]);
        add_length += gr.nodes[path.back()]->s.length();
      }                                                                                         
      path.push_back(next);  
      add_length += gr.nodes[path.back()]->s.length();
      double p = exp(-add_length / 1000.0);
      uniform_real_distribution<double> dist(0.0, 1.0);
      double samp = dist(generator);
      if (samp > p) {
        break;
      }
    }
    if (path_ends.count(path.back())) {
      join = path_ends[path.back()][rand()%path_ends[path.back()].size()];
      found = true;
    }
    if (rand() % 5 == 0) {
      found = true;
    }
  }
  if (!found) {
    return false;
  }
  int pt = path.size() - 1;

  /*      geometric_distribution<int> dist(0.2);
          int extend_limit = 2 + dist(generator);

          for (int extend = 0; extend < extend_limit; extend++) {
          if (path_poses.count(path.back())) {
          vector<pair<int, int> > choices;
          for (int i = 0; i < path_poses[path.back()].size(); i++) {
          if (path_poses[path.back()][i].first != rp) {
          choices.push_back(path_poses[path.back()][i]);
          }
          }
          if (!choices.empty()) {
          inner_join = choices[rand()%choices.size()];
          found = true;
          break;
          }
          }
          if (gr.nodes[path.back()]->s.length() > threshold && extend > 0) {
          break;
          }
          Node*next = gr.nodes[path.back()]->SampleNext();
          if (!next) {
          break;
          }
          path.push_back(next->id);
          }
          if (!found) {
          continue;
          }*/

  if (join != 0) {
    vector<int> join_path;
    int join_num;
    if (join < 0) {
      join_num = (-join)-1;
      join_path = new_paths[(-join)-1];
      reverse(join_path.begin(), join_path.end());
      for (int i = 0; i < join_path.size(); i++) {
        if (join_path[i] >= 0)
          join_path[i] ^= 1;
      }
    } else {
      join_num = join-1;
      join_path = new_paths[join-1];
    }
    assert(path.back() == join_path[0]);
    if (join_num != rp) {
      for (int i = 1; i < join_path.size(); i++) {
        path.push_back(join_path[i]);
      }
    } else {
      printf("repeat join\n");
    }

    new_paths.erase(new_paths.begin() + max(join_num, rp));
    if (join_num != rp) {
      new_paths.erase(new_paths.begin() + min(join_num, rp));
    }
    new_paths.push_back(path);
  } else {
    /*          if (path.back() == paths[inner_join.first][inner_join.second]) {
                for (int i = inner_join.second + 1; i < paths[inner_join.first].size(); i++) {
                path.push_back(paths[inner_join.first][i]);
                }
                new_paths[inner_join.first].resize(inner_join.second+1);
                new_paths.erase(new_paths.begin() + rp);
                new_paths.push_back(path);
                } else if (path.back()^1 == paths[inner_join.first][inner_join.second]) {
                for (int i = inner_join.second - 1; i >= 0; i--) {
                path.push_back(paths[inner_join.first][i]^1);
                }
                new_paths[inner_join.first].erase(
                new_paths[inner_join.first].begin(),
                new_paths[inner_join.first].begin() + inner_join.second);
                new_paths.erase(new_paths.begin() + rp);
                new_paths.push_back(path);
                } else {
                assert(false);
                }*/
    new_paths.erase(new_paths.begin() + rp);
    new_paths.push_back(path);
    //          printf("path join %d %d %d\n", rp, inner_join.first, inner_join.second);
  }
  if (!found) {
    return false;
  }
  vector<vector<int>> paths2 = new_paths;
  if (LocalChange2(paths2, gr, threshold, paths2.size() - 1, ps, pt, prob_calc)) {
    new_paths = paths2;
    printf("local 2 ok\n");
  } else {
    printf("local 2 fail\n");
  }
  return true;
}

int SamplePathByLength(vector<vector<int> >& paths, Graph& gr) {
  vector<int> lens(paths.size());
  int ss = 0;
  for (int i = 0; i < paths.size(); i++) {
    for (int j = 0; j < paths[i].size(); j++) {
      if (paths[i][j] >= 0) {
        lens[i] += gr.nodes[paths[i][j]]->s.length();
        ss += gr.nodes[paths[i][j]]->s.length();
      } else {
        lens[i] += -paths[i][j];
        ss += -paths[i][j];
      }
    }
  }
  int r = rand()%ss;
  ss = 0;
  for (int i = 0; i < lens.size(); i++) {
    ss += lens[i];
    if (r < ss)
      return i;
  }
  return paths.size()-1;
}

bool FixGapLength(vector<vector<int> >& paths, int path_id, int gap_pos,
                  ProbCalculator& prob_calc, int prev_len) {
  int cur_length = -paths[path_id][gap_pos];
  printf("fix len %d %d %d\n", path_id, gap_pos, cur_length); 
  assert(cur_length > 0);
  int max_move = cur_length - 1;
  if (prev_len != -1) {
    max_move = min(max_move, abs(prev_len - cur_length));
  }
  double cur_p = prob_calc.CalcProb(paths);
  for (int i = min(20, max_move); i >= 1; i--) {
    paths[path_id][gap_pos] = -(cur_length - i);
    double minus_p = prob_calc.CalcProb(paths);
    if (minus_p > cur_p) {
      return FixGapLength(paths, path_id, gap_pos, prob_calc, cur_length);
    }
    paths[path_id][gap_pos] = -(cur_length + i);
    double plus_p = prob_calc.CalcProb(paths);
    if (plus_p > cur_p) {
      return FixGapLength(paths, path_id, gap_pos, prob_calc, cur_length);
    }
    paths[path_id][gap_pos] = -(cur_length);
  }
  return true;
}

bool ExtendPathsAdv(vector<vector<int> >& paths, Graph&gr, int threshold,
                    PacbioReadSet& rs, int kmer, ProbCalculator& prob_calc) {
  printf("extend adv\n");
  int rp = SamplePathByLength(paths, gr);
  printf("rp %d\n", rp);
  vector<int> path = paths[rp];
  int rev = rand()%2;
  if (rev == 1) {
    for (int i = 0; i < path.size(); i++) {
      if (path[i] >= 0)
        path[i] ^= 1;
    }
    reverse(path.begin(), path.end());
  }
  paths.erase(paths.begin()+rp);

  int tl1;
  unordered_set<int> path_v(path.begin(), path.end());
  for (auto &e: path)
    path_v.insert(e^1);
  vector<pair<int, int> > cands;
  bool only_out = true;
  if (rand() % 5 == 0) only_out = false;
  bool allow_gaps = false;
  if (rand() % 5 == 0) allow_gaps = true;

  for (auto &r: rs.anchors_end_[path.back()]) {
    for (auto &x: rs.anchors_reverse_[r]) {
      if (gr.nodes[x]->s.length() > threshold)
        cands.push_back(make_pair(x, r));
    }
  }
  printf("cands len %d\n", cands.size());

  unordered_map<int, vector<int> > path_ends;
  for (int i = 0; i < paths.size(); i++) {
    path_ends[paths[i][0]].push_back(i+1);
    path_ends[paths[i].back()^1].push_back(-(i+1));
  }

  if (cands.empty()) return false;
  pair<int, int> cand = cands[rand()%cands.size()];
  int next = cand.first;
  bool gap = false;
  int gap_len = 0;
  if (gr.reach_limit_[path.back()].count(next) == 0) {
    gap = true;
  } else if (allow_gaps && rand() % 2 == 0) {
    gap = true;
  }
  if (gap) {
    gap_len = rs.GetGap(gr, path.back(), next, cand.second);
    printf("force gap %d\n", gap_len);
    if (gap_len < 0) return false;
  }

  int ps = path.size() - 1;
  int s = path.back();
  int gap_pos = -1;
  if (gap) {
    gap_pos = path.size();
    path.push_back(-gap_len);
    path.push_back(next);
  } else {
    for (int i = 0; i < gr.reach_limit_[s][next].size(); i++) { 
      path.push_back(gr.reach_limit_[s][next][i]);
    }                                                                                         
    path.push_back(next);
  }
  int pt = path.size() - 1;

  int join = 0;
  bool found = false;
  if (path_ends.count(path.back())) {
    join = path_ends[path.back()][rand()%path_ends[path.back()].size()];
    found = true;
  }
  if (rand() % 5 == 0) {
    found = true;
  }
  if (!found) {
    return false;
  }
  if (join != 0) {
    vector<int> join_path;
    int join_num;
    if (join < 0) {
      join_num = (-join)-1;
      join_path = paths[(-join)-1];
      reverse(join_path.begin(), join_path.end());
      for (int i = 0; i < join_path.size(); i++) {
        if (join_path[i] >= 0)
          join_path[i] ^= 1;
      }
    } else {
      join_num = join-1;
      join_path = paths[join-1];
    }
    assert(path.back() == join_path[0]);
    for (int i = 1; i < join_path.size(); i++) {
      path.push_back(join_path[i]);
    }

    paths.erase(paths.begin() + join_num);
    paths.push_back(path);
  } else {
    paths.push_back(path);
  }
  vector<vector<int> > paths2 = paths;
  if (gap) {
    printf("adv done gap\n");
/*    FixGapLength(paths, paths.size() - 1, gap_pos, prob_calc, -1);  
    if (paths.back()[gap_pos] == -1) {
      printf("bad -1 gap\n");
      return false;
    }*/
  } else {
    printf("adv done normal\n");
    if (LocalChange2(paths2, gr, threshold, paths.size() - 1, ps, pt, prob_calc)) {
      paths = paths2;
      printf("local 2 ok\n");
    } else {
      printf("local 2 fail\n");
    }
  }
  return true;
}

bool ExtendPathsAdv(vector<vector<int> >& paths, Graph&gr, int threshold,
                    ReadSet& rs1, ReadSet& rs2, int kmer, ProbCalculator& prob_calc) {
  printf("extend adv\n");
  int rp = SamplePathByLength(paths, gr);
  printf("rp %d\n", rp);
  vector<int> path = paths[rp];
  int rev = rand()%2;
  if (rev == 1) {
    for (int i = 0; i < path.size(); i++) {
      if (path[i] >= 0)
        path[i] ^= 1;
    }
    reverse(path.begin(), path.end());
  }
  paths.erase(paths.begin()+rp);
  unordered_map<int, vector<int> > read_poses;
  unordered_map<int, vector<int> > read_poses_1;
  for (int i = 0; i < gr.nodes.size(); i++) {
    if (gr.nodes[i]->s.length() > threshold) {
      vector<int> path({i});
      int tl1;
      vector<vector<pair<int, pair<int, int> > > >& positions1 = 
        rs2.GetPositions(gr, path, tl1);
      for (int j = 0; j < rs2.GetNumberOfReads(); j++) {
        if (!positions1[j].empty()) {
          read_poses[j].push_back(i);
          if (positions1[j][0].second.second == 1)
            read_poses_1[j].push_back(i);
        }
      }
    }
  }

  int tl1;
  unordered_set<int> path_v(path.begin(), path.end());
  for (auto &e: path)
    path_v.insert(e^1);
  vector<vector<pair<int, pair<int, int> > > >& positions1 = 
    rs1.GetPositions(gr, path, tl1);
  vector<int> cands;
  bool only_out = true;
  if (rand() % 5 == 0) only_out = false;
  bool allow_gaps = false;
  if (rand() % 5 == 0) allow_gaps = true;

  for (int i = 0; i < rs1.GetNumberOfReads(); i++) {
    if (positions1[i].empty()) continue;
    if (positions1[i][0].second.second != 0) continue;
    for (int j = 0; j < read_poses_1[i].size(); j++) {
      if (path_v.count(read_poses_1[i][j]) && only_out) continue;
      if (gr.reach_limit_[path.back()].count(read_poses_1[i][j]) != 0 || allow_gaps) {
        cands.push_back(read_poses_1[i][j]);
      }
    } 
  }
  unordered_map<int, vector<int> > path_ends;
  for (int i = 0; i < paths.size(); i++) {
    path_ends[paths[i][0]].push_back(i+1);
    path_ends[paths[i].back()^1].push_back(-(i+1));
  }

  if (cands.empty()) return false;
  int next = cands[rand()%cands.size()];
  bool gap = false;
  if (gr.reach_limit_[path.back()].count(next) == 0) {
    gap = true;
  } else if (allow_gaps && rand() % 2 == 0) {
    printf("force gap\n");
    gap = true;
  }

  int ps = path.size() - 1;
  int s = path.back();
  int gap_pos = -1;
  if (gap) {
    gap_pos = path.size();
    path.push_back(-21);
    path.push_back(next);
  } else {
    for (int i = 0; i < gr.reach_limit_[s][next].size(); i++) { 
      path.push_back(gr.reach_limit_[s][next][i]);
    }                                                                                         
    path.push_back(next);
  }
  int pt = path.size() - 1;

  int join = 0;
  bool found = false;
  if (path_ends.count(path.back())) {
    join = path_ends[path.back()][rand()%path_ends[path.back()].size()];
    found = true;
  }
  if (rand() % 5 == 0) {
    found = true;
  }
  if (!found) {
    return false;
  }
  if (join != 0) {
    vector<int> join_path;
    int join_num;
    if (join < 0) {
      join_num = (-join)-1;
      join_path = paths[(-join)-1];
      reverse(join_path.begin(), join_path.end());
      for (int i = 0; i < join_path.size(); i++) {
        if (join_path[i] >= 0)
          join_path[i] ^= 1;
      }
    } else {
      join_num = join-1;
      join_path = paths[join-1];
    }
    assert(path.back() == join_path[0]);
    for (int i = 1; i < join_path.size(); i++) {
      path.push_back(join_path[i]);
    }

    paths.erase(paths.begin() + join_num);
    paths.push_back(path);
  } else {
    paths.push_back(path);
  }
  vector<vector<int> > paths2 = paths;
  if (gap) {
    printf("adv done gap\n");
    FixGapLength(paths, paths.size() - 1, gap_pos, prob_calc, -1);  
    if (paths.back()[gap_pos] == -1) {
      printf("bad -1 gap\n");
      return false;
    }
  } else {
    printf("adv done normal\n");
    if (LocalChange2(paths2, gr, threshold, paths.size() - 1, ps, pt, prob_calc)) {
      paths = paths2;
      printf("local 2 ok\n");
    } else {
      printf("local 2 fail\n");
    }
  }
  return true;
}

bool FixGapLength(vector<vector<int> >& paths, ProbCalculator& prob_calc) {
  vector<pair<int, int> > opts;
  for (int i = 0; i < paths.size(); i++) {
    for (int j = 0; j < paths[i].size(); j++) {
      if (paths[i][j] < 0) {
        opts.push_back(make_pair(i, j));
      }
    }
  }
  if (opts.empty()) return false;
  pair<int, int> opt = opts[rand()%opts.size()];
  return FixGapLength(paths, opt.first, opt.second, prob_calc, -1);
}

bool SplitOnNode(int node, vector<vector<int>>& paths) {
  vector<vector<int>>paths2(paths);
  vector<vector<int>> paths_with_node;
  for (int i = paths2.size()-1; i>=0; i--) {
    bool has_node = false;
    for (int j = 0; j < paths2[i].size(); j++) {
      if ((paths2[i][j]/2)*2 == node) {
        has_node = true;
        break;
      }
    }
    if (has_node) {
      paths_with_node.push_back(paths2[i]);
      swap(paths2[paths2.size()-1], paths2[i]);
      paths2.pop_back();
    }
  }

  for (auto &p: paths_with_node) {
    int last = 0;
    for (int i = 1; i < p.size(); i++) {
      if ((p[i]/2)*2 == node) {
        vector<int> pp(p.begin()+last, p.begin()+i+1);
        paths2.push_back(pp);
        last = i;
      }
    }
    if (last != p.size()-1) {
      vector<int> pp(p.begin()+last, p.end());
      paths2.push_back(pp);      
    }
  }
  paths = paths2;
}

void FixRepForNode2(vector<vector<int>>& paths, Graph&gr, int threshold, bool disjoin_similar,
                   int node, ProbCalculator& prob_calc) {
  vector<pair<int, int> > poses;
  vector<pair<int, pair<int, int>>> doubles;
  vector<pair<int, pair<int, int>>> pals;
  for (int i = 0; i < paths.size(); i++) {
    int lp = -1;
    vector<int> cur_poses;
    for (int j = 0; j < paths[i].size(); j++) {
      if (paths[i][j] < 0) continue;
      if ((paths[i][j]/2)*2 == node) {
        poses.push_back(make_pair(i, j));
        if (lp != -1 && paths[i][j] == paths[i][lp]) {
          doubles.push_back(make_pair(i, make_pair(lp, j)));
        }
        lp = j;
        for (int k = 0; k < cur_poses.size(); k++) {
          if (paths[i][j] != paths[i][cur_poses[k]]) {
            pals.push_back(make_pair(i, make_pair(cur_poses[k], j)));
          }
        }
        cur_poses.push_back(j);
      }
    }
  }
  bool better = false;
  printf("fix rep node2 %d %d %d %d\n", node, gr.nodes[node]->s.length(), poses.size(), doubles.size());
  double cur_score = prob_calc.CalcProb(paths);
  set<pair<int, int> > disjoint;
  for (int i = 0; i < poses.size(); i++) {
    for (int j = 0; j < i; j++) {
      vector<vector<int>> paths2(paths);
      if (poses[i].first != poses[j].first) {
        vector<int> p1 = paths[poses[i].first];
        vector<int> p2 = paths[poses[j].first];
        vector<int> pp1, pp2;
        if (p1[poses[i].second] == p2[poses[j].second]) {
          pp1 = vector<int>(p1.begin(), p1.begin()+poses[i].second);
          pp1.insert(pp1.end(), p2.begin()+poses[j].second, p2.end());
          pp2 = vector<int>(p2.begin(), p2.begin()+poses[j].second);
          pp2.insert(pp2.end(), p1.begin()+poses[i].second, p1.end());
        } else {
          vector<int> s1(p1.begin(), p1.begin()+poses[i].second+1);
          vector<int> e1(p1.begin()+poses[i].second+1, p1.end());
          vector<int> s2(p2.begin(), p2.begin()+poses[j].second);
          vector<int> e2(p2.begin()+poses[j].second, p2.end());
          ReversePath(s2); ReversePath(e2);
          pp1.insert(pp1.end(), s1.begin(), s1.end());
          pp1.insert(pp1.end(), s2.begin(), s2.end());
          pp2.insert(pp2.end(), e2.begin(), e2.end());
          pp2.insert(pp2.end(), e1.begin(), e1.end());
        }
        paths2[poses[i].first] = pp1;
        paths2[poses[j].first] = pp2;
        if (paths2[max(poses[i].first, poses[j].first)].size() <= 1) {
          paths2.erase(paths2.begin()+max(poses[i].first, poses[j].first));
        }
        if (paths2[min(poses[i].first, poses[j].first)].size() <= 1) {
          paths2.erase(paths2.begin()+min(poses[i].first, poses[j].first));
        }
        double score = prob_calc.CalcProb(paths2);
        printf("scr %lf %lf\n", score, cur_score);
        if (fabs(score - cur_score) < 0.001) {
          if (disjoin_similar) {
            disjoint.insert(poses[i]);
            disjoint.insert(poses[j]);
          }
        }
        if (score > cur_score) {
          paths = paths2;
          FixRepForNode2(paths, gr, threshold, disjoin_similar, node, prob_calc);
          return;
        }
      }
    }
  }
  for (int i = 0; i < poses.size(); i++) {
    for (int j = 0; j < doubles.size(); j++) {
      vector<vector<int>> paths2(paths);
      if (poses[i].first != doubles[j].first) {
        printf("dj %d %d\n", poses[i].first, doubles[j].first);
        vector<int> p1(paths[poses[i].first].begin(), paths[poses[i].first].begin()+poses[i].second);
        vector<int> p2(paths[doubles[j].first].begin(), paths[doubles[j].first].begin()+doubles[j].second.first);
        p2.insert(p2.end(), paths[doubles[j].first].begin()+doubles[j].second.second,
                  paths[doubles[j].first].end());
        vector<int> pj(paths[doubles[j].first].begin()+doubles[j].second.first,
                       paths[doubles[j].first].begin()+doubles[j].second.second+1);
        if (pj[0] != paths[poses[i].first][poses[i].second]) {
          ReversePath(pj);
        }
        p1.insert(p1.end(), pj.begin(), pj.end());
        p1.insert(p1.end(), paths[poses[i].first].begin()+poses[i].second+1,
            paths[poses[i].first].end());
        paths2[poses[i].first] = p1;
        paths2[doubles[j].first] = p2;
      } else {
        printf("hdj %d\n", poses[i].first);
        vector<int> pj(paths[doubles[j].first].begin()+doubles[j].second.first,
                       paths[doubles[j].first].begin()+doubles[j].second.second);
        if (pj[0] != paths[poses[i].first][poses[i].second]) {
          ReversePath(pj);
          pj.insert(pj.begin(), pj.back());
          pj.pop_back();
        }
        if (poses[i].second < doubles[j].second.first) {
          vector<int> p1(paths[poses[i].first]);
          p1.erase(p1.begin()+doubles[j].second.first,
                   p1.begin()+doubles[j].second.second);
          p1.insert(p1.begin()+poses[i].second,
                    pj.begin(), pj.end());
          paths2[poses[i].first] = p1;
        } else if (poses[i].second > doubles[j].second.second) {
          vector<int> p1(paths[poses[i].first]);
          p1.insert(p1.begin()+poses[i].second,
                    pj.begin(), pj.end());
          p1.erase(p1.begin()+doubles[j].second.first,
                   p1.begin()+doubles[j].second.second);
          paths2[poses[i].first] = p1;
        } else {
          continue;
        }
      }
/*        printf("before:");
        printf("(");
        for (int k = 0; k < paths[poses[i].first].size(); k++) printf("%d ", paths[poses[i].first][k]);
        printf(") (");
        for (int k = 0; k < paths[doubles[j].first].size(); k++) printf("%d ", paths[doubles[j].first][k]);
        printf(")\nafter: (");
        for (int k = 0; k < p1.size(); k++) printf("%d ", p1[k]);
        printf(") (");
        for (int k = 0; k < p2.size(); k++) printf("%d ", p2[k]);
        printf(")\n\n");*/
      if (paths2[doubles[j].first].size() <= 1) {
        paths2.erase(paths2.begin()+doubles[j].first);
      }
      double score = prob_calc.CalcProb(paths2);
      printf("scr %lf %lf %d %d\n", score, cur_score, (int)paths.size(),
          (int)paths2.size());
      if (fabs(score - cur_score) < 0.002) {
        if (disjoin_similar) {
          disjoint.insert(poses[i]);
          disjoint.insert(make_pair(doubles[j].first, doubles[j].second.first));
          disjoint.insert(make_pair(doubles[j].first, doubles[j].second.second));
        }
        printf("similar\n");
      }
      if (score > cur_score) {
        paths = paths2;
        FixRepForNode2(paths, gr, threshold, disjoin_similar, node, prob_calc);
        return;
      }
    }
  }
  for (auto &pal: pals) {
    vector<vector<int>> paths2(paths);
    vector<int> p(paths2[pal.first].begin()+pal.second.first,
                  paths2[pal.first].begin()+pal.second.second+1);
    ReversePath(p);
    for (int i = 0; i < p.size(); i++) {
      paths2[pal.first][pal.second.first+i] = p[i];
    }
    double score = prob_calc.CalcProb(paths2);
    printf("scr %lf %lf %d %d\n", score, cur_score, (int)paths.size(),
        (int)paths2.size());
    if (fabs(score - cur_score) < 0.002) {
      if (disjoin_similar) {
        disjoint.insert(make_pair(pal.first, pal.second.first));
        disjoint.insert(make_pair(pal.first, pal.second.second));
      }
      printf("similar\n");
    }
    if (score > cur_score) {
      paths = paths2;
      FixRepForNode2(paths, gr, threshold, disjoin_similar, node, prob_calc);
      return;
    }
  }
  if (disjoin_similar) {
    for (auto it = disjoint.rbegin(); it != disjoint.rend(); it++) {
      printf("disjoin %d %d %d\n", it->first, it->second, paths[it->first].size());
      printf("before: ");
      for (int i = 0; i < paths[it->first].size(); i++) printf("%d ", paths[it->first][i]);
      printf("\n");
      paths.push_back(vector<int>(paths[it->first].begin()+it->second, paths[it->first].end()));
      paths[it->first].erase(paths[it->first].begin()+it->second+1, paths[it->first].end());
      printf("after: ");
      for (int i = 0; i < paths[it->first].size(); i++) printf("%d ", paths[it->first][i]);
      printf(" second: ");
      for (int i = 0; i < paths.back().size(); i++) printf("%d ", paths.back()[i]);
      printf("\n");
      if (paths[it->first].empty()) {
        paths.erase(paths.begin()+it->first);
      }
    }
  }
}

bool FixBigReps(vector<vector<int>>& paths, Graph&gr, int threshold, bool disjoin_similar,
                ProbCalculator& prob_calc) {
  unordered_map<int, int> counts;
  for (int i = 0; i < paths.size(); i++) {
    for (int j = 0; j < paths[i].size(); j++) {
      if (paths[i][j] < 0) continue;
      if (gr.nodes[paths[i][j]]->s.length() > threshold) {
        counts[(paths[i][j]/2)*2]++;
      }
    }
  }
  vector<int> rr;
  for (auto &e: counts) {
    if (e.second >= 2) {
      printf("rep %d %d\n", e.first, e.second);
      rr.push_back(e.first);
    }
  }
  for (auto &e: rr) {
    FixRepForNode2(paths, gr, threshold, disjoin_similar, e, prob_calc); 
  }
  return true;
}

bool FixSomeBigReps(vector<vector<int>>& paths, Graph&gr, int threshold, bool disjoin_similar,
                ProbCalculator& prob_calc) {
  unordered_map<int, int> counts;
  for (int i = 0; i < paths.size(); i++) {
    for (int j = 0; j < paths[i].size(); j++) {
      if (paths[i][j] < 0) continue;
      if (gr.nodes[paths[i][j]]->s.length() > threshold) {
        counts[(paths[i][j]/2)*2]++;
      }
    }
  }
  vector<int> rr;
  for (auto &e: counts) {
    if (e.second >= 2) {
      printf("rep %d %d\n", e.first, e.second);
      rr.push_back(e.first);
    }
  }
  if (rr.empty()) return false;
  int e = rr[rand()%rr.size()];
  FixRepForNode2(paths, gr, threshold, disjoin_similar, e, prob_calc); 
  return true;
}

bool FixRepForNode(int node, vector<vector<int>>& paths, int threshold, Graph& gr,
                   ProbCalculator& prob_calc) {
  vector<vector<int>>paths2(paths);
  vector<vector<int>> paths_with_node;
  vector<int> fp, bfp, afp;
  for (int i = paths2.size()-1; i>=0; i--) {
    bool has_node = false;
    for (int j = 0; j < paths2[i].size(); j++) {
      if ((paths2[i][j]/2)*2 == node) {
        has_node = true;
        break;
      }
    }
    int nl = 0;
    for (int j = 0; j < paths2[i].size(); j++) {
      if (paths2[i][j] < 0) nl += -paths2[i][j];
      else nl += gr.nodes[paths2[i][j]]->s.length();
    }
    if (has_node) {
      fp.push_back(nl);
      paths_with_node.push_back(paths2[i]);
      swap(paths2[paths2.size()-1], paths2[i]);
      paths2.pop_back();
    }
  }

  printf("fix rep node %d %d\n", node, gr.nodes[node]->s.length());
  vector<vector<int>> before;
  vector<vector<int>> after;
  for (int i = 0; i < paths_with_node.size(); i++) {
    int last = -1; bool last_inv = false;
    for (int j = 0; j < paths_with_node[i].size(); j++) {
      if ((paths_with_node[i][j]/2)*2 == node) {
        if (last != -1) {
          printf("self rep\n");
          return false;
        } else {
          if (paths_with_node[i][j] == node) {
            last_inv = false;
            before.push_back(vector<int>(paths_with_node[i].begin()+(last+1),
                                         paths_with_node[i].begin()+j));
            bfp.push_back(fp[i]);
          } else {
            vector<int> pp(paths_with_node[i].begin()+(last+1),
                           paths_with_node[i].begin()+j);
            ReversePath(pp);
            after.push_back(pp);
            afp.push_back(fp[i]);
            last_inv = true;
          }
        }
        last = j;
      }
    }
    assert(last != -1);
    if (!last_inv) {
      after.push_back(vector<int>(paths_with_node[i].begin()+(last+1),
                                  paths_with_node[i].end()));
      afp.push_back(fp[i]);
    } else {
      vector<int> pp(paths_with_node[i].begin()+(last+1),
                     paths_with_node[i].end());
      ReversePath(pp);
      before.push_back(pp);
      bfp.push_back(fp[i]);
    }
  }
  printf("ok bef %d aft %d\n", (int)before.size(), (int)after.size());

  vector<int> opts;
  for (int i = 0; i < after.size(); i++) {
    opts.push_back(i);
  }

  vector<int> best_opts;
  double best_score = -1000000;
  do {
    printf("opt begin\n");
    vector<vector<int>> paths3(paths2);
    for (int i = 0; i < opts.size(); i++) {
      for (int j = (int)before[i].size()-1; j >= 0; j--) {
        if (before[i][j] < 0) continue;
        if (gr.nodes[before[i][j]]->s.length() > threshold) {
          printf(" %d(%d) ", before[i][j], bfp[i]);
          break;
        }
      }
      printf(" %d ", node);
      vector<int> pp(before[i]);
      pp.push_back(node);
      pp.insert(pp.end(), after[opts[i]].begin(), after[opts[i]].end());
      for (int j = 0; j < after[opts[i]].size(); j++) {
        if (after[opts[i]][j] < 0) continue;
        if (gr.nodes[after[opts[i]][j]]->s.length() > threshold) {
          printf("%d(%d)", after[opts[i]][j], afp[opts[i]]);
          break;
        }
      }
      printf("\n");
      if (pp.size() > 1)
        paths3.push_back(pp);
    }
    double score = prob_calc.CalcProb(paths3);
    printf("score %lf\n", score);
    if (score > best_score) {
      best_score = score;
      best_opts = opts;
    }
  } while(next_permutation(opts.begin(), opts.end()));
  vector<vector<int>> paths3(paths2);
  for (int i = 0; i < best_opts.size(); i++) {
    vector<int> pp(before[i]);
    pp.push_back(node);
    pp.insert(pp.end(), after[best_opts[i]].begin(), after[best_opts[i]].end());
    if (pp.size() > 1)
      paths3.push_back(pp);
  }
  paths = paths3;
  return true;
}

