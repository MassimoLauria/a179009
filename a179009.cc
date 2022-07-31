#include <vector>
#include <iostream>
#include <iomanip>      // std::setprecision
#include <sstream>
#include <limits>       // std::numeric_limits
#include <cmath>
#include <string>
#include <algorithm>    // std::fill
#include <numeric>      // std::accumulate
#include <thread>
#include <omp.h>

using namespace std;

constexpr int NOREFINEMENTS = 0;
constexpr int NONINCREASING = -1;
constexpr int RANGEISSUE = -2;

constexpr int MODE_COUNT = 0;
constexpr int MODE_LIST = 1;
constexpr int MODE_EXPLORE = 2;

using int64 = long long int;

struct result {
  int64 discovered;
  int64 loop;
  vector<int64> sequences;
  vector<int> maximum;
  vector<int> maxcount;
public:
  void init(int n) {
    if (n>=0) {
      discovered = 0;
      loop = 0;
      sequences.resize(n+1);
      maximum.resize(n+1);
      maxcount.resize(n+1);
      fill(sequences.begin(), sequences.end(),0);
      fill(maximum.begin()  , maximum.end()  ,0);
      fill(maxcount.begin() , maxcount.end() ,0);
    }
  }
};


vector<result> Workers;

int setup_workers(int goal, unsigned int max_threads = 0) {
  // max_thread == 0 means no limit.

  //may return 0 when not able to detect
  auto threads = thread::hardware_concurrency();
  if (threads==0) threads = 1;

  // Cap  the number of threads
  if (max_threads!=0 and max_threads < threads) threads = max_threads;

  /* Init partial statistics (one per thread)

     Each thread counts the numbers of sequences and the relative
     data. The statistic vectors are reversed, in the sense that while
     working toward a number N, the stats for value M are at position
     N-M.
   */
  cerr << "# Cores: "<<threads<<endl;
  omp_set_num_threads(threads);
  Workers.resize(threads);
  for(unsigned int i=0;i<threads;i++)
    Workers[i].init(goal);
  return threads;
}

result collect_from_workers() {

  int n = Workers[0].sequences.size() - 1;
  int threads = Workers.size();
  result solutions;
  solutions.init(n);

  for(int i=0;i<n+1;i++){
    for(int j=0;j<threads;j++) {
      solutions.sequences[i] += Workers[j].sequences[i];
      if (solutions.maximum[i] < Workers[j].maximum[i]) {
        solutions.maximum[i]  = Workers[j].maximum[i];
        solutions.maxcount[i] = Workers[j].maxcount[i];
      } else if (solutions.maximum[i] == Workers[j].maximum[i]) {
        solutions.maxcount[i] += Workers[j].maxcount[i];
      }
    }
  }
  for(int j=0;j<threads;j++) {
    solutions.discovered += Workers[j].discovered;
    solutions.loop       += Workers[j].loop;
  }
  return solutions;
}


class residues {
private:
  int excluded_classes;
  int num_classes;
  vector<int> representants;

  static int infinity;

public:
  residues(residues &x) {
    excluded_classes = x.excluded_classes;
    num_classes = x.num_classes;
    representants = x.representants;
  }

  residues(int m) : excluded_classes{0},
                    num_classes{m},
                    representants(m) {
    if(m<1) {
      cerr<<"Internal error: can't create residues of size "<<m<<endl;
      exit(-1);
    }
    /* Set the representants to infinity (i.e. a very large value in
       the same residue class. Infinity must have the properties:

       1. infinity % m == 0
       2. 2*infinity + 2*m < numeric_limits<int>::max()

    this allows to simplify the update method.
    */
    for(int i=0;i<m;i++) {
      representants[i]=infinity + i;
    }
  }

  int size() {
    return num_classes;
  }

  void update(int value) {
    int rvalue = value % num_classes;
    int update;
    int rupdate;

    if (representants[rvalue] < value) {
      update = representants[rvalue] + value;
      rupdate = update % num_classes;
      representants[rupdate] = min(representants[rupdate], update);
      return;
    }

    for (int i=1;i<num_classes;i++) {
      update = representants[i] + value;
      rupdate = update % num_classes;
      representants[rupdate] = min(representants[rupdate], update);
    }
    representants[rvalue] = value;
    excluded_classes++;

  }

  bool saturated() {
    return excluded_classes >= num_classes;
  }

  int upper_limit() {
    return (infinity - 1) / 2;
  }

  bool isrefinable(int value) {
    int d = representants[value % num_classes];
    return d <= value;
  }

};

struct searchpoint {
  residues *seq;
  int next;
};

vector<searchpoint> SearchSpacePartitions;

/* Set the representants to infinity (i.e. a very large value in
   the same residue class. Infinity must have the properties:

   1. infinity % m == 0
   2. 2*infinity + 2*m < numeric_limits<int>::max()

   this allows to simplify the update method.
*/
int residues::infinity = numeric_limits<int>::max() / 4;

/* Tests whether the sequence is refinable
   It takes for granted that the sequence
   - is increasing
   - is non empty
   - contains only strictly positive numbers

   returns
   - 0 no refinable elements
   - -1 error
   - N>0 a refinable element
 */
int is_refinable(vector<int> &sequence) {
  int current_value;
  typedef vector<int>::size_type index;
  index pointer;
  index length=sequence.size();

  if (length==0) return NOREFINEMENTS;

  // the sequence must be positive and monotone increasing
  if (sequence[0] < 1) return NONINCREASING;
  for (index i=1; i<length; ++i) {
    if (sequence[i-1]>=sequence[i]) {
      return NONINCREASING;
    }
  }

  /* Look for the minimum excluded element */
  current_value = 1;
  pointer = 0;
  while (pointer<length) {
    if (sequence[pointer] > current_value) {
      break;
    }
    ++pointer;
    ++current_value;
  }
  if (pointer>=length)
    return NOREFINEMENTS;    /* No excluded elements */

  residues  watched(current_value);

  if (watched.upper_limit()<= *sequence.crbegin()) {
    return RANGEISSUE;
  }

  ++current_value;
  while(pointer<length) {

    /* Analysis of holes */
    while (current_value < sequence[pointer]) {
      watched.update(current_value);

      if (watched.saturated()) {
        return sequence[pointer];
      }
      ++current_value;
    }

    /* Analysis of sequence values */
    if (watched.isrefinable(current_value)) {
      return current_value;
    }

    ++current_value;
    ++pointer;
  }

  return NOREFINEMENTS;
}


/* Reads integers from a string and append them
   to the vector */
void append_data(vector<int> &V, string data) {
  int    value_read;
  string text_read;
  stringstream buffer;
  buffer << data;

  while (!buffer.eof()) {
    text_read.clear();
    buffer >> text_read;  /* read as text first */
    if (text_read.length()==0)  /* no more numbers */
      return;
    try {
      value_read = stoi(text_read);
      V.push_back(value_read);
    }
    catch (const std::invalid_argument& ia) {
	  cerr << "Invalid input: " << ia.what() << '\n';
      exit(-1);
    }
  }
    return;
}

void emit_sequence(residues &watched, int upper_bound) {
  int i=1;
  int m=watched.size();
  while(i<m) {cout<< i++ <<" ";}
  i=m+1;
  while (i<upper_bound) {
    if (!watched.isrefinable(i)) { cout<<i<<" ";}
    i++;
  }
}

void save_searchpoint(residues &watched, int upper_bound) {
  auto *p = new residues(watched);
  SearchSpacePartitions.push_back({p,upper_bound});
}

int compute_sum(residues &watched, int upper_bound) {
  int i=1;
  int m=watched.size();
  int total = m*(m-1) / 2;
  i=m+1;
  while (i<upper_bound) {
    if (!watched.isrefinable(i)) total += i;
    i++;
  }
  return total;
}


void partitions_recursive(int n,
                          int sum,
                          int current,
                          residues watched,
                          int mode = MODE_COUNT) {
  /* Guarantees:
     - n - sum >= 0

     The state of watched/current is the same we would obtain checking
     the refinability of the sequence in head, going forward until
     current-1, as in the function is_refinable.
  */

  int th = omp_get_thread_num();

  while (current <= (n - sum)) {

    Workers[th].loop++;

    if (!watched.isrefinable(current)) {

      // Pick it if unrefinable
      Workers[th].discovered++;
      Workers[th].sequences[sum+current]++;
      int maximum   = Workers[th].maximum[sum+current];
      int maxcount  = Workers[th].maxcount[sum+current];
      if (maximum < current) {
        maximum  = current;
        maxcount = 1;
      } else if (maximum == current) {
        maxcount += 1;
      }
      Workers[th].maximum[sum+current] = maximum;
      Workers[th].maxcount[sum+current]= maxcount;

      partitions_recursive(n,
                           sum+current,
                           current+1,
                           watched,
                           mode);
      }

    // Don't pick it either ways
    watched.update(current);
    if (watched.saturated()) return;
    current++;
  }

  // end of the recursion
  if (mode == MODE_LIST and sum==n) {
    emit_sequence(watched, current);
    cout<<endl;
  } else if (mode == MODE_EXPLORE) {
    save_searchpoint(watched, current);
  }
}


result partitions(int n,
                  int mode = MODE_COUNT) {
  if (n==0) {
    result S;
    S.init(0);
    S.discovered = 1;
    S.loop = 0;
    S.sequences[0] = 1;
    S.maximum[0]   = 0;
    S.maxcount[0]  = 1;
    return S;
  }

  setup_workers(n, mode == MODE_COUNT? 0 : 1);

  /* The shortest sequence 1 2 3... H, with sum > n */
  int H=0;
  while (H*(H+1) <= 2*n) ++H;

  /* Generates all prefixes with the minimum excluded in decreasing
     position.

     No excluded : if n = 1 + 2 + ... + K, then in the first cycle we
       have that sum==n and that the minimum excluded is at value
       (K+1). This essentially represents the sequence with no
       excluded elements. This gives a partition only if n is
       a triangular number.

     Min excluded is 1: in this case no other element can be excluded
       in a non refinable partition. Then n = 2 + 3 + 4 + ... + K and
       therefore n+1 is a triangular number.
   */

  for(int i=H; i>0; --i) {
    // sequence 1 2 3 4 ... i-1
    int sum = (i*(i-1))/2;
    // int ThrId = omp_get_thread_num();
    Workers[0].discovered += 1;
    Workers[0].sequences[sum] = 1;
    Workers[0].maximum[sum]   = i-1;
    Workers[0].maxcount[sum]  = 1;
    // observe min excluded at position i
    residues watched(i);
    partitions_recursive(min(150+sum,n), sum, i+1,
                         watched,
                         MODE_EXPLORE);
  }

  if (mode == MODE_EXPLORE) {
    return collect_from_workers();
  }

  /* Compose together the solutions found in each thread. */
  int Parts=SearchSpacePartitions.size();
  cerr<<"# Search space split in "<<Parts<<" parts"<<endl;
  #pragma omp parallel for schedule(dynamic)
  for(int i=0;i<Parts;i++) {
    residues *watched = SearchSpacePartitions[i].seq;
    int next = SearchSpacePartitions[i].next;
    int sum = compute_sum(*watched, next);
    //emit_sequence(*watched, next);
    partitions_recursive(n, sum, next, *watched, mode);
  }

  for(int i=0;i<Parts;i++) {
    delete SearchSpacePartitions[i].seq;
  }

  return collect_from_workers();
}


void batch() {
  vector<int> sequence;
  int res;
  string temp;

  while(!cin.eof()) {

    /* Get another input line */
    sequence.clear();
    getline(cin,temp);
    append_data(sequence,temp);
    if (sequence.size() < 1) continue;

    /* Test the sequence */
    res = is_refinable(sequence);
    for(const auto v : sequence) cout<<v<<" ";
    if (res > 0) cout<<"# refinable at "<<res;
    else if (res == NOREFINEMENTS) cout<<"# not refinable";
    else if (res == NONINCREASING) cout<<"# invalid sequence";
    else if (res == RANGEISSUE) cout<<"# too large values";
    else cout<<"# error";
    cout<<endl;
  }
}

void usage_exit(char const* prgname) {
    cerr << "Usage: " << prgname << " [N | -e N | N M]"<< endl;
    cerr << endl;
    cerr << " - no arguments: batch mode, check sequences in stdin"<< endl;
    cerr << " - one arguments: computes a(N)"<< endl;
    cerr << " - two arguments (N<=M): computes a(N)...a(M)"<< endl;
    cerr << " - '-e N' enumerates non refinable partitions of N"<< endl;
    cerr << " - '-s N' partition search space around partitions of N"<< endl;
    cerr << endl;
    cerr << "See: https://oeis.org/A179009"<< endl;
    cerr << "     https://oeis.org/A179009/b179009.txt"<< endl;
    exit(-1);
}

string human_duration(double duration) {
  stringstream format;

  int64 totalms = round(duration);
  int64 ms = totalms % 1000;
  int64 ss = (totalms / 1000) % 60;
  int64 mm = (totalms / (1000*60)) % 60;
  int64 hh = (totalms / (1000*60*60)) % 24;
  int64 gg = totalms / (1000*60*60*24);

  if (gg>0) format<<gg<<"g";
  if (gg>0 or hh>0) format<<hh<<"h ";
  if (gg>0 or hh>0 or mm>0) format<<mm<<"m ";
  if (gg>0 or hh>0 or mm>0 or ss>0) format<<ss<<"s ";
  format<<ms<<"ms";
  return format.str();
}

int main(int argc,char **argv) {
  int    mode = MODE_COUNT;
  int    pos_argc = argc-1;
  char** pos_argv = argv+1;
  string enumerateflag{"-e"};
  string splitflag{"-s"};

  // Process command line args
  if (argc > 1) {
    if (enumerateflag == argv[1]) {
      mode = MODE_LIST;
      pos_argc -=1;
      pos_argv +=1;
    }

    if (splitflag == argv[1]) {
      mode = MODE_EXPLORE;
      pos_argc -=1;
      pos_argv +=1;
    }
  }

  // Batch verification
  if (pos_argc == 0) {
    cerr << "# Batch verification of sequences."<<endl;
    batch();
    exit(0);
  }

  // split and enumeration accept only one argument.
  if ((mode == MODE_LIST or mode == MODE_EXPLORE )
      and pos_argc>1) {
    usage_exit(argv[0]);
  }

  // Parse limits for counting
  int low, high;
  try {
    low = stoi(string(pos_argv[0]));
    high = low;
    if (pos_argc >= 2) {
      high = stoi(string(pos_argv[1]));
    }
  }
  catch (const std::invalid_argument& ia) {
    usage_exit(argv[0]);
  }

  // A bit of info
  cerr << "# Email: "<<"Massimo Lauria <massimo.lauria@uniroma1.it>"<<endl;
  cerr << "# N    : "<<low;
  if (high != low) cerr << "..."<<high;
  cerr << endl;
  cerr << "# Mode : ";
  switch(mode) {
  case MODE_LIST:
    cerr<<"enumeration"; break;
  case MODE_EXPLORE:
    cerr<<"split"; break;
  case MODE_COUNT:
    cerr<<"counting"; break;
  default:
    cerr<<"counting";
  }
  cerr <<endl;

  // enumeration mode: list all partitions and exit
  if (mode == MODE_LIST) {
    partitions(high, mode);
    return 0;
  }

  if (mode == MODE_EXPLORE) {
    partitions(high, mode);
    int Parts=SearchSpacePartitions.size();
    for(int i=0;i<Parts;i++) {
      residues *watched = SearchSpacePartitions[i].seq;
      int next = SearchSpacePartitions[i].next;
      emit_sequence(*watched, next);
      cout<<"* "<<next<<endl;
    }
    return 0;
  }


  // counting mode: count partitions and other stats
  auto start = chrono::steady_clock::now();
  result S = partitions(high);
  auto end = chrono::steady_clock::now();
  double duration = chrono::duration <double, milli> (end-start).count();

  cerr << "# Time : "<< duration << " ms" << endl;
  cerr << "# TimeH: "<< human_duration(duration)<<endl;
  cerr<<"# Discovered "<<S.discovered<<endl;
  cerr<<"# Loop       "<<S.loop<<endl;
  int threads = Workers.size();
  for (int i=0; i < threads; i++) {
    cerr<<"# Discovered["<<i<<"] "<<Workers[i].discovered<<endl;
    cerr<<"# Loop      ["<<i<<"] "<<Workers[i].loop<<endl;
  }
  for (int n=low; n <= high; n++) {
    cout<<n
        <<" "<<S.sequences[n]
        <<" "<<S.maximum[n]
        <<" "<<S.maxcount[n]
        <<endl;
    }
}
