// lrc_for_nine_runners.cpp
// Implements a three-step hierarchical lifting sieve to verify the Lonely Runner Conjecture (LRC) for k+1 = 9 runners.
// The method lifts covers over moduli n = 1 -> 3 -> 9 using bitset representations and subset-GCD sieving.
// Parallelization is used over seed tuples for efficiency.

#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <cstdint>
#include <climits>
#include <functional>
#include <chrono>
#include <thread>
#include <mutex>
#include <atomic>

using namespace std;
using u64 = uint64_t;

#ifndef PRIME
#define PRIME 151 
#endif
#ifndef K
#define K 8       // k = 8 corresponds to k+1 = 9 runners
#endif

constexpr int P_CONST = PRIME;
constexpr int K_CONST = K;

// Fallback GCD for subset-gcd checking (used in sieve condition)
static inline long long gcd_fallback(long long a, long long b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) { long long t = b; b = a % b; a = t; }
    return a;
}

// Compressed bit vector for fast modular coverage representation
struct WordBitset {
    int nbits{}, nwords{};
    vector<u64> w;
    WordBitset() = default;
    explicit WordBitset(int bits){ reset(bits); }
    void reset(int bits){
        nbits = bits;
        nwords = (bits + 63) >> 6;
        w.assign(nwords, 0ULL);
    }
    inline void setBit(int pos){ w[pos >> 6] |= (1ULL << (pos & 63)); }
    inline bool testBit(int pos) const { return ((w[pos >> 6] >> (pos & 63)) & 1ULL) != 0; }
    inline long long count() const {
        long long s = 0;
        for (u64 x : w) s += __builtin_popcountll(x);
        return s;
    }
    inline void orWith(const WordBitset &o){
        int m = min(nwords, o.nwords);
        for (int i = 0; i < m; ++i) w[i] |= o.w[i];
    }
};

// Context struct encapsulates the lifting level (mod Q = np) and bitset precomputations
struct Context {
    int p{}, k{}, n{}, Q{};
    int maxIndex{}, bitlen{}, nwords{};
    vector<WordBitset> vec;  // For each index i, vec[i] is the bitset for coverage testing
};

// Precompute coverage bitsets for a given level n
static Context make_context(int p, int k, int n, bool fullRange){
    Context C{};
    C.p = p; C.k = k; C.n = n; C.Q = n * p;
    C.maxIndex = fullRange ? C.Q - 1 : C.Q / 2;
    C.bitlen = fullRange ? C.Q : C.Q / 2;
    C.nwords = (C.bitlen + 63) >> 6;
    C.vec.resize(C.maxIndex + 1, WordBitset(C.bitlen));

    // For each i in [0, Q) compute modular coverage bitset
    for (int i = 0; i <= C.maxIndex; ++i) {
        WordBitset &B = C.vec[i];
        for (int t = 1; t <= C.bitlen; ++t) {
            int pos = C.bitlen - t;
            int rem = (int)((1LL * t * i) % C.Q);
            // Check if i is not lonely at time t/Q (mod 1); condition per LRC definition
            bool cond = (1LL * rem * (C.k + 1) < C.Q) || (1LL * (C.Q - rem) * (C.k + 1) < C.Q);
            if (cond) B.setBit(pos);
        }
    }
    return C;
}

// Step 1: Brute-force search for all k-tuples (mod p) whose union covers all unsafe times
static vector<vector<int>> find_all_covers(const Context &C){
    const int k = C.k;
    const int p = C.p;
    const int maxI = C.maxIndex;
    const int bitlen = C.bitlen;
    const auto &vec = C.vec;

    vector<int> elems(k, 0);
    vector<WordBitset> covered(k+1, WordBitset(bitlen));
    vector<vector<char>> eliminated(k+1, vector<char>(maxI+1, 0));
    vector<vector<int>> remaining(k+1, vector<int>(bitlen, 0));
    unordered_set<string> seen;
    vector<vector<int>> solutions;
    solutions.reserve(1024);

    // init remaining[0]
    for (int i = 0; i <= maxI; ++i) {
        if (i % p == 0) continue;
        for (int pos = 0; pos < bitlen; ++pos) if (vec[i].testBit(pos)) remaining[0][pos]++;
    }

    function<void(int)> dfs = [&](int used){
        if (used == k) {
            if (covered[used].count() == bitlen) {
                vector<int> tmp(elems.begin(), elems.begin()+used);
                sort(tmp.begin(), tmp.end());
                string key;
                for (size_t i = 0; i < tmp.size(); ++i) { if (i) key.push_back(','); key += to_string(tmp[i]); }
                if (seen.insert(key).second) solutions.push_back(tmp);
            }
            return;
        }

        int nextToCover = -1, best = INT_MAX;
        for (int pos = 0; pos < bitlen; ++pos) {
            if (!covered[used].testBit(pos) && remaining[used][pos] < best) {
                best = remaining[used][pos]; nextToCover = pos;
            }
        }
        int totalToCover = bitlen - (int)covered[used].count();

        // branch-and-bound heuristic
        if (used >= C.k-3 && nextToCover != -1) {
            int largestAvailable = 0;
            int bestCovering_next = 0;
            vector<u64> nextC(C.nwords);
            for (int w = 0; w < C.nwords; ++w) nextC[w] = ~covered[used].w[w];
            int excess = C.nwords*64 - bitlen;
            if (excess > 0) nextC[C.nwords-1] &= (~u64(0)) >> excess;
            nextC[nextToCover>>6] &= ~(1ULL << (nextToCover & 63));
            for (int i = 0; i <= maxI; ++i) {
                if (i % p == 0) continue;
                if (eliminated[used][i]) continue;
                int cnt = 0;
                for (int w = 0; w < C.nwords; ++w) cnt += __builtin_popcountll(vec[i].w[w] & nextC[w]);
                if (cnt > largestAvailable) largestAvailable = cnt;
                if (vec[i].testBit(nextToCover) && cnt+1 > bestCovering_next) bestCovering_next = cnt+1;
            }
            if (totalToCover > (k-used-1) * largestAvailable + bestCovering_next) return;
        }

        eliminated[used+1] = eliminated[used];
        remaining[used+1] = remaining[used];

        for (int i = 0; i <= maxI; ++i) {
            if (i % p == 0) continue;
            if (eliminated[used][i]) continue;
            if (nextToCover == -1 || vec[i].testBit(nextToCover)) {
                elems[used] = i;
                bool extra_ok = true;
                if (used+1 >= k-1) {
                    if (C.n==3 || C.n==5 || C.n==7) {
                        int cnt=0; for(int z=0; z<used+1; ++z) if(elems[z] % C.n == 0) ++cnt;
                        if (cnt >= k-1) extra_ok = false;
                    } else if (C.n==4) {
                        int cnt=0; for(int z=0; z<used+1; ++z) if(elems[z] % 2 == 0) ++cnt;
                        if (cnt >= k-1) extra_ok = false;
                    } else if (C.n==6) {
                        int c2=0,c3=0; for(int z=0; z<used+1; ++z){ if(elems[z]%2==0) ++c2; if(elems[z]%3==0) ++c3; }
                        if (c2>=k-1 || c3>=k-1) extra_ok=false;
                    } else if (C.n==8) {
                        int c2=0; for(int z=0; z<used+1; ++z) if(elems[z]%2==0) ++c2; if(c2>=k-1) extra_ok=false;
                    } else if (C.n==9) {
                        int c3=0; for(int z=0; z<used+1; ++z) if(elems[z]%3==0) ++c3; if(c3>=k-1) extra_ok=false;
                    }
                }
                if (extra_ok) {
                    auto elim_snap = eliminated[used+1];
                    auto rem_snap = remaining[used+1];
                    covered[used+1] = covered[used];
                    covered[used+1].orWith(vec[i]);
                    dfs(used+1);
                    eliminated[used+1] = move(elim_snap);
                    remaining[used+1] = move(rem_snap);
                }
                eliminated[used+1][i] = 1;
                for (int pos = 0; pos < bitlen; ++pos) if (vec[i].testBit(pos)) remaining[used+1][pos]--;
            }
        }
    };

    dfs(0);
    return solutions;
}

// Step 2 & 3: Lifting seeds from prior level to next level (n -> m*n)
// This function performs parallel lifting over the seed list and applies subset-GCD sieve and coverage test
static vector<vector<int>> find_lifted_covers_parallel(const Context &C, const vector<vector<int>> &seeds, int multiplier, int m) {
    const int k = C.k;
    const int maxI = C.maxIndex;
    const int bitlen = C.bitlen;
    const int nwords = C.nwords;
    const auto &vec = C.vec;

    vector<u64> fullmask(nwords, ~u64(0));
    int excess = nwords*64 - bitlen;
    if (excess > 0) fullmask[nwords-1] = (~u64(0)) >> excess;

    size_t N = seeds.size();
    if (N == 0) return {};

    unsigned int hw = thread::hardware_concurrency();
    unsigned int nthreads = hw ? hw : 4u;
    if (nthreads > N) nthreads = (unsigned int)N;

    vector<vector<vector<int>>> thread_results(nthreads);
    vector<unordered_set<string>> thread_seen(nthreads);
    vector<thread> threads;
    mutex merge_mutex;

    auto worker = [&](size_t lo, size_t hi, unsigned tid){
        auto &local_results = thread_results[tid];
        auto &local_seen = thread_seen[tid];
        local_results.reserve(max<size_t>(1, (hi-lo)/4));
        for (size_t si = lo; si < hi; ++si) {
            const auto &s = seeds[si];
            // build candidate lists
            vector<vector<int>> cand(k);
            bool any_empty = false;
            for (int j = 0; j < k; ++j) {
                for (int a = 0; a < m; ++a) {
                    long long val = (long long)s[j] + (long long)a * multiplier;
                    if (val >= 0 && val <= maxI) cand[j].push_back((int)val);
                }
                if (cand[j].empty()) { any_empty = true; break; }
            }
            if (any_empty) continue;

            // order positions by candidate size
            vector<int> order(k);
            iota(order.begin(), order.end(), 0);
            sort(order.begin(), order.end(), [&](int A, int B){ return cand[A].size() < cand[B].size(); });

            vector<int> idx(k, -1);

            function<void(int)> dfs = [&](int depth){
                if (depth == k) {
                    // construct final_idx in natural order
                    vector<int> final_idx(k);
                    for (int t = 0; t < k; ++t) final_idx[order[t]] = idx[order[t]];

                    // subset (k-1)-gcd check with C.n
                    vector<long long> pref(k), suf(k);
                    pref[0] = final_idx[0];
                    for (int i = 1; i < k; ++i) pref[i] = gcd_fallback(pref[i-1], final_idx[i]);
                    suf[k-1] = final_idx[k-1];
                    for (int i = k-2; i >= 0; --i) suf[i] = gcd_fallback(suf[i+1], final_idx[i]);
                    for (int removed = 0; removed < k; ++removed) {
                        long long g;
                        if (removed == 0) g = suf[1];
                        else if (removed == k-1) g = pref[k-2];
                        else g = gcd_fallback(pref[removed-1], suf[removed+1]);
                        g = gcd_fallback(g, C.n);
                        if (g != 1) return;
                    }

                    // wordwise OR
                    vector<u64> acc(nwords);
                    for (int t = 0; t < k; ++t) {
                        const auto &wv = vec[ final_idx[t] ].w;
                        for (int w = 0; w < nwords; ++w) acc[w] |= wv[w];
                    }
                    for (int w = 0; w < nwords; ++w) if (acc[w] != fullmask[w]) return;

                    // success: canonical sorted tuple key
                    vector<int> out = final_idx;
                    sort(out.begin(), out.end());
                    string key;
                    key.reserve(out.size()*6);
                    for (size_t i = 0; i < out.size(); ++i) { if (i) key.push_back(','); key += to_string(out[i]); }
                    if (local_seen.insert(key).second) local_results.push_back(move(out));
                    return;
                }
                int pos = order[depth];
                for (int candidate : cand[pos]) {
                    idx[pos] = candidate;
                    dfs(depth+1);
                }
                idx[pos] = -1;
            }; // end dfs

            dfs(0);
        } // end seeds loop
    }; // end worker

    // partition seeds
    size_t chunk = (N + nthreads - 1) / nthreads;
    size_t start = 0;
    for (unsigned int t = 0; t < nthreads && start < N; ++t) {
        size_t end = min(N, start + chunk);
        threads.emplace_back(worker, start, end, t);
        start = end;
    }
    for (auto &th : threads) th.join();

    // merge thread results into single vector, dedup across threads using global set
    unordered_set<string> global_seen;
    vector<vector<int>> results;
    results.reserve(1024);
    for (unsigned int t = 0; t < nthreads; ++t) {
        for (auto &tuple : thread_results[t]) {
            string key;
            key.reserve(tuple.size()*6);
            for (size_t i = 0; i < tuple.size(); ++i) { if (i) key.push_back(','); key += to_string(tuple[i]); }
            if (global_seen.insert(key).second) results.push_back(tuple);
        }
    }

    return results;
}

// Main driver: constructs and applies the lifting sieve over levels n = 1, 3, 9
int main(){
    const int p = P_CONST;
    const int k = K_CONST;
    using clk = chrono::high_resolution_clock;

    cerr << "Parameters: p = " << p << ", k = " << k << "  (k+1 should be 9 in your setting)\n";

    // Step 1: Compute seed covers S at level n = 1 (half-range mod p)
    auto t1s = clk::now();
    Context C1 = make_context(p, k, 1, false);
    auto S = find_all_covers(C1);
    auto t1e = clk::now();
    double time1 = chrono::duration<double>(t1e - t1s).count();
    cout << "Step1 (n=1): S size = " << S.size() << "  (time " << time1 << " s)\n";
    if (S.empty()) { cout << "S empty; terminating.\n"; return 0; }

    // Step 2: Lift each seed from S using multiplier p, m = 3 (to n = 3)
    auto t2s = clk::now();
    Context C2 = make_context(p, k, 3, true);
    auto T = find_lifted_covers_parallel(C2, S, p, 3);
    auto t2e = clk::now();
    double time2 = chrono::duration<double>(t2e - t2s).count();
    cout << "Step2 (n=3): T size = " << T.size() << "  (time " << time2 << " s)\n";
    if (T.empty()) { cout << "T empty; terminating.\n"; return 0; }

    // Step 3: Lift each seed from T using multiplier 3p, m = 3 (to n = 9)
    auto t3s = clk::now();
    Context C3 = make_context(p, k, 9, true);
    auto U = find_lifted_covers_parallel(C3, T, 3*p, 3);
    auto t3e = clk::now();
    double time3 = chrono::duration<double>(t3e - t3s).count();
    cout << "Step3 (n=9): U size = " << U.size() << "  (time " << time3 << " s)\n";

    // Final result: if U non-empty, LRC verified for this p
    if (!U.empty()) {
        cout << "Example in U: (";
        for (size_t i = 0; i < U.front().size(); ++i) { if (i) cout << ","; cout << U.front()[i]; }
        cout << ")\n";
    } else {
        cout << "U empty; there is no cover.\n";
    }

    cout << "Total times: step1=" << time1 << " s, step2=" << time2 << " s, step3=" << time3 << " s\n";
    return 0;
}