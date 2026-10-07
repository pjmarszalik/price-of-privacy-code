#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
using namespace std;

/*
Exact local-view / forced-edge normal-form checker for three-state graph-CMSS.

The program uses only finite partitions and exact integer arithmetic.  It does
not assume linear shares, binary share alphabets, or an XOR reconstruction rule.
For the K4 proof, only the uniform stratum U=(1,1,1) and the exceptional stratum
Q=(2,1,1) are logically needed after the source-weight collapse theorem.

Modes:
  ./normal_form_search --proof        required O2--O5 U/Q UNSAT checks
  ./normal_form_search --k4           all six K4 orbit representatives in U/Q
  ./normal_form_search --calibration  four published K3 calibration domains, all strata
  ./normal_form_search --paper        --k4 followed by --calibration
  ./normal_form_search NAME STRATUM    one instance, e.g. O2 U or F11 Q
*/

using Domain = vector<vector<uint8_t>>;

struct VConfig {
    int k = 0;
    vector<array<uint8_t, 3>> lab;  // input index -> source state -> local label
    uint16_t full = 0;              // all unordered pairs of distinct labels
};

struct EdgeMatrix {
    int na = 0, nb = 0;
    vector<uint16_t> sepA, sepB;
    uint16_t A(int a, int b) const { return sepA[(size_t)a * nb + b]; }
    uint16_t B(int a, int b) const { return sepB[(size_t)a * nb + b]; }
};

struct DSU {
    vector<int> p, sz;
    explicit DSU(int n) : p(n), sz(n, 1) { iota(p.begin(), p.end(), 0); }
    int find(int x) {
        while (p[x] != x) {
            p[x] = p[p[x]];
            x = p[x];
        }
        return x;
    }
    void unite(int a, int b) {
        a = find(a);
        b = find(b);
        if (a == b) return;
        if (sz[a] < sz[b]) swap(a, b);
        p[b] = a;
        sz[a] += sz[b];
    }
};

static vector<array<uint8_t, 3>> canonical_maps() {
    // The five set partitions of a three-point source, with canonical labels.
    return {{{0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1}, {0, 1, 2}}};
}

static vector<array<uint8_t, 3>> allowed_maps(
    const array<uint8_t, 3>& first, const array<int, 3>& weights) {
    const int k = 1 + max({(int)first[0], (int)first[1], (int)first[2]});
    vector<int> target_mass(k, 0);
    for (int r = 0; r < 3; ++r) target_mass[first[r]] += weights[r];

    vector<array<uint8_t, 3>> out;
    for (int a = 0; a < k; ++a) {
        for (int b = 0; b < k; ++b) {
            for (int c = 0; c < k; ++c) {
                array<uint8_t, 3> f = {(uint8_t)a, (uint8_t)b, (uint8_t)c};
                vector<int> mass(k, 0), count(k, 0);
                for (int r = 0; r < 3; ++r) {
                    mass[f[r]] += weights[r];
                    count[f[r]]++;
                }
                bool ok = true;
                for (int j = 0; j < k; ++j) {
                    if (count[j] == 0 || mass[j] != target_mass[j]) ok = false;
                }
                if (ok) out.push_back(f);
            }
        }
    }
    return out;
}

struct ClassConfig {
    vector<array<uint8_t, 3>> maps;
    int k = 0;
};

static void rec_class(int pos,
                      const vector<int>& indices,
                      const vector<vector<array<uint8_t, 3>>>& choices,
                      ClassConfig& cur,
                      vector<ClassConfig>& out) {
    if (pos == (int)indices.size()) {
        out.push_back(cur);
        return;
    }
    const int idx = indices[pos];
    for (const auto& f : choices[pos]) {
        cur.maps[idx] = f;
        rec_class(pos + 1, indices, choices, cur, out);
    }
}

static vector<ClassConfig> gen_class(const vector<int>& indices,
                                     int domain_size,
                                     const array<int, 3>& weights) {
    if (indices.empty()) {
        ClassConfig c;
        c.maps.assign(domain_size, {255, 255, 255});
        return {c};
    }

    vector<ClassConfig> out;
    for (const auto& first : canonical_maps()) {
        const int k = 1 + max({(int)first[0], (int)first[1], (int)first[2]});
        ClassConfig cur;
        cur.maps.assign(domain_size, {255, 255, 255});
        cur.maps[indices[0]] = first;
        cur.k = k;

        vector<int> rest;
        vector<vector<array<uint8_t, 3>>> choices;
        const auto allowed = allowed_maps(first, weights);
        for (size_t j = 1; j < indices.size(); ++j) {
            rest.push_back(indices[j]);
            choices.push_back(allowed);
        }
        if (rest.empty()) {
            out.push_back(cur);
        } else {
            rec_class(0, rest, choices, cur, out);
        }
    }
    return out;
}

static vector<VConfig> gen_vertex(const Domain& D,
                                  int vertex,
                                  const array<int, 3>& weights) {
    const int d = (int)D.size();
    vector<int> I0, I1;
    for (int i = 0; i < d; ++i) (D[i][vertex] ? I1 : I0).push_back(i);

    const auto C0 = gen_class(I0, d, weights);
    const auto C1 = gen_class(I1, d, weights);
    vector<VConfig> out;
    out.reserve(C0.size() * C1.size());

    for (const auto& a : C0) {
        for (const auto& b : C1) {
            VConfig c;
            c.k = a.k + b.k;
            c.lab.resize(d);
            for (int i = 0; i < d; ++i) {
                if (a.maps[i][0] != 255) {
                    c.lab[i] = a.maps[i];
                } else {
                    for (int r = 0; r < 3; ++r) c.lab[i][r] = b.maps[i][r] + a.k;
                }
            }
            const int pairs = c.k * (c.k - 1) / 2;
            if (pairs > 16) throw runtime_error("local-label pair mask overflow");
            c.full = pairs == 16 ? 0xFFFFu
                                 : (pairs == 0 ? 0u : (uint16_t)((1u << pairs) - 1u));
            out.push_back(std::move(c));
        }
    }
    return out;
}

static vector<vector<int>> pair_index(int k) {
    vector<vector<int>> ix(k, vector<int>(k, -1));
    int p = 0;
    for (int a = 0; a < k; ++a) {
        for (int b = a + 1; b < k; ++b) ix[a][b] = ix[b][a] = p++;
    }
    return ix;
}

static EdgeMatrix build_edge(const vector<VConfig>& A,
                             const vector<VConfig>& B,
                             int domain_size) {
    EdgeMatrix M;
    M.na = (int)A.size();
    M.nb = (int)B.size();
    M.sepA.resize((size_t)M.na * M.nb);
    M.sepB.resize((size_t)M.na * M.nb);

    vector<vector<vector<int>>> groupsA(M.na), groupsB(M.nb);
    for (int ia = 0; ia < M.na; ++ia) {
        groupsA[ia].assign(A[ia].k, {});
        for (int x = 0; x < domain_size; ++x) {
            for (int r = 0; r < 3; ++r) groupsA[ia][A[ia].lab[x][r]].push_back(3 * x + r);
        }
    }
    for (int ib = 0; ib < M.nb; ++ib) {
        groupsB[ib].assign(B[ib].k, {});
        for (int x = 0; x < domain_size; ++x) {
            for (int r = 0; r < 3; ++r) groupsB[ib][B[ib].lab[x][r]].push_back(3 * x + r);
        }
    }

    vector<vector<vector<int>>> pixA(M.na), pixB(M.nb);
    for (int ia = 0; ia < M.na; ++ia) pixA[ia] = pair_index(A[ia].k);
    for (int ib = 0; ib < M.nb; ++ib) pixB[ib] = pair_index(B[ib].k);

    for (int ia = 0; ia < M.na; ++ia) {
        for (int ib = 0; ib < M.nb; ++ib) {
            DSU uf(3 * domain_size);
            for (const auto& g : groupsA[ia]) {
                for (size_t j = 1; j < g.size(); ++j) uf.unite(g[0], g[j]);
            }
            for (const auto& g : groupsB[ib]) {
                for (size_t j = 1; j < g.size(); ++j) uf.unite(g[0], g[j]);
            }

            vector<int> compA(A[ia].k), compB(B[ib].k);
            for (int l = 0; l < A[ia].k; ++l) compA[l] = uf.find(groupsA[ia][l][0]);
            for (int l = 0; l < B[ib].k; ++l) compB[l] = uf.find(groupsB[ib][l][0]);

            uint16_t sepA = 0, sepB = 0;
            for (int x = 0; x < A[ia].k; ++x) {
                for (int y = x + 1; y < A[ia].k; ++y) {
                    if (compA[x] != compA[y]) sepA |= (uint16_t)(1u << pixA[ia][x][y]);
                }
            }
            for (int x = 0; x < B[ib].k; ++x) {
                for (int y = x + 1; y < B[ib].k; ++y) {
                    if (compB[x] != compB[y]) sepB |= (uint16_t)(1u << pixB[ib][x][y]);
                }
            }
            M.sepA[(size_t)ia * M.nb + ib] = sepA;
            M.sepB[(size_t)ia * M.nb + ib] = sepB;
        }
    }
    return M;
}

struct CoverTable {
    int na = 0, nb = 0, words = 0;
    vector<size_t> offset;
    vector<int> pair_count;
    vector<uint64_t> data;
    const uint64_t* ptr(int a, uint16_t missing) const {
        return &data[offset[a] + (size_t)missing * words];
    }
};

static CoverTable build_cover(const vector<VConfig>& A,
                              int nb,
                              const vector<uint16_t>& sepA) {
    CoverTable T;
    T.na = (int)A.size();
    T.nb = nb;
    T.words = (nb + 63) / 64;
    T.offset.resize(T.na);
    T.pair_count.resize(T.na);

    size_t total = 0;
    for (int a = 0; a < T.na; ++a) {
        const int p = A[a].k * (A[a].k - 1) / 2;
        T.pair_count[a] = p;
        T.offset[a] = total;
        total += ((size_t)1 << p) * T.words;
    }
    T.data.assign(total, 0);

    for (int a = 0; a < T.na; ++a) {
        const int p = T.pair_count[a];
        const size_t masks = (size_t)1 << p;
        const size_t base = T.offset[a];
        for (int b = 0; b < nb; ++b) {
            const uint16_t s = sepA[(size_t)a * nb + b];
            T.data[base + (size_t)s * T.words + b / 64] |= 1ULL << (b % 64);
        }
        // Superset zeta transform: entry[missing] is the union of exact supersets.
        for (int bit = 0; bit < p; ++bit) {
            for (size_t m = 0; m < masks; ++m) {
                if ((m >> bit) & 1u) continue;
                const size_t dst = base + m * T.words;
                const size_t src = base + (m | (1ULL << bit)) * T.words;
                for (int w = 0; w < T.words; ++w) T.data[dst + w] |= T.data[src + w];
            }
        }
    }
    return T;
}

struct Result {
    bool sat = false;
    vector<int> counts;
    vector<int> witness;
    unsigned long long triples = 0;
};

static Result solve4(const Domain& D, const array<int, 3>& weights) {
    Result out;
    vector<vector<VConfig>> C(4);
    for (int v = 0; v < 4; ++v) {
        C[v] = gen_vertex(D, v, weights);
        out.counts.push_back((int)C[v].size());
    }

    const vector<pair<int, int>> endpoints = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
    array<EdgeMatrix, 6> E;
    for (int k = 0; k < 6; ++k) {
        E[k] = build_edge(C[endpoints[k].first], C[endpoints[k].second], (int)D.size());
    }
    const auto& e01 = E[0];
    const auto& e02 = E[1];
    const auto& e03 = E[2];
    const auto& e12 = E[3];
    const auto& e13 = E[4];
    const auto& e23 = E[5];

    const auto cv03 = build_cover(C[0], (int)C[3].size(), e03.sepA);
    const auto cv13 = build_cover(C[1], (int)C[3].size(), e13.sepA);
    const auto cv23 = build_cover(C[2], (int)C[3].size(), e23.sepA);
    const int words = cv03.words;
    vector<uint64_t> candidates(words);

    for (int a = 0; a < (int)C[0].size(); ++a) {
        for (int b = 0; b < (int)C[1].size(); ++b) {
            const uint16_t s0ab = e01.A(a, b);
            const uint16_t s1ab = e01.B(a, b);
            for (int c = 0; c < (int)C[2].size(); ++c) {
                ++out.triples;
                const uint16_t miss0 = C[0][a].full & ~(uint16_t)(s0ab | e02.A(a, c));
                const uint16_t miss1 = C[1][b].full & ~(uint16_t)(s1ab | e12.A(b, c));
                const uint16_t miss2 = C[2][c].full & ~(uint16_t)(e02.B(a, c) | e12.B(b, c));
                const uint64_t* p0 = cv03.ptr(a, miss0);
                const uint64_t* p1 = cv13.ptr(b, miss1);
                const uint64_t* p2 = cv23.ptr(c, miss2);

                bool any = false;
                for (int w = 0; w < words; ++w) {
                    candidates[w] = p0[w] & p1[w] & p2[w];
                    any = any || candidates[w] != 0;
                }
                if (!any) continue;

                for (int w = 0; w < words; ++w) {
                    uint64_t z = candidates[w];
                    while (z) {
                        const int bit = __builtin_ctzll(z);
                        const int d = w * 64 + bit;
                        z &= z - 1;
                        if (d >= (int)C[3].size()) continue;
                        const uint16_t s3 = (uint16_t)(e03.B(a, d) | e13.B(b, d) | e23.B(c, d));
                        if (s3 == C[3][d].full) {
                            out.sat = true;
                            out.witness = {a, b, c, d};
                            return out;
                        }
                    }
                }
            }
        }
    }
    return out;
}

static Result solve3(const Domain& D, const array<int, 3>& weights) {
    Result out;
    vector<vector<VConfig>> C(3);
    for (int v = 0; v < 3; ++v) {
        C[v] = gen_vertex(D, v, weights);
        out.counts.push_back((int)C[v].size());
    }

    const auto e01 = build_edge(C[0], C[1], (int)D.size());
    const auto e02 = build_edge(C[0], C[2], (int)D.size());
    const auto e12 = build_edge(C[1], C[2], (int)D.size());

    for (int a = 0; a < (int)C[0].size(); ++a) {
        for (int b = 0; b < (int)C[1].size(); ++b) {
            const uint16_t s0ab = e01.A(a, b);
            const uint16_t s1ab = e01.B(a, b);
            for (int c = 0; c < (int)C[2].size(); ++c) {
                ++out.triples;
                if ((uint16_t)(s0ab | e02.A(a, c)) != C[0][a].full) continue;
                if ((uint16_t)(s1ab | e12.A(b, c)) != C[1][b].full) continue;
                if ((uint16_t)(e02.B(a, c) | e12.B(b, c)) != C[2][c].full) continue;
                out.sat = true;
                out.witness = {a, b, c};
                return out;
            }
        }
    }
    return out;
}

static Domain parse_domain(const vector<string>& points) {
    if (points.empty()) throw runtime_error("empty domain");
    const size_t n = points.front().size();
    Domain D;
    for (const string& s : points) {
        if (s.size() != n) throw runtime_error("inconsistent point length");
        vector<uint8_t> x(n);
        for (size_t i = 0; i < n; ++i) {
            if (s[i] != '0' && s[i] != '1') throw runtime_error("non-binary domain point");
            x[i] = (uint8_t)(s[i] - '0');
        }
        D.push_back(std::move(x));
    }
    return D;
}

static const map<string, vector<string>>& instances() {
    static const map<string, vector<string>> data = {
        // Representatives used in the paper.  In particular O6 is T6.
        {"O1", {"0000", "0001", "0010", "0100"}},
        {"O2", {"0000", "0001", "0011", "0100"}},
        {"O3", {"0000", "0010", "0011", "0101"}},
        {"O4", {"0000", "0011", "0101", "1000"}},
        {"O5", {"0000", "0101", "1000", "1011"}},
        {"O6", {"0000", "0111", "1001", "1010"}},
        // Published K3 calibration instances.
        {"PARITY", {"000", "011", "101", "110"}},
        {"F11", {"000", "001", "010", "100"}},
        {"F13", {"000", "001", "010", "111"}},
        {"F14", {"000", "010", "100", "101"}},
    };
    return data;
}

static const map<string, array<int, 3>>& strata() {
    static const map<string, array<int, 3>> data = {
        {"G", {1, 2, 4}}, {"E", {1, 1, 3}}, {"U", {1, 1, 1}},
        {"H", {3, 2, 1}}, {"Q", {2, 1, 1}},
    };
    return data;
}

static Result run_case(const string& name, const string& stratum) {
    const auto itD = instances().find(name);
    const auto itW = strata().find(stratum);
    if (itD == instances().end()) throw runtime_error("unknown instance: " + name);
    if (itW == strata().end()) throw runtime_error("unknown stratum: " + stratum);
    const Domain D = parse_domain(itD->second);
    if (D.front().size() == 4) return solve4(D, itW->second);
    if (D.front().size() == 3) return solve3(D, itW->second);
    throw runtime_error("only K3 and K4 instances are supported");
}

static string join(const vector<int>& xs) {
    string out;
    for (size_t i = 0; i < xs.size(); ++i) {
        if (i) out += ',';
        out += to_string(xs[i]);
    }
    return out;
}

static void print_result(const string& name, const string& stratum, const Result& r) {
    cout << name << ' ' << stratum << ' ' << (r.sat ? "SAT" : "UNSAT")
         << " counts=" << join(r.counts) << " triples=" << r.triples;
    if (r.sat) cout << " witness=" << join(r.witness);
    cout << '\n';
}

static bool expected_k4(const string& name, const string& stratum) {
    return name == "O6" && stratum == "U";
}

static bool expected_k3(const string& name, const string& stratum) {
    if (name != "PARITY") return false;
    return stratum == "H" || stratum == "Q";
}

static void checked_case(const string& name,
                         const string& stratum,
                         bool expected_sat) {
    const Result r = run_case(name, stratum);
    print_result(name, stratum, r);
    if (r.sat != expected_sat) {
        throw runtime_error("unexpected result for " + name + " " + stratum);
    }
}

static void run_k4(bool proof_only) {
    cout << "# K4 exact normal-form evaluation\n";
    const vector<string> names = proof_only
        ? vector<string>{"O2", "O3", "O4", "O5"}
        : vector<string>{"O1", "O2", "O3", "O4", "O5", "O6"};
    for (const string& name : names) {
        for (const string& stratum : {string("U"), string("Q")}) {
            checked_case(name, stratum, expected_k4(name, stratum));
        }
    }
    cout << "K4 CHECKS PASS\n";
}

static void run_calibration() {
    cout << "# Published K3 external calibration\n";
    for (const string& name : {string("PARITY"), string("F11"), string("F13"), string("F14")}) {
        for (const string& stratum : {string("G"), string("E"), string("U"), string("H"), string("Q")}) {
            checked_case(name, stratum, expected_k3(name, stratum));
        }
    }
    cout << "K3 CALIBRATION PASS\n";
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    try {
        if (argc == 1 || (argc == 2 && string(argv[1]) == "--paper")) {
            run_k4(false);
            run_calibration();
            cout << "ALL CHECKS PASS\n";
            return 0;
        }
        if (argc == 2 && string(argv[1]) == "--proof") {
            run_k4(true);
            return 0;
        }
        if (argc == 2 && string(argv[1]) == "--k4") {
            run_k4(false);
            return 0;
        }
        if (argc == 2 && string(argv[1]) == "--calibration") {
            run_calibration();
            return 0;
        }
        if (argc == 3) {
            const string name = argv[1], stratum = argv[2];
            const Result r = run_case(name, stratum);
            print_result(name, stratum, r);
            return 0;
        }
        cerr << "usage: " << argv[0]
             << " [--paper|--proof|--k4|--calibration|INSTANCE STRATUM]\n";
        return 2;
    } catch (const exception& e) {
        cerr << "error: " << e.what() << '\n';
        return 1;
    }
}
