//#include <iostream>
//#include <fstream>
//#include <string>
//#include <vector>
//#include <unordered_set>
//#include <unordered_map>
//#include <random>
//#include <algorithm>
//#include <ctime>
//#include <bitset>
//
//// --- 간단 Bloom Filter 구현 --- //
//struct Bloom {
//    std::vector<bool> bitset;
//    int size, hash_num;
//
//    Bloom(size_t n, int k=3) : size(n), hash_num(k), bitset(n, false) {}
//
//    size_t hash(const std::string &s, int seed) const {
//        size_t h = seed;
//        for (char c : s) h = h*101 + c;
//        return h % size;
//    }
//    void add(const std::string &s) {
//        for (int i=0; i<hash_num; ++i) bitset[hash(s, i)] = true;
//    }
//    bool possibly_contains(const std::string &s) const {
//        for (int i=0; i<hash_num; ++i)
//            if (!bitset[hash(s, i)]) return false;
//        return true;
//    }
//};
//
//// --- 간단 Band DP(Levenshtein) 구현 --- //
//int banded_edit_distance(const std::string& a, const std::string& b, int band) {
//    int m = a.size(), n = b.size();
//    const int INF = m+n;
//    int D = band;
//    std::vector<std::vector<int>> dp(m+1, std::vector<int>(n+1, INF));
//    dp[0][0] = 0;
//    for (int i = 0; i <= m; ++i) {
//        for (int j = std::max(0, i-D); j <= std::min(n, i+D); ++j) {
//            if (i>0 && j>0) dp[i][j] = std::min(dp[i][j], dp[i-1][j-1]+(a[i-1]!=b[j-1]));
//            if (i>0) dp[i][j] = std::min(dp[i][j], dp[i-1][j]+1);
//            if (j>0) dp[i][j] = std::min(dp[i][j], dp[i][j-1]+1);
//        }
//    }
//    return dp[m][n];
//}
//
//// --- LSH hashing (상위 b=10bit만 bucket) --- //
//size_t lsh_hash(const std::string &s, int lsh_bit=10) {
//    size_t h = 5381;
//    for (char c : s) h = ((h << 5) + h) + c;
//    return h >> (sizeof(size_t)*8 - lsh_bit); // 상위 비트만 사용
//}
//
//// --- 레퍼런스 fasta에서 ATGC만 추출 --- //
//std::string extract_atgc(const std::string& file) {
//    std::ifstream fin(file);
//    if (!fin) { std::cerr << "파일 오픈 실패!\n"; exit(1);}
//    std::string seq, line;
//    while (getline(fin, line)) {
//        if (line.empty() || line[0]=='>') continue;
//        for (char c : line) {
//            char uc = toupper(c);
//            if (uc=='A'||uc=='T'||uc=='G'||uc=='C') seq += uc;
//        }
//    }
//    return seq;
//}
//
//// 랜덤 쿼리(read) M개 추출
//std::vector<std::string> pick_random_reads(const std::string &ref, int M, int L) {
//    std::vector<std::string> reads;
//    std::mt19937 gen(time(0));
//    std::uniform_int_distribution<> dist(0, ref.size()-L-1);
//    for (int i=0; i<M; ++i)
//        reads.push_back(ref.substr(dist(gen), L));
//    return reads;
//}
//
//// --- 전체 파이프라인 --- //
//int main(int argc, char** argv) {
//    // --- 파라미터 설정 --- //
//    const std::string fasta = "/Users/kosora/Documents/Xcode/C++/AlgorithmProject/AlgorithmProject/Homo_sapiens.GRCh38.dna.alt.fa";
//    int M = 20000;        // 쿼리 개수
//    int L = 100;               // read 길이
//    int k = 15;                // seed 길이
//    int D = 5;                 // 오차 허용
//    int BLOOM_SIZE = 40000000; // 대략 필요용량에 맞춤
//    int LSH_BIT = 14;          // LSH bucket(상위 몇 비트 사용할지)
//
//    // --- 레퍼런스 로딩 --- //
//    std::cout << "Reference (ATGC) 추출 중...\n";
//    std::string N = extract_atgc(fasta);
//    std::cout << "Reference length: " << N.size() << "\n";
//
//    // --- Bloom Filter 생성 및 삽입 --- //
//    Bloom bf(BLOOM_SIZE, 4);
//    std::cout << "Bloom filter 생성 중...\n";
//    for (size_t i=0; i<=N.size()-k; ++i) bf.add(N.substr(i, k));
//
//    // --- LSH 버킷화 --- //
//    std::unordered_map<size_t, std::vector<size_t>> lsh_table; // bucket hash → positions
//    std::cout << "LSH bucket화 중...\n";
//    for (size_t i=0; i<=N.size()-k; ++i) {
//        std::string seed = N.substr(i, k);
//        size_t bkt = lsh_hash(seed, LSH_BIT);
//        lsh_table[bkt].push_back(i);
//    }
//
//    // --- 쿼리(read) 준비 --- //
//    std::vector<std::string> reads = pick_random_reads(N, M, L);
//
//    // --- 메인 탐색 (쿼리별 진행) --- //
//    std::cout << "쿼리 처리 시작...\n";
//    for (int qid=0; qid<M; ++qid) {
//        const std::string &read = reads[qid];
//        bool matched = false;
//        // 1. Read에서 seed(k-mer) 전수 추출
//        for (int j=0; j<=L-k; ++j) {
//            std::string seed = read.substr(j, k);
//            // 2. Bloom filter에 k-mer가 없으면 skip (후보 탈락)
//            if (!bf.possibly_contains(seed)) continue;
//            // 3. LSH bucket에 동일한 hash bucket에 들어있는 레퍼런스 seed 확인
//            size_t bkt = lsh_hash(seed, LSH_BIT);
//            if (lsh_table.find(bkt) == lsh_table.end()) continue;
//            // 4. Bucket 내 후보 position마다 Band DP로 check
//            for (size_t refpos : lsh_table[bkt]) {
//                // 전체 read(L)와 ref substring 비교
//                if (refpos + L > N.size()) continue;
//                std::string refseg = N.substr(refpos, L);
//                int edist = banded_edit_distance(read, refseg, D);
//                if (edist <= D) {
//                    std::cout << "쿼리#" << qid << " (" << read << ")"
//                              << " 매칭포지션: " << refpos << " EDIT_DIST=" << edist << "\n";
//                    matched = true; // 하나라도 찾으면 break
//                    break;
//                }
//            }
//            if (matched) break; // 하나라도 찾으면 read 처리 종료
//        }
//        if (!matched) std::cout << "쿼리#" << qid << " (" << read << ") 매칭 없음\n";
//    }
//
//    return 0;
//}
