#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <random>
#include <algorithm>
#include <ctime>
#include <bitset>
#include <cmath>

//파일로만 출력 콘솔로 매번 출력 X -> O
//LSH 버킷 수 문제 100개 ? ->
//bloom + (LSH + BANDDP) 분할 적용 코드 -> O
//bloom filter 개선 -> O

using namespace std;

const int d = 256;  // ASCII 문자 수
const int q = 10000019;//큰 N <  mod 값
//const int q = 99991; //작은 N < mod 값
int patternLength = 10; //자르는 길이
const int error = 20; //허용 오차갯수 D

// ==================== Rabin-Karp Hash 기반 ====================
unordered_map<int, vector<int>> buildHashTable(const string& T, int m) {
    unordered_map<int, vector<int>> hashTable;
    int n = T.length();
    int D = 1;

    for (int i = 0; i < m - 1; i++) {
        D = (D * d) % q;
    }

    int hash = 0;
    for (int i = 0; i < m; i++) {
        hash = (d * hash + T[i]) % q;
    }

    hashTable[hash].push_back(0);

    for (int i = 1; i <= n - m; i++) {
        hash = (d * (hash - T[i - 1] * D) + T[i + m - 1]) % q;
        if (hash < 0)
            hash += q;
        hashTable[hash].push_back(i);
    }

    return hashTable;
}

// 염기서열 읽어서 sequence에 저장 반환
string readSequence(const string& fileName) {
    ifstream file(fileName);
    string line, sequence = "";
    while (getline(file, line)) {
        if (line.empty() || line[0] == '>' || line[0] == 'N') continue;
        //원본 파일에는 space, tab같은 이스케이프 문자 있어서 처리
        for (char c : line) {
            char uc = toupper(c);
            if (uc == 'A' || uc == 'T' || uc == 'G' || uc == 'C')
                sequence += uc;
        }
    }
    return sequence;
}

// 해시 테이블 저장
void saveHashTable(const unordered_map<int, vector<int>>& hashTable, const string& fileName) {
    ofstream file(fileName);
    for (auto iter = hashTable.begin(); iter != hashTable.end(); iter++) {
        int hash = iter->first;
        const vector<int>& idx = iter->second;
        file << "Hash: " << hash << " -> Positions: ";
        for (int i = 0; i < idx.size(); i++) {
            file << idx[i] << " ";
        }
        file << "\n";
    }
    file.close();
}

// read 파일 읽기
vector<string> readReads(const string& filename) {
    ifstream file(filename);
    vector<string> reads;
    string line;
    while (getline(file, line)) {
        if (!line.empty()) reads.push_back(line);
    }
    return reads;
}

// mismatch 개수 계산
int countError(const string& subSequence, const string& subRead) {
    int count = 0;
    for (int i = 0; i < subSequence.size(); i++) {
        if (subSequence[i] != subRead[i]) count++;
        if (count > error) return -1;
    }
    return count;
}

// ===================== [추가] Bloom Filter/LSH/BandDP 구조체 =====================

const size_t LSH_BUCKET_LIMIT = 100;
/*
 Bloom Filter 설명 (부족하면 말씀해주세요)
 후보군(매칭 가능성 있는 시드)만 남기고, 나머지는 “비교” 자체를 건너뛸 수 있음
 
 배열크기 크게 늘려보고 리드가 있나 없나만 판단하는거 아닌가 ?
 */

struct FastBloom {
    vector<uint64_t> bits;
    size_t size, hash_num;
    FastBloom(size_t n_bits, int k_hash) : size(n_bits), hash_num(k_hash), bits((n_bits+63)/64, 0) {}
    size_t hash(const string &s, int seed) const {
        uint64_t h = 0xcbf29ce484222325 ^ seed;
        for (char c : s) h = (h^c)*0x100000001b3;
        return h % size;
    }
    void add(const string &s) {
        for (int i=0; i<hash_num; ++i) {
            size_t idx = hash(s, i);
            bits[idx/64] |= (1ULL<<(idx%64));
        }
    }
    bool possibly_contains(const string &s) const {
        for (int i=0; i<hash_num; ++i) {
            size_t idx = hash(s, i);
            if (!(bits[idx/64] & (1ULL<<(idx%64)))) return false;
        }
        return true;
    }
};
//Bloom filter 신버전 -> size 최적화 설정
size_t estimate_bloom_size(size_t n_items, double fpr=0.02) {
    return static_cast<size_t>(-(double)n_items * log(fpr) / pow(log(2),2));
}
int recommend_hash_num(size_t bloom_size, size_t n_items) {
    return max(2, int((double)bloom_size / n_items * log(2) + 0.5));
}


/*
 BandDP 추가설명 -> 편집 거리를 임의로 지정해서 오차 검증 빠르게 하도록 -> 시간단축
 read(길이 L) 전체와 reference의 (후보) substring 전체를 “최대 D개 편집거리(오차)만 허용”해서 매칭 여부를 판정
 D만큼의 대각선 주변(2D+1 폭)의 셀만 dp 수행 -> (O(D*L) 형태로 논문에서 시간복잡도 나왔습니다. )
 오차 D 이하인 후보만 “최종 매칭”으로 채택
 */

int banded_edit_distance(const string& a, const string& b, int band) {
    int m = a.size(), n = b.size();
    const int INF = m + n;
    int D = band;
    vector<vector<int>> dp(m + 1, vector<int>(n + 1, INF));
    dp[0][0] = 0;
    for (int i = 0; i <= m; ++i) {
        int min_val = INF;
        for (int j = max(0, i - D); j <= min(n, i + D); ++j) {
            if (i > 0 && j > 0) dp[i][j] = min(dp[i][j], dp[i - 1][j - 1] + (a[i - 1] != b[j - 1]));
            if (i > 0) dp[i][j] = min(dp[i][j], dp[i - 1][j] + 1);
            if (j > 0) dp[i][j] = min(dp[i][j], dp[i][j - 1] + 1);
            min_val = min(min_val, dp[i][j]);
        }
        if (min_val > D) return INF;
    }
    return dp[m][n];
}

size_t lsh_hash(const string &s, int lsh_bit=10) {
    size_t h = 5381; //임의로 지정 << 벤치마킹 때에서
    for (char c : s) h = ((h << 5) + h) + c;
    return h >> (sizeof(size_t)*8 - lsh_bit);
}

// main
int main() {
    // [1] 실행 시간 측정 시작
    clock_t start = clock();

    // [2] 결과 저장을 위한 테이블/버퍼
    unordered_map<int, vector<int>> resultTable; // read별 매칭 포지션 기록
    vector<string> result_buf;                   // 파일로 쓸 결과 문자열 버퍼

    // [3] 입력(reference/read) 파일명 설정 및 데이터 로딩
    string inputFile = "/Users/kosora/Documents/Xcode/C++/AlgorithmProject/AlgorithmProject/smallinput.txt";
    string readFile = "/Users/kosora/Documents/Xcode/C++/AlgorithmProject/AlgorithmProject/smallread.txt";
    string sequence = readSequence(inputFile);      // 레퍼런스 시퀀스 로딩
    vector<string> reads = readReads(readFile);     // read 리스트 로딩

    // [4] Bloom+LSH+BandDP 초기화
    int k = 15;            // seed(k-mer) 길이
    int lsh_bit = 10;      // LSH 해시 상위 비트 수(버킷 수=1024)
    int D = error;         // 허용 오차(편집거리)

    // --- Bloom filter 최적크기/개선 버전 사용 ---
    size_t bloom_size = estimate_bloom_size(sequence.size(), 0.02);
    int hash_num = recommend_hash_num(bloom_size, sequence.size());
    FastBloom bf(bloom_size, hash_num);

    // --- Reference 전체에서 seed를 Bloom filter에 추가 ---
    for (size_t i = 0; i + k <= sequence.size(); ++i) {
        bf.add(sequence.substr(i, k));
    }

    // --- Reference의 모든 seed를 LSH bucket에 저장해서 후보 압축 ---
    #ifndef LSH_BUCKET_LIMIT
    #define LSH_BUCKET_LIMIT 1000    // 필요시 조정 (버킷 후보군 제한)
    #endif
    unordered_map<size_t, vector<size_t>> lsh_table;
    for (size_t i = 0; i + k <= sequence.size(); ++i) {
        string seed = sequence.substr(i, k);
        size_t bkt = lsh_hash(seed, lsh_bit);
        if (lsh_table[bkt].size() < LSH_BUCKET_LIMIT)
            lsh_table[bkt].push_back(i);
    }

    // [5] 각 read별로 Bloom+LSH+BandDP 매칭 시도
    int readIdx = 0;
    for (const string& read : reads) {
        // 방어: read가 너무 짧으면 스킵 (seed k, 오차 D보다 짧을 때)
        if (read.size() < k || read.size() < D) {
            ++readIdx;
            continue;
        }
        // BandDP를 통과한 후보 위치만 resultTable에 기록
        bool found = false;
        for (int j = 0; j + k <= read.size(); ++j) {
            string seed = read.substr(j, k);
            if (!bf.possibly_contains(seed)) continue;             // Bloom filter 1차 필터링
            size_t bkt = lsh_hash(seed, lsh_bit);
            if (lsh_table.find(bkt) == lsh_table.end()) continue;  // LSH bucket 필터링
            // 버킷에 있는 모든 후보(reference 위치)에서 BandDP 실행
            for (size_t refpos : lsh_table[bkt]) {
                if (refpos + read.size() > sequence.size()) continue;
                string refseg = sequence.substr(refpos, read.size());
                // BandDP를 D/2 band로 먼저 시도, 실패시 D로 확장(속도up)
                int edist = banded_edit_distance(read, refseg, D / 2);
                if (edist > D) {
                    edist = banded_edit_distance(read, refseg, D);
                }
                // 편집거리 D 이하라면 최종 매칭!
                if (edist <= D) {
                    resultTable[readIdx].push_back(refpos);
                    found = true;
                    break; // 첫 매칭만 기록, 여러 위치 원하면 break 제거 가능
                }
            }
            if (found) break;
        }
        readIdx++;
    }

    // [6] 결과를 버퍼에 모은 후 한번에 파일로 저장 (디스크 I/O 최적화)
    for (const auto& [readIdx, positions] : resultTable) {
        const string& read = reads[readIdx];
        string line = "read index: " + to_string(readIdx) + " -> positions: ";
        for (int pos : positions) {
            string refseg = sequence.substr(pos, read.size());
            int edist = banded_edit_distance(read, refseg, error);
            line += "(" + to_string(pos) + ",오차=" + to_string(edist) + ") ";
        }
        line += "\n";
        result_buf.push_back(line);
    }

    // result.txt 파일에 결과 기록
    ofstream file("result.txt");
    for (const string& line : result_buf) {
        file << line;
    }
    file.close();

    // [마지막] 수행시간 출력
    clock_t finish = clock();
    double duration = double(finish - start) / CLOCKS_PER_SEC;
    cout << "전체 프로세스 완료!\n";
    cout << "총 소요 시간: " << duration << "초\n";
    return 0;
}
