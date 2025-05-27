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

using namespace std;

const int d = 256;  // ASCII 문자 수
const int q = 10000019;//큰 N
//const int q = 99991;//작은 N
int patternLength = 10;//자르는 길이
const int error = 5;//허용 오차갯수 D

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
/*
 Bloom Filter 설명 (부족하면 말씀해주세요)
 후보군(매칭 가능성 있는 시드)만 남기고, 나머지는 “비교” 자체를 건너뛸 수 있음
 */

struct Bloom {
    vector<bool> bitset;
    int size, hash_num;

    Bloom(size_t n, int k=3) : size(n), hash_num(k), bitset(n, false) {}

    size_t hash(const string &s, int seed) const {
        size_t h = seed;
        for (char c : s) h = h*101 + c;
        return h % size;
    }
    //  미리 BloomFilter에 미리 저장해놓고. 000.add(seed)
    void add(const string &s) {
        for (int i=0; i<hash_num; ++i) bitset[hash(s, i)] = true;
    }
    //read의 seed를 검사한다.
    bool possibly_contains(const string &s) const {
        for (int i=0; i<hash_num; ++i)
            //없을 경우 "!bitset[hash(s, i)]" false 후 skip
            if (!bitset[hash(s, i)]) return false;
        //있다면 후보군으로 올리기
        return true;
    }
};
/*
 BandDP 추가설명 -> 편집 거리를 임의로 지정해서 오차 검증 빠르게 하도록 -> 시간단축
 read(길이 L) 전체와 reference의 (후보) substring 전체를 “최대 D개 편집거리(오차)만 허용”해서 매칭 여부를 판정
 D만큼의 대각선 주변(2D+1 폭)의 셀만 dp 수행 -> (O(D*L) 형태로 논문에서 시간복잡도 나왔습니다. )
 오차 D 이하인 후보만 “최종 매칭”으로 채택
 */

int banded_edit_distance(const string& a, const string& b, int band) {
    int m = a.size(), n = b.size();
    const int INF = m+n;
    int D = band;
    vector<vector<int>> dp(m+1, vector<int>(n+1, INF));
    dp[0][0] = 0;
    for (int i = 0; i <= m; ++i) {
        for (int j = max(0, i-D); j <= min(n, i+D); ++j) {
            if (i>0 && j>0) dp[i][j] = min(dp[i][j], dp[i-1][j-1]+(a[i-1]!=b[j-1]));
            if (i>0) dp[i][j] = min(dp[i][j], dp[i-1][j]+1);
            if (j>0) dp[i][j] = min(dp[i][j], dp[i][j-1]+1);
        }
    }
    return dp[m][n];
}

size_t lsh_hash(const string &s, int lsh_bit=10) {
    size_t h = 5381; //임의로 지정
    for (char c : s) h = ((h << 5) + h) + c;
    return h >> (sizeof(size_t)*8 - lsh_bit);
}

// main
int main() {
    // ---- [1] 기본 read/해시테이블 구조 ----
    //절대경로
    string inputFile = "/Users/kosora/Documents/Xcode/C++/AlgorithmProject/AlgorithmProject/Homo_sapiens.GRCh38.dna.alt.fa";
//    string inputFile = "/Users/kosora/Documents/Xcode/C++/AlgorithmProject/AlgorithmProject/smallinput.txt";
    string outputFile = "hash_table.txt";
    string readFile = "/Users/kosora/Documents/Xcode/C++/AlgorithmProject/AlgorithmProject/smallread.txt"; //read 목록파일

    string sequence = readSequence(inputFile);
    unordered_map<int, vector<int>> hashTable = buildHashTable(sequence, patternLength);
    cout << "Rabin-Karp 기반 해시테이블 Build 완료!\n";
    saveHashTable(hashTable, outputFile);
    cout << "해시테이블 파일 저장 완료!\n";

    vector<string> reads = readReads(readFile);
    unordered_map<int, vector<int>> resultTable;
    int readIdx = 0;

    // ---- [2] 확장 파트: BloomFilter+LSH+BandDP 기반 후보군 탐색 ----
    cout << "Bloom+LSH 후보군 기반 추가탐색...\n";
    int k = 15; // seed 길이
    int bloom_size = 1000000;
    int lsh_bit = 10;
    int D = error;

    Bloom bf(bloom_size, 4);

    // Bloom Filter 삽입
    for (size_t i=0; i+ k <= sequence.size(); ++i) {
        bf.add(sequence.substr(i, k));
    }
    // Reference의 모든 seed(부분문자열)를 LSH bucket에 저장
    /*
     LSH 설명 (부족하면 말씀해주세요) -> 유사한 seed들 후보 압축해서 시간 단축.
     seed 해시값의 상위 몇 비트만 “bucket(버킷)”에 저장하여, 비슷한 k-mer(seed)들이 같은 bucket에 모으기
     모든 seed를 hash→bucket에 위치 저장
     (예: lsh_table[버킷] = {seq위치1, seq위치2, ...})
     */
    unordered_map<size_t, vector<size_t>> lsh_table;
    for (size_t i=0; i+k <= sequence.size(); ++i) {
        string seed = sequence.substr(i, k);
        size_t bkt = lsh_hash(seed, lsh_bit);
        // (최대 후보군을 너무 늘리지 않게, 한 bucket당 100개까지만 저장)
        if (lsh_table[bkt].size()<100)
            lsh_table[bkt].push_back(i);
    }

    // ---- [3] 기존 read마다 Rabin-Karp 검색 + [확장] Bloom/LSH/BandDP 추가 ----
    for (const string& read : reads) {
        // 방어: 길이 k, 100보다 짧으면 바로 skip!
            if (read.size() < 100 || read.size() < k) {
                cout << "[SKIP] readIdx=" << readIdx << ", 길이=" << read.size() << " (분석 제외)\n";
                ++readIdx;
                continue;
            }
        bool found = false;

        // --- (a) 기존 해시테이블 기반 탐색 ---
        for (int i = 0; i < 10; i++) {
            if (i*10 + 10 >= read.size()) continue; //size() failed: string index out of bounds 오류해결
            string sub = read.substr(i * 10, 10);
            int hash = 0, Dval = 1;
            for (int j = 0; j < 9; j++) Dval = (Dval * d) % q;
            for (int j = 0; j < 10; j++) hash = (d * hash + sub[j]) % q;

            if (hashTable.find(hash) != hashTable.end()) {
                for (int idx : hashTable[hash]) {
                    int start = idx - i * 10;
                    if (start < 0 || start + 100 > sequence.size()) continue;

                    string candidate = sequence.substr(start, 100);
                    if (countError(read, candidate) != -1) {
                        resultTable[readIdx].push_back(start);
                        found = true;
                    }
                }
            }
            if(found) break;
        }

        // --- (b) 추가 Bloom+LSH+BandDP 기반 탐색 ----
        // Read의 각 seed를 LSH bucket에 해싱해서 후보군 추출
        if (!found) {
            for (int j=0; j+k<=read.size(); ++j) { //out of bound 오류 해결
                string seed = read.substr(j, k);
                if (!bf.possibly_contains(seed)) continue;
                size_t bkt = lsh_hash(seed, lsh_bit);
                if (lsh_table.find(bkt) == lsh_table.end()) continue;
                for (size_t refpos : lsh_table[bkt]) {
                    if (refpos + read.size() > sequence.size()) continue;
                    string refseg = sequence.substr(refpos, read.size());
                    int edist = banded_edit_distance(read, refseg, D); // BandDP 적용
                    if (edist <= D) {
                        // 편집거리 D이내면 -> 여기서 매칭 성공
                        cout << "확장: read#" << readIdx << " (" << read.substr(0,15) << "...)"
                             << " @Pos " << refpos << " EditDist=" << edist << "\n";
                        found = true;
                        break;
                    }
                }
                if(found) break;
            }
        }
        readIdx++;
    }

    //파일 출력
    ofstream file("result.txt");
    for (auto& result : resultTable) {
        file << "Read Index: " << result.first << " -> Positions: ";
        for (int idx : result.second) {
            file << idx << " ";
        }
        file << "\n";
    }
    file.close();

    cout << "전체 프로세스 완료!\n";
    return 0;
}
