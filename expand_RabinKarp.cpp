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

using namespace std;

const int d = 256;  // ASCII 문자 수
const int q = 10000000;//큰 N
const int error = 2; //허용 오차갯수 D

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

// [참조 내용 + 개선 아이디어 (비트방식으로 저장)] Bloom Filter
/*
 Bloom Filter 설명 (부족하면 말씀해주세요)
 후보군(매칭 가능성 있는 시드)만 남기고, 나머지는 “비교” 자체를 건너뛸 수 있음

 배열크기 크게 늘려보고 리드가 있나 없나만 판단하는거 아닌가 ?
 메모리를 아끼는 방법 -> 비트 배열로 해결해봄!
 */
struct FastBloom {
    vector<uint64_t> bits; // filter 저장 역할
    size_t size, hash_num; // 비트 배열 크기, 사용할 해시 함수 개수

    FastBloom(size_t n_bits, int k_hash)
        : size(n_bits), hash_num(k_hash), bits((n_bits+63)/64, 0) {}

    // 단일 seed 문자열 + 해시 시드(seed)로 bit index 계산용 hash 함수
    size_t hash(const string &s, int seed) const {
        // 해시값은 참조를 통해 구했습니다.
        uint64_t h = 0xcbf29ce484222325 ^ seed;
        for (char c : s) h = (h ^ c)*0x100000001b3;
        return h % size;
    }

    // 블룸필터에 시드 삽입
    void add(const string &s) {
        // 여러 해시함수를 돌려 각각의 비트 위치를 켜줌
        for (int i=0; i<hash_num; ++i) {
            size_t idx = hash(s, i);
            bits[idx/64] |= (1ULL<<(idx%64));
        }
    }
    bool possibly_contains(const string &s) const {
        for (int i=0; i<hash_num; ++i) {
            size_t idx = hash(s, i);
            //0 있으면 false로 존재 가능여부 파단
            if (!(bits[idx/64] & (1ULL<<(idx%64)))) return false;
        }
        // 모두 1이면 true
        return true;
    }
};


//Bloom filter 사이즈 설정
size_t estimate_bloom_size(size_t n_items, double fpr=0.02) {
    return static_cast<size_t>(-(double)n_items * log(fpr) / pow(log(2),2));
}
int recommend_hash_num(size_t bloom_size, size_t n_items) {
    return max(2, int((double)bloom_size / n_items * log(2) + 0.5));
}

// [참조 내용 + 개선 아이디어 (오차를 이용하고자)] BandDP
/*
    오차가 주어지고 시작하니까 오차를 이용해서 비교시간을 줄이고자함. <- 아이디어
    D를 기준으로 그 근방의 값만 비교하자
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
// [참조 내용] LSH
size_t lsh_hash(const string &s, int lsh_bit=10) {
    size_t h = 5381; //임의로 지정 << 벤치마킹 때에서
    for (char c : s) h = ((h << 5) + h) + c;
    return h >> (sizeof(size_t)*8 - lsh_bit);
}

// main
int main() {
    // 실행 시간 측정 시작
    clock_t start = clock();

    // 결과 저장을 위한 테이블/버퍼
    unordered_map<int, vector<int>> resultTable; // read별 매칭 포지션 기록
    vector<string> result_buf; // 파일로 쓸 결과 문자열 버퍼

    // 입력(reference/read) 파일명 설정 및 데이터 로딩
    // Xcode 환경에서 경우 절대경로로 호출
    string inputFile = "/Users/kosora/Documents/Xcode/C++/AlgorithmProject/AlgorithmProject/smallinput.txt";
    string readFile = "/Users/kosora/Documents/Xcode/C++/AlgorithmProject/AlgorithmProject/smallread.txt";
    string sequence = readSequence(inputFile);
    vector<string> reads = readReads(readFile);

    // 변수들 초기화
    int k = 15;            // seed(k-mer) 길이
    int lsh_bit = 10;      // LSH 해시 상위 비트 수(버킷 수=1024)
    int D = error;         // 허용 오차(편집거리)

    // Bloom 필터 적용
    size_t bloom_size = estimate_bloom_size(sequence.size(), 0.02);
    int hash_num = recommend_hash_num(bloom_size, sequence.size());
    FastBloom bf(bloom_size, hash_num);

    for (size_t i = 0; i + k <= sequence.size(); ++i) {
        bf.add(sequence.substr(i, k));
    }

    // LSH 버킷수 설정
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

    // 각 read별로 Bloom+LSH+BandDP 매칭 시도
    int readIdx = 0;
    for (const string& read : reads) {
        // read가 너무 짧으면 스킵 (seed k, 오차 D보다 짧을 때)
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
                
                // BandDP를 D/2 band로 먼저 시도, 실패시 D로 확장
                // 과제를 해결하는데 있어 D가 1이라면 기존 라빈카프 countError 함수방식과 동일.
                int edist = banded_edit_distance(read, refseg, D / 2);
                if (edist > D) {
                    edist = banded_edit_distance(read, refseg, D);
                }

                if (edist <= D) {
                    resultTable[readIdx].push_back(refpos);
                    found = true;
                    break;
                }
            }
            if (found) break;
        }
        readIdx++;
    }

    // 결과 출력
    for (const auto& [readIdx, positions] : resultTable) {
        const string& read = reads[readIdx];
        string line = "read index: " + to_string(readIdx) + " -> positions: ";
        for (int pos : positions) {
            string refseg = sequence.substr(pos, read.size());
            int edist = banded_edit_distance(read, refseg, D); // D-band로 오차 계산
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

    // 수행시간 출력
    clock_t finish = clock();
    double duration = double(finish - start) / CLOCKS_PER_SEC;
    cout << "전체 프로세스 완료!\n";
    cout << "총 소요 시간: " << duration << "초\n";
    
    return 0;
}
