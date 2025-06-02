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
//#include <cmath>
//#include <time.h>
//
////파일로만 출력 콘솔로 매번 출력 X -> O
////LSH 버킷 수 문제 100개 ? ->
////bloom + (LSH + BANDDP) 분할 적용 코드 -> O
////bloom filter 개선 -> O
//
//using namespace std;
//
//const int d = 256;  // ASCII 문자 수
////const int q = 10000019;//큰 N <  mod 값
//const int q = 99991; //작은 N < mod 값
//int patternLength = 10; //자르는 길이
//const int error = 20; //허용 오차갯수 D
//
//clock_t start, finish;
//double duration;
//// ==================== Rabin-Karp Hash 기반 ====================
//unordered_map<int, vector<int>> buildHashTable(const string& T, int m) {
//    unordered_map<int, vector<int>> hashTable;
//    int n = T.length();
//    int D = 1;
//
//    for (int i = 0; i < m - 1; i++) {
//        D = (D * d) % q;
//    }
//
//    int hash = 0;
//    for (int i = 0; i < m; i++) {
//        hash = (d * hash + T[i]) % q;
//    }
//
//    hashTable[hash].push_back(0);
//
//    for (int i = 1; i <= n - m; i++) {
//        hash = (d * (hash - T[i - 1] * D) + T[i + m - 1]) % q;
//        if (hash < 0)
//            hash += q;
//        hashTable[hash].push_back(i);
//    }
//
//    return hashTable;
//}
//
//// 염기서열 읽어서 sequence에 저장 반환
//string readSequence(const string& fileName) {
//    ifstream file(fileName);
//    string line, sequence = "";
//    while (getline(file, line)) {
//        if (line.empty() || line[0] == '>' || line[0] == 'N') continue;
//        //원본 파일에는 space, tab같은 이스케이프 문자 있어서 처리
//        for (char c : line) {
//            char uc = toupper(c);
//            if (uc == 'A' || uc == 'T' || uc == 'G' || uc == 'C')
//                sequence += uc;
//        }
//    }
//    return sequence;
//}
//
//// 해시 테이블 저장
//void saveHashTable(const unordered_map<int, vector<int>>& hashTable, const string& fileName) {
//    ofstream file(fileName);
//    for (auto iter = hashTable.begin(); iter != hashTable.end(); iter++) {
//        int hash = iter->first;
//        const vector<int>& idx = iter->second;
//        file << "Hash: " << hash << " -> Positions: ";
//        for (int i = 0; i < idx.size(); i++) {
//            file << idx[i] << " ";
//        }
//        file << "\n";
//    }
//    file.close();
//}
//
//// read 파일 읽기
//vector<string> readReads(const string& filename) {
//    ifstream file(filename);
//    vector<string> reads;
//    string line;
//    while (getline(file, line)) {
//        if (!line.empty()) reads.push_back(line);
//    }
//    return reads;
//}
//
//// mismatch 개수 계산
//int countError(const string& subSequence, const string& subRead) {
//    int count = 0;
//    for (int i = 0; i < subSequence.size(); i++) {
//        if (subSequence[i] != subRead[i]) count++;
//        if (count > error) return -1;
//    }
//    return count;
//}
//
//
//
//// main
//int main() {
//    // ---- [1] 기본 read/해시테이블 구조 ----
//    //절대경로
//    //string inputFile = "/Users/kosora/Documents/Xcode/C++/AlgorithmProject/AlgorithmProject/Homo_sapiens.GRCh38.dna.alt.fa";
//    string inputFile = "/Users/kosora/Documents/Xcode/C++/AlgorithmProject/AlgorithmProject/smallinput.txt";
//    string outputFile = "hash_table.txt";
//    string readFile = "/Users/kosora/Documents/Xcode/C++/AlgorithmProject/AlgorithmProject/smallread.txt"; //read 목록파일
//    start = clock();
//    string sequence = readSequence(inputFile);
//    unordered_map<int, vector<int>> hashTable = buildHashTable(sequence, patternLength);
//    cout << "Rabin-Karp 기반 해시테이블 Build 완료!\n";
//
//    saveHashTable(hashTable, outputFile);
//    cout << "해시테이블 파일 저장 완료!\n";
//
//    vector<string> reads = readReads(readFile);
//    unordered_map<int, vector<int>> resultTable;
//    unordered_map<int, vector<int>> resultTable2;
//    int readIdx = 0;
//
//    
//
//    
//
//    // ---- [3] 기존 read마다 Rabin-Karp 검색 + [확장] Bloom/LSH/BandDP 비교 ----
//    for (const string& read : reads) {
//        // 방어: 길이 k보다 짧으면 바로 skip!
//        if (read.size() < 100) {
//            ++readIdx;
//            continue;
//        }
//        bool found = false;
//
//        // --- 기존 해시테이블 기반 탐색 ---
//        for (int i = 0; i < 10; i++) {
//            if (i * 10 + 10 >= read.size()) continue; //size() failed: string index out of bounds 오류해결
//            string sub = read.substr(i * 10, 10);
//            int hash = 0, Dval = 1;
//            for (int j = 0; j < 9; j++) Dval = (Dval * d) % q;
//            for (int j = 0; j < 10; j++) hash = (d * hash + sub[j]) % q;
//            if (hash < 0)
//                hash += q;
//
//            if (hashTable.find(hash) != hashTable.end()) {
//                for (int idx : hashTable[hash]) {
//                    int start = idx - i * 10;
//                    if (start < 0 || start + 100 > sequence.size()) continue;
//
//                    string candidate = sequence.substr(start, 100);
//                    if (countError(read, candidate) != -1) {
//                        resultTable[readIdx].push_back(start);
//                        found = true;
//                    }
//                }
//            }
//            if (found) break;
//        }
//
//        readIdx++;
//    }
//    finish = clock();
//    duration = (double)(finish - start) / CLOCKS_PER_SEC;
//    cout << duration << "second" << endl;
//
//
//    ofstream file("result.txt");
//    for (auto& result : resultTable) {
//        file << "read index: " << result.first << " -> positions: ";
//        for (int idx : result.second) {
//            file << idx << " ";
//        }
//        file << "\n";
//    }
//    file.close();
//    cout << "전체 프로세스 완료!\n";
//
//    return 0;
//}
