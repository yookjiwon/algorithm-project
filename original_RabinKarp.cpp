#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <map>
#include <vector>
#include <time.h>
using namespace std;

const int d = 256;  // ASCII 문자 수
const int q = 99991;//작은 N
int patternLength = 10;//자르는 길이  
const int error = 45;//허용 오차
const int readLength = 500; //read 길이
clock_t start, finish;


unordered_map<int, vector<int>> buildHashTable(const string& T, int m) {
    unordered_map<int, vector<int>> hashTable;
    int n = T.length();
    int D = 1;

    int h = 1;
    for (int i = 0; i < patternLength - 1; i++)
        h = (h * d) % q; // d^(m-1)

    int hash = 0;
    for (int i = 0; i < patternLength; i++)
        hash = (d * hash + T[i]) % q;

    hashTable[hash].push_back(0);
    cout << hash << endl;

    for (int i = 1; i <= n - patternLength; i++) {
        hash = (d * (hash - T[i - 1] * h % q + q) + T[i + patternLength - 1]) % q;
        hashTable[hash].push_back(i);
        //cout << i << " " << hash << endl;
    }


    return hashTable;
}

// 염기서열 읽어서 sequence에 저장 반환
string readSequence(const string& fileName) {
    ifstream file(fileName);

    string line, sequence = "";
    while (getline(file, line)) {
        if (line.empty() || line[0] == '>' || line[0] == 'N') continue;
        sequence += line;
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
        if (!line.empty()) {
            reads.push_back(line);
        }
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

int main() {
    //string inputFile = "Homo_sapiens.GRCh38.dna.alt.fa";
    string inputFile = "10M.txt";
    string readFile = "10M_reads.txt";
    string outputFile = "hash_table.txt";
    

    start = clock();

    //입력 자료 읽기
    string sequence = readSequence(inputFile);
    unordered_map<int, vector<int>> hashTable = buildHashTable(sequence, patternLength);
    cout << "테이블 생성 완료" << endl;

    
    //saveHashTable(hashTable, outputFile);
    //cout << "해시 테이블 저장 완료" << endl;


    vector<string> reads = readReads(readFile);

    map<int, vector<int>> resultTable;

    int readIdx = 0;

    //해시 테이블 기반으로 실제 값 비교
    for (const string& read : reads) {
        cout << readIdx << "번쨰" << endl;
        bool found = false;
        for (int i = 0; i <= read.size() - patternLength; i += patternLength) {
            cout << i << endl;
            string sub = read.substr(i, patternLength);
            int hash = 0;
            for (int i = 0; i < patternLength; i++) {
                hash = (d * hash + sub[i]) % q;
            }
            cout<< hash << endl;
            if (hashTable.find(hash) != hashTable.end()) {
                for (int idx : hashTable[hash]) {
                    int start = idx - i;
                    if (start < 0 || start + 100 > sequence.size()) continue;

                    string candidate = sequence.substr(start, readLength);

                    if (countError(read, candidate) != -1) {
                        resultTable[readIdx].push_back(start);
                        found = true;
                        break;
                    }
                }
            }
            if (found) break;

        }
        readIdx++;
    }

    // 결과 저장
    ofstream file("result.txt");
    for (auto& result : resultTable) {
        file << "Read Index: " << result.first << " -> Positions: ";
        for (int idx : result.second) {
            file << idx << " ";
        }
        file << "\n";
    }
    file.close();
    
    finish = clock();
    double timeTaken = double(finish - start) / CLOCKS_PER_SEC;
    cout << "총 소요 시간: " << timeTaken << "초" << endl;

    return 0;
}