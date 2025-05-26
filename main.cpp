#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;

const int d = 256;  // ASCII 문자 수
//const int q = 10000019;//큰 N
const int q = 99991;//작은 N
int patternLength = 10;//자르는 길이  
const int error = 5;//허용 오차


// Rabin-Karp 해시 테이블 생성
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
        if (subSequence[i] != subSequence[i]) count++;
        if (count > error) return -1;
    }
    return count;
}

int main() {
    //string inputFile = "Homo_sapiens.GRCh38.dna.alt.fa";
    string inputFile = "input.txt";//과제에서 쓴 자료
    string outputFile = "hash_table.txt";                        

    //입력 자료 읽기
    string sequence = readSequence(inputFile);
    unordered_map<int, vector<int>> hashTable = buildHashTable(sequence, patternLength);
    cout << "테이블 생성 완료" << endl;
    
    //생략가능- 자료크기 크면 너무 오래걸림
    saveHashTable(hashTable, outputFile);
    cout << "해시 테이블 저장 완료" << endl;

    //read길이는 100으로 설정
    string readFile = "reads_1000.txt";
    //string readFile = "smallread.txt";

    vector<string> reads = readReads(readFile);

    unordered_map<int, vector<int>> resultTable;

    int readIdx = 0;

    //해시 테이블 기반으로 실제 값 비교
    for (const string& read : reads) {
        bool found = false;
        for (int i = 0; i < 10; i++) {
            string sub = read.substr(i * 10, 10);
            int hash = 0, D = 1;
            for (int j = 0; j < 9; j++) D = (D * d) % q;
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

    return 0;
}