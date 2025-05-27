#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;

const int d = 256;  // ASCII ���� ��
//const int q = 10000019;//ū N
const int q = 99991;//���� N
int patternLength = 10;//�ڸ��� ����  
const int error = 5;//��� ����


// Rabin-Karp �ؽ� ���̺� ����
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

// ���⼭�� �о sequence�� ���� ��ȯ
string readSequence(const string& fileName) {
    ifstream file(fileName);

    string line, sequence = "";
    while (getline(file, line)) {
        if (line.empty() || line[0] == '>' || line[0] == 'N') continue;
        sequence += line;
    }
    return sequence;
}

// �ؽ� ���̺� ����
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

// read ���� �б�
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

// mismatch ���� ���
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
    string inputFile = "input.txt";//�������� �� �ڷ�
    string outputFile = "hash_table.txt";                        

    //�Է� �ڷ� �б�
    string sequence = readSequence(inputFile);
    unordered_map<int, vector<int>> hashTable = buildHashTable(sequence, patternLength);
    cout << "���̺� ���� �Ϸ�" << endl;
    
    //��������- �ڷ�ũ�� ũ�� �ʹ� �����ɸ�
    saveHashTable(hashTable, outputFile);
    cout << "�ؽ� ���̺� ���� �Ϸ�" << endl;

    //read���̴� 100���� ����
    string readFile = "reads_1000.txt";
    //string readFile = "smallread.txt";

    vector<string> reads = readReads(readFile);

    unordered_map<int, vector<int>> resultTable;

    int readIdx = 0;

    //�ؽ� ���̺� ������� ���� �� ��
    for (const string& read : reads) {
        bool found = false;
        for (int i = 0; i < patternLength; i++) {
            string sub = read.substr(i * patternLength, patternLength);
            int hash = 0, D = 1;
            for (int j = 0; j < patternLength - 1; j++) D = (D * d) % q;
            for (int j = 0; j < patternLength; j++) hash = (d * hash + sub[j]) % q;

            if (hashTable.find(hash) != hashTable.end()) {
                for (int idx : hashTable[hash]) {
                    int start = idx - i * patternLength;
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

    // ��� ����
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