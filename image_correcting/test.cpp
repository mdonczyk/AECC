#include <iostream>
#include <random>
#include <fstream>
#include <bitset>
#include <memory>
using namespace std;

#define k 51
#define HEADER_BYTES 30
#define ERROR_PROBABILITY 1000 // 100 is 1% error probability, 1000 is 0.1% ... 100000 is 0.001%

void introduceErrors(const int &probability, char &buffer,const int &bufferSize) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> distribution(0, probability);
    for (int j = bufferSize-1; j >= 0; j--) {
        int randnum = distribution(gen);
        if (randnum < 1) {
            buffer ^= (1<<j);
        }
    }
}

vector <bitset <k>> bits_to_bitsets (const vector <bool> &buffer_bits) {
    // put each bit from bool vector to a bitset and then the full bitset to a vector:
    vector <bitset <k>> vector_of_bits;
    bitset <k> bits;
    int it = 0;
    for (int i=0; i<buffer_bits.size(); i++) {
        bits[it] = buffer_bits[i];
        it++;
        if (it == k || i == buffer_bits.size()-1) {
            vector_of_bits.push_back(bits);
            it = 0;
            bits = 0;
        }
    }
    return vector_of_bits;
}

vector <bool> bytes_to_bits (const char *buffer, const int fileSize) {
    // put each bit of the char stream into a bool vector:
    vector <bool> buffer_bits;
     for (int i = HEADER_BYTES; i <fileSize; i++) {
        int mask = 0;
        char temp_char = 0;
            for (int j=0; j<8; j++) {
            bool is_set = buffer[i] & (1 << mask);
            buffer_bits.push_back(is_set);
            mask ++;
        }
    }
    return buffer_bits;
}

vector <bool> bitset_to_bits (const vector <bitset <k>> &vector_of_bits) {
// put all bits from bitset into a bool vector:
    vector <bool> recovered_bits;
    int it = 0;
    bool break_out = false;
    for (auto const &bits : vector_of_bits) {
        for (int j=0; j<k; j++) {
            recovered_bits.push_back(bits[j]);
            it ++;
        }
    }
    return recovered_bits;
}

vector <char> bits_to_bytes (const vector <bool> &recovered_bits, const int fileSize) {
// put bits from buffer_bits into vector of chars:
    int it = 0;
    vector <char> charstream;
    for (int i = HEADER_BYTES; i <fileSize; i++) {
        char temp_char = 0;
            for (int j=0; j<8; j++) {
            temp_char ^= (recovered_bits[it] << j);
            it ++;
        }
        introduceErrors(ERROR_PROBABILITY, temp_char, 8);
        charstream.push_back(temp_char);
    }
    return charstream;
}

int main() {
    ifstream file1("image.bmp", ios::binary);
    ofstream file2;
    remove("image2.bmp");
    file2.open ("image2.bmp", ios::out | ios::app | ios::binary);
    // find the number of bytes that are in the file:
    file1.seekg(0, ios::end);
    int fileSize = file1.tellg();
    file1.seekg(0, ios::beg);
    // allocate enough memory
    char *buffer = new char[fileSize];
    // read the bytes to buffer:
    file1.read(buffer, fileSize);
    // write all header bytes to new file without modification:
    for (int i = 0; i < HEADER_BYTES; i++) {
        file2 << buffer[i];
    }

    // put each bit of the char stream into a bool vector:
    vector <bool> buffer_bits = bytes_to_bits(buffer, fileSize);

    // put each bit from bool vector to a bitset and then the full bitset to a vector:
    vector <bitset <k>> vector_of_bits = bits_to_bitsets(buffer_bits);

    // put all bits from bitset into a bool vector:
    vector <bool> recovered_bits = bitset_to_bits(vector_of_bits);

    // put bits from buffer_bits into vector of chars:
    vector <char> charstream = bits_to_bytes(recovered_bits, fileSize);

    // put the rest of modified bytes to new file:
    for (auto temp_char : charstream) {
        file2 << temp_char;
    }

    // delete allocated memory:
    delete[] buffer;
    // close opened files:
    file1.close();
    file2.close();
    return 0;
}