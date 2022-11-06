// ------------------------------------------------------------------------
//
//        File: bch4836.c 
// Description: Encoder/Decoder for a (48, 36, 5) binary BCH code 
//              This programs illustrates the use of PGZ decoder for t=2
//
// ------------------------------------------------------------------------
// This program is complementary material for the book:
//
// R.H. Morelos-Zaragoza, The Art of Error Correcting Coding, Wiley, 2002.
//
// ISBN 0471 49581 6
//
// This and other programs are available at http://the-art-of-ecc.com
//
// You may use this program for academic and personal purposes only. 
// If this program is used to perform simulations whose results are 
// published in a journal or book, please refer to the book above.
//
// The use of this program in a commercial product requires explicit
// written permission from the author. The author is not responsible or 
// liable for damage or loss that may be caused by the use of this program. 
//
// Copyright (c) 2002. Robert H. Morelos-Zaragoza. All rights reserved.
// ------------------------------------------------------------------------

#include <iostream>
#include <algorithm>
#include <vector>

using namespace std;

int m = 6, n = 63, k = 36, t = 2, d = 5;
int	length = 48;
int p[7]; //irreducible polynomial 
int alpha_to[64], index_of[64], g[13];
int recd[48], _data[36], bb[13];
int numerr, errpos[64], decerror = 0;
int seed;

void read_p(){
// Primitive polynomial of degree 6 
    p[0] = p[1] = p[6] = 1;
    p[2] = p[3] = p[4] = p[5] =0;
}

void generate_gf() {
/*
	generate GF(pow(2, m)) ((2**m)) from the irreducible polynomial p(X) in p[0]..p[m]
	lookup tables:  index->polynomial form   alpha_to[] contains j=alpha**i;
	polynomial form -> index form  index_of[j=alpha**i] = i alpha=2 is the
	primitive element of GF(2**m) 
*/
	int mask = 1;
	alpha_to[m] = 0; //m = 6
	for (int i = 0; i < m; i++) {
		alpha_to[i] = mask;
		index_of[alpha_to[i]] = i;
		if (p[i] != 0)
			alpha_to[m] ^= mask;
		mask <<= 1;
	}
	index_of[alpha_to[m]] = m;
	mask >>= 1;
	for (int i = m + 1; i < n; i++) {
		if (alpha_to[i - 1] >= mask)
		  alpha_to[i] = alpha_to[m] ^ ((alpha_to[i - 1] ^ mask) << 1);
		else
		  alpha_to[i] = alpha_to[i - 1] << 1;
		index_of[alpha_to[i]] = i;
	}
	index_of[0] = -1;
}



void gen_poly() {
/*  
	Compute generator polynomial of BCH code of length = 48, redundancy = 12
*/
	vector <vector <int>> cycle_coset;
	vector <int> allnumbers, zeros, temp_coset_index;
	vector <int>::iterator it;

	// Generate cycle sets modulo 63
	for (int i = 0, j; i <= 31; i++){
		for (j = 0; ; j++){
			if (j == 0)
				temp_coset_index.push_back(i);
			else
				temp_coset_index.push_back((temp_coset_index[j - 1] *2 ) % n);
			int last_element = temp_coset_index[j];
			it = (find(allnumbers.begin(), allnumbers.end(), last_element));
			if (it != allnumbers.end()){
				temp_coset_index.clear();
				break;
			}
			allnumbers.push_back(last_element);
			if (temp_coset_index[0] == (temp_coset_index[j] * 2 ) % n){
				cycle_coset.push_back(temp_coset_index);
				temp_coset_index.clear();
				break;
			}
		}
	}
	allnumbers.clear();

	int rdncy = 0, size = 0;	
	// Search for roots 1, 2, ..., d-1 in cycle sets (d = 5)
	for (const auto& index : cycle_coset){
		for (const auto& second_index : index){
			for (int root = 1; root < d; root++)
				if (root == second_index)
					size = index.size();
		}
		rdncy += size;
		if(size != 0) {
		//populate zeros with cosets that have roots 1 - d-1
			for (const auto& second_index : index)
				zeros.push_back(second_index);
		}
		size = 0;
	}

	cout<<"This is a ("<<length<<","<<k<<","<<d<<") binary BCH code"<<endl;
	// Compute generator polynomial 
	g[0] = alpha_to[zeros[0]];
	g[1] = 1;		// g(x) = (X + zeros[0]) initially
	for (int i = 2; i <= rdncy; i++) {
	  g[i] = 1;
	  for (int j = i - 1; j > 0; j--)
	    if (g[j] != 0)
	      g[j] = g[j - 1] ^ alpha_to[(index_of[g[j]] + zeros[i-1]) % n];
	    else
	      g[j] = g[j - 1];
	  g[0] = alpha_to[(index_of[g[0]] + zeros[i-1]) % n];
	}
	cout<<"g(x) = ";
	for (int i = 0; i <= rdncy; i++) {
	  cout<<g[i];
	}
	cout<<endl;
}


void encode_bch() {
/* 
	Calculate redundant bits bb[], codeword is c(X) = _data(X)*X**(n-k)+ bb(X)
*/
	int i, j;
	int feedback;
	for (i = 0; i < length - k; i++)
		bb[i] = 0;
	for (i = k - 1; i >= 0; i--) {
		feedback = _data[i] ^ bb[length - k - 1];
		if (feedback != 0) {
			for (j = length - k - 1; j > 0; j--)
				if (g[j] != 0)
					bb[j] = bb[j - 1] ^ feedback;
				else
					bb[j] = bb[j - 1];
			bb[0] = g[0] && feedback;
		} else {
			for (j = length - k - 1; j > 0; j--)
				bb[j] = bb[j - 1];
			bb[0] = 0;
		}
	}
}


void decode_bch() {
/*
	We do not need the Berlekamp algorithm to decode.
	We solve before hand two equations in two variables.
*/
	int i, j, q;
	int elp[3] = {0}, s[5] = {0}, s3;
	int count = 0, syn_error = 0;
	int loc[3] = {0}, err[3] = {0}, reg[3] = {0};
	int	aux;
	// first form the syndromes 
	cout<<endl<<"s[] = (";
	for (i = 1; i <= 4; i++) {
		s[i] = 0;
		for (j = 0; j < length; j++)
			if (recd[j] != 0)
				s[i] ^= alpha_to[(i * j) % n];
		if (s[i] != 0)
			syn_error = 1;	/* set flag if non-zero syndrome 
								NOTE: If only error detection is needed,
								then exit the program here...
							*/
		// convert syndrome from polynomial form to index form  
		s[i] = index_of[s[i]];
		if (i<4)
			cout<<s[i]<<"  ";
		else
			cout<<s[i]<<")"<<endl;
	}
	if (syn_error) {	// If there are errors, try to correct them 
		if (s[1] != -1) {
			s3 = (s[1] * 3) % n;
			if ( s[3] == s3 )  // Was it a single error ? 
				{
				cout<<"One error at "<<s[1];
				recd[s[1]] ^= 1;		// Yes: Correct it 
				}
			else {				/* Assume two errors occurred and solve
									for the coefficients of sigma(x), the
									error locator polynomail
								*/
                if (s[3] != -1)
                  aux = alpha_to[s3] ^ alpha_to[s[3]];
                else
                  aux = alpha_to[s3];

				elp[0] = 0;
				elp[1] = (s[2] - index_of[aux] + n) % n;
				elp[2] = (s[1] - index_of[aux] + n) % n;
				cout<<"Sigma(x) = ";
				for (i = 0; i <= 2; i++)
					cout<<elp[i]<<" ";
				cout<<endl<<"Roots: ";
				// find roots of the error location polynomial 
				for (i = 1; i <= 2; i++)
					reg[i] = elp[i];
				count = 0;
				for (i = 1; i <= 63; i++) { // Chien search 
					q = 1;
					for (j = 1; j <= 2; j++)
						if (reg[j] != -1) {
							reg[j] = (reg[j] + j) % n;
							q ^= alpha_to[reg[j]];
						}
					if (!q) {	// store error location number indices 
						loc[count] = i % n;
						count++;
						cout<<(i%n)<<" ";
					}
				}
				if (count == 2)	
				// no. roots = degree of elp hence 2 errors 
					for (i = 0; i < 2; i++)
						recd[loc[i]] ^= 1;
				else	// Cannot solve: Error detection 
					cout<<endl<<"Incomplete decoding";
				}
			}
		else if (s[2] != -1) // Error detection 
			cout<<endl<<"Incomplete decoding";
	}
}

bool check_input(int input, int boundry_1, int boundry_2){
	if (input < boundry_1 || input > boundry_2){
		return true;
	}
	return false;
}

int main() {
	int i;
	read_p();		// read generator polynomial g(x) 
	generate_gf();	// generate the Galois Field GF(2**m) 
	gen_poly();		// Compute the generator polynomial of BCH code 

	seed = 1;
	srandom(seed);
	// Randomly generate _data 
	for (i = 0; i < k; i++)
		_data[i] = (random() & 67108864) >> 26;

	// ENCODE 
	encode_bch();	// encode _data 
 
	for (i = 0; i < length - k; i++)
		recd[i] = bb[i];	// first (length-k) bits are redundancy 
	for (i = 0; i < k; i++)
		recd[i + length - k] = _data[i];	// last k bits are _data 
	cout<<"c(x) = ";
	for (i = 0; i < length; i++) {
		cout<<recd[i];
	}

	// ERRORS
    cout<<endl<<"Enter the number of errors (this program can fix only two errors): "<<endl;
    cin>>numerr;
	while (check_input(numerr, 1, 10)){
		cout<<"Wrong number of errors:"<<endl;
		cin>>numerr;
	}

	cout<<"Enter the position of errors: "<<endl;
	for (i = 0; i < numerr; i++) {
		cin>>errpos[i];
		while (check_input(errpos[i], 0, 36)){
			cout<<"Wrong error position (choose a number between 0 and 36):"<<endl;
			cin>>errpos[i];
		} 
		recd[errpos[i]] ^= 1;
	}
	cout<<"r(x) = ";
	for (i = 0; i < length; i++)
		cout<<recd[i];

    // DECODE
	decode_bch();
	/*
	  print out original and decoded _data
	*/
	cout<<endl<<"Results: "<<endl;
	cout<<"Original _data  = ";
	for (i = 0; i < k; i++)
		cout<<_data[i];
	cout<<endl<<"Recovered _data = ";
	for (i = length - k; i < length; i++)
		cout<<recd[i];
	/* decoding errors: we compare only the _data portion */
	for (i = length - k; i < length; i++)
		if (_data[i - length + k] != recd[i])
			decerror++;
	if (decerror)
		cout<<endl<<decerror<<" Message decoding errors"<<endl;
	else
		cout<<endl<<"Succesful decoding"<<endl;
    return 0;
}