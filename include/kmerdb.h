#include<kmer_general.h>
#include<iostream>
#include<fstream>
#include<vector>
#include<string>
using namespace std;

KmerUint64Hash1* read_kmer_db_hash(string  kmer_database);

uint64_t find_kmer_db_hash(KmerUint64Hash1 *kmers , uint64_t kmer, uint32_t kmer_size);

vector<uint64_t>* read_kmer_db_vector(string  kmer_database);

uint64_t find_kmer_db_vector(vector<uint64_t>*kmers , uint64_t kmer, uint32_t kmer_size); 