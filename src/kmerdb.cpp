#include<kmerdb.h>
#include<algorithm>
inline tuple<uint64_t, uint64_t> is_canonized_kmer_representation_flag(const uint64_t& k, const uint32& k_len) {
	uint64_t k_rc = kmer_reverse_complement(k, k_len); 
	if(k<k_rc) 
		return make_tuple(k,    0x4000000000000000);
	else
		return make_tuple(k_rc, 0x8000000000000000);
}

KmerUint64Hash1* read_kmer_db_hash(string  kmer_database){
    // KmerUint64Hash1 kmerss(10000);
    KmerUint64Hash1 * kmerss = new KmerUint64Hash1(1000);
    (*kmerss).set_empty_key(NULL_KEY);

    ifstream kmer_database_handle(kmer_database, ios::binary); 
    if(kmer_database_handle.is_open()) {
        // cout<<333<<endl;
        uint64_t tmp_kmer_inf = 0;
        uint64_t c=0;
        while (1){
            kmer_database_handle.read(reinterpret_cast<char *>(&tmp_kmer_inf), sizeof(tmp_kmer_inf));
            // cout << bits2kmer31(tmp_kmer_inf& 0x3FFFFFFFFFFFFFFF , 31) <<endl;

            if(!kmer_database_handle.eof()){
                c++;
                // if (c%10000 == 0){
                //     cerr << c << endl;
                // }
                (*kmerss).insert(KmerUint64Hash1::value_type(tmp_kmer_inf& 0x3FFFFFFFFFFFFFFF,tmp_kmer_inf>>62));
            }
            else{
                break;
            }
        }
    }
    cerr << "total kmers: "<< (*kmerss).size()<<endl; 
    kmer_database_handle.close();
    return kmerss;
}

uint64_t find_kmer_db_hash(KmerUint64Hash1 *kmers , uint64_t kmer, uint32_t kmer_size){
    uint64_t k, flag;
    uint64_t returnv = 0;
    tie(k, flag)  = is_canonized_kmer_representation_flag(kmer, kmer_size);
    KmerUint64Hash1::iterator it_hash;
    it_hash = (*kmers).find(k);
    if(it_hash != (*kmers).end()) {
        // returnv = (it_hash->second) & (flag>>62);// return true only the kmer is same
        returnv = 1;  // return true if find Inverse complementary sequence or himself
    }
    // cerr<<returnv<<endl;
    return returnv;
}

vector<uint64_t>* read_kmer_db_vector(string  kmer_database){
    ifstream kmers_handle(kmer_database, ios::binary | ios::ate); 
    uint64_t kmers_count = 0;
    if(kmers_handle.is_open()) {
		uint64_t left_in_file = kmers_handle.tellg();
        kmers_count = left_in_file/sizeof(uint64_t);
        cerr << "total kmers: "<< kmers_count <<endl; 
        kmers_handle.seekg(0, ios::beg);
        // vector<uint64_t> * kmers_db = new vector<uint64_t> (kmers_count);
        vector<uint64_t> * kmers_db = new vector<uint64_t> ();
        cerr << "calloc mem " <<endl; 
        cerr << "total kmers: "<< kmers_db->size() <<endl; 
        uint64_t tmp_kmer_inf = 0;
        uint64_t c = 0;
        while(1){
            kmers_handle.read(reinterpret_cast<char *>(&tmp_kmer_inf), sizeof(tmp_kmer_inf));
            // cout << bits2kmer31(tmp_kmer_inf& 0x3FFFFFFFFFFFFFFF , 31) <<endl;

            if(!kmers_handle.eof()){
                // if (c%10000 == 0){
                //     cerr << c << endl;
                // }
                // (*kmers_db)[c]=tmp_kmer_inf;

                kmers_db->push_back(tmp_kmer_inf);
                c++;
            }
            else{
                break;
            }
        }
        cerr << "total kmers: "<< kmers_db->size()<<endl; 
        kmers_handle.close();
        return kmers_db;
    }else{
        cerr<< "cant open kmers_with_strand"<<endl;
        exit(-1);
    }

}

int kmer_comp(const void * a, const void * b){
    uint64_t ka = *(uint64_t * )a;
    uint64_t kb = *(uint64_t * )b;

    return ( (ka&0x3FFFFFFFFFFFFFFF) > (kb&0x3FFFFFFFFFFFFFFF) ) - ( (ka&0x3FFFFFFFFFFFFFFF) < (kb&0x3FFFFFFFFFFFFFFF) ) ;
}



uint64_t find_kmer_db_vector(vector<uint64_t> *kmers , uint64_t kmer, uint32_t kmer_size){
    uint64_t k, flag;
    tie(k, flag)  = is_canonized_kmer_representation_flag(kmer, kmer_size);
    
    uint64_t ret =(uint64_t ) binary_search(kmers->begin(), kmers->end(), k, 
        [](const uint64_t & a, const uint64_t & b) -> bool
		{return (a&0x3FFFFFFFFFFFFFFF) < (b&0x3FFFFFFFFFFFFFFF);});

    return ret;
}