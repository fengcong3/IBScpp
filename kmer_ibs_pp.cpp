///
///      @file  filter_kmers.cpp
///     @brief  Given a list of k-mers and a k-mers table output their presence/absence pattern
///
///    @author  Yoav Voichek (YV), yoav.voichek@tuebingen.mpg.de
///
///  @internal
///    Created  01/14/19
///   Compiler  gcc/g++
///    Company  Max Planck Institute for Developmental Biology Dep 6
///  Copyright  Copyright (c) 2019, Yoav Voichek
///
/// This source code is released for free distribution under the terms of the
/// GNU General Public License as published by the Free Software Foundation.
///=====================================================================================
///

#include <fstream>
#include <iostream>
#include <utility>
#include <vector>
#include <pthread.h>
#include <stdlib.h>

#include <cxxopts/include/cxxopts.hpp>

#include "kmer_general.h" // read DBs

using namespace std;

struct thread_data{
   string   kmer_talbe;
   uint64_t offset;
   uint64_t words_per_kmer;
   uint64_t kmer_number;
   uint32_t sample_size;
};

void * partial_matrix(void *threadarg)
{
	//返回这一一部分的kmer情况
	struct thread_data *my_data;
	my_data = (struct thread_data *) threadarg;
	string   kmer_talbe = my_data->kmer_talbe;
	uint64_t offset = my_data->offset;
	uint64_t words_per_kmer= my_data->words_per_kmer;
	uint64_t kmer_number= my_data->kmer_number;
	uint32_t sample_size= my_data->sample_size;
	//分配矩阵空间
	uint64_t ** p_matrix=(uint64_t **)calloc(sample_size,sizeof(uint64_t *));
	for(uint32_t i=0;i<sample_size;i++){
		p_matrix[i] = (uint64_t *)calloc(sample_size,sizeof(uint64_t));
	} 

	// cout<<kmer_talbe + string(".table")<<endl;

	ifstream table_handle(kmer_talbe , ios::binary | ios::ate); 
	if(table_handle.is_open()) {
		table_handle.seekg(offset, ios::beg); // go to offset
		uint64_t i_kt(0); // index kmers table
		vector<uint64_t> buffer(words_per_kmer+1);
		vector<uint64_t> pattern(sample_size);

		bool advance_row = true;
		while(i_kt < kmer_number) {
			if(advance_row) { 
				table_handle.read(reinterpret_cast<char *>(buffer.data()), sizeof(uint64_t)*buffer.size());
				//cout << bits2kmer31(buffer[0], 31)<<"\n";
			}
			for(size_t col_index=0; col_index<sample_size; col_index++) {
				uint64_t new_bit = buffer[(col_index >> 6) + 1] >> (col_index&(WLEN-1))  & 1ull; //得到这个样本是0 or 1
				// fout << "\t" << new_bit;
				pattern[col_index]=new_bit; //记录下来

				for(size_t sample_index1=0;sample_index1<=col_index;sample_index1++){
					// printf("%d\t%d\n",pattern[col_index],pattern[sample_index1]);
					p_matrix[col_index][sample_index1] += !(pattern[col_index] ^ pattern[sample_index1]);
				}
			}
			i_kt++;
		}


		table_handle.close();
	}else {
		cerr << "Can't open table file" << endl;
	}

	return (void *)p_matrix;

}


int main(int argc, char* argv[]) {
	cxxopts::Options options("kmer ibs", "Output the same Kmer number between two samples");
	try
	{
		options.add_options()
			("t,kmers_table", "k-mers table path", cxxopts::value<string>())
			("o,output", "output file", cxxopts::value<string>())
			("p,threads","threads",cxxopts::value<uint32_t>())
			("help", "print help")
			;
		auto result = options.parse(argc, argv);
		if (result.count("help"))
		{
			cerr << options.help() << endl;
			exit(0);
		}
		vector<string> required_parametrs({"kmers_table",  "output","threads"});
		for(size_t i=0; i<required_parametrs.size(); i++) {
			if(result.count(required_parametrs[i]) == 0) {
				cerr << required_parametrs[i] << " is a required parameter" << endl;
				cerr << options.help() << endl;
				exit(1);
			}
		}
		// Load parameters
		string fn_kmers_table(result["kmers_table"].as<string>());
		string fn_output_file(result["output"].as<string>());
		uint32_t tra_threads(result["threads"].as<uint32_t>());

		// Check if all input files exist
		vector<string> required_files({fn_kmers_table+".names",fn_kmers_table+".table"});
		for(size_t i=0; i<required_files.size(); i++) {
			if(!is_file_exist(required_files[i])) {
				cerr << "Couldn't find file: " << required_files[i] << endl;
				exit(1);
			}
		}
		/****************************************************************************************************/
		/* END of parsing and checking input parameters */
		/****************************************************************************************************/

		
		uint32_t kmer_len=31;
		
		vector<string> accession_names = load_kmers_talbe_column_names(fn_kmers_table);
		uint64_t words_per_kmer = (accession_names.size() +  WLEN - 1) / WLEN;  // 一个kmer 的p/a pattern占13个words，剩余补充为0

		ifstream table_handle(fn_kmers_table + string(".table"), ios::binary | ios::ate);  //打开table file
		if(table_handle.is_open()) {
			size_t left_in_file = table_handle.tellg();
			if(left_in_file <= (4 + 8 + 4)) {
				cerr << "table file is too small" << endl;
				return 1;
			}
			table_handle.seekg(0, ios::beg); // go to file begining

			uint32 prefix, file_kmer_len;
			uint64_t file_accession_number;
			table_handle.read(reinterpret_cast<char *>(&prefix), sizeof(prefix));   // 存储prefix  4字节
			table_handle.read(reinterpret_cast<char *>(&file_accession_number), sizeof(file_accession_number)); //存储样本数量  8字节
			table_handle.read(reinterpret_cast<char *>(&file_kmer_len), sizeof(file_kmer_len));  //存储kmer的长度 ，4字节
			left_in_file -= (sizeof(prefix) + sizeof(file_accession_number) + sizeof(file_kmer_len));

			if(prefix!=0xDDCCBBAA)  //判断头部信息是否正确
				throw std::logic_error("Incorrect prefix");
			if(file_accession_number != accession_names.size() )
				throw std::logic_error("number of accession in file not as defined in class");
			if(file_kmer_len != kmer_len)
				throw std::logic_error("kmer length in table and in list are not the same");

			// From the size of the file we can calculate the number of k_mers
			size_t size_per_kmer = sizeof(uint64_t) * (1 + words_per_kmer);
			// printf("%ull\t%ull\n",left_in_file,size_per_kmer);
			if((left_in_file % size_per_kmer) != 0)
				throw std::logic_error("size of file not valid");

			uint64_t kmer_number = left_in_file / size_per_kmer;
			cerr << "We have " << kmer_number << endl;
			//给每个线程分配任务
			vector<uint64_t> job_count(tra_threads);
			for(size_t i=0;i<job_count.size()-1;i++){
				job_count[i]=(int)(kmer_number/tra_threads);
				cout << "\t" << job_count[i]; //
			}
			job_count[job_count.size()-1] = kmer_number-(job_count[0]*(job_count.size()-1));
			cout << "\t" << job_count[job_count.size()-1] << "\n";//

			//每个线程的偏移量
			vector<uint64_t> job_offset(tra_threads);
			job_offset[0]=sizeof(file_kmer_len)+sizeof(file_accession_number)+sizeof(prefix);
			cout << "\t" << job_offset[0] ; //
			for(size_t i=1;i<job_offset.size();i++){
				job_offset[i]=job_offset[i-1] + (sizeof(uint64_t)*(words_per_kmer+1))*job_count[i-1];
				cout << "\t" << job_offset[i]; //
			}
			cout <<"\n"<<endl;

			//创建每个线程
			struct thread_data * td = (struct thread_data *)calloc(tra_threads,sizeof(thread_data));

			pthread_t * tids=(pthread_t *)calloc(tra_threads,sizeof(pthread_t));
			for(size_t i = 0; i < tra_threads; ++i)
			{
				(td+i)->kmer_number = job_count[i];
				(td+i)->kmer_talbe = fn_kmers_table + string(".table");
				(td+i)->offset = job_offset[i];
				(td+i)->sample_size = accession_names.size();
				(td+i)->words_per_kmer = words_per_kmer;
				//参数依次是：创建的线程id，线程参数，调用的函数，传入的函数参数
				int ret = pthread_create(tids+i, NULL, partial_matrix,(void *)(td+i));
				if (ret != 0)
				{
					cerr << "pthread_create error: error_code=" << ret << endl;
				}
			}

			//thread join
			vector<void *> thread_res(tra_threads);
			for(size_t i = 0; i < tra_threads; ++i)
			{
				int ret = pthread_join ( *(tids+i), &thread_res[i] );
				if (ret != 0)
				{
					cerr << "pthread_join error: error_code=" << ret << endl;
				}

			}
			

			uint64_t ibs_matrix[accession_names.size()][accession_names.size()]={0};
			//deal the res
			for(size_t i = 0; i < tra_threads; ++i)
			{
				for(size_t j =0;j < accession_names.size();j++){
					for(size_t k=0;k<accession_names.size();k++){
						ibs_matrix[j][k] += ((uint64_t **)thread_res[i])[j][k];
					}
				}
			}

			
			


			// Open output file
			ofstream fout(fn_output_file);
			if(!fout.is_open()) {
				cerr << "can't open output file " << endl;
				return 1;
			}
			fout << "Samples"; // output header
			for(size_t i=0; i<accession_names.size(); i++)
				fout << "\t" << accession_names[i];
			fout << "\n";

			// // Start reading files
			// uint64_t i_kt(0); // index kmers table
			// vector<uint64_t> buffer(words_per_kmer+1);
			// vector<uint64_t> pattern(accession_names.size());

			// uint64_t ibs_matrix[accession_names.size()][accession_names.size()]={0}; //记录两两间的距离


			// bool advance_row = true;
			// while(i_kt < kmer_number) {
			// 	if(advance_row) { 
			// 		table_handle.read(reinterpret_cast<char *>(buffer.data()), sizeof(uint64_t)*buffer.size());
			// 		advance_row = false;
			// 	}
			// 	for(size_t col_index=0; col_index<accession_names.size(); col_index++) {
			// 		uint64_t new_bit = buffer[(col_index >> 6) + 1] >> (col_index&(WLEN-1))  & 1ull; //得到这个样本是0 or 1
			// 		// fout << "\t" << new_bit;
			// 		pattern[col_index]=new_bit; //记录下来

			// 		for(size_t sample_index1=0;sample_index1<=col_index;sample_index1++){
			// 			// printf("%d\t%d\n",pattern[col_index],pattern[sample_index1]);
			// 			ibs_matrix[col_index][sample_index1] += !(pattern[col_index] ^ pattern[sample_index1]);
			// 		}
			// 	}
			// 	i_kt++;
				
			// }

			//最终输出ibs matrix
			for(size_t col_index=0; col_index<accession_names.size(); col_index++) {
				fout << accession_names[col_index];
				for (size_t row_index=0; row_index<accession_names.size(); row_index++){
					fout << "\t" << ibs_matrix[col_index][row_index];
				}
				fout << "\n";
			}

			fout.close();
			table_handle.close();
		} else {
			cerr << "Can't open table file" << endl;
			return 1;
		}
	} catch (const cxxopts::OptionException& e)
	{
		cerr << "error parsing options: " << e.what() << endl;
		cerr << options.help() << endl;
		exit(1);
	}
	return 0;
}
