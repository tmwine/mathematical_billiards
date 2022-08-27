/*

	modification of MultAlpha_Fast_File:
		while MultAlpha_Fast_File only took a file as a single (potentially very long) symbol string, this (_Batch) allows batching of a text file with symbol strings on each line, developing an intermediary array of the log entrywise norm for each line's symbol string; the array is saved in an output file, one log||.|| entry on each line
		note, this expects only symbols in accepted range (numbers 1 to k_sym) on each line of the file
		note, this ~latches onto the input symbol file (eg <path>/string_symbols.txt) and saves the log||.|| output doubles, one on each line, in a newly-made file <path>/string_symbols_output.txt)

	fast matrix multiplication, for alphabets with > 2 symbols

	this expects command line inputs: <Mi*P data filename> <symbol string file>
	<Mi*P data filename> is the appropriate .bin file containing the relevant Mi*P set
	<symbol string> should be a file with a string of numbers on each line, each in the range [1,...,k_sym]; nb: this expects the values in symbol_string to be in the valid range--it does not check this beforehand--if there is an error, it will give an error message
	nb: as of latest, this assumes a consecutive mapping between chars '1','2',...'9' and their representative (byte) values--this to quickly convert a char decimal digit into its value (through char-'1')

	this reads in a standard format M_i*P matrix data file:
	<k_sym>, <di_vals>, <hk>, <delta>, <p_and_s>, <nbrows>, <mxsz_vals>, <MiP>
	k_sym = number of symbols
	di_vals = shift values (d_1,...,d_k)
	hk = simplex height
	delta = pad/splice penalty parameter, for degree of regularity/randomness 
	p_and_s = 0 for pad only, 1 for pad and splice, 2 for splice only
	nbrows = number of block rows (=1 for binary alphabets)
	mxsz_vals = sizes of the matrices in respective block rows (array of elements nbrows long)
	MiP = matrix data; it's stored in order MiP_1,1,...,MiP_1,k; MiP_2,1,...,MiP_2,k; ...; ie the first block row, followed by the second block row, and so on
	nb: the individual MiP matrix entries are expected to be written in row-major order (that is, (r1,c1),(r1,c2),...,(r2,c1),(r2,c2),...,...)

	nb: this should compile on C++11 or later; eg $g++ -std=c++11 my_code.cpp -o my_code

*/


using namespace std; // this is generally bad practice

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <string.h>
#include <sys/time.h>
#include <iomanip> // for setprecision (at least)
#include <climits> // for checking bits in a byte
#include <fstream> // for string file access
#include <vector> // for storing output log||.|| values
#include <iterator> // for saving output log||.|| values


#define PAD_ONLY 0
#define PAD_AND_SPLICE 1
#define SPLICE_ONLY 2


void cout_vec(double *vec, int s); // for debugging



int main(int argc,char *argv[])
{

	if (argc != 3) {
		cout << "Problem with input arguments." << endl;
		exit(1);
	}

	char *MiPfile = argv[1];
	char *symbol_string_file = argv[2];

	// check C++ here is using 8 bit bytes, and 8 bytes to a "double"
	if (CHAR_BIT!=8) {
		cout << "Number of bits in a byte here is " << CHAR_BIT << ", not 8. This will likely produce problems with reading the Octave-generated M_i*P bin data file, which renders doubles with 64 bits. Exiting." << endl;
		exit(1);
	} else if (sizeof(double)!=8) {
		cout << "C++ here defines type double as, " << sizeof(double) << " bytes. This will likely produce problems reading the Octave-generated M_i*P bin data file, which renders doubles with 64 bits. Exiting." << endl;
	}



	//***
	// read in .bin file header
	//***

	// get prelim parameters for matrix set
	// the matrix data file has (in order) {matrices' size, number of symbols, p_1,...,p_k, M_1,...,M_k}

	FILE *pfile = NULL;
	if ((pfile = fopen (MiPfile, "rb")) == (FILE *)0) {
		printf("\nUnable to open file.\n");
		return 1;
	} 

	double tmp;
	long ii,jj,kk,mm,rr,cc;

	// get number of symbols
	int k_sym; // number of symbols
	if (fread(&tmp, sizeof(double), 1, pfile) < 1){
		printf("\nUnable to read from opened file.\n");
		fclose(pfile);
		return 1;
	}
	k_sym = (int)round(tmp);

	// allocate memory for the di_vals shift amount array (d_1,...,d_k), and retrieve from file
	long *di_vals; // 1D array of integer shift amounts, (d_1,...,d_k)
	if (!(di_vals = new long[k_sym])) { // assign pointer to appropriate allocated spot in memory
		cout << "Error: out of memory." << endl;
		exit(1);
	}
	for (ii=0; ii<k_sym; ii++) {
		if (fread(&tmp, sizeof(double), 1, pfile) < 1){
			printf("\nUnable to read from opened file.\n");
			fclose(pfile);
			return 1;
		}
		di_vals[ii] = (long)round(tmp);
	}

	// get simplex height (=matrix size for binary sequences, k=2)
	long hk; // simplex height, in units of 1/(k-1); for binary symbol set, hk is just the matrix size
	if (fread(&tmp, sizeof(double), 1, pfile) < 1){
		printf("\nUnable to read from opened file.\n");
		fclose(pfile);
		return 1;
	}
	hk = (long)round(tmp);

	// get delta
	double delta; // matrix regularity parameter
	if (fread(&tmp, sizeof(double), 1, pfile) < 1){
		printf("\nUnable to read from opened file.\n");
		fclose(pfile);
		return 1;
	}
	delta = tmp;

	// get regularity type (pad only, pad and splice, or splice only)
	int p_and_s; // type of regularity: 0 for pad only; 1 for pad and splice; 2 for splice only
	if (fread(&tmp, sizeof(double), 1, pfile) < 1){
		printf("\nUnable to read from opened file.\n");
		fclose(pfile);
		return 1;
	}
	p_and_s = (int)round(tmp);

	// get number of block rows in the M_i*P set
	long nbrows; // number of block rows in the M_i*P set
	if (fread(&tmp, sizeof(double), 1, pfile) < 1){
		printf("\nUnable to read from opened file.\n");
		fclose(pfile);
		return 1;
	}
	nbrows = (long)round(tmp);

	// allocate memory for mxsz_vals and read in matrix block size values
	long *mxsz_vals; // 1D array holding the block matrix sizes (nbrows total)
	if (!(mxsz_vals = new long[nbrows])) { // assign pointer to appropriate allocated spot in memory
		cout << "Error: out of memory." << endl;
		exit(1);
	}
	long max_mxsz = 0;
	for (ii=0; ii<nbrows; ii++) {
		if (fread(&tmp, sizeof(double), 1, pfile) < 1){
			printf("\nUnable to read from opened file.\n");
			fclose(pfile);
			return 1;
		}
		mxsz_vals[ii] = (long)round(tmp);
		if (mxsz_vals[ii]>max_mxsz) {
			max_mxsz = mxsz_vals[ii];
		}
	}

	//allocate memory for MiP arrays, to hold matrices; there will be k matrices for the nth block row, each of the k matrices being of dimension mxsz_vals[n-1]
	double **MiP; // 2D array for the matrix entry data; given block row i, symbol #j, matrix row r, column c, format is: MiP[i*k_sym+j][r*mxsz_vals[i]+c]
	if (!(MiP = new double*[nbrows*k_sym])) {
		cout << "Error: out of memory." << endl;
		exit(1);
	}
	for (ii=0; ii<nbrows; ii++) {
		for (jj=0; jj<k_sym; jj++) {
			if (!(MiP[ii*k_sym + jj] = new double[mxsz_vals[ii]*mxsz_vals[ii]])) { // assign pointer to appropriate allocated spot in memory
				cout << "Error: out of memory." << endl;
				exit(1);
			}
		}
	}



	// ***
	// read in .bin file matrices (all block rows)
	// ***

	for (ii=0; ii<nbrows; ii++) {
		for (jj=0; jj<k_sym; jj++) {
			for (kk=0;kk<(mxsz_vals[ii]*mxsz_vals[ii]);kk++) {
	  			if (fread(&MiP[ii*k_sym+jj][kk], sizeof(double), 1, pfile) < 1) {
					printf("\nUnable to read from opened file.\n");
					return 1;
				}
			}
		}
	}

	fclose(pfile);



	//***
	// show header info
	//***

	cout << "Matrix data read:" << endl;
	cout << k_sym << " symbols; simplex height=" << hk << "; delta=" << delta;
	if (p_and_s==PAD_ONLY) {
		cout << "; pad only";
	} else if (p_and_s==PAD_AND_SPLICE) {
		cout << "; pad and splice";
	} else if (p_and_s==SPLICE_ONLY) {
		cout << "; splice only";
	} else {
		cout << "; p_and_s unknown";
	}
	cout << "; " << nbrows << " block rows." << endl;
	cout << "shift amounts: ";
	for (ii=0; ii<k_sym; ii++) {
		cout << di_vals[ii] << " ";
	}
	cout << endl;


	// DEBUG / for displaying / checking matrix loading
	/*
	ii = 0; // for which block row
	jj = 0; // for which symbol
	long mxsz_tmp = mxsz_vals[ii];
		for (rr=0;rr<mxsz_tmp;rr++) { // row number
			for (cc=0;cc<mxsz_tmp;cc++) { // column number
				cout << setprecision(8) << MiP[ii*k_sym+jj][mxsz_tmp*rr+cc] << " ";
			}
			cout << endl;
		}
		cout << endl;
	*/
	//

	
	//***
	// memory initializations for multiplication loop
	//***

	// create the array of vectors; a v_0 vector (array) will be assigned to each block row; addressing will be eg v_0[block row #, in [1,nbrows]][element number, in [1,mxsz_vals[block row #]]]
	double **v_0;
	if (!(v_0 = new double*[nbrows])) {
		cout << "Error: out of memory." << endl;
		exit(1);
	}
	for (ii=0; ii<nbrows; ii++) {
		if (!(v_0[ii] = new double[mxsz_vals[ii]])) { // assign pointer to appropriate allocated spot in memory
			cout << "Error: out of memory." << endl;
			exit(1);
		}
	}
	// temporary vector; it only needs to be as long as the maximum matrix size, over all the block rows
	double *v_tmp;
	if (!(v_tmp = new double[max_mxsz])) {
		cout << "Error: out of memory." << endl;
		exit(1);
	}
	// vector (array) to store the cumulative log-of-norms for each block row's v_0
	double *log_vec_norm;
	if (!(log_vec_norm = new double[nbrows])) {
		cout << "Error: out of memory." << endl;
		exit(1);
	}



	//***
	// matrix multiplication routine
	//***

	long r_1,c_1;
	double norm_v_0;
	int sym_val;
	bool bad_char = false;
	string symbol_string;
	vector<double> out_vec;


	ifstream string_file(symbol_string_file);	// try to create and open ifstream object "string_file"
	if (!string_file) {
		cout << "failed to open file " << symbol_string_file << "\n";
		exit(1);
	}

	// FETCH ANOTHER SYMBOL STRING LINE FROM THE FILE
	while (getline(string_file, symbol_string)) {
	

	// initialize entrywise norm accumulation vector
	long str_len = symbol_string.length();
	for (ii=0; ii<nbrows; ii++) {
		for (r_1 = 0; r_1<mxsz_vals[ii]; r_1++) {
			v_0[ii][r_1] = 1.0;
		}
		log_vec_norm[ii] = 0.0;
	}

	// loop through symbols in symbol_string (outer loop) and loop through block rows (inner loop)
	// this keeps track of the norm of M_iP...M_iPv_0, where v_0 is initially [1...1]'
	// given a sequence eg 12452, if we "convert" this to a matrix product, we're acting with 2 first, then 5, etc, so the symbol index counts "down" from the right-hand end of the array
	for (kk=str_len-1;kk>=0;kk--) {

		sym_val = symbol_string[kk]-'1';

		if (sym_val<0 || sym_val>=k_sym) {
			bad_char = true;
			break;
		}

		// perform the matrix multiplication (in the chain M_i P M_i P ... M_i P [1,...,1]') one block row at a time
		for (mm=0; mm<nbrows; mm++) {

			for (r_1=0; r_1<mxsz_vals[mm]; r_1++) {
				v_tmp[r_1] = 0.0;
				for (c_1=0; c_1<mxsz_vals[mm]; c_1++) {
					v_tmp[r_1] += MiP[mm*k_sym+sym_val][mxsz_vals[mm]*r_1+c_1]*v_0[mm][c_1];
				}
			}

			norm_v_0 = 0.0;
			for (r_1=0; r_1<mxsz_vals[mm]; r_1++) {
				norm_v_0 += v_tmp[r_1]*v_tmp[r_1];
			}
			norm_v_0 = sqrt(norm_v_0);
			for (r_1=0; r_1<mxsz_vals[mm]; r_1++) {
				v_0[mm][r_1] = v_tmp[r_1]/norm_v_0;
			}

			log_vec_norm[mm] += log(norm_v_0);

		} // end loop over block rows


	} // end loop over symbols in symbol string



	//***
	// sum over block rows and output results
	//***

	// at this point, v_0[mm] stores the normed result of the mmth block row of M_i P M_i P ... M_i P [1,...,1]', with log of the norm in log_vec_norm[mm] (could be used for product chaining, eg)

	double log_en_nrm;
	double max_br_nrm;

	if (bad_char) {
		cout << "symbol out of range in input string; multiplication failed" << endl;
	} else {
		double tmp_sum;
		for (mm=0; mm<nbrows; mm++) {
			tmp_sum = 0.0;
			for (r_1=0; r_1<mxsz_vals[mm]; r_1++) {
				tmp_sum += v_0[mm][r_1];
			}
			log_vec_norm[mm] += log(tmp_sum); // this now holds the log(||.||) value for block row mm
			if (mm==0) {
				max_br_nrm = log_vec_norm[0];
			} else {
				if (log_vec_norm[mm]>max_br_nrm) {
					max_br_nrm = log_vec_norm[mm];
				}
			}
		}
		tmp_sum = 0.0;
		for (mm=0; mm<nbrows; mm++) { // compute log(e^(log_vec_norm[0])+e^(log_vec_norm[1])+...)
			tmp_sum += exp(log_vec_norm[mm]-max_br_nrm);
		}
		log_en_nrm = max_br_nrm + log(tmp_sum);
	}

	out_vec.push_back(log_en_nrm);
	//printf("%.14f\n",log_en_nrm);


	// END LOOP OVER SYMBOL STRING FILE LINES
	}

	printf("%ld total file lines processed\n", out_vec.size());

	string file_out(symbol_string_file);
	string str_trm;
	str_trm = file_out.substr(file_out.size()-4);
	if (str_trm==".txt") {
		file_out = file_out.substr(0,file_out.size()-4)+"_output.txt";
	} else {
		file_out = file_out+"_output";
	}

	ofstream log_norm_out_file(file_out);
	log_norm_out_file.precision(15);
	copy(out_vec.begin(), out_vec.end(), ostream_iterator<double>(log_norm_out_file,"\n"));


	//***
	// release memory and exit
	//***

	delete [] di_vals;
	di_vals = NULL;
	delete [] mxsz_vals;
	mxsz_vals = NULL;
	for (ii=0; ii<k_sym; ii++) {
		delete [] MiP[ii];
	}
	delete [] MiP;
	MiP = NULL;
	delete [] v_tmp;
	v_tmp = NULL;
	delete [] log_vec_norm;
	log_vec_norm = NULL;
	for (ii=0; ii<nbrows; ii++) {
		delete [] v_0[ii];
	}
	delete [] v_0;
	v_0 = NULL;

	if (!bad_char) {
		return 0;
	} else {
		return 1;
	}

}



// for debugging
void cout_vec(double *vec, int s)
{
	int j;
	for (j=0; j<s; j++) {
		cout << vec[j] << " ";
	}
	cout << endl;
}







