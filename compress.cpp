# include <iostream>
# include <sys/timeb.h>
# include <cstring>
# include <vector>
# include <windows.h>
# include <stdio.h>
# include <stdlib.h>
# include <cmath>
# include "decompress.h"
# define _CRT_SECURE_NO_WARNINGS

using namespace std;

typedef struct
{
	int begin;
	int length;
} POSITION_RANGE;             //for N character and lowercase character

typedef struct
{
	int pos;
	int ch;
} POSITION_SPE_CHA;         //for special character

typedef struct {
	int pos;
	int length;
	string misStr;
}MatchEntry;                  //for first match result



const int MAX_SEQ_NUM = 2000;//maximum sequence number
const int MAX_CHA_NUM = 1 << 28;//maximum length of a chromosome
const int LINE_CHA_NUM =200;    //maximum character number of one line
const int kMerLen = 14; //the length of k-mer
const int kmer_bit_num = 2 * kMerLen; //bit numbers of k-mer
const int hashTableLen = 1 << kmer_bit_num; // length of hash table
const int VEC_SIZE = 1 <<20; //length for other character arrays
const int min_rep_len = 15;   //minimum replace length, matched string length exceeds min_rep_len, saved as matched information

string identifier;
int lineWidth, ref_code_len, seq_code_len, ref_low_len, seq_low_len, diff_low_len, nCha_len, spe_cha_len, seqNumber, seqBucketLen;
int percent, sec_seq_num; //the percentage and reference sequence number used for the second compress 
char *ref_code, *seq_code;
char *dismatched_str; //mismatched subsequence
int *refLoc; //reference hash location
int *refBucket; //reference hash bucket
int *seqBucket;        //sequence hash bucket
int *low_loc;  //lowercase tuple location

POSITION_RANGE *ref_low, *seq_low, *diff_low, *nCha;
POSITION_SPE_CHA *spe_cha;


vector<char *> seqName;
vector <string> identifier_vec;
vector <int> lineWidth_vec;
vector <POSITION_RANGE> ref_nCha, seq_nCha, lowCha; //N character vector and lowercase character vector;
vector <MatchEntry> matchResult;                   //store the first match result of one sequence
vector <MatchEntry> misMatchEntry;                //store the mismatched match entity of the second match
vector < vector <MatchEntry> > matchResult_vec;   //store the match result of all sequences of the first match
vector <int *> seqBucket_vec;                    //store the hash bucket of all sequences
vector < vector<int> > seqLoc_vec;               //store the collision elements of all sequences when creating hash table


inline void initial()   // allocate memory
{
	ref_code = new char[MAX_CHA_NUM];
	seq_code = new char[MAX_CHA_NUM];
	refBucket = new int[hashTableLen];
	refLoc = new int[MAX_CHA_NUM];
	ref_low = new POSITION_RANGE[VEC_SIZE];
	seq_low = new POSITION_RANGE[VEC_SIZE];
	diff_low = new POSITION_RANGE[VEC_SIZE];
	low_loc = new int[VEC_SIZE];
	nCha = new POSITION_RANGE[VEC_SIZE];
	spe_cha = new POSITION_SPE_CHA[VEC_SIZE];
	identifier_vec.reserve(seqNumber);
	lineWidth_vec.reserve(seqNumber);
	seqBucket_vec.reserve(seqNumber);
	seqLoc_vec.reserve(seqNumber);
}

// encoding ACGT - replaces each letter with specified number
int integerCoding(char ch) {
    if (ch == 'A') {
        return 0;
    } else if (ch == 'C') {
        return 1;
    } else if (ch == 'G') {
        return 2;
    } else if (ch == 'T') {
        return 3;
    } else {
        return 4;
    }
}


//read the input file, get the filename of every sequence and the sequence amount
int readFile(char *filename)
{
	FILE* fp = fopen(filename, "r");
	if (NULL == fp) {
		printf("Error: failed to open file %s\n", filename);
		exit(-1);
	}
	char *temp_name = new char[LINE_CHA_NUM];  //the length of filename
	while (fscanf(fp,"%s", temp_name) != EOF)
	{
		seqName.push_back(temp_name);
		temp_name = new char[LINE_CHA_NUM];
	}
	seqNumber = seqName.size();
	delete temp_name;
	return seqNumber;
}

//finds the first matches between the target sequence and the reference sequence using the k-mer hash table
void codeFirstMatch(char *tar_seq_code, int tar_seq_len, vector<MatchEntry> &matchResult)
{
	int pre_pos = 0;
	int step_len = tar_seq_len - kMerLen + 1;
	int max_length, max_k;
	int i, id, k, ref_idx, tar_idx, length, cur_pos, tar_value;
	string mismatched_str;
	mismatched_str.reserve(10240);

	MatchEntry me;
	matchResult.reserve(VEC_SIZE);

	for (int i = 0; i < step_len; i++) 
{
    unsigned int tar_value = 0;

    // Calculate the hash value of the k-mer starting at position i
    for (int k = kMerLen - 1; k >= 0; k--) 
    {
        tar_value <<= 2;
        tar_value += integerCoding(tar_seq_code[i + k]);
    }

    int id = refBucket[tar_value];
    if (id > -1) 
    {
        int max_length = -1;
        int max_k = -1;

        // Traverse the linked list to find the longest match
        for (int k = id; k != -1; k = refLoc[k]) 
        {
            int ref_idx = k + kMerLen;
            int tar_idx = i + kMerLen;
            int length = kMerLen;

            // Continue matching as long as characters in ref_code and tar_seq_code match
            while (ref_idx < ref_code_len && tar_idx < tar_seq_len && ref_code[ref_idx] == tar_seq_code[tar_idx]) 
            {
                ref_idx++;
                tar_idx++;
                length++;
            }

            // Update max_length and max_k if a longer match is found
            if (length >= min_rep_len && length > max_length)
            {
                max_length = length;
                max_k = k;
            }
        }

        // If a valid match is found, save the matched information
        if (max_length > -1) 
        {
            int cur_pos = max_k - pre_pos;  // Delta-coding for cur_pos
            me.pos = cur_pos;
            me.length = max_length - min_rep_len;
            me.misStr = mismatched_str;
            matchResult.push_back(me);

            i += max_length;  // Skip the length of the matched segment
            pre_pos = max_k + max_length;
            mismatched_str = "";
            
            if (i < tar_seq_len) 
                mismatched_str += '0' + integerCoding(tar_seq_code[i]);  // Store the integer code of mismatched nucleotides
            
            continue;
        }
    }

    // Append the integer code of the current nucleotide to mismatched_str
    mismatched_str += '0' + integerCoding(tar_seq_code[i]);
}

	if (i < tar_seq_len)
	{
		for (; i < tar_seq_len; i++)
			mismatched_str += '0' + integerCoding(tar_seq_code[i]);
		me.pos = 0;
		me.length = -min_rep_len;                //no match information, not 0 ,is -min_rep_len;
		me.misStr = mismatched_str;
		matchResult.push_back(me);
	}
}

//Reads a reference sequence file, extracts the sequence, converts lowercase characters to uppercase, and records their positions.
void referenceSequenceExtraction(char *str_referenceName)
{
	int _seq_code_len = 0, _ref_low_len = 1, letters_len = 0;//record lowercase from 1, diff_lowercase_loc[i]=0 means mismatching
	char temp_cha;
	bool flag = true;
	char cha[LINE_CHA_NUM];      //the content of one line

	printf("The reference file is %s.\n", str_referenceName);
	FILE* fp = fopen(str_referenceName, "r");
	if (NULL == fp)
	{
		printf("Error: fail to open reference file %s.\n", str_referenceName);
		exit(-1);
	}

	fgets(cha, LINE_CHA_NUM, fp);
	int c; // Use int to accommodate EOF

	while ((c = getc(fp)) != EOF)
	{
    	if (islower(c))
    	{
        	if (flag) //previous is upper case
        	{
            	flag = false; //change status of flag
            	ref_low[_ref_low_len].begin = letters_len;
            	letters_len = 0;
        	}
        	c = toupper(c);
    	}
    	else    //this case is upper case
    	{
        	if (!flag)  //previous is lower case
        	{
            	flag = true;
            	ref_low[_ref_low_len++].length = letters_len;
            	letters_len = 0;
        	}
    	}
    	if (c == 'A' || c == 'C' || c == 'G' || c == 'T')
        	ref_code[_seq_code_len++] = c;
    	letters_len++;
	}
	if (!flag)  //if flag=false, don't forget record the length
		ref_low[_ref_low_len++].length = letters_len;

	fclose(fp);
	ref_code_len = _seq_code_len;
    ref_low_len = _ref_low_len - 1;
}

//Similar to referenceSequenceExtraction, but for target sequences, including special characters and 'N' characters
void targetSequenceExtraction(char *str_sequenceName, int &_seq_code_len, int &_seq_low_len, int &_nCha_len, int &_spe_cha_len)
{
	FILE *fp = fopen(str_sequenceName, "r");
	if (NULL == fp) {
		printf("Error: fail to open sequence file %s.\n", str_sequenceName);
		return;
	}

	_seq_code_len = _seq_low_len = _nCha_len = _spe_cha_len = 0;
	int letters_len = 0, n_letters_len = 0;
	bool flag = true, n_flag = false;
	char cha[LINE_CHA_NUM];      //the content of one line
	char temp_cha;

	//get the identifier
	fgets(cha, LINE_CHA_NUM, fp);
	identifier = cha;
	identifier_vec.push_back(identifier);

	//get the lineWidth
	if (fscanf(fp, "%s", cha) != EOF)
		lineWidth = strlen(cha);
	lineWidth_vec.push_back(lineWidth);
	fseek(fp, -1L * (lineWidth+1), SEEK_CUR);//set the 'fp' to the beginning of the file again

	int c; // Use int to accommodate EOF

while ((c = fgetc(fp)) != EOF)
{
    if (islower(c))
    {
        if (flag) //previous is upper case
        {
            flag = false;
            seq_low[_seq_low_len].begin = letters_len;
            letters_len = 0;
        }
        c = toupper(c);
    }
    else if (isupper(c))
    {
        if (!flag)
        {
            flag = true;
            seq_low[_seq_low_len++].length = letters_len;
            letters_len = 0;
        }
    }

    letters_len++;

    if (c == 'A' || c == 'C' || c == 'G' || c == 'T')
    {
        seq_code[_seq_code_len++] = c;
    }
    else if (c != 'N')
    {
        spe_cha[_spe_cha_len].pos = _seq_code_len;
        spe_cha[_spe_cha_len++].ch = c - 'A';
    }

    if (!n_flag)
    {
        if (c == 'N')
        {
            nCha[_nCha_len].begin = n_letters_len;
            n_letters_len = 0;
            n_flag = true;
        }
    }
    else
    {
        if (c != 'N')
        {
            nCha[_nCha_len++].length = n_letters_len;
            n_letters_len = 0;
            n_flag = false;
        }
    }
    n_letters_len++;
}


	if (!flag)
		seq_low[_seq_low_len++].length = letters_len;

	if (n_flag)
		nCha[_nCha_len++].length = n_letters_len;

	for (int i = _spe_cha_len - 1; i > 0; i--)
		spe_cha[i].pos -= spe_cha[i-1].pos;

	fclose(fp);
}

// construction of k-mer hashing table and linked list
//this is used for efficient matching
void kMerHashingConstruct()//complete
{
	//initialize the point array
	for (int i = 0; i < hashTableLen; i++)
		refBucket[i] = -1;
	unsigned int value = 0;
	int step_len = ref_code_len - kMerLen + 1;

// Calculate the value of the first k-mer
for (int k = 0; k < kMerLen; k++) {
    value <<= 2;
    value += integerCoding(ref_code[k]);
}
refLoc[0] = refBucket[value];
refBucket[value] = 0;

int shift_bit_num = (kMerLen - 1) * 2;

// Calculate the value of the following k-mers using the last k-mer
for (int i = 1; i < step_len; i++) {
    value = (value & ((1 << shift_bit_num) - 1)) << 2; // Shift out the leftmost base
    value += integerCoding(ref_code[i + kMerLen - 1]);
    refLoc[i] = refBucket[value]; // refLoc[i] record the list of same values
    refBucket[value] = i;
}

}

//matches lowercase character information between the reference and target sequences
void seqLowercaseMatching(int _seq_low_len, int &_diff_low_len)
{
	int start_position = 1;
_diff_low_len = 0;

// Initialize the diff_low_loc array to zero
memset(low_loc, 0, sizeof(int) * _seq_low_len);

for (int i = 0; i < _seq_low_len; i++) {
    bool found = false;

    // Search from start_position to the end
    for (int j = start_position; j < ref_low_len; j++) {
        if (seq_low[i].begin == ref_low[j].begin && seq_low[i].length == ref_low[j].length) {
            low_loc[i] = j;
            start_position = j + 1;
            found = true;
            break;
        }
    }

    // If not found, search from start_position to the beginning
    if (!found) {
        for (int j = start_position - 1; j > 0; j--) {
            if (seq_low[i].begin == ref_low[j].begin && seq_low[i].length == ref_low[j].length) {
                low_loc[i] = j;
                start_position = j + 1;
                found = true;
                break;
            }
        }
    }

    // Record the mismatched information if still not found
    if (!found) {
        diff_low[_diff_low_len].begin = seq_low[i].begin;
        diff_low[_diff_low_len++].length = seq_low[i].length;
    }
}

}

//compresses data using run-length encoding and writes it to a file
void runLengthCoding(FILE *fp, int *vec, int length, int tolerance)
{
	vector<int> code;
	if (length > 0) {
        code.push_back(vec[0]);
        int count = 1;
        int i = 1;

        // Use a while loop to iterate through the vector
        while (i < length) {
            if (vec[i] - vec[i - 1] == tolerance) {
                count++;
            } else {
                code.push_back(count);
                code.push_back(vec[i]);
                count = 1;
            }
            i++;
        }

        // Push the final count
        code.push_back(count);
    }

    // Write the length of the code vector to the file
    int code_len = code.size();
    fprintf(fp, "%d ", code_len);

    // Write each element of the code vector to the file
    for (int i = 0; i < code_len; i++) {
        fprintf(fp, "%d ", code[i]);
    }
}

//writes a match entry to a file
void saveMatchEntry(FILE *fp, MatchEntry &_me)
{
	if (!_me.misStr.empty())
		fprintf(fp, "%s\n", _me.misStr.c_str());
	fprintf(fp, "%d %d\n", _me.pos, _me.length);
}

//saves the results of the first match to a file
void saveFirstMatchResult(FILE *fp, vector <MatchEntry> &_mr)
{
	for (unsigned int i = 0; i < _mr.size(); i++)
		saveMatchEntry(fp, _mr[i]);
    fprintf(fp, "\n");//The end flag of the first target sequence.
}

//saves identifier data using run-length encoding
void saveIdentifierData(FILE *fp, vector<string> &vec)
{
	vector<string> _vec;
	vector<int> code;

	_vec.push_back(vec[0]);
	string pre_str = vec[0];
	int cnt = 1;
	unsigned int i = 1;

    while (i < vec.size()) {
        if (vec[i] == pre_str) {
            cnt++;
        } else {
            code.push_back(cnt);
            _vec.push_back(vec[i]);
            pre_str = vec[i];
            cnt = 1;
        }
        i++;
    }
	code.push_back(cnt);
	int code_len = code.size();
	fprintf(fp, " %d", code_len);
	for (int i = 0; i < code_len; i++)
		fprintf(fp, " %d", code[i]);
	fprintf(fp, "\n");
	for (int i = 0; i < code_len; i++)
		fprintf(fp, "%s", _vec[i].c_str());
}

//writes position range data to a file
void savePositionRangeData(FILE *fp, int _vec_len, POSITION_RANGE *_vec)
{
	fprintf(fp, "%d ", _vec_len);
	for (int i = 0; i < _vec_len; i++)
		fprintf(fp, "%d %d ", _vec[i].begin, _vec[i].length);
}

//saves special character data to a file
void saveSpeChaData(FILE *filePointer, int length, POSITION_SPE_CHA *vector)
{
	// Initialize an array to mark special characters encountered
	int mark[26], temp;
	int count = 0;
	while (count < 26)
	{
		mark[count] = -1;
		count++;
	}
	// Create a vector to store unique special characters encountered
	std::vector<int> storage;
	count = 0;
	// Iterate through the given vector
	while (count < length)
	{
		// Write the position of the special character to the file
		fprintf(filePointer, "%d ", vector[count].pos);
		temp = vector[count].ch;
		// If the special character is encountered for the first time
		if (mark[temp] == -1)
		{
			// Store the special character in the vector and mark its index
			storage.push_back(temp);
			mark[temp] = storage.size() - 1;   //storage vector stores all types of special characters and mark array constructs the index
		}
		count++;
	}

	// Write the number of unique special characters to the file
	int size = storage.size();       
	fprintf(filePointer, "%d ", size);             //size is the type amount of special characters
	count = 0;
	// Write the unique special characters to the file
	while (count < size)
	{
		fprintf(filePointer, "%d ", storage[count]);
		count++;
	}

	// If there is more than one type of special character
	if (size != 1)
	{
		// Calculate the number of bits required to represent each special character
		unsigned int bit_num = ceil(log(size) / log(2));//the bit number of representing a special character
		// Calculate the number of characters that can be represented in 4 bytes
		unsigned int v_num = floor(32.0 / bit_num);     //the number of characters can be represented in 4 bytes
		count = 0;
		// Encode special characters and write to the file
		while (count < length)
		{
			unsigned int value = 0;
			unsigned int j = 0;
			while (j < v_num && count < length)
			{
				value <<= bit_num;
				value += mark[vector[count].ch];
				j++;
				count++;
			}
			fprintf(filePointer, "%u ", value);
		}
	}
}

//saves various sequence-related data to a file, including lowercase, 'N' characters, and special characters
void saveOtherData(FILE *fp, int _seq_low_len, int _nCha_len, int _spe_cha_len)
{
	int flag = 0;
	(_seq_low_len > 0 && ref_low_len > 0) ? (
		seqLowercaseMatching(_seq_low_len, diff_low_len),
		(2 * diff_low_len < _seq_low_len) ? (
			flag = 1,
			fprintf(fp, "%d ", flag),
			runLengthCoding(fp, low_loc, _seq_low_len, 1),
			savePositionRangeData(fp, diff_low_len, diff_low)
		) : (
			fprintf(fp, "%d ", flag),
			savePositionRangeData(fp, _seq_low_len, seq_low)
		)
	) : (
		fprintf(fp, "%d ", flag),
		savePositionRangeData(fp, _seq_low_len, seq_low)
	);

	savePositionRangeData(fp, _nCha_len, nCha);

	fprintf(fp, "%d ", _spe_cha_len);
	if (_spe_cha_len > 0)
	{
		saveSpeChaData(fp, _spe_cha_len, spe_cha);
	}
	fprintf(fp, "\n"); // The end of other data
}

//computes a hash value for a match entry for use in the second matching phase
int getHashValue(MatchEntry &_me)
{
	int result = 0;
	for (unsigned int i = 0; i < _me.misStr.size(); i++)
		result += _me.misStr[i] * 92083;
	result += _me.pos * 69061 + _me.length * 51787;
	return result % seqBucketLen;
}

//get the nearest prime larger than number as the length of hash bucket
int getNextPrime(const int startingNumber)
{
	int current = startingNumber + 1;
	bool isPrime = false;
	while (!isPrime)
	{
		isPrime = true;
		for (int divisor = 2; divisor < sqrt(startingNumber) + 1; divisor++)
		{
			if (current % divisor == 0)
			{
				isPrime = false;
				break;
			}
		}
		if (!isPrime) 
			current++;
	}
	return current;
}

//compares two match entries to see if they are the same
bool compareMatchEntry(MatchEntry &ref, MatchEntry &tar)//complete
{
	if (ref.pos == tar.pos && ref.length == tar.length && ref.misStr == tar.misStr)
		return true;
	else
		return false;
}

//computes the length of matching entries between reference and target sequences
int getMatchLength(vector <MatchEntry> &ref_me, unsigned int ref_idx, vector <MatchEntry> &tar_me, unsigned int tar_idx)
{
	int length = 0;
	while (ref_idx < ref_me.size() && tar_idx < tar_me.size() && compareMatchEntry(ref_me[ref_idx++], tar_me[tar_idx++]))
		length++;
	return length;

}


//constructs a hash table for the second matching phase using match entries
void matchResultHashConstruct(vector<MatchEntry> &_mr)
{
    int hashValue1, hashValue2, hashValue;
    vector<int> seqLoc;
    seqLoc.reserve(VEC_SIZE);
    int *seqBucket = new int[seqBucketLen];
    // Initialize sequence bucket array
    for (int i = 0; i < seqBucketLen; i++)
        seqBucket[i] = -1;

    // Construct hash table
    hashValue1 = getHashValue(_mr[0]);
    hashValue2 = (_mr.size() < 2) ? 0 : getHashValue(_mr[1]);
    hashValue = abs(hashValue1 + hashValue2) % seqBucketLen;
    seqLoc.push_back(seqBucket[hashValue]);
    seqBucket[hashValue] = 0;

    for (unsigned int i = 1; i < _mr.size() - 1; i++)
    {
        hashValue1 = hashValue2;
        hashValue2 = getHashValue(_mr[i + 1]);
        hashValue = abs(hashValue1 + hashValue2) % seqBucketLen;
        seqLoc.push_back(seqBucket[hashValue]);
        seqBucket[hashValue] = i;
    }

    seqLoc_vec.push_back(seqLoc);
    seqBucket_vec.push_back(seqBucket);
    seqLoc.clear();
}

//performs a second match between sequences using the constructed hash tables and writes the results to a file
void codeSecondMatch(FILE *fp, vector<MatchEntry> &_mr, int seqNum)
{
    int hashValue;
    int pre_seq_id = 1;
    int max_pos = 0, pre_pos = 0, delta_pos, length, max_length = 0, delta_length, seq_id = 0, delta_seq_id;
    int id, pos, secondMatchTotalLength = 0;
    unsigned int i;

    // Iterate through match entries
    for (i = 0; i < _mr.size() - 1; i++)
    {
        // Compute hash value based on match entries
        if (_mr.size() < 2)
            hashValue = abs(getHashValue(_mr[i])) % seqBucketLen;
        else
            hashValue = abs(getHashValue(_mr[i]) + getHashValue(_mr[i + 1])) % seqBucketLen;

        max_length = 0;
        // Check each sequence for potential matches
        for (int m = 0; m < min(seqNum - 1, sec_seq_num); m++)
        {
            id = seqBucket_vec[m][hashValue];
            // Traverse linked list to find the best match
            while (id != -1)
            {
                pos = id;
                length = getMatchLength(matchResult_vec[m], pos, _mr, i);
                if (length > 1 && length > max_length)
                {
                    seq_id = m + 1;
                    max_pos = pos;
                    max_length = length;
                }
                id = seqLoc_vec[m][pos];
            }
        }

        // If a match is found, process it
        if (max_length > 0)
        {
            // Delta encoding for sequence ID, position, and length
            delta_seq_id = seq_id - pre_seq_id;
            delta_length = max_length - 2;
            delta_pos = max_pos - pre_pos;

            // Update previous sequence ID and position
            pre_seq_id = seq_id;
            pre_pos = max_pos + max_length;
            secondMatchTotalLength += max_length;

            // Save mismatched entries first
            if (!misMatchEntry.empty())
            {
                for (unsigned int j = 0; j < misMatchEntry.size(); j++)
                    saveMatchEntry(fp, misMatchEntry[j]);
                misMatchEntry.clear();
            }

            // Save matched entry
            fprintf(fp, "%d %d %d\n", delta_seq_id, delta_pos, delta_length);
            i += max_length - 1;
        }
        else
        {
            misMatchEntry.push_back(_mr[i]);
        }
    }

    // Handle the last entry if not processed
    if (i == _mr.size() - 1)
    {
        misMatchEntry.push_back(_mr[i]);
    }

    // Save remaining mismatched entries
    if (!misMatchEntry.empty())
    {
        for (unsigned int j = 0; j < misMatchEntry.size(); j++)
            saveMatchEntry(fp, misMatchEntry[j]);
        misMatchEntry.clear();
    }

    // Mark the end of the target sequence
    fprintf(fp, "\n");
}

//performs the main sequence compression logic
void codeMatch(FILE *fp)
{
	vector<MatchEntry> mr;
sec_seq_num = ceil(percent * seqNumber / 100.0);  // Calculate the number of reference sequences for the second match
seqBucketLen = getNextPrime(VEC_SIZE);

for (int i = 1; i < seqNumber; i++)
{
    // Extract target sequence information
    targetSequenceExtraction(seqName[i], seq_code_len, seq_low_len, nCha_len, spe_cha_len);

    // Save additional data related to the sequence
    saveOtherData(fp, seq_low_len, nCha_len, spe_cha_len);

    // Perform the first match encoding
    codeFirstMatch(seq_code, seq_code_len, mr);

    if (i <= sec_seq_num)
    {
        // Store match results and construct hash for sequences within the secondary sequence number limit
        matchResult_vec.push_back(mr);
        matchResultHashConstruct(mr);
    }

    if (i == 1)
    {
        // Directly save the match results for the first sequence to be compressed
        saveFirstMatchResult(fp, mr);
    }
    else
    {
        // Perform second match encoding for subsequent sequences
        codeSecondMatch(fp, mr, i);
    }

    // Display compression status
    printf("Compressed sequence %s (%d/%d).\n", seqName[i], i, seqNumber - 1);

    // Clear match entries for the next iteration
    mr.clear();
}
}

//releases allocated memory
void compressClear()
{
    // Deallocate dynamically allocated arrays
    delete[] ref_code;
    delete[] seq_code;
    delete[] refLoc;
    delete[] refBucket;
    delete[] ref_low;
    delete[] seq_low;
    delete[] diff_low;
    delete[] low_loc;
    delete[] nCha;
    delete[] spe_cha;

    // Clear the sequence bucket vector
    for (auto it = seqBucket_vec.begin(); it != seqBucket_vec.end(); ++it)
    {
        if (*it != nullptr)
        {
            delete[] *it;
            *it = nullptr;
        }
    }
    seqBucket_vec.clear();

    // Clear the sequence name vector
    for (auto it = seqName.begin(); it != seqName.end(); ++it)
    {
        if (*it != nullptr)
        {
            delete[] *it;
            *it = nullptr;
        }
    }
    seqName.clear();
}

//extracts the filename without the path and extension
void extractFileName(char *srcFileName, char *destFileName)
{
    int len = strlen(srcFileName);
    int start = 0, end = len;
    int k = 0;
    bool dotFound = false;

    // Traverse the source filename from the end to determine start and end positions
    for (int i = len - 1; i >= 0; i--) 
    {
        char ch = srcFileName[i];

        if (ch == '.' && !dotFound)
        {
            end = i;
            dotFound = true;
        }
        else if (ch == '/')
        {
            start = i + 1;
            break;
        }
    }

    // Copy the relevant portion of the source filename to the destination
    for (int j = start; j < end; j++)
    {
        destFileName[k++] = srcFileName[j];
    }

    // Null-terminate the destination filename
    destFileName[k] = '\0';
}

//manages the entire compression process
void compress(char *filename)
{
	
	if (seqNumber > 1)
	{
		printf("Info: Compressing...Please wait for a moment.\n");
		initial();
		char temp_filename[100], resultFilename1[100], resultFilename2[100], cmd[100];
		extractFileName(filename, temp_filename);
		sprintf(resultFilename1, "%s.hrcm", temp_filename);//stores the match result of all to-be-compressed sequences
		sprintf(resultFilename2, "%s.desc", temp_filename);//stores the identifier and line width of all to-be-compressed sequences

		referenceSequenceExtraction(seqName[0]);//reference sequence information extraction
		kMerHashingConstruct();//construct the hash index for the first match based on k-mer

		FILE *fp1 = fopen(resultFilename1, "w");
		if (NULL == fp1) {
			printf("Error: fail to open %s.\n", resultFilename1);
			exit(-1);
		}
		codeMatch(fp1);
		fclose(fp1);

		FILE *fp2 = fopen(resultFilename2, "w");
		if (NULL == fp2) {
			printf("Error: fail to open %s.\n", resultFilename2);
			exit(-1);
		}
		runLengthCoding(fp2, &lineWidth_vec[0], seqNumber - 1, 0);//save lineWidth data
		saveIdentifierData(fp2, identifier_vec);//save identifier data
		fclose(fp2);
		sprintf(cmd, "./7za a %s.7z %s %s -m0=PPMd", temp_filename, resultFilename1, resultFilename2);
		sprintf(cmd, "rm -f %s %s", resultFilename1, resultFilename2);
		compressClear();
	}
	else
		printf("Error: There is not any to-be-compressed sequence, nothing to be done.\n");
	
}

void show_usage() {
	cout << "Usage: hrcm {compress | decompress}  -r {ref-file-path}{ [-t] {tar-file-path}|[-f] {filename} [percent]}\n";
	cout << "  {compress | decompress} is mode,  choose one of them according to requirement, required\n";
	cout << "  -r is the reference, the {ref-file-path} followed, required\n";
	cout << "  -t is the target, a single to-be-compressed file path {tar-file-path} followed, optional\n";
	cout << "  -f is the alternative option of -t, a set of to-be-compressed file paths included in {filename}, optional\n";
    cout << "  [percent] is the percentage of the second-level matching, default is 10, means 10% of sequences will be used for the second-level matching, optional when -f, illegal when -t\n";
	cout << "Examples:\n";
	cout << "  hrcm compress -r hg17_chr22.fa -t hg18_chr22.fa\n";
	cout << "  hrcm decompress -r hg17_chr22.fa -t hg18_chr22.7z\n";
	cout << "  hrcm compress -r hg17_chr22.fa -f filename.txt 20\n";
	cout << "  hrcm decompress -r hg17_chr22.fa -f filename.txt 20\n";
}

int main(int argc, char *argv[]) {
	bool arg_flag = true, ref_flag = false, tar_flag = false, tar_set_flag = false, compressed = false, decompressed = false;
	char *mode = NULL, *ref_file = NULL, *tar_file = NULL, *per = NULL;
	int oc;

	struct  timeval  start;
	struct  timeval  end;
	unsigned long timer;
	//gettimeofday(&start,NULL);//first argument get the time result, second argument get the timezone

	if (argc < 6 || argc > 7)
	{
		show_usage();
		return 0;
	}
	mode = argv[1];                         //mode = compress or decompress
	percent = 10; //default is 10
	if (argc == 7)
	{
		per = argv[6];
		percent = atoi(per); 
	}

//manually parsing command line arguments using loop
for (int i = 1; i < argc; ++i) 
{
    if (argv[i][0] == '-') 
    {
        switch (argv[i][1]) 
        {
            case 'r':
                if (i + 1 < argc) 
                {
                    ref_file = argv[++i];
                    ref_flag = true;
                }
                else 
                {
                    arg_flag = false;
                }
                break;
            case 't':
                if (i + 1 < argc) 
                {
                    tar_file = argv[++i];
                    tar_flag = true;
                }
                else 
                {
                    arg_flag = false;
                }
                break;
            case 'f':
                if (i + 1 < argc) 
                {
                    tar_file = argv[++i];
                    tar_set_flag = true;
                }
                else 
                {
                    arg_flag = false;
                }
                break;
            default:
                arg_flag = false;
                break;
        }
    }
    else 
    {
        arg_flag = false;
    }
}

	if (!(tar_flag ^ tar_set_flag))
		arg_flag = false;
	if (arg_flag && ref_flag)
	{
		char *temp_name = new char[LINE_CHA_NUM];
		strcpy(temp_name, ref_file);
		seqName.push_back(temp_name);
		if (tar_flag)
		{
			temp_name = new char[LINE_CHA_NUM];
			strcpy(temp_name, tar_file);
			seqName.push_back(temp_name);
			seqNumber = 2;
		}
		else
			readFile(tar_file);		
		if (strcmp(mode, "compress") == 0)
		{
			compress(tar_file);
			compressed = true;
		}
		if (strcmp(mode, "decompress") == 0)
		{
			decompress(tar_file);
			decompressed = true;
		}

	}
	if (!(compressed||decompressed))
	{
		show_usage();
		return 0;
	}
	//gettimeofday(&end,NULL);
	timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
	if (compressed)
		printf("Total compression time = %lf ms; %lf s\n", timer/1000.0, timer/1000.0/1000.0);
	else
		printf("Total decompression time = %lf ms; %lf s\n", timer/1000.0, timer/1000.0/1000.0);
}

