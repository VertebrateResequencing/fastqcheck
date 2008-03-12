#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "scoredsequence.h"

using namespace std;

class FastqReader {
public:

  string filename;
  ifstream fastq_handle;
  int quality_conversion; // This value is subtracted for the ascii character value to obtain the quality score.

  FastqReader(string filename_in) : filename(filename_in) {
    quality_conversion = 33;
  }

  void qualityoffset(int qualityoffset_in) {
    quality_conversion = qualityoffset_in;
  }

  vector<ScoredSequence> parse_all() {
    vector<ScoredSequence> v;
    
    open();

    bool eof = false;
    for(;!eof;) {
      ScoredSequence s = next_line(eof);
      if(eof == false) v.push_back(s);
    }

    close();

    return v;
  }

  void open() {
    fastq_handle.open(filename.c_str());
  }

  // reads a line in to a ScoredSequence
  ScoredSequence next_line(bool &eof) {
    string line;
    ScoredSequence s;
    
    // ID line (strip leading character)
    getline(fastq_handle,line);
    if(fastq_handle.eof()) {eof = true; return s;}

    s.id = line.substr(1);
    
    // Sequence line
    getline(fastq_handle,line);
    if(fastq_handle.eof()) {eof = true; return s;}
    s.sequence = line;
    
    // Another ID line, which I will ignore
    getline(fastq_handle,line);
    if(fastq_handle.eof()) {eof = true; return s;}

    // Quality line
    getline(fastq_handle,line);
    if(fastq_handle.eof()) {eof = true; return s;}
   
    for(int n=0;n<line.length();n++) {
      s.quality.push_back((int) line[n]-quality_conversion);
    }
  }

  void close() {
    fastq_handle.close();
  }

};
