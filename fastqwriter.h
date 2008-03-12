#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "scoredsequence.h"

using namespace std;

class FastqWriter {
public:

  string filename;
  ofstream fastq_handle;
  int quality_conversion; // This value is subtracted for the ascii character value to obtain the quality score.
  bool textmode;

  FastqWriter(string filename_in) : filename(filename_in), textmode(false) {
    quality_conversion = 33;
  }

  void set_textmode() {
    textmode = true;
  }

  void qualityoffset(int qualityoffset_in) {
    quality_conversion = qualityoffset_in;
  }


  void write(const ScoredSequence &s) {
    fastq_handle << "@" << s.id << endl;
    fastq_handle << s.sequence << endl;

    fastq_handle << "+" << s.id << endl;
    if(!textmode) fastq_handle << quality_int_to_ascii(s.quality) << endl;
    else {
      for(ScoredSequence::quality_type::const_iterator i = s.quality.begin();i != s.quality.end();i++) {
        fastq_handle << *i << " ";
      }
      fastq_handle << endl;
    }
  }

  void open() {
    fastq_handle.open(filename.c_str());
  }

  void close() {
    fastq_handle.close();
  }

  string quality_int_to_ascii(const vector<int> qual_in) {
    string s;
    for(vector<int>::const_iterator i = qual_in.begin();i != qual_in.end();i++) {
      s.push_back((char) (*i)+quality_conversion);
    }

    return s;
  }

};
