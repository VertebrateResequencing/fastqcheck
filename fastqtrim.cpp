#include "fastqreader.h"
#include "fastqwriter.h"

#include <string>

using namespace std;

int main(int argc,char **argv) {

  if(argc < 4) {
    cout << "fastqtrim <fastq file> <trimmed file> <max length>" << endl;
   return 0;
  }

  FastqReader fastq_in(argv[1]); 
  FastqWriter fastq_out(argv[2]);

  fastq_in.open();
  fastq_out.open();

  int max_length  = atoi(argv[3]);
  cout << "Input file  : " << argv[1] << endl;
  cout << "Output file : " << argv[2] << endl;
  cout << "Max length  : " << max_length << endl;

  bool eof=false;
  
  for(;!eof;) {
    ScoredSequence s = fastq_in.next_line(eof);
   
    // If read is ok, trim it and write it out.
    if(!eof) {
      s.sequence = s.sequence.substr(0,max_length);
      s.quality.erase(s.quality.begin()+max_length,s.quality.end());

      fastq_out.write(s);
    }
  }

  fastq_in.close();
  fastq_out.close();
}
