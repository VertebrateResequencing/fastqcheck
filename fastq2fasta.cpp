#include "fastqreader.h"
#include "fastqwriter.h"

#include <string>

using namespace std;

int main(int argc,char **argv) {

  if(argc < 3) {
    cout << "fastqreader <fastq file> <fasta file>" << endl;
   return 0;
  }

  FastqReader fastq_in(argv[1]); 

  fastq_in.open();
  ofstream outfile(argv[2]);

  cout << "Input file  : " << argv[1] << endl;
  cout << "Output file : " << argv[2] << endl;

  bool eof=false;
  
  for(;!eof;) {
    ScoredSequence s = fastq_in.next_line(eof);
   
    if(!eof) {
      outfile << ">" << s.id << endl;

      outfile << s.sequence;

      outfile << endl;
    }
  }

  fastq_in.close();
  outfile.close();
}
