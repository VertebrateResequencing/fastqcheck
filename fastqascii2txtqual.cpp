#include "fastqreader.h"
#include "fastqwriter.h"

#include <string>

using namespace std;

int main(int argc,char **argv) {

  if(argc < 3) {
    cout << "fastqascii2txtqual <fastq file> <fastq txt qual file>" << endl;
    cout << "Unpacks quality values" << endl;
    cout << "Optional [quality offset] argument, allows you to select the offset from the ascii value (i.e. ASCII-33 = quality value)" << endl;
    cout << "the default is 33 (Phred uses 33, Solexa use 64)." << endl;
    return 0;
  }

  FastqReader fastq_in(argv[1]); 
  FastqWriter fastq_out(argv[2]);

  int qualityoffset = -1;
  if(argc > 3) qualityoffset = atoi(argv[3]);

  if(qualityoffset != -1) {
    fastq_in.qualityoffset(qualityoffset);
    fastq_out.qualityoffset(qualityoffset);
  }

  fastq_in.open();
  fastq_out.open();

  cout << "Input file  : " << argv[1] << endl;
  cout << "Output file : " << argv[2] << endl;

  bool eof=false;
 
  fastq_out.set_textmode();

  for(;!eof;) {
    ScoredSequence s = fastq_in.next_line(eof);
    
    fastq_out.write(s);
  }

  fastq_in.close();
  fastq_out.close();
}
