#include "fastqreader.h"

#include <string>

using namespace std;


const int max_duplicates = 200;

int main(int argc,char **argv) {
  

  if(argc < 2) {
    cout << "fastq2duplicatecheck <fastq in file>" << endl;
    return 0;
  }

  FastqReader fastq_in(argv[1]); 

  fastq_in.open();

  cout << "Input file  : " << argv[1] << endl;

  bool eof=false;
 
  vector<string> sequences;
  for(;!eof;) {
    ScoredSequence s = fastq_in.next_line(eof);
    sequences.push_back(s.sequence);
  }

  fastq_in.close();

  sort(sequences.begin(),sequences.end());

  vector<int> occur_count(max_duplicates,0);
  int greater_than_max_count = 0;
  int greater_than_max_read_count = 0;

  int occurs=0;
  for(int i=0;i<(sequences.size()-1);i++) {
    if(sequences[i] == sequences[i+1]) {
      occurs++;
    } else {
      if(occurs < max_duplicates) 
        occur_count[occurs]++;
      else {
        greater_than_max_count++;
        greater_than_max_read_count += occurs+1;
      }
      occurs=0;
    }
  }

  if(sequences[sequences.size()-1] != sequences[sequences.size()-2]) {
    occur_count[0]++;
  }

  // Dump occurances
  cout << "Raw occurence (number of reads that occur a given number of times)" << endl;
  for(int i=0;i<occur_count.size();i++) {
    cout << i+1 << " " << occur_count[i] << endl;
  }
  cout << ">" << max_duplicates << " " << greater_than_max_count << endl;
  
  cout << "Cumulative occurences (number of reads)" << endl;
  int cumulative=greater_than_max_read_count;
  for(int i=occur_count.size()-1;i>=0;i--) {
    cumulative += occur_count[i]*(i+1);
    cout << ">=" << i+1 << " " << cumulative << " " << (static_cast<double>(cumulative)/static_cast<double>(sequences.size()))*100 <<  endl;
  }

}
