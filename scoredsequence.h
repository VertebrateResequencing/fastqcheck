#ifndef SCOREDSEQUENCE_H
#define SCOREDSEQUENCE_H

#include <string>
#include <vector>

using namespace std;

/// Represents a sequence with associated quality score
class ScoredSequence {
public:

  typedef vector<int> quality_type;
  typedef string      sequence_type;
  typedef string      id_type;

  id_type       id;
  sequence_type sequence;
  quality_type  quality;


  void set_quality(int new_quality) {
    quality.clear();
    quality.insert(quality.begin(),sequence.length(),new_quality);
  }
};

#endif
