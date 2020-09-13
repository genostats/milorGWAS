#include <string>

#ifndef STRINGSTREAM_LITE
#define STRINGSTREAM_LITE

// to stream strings delimited by single characters.
// the string is modified in-place, make a copy if needed !
// !! don't deallocate the string !!
class stringstream_lite {
  public:

  char * a;
  char delim;
  char * token;
  bool eof;
  int token_length;

  stringstream_lite(char * a_, char delim_) : a(a_), delim(delim_), eof(*a == 0) {}

  stringstream_lite(std::string & s, char delim_) : a(&s[0]), delim(delim_), eof(*a == 0) {}

  // insère un 0 à la fin du token, renvoie sa longueur,
  // positionne a au début du token suivant (sauf si fin de chaine)
  int next_token() {
    token = a;
    if(*a == 0) {
      eof = true;
      return 0; 
    }
    while(*a != delim && *a != 0) a++;
    token_length = (a-token);
    if(*a == delim) {
      *a = 0;
      a++;
    } 
    return token_length;
  }

  stringstream_lite & operator>>(std::string & s) {
    token_length = next_token();
    if(token_length > 0) {
      s.assign(token);
    }
    return *this;
  }

  stringstream_lite & operator>>(int & x) {
    token_length = next_token();
    if(token_length > 0) {
      x = atoi(token);
    }
    return *this;
  } 

  stringstream_lite & operator>>(double & x) {
    token_length = next_token();
    if(token_length > 0) {
      x = atof(token);
    }
    return *this;
  } 

  operator bool() { 
    return !eof;
  }
};

#endif
