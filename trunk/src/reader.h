#ifndef __READER_H__
#define __READER_H__

#include <limits> // for std::numeric_limits
#include <string> // for std::string
#include <iostream> // I/O
#include <fstream>// file I/O
#include <cassert>
#include <cstdlib>

//------------------------------------------------------------------------------
// Class to read input file
// Main job of this is to skip comments !!!
//------------------------------------------------------------------------------
class Reader
{
   public:
      Reader (char* file)
      {
         f.open (file);
         assert (f.is_open());
      };
      ~Reader()
      {
         f.close();
      };
      void operator>> (std::string &s);
      void operator>> (int &i);
      void operator>> (unsigned int &i);
      void operator>> (double &d);
      void getline (std::string &s);
      bool bos ();
      bool eos ();
      void entry (std::string);
      void begin_section (std::string);
      void end_section ();

   private:
      std::ifstream f;
      void skipComment ();

};

//------------------------------------------------------------------------------
// Jump over comments. These begin with // or are enclosed b/w /* ... */
//------------------------------------------------------------------------------
inline
void Reader::skipComment ()
{
   char c;
   const int c_MAX_LINESIZE=std::numeric_limits<int>::max();
   bool has_comment = false;
   
   if(!f) return;
   
   while(isspace(c=f.get())||c==','||c==';');
   
   if(c=='#'||c=='%'||(c=='/' && f.peek()=='/'))
   {
      //skip the rest of the line
      f.ignore(c_MAX_LINESIZE, '\n');
      has_comment = true;
   }
   else if (c=='/' && f.peek()=='*')
   {
      //skip everything in the comment block
      c=f.get();          //skip the first '*'
      char last='\0';
      while(!(f.eof())&&f.good())
      {
         c=f.get();
         if(c=='/'&&last=='*')break;
         else last=c;
      }
      has_comment = true;
   }
   else if(c!=EOF)
   {
      f.putback(c);
   }

   // If comment line was found, then there may be more comments
   // We skip until no comment lines are found
   if(has_comment)
      skipComment ();
}

//------------------------------------------------------------------------------
// Read a string
//------------------------------------------------------------------------------
inline
void Reader::operator>> (std::string &s)
{
   skipComment ();
   f >> s;
}

//------------------------------------------------------------------------------
// Read an integer
//------------------------------------------------------------------------------
inline
void Reader::operator>> (int &i)
{
   skipComment ();
   f >> i;
}

//------------------------------------------------------------------------------
// Read an unsigned integer
//------------------------------------------------------------------------------
inline
void Reader::operator>> (unsigned int &i)
{
   skipComment ();
   f >> i;
}

//------------------------------------------------------------------------------
// Read a double
//------------------------------------------------------------------------------
inline
void Reader::operator>> (double &d)
{
   skipComment ();
   f >> d;
}

//------------------------------------------------------------------------------
// Read whole line
//------------------------------------------------------------------------------
inline
void Reader::getline (std::string &s)
{
   skipComment ();
   std::getline(f, s);
}

//------------------------------------------------------------------------------
// Return true if end of section string "}" is encountered
// otherwise, put undo read and return false
//------------------------------------------------------------------------------
inline
bool Reader::eos ()
{
   std::streampos curr_pos = f.tellg ();

   std::string input;
   *this >> input;

   if(input != "}")
   {
      f.seekg (curr_pos);
      return false;
   }
   else
      return true;
}

//------------------------------------------------------------------------------
// Return true if beginning of section string "{" is encountered
// otherwise, put undo read and return false
//------------------------------------------------------------------------------
inline
bool Reader::bos ()
{
   std::streampos curr_pos = f.tellg ();

   std::string input;
   *this >> input;

   if(input != "{")
   {
      f.seekg (curr_pos);
      return false;
   }
   else
      return true;
}

//------------------------------------------------------------------------------
// Check that currently read string matches str1
//------------------------------------------------------------------------------
inline
void Reader::entry (std::string str1)
{
   std::string str2;

   *this >> str2;

   if(str1 != str2)
   {
      std::cout << "   Expecting " << str1 << " but found " << str2 << std::endl;
      abort ();
   }
}

//------------------------------------------------------------------------------
// Check that section "str" is reached
//------------------------------------------------------------------------------
inline
void Reader::begin_section (std::string str)
{
   this->entry (str);
   this->entry ("{");
}

//------------------------------------------------------------------------------
// Check that end of section is reached
//------------------------------------------------------------------------------
inline
void Reader::end_section ()
{
   this->entry ("}");
}

#endif
