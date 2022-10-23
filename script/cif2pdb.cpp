/* Compile this program by: 
 * $ g++ -O3 cif2pdb.cpp -o cif2pdb
 */

const char* docstring=""
"cif2pdb input.cif prefix\n"
"    convert PDBx/mmCIF format input file 'input.cif' to PDB format output prefix*\n"
"\n"
"Ouput:\n"
"    prefix[chainID].pdb        - protein and nucleic acid\n"
"    prefix_res_chainID_idx.pdb - ligand\n"
"    prefix.txt                 - citation title\n"
"                                 pubmed ID\n"
"                                 fasta sequence\n"
"                                 ligand binding residues\n"
"                                 [1] ligand filename\n"
"                                 [2] receptor chainID\n"
"                                 [3] idx\n"
"                                 [4] original index of residue in contact\n"
"                                 [5] renumber index of residue in contact\n"
;

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <cstdlib>
using namespace std;

/* PDBParser */
struct AtomUnit
{
    string name; // atom name
    vector<double> xyz;
    double bfactor;
    string element;
};

struct ResidueUnit
{
    int resi;
    char icode;
    int het; // 1 - ATOM; 2 - modified residue; 3 - ligand
    string resn;
    vector<AtomUnit> atoms;
};

struct ChainUnit
{
    string asym_id;
    int mol_type;
    vector<ResidueUnit> residues;
};

struct ModelUnit
{
    vector<ChainUnit> chains;
};

void deepClean(AtomUnit &atom)
{
    atom.name.clear();
    atom.xyz.clear();
    atom.element.clear();
}

void deepClean(ResidueUnit &residue)
{
    size_t a;
    for (a=0;a<residue.atoms.size();a++) deepClean(residue.atoms[a]);
    residue.atoms.clear();
    residue.resn.clear();
}

void deepClean(ChainUnit &chain)
{
    size_t r;
    for (r=0;r<chain.residues.size();r++) deepClean(chain.residues[r]);
    chain.residues.clear();
    chain.asym_id.clear();
}

void deepClean(ModelUnit &pep)
{
    size_t c;
    for (c=0;c<pep.chains.size();c++) deepClean(pep.chains[c]);
    pep.chains.clear();
}

/* StringTools START */
string Upper(const string &inputString)
{
    string result=inputString;
    transform(result.begin(), result.end(), result.begin(), ::toupper);
    return result;
}

string Lower(const string &inputString)
{
    string result=inputString;
    transform(result.begin(), result.end(), result.begin(), ::tolower);
    return result;
}

string Join(const string sep, const vector<string>& string_vec,
    const int joinFrom=0)
{
    if (string_vec.size()<=joinFrom) return "";
    string joined_str=string_vec[joinFrom];
    for (int s=joinFrom+1;s<string_vec.size();s++)
        joined_str+=sep+string_vec[s];
    return joined_str;
}

/* split a long string into vectors by whitespace 
 * line          - input string
 * line_vec      - output vector 
 * delimiter     - delimiter */
void Split(const string &line, vector<string> &line_vec,
    const char delimiter=' ',const bool ignore_quotation=false)
{
    bool within_word = false;
    bool within_quotation = false;
    for (size_t pos=0;pos<line.size();pos++)
    {
        if (ignore_quotation==false && (line[pos]=='"' || line[pos]=='\''))
        {
            if (within_quotation) within_quotation=false;
            else within_quotation=true;
        }
        else if (line[pos]=='\n' || line[pos]=='\r')
        {
            within_quotation=false;
        }
        if (line[pos]==delimiter && within_quotation==false)
        {
            within_word = false;
            continue;
        }
        if (!within_word)
        {
            within_word = true;
            line_vec.push_back("");
        }
        line_vec.back()+=line[pos];
    }
}

void clear_line_vec(vector<string> &line_vec)
{
    int i;
    for (i=0;i<line_vec.size();i++) line_vec[i].clear();
    line_vec.clear();
}

string Basename(const string &inputString)
{
    string result="";
    int i;
    vector<string> line_vec;
    Split(inputString,line_vec,'/');
    if (line_vec.size())
    {
        result=line_vec.back();
        clear_line_vec(line_vec);
        Split(result,line_vec,'\\');
        result=line_vec.back();
        clear_line_vec(line_vec);
    }
    return result;
}

string Trim(const string &inputString,const string &char_list=" \n\r\t")
{
    string result = inputString;
    int idxBegin = inputString.find_first_not_of(char_list);
    int idxEnd = inputString.find_last_not_of(char_list);
    if (idxBegin >= 0 && idxEnd >= 0)
        result = inputString.substr(idxBegin, idxEnd + 1 - idxBegin);
    else result = "";
    return result;
}

string lstrip(const string &inputString,const string &char_list=" \n\r\t")
{
    string result = inputString;
    int idxBegin = inputString.find_first_not_of(char_list);
    if (idxBegin >= 0) result = inputString.substr(idxBegin);
    else result = "";
    return result;
}

string rstrip(const string &inputString,const string &char_list=" \n\r\t")
{
    string result=inputString;
    int idxEnd = inputString.find_last_not_of(char_list);
    if (idxEnd >= 0) result = inputString.substr(0, idxEnd + 1);
    else result = "";
    return result;
}

inline bool StartsWith(const string &longString, const string &shortString)
{
    return (longString.size()>=shortString.size() &&
            longString.substr(0,shortString.size())==shortString);
}

inline bool EndsWith(const string &longString, const string &shortString)
{
    return (longString.size()>=shortString.size() &&
            longString.substr(longString.size()-shortString.size(),
                shortString.size())==shortString);
}

inline string formatString(const string &inputString,const int width=8, 
    const int digit=3)
{
    string result=Trim(inputString," ");
    if (StartsWith(result,"00")) result='0'+lstrip(result,"0");
    size_t found=result.find_first_of('.');
    int i;
    if (found==string::npos)
    {
        result+='.';
        found=result.find_first_of('.');
    }
    int curWidth=result.size();
    if (curWidth<found+digit+1)
        for (i=0;i<((found+digit+1)-curWidth);i++) result+='0';
    else if (curWidth>found+digit+1)
    {
        long int extra_prod=1;
        for (i=0;i+found+1<curWidth;i++) extra_prod*=10;
        long int first=atoi((result.substr(0,found)).c_str())*extra_prod;
        long int second=atoi((lstrip(result.substr(found+1),"0")).c_str());
        if (result[0]=='-') second=-second;
        stringstream buf;
        buf<<fixed<<setprecision(digit)<<(first+second+.5)/extra_prod;
        result=buf.str();
        buf.str(string());
    }
    // -0.000
    if (StartsWith(result,"-0."))
    {
        bool allzero=true;
        for (i=found+1;i<result.size();i++)
        {
            if (result[i]!='0')
            {
                allzero=false;
                break;
            }
        }
        if (allzero) result=result.substr(1);
    }
    if (width)
    {
        curWidth=result.size();
        if (curWidth>width)
        {
            result=result.substr(0,width);
            //result=result.substr(result.size()-width);
        }
        else if (curWidth<width)
            for (i=0;i<width-curWidth;i++) result=' '+result;
    }
    return result;
}

inline double ReadDouble(const string &inputString,const int digit=3)
{
    size_t found=inputString.find_first_of('.');
    double value=atof(inputString.c_str());
    if (found!=string::npos)
    {
        string right=inputString.substr(found+1,digit);
        while (right.size()<digit) right+='0';
        int decimal=atoi(right.c_str());
        int prod=1;
        int i;
        for (i=0;i<digit;i++) prod*=10;
        if (value>0) value= double(decimal)/prod;
        else         value=-double(decimal)/prod;
        value+=atol(inputString.substr(0,found).c_str());
    }
    else value=atol(inputString.c_str());
    return value;
}

inline string writeDouble(const double value, 
    const int width=8, const int digit=3)
{
    stringstream buf;
    buf<<setw(width)<<setiosflags(ios::fixed)<<setprecision(digit)<<value;
    string result=buf.str().substr(0,width);
    buf.str(string());
    return result;
}

/* StringTools END */
/* pstream START */

// PStreams - POSIX Process I/O for C++

//        Copyright (C) 2001 - 2017 Jonathan Wakely
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)
//

/* do not compile on windows, which does not have cygwin */
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) && !defined(__CYGWIN__)
#define NO_PSTREAM
#else

#ifndef REDI_PSTREAM_H_SEEN
#define REDI_PSTREAM_H_SEEN

#include <ios>
#include <streambuf>
#include <istream>
#include <ostream>
#include <string>
#include <vector>
#include <algorithm>    // for min()
#include <cerrno>       // for errno
#include <cstddef>      // for size_t, NULL
#include <cstdlib>      // for exit()
#include <sys/types.h>  // for pid_t
#include <sys/wait.h>   // for waitpid()
#include <sys/ioctl.h>  // for ioctl() and FIONREAD
#if defined(__sun)
# include <sys/filio.h> // for FIONREAD on Solaris 2.5
#endif
#include <unistd.h>     // for pipe() fork() exec() and filedes functions
#include <signal.h>     // for kill()
#include <fcntl.h>      // for fcntl()
#if REDI_EVISCERATE_PSTREAMS
# include <stdio.h>     // for FILE, fdopen()
#endif


/// The library version.
#define PSTREAMS_VERSION 0x0101   // 1.0.1

/**
 *  @namespace redi
 *  @brief  All PStreams classes are declared in namespace redi.
 *
 *  Like the standard iostreams, PStreams is a set of class templates,
 *  taking a character type and traits type. As with the standard streams
 *  they are most likely to be used with @c char and the default
 *  traits type, so typedefs for this most common case are provided.
 *
 *  The @c pstream_common class template is not intended to be used directly,
 *  it is used internally to provide the common functionality for the
 *  other stream classes.
 */
namespace redi
{
  /// Common base class providing constants and typenames.
  struct pstreams
  {
    /// Type used to specify how to connect to the process.
    typedef std::ios_base::openmode           pmode;

    /// Type used to hold the arguments for a command.
    typedef std::vector<std::string>          argv_type;

    /// Type used for file descriptors.
    typedef int                               fd_type;

    static const pmode pstdin  = std::ios_base::out; ///< Write to stdin
    static const pmode pstdout = std::ios_base::in;  ///< Read from stdout
    static const pmode pstderr = std::ios_base::app; ///< Read from stderr

    /// Create a new process group for the child process.
    static const pmode newpg   = std::ios_base::trunc;

  protected:
    enum { bufsz = 32 };  ///< Size of pstreambuf buffers.
    enum { pbsz  = 2 };   ///< Number of putback characters kept.
  };

  /// Class template for stream buffer.
  template <typename CharT, typename Traits = std::char_traits<CharT> >
    class basic_pstreambuf
    : public std::basic_streambuf<CharT, Traits>
    , public pstreams
    {
    public:
      // Type definitions for dependent types
      typedef CharT                             char_type;
      typedef Traits                            traits_type;
      typedef typename traits_type::int_type    int_type;
      typedef typename traits_type::off_type    off_type;
      typedef typename traits_type::pos_type    pos_type;
      /** @deprecated use pstreams::fd_type instead. */
      typedef fd_type                           fd_t;

      /// Default constructor.
      basic_pstreambuf();

      /// Constructor that initialises the buffer with @a cmd.
      basic_pstreambuf(const std::string& cmd, pmode mode);

      /// Constructor that initialises the buffer with @a file and @a argv.
      basic_pstreambuf( const std::string& file,
                        const argv_type& argv,
                        pmode mode );

      /// Destructor.
      ~basic_pstreambuf();

      /// Initialise the stream buffer with @a cmd.
      basic_pstreambuf*
      open(const std::string& cmd, pmode mode);

      /// Initialise the stream buffer with @a file and @a argv.
      basic_pstreambuf*
      open(const std::string& file, const argv_type& argv, pmode mode);

      /// Close the stream buffer and wait for the process to exit.
      basic_pstreambuf*
      close();

      /// Send a signal to the process.
      basic_pstreambuf*
      kill(int signal = SIGTERM);

      /// Send a signal to the process' process group.
      basic_pstreambuf*
      killpg(int signal = SIGTERM);

      /// Close the pipe connected to the process' stdin.
      void
      peof();

      /// Change active input source.
      bool
      read_err(bool readerr = true);

      /// Report whether the stream buffer has been initialised.
      bool
      is_open() const;

      /// Report whether the process has exited.
      bool
      exited();

#if REDI_EVISCERATE_PSTREAMS
      /// Obtain FILE pointers for each of the process' standard streams.
      std::size_t
      fopen(FILE*& in, FILE*& out, FILE*& err);
#endif

      /// Return the exit status of the process.
      int
      status() const;

      /// Return the error number (errno) for the most recent failed operation.
      int
      error() const;

    protected:
      /// Transfer characters to the pipe when character buffer overflows.
      int_type
      overflow(int_type c);

      /// Transfer characters from the pipe when the character buffer is empty.
      int_type
      underflow();

      /// Make a character available to be returned by the next extraction.
      int_type
      pbackfail(int_type c = traits_type::eof());

      /// Write any buffered characters to the stream.
      int
      sync();

      /// Insert multiple characters into the pipe.
      std::streamsize
      xsputn(const char_type* s, std::streamsize n);

      /// Insert a sequence of characters into the pipe.
      std::streamsize
      write(const char_type* s, std::streamsize n);

      /// Extract a sequence of characters from the pipe.
      std::streamsize
      read(char_type* s, std::streamsize n);

      /// Report how many characters can be read from active input without blocking.
      std::streamsize
      showmanyc();

    protected:
      /// Enumerated type to indicate whether stdout or stderr is to be read.
      enum buf_read_src { rsrc_out = 0, rsrc_err = 1 };

      /// Initialise pipes and fork process.
      pid_t
      fork(pmode mode);

      /// Wait for the child process to exit.
      int
      wait(bool nohang = false);

      /// Return the file descriptor for the output pipe.
      fd_type&
      wpipe();

      /// Return the file descriptor for the active input pipe.
      fd_type&
      rpipe();

      /// Return the file descriptor for the specified input pipe.
      fd_type&
      rpipe(buf_read_src which);

      void
      create_buffers(pmode mode);

      void
      destroy_buffers(pmode mode);

      /// Writes buffered characters to the process' stdin pipe.
      bool
      empty_buffer();

      bool
      fill_buffer(bool non_blocking = false);

      /// Return the active input buffer.
      char_type*
      rbuffer();

      buf_read_src
      switch_read_buffer(buf_read_src);

    private:
      basic_pstreambuf(const basic_pstreambuf&);
      basic_pstreambuf& operator=(const basic_pstreambuf&);

      void
      init_rbuffers();

      pid_t         ppid_;        // pid of process
      fd_type       wpipe_;       // pipe used to write to process' stdin
      fd_type       rpipe_[2];    // two pipes to read from, stdout and stderr
      char_type*    wbuffer_;
      char_type*    rbuffer_[2];
      char_type*    rbufstate_[3];
      /// Index into rpipe_[] to indicate active source for read operations.
      buf_read_src  rsrc_;
      int           status_;      // hold exit status of child process
      int           error_;       // hold errno if fork() or exec() fails
    };

  /// Class template for common base class.
  template <typename CharT, typename Traits = std::char_traits<CharT> >
    class pstream_common
    : virtual public std::basic_ios<CharT, Traits>
    , virtual public pstreams
    {
    protected:
      typedef basic_pstreambuf<CharT, Traits>       streambuf_type;

      typedef pstreams::pmode                       pmode;
      typedef pstreams::argv_type                   argv_type;

      /// Default constructor.
      pstream_common();

      /// Constructor that initialises the stream by starting a process.
      pstream_common(const std::string& cmd, pmode mode);

      /// Constructor that initialises the stream by starting a process.
      pstream_common(const std::string& file, const argv_type& argv, pmode mode);

      /// Pure virtual destructor.
      virtual
      ~pstream_common() = 0;

      /// Start a process.
      void
      do_open(const std::string& cmd, pmode mode);

      /// Start a process.
      void
      do_open(const std::string& file, const argv_type& argv, pmode mode);

    public:
      /// Close the pipe.
      void
      close();

      /// Report whether the stream's buffer has been initialised.
      bool
      is_open() const;

      /// Return the command used to initialise the stream.
      const std::string&
      command() const;

      /// Return a pointer to the stream buffer.
      streambuf_type*
      rdbuf() const;

#if REDI_EVISCERATE_PSTREAMS
      /// Obtain FILE pointers for each of the process' standard streams.
      std::size_t
      fopen(FILE*& in, FILE*& out, FILE*& err);
#endif

    protected:
      std::string       command_; ///< The command used to start the process.
      streambuf_type    buf_;     ///< The stream buffer.
    };


  /**
   * @class basic_ipstream
   * @brief Class template for Input PStreams.
   *
   * Reading from an ipstream reads the command's standard output and/or
   * standard error (depending on how the ipstream is opened)
   * and the command's standard input is the same as that of the process
   * that created the object, unless altered by the command itself.
   */

  template <typename CharT, typename Traits = std::char_traits<CharT> >
    class basic_ipstream
    : public std::basic_istream<CharT, Traits>
    , public pstream_common<CharT, Traits>
    , virtual public pstreams
    {
      typedef std::basic_istream<CharT, Traits>     istream_type;
      typedef pstream_common<CharT, Traits>         pbase_type;

      using pbase_type::buf_;  // declare name in this scope

      // Ensure a basic_ipstream will read from at least one pipe
      pmode readable(pmode mode)
      {
        if (!(mode & (pstdout|pstderr)))
          mode |= pstdout;
        return mode;
      }

    public:
      /// Type used to specify how to connect to the process.
      typedef typename pbase_type::pmode            pmode;

      /// Type used to hold the arguments for a command.
      typedef typename pbase_type::argv_type        argv_type;

      /// Default constructor, creates an uninitialised stream.
      basic_ipstream()
      : istream_type(NULL), pbase_type()
      { }

      /**
       * @brief Constructor that initialises the stream by starting a process.
       *
       * Initialises the stream buffer by calling do_open() with the supplied
       * arguments.
       *
       * @param cmd   a string containing a shell command.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, pmode)
       */
      explicit
      basic_ipstream(const std::string& cmd, pmode mode = pstdout)
      : istream_type(NULL), pbase_type(cmd, readable(mode))
      { }

      /**
       * @brief Constructor that initialises the stream by starting a process.
       *
       * Initialises the stream buffer by calling do_open() with the supplied
       * arguments.
       *
       * @param file  a string containing the pathname of a program to execute.
       * @param argv  a vector of argument strings passed to the new program.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, const argv_type&, pmode)
       */
      basic_ipstream( const std::string& file,
                      const argv_type& argv,
                      pmode mode = pstdout )
      : istream_type(NULL), pbase_type(file, argv, readable(mode))
      { }

      /**
       * @brief Constructor that initialises the stream by starting a process.
       *
       * Initialises the stream buffer by calling
       * @c do_open(argv[0],argv,mode|pstdout)
       *
       * @param argv  a vector of argument strings passed to the new program.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, const argv_type&, pmode)
       */
      explicit
      basic_ipstream(const argv_type& argv, pmode mode = pstdout)
      : istream_type(NULL), pbase_type(argv.at(0), argv, readable(mode))
      { }

#if __cplusplus >= 201103L
      template<typename T>
        explicit
        basic_ipstream(std::initializer_list<T> args, pmode mode = pstdout)
        : basic_ipstream(argv_type(args.begin(), args.end()), mode)
        { }
#endif

      /**
       * @brief Destructor.
       *
       * Closes the stream and waits for the child to exit.
       */
      ~basic_ipstream()
      { }

      /**
       * @brief Start a process.
       *
       * Calls do_open( @a cmd , @a mode|pstdout ).
       *
       * @param cmd   a string containing a shell command.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, pmode)
       */
      void
      open(const std::string& cmd, pmode mode = pstdout)
      {
        this->do_open(cmd, readable(mode));
      }

      /**
       * @brief Start a process.
       *
       * Calls do_open( @a file , @a argv , @a mode|pstdout ).
       *
       * @param file  a string containing the pathname of a program to execute.
       * @param argv  a vector of argument strings passed to the new program.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, const argv_type&, pmode)
       */
      void
      open( const std::string& file,
            const argv_type& argv,
            pmode mode = pstdout )
      {
        this->do_open(file, argv, readable(mode));
      }

      /**
       * @brief Set streambuf to read from process' @c stdout.
       * @return  @c *this
       */
      basic_ipstream&
      out()
      {
        this->buf_.read_err(false);
        return *this;
      }

    };


  /**
   * @class basic_opstream
   * @brief Class template for Output PStreams.
   *
   * Writing to an open opstream writes to the standard input of the command;
   * the command's standard output is the same as that of the process that
   * created the pstream object, unless altered by the command itself.
   */

  template <typename CharT, typename Traits = std::char_traits<CharT> >
    class basic_opstream
    : public std::basic_ostream<CharT, Traits>
    , public pstream_common<CharT, Traits>
    , virtual public pstreams
    {
      typedef std::basic_ostream<CharT, Traits>     ostream_type;
      typedef pstream_common<CharT, Traits>         pbase_type;

      using pbase_type::buf_;  // declare name in this scope

    public:
      /// Type used to specify how to connect to the process.
      typedef typename pbase_type::pmode            pmode;

      /// Type used to hold the arguments for a command.
      typedef typename pbase_type::argv_type        argv_type;

      /// Default constructor, creates an uninitialised stream.
      basic_opstream()
      : ostream_type(NULL), pbase_type()
      { }

      /**
       * @brief Constructor that initialises the stream by starting a process.
       *
       * Initialises the stream buffer by calling do_open() with the supplied
       * arguments.
       *
       * @param cmd   a string containing a shell command.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, pmode)
       */
      explicit
      basic_opstream(const std::string& cmd, pmode mode = pstdin)
      : ostream_type(NULL), pbase_type(cmd, mode|pstdin)
      { }

      /**
       * @brief Constructor that initialises the stream by starting a process.
       *
       * Initialises the stream buffer by calling do_open() with the supplied
       * arguments.
       *
       * @param file  a string containing the pathname of a program to execute.
       * @param argv  a vector of argument strings passed to the new program.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, const argv_type&, pmode)
       */
      basic_opstream( const std::string& file,
                      const argv_type& argv,
                      pmode mode = pstdin )
      : ostream_type(NULL), pbase_type(file, argv, mode|pstdin)
      { }

      /**
       * @brief Constructor that initialises the stream by starting a process.
       *
       * Initialises the stream buffer by calling
       * @c do_open(argv[0],argv,mode|pstdin)
       *
       * @param argv  a vector of argument strings passed to the new program.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, const argv_type&, pmode)
       */
      explicit
      basic_opstream(const argv_type& argv, pmode mode = pstdin)
      : ostream_type(NULL), pbase_type(argv.at(0), argv, mode|pstdin)
      { }

#if __cplusplus >= 201103L
      /**
       * @brief Constructor that initialises the stream by starting a process.
       *
       * @param args  a list of argument strings passed to the new program.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, const argv_type&, pmode)
       */
      template<typename T>
        explicit
        basic_opstream(std::initializer_list<T> args, pmode mode = pstdin)
        : basic_opstream(argv_type(args.begin(), args.end()), mode)
        { }
#endif

      /**
       * @brief Destructor
       *
       * Closes the stream and waits for the child to exit.
       */
      ~basic_opstream() { }

      /**
       * @brief Start a process.
       *
       * Calls do_open( @a cmd , @a mode|pstdin ).
       *
       * @param cmd   a string containing a shell command.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, pmode)
       */
      void
      open(const std::string& cmd, pmode mode = pstdin)
      {
        this->do_open(cmd, mode|pstdin);
      }

      /**
       * @brief Start a process.
       *
       * Calls do_open( @a file , @a argv , @a mode|pstdin ).
       *
       * @param file  a string containing the pathname of a program to execute.
       * @param argv  a vector of argument strings passed to the new program.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, const argv_type&, pmode)
       */
      void
      open( const std::string& file,
            const argv_type& argv,
            pmode mode = pstdin)
      {
        this->do_open(file, argv, mode|pstdin);
      }
    };


  /**
   * @class basic_pstream
   * @brief Class template for Bidirectional PStreams.
   *
   * Writing to a pstream opened with @c pmode @c pstdin writes to the
   * standard input of the command.
   * Reading from a pstream opened with @c pmode @c pstdout and/or @c pstderr
   * reads the command's standard output and/or standard error.
   * Any of the process' @c stdin, @c stdout or @c stderr that is not
   * connected to the pstream (as specified by the @c pmode)
   * will be the same as the process that created the pstream object,
   * unless altered by the command itself.
   */
  template <typename CharT, typename Traits = std::char_traits<CharT> >
    class basic_pstream
    : public std::basic_iostream<CharT, Traits>
    , public pstream_common<CharT, Traits>
    , virtual public pstreams
    {
      typedef std::basic_iostream<CharT, Traits>    iostream_type;
      typedef pstream_common<CharT, Traits>         pbase_type;

      using pbase_type::buf_;  // declare name in this scope

    public:
      /// Type used to specify how to connect to the process.
      typedef typename pbase_type::pmode            pmode;

      /// Type used to hold the arguments for a command.
      typedef typename pbase_type::argv_type        argv_type;

      /// Default constructor, creates an uninitialised stream.
      basic_pstream()
      : iostream_type(NULL), pbase_type()
      { }

      /**
       * @brief Constructor that initialises the stream by starting a process.
       *
       * Initialises the stream buffer by calling do_open() with the supplied
       * arguments.
       *
       * @param cmd   a string containing a shell command.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, pmode)
       */
      explicit
      basic_pstream(const std::string& cmd, pmode mode = pstdout|pstdin)
      : iostream_type(NULL), pbase_type(cmd, mode)
      { }

      /**
       * @brief Constructor that initialises the stream by starting a process.
       *
       * Initialises the stream buffer by calling do_open() with the supplied
       * arguments.
       *
       * @param file  a string containing the pathname of a program to execute.
       * @param argv  a vector of argument strings passed to the new program.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, const argv_type&, pmode)
       */
      basic_pstream( const std::string& file,
                     const argv_type& argv,
                     pmode mode = pstdout|pstdin )
      : iostream_type(NULL), pbase_type(file, argv, mode)
      { }

      /**
       * @brief Constructor that initialises the stream by starting a process.
       *
       * Initialises the stream buffer by calling
       * @c do_open(argv[0],argv,mode)
       *
       * @param argv  a vector of argument strings passed to the new program.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, const argv_type&, pmode)
       */
      explicit
      basic_pstream(const argv_type& argv, pmode mode = pstdout|pstdin)
      : iostream_type(NULL), pbase_type(argv.at(0), argv, mode)
      { }

#if __cplusplus >= 201103L
      /**
       * @brief Constructor that initialises the stream by starting a process.
       *
       * @param l     a list of argument strings passed to the new program.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, const argv_type&, pmode)
       */
      template<typename T>
        explicit
        basic_pstream(std::initializer_list<T> l, pmode mode = pstdout|pstdin)
        : basic_pstream(argv_type(l.begin(), l.end()), mode)
        { }
#endif

      /**
       * @brief Destructor
       *
       * Closes the stream and waits for the child to exit.
       */
      ~basic_pstream() { }

      /**
       * @brief Start a process.
       *
       * Calls do_open( @a cnd , @a mode ).
       *
       * @param cmd   a string containing a shell command.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, pmode)
       */
      void
      open(const std::string& cmd, pmode mode = pstdout|pstdin)
      {
        this->do_open(cmd, mode);
      }

      /**
       * @brief Start a process.
       *
       * Calls do_open( @a file , @a argv , @a mode ).
       *
       * @param file  a string containing the pathname of a program to execute.
       * @param argv  a vector of argument strings passed to the new program.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, const argv_type&, pmode)
       */
      void
      open( const std::string& file,
            const argv_type& argv,
            pmode mode = pstdout|pstdin )
      {
        this->do_open(file, argv, mode);
      }

      /**
       * @brief Set streambuf to read from process' @c stdout.
       * @return  @c *this
       */
      basic_pstream&
      out()
      {
        this->buf_.read_err(false);
        return *this;
      }
    };


  /**
   * @class basic_rpstream
   * @brief Class template for Restricted PStreams.
   *
   * Writing to an rpstream opened with @c pmode @c pstdin writes to the
   * standard input of the command.
   * It is not possible to read directly from an rpstream object, to use
   * an rpstream as in istream you must call either basic_rpstream::out()
   * or basic_rpstream::err(). This is to prevent accidental reads from
   * the wrong input source. If the rpstream was not opened with @c pmode
   * @c pstderr then the class cannot read the process' @c stderr, and
   * basic_rpstream::err() will return an istream that reads from the
   * process' @c stdout, and vice versa.
   * Reading from an rpstream opened with @c pmode @c pstdout and/or
   * @c pstderr reads the command's standard output and/or standard error.
   * Any of the process' @c stdin, @c stdout or @c stderr that is not
   * connected to the pstream (as specified by the @c pmode)
   * will be the same as the process that created the pstream object,
   * unless altered by the command itself.
   */

  template <typename CharT, typename Traits = std::char_traits<CharT> >
    class basic_rpstream
    : public std::basic_ostream<CharT, Traits>
    , private std::basic_istream<CharT, Traits>
    , private pstream_common<CharT, Traits>
    , virtual public pstreams
    {
      typedef std::basic_ostream<CharT, Traits>     ostream_type;
      typedef std::basic_istream<CharT, Traits>     istream_type;
      typedef pstream_common<CharT, Traits>         pbase_type;

      using pbase_type::buf_;  // declare name in this scope

    public:
      /// Type used to specify how to connect to the process.
      typedef typename pbase_type::pmode            pmode;

      /// Type used to hold the arguments for a command.
      typedef typename pbase_type::argv_type        argv_type;

      /// Default constructor, creates an uninitialised stream.
      basic_rpstream()
      : ostream_type(NULL), istream_type(NULL), pbase_type()
      { }

      /**
       * @brief  Constructor that initialises the stream by starting a process.
       *
       * Initialises the stream buffer by calling do_open() with the supplied
       * arguments.
       *
       * @param cmd   a string containing a shell command.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, pmode)
       */
      explicit
      basic_rpstream(const std::string& cmd, pmode mode = pstdout|pstdin)
      : ostream_type(NULL) , istream_type(NULL) , pbase_type(cmd, mode)
      { }

      /**
       * @brief  Constructor that initialises the stream by starting a process.
       *
       * Initialises the stream buffer by calling do_open() with the supplied
       * arguments.
       *
       * @param file a string containing the pathname of a program to execute.
       * @param argv a vector of argument strings passed to the new program.
       * @param mode the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, const argv_type&, pmode)
       */
      basic_rpstream( const std::string& file,
                      const argv_type& argv,
                      pmode mode = pstdout|pstdin )
      : ostream_type(NULL), istream_type(NULL), pbase_type(file, argv, mode)
      { }

      /**
       * @brief Constructor that initialises the stream by starting a process.
       *
       * Initialises the stream buffer by calling
       * @c do_open(argv[0],argv,mode)
       *
       * @param argv  a vector of argument strings passed to the new program.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, const argv_type&, pmode)
       */
      explicit
      basic_rpstream(const argv_type& argv, pmode mode = pstdout|pstdin)
      : ostream_type(NULL), istream_type(NULL),
        pbase_type(argv.at(0), argv, mode)
      { }

#if __cplusplus >= 201103L
      /**
       * @brief Constructor that initialises the stream by starting a process.
       *
       * @param l     a list of argument strings passed to the new program.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, const argv_type&, pmode)
       */
      template<typename T>
        explicit
        basic_rpstream(std::initializer_list<T> l, pmode mode = pstdout|pstdin)
        : basic_rpstream(argv_type(l.begin(), l.end()), mode)
        { }
#endif

      /// Destructor
      ~basic_rpstream() { }

      /**
       * @brief  Start a process.
       *
       * Calls do_open( @a cmd , @a mode ).
       *
       * @param cmd   a string containing a shell command.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, pmode)
       */
      void
      open(const std::string& cmd, pmode mode = pstdout|pstdin)
      {
        this->do_open(cmd, mode);
      }

      /**
       * @brief  Start a process.
       *
       * Calls do_open( @a file , @a argv , @a mode ).
       *
       * @param file a string containing the pathname of a program to execute.
       * @param argv a vector of argument strings passed to the new program.
       * @param mode the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, const argv_type&, pmode)
       */
      void
      open( const std::string& file,
            const argv_type& argv,
            pmode mode = pstdout|pstdin )
      {
        this->do_open(file, argv, mode);
      }

      /**
       * @brief  Obtain a reference to the istream that reads
       *         the process' @c stdout.
       * @return @c *this
       */
      istream_type&
      out()
      {
        this->buf_.read_err(false);
        return *this;
      }
    };


  /// Type definition for common template specialisation.
  typedef basic_pstreambuf<char> pstreambuf;
  /// Type definition for common template specialisation.
  typedef basic_ipstream<char> ipstream;
  /// Type definition for common template specialisation.
  typedef basic_opstream<char> opstream;
  /// Type definition for common template specialisation.
  typedef basic_pstream<char> pstream;
  /// Type definition for common template specialisation.
  typedef basic_rpstream<char> rpstream;


  /**
   * When inserted into an output pstream the manipulator calls
   * basic_pstreambuf<C,T>::peof() to close the output pipe,
   * causing the child process to receive the end-of-file indicator
   * on subsequent reads from its @c stdin stream.
   *
   * @brief   Manipulator to close the pipe connected to the process' stdin.
   * @param   s  An output PStream class.
   * @return  The stream object the manipulator was invoked on.
   * @see     basic_pstreambuf<C,T>::peof()
   * @relates basic_opstream basic_pstream basic_rpstream
   */
  template <typename C, typename T>
    inline std::basic_ostream<C,T>&
    peof(std::basic_ostream<C,T>& s)
    {
      typedef basic_pstreambuf<C,T> pstreambuf_type;
      if (pstreambuf_type* p = dynamic_cast<pstreambuf_type*>(s.rdbuf()))
        p->peof();
      return s;
    }


  /*
   * member definitions for pstreambuf
   */


  /**
   * @class basic_pstreambuf
   * Provides underlying streambuf functionality for the PStreams classes.
   */

  /** Creates an uninitialised stream buffer. */
  template <typename C, typename T>
    inline
    basic_pstreambuf<C,T>::basic_pstreambuf()
    : ppid_(-1)   // initialise to -1 to indicate no process run yet.
    , wpipe_(-1)
    , wbuffer_(NULL)
    , rsrc_(rsrc_out)
    , status_(-1)
    , error_(0)
    {
      init_rbuffers();
    }

  /**
   * Initialises the stream buffer by calling open() with the supplied
   * arguments.
   *
   * @param cmd   a string containing a shell command.
   * @param mode  the I/O mode to use when opening the pipe.
   * @see   open()
   */
  template <typename C, typename T>
    inline
    basic_pstreambuf<C,T>::basic_pstreambuf(const std::string& cmd, pmode mode)
    : ppid_(-1)   // initialise to -1 to indicate no process run yet.
    , wpipe_(-1)
    , wbuffer_(NULL)
    , rsrc_(rsrc_out)
    , status_(-1)
    , error_(0)
    {
      init_rbuffers();
      open(cmd, mode);
    }

  /**
   * Initialises the stream buffer by calling open() with the supplied
   * arguments.
   *
   * @param file  a string containing the name of a program to execute.
   * @param argv  a vector of argument strings passsed to the new program.
   * @param mode  the I/O mode to use when opening the pipe.
   * @see   open()
   */
  template <typename C, typename T>
    inline
    basic_pstreambuf<C,T>::basic_pstreambuf( const std::string& file,
                                             const argv_type& argv,
                                             pmode mode )
    : ppid_(-1)   // initialise to -1 to indicate no process run yet.
    , wpipe_(-1)
    , wbuffer_(NULL)
    , rsrc_(rsrc_out)
    , status_(-1)
    , error_(0)
    {
      init_rbuffers();
      open(file, argv, mode);
    }

  /**
   * Closes the stream by calling close().
   * @see close()
   */
  template <typename C, typename T>
    inline
    basic_pstreambuf<C,T>::~basic_pstreambuf()
    {
      close();
    }

  /**
   * Starts a new process by passing @a command to the shell (/bin/sh)
   * and opens pipes to the process with the specified @a mode.
   *
   * If @a mode contains @c pstdout the initial read source will be
   * the child process' stdout, otherwise if @a mode  contains @c pstderr
   * the initial read source will be the child's stderr.
   *
   * Will duplicate the actions of  the  shell  in searching for an
   * executable file if the specified file name does not contain a slash (/)
   * character.
   *
   * @warning
   * There is no way to tell whether the shell command succeeded, this
   * function will always succeed unless resource limits (such as
   * memory usage, or number of processes or open files) are exceeded.
   * This means is_open() will return true even if @a command cannot
   * be executed.
   * Use pstreambuf::open(const std::string&, const argv_type&, pmode)
   * if you need to know whether the command failed to execute.
   *
   * @param   command  a string containing a shell command.
   * @param   mode     a bitwise OR of one or more of @c out, @c in, @c err.
   * @return  NULL if the shell could not be started or the
   *          pipes could not be opened, @c this otherwise.
   * @see     <b>execl</b>(3)
   */
  template <typename C, typename T>
    basic_pstreambuf<C,T>*
    basic_pstreambuf<C,T>::open(const std::string& command, pmode mode)
    {
      const char * shell_path = "/bin/sh";
#if 0
      const std::string argv[] = { "sh", "-c", command };
      return this->open(shell_path, argv_type(argv, argv+3), mode);
#else
      basic_pstreambuf<C,T>* ret = NULL;

      if (!is_open())
      {
        switch(fork(mode))
        {
        case 0 :
          // this is the new process, exec command
          ::execl(shell_path, "sh", "-c", command.c_str(), (char*)NULL);

          // can only reach this point if exec() failed

          // parent can get exit code from waitpid()
          ::_exit(errno);
          // using std::exit() would make static dtors run twice

        case -1 :
          // couldn't fork, error already handled in pstreambuf::fork()
          break;

        default :
          // this is the parent process
          // activate buffers
          create_buffers(mode);
          ret = this;
        }
      }
      return ret;
#endif
    }

  /**
   * @brief  Helper function to close a file descriptor.
   *
   * Inspects @a fd and calls <b>close</b>(3) if it has a non-negative value.
   *
   * @param   fd  a file descriptor.
   * @relates basic_pstreambuf
   */
  inline void
  close_fd(pstreams::fd_type& fd)
  {
    if (fd >= 0 && ::close(fd) == 0)
      fd = -1;
  }

  /**
   * @brief  Helper function to close an array of file descriptors.
   *
   * Calls @c close_fd() on each member of the array.
   * The length of the array is determined automatically by
   * template argument deduction to avoid errors.
   *
   * @param   fds  an array of file descriptors.
   * @relates basic_pstreambuf
   */
  template <int N>
    inline void
    close_fd_array(pstreams::fd_type (&fds)[N])
    {
      for (std::size_t i = 0; i < N; ++i)
        close_fd(fds[i]);
    }

  /**
   * Starts a new process by executing @a file with the arguments in
   * @a argv and opens pipes to the process with the specified @a mode.
   *
   * By convention @c argv[0] should be the file name of the file being
   * executed.
   *
   * If @a mode contains @c pstdout the initial read source will be
   * the child process' stdout, otherwise if @a mode  contains @c pstderr
   * the initial read source will be the child's stderr.
   *
   * Will duplicate the actions of  the  shell  in searching for an
   * executable file if the specified file name does not contain a slash (/)
   * character.
   *
   * Iff @a file is successfully executed then is_open() will return true.
   * Otherwise, pstreambuf::error() can be used to obtain the value of
   * @c errno that was set by <b>execvp</b>(3) in the child process.
   *
   * The exit status of the new process will be returned by
   * pstreambuf::status() after pstreambuf::exited() returns true.
   *
   * @param   file  a string containing the pathname of a program to execute.
   * @param   argv  a vector of argument strings passed to the new program.
   * @param   mode  a bitwise OR of one or more of @c out, @c in and @c err.
   * @return  NULL if a pipe could not be opened or if the program could
   *          not be executed, @c this otherwise.
   * @see     <b>execvp</b>(3)
   */
  template <typename C, typename T>
    basic_pstreambuf<C,T>*
    basic_pstreambuf<C,T>::open( const std::string& file,
                                 const argv_type& argv,
                                 pmode mode )
    {
      basic_pstreambuf<C,T>* ret = NULL;

      if (!is_open())
      {
        // constants for read/write ends of pipe
        enum { RD, WR };

        // open another pipe and set close-on-exec
        fd_type ck_exec[] = { -1, -1 };
        if (-1 == ::pipe(ck_exec)
            || -1 == ::fcntl(ck_exec[RD], F_SETFD, FD_CLOEXEC)
            || -1 == ::fcntl(ck_exec[WR], F_SETFD, FD_CLOEXEC))
        {
          error_ = errno;
          close_fd_array(ck_exec);
        }
        else
        {
          switch(fork(mode))
          {
          case 0 :
            // this is the new process, exec command
            {
              char** arg_v = new char*[argv.size()+1];
              for (std::size_t i = 0; i < argv.size(); ++i)
              {
                const std::string& src = argv[i];
                char*& dest = arg_v[i];
                dest = new char[src.size()+1];
                dest[ src.copy(dest, src.size()) ] = '\0';
              }
              arg_v[argv.size()] = NULL;

              ::execvp(file.c_str(), arg_v);

              // can only reach this point if exec() failed

              // parent can get error code from ck_exec pipe
              error_ = errno;

              while (::write(ck_exec[WR], &error_, sizeof(error_)) == -1
                  && errno == EINTR)
              { }

              ::close(ck_exec[WR]);
              ::close(ck_exec[RD]);

              ::_exit(error_);
              // using std::exit() would make static dtors run twice
            }

          case -1 :
            // couldn't fork, error already handled in pstreambuf::fork()
            close_fd_array(ck_exec);
            break;

          default :
            // this is the parent process

            // check child called exec() successfully
            ::close(ck_exec[WR]);
            switch (::read(ck_exec[RD], &error_, sizeof(error_)))
            {
            case 0:
              // activate buffers
              create_buffers(mode);
              ret = this;
              break;
            case -1:
              error_ = errno;
              break;
            default:
              // error_ contains error code from child
              // call wait() to clean up and set ppid_ to 0
              this->wait();
              break;
            }
            ::close(ck_exec[RD]);
          }
        }
      }
      return ret;
    }

  /**
   * Creates pipes as specified by @a mode and calls @c fork() to create
   * a new process. If the fork is successful the parent process stores
   * the child's PID and the opened pipes and the child process replaces
   * its standard streams with the opened pipes.
   *
   * If an error occurs the error code will be set to one of the possible
   * errors for @c pipe() or @c fork().
   * See your system's documentation for these error codes.
   *
   * @param   mode  an OR of pmodes specifying which of the child's
   *                standard streams to connect to.
   * @return  On success the PID of the child is returned in the parent's
   *          context and zero is returned in the child's context.
   *          On error -1 is returned and the error code is set appropriately.
   */
  template <typename C, typename T>
    pid_t
    basic_pstreambuf<C,T>::fork(pmode mode)
    {
      pid_t pid = -1;

      // Three pairs of file descriptors, for pipes connected to the
      // process' stdin, stdout and stderr
      // (stored in a single array so close_fd_array() can close all at once)
      fd_type fd[] = { -1, -1, -1, -1, -1, -1 };
      fd_type* const pin = fd;
      fd_type* const pout = fd+2;
      fd_type* const perr = fd+4;

      // constants for read/write ends of pipe
      enum { RD, WR };

      // N.B.
      // For the pstreambuf pin is an output stream and
      // pout and perr are input streams.

      if (!error_ && mode&pstdin && ::pipe(pin))
        error_ = errno;

      if (!error_ && mode&pstdout && ::pipe(pout))
        error_ = errno;

      if (!error_ && mode&pstderr && ::pipe(perr))
        error_ = errno;

      if (!error_)
      {
        pid = ::fork();
        switch (pid)
        {
          case 0 :
          {
            // this is the new process

            // for each open pipe close one end and redirect the
            // respective standard stream to the other end

            if (*pin >= 0)
            {
              ::close(pin[WR]);
              ::dup2(pin[RD], STDIN_FILENO);
              ::close(pin[RD]);
            }
            if (*pout >= 0)
            {
              ::close(pout[RD]);
              ::dup2(pout[WR], STDOUT_FILENO);
              ::close(pout[WR]);
            }
            if (*perr >= 0)
            {
              ::close(perr[RD]);
              ::dup2(perr[WR], STDERR_FILENO);
              ::close(perr[WR]);
            }

#ifdef _POSIX_JOB_CONTROL
            if (mode&newpg)
              ::setpgid(0, 0); // Change to a new process group
#endif

            break;
          }
          case -1 :
          {
            // couldn't fork for some reason
            error_ = errno;
            // close any open pipes
            close_fd_array(fd);
            break;
          }
          default :
          {
            // this is the parent process, store process' pid
            ppid_ = pid;

            // store one end of open pipes and close other end
            if (*pin >= 0)
            {
              wpipe_ = pin[WR];
              ::close(pin[RD]);
            }
            if (*pout >= 0)
            {
              rpipe_[rsrc_out] = pout[RD];
              ::close(pout[WR]);
            }
            if (*perr >= 0)
            {
              rpipe_[rsrc_err] = perr[RD];
              ::close(perr[WR]);
            }
          }
        }
      }
      else
      {
        // close any pipes we opened before failure
        close_fd_array(fd);
      }
      return pid;
    }

  /**
   * Closes all pipes and calls wait() to wait for the process to finish.
   * If an error occurs the error code will be set to one of the possible
   * errors for @c waitpid().
   * See your system's documentation for these errors.
   *
   * @return  @c this on successful close or @c NULL if there is no
   *          process to close or if an error occurs.
   */
  template <typename C, typename T>
    basic_pstreambuf<C,T>*
    basic_pstreambuf<C,T>::close()
    {
      const bool running = is_open();

      sync(); // this might call wait() and reap the child process

      // rather than trying to work out whether or not we need to clean up
      // just do it anyway, all cleanup functions are safe to call twice.

      destroy_buffers(pstdin|pstdout|pstderr);

      // close pipes before wait() so child gets EOF/SIGPIPE
      close_fd(wpipe_);
      close_fd_array(rpipe_);

      do
      {
        error_ = 0;
      } while (wait() == -1 && error() == EINTR);

      return running ? this : NULL;
    }

  /**
   *  Called on construction to initialise the arrays used for reading.
   */
  template <typename C, typename T>
    inline void
    basic_pstreambuf<C,T>::init_rbuffers()
    {
      rpipe_[rsrc_out] = rpipe_[rsrc_err] = -1;
      rbuffer_[rsrc_out] = rbuffer_[rsrc_err] = NULL;
      rbufstate_[0] = rbufstate_[1] = rbufstate_[2] = NULL;
    }

  template <typename C, typename T>
    void
    basic_pstreambuf<C,T>::create_buffers(pmode mode)
    {
      if (mode & pstdin)
      {
        delete[] wbuffer_;
        wbuffer_ = new char_type[bufsz];
        this->setp(wbuffer_, wbuffer_ + bufsz);
      }
      if (mode & pstdout)
      {
        delete[] rbuffer_[rsrc_out];
        rbuffer_[rsrc_out] = new char_type[bufsz];
        rsrc_ = rsrc_out;
        this->setg(rbuffer_[rsrc_out] + pbsz, rbuffer_[rsrc_out] + pbsz,
            rbuffer_[rsrc_out] + pbsz);
      }
      if (mode & pstderr)
      {
        delete[] rbuffer_[rsrc_err];
        rbuffer_[rsrc_err] = new char_type[bufsz];
        if (!(mode & pstdout))
        {
          rsrc_ = rsrc_err;
          this->setg(rbuffer_[rsrc_err] + pbsz, rbuffer_[rsrc_err] + pbsz,
              rbuffer_[rsrc_err] + pbsz);
        }
      }
    }

  template <typename C, typename T>
    void
    basic_pstreambuf<C,T>::destroy_buffers(pmode mode)
    {
      if (mode & pstdin)
      {
        this->setp(NULL, NULL);
        delete[] wbuffer_;
        wbuffer_ = NULL;
      }
      if (mode & pstdout)
      {
        if (rsrc_ == rsrc_out)
          this->setg(NULL, NULL, NULL);
        delete[] rbuffer_[rsrc_out];
        rbuffer_[rsrc_out] = NULL;
      }
      if (mode & pstderr)
      {
        if (rsrc_ == rsrc_err)
          this->setg(NULL, NULL, NULL);
        delete[] rbuffer_[rsrc_err];
        rbuffer_[rsrc_err] = NULL;
      }
    }

  template <typename C, typename T>
    typename basic_pstreambuf<C,T>::buf_read_src
    basic_pstreambuf<C,T>::switch_read_buffer(buf_read_src src)
    {
      if (rsrc_ != src)
      {
        char_type* tmpbufstate[] = {this->eback(), this->gptr(), this->egptr()};
        this->setg(rbufstate_[0], rbufstate_[1], rbufstate_[2]);
        for (std::size_t i = 0; i < 3; ++i)
          rbufstate_[i] = tmpbufstate[i];
        rsrc_ = src;
      }
      return rsrc_;
    }

  /**
   * Suspends execution and waits for the associated process to exit, or
   * until a signal is delivered whose action is to terminate the current
   * process or to call a signal handling function. If the process has
   * already exited (i.e. it is a "zombie" process) then wait() returns
   * immediately.  Waiting for the child process causes all its system
   * resources to be freed.
   *
   * error() will return EINTR if wait() is interrupted by a signal.
   *
   * @param   nohang  true to return immediately if the process has not exited.
   * @return  1 if the process has exited and wait() has not yet been called.
   *          0 if @a nohang is true and the process has not exited yet.
   *          -1 if no process has been started or if an error occurs,
   *          in which case the error can be found using error().
   */
  template <typename C, typename T>
    int
    basic_pstreambuf<C,T>::wait(bool nohang)
    {
      int child_exited = -1;
      if (is_open())
      {
        int exit_status;
        switch(::waitpid(ppid_, &exit_status, nohang ? WNOHANG : 0))
        {
          case 0 :
            // nohang was true and process has not exited
            child_exited = 0;
            break;
          case -1 :
            error_ = errno;
            break;
          default :
            // process has exited
            ppid_ = 0;
            status_ = exit_status;
            child_exited = 1;
            // Close wpipe, would get SIGPIPE if we used it.
            destroy_buffers(pstdin);
            close_fd(wpipe_);
            // Must free read buffers and pipes on destruction
            // or next call to open()/close()
            break;
        }
      }
      return child_exited;
    }

  /**
   * Sends the specified signal to the process.  A signal can be used to
   * terminate a child process that would not exit otherwise.
   *
   * If an error occurs the error code will be set to one of the possible
   * errors for @c kill().  See your system's documentation for these errors.
   *
   * @param   signal  A signal to send to the child process.
   * @return  @c this or @c NULL if @c kill() fails.
   */
  template <typename C, typename T>
    inline basic_pstreambuf<C,T>*
    basic_pstreambuf<C,T>::kill(int signal)
    {
      basic_pstreambuf<C,T>* ret = NULL;
      if (is_open())
      {
        if (::kill(ppid_, signal))
          error_ = errno;
        else
        {
#if 0
          // TODO call exited() to check for exit and clean up? leave to user?
          if (signal==SIGTERM || signal==SIGKILL)
            this->exited();
#endif
          ret = this;
        }
      }
      return ret;
    }

  /**
   * Sends the specified signal to the process group of the child process.
   * A signal can be used to terminate a child process that would not exit
   * otherwise, or to kill the process and its own children.
   *
   * If an error occurs the error code will be set to one of the possible
   * errors for @c getpgid() or @c kill().  See your system's documentation
   * for these errors. If the child is in the current process group then
   * NULL will be returned and the error code set to EPERM.
   *
   * @param   signal  A signal to send to the child process.
   * @return  @c this on success or @c NULL on failure.
   */
  template <typename C, typename T>
    inline basic_pstreambuf<C,T>*
    basic_pstreambuf<C,T>::killpg(int signal)
    {
      basic_pstreambuf<C,T>* ret = NULL;
#ifdef _POSIX_JOB_CONTROL
      if (is_open())
      {
        pid_t pgid = ::getpgid(ppid_);
        if (pgid == -1)
          error_ = errno;
        else if (pgid == ::getpgrp())
          error_ = EPERM;  // Don't commit suicide
        else if (::killpg(pgid, signal))
          error_ = errno;
        else
          ret = this;
      }
#else
      error_ = ENOTSUP;
#endif
      return ret;
    }

  /**
   *  This function can call pstreambuf::wait() and so may change the
   *  object's state if the child process has already exited.
   *
   *  @return  True if the associated process has exited, false otherwise.
   *  @see     basic_pstreambuf<C,T>::wait()
   */
  template <typename C, typename T>
    inline bool
    basic_pstreambuf<C,T>::exited()
    {
      return ppid_ == 0 || wait(true)==1;
    }

  /**
   *  @return  The error code of the most recently failed operation, or zero.
   */
  template <typename C, typename T>
    inline int
    basic_pstreambuf<C,T>::error() const
    {
      return error_;
    }

  /**
   *  Closes the output pipe, causing the child process to receive the
   *  end-of-file indicator on subsequent reads from its @c stdin stream.
   */
  template <typename C, typename T>
    inline void
    basic_pstreambuf<C,T>::peof()
    {
      sync();
      destroy_buffers(pstdin);
      close_fd(wpipe_);
    }

  /**
   * Unlike pstreambuf::exited(), this function will not call wait() and
   * so will not change the object's state.  This means that once a child
   * process is executed successfully this function will continue to
   * return true even after the process exits (until wait() is called.)
   *
   * @return  true if a previous call to open() succeeded and wait() has
   *          not been called and determined that the process has exited,
   *          false otherwise.
   */
  template <typename C, typename T>
    inline bool
    basic_pstreambuf<C,T>::is_open() const
    {
      return ppid_ > 0;
    }

  /**
   * Toggle the stream used for reading. If @a readerr is @c true then the
   * process' @c stderr output will be used for subsequent extractions, if
   * @a readerr is false the the process' stdout will be used.
   * @param   readerr  @c true to read @c stderr, @c false to read @c stdout.
   * @return  @c true if the requested stream is open and will be used for
   *          subsequent extractions, @c false otherwise.
   */
  template <typename C, typename T>
    inline bool
    basic_pstreambuf<C,T>::read_err(bool readerr)
    {
      buf_read_src src = readerr ? rsrc_err : rsrc_out;
      if (rpipe_[src]>=0)
      {
        switch_read_buffer(src);
        return true;
      }
      return false;
    }

  /**
   * Called when the internal character buffer is not present or is full,
   * to transfer the buffer contents to the pipe.
   *
   * @param   c  a character to be written to the pipe.
   * @return  @c traits_type::eof() if an error occurs, otherwise if @a c
   *          is not equal to @c traits_type::eof() it will be buffered and
   *          a value other than @c traits_type::eof() returned to indicate
   *          success.
   */
  template <typename C, typename T>
    typename basic_pstreambuf<C,T>::int_type
    basic_pstreambuf<C,T>::overflow(int_type c)
    {
      if (!empty_buffer())
        return traits_type::eof();
      else if (!traits_type::eq_int_type(c, traits_type::eof()))
        return this->sputc(c);
      else
        return traits_type::not_eof(c);
    }


  template <typename C, typename T>
    int
    basic_pstreambuf<C,T>::sync()
    {
      return !exited() && empty_buffer() ? 0 : -1;
    }

  /**
   * @param   s  character buffer.
   * @param   n  buffer length.
   * @return  the number of characters written.
   */
  template <typename C, typename T>
    std::streamsize
    basic_pstreambuf<C,T>::xsputn(const char_type* s, std::streamsize n)
    {
      std::streamsize done = 0;
      while (done < n)
      {
        if (std::streamsize nbuf = this->epptr() - this->pptr())
        {
          nbuf = std::min(nbuf, n - done);
          traits_type::copy(this->pptr(), s + done, nbuf);
          this->pbump(nbuf);
          done += nbuf;
        }
        else if (!empty_buffer())
          break;
      }
      return done;
    }

  /**
   * @return  true if the buffer was emptied, false otherwise.
   */
  template <typename C, typename T>
    bool
    basic_pstreambuf<C,T>::empty_buffer()
    {
      const std::streamsize count = this->pptr() - this->pbase();
      if (count > 0)
      {
        const std::streamsize written = this->write(this->wbuffer_, count);
        if (written > 0)
        {
          if (const std::streamsize unwritten = count - written)
            traits_type::move(this->pbase(), this->pbase()+written, unwritten);
          this->pbump(-written);
          return true;
        }
      }
      return false;
    }

  /**
   * Called when the internal character buffer is is empty, to re-fill it
   * from the pipe.
   *
   * @return The first available character in the buffer,
   * or @c traits_type::eof() in case of failure.
   */
  template <typename C, typename T>
    typename basic_pstreambuf<C,T>::int_type
    basic_pstreambuf<C,T>::underflow()
    {
      if (this->gptr() < this->egptr() || fill_buffer())
        return traits_type::to_int_type(*this->gptr());
      else
        return traits_type::eof();
    }

  /**
   * Attempts to make @a c available as the next character to be read by
   * @c sgetc().
   *
   * @param   c   a character to make available for extraction.
   * @return  @a c if the character can be made available,
   *          @c traits_type::eof() otherwise.
   */
  template <typename C, typename T>
    typename basic_pstreambuf<C,T>::int_type
    basic_pstreambuf<C,T>::pbackfail(int_type c)
    {
      if (this->gptr() != this->eback())
      {
        this->gbump(-1);
        if (!traits_type::eq_int_type(c, traits_type::eof()))
          *this->gptr() = traits_type::to_char_type(c);
        return traits_type::not_eof(c);
      }
      else
         return traits_type::eof();
    }

  template <typename C, typename T>
    std::streamsize
    basic_pstreambuf<C,T>::showmanyc()
    {
      int avail = 0;
      if (sizeof(char_type) == 1)
        avail = fill_buffer(true) ? this->egptr() - this->gptr() : -1;
#ifdef FIONREAD
      else
      {
        if (::ioctl(rpipe(), FIONREAD, &avail) == -1)
          avail = -1;
        else if (avail)
          avail /= sizeof(char_type);
      }
#endif
      return std::streamsize(avail);
    }

  /**
   * @return  true if the buffer was filled, false otherwise.
   */
  template <typename C, typename T>
    bool
    basic_pstreambuf<C,T>::fill_buffer(bool non_blocking)
    {
      const std::streamsize pb1 = this->gptr() - this->eback();
      const std::streamsize pb2 = pbsz;
      const std::streamsize npb = std::min(pb1, pb2);

      char_type* const rbuf = rbuffer();

      if (npb)
        traits_type::move(rbuf + pbsz - npb, this->gptr() - npb, npb);

      std::streamsize rc = -1;

      if (non_blocking)
      {
        const int flags = ::fcntl(rpipe(), F_GETFL);
        if (flags != -1)
        {
          const bool blocking = !(flags & O_NONBLOCK);
          if (blocking)
            ::fcntl(rpipe(), F_SETFL, flags | O_NONBLOCK);  // set non-blocking

          error_ = 0;
          rc = read(rbuf + pbsz, bufsz - pbsz);

          if (rc == -1 && error_ == EAGAIN)  // nothing available
            rc = 0;
          else if (rc == 0)  // EOF
            rc = -1;

          if (blocking)
            ::fcntl(rpipe(), F_SETFL, flags); // restore
        }
      }
      else
        rc = read(rbuf + pbsz, bufsz - pbsz);

      if (rc > 0 || (rc == 0 && non_blocking))
      {
        this->setg( rbuf + pbsz - npb,
                    rbuf + pbsz,
                    rbuf + pbsz + rc );
        return true;
      }
      else
      {
        this->setg(NULL, NULL, NULL);
        return false;
      }
    }

  /**
   * Writes up to @a n characters to the pipe from the buffer @a s.
   *
   * @param   s  character buffer.
   * @param   n  buffer length.
   * @return  the number of characters written.
   */
  template <typename C, typename T>
    inline std::streamsize
    basic_pstreambuf<C,T>::write(const char_type* s, std::streamsize n)
    {
      std::streamsize nwritten = 0;
      if (wpipe() >= 0)
      {
        nwritten = ::write(wpipe(), s, n * sizeof(char_type));
        if (nwritten == -1)
          error_ = errno;
        else
          nwritten /= sizeof(char_type);
      }
      return nwritten;
    }

  /**
   * Reads up to @a n characters from the pipe to the buffer @a s.
   *
   * @param   s  character buffer.
   * @param   n  buffer length.
   * @return  the number of characters read.
   */
  template <typename C, typename T>
    inline std::streamsize
    basic_pstreambuf<C,T>::read(char_type* s, std::streamsize n)
    {
      std::streamsize nread = 0;
      if (rpipe() >= 0)
      {
        nread = ::read(rpipe(), s, n * sizeof(char_type));
        if (nread == -1)
          error_ = errno;
        else
          nread /= sizeof(char_type);
      }
      return nread;
    }

  /** @return a reference to the output file descriptor */
  template <typename C, typename T>
    inline pstreams::fd_type&
    basic_pstreambuf<C,T>::wpipe()
    {
      return wpipe_;
    }

  /** @return a reference to the active input file descriptor */
  template <typename C, typename T>
    inline pstreams::fd_type&
    basic_pstreambuf<C,T>::rpipe()
    {
      return rpipe_[rsrc_];
    }

  /** @return a reference to the specified input file descriptor */
  template <typename C, typename T>
    inline pstreams::fd_type&
    basic_pstreambuf<C,T>::rpipe(buf_read_src which)
    {
      return rpipe_[which];
    }

  /** @return a pointer to the start of the active input buffer area. */
  template <typename C, typename T>
    inline typename basic_pstreambuf<C,T>::char_type*
    basic_pstreambuf<C,T>::rbuffer()
    {
      return rbuffer_[rsrc_];
    }


  /*
   * member definitions for pstream_common
   */

  /**
   * @class pstream_common
   * Abstract Base Class providing common functionality for basic_ipstream,
   * basic_opstream and basic_pstream.
   * pstream_common manages the basic_pstreambuf stream buffer that is used
   * by the derived classes to initialise an iostream class.
   */

  /** Creates an uninitialised stream. */
  template <typename C, typename T>
    inline
    pstream_common<C,T>::pstream_common()
    : std::basic_ios<C,T>(NULL)
    , command_()
    , buf_()
    {
      this->std::basic_ios<C,T>::rdbuf(&buf_);
    }

  /**
   * Initialises the stream buffer by calling
   * do_open( @a command , @a mode )
   *
   * @param cmd   a string containing a shell command.
   * @param mode  the I/O mode to use when opening the pipe.
   * @see   do_open(const std::string&, pmode)
   */
  template <typename C, typename T>
    inline
    pstream_common<C,T>::pstream_common(const std::string& cmd, pmode mode)
    : std::basic_ios<C,T>(NULL)
    , command_(cmd)
    , buf_()
    {
      this->std::basic_ios<C,T>::rdbuf(&buf_);
      do_open(cmd, mode);
    }

  /**
   * Initialises the stream buffer by calling
   * do_open( @a file , @a argv , @a mode )
   *
   * @param file  a string containing the pathname of a program to execute.
   * @param argv  a vector of argument strings passed to the new program.
   * @param mode  the I/O mode to use when opening the pipe.
   * @see do_open(const std::string&, const argv_type&, pmode)
   */
  template <typename C, typename T>
    inline
    pstream_common<C,T>::pstream_common( const std::string& file,
                                         const argv_type& argv,
                                         pmode mode )
    : std::basic_ios<C,T>(NULL)
    , command_(file)
    , buf_()
    {
      this->std::basic_ios<C,T>::rdbuf(&buf_);
      do_open(file, argv, mode);
    }

  /**
   * This is a pure virtual function to make @c pstream_common abstract.
   * Because it is the destructor it will be called by derived classes
   * and so must be defined.  It is also protected, to discourage use of
   * the PStreams classes through pointers or references to the base class.
   *
   * @sa If defining a pure virtual seems odd you should read
   * http://www.gotw.ca/gotw/031.htm (and the rest of the site as well!)
   */
  template <typename C, typename T>
    inline
    pstream_common<C,T>::~pstream_common()
    {
    }

  /**
   * Calls rdbuf()->open( @a command , @a mode )
   * and sets @c failbit on error.
   *
   * @param cmd   a string containing a shell command.
   * @param mode  the I/O mode to use when opening the pipe.
   * @see   basic_pstreambuf::open(const std::string&, pmode)
   */
  template <typename C, typename T>
    inline void
    pstream_common<C,T>::do_open(const std::string& cmd, pmode mode)
    {
      if (!buf_.open((command_=cmd), mode))
        this->setstate(std::ios_base::failbit);
    }

  /**
   * Calls rdbuf()->open( @a file, @a  argv, @a mode )
   * and sets @c failbit on error.
   *
   * @param file  a string containing the pathname of a program to execute.
   * @param argv  a vector of argument strings passed to the new program.
   * @param mode  the I/O mode to use when opening the pipe.
   * @see   basic_pstreambuf::open(const std::string&, const argv_type&, pmode)
   */
  template <typename C, typename T>
    inline void
    pstream_common<C,T>::do_open( const std::string& file,
                                  const argv_type& argv,
                                  pmode mode )
    {
      if (!buf_.open((command_=file), argv, mode))
        this->setstate(std::ios_base::failbit);
    }

  /** Calls rdbuf->close() and sets @c failbit on error. */
  template <typename C, typename T>
    inline void
    pstream_common<C,T>::close()
    {
      if (!buf_.close())
        this->setstate(std::ios_base::failbit);
    }

  /**
   * @return  rdbuf()->is_open().
   * @see     basic_pstreambuf::is_open()
   */
  template <typename C, typename T>
    inline bool
    pstream_common<C,T>::is_open() const
    {
      return buf_.is_open();
    }

  /** @return a pointer to the private stream buffer member. */
  // TODO  document behaviour if buffer replaced.
  template <typename C, typename T>
    inline typename pstream_common<C,T>::streambuf_type*
    pstream_common<C,T>::rdbuf() const
    {
      return const_cast<streambuf_type*>(&buf_);
    }


#if REDI_EVISCERATE_PSTREAMS
  /**
   * @def REDI_EVISCERATE_PSTREAMS
   * If this macro has a non-zero value then certain internals of the
   * @c basic_pstreambuf template class are exposed. In general this is
   * a Bad Thing, as the internal implementation is largely undocumented
   * and may be subject to change at any time, so this feature is only
   * provided because it might make PStreams useful in situations where
   * it is necessary to do Bad Things.
   */

  /**
   * @warning  This function exposes the internals of the stream buffer and
   *           should be used with caution. It is the caller's responsibility
   *           to flush streams etc. in order to clear any buffered data.
   *           The POSIX.1 function <b>fdopen</b>(3) is used to obtain the
   *           @c FILE pointers from the streambuf's private file descriptor
   *           members so consult your system's documentation for
   *           <b>fdopen</b>(3).
   *
   * @param   in    A FILE* that will refer to the process' stdin.
   * @param   out   A FILE* that will refer to the process' stdout.
   * @param   err   A FILE* that will refer to the process' stderr.
   * @return  An OR of zero or more of @c pstdin, @c pstdout, @c pstderr.
   *
   * For each open stream shared with the child process a @c FILE* is
   * obtained and assigned to the corresponding parameter. For closed
   * streams @c NULL is assigned to the parameter.
   * The return value can be tested to see which parameters should be
   * @c !NULL by masking with the corresponding @c pmode value.
   *
   * @see <b>fdopen</b>(3)
   */
  template <typename C, typename T>
    std::size_t
    basic_pstreambuf<C,T>::fopen(FILE*& in, FILE*& out, FILE*& err)
    {
      in = out = err = NULL;
      std::size_t open_files = 0;
      if (wpipe() > -1)
      {
        if ((in = ::fdopen(wpipe(), "w")))
        {
            open_files |= pstdin;
        }
      }
      if (rpipe(rsrc_out) > -1)
      {
        if ((out = ::fdopen(rpipe(rsrc_out), "r")))
        {
            open_files |= pstdout;
        }
      }
      if (rpipe(rsrc_err) > -1)
      {
        if ((err = ::fdopen(rpipe(rsrc_err), "r")))
        {
            open_files |= pstderr;
        }
      }
      return open_files;
    }

  /**
   *  @warning This function exposes the internals of the stream buffer and
   *  should be used with caution.
   *
   *  @param  in   A FILE* that will refer to the process' stdin.
   *  @param  out  A FILE* that will refer to the process' stdout.
   *  @param  err  A FILE* that will refer to the process' stderr.
   *  @return A bitwise-or of zero or more of @c pstdin, @c pstdout, @c pstderr.
   *  @see    basic_pstreambuf::fopen()
   */
  template <typename C, typename T>
    inline std::size_t
    pstream_common<C,T>::fopen(FILE*& fin, FILE*& fout, FILE*& ferr)
    {
      return buf_.fopen(fin, fout, ferr);
    }

#endif // REDI_EVISCERATE_PSTREAMS


} // namespace redi


#endif  // REDI_PSTREAM_H_SEEN
#endif  // WIN32

/* pstream END */
/* main START */

inline string formatANISOU(const string &inputString)
{
    string result=Trim(inputString," ");
    size_t found=result.find_first_of('.');
    stringstream buf;
    if (found==string::npos)
    {
        result+=".0000";
        found=result.find_first_of('.');
    }
    string post_decimal=result.substr(found+1);
    int i;
    long int second=atoi(post_decimal.c_str());
    if (result[0]=='-') second=-second;
    if (post_decimal.size()>4)
    {
        double anisou_dbl=atof(result.c_str())*10000;
        long int anisou_int=round(anisou_dbl);
        if (rstrip(post_decimal.substr(4),"0")=="5" &&
            (post_decimal[3]-'0') % 2 ==0)  anisou_int=(long int)(anisou_dbl);
        buf<<right<<setw(7)<<anisou_int<<flush;
        result=buf.str();
        buf.str(string());
        post_decimal.clear();
        return result.substr(0,7);
    }
    else if (post_decimal.size()<4)
        for (i=0;i<4-post_decimal.size();i++) second*=10;
    
    string pre_decimal=result.substr(0,found);
    long int first=atoi(pre_decimal.c_str());
    buf<<right<<setw(7)<<first*10000+second<<flush;
    result=buf.str();
    buf.str(string());
    //if (result.size()>7) result=result.substr(7-result.size());
    if (result.size()>7) result=result.substr(0,7);
    pre_decimal.clear();
    post_decimal.clear();
    return result;
}

int read_semi_colon(vector<string> &line_vec, const int fields, int l,
    const vector<string> &lines, vector<string> &line_append_vec, string &line,
    const bool ignore_quotation=false,const bool add_space=false)
{
    int i;
    if (line_vec.size() && StartsWith(line_vec[0],";"))
    {
        l--;
        clear_line_vec(line_vec);
    }
    while (line_vec.size()<fields)
    {
        l++;
        if (StartsWith(lines[l],";"))
        {
            line="";
            while (l<lines.size())
            {
                if (StartsWith(lines[l],";"))
                {
                    if (Trim(lines[l])==";") break;
                    else line+=lines[l].substr(1);
                }
                else line+=lines[l];
                if (add_space) line+=' ';
                l++;
            }
            line_vec.push_back(line);
        }
        else
        {
            Split(lines[l],line_append_vec,' ',ignore_quotation);
            for (i=0;i<line_append_vec.size();i++)
                line_vec.push_back(line_append_vec[i]);
            clear_line_vec(line_append_vec);
        }
    }
    return l;
}

inline char aa3to1(const string resn)
{
    if (resn[0]==' ') return tolower(resn[2]);
    else if (resn=="PSU") return 'u';

    // 20 standard amino acid + MSE
    else if (resn=="ALA") return 'A';
    else if (resn=="CYS") return 'C';
    else if (resn=="ASP") return 'D';
    else if (resn=="GLU") return 'E';
    else if (resn=="PHE") return 'F';
    else if (resn=="GLY") return 'G';
    else if (resn=="HIS") return 'H';
    else if (resn=="ILE") return 'I';
    else if (resn=="LYS") return 'K';
    else if (resn=="LEU") return 'L';
    else if (resn=="MET") return 'M';
    else if (resn=="ASN") return 'N';
    else if (resn=="PRO") return 'P';
    else if (resn=="GLN") return 'Q';
    else if (resn=="ARG") return 'R';
    else if (resn=="SER") return 'S';
    else if (resn=="THR") return 'T';
    else if (resn=="VAL") return 'V'; 
    else if (resn=="TRP") return 'W';
    else if (resn=="TYR") return 'Y';

    if (resn=="MSE") return 'M';

    // non-standard amino acid with known parent
    if (resn=="CHG"||resn=="HAC"||resn=="AYA"||resn=="TIH"||resn=="BNN"||
        resn=="ALM"||resn=="TPQ"||resn=="MAA"||resn=="PRR"||resn=="FLA"||
        resn=="AIB"||resn=="DAL"||resn=="CSD"||resn=="DHA"||resn=="DNP") 
        return 'A';
    else if (resn=="PR3"||resn=="CCS"||resn=="C6C"||resn=="SMC"||resn=="BCS"||
             resn=="SCY"||resn=="DCY"||resn=="SCS"||resn=="CME"||resn=="CY1"||
             resn=="CYQ"||resn=="CEA"||resn=="CYG"||resn=="BUC"||resn=="PEC"||
             resn=="CYM"||resn=="CY3"||resn=="CSO"||resn=="SOC"||resn=="CSX"||
             resn=="CSW"||resn=="EFC"||resn=="CSP"||resn=="CSS"||resn=="SCH"||
             resn=="OCS"||resn=="SHC"||resn=="C5C") return 'C';
    else if (resn=="DGL"||resn=="GGL"||resn=="CGU"||resn=="GMA"||resn=="5HP"||
             resn=="PCA") return 'E';
    else if (resn=="ASQ"||resn=="ASB"||resn=="ASA"||resn=="ASK"||resn=="ASL"||
             resn=="2AS"||resn=="DAS"||resn=="DSP"||resn=="BHD") return 'D';
    else if (resn=="PHI"||resn=="PHL"||resn=="DPN"||resn=="DAH"||resn=="HPQ")
        return 'F';
    else if (resn=="GLZ"||resn=="SAR"||resn=="GSC"||resn=="GL3"||resn=="MSA"||
             resn=="MPQ"||resn=="NMC") return 'G';
    else if (resn=="NEM"||resn=="NEP"||resn=="HSD"||resn=="HSP"||resn=="MHS"||
             resn=="3AH"||resn=="HIC"||resn=="HIP"||resn=="DHI"||resn=="HSE") 
        return 'H';
    else if (resn=="IIL"||resn=="DIL") return 'I';
    else if (resn=="DLY"||resn=="LYZ"||resn=="SHR"||resn=="ALY"||resn=="TRG"||
             resn=="LYM"||resn=="LLY"||resn=="KCX") return 'K';
    else if (resn=="NLE"||resn=="CLE"||resn=="NLP"||resn=="DLE"||resn=="BUG"||
             resn=="NLN"||resn=="MLE") return 'L';
    else if (resn=="FME"||resn=="CXM"||resn=="OMT") return 'M';
    else if (resn=="MEN") return 'N';
    else if (resn=="DPR"||resn=="HYP") return 'P';
    else if (resn=="DGN") return 'Q';
    else if (resn=="AGM"||resn=="ACL"||resn=="DAR"||resn=="HAR"||resn=="HMR"||
             resn=="ARM") return 'R';
    else if (resn=="OAS"||resn=="MIS"||resn=="SAC"||resn=="SEL"||resn=="SVA"||
             resn=="SET"||resn=="DSN"||resn=="SEP") return 'S';
    else if (resn=="DTH"||resn=="TPO"||resn=="ALO"||resn=="BMT") return 'T';
    else if (resn=="DVA"||resn=="MVA"||resn=="DIV") return 'V';
    else if (resn=="LTR"||resn=="DTR"||resn=="TRO"||resn=="TPL"||resn=="HTR") 
        return 'W';
    else if (resn=="PAQ"||resn=="STY"||resn=="TYQ"||resn=="IYR"||resn=="TYY"||
             resn=="DTY"||resn=="TYB"||resn=="PTR"||resn=="TYS") return 'Y';
    
    // undeterminted amino acid
    else if (resn=="ASX") return 'B'; // or D or N
    else if (resn=="GLX") return 'Z'; // or Q or E
    else if (resn=="SEC") return 'U';
    else if (resn=="PYL") return 'O';
    return 'X';
}

int cif2fasta(const string &infile, string &pdbid, const int do_upper,
    const int do_gzip, const vector<string> &outputChain_vec)
{

    stringstream buf;
    if (infile=="-") buf<<cin.rdbuf();
#if defined(REDI_PSTREAM_H_SEEN)
    else if (EndsWith(infile,".gz"))
    {
        redi::ipstream fp_gz; // if file is compressed
        fp_gz.open("gunzip -c "+infile);
        buf<<fp_gz.rdbuf();
        fp_gz.close();
    }
#endif
    else
    {
        ifstream fp;
        fp.open(infile.c_str(),ios::in); //ifstream fp(filename,ios::in);
        buf<<fp.rdbuf();
        fp.close();
    }
    vector<string> lines;
    Split(buf.str(),lines,'\n',true); 
    buf.str(string());
    if (lines.size()<=1)
    {
        cerr<<"ERROR! Empty structure "<<infile<<endl;
        vector<string>().swap(lines);
        return 0;
    }

    /* parse ATOM/HETATM */
    map<string,int> _atom_site;
    string comp_id    ="UNK";   // auth_comp_id, label_comp_id (residue name)
    string asym_prev  ="";
    string asym_id    =" ";     // auth_asym_id, label_asym_id (chain ID)
    string seq_id     ="";      // label_seq_id, auth_seq_id (residue index)
    string seq_prev   ="";
    string pdbx_PDB_ins_code="";// (insertion code)
    string pdbx_PDB_model_num="1"; // model index

    string sequence="";
    vector<string> chainID_vec;
    vector<string> sequence_vec;
    vector<size_t> mol_type_vec(3,0); // protein, dna, rna
    vector<vector<size_t> >mol_type_mat;
    
    size_t l;
    int i,j;
    string line;
    vector<string> line_vec;
    for (l=0;l<lines.size();l++)
    {
        line=lines[l];
        Split(line,line_vec,' ',true);
        if (line_vec.size()==0) continue;
        else if (line_vec.size() && line_vec[0]=="#")
            _atom_site.clear();
        else if (pdbid.size()==0 && l==0 && StartsWith(line,"data_"))
            pdbid=Lower(line.substr(5));
        else if (pdbid.size()==0 && line_vec.size()>1 && line_vec[0]=="_entry.id")
            pdbid=Lower(line_vec[1]);
        else if (StartsWith(line,"_atom_site."))
        {
            line=line_vec[0];
            clear_line_vec(line_vec);
            Split(line,line_vec,'.');
            if (line_vec.size()>1)
            {
                line=line_vec[1];
                j=_atom_site.size();
                _atom_site[line]=j;
            }
        }
        else if (_atom_site.size())
        {
            if (_atom_site.count("pdbx_PDB_model_num"))
            {
                pdbx_PDB_model_num=line_vec[_atom_site["pdbx_PDB_model_num"]];
                if (pdbx_PDB_model_num!="." &&
                    pdbx_PDB_model_num=="?" && pdbx_PDB_model_num!="1")
                {
                    clear_line_vec(line_vec);
                    continue;
                }
            }
            if (_atom_site.count("label_seq_id") && 
                line_vec[_atom_site["label_seq_id"]]==".")
            {                
                clear_line_vec(line_vec);
                continue;
            }
            
            if (_atom_site.count("auth_asym_id"))
                asym_id=line_vec[_atom_site["auth_asym_id"]];
            else if (_atom_site.count("label_asym_id"))
                asym_id=line_vec[_atom_site["label_asym_id"]];
            else if (_atom_site.count("pdbx_auth_asym_id"))
                asym_id=line_vec[_atom_site["pdbx_auth_asym_id"]];
            else if (_atom_site.count("pdbx_label_asym_id"))
                asym_id=line_vec[_atom_site["pdbx_label_asym_id"]];
            if (asym_id=="." || asym_id=="?") asym_id="_";
            if (outputChain_vec.size() && find(outputChain_vec.begin(),
                outputChain_vec.end(), asym_id)==outputChain_vec.end())
            {
                clear_line_vec(line_vec);
                continue;
            }

            if (_atom_site.count("auth_seq_id"))
                seq_id=line_vec[_atom_site["auth_seq_id"]];
            else if (_atom_site.count("label_seq_id"))
                seq_id=line_vec[_atom_site["label_seq_id"]];
            else if (_atom_site.count("pdbx_auth_seq_id"))
                seq_id=line_vec[_atom_site["pdbx_auth_seq_id"]];
            else if (_atom_site.count("pdbx_label_seq_id"))
                seq_id=line_vec[_atom_site["pdbx_label_seq_id"]];

            if (_atom_site.count("pdbx_PDB_ins_code"))
            {
                pdbx_PDB_ins_code=line_vec[_atom_site["pdbx_PDB_ins_code"]];
                if (pdbx_PDB_ins_code!="." && pdbx_PDB_ins_code!="?")
                    seq_id+=pdbx_PDB_ins_code;
            }

            if (asym_prev==asym_id && seq_prev==seq_id)
            {
                clear_line_vec(line_vec);
                continue;
            }
            
            if (_atom_site.count("auth_comp_id"))
            {
                comp_id=line_vec[_atom_site["auth_comp_id"]];
                if (comp_id.size()>3 && _atom_site.count("label_comp_id"))
                    comp_id=line_vec[_atom_site["label_comp_id"]];
            }
            else if (_atom_site.count("label_comp_id"))
                comp_id=line_vec[_atom_site["label_comp_id"]];
            else if (_atom_site.count("pdbx_auth_comp_id"))
            {
                comp_id=line_vec[_atom_site["pdbx_auth_comp_id"]];
                if (comp_id.size()>3 && _atom_site.count("pdbx_label_comp_id"))
                    comp_id=line_vec[_atom_site["pdbx_label_comp_id"]];
            }
            else if (_atom_site.count("pdbx_label_comp_id"))
                comp_id=line_vec[_atom_site["pdbx_label_comp_id"]];
            while (comp_id.size()<3) comp_id=' '+comp_id;

            if (asym_prev!=asym_id)
            {
                if (asym_prev.size())
                {
                    chainID_vec.push_back(asym_prev);
                    sequence_vec.push_back(sequence);
                    mol_type_mat.push_back(mol_type_vec);
                    mol_type_vec[0]=mol_type_vec[1]=mol_type_vec[2]=0;
                    sequence="";
                }

                asym_prev=asym_id;
            }
            sequence+=aa3to1(comp_id);
            if (sequence[sequence.size()-1]!='X')
            {
                if      (comp_id.substr(0,2)==" D") mol_type_vec[1]++;
                else if (comp_id.substr(0,2)=="  ") mol_type_vec[2]++;
                else mol_type_vec[0]++;
            }
            seq_prev=seq_id;
        }

        /* clean up */
        clear_line_vec(line_vec);
        lines[l].clear();
    }
    lines.clear();
    if (sequence.size())
    {
        chainID_vec.push_back(asym_prev);
        sequence_vec.push_back(sequence);
        mol_type_mat.push_back(mol_type_vec);
    }

    if (pdbid.size()==0)
    {
        cerr<<"ERROR: no PDB ID in "<<infile<<'\n'
            <<"PDB ID can be specified by option -p=xxxx"<<endl;
        return -1;
    }

    size_t seqNum=chainID_vec.size();
    for (l=0;l<seqNum;l++)
    {
        buf<<'>'<<pdbid<<':'<<chainID_vec[l]<<'\t';
        if (mol_type_mat[l][0]>=mol_type_mat[l][1] && 
            mol_type_mat[l][1]>=mol_type_mat[l][2])
        {
            sequence=sequence_vec[l];
            if (do_upper==2) sequence=Upper(sequence);
            else if (do_upper==0) sequence=Lower(sequence);
            buf<<"PROTEIN\t"<<sequence.size()<<'\n'<<sequence<<'\n';
        }
        else
        {
            sequence="";
            for (j=0;j<sequence_vec[l].size();j++)
            {
                if (sequence_vec[l][j]=='X') sequence+='n';
                else sequence+=sequence_vec[l][j];
            }
            if (do_upper==2) sequence=Upper(sequence);
            else if (do_upper==0) sequence=Lower(sequence);
            if (mol_type_mat[l][1]>=mol_type_mat[l][2])
                 buf<<"DNA\t"<<sequence.size()<<'\n'<<sequence<<'\n';
            else buf<<"RNA\t"<<sequence.size()<<'\n'<<sequence<<'\n';
        }
    }
    buf<<flush;
    
    ofstream fout;
    string filename=pdbid+".fasta";
    fout.open(filename.c_str());
    fout<<buf.str()<<flush;
    buf.str(string());
    fout.close();
    cout<<filename<<endl;
    
    if (do_gzip)
    {
        line="gzip -f "+filename;
        j=system(line.c_str());
    }

    /* clean up */
    map<string,int> ().swap(_atom_site);
    string ().swap(filename);
    
    vector<string>().swap(chainID_vec);
    vector<string>().swap(sequence_vec);
    vector<size_t>().swap(mol_type_vec);
    vector<vector<size_t> >().swap(mol_type_mat);

    vector<string>().swap(lines);
    vector<string>().swap(line_vec);
    
    comp_id.clear();
    asym_prev.clear();
    asym_id.clear();
    seq_prev.clear();
    seq_id.clear();
    pdbx_PDB_ins_code.clear();
    pdbx_PDB_model_num.clear();
    line.clear();
    return seqNum;
}

int cif2pdb(const string &infile, string &pdbid)
{

    stringstream buf;
    if (infile=="-") buf<<cin.rdbuf();
#if defined(REDI_PSTREAM_H_SEEN)
    else if (EndsWith(infile,".gz"))
    {
        redi::ipstream fp_gz; // if file is compressed
        fp_gz.open("gunzip -c "+infile);
        buf<<fp_gz.rdbuf();
        fp_gz.close();
    }
#endif
    else
    {
        ifstream fp;
        fp.open(infile.c_str(),ios::in); //ifstream fp(filename,ios::in);
        buf<<fp.rdbuf();
        fp.close();
    }
    vector<string> lines;
    Split(buf.str(),lines,'\n',true); 
    buf.str(string());
    if (lines.size()<=1)
    {
        cerr<<"ERROR! Empty structure "<<infile<<endl;
        vector<string>().swap(lines);
        return 0;
    }

    /* parse PDB ID
     * HEADER, AUTHOR, JRNL, CRYST1, SCALEn */
    bool loop_=false;

    map<string,int> _citation;
    map<string,int> _atom_site;
    map<string,int> _struct_ref_seq;
    string _citation_title="";
    string _citation_pdbx_database_id_PubMed="";

    string group_PDB  =""; // (ATOM/HETATM)
    string type_symbol="";    // (element symbol)
    string atom_id    ="";   // auth_atom_id, label_atom_id (atom name)
    string comp_id    ="";  // auth_comp_id, label_comp_id (residue name)
    string asym_id    =""; // auth_asym_id, label_asym_id (chain ID)
    int    seq_id     =0; // label_seq_id, auth_seq_id (residue index)
    string pdbx_PDB_ins_code="";// (insertion code)
    double Cartn_x    =0; 
    double Cartn_y    =0; 
    double Cartn_z    =0; 
    double B_iso_or_equiv=0;   // Bfactor
    string alt_id     =" ";    // auth_alt_id, label_alt_id
    string pdbx_PDB_model_num=""; // model index

    string asym_prev    ="";
    string ins_code_prev="";
    int    seq_id_prev  =0;
    int    het          =0;

    vector<pair<string,string> > dbref_vec;
    string pdbx_db_accession;

    AtomUnit atom;
    atom.xyz.assign(3,0);
    atom.name="";

    ResidueUnit residue;
    residue.resn="";
    residue.resi=0;
    residue.icode=0;

    ChainUnit chain;
    chain.asym_id="";
    
    ModelUnit pep;


    size_t l;
    size_t a,c,r;
    int i,j;
    string line;
    vector<string> line_vec;
    vector<string> line_append_vec;
    for (l=0;l<lines.size();l++)
    {
        line=lines[l];
        if (_atom_site.size() && !StartsWith(line,"_atom_site"))
             Split(line,line_vec,' ',true);
        else Split(line,line_vec,' ');
        //cout<<"["<<l<<"] "<<line<<endl;
        if (line_vec.size()==0) continue;
        else if (line_vec.size()==1 && line_vec[0]=="#")
        {
            if (pdbx_db_accession.size() && pdbx_db_accession!="?")
            {
                dbref_vec.push_back(make_pair(asym_id,pdbx_db_accession));
                pdbx_db_accession.clear();
            }
            _citation.clear();
            _atom_site.clear();
            _struct_ref_seq.clear();
            loop_=false;
        }
        else if (line_vec.size()==1 && line_vec[0]=="loop_")
            loop_=true;
        else if (pdbid.size()==0 && l==0 && StartsWith(line,"data_"))
            pdbid=Lower(line.substr(5));
        else if (pdbid.size()==0 && line_vec.size()>1 && line_vec[0]=="_entry.id")
            pdbid=Lower(line_vec[1]);
        else if (StartsWith(line,"_struct_ref_seq."))
        {
            if (loop_)
            {
                j=_struct_ref_seq.size();
                line=line_vec[0];
                clear_line_vec(line_vec);
                Split(line,line_vec,'.');
                if (line_vec.size()>1)
                {
                    line=line_vec[1];
                    _struct_ref_seq[line]=j;
                }
            }
            else
            {
                l=read_semi_colon(line_vec, 2, l, lines, line_append_vec, line);
                if (line_vec[0]=="_struct_ref_seq.pdbx_strand_id")
                    asym_id=line_vec[1];
                else if (line_vec[0]=="_struct_ref_seq.pdbx_db_accession")
                    pdbx_db_accession=line_vec[1];
            }
        }
        else if (_struct_ref_seq.count("pdbx_strand_id") &&
            _struct_ref_seq.count("pdbx_db_accession"))
        {
            l=read_semi_colon(line_vec, _struct_ref_seq.size(),
                l, lines, line_append_vec, line);
            asym_id=line_vec[_struct_ref_seq["pdbx_strand_id"]];
            pdbx_db_accession=line_vec[_struct_ref_seq["pdbx_db_accession"]];
            if (pdbx_db_accession!="?")
                dbref_vec.push_back(make_pair(asym_id,pdbx_db_accession));
        }
        else if (StartsWith(line,"_citation."))
        {
            if (loop_)
            {
                line=line_vec[0];
                clear_line_vec(line_vec);
                Split(line,line_vec,'.');
                if (line_vec.size()>1)
                {
                    j=_citation.size();
                    line=line_vec[1];
                    _citation[line]=j;
                }
            }
            else
            {
                l=read_semi_colon(line_vec, 2, l, lines, line_append_vec, line,false,true);
                if      (line_vec[0]=="_citation.title")
                    _citation_title=Trim(line_vec[1],"'\"");
                else if (line_vec[0]=="_citation.pdbx_database_id_PubMed")
                    _citation_pdbx_database_id_PubMed=line_vec[1];
            }
        }
        else if (_citation.size() && (_citation.count("title")==0 ||
                 line_vec[_citation["id"]]=="pdbx_database_id_PubMed"))
        {
            l=read_semi_colon(line_vec, _citation.size(),
                l, lines, line_append_vec, line,false,true);
            if (_citation.count("title"))
                _citation_title=Trim(line_vec[_citation["title"]],"'\"");
            if (_citation.count("pdbx_database_id_PubMed"))
                _citation_pdbx_database_id_PubMed=line_vec[
                _citation["pdbx_database_id_PubMed"]];
        }
        else if (StartsWith(line,"_atom_site."))
        {
            line=line_vec[0];
            clear_line_vec(line_vec);
            Split(line,line_vec,'.');
            if (line_vec.size()>1)
            {
                line=line_vec[1];
                j=_atom_site.size();
                _atom_site[line]=j;
            }
        }
        else if (_atom_site.size() && _atom_site.count("Cartn_x") &&
            _atom_site.count("Cartn_y") && _atom_site.count("Cartn_z"))
        {
            if (_atom_site.count("pdbx_PDB_model_num"))
            {
                pdbx_PDB_model_num=line_vec[_atom_site["pdbx_PDB_model_num"]];
                if (pdbx_PDB_model_num!="." && pdbx_PDB_model_num!="?" &&
                    pdbx_PDB_model_num!="1")
                {
                    clear_line_vec(line_vec);
                    continue;
                }
            }
            
            if (_atom_site.count("label_comp_id"))
                comp_id=line_vec[_atom_site["label_comp_id"]];
            else if (_atom_site.count("auth_comp_id"))
                comp_id=line_vec[_atom_site["auth_comp_id"]];
            if (comp_id=="HOH")
            {
                clear_line_vec(line_vec);
                continue;
            }

            if (_atom_site.count("label_atom_id"))
                atom_id=line_vec[_atom_site["label_atom_id"]];
            else if (_atom_site.count("auth_atom_id"))
                atom_id=line_vec[_atom_site["auth_atom_id"]];
            if ((atom_id[0]=='"'  && atom_id[atom_id.size()-1]=='"')||
                (atom_id[0]=='\'' && atom_id[atom_id.size()-1]=='\''))
                atom_id=atom_id.substr(1,atom_id.size()-2);
            atom_id=atom_id.substr(0,4);
            
            if (_atom_site.count("type_symbol"))
                type_symbol=line_vec[_atom_site["type_symbol"]].substr(0,2);
            else type_symbol=lstrip(atom_id,"1234567890 ")[0];
            if (type_symbol=="H")
            {
                clear_line_vec(line_vec);
                continue;
            }

            if (type_symbol.size()==2) while (atom_id.size()<4) atom_id+=' ';
            else
            {
                if (atom_id.size()==1) atom_id+=' ';
                if (atom_id.size()==2) atom_id+=' ';
                if (atom_id.size()==3) atom_id=' '+atom_id;
            }
            if (type_symbol.size()==1) type_symbol=' '+type_symbol;
            
            if (_atom_site.count("group_PDB"))
                group_PDB=line_vec[_atom_site["group_PDB"]];

            if (_atom_site.count("auth_asym_id"))
                asym_id=line_vec[_atom_site["auth_asym_id"]];
            else if (_atom_site.count("label_asym_id"))
                asym_id=line_vec[_atom_site["label_asym_id"]];
            if (asym_id=="." || asym_id=="?") asym_id=" ";

            if (_atom_site.count("auth_seq_id"))
                seq_id=atoi(line_vec[_atom_site["auth_seq_id"]].c_str());
            else if (_atom_site.count("label_seq_id"))
                seq_id=atoi(line_vec[_atom_site["label_seq_id"]].c_str());

            if (_atom_site.count("pdbx_PDB_ins_code"))
                pdbx_PDB_ins_code=line_vec[_atom_site["pdbx_PDB_ins_code"]];
            if (pdbx_PDB_ins_code=="." || pdbx_PDB_ins_code=="?")
                pdbx_PDB_ins_code=" ";
            else pdbx_PDB_ins_code=pdbx_PDB_ins_code[0];

            Cartn_x=ReadDouble(line_vec[_atom_site["Cartn_x"]],3);
            Cartn_y=ReadDouble(line_vec[_atom_site["Cartn_y"]],3);
            Cartn_z=ReadDouble(line_vec[_atom_site["Cartn_z"]],3);

            if (_atom_site.count("B_iso_or_equiv"))
                B_iso_or_equiv=ReadDouble(line_vec[_atom_site["B_iso_or_equiv"]],2);

            if (_atom_site.count("auth_alt_id"))
                alt_id=line_vec[_atom_site["auth_alt_id"]];
            else if (_atom_site.count("label_alt_id"))
                alt_id=line_vec[_atom_site["label_alt_id"]];
            else if (_atom_site.count("pdbx_auth_alt_id"))
                alt_id=line_vec[_atom_site["pdbx_auth_alt_id"]];
            else if (_atom_site.count("pdbx_label_alt_id"))
                alt_id=line_vec[_atom_site["pdbx_label_alt_id"]];
            if (alt_id=="." || alt_id=="?") alt_id=" ";
            else if (asym_prev==asym_id && seq_id_prev==seq_id &&
                ins_code_prev==pdbx_PDB_ins_code)
            {
                bool found=(comp_id!=residue.resn);
                for (a=0;a<residue.atoms.size();a++)
                {
                    if (residue.atoms[a].name==atom_id)
                    {
                        found=true;
                        break;
                    }
                }
                if (found)
                {
                    clear_line_vec(line_vec);
                    continue;
                }
            }
            
            if (asym_id!=asym_prev || seq_id!=seq_id_prev || 
                pdbx_PDB_ins_code!=ins_code_prev)
            {
                if (residue.atoms.size())
                {
                    chain.residues.push_back(residue);
                    deepClean(residue);
                }
                if (asym_id!=asym_prev)
                {
                    //cout<<infile<<":"<<asym_id<<endl;
                    if (chain.residues.size())
                    {
                        pep.chains.push_back(chain);
                        deepClean(chain);
                    }
                    chain.asym_id=asym_id;
                }
                residue.resi=seq_id;
                residue.icode=pdbx_PDB_ins_code[0];
                residue.resn=comp_id;
                if (group_PDB=="ATOM") residue.het=1;
                else if (group_PDB=="HETATM")
                {
                    residue.het=3;
                    if (_atom_site.count("label_seq_id") && 
                        line_vec[_atom_site["label_seq_id"]]!=".")
                        residue.het=2;
                }
            }

            atom.name=atom_id;
            atom.element=type_symbol;
            atom.xyz[0]=Cartn_x;
            atom.xyz[1]=Cartn_y;
            atom.xyz[2]=Cartn_z;
            atom.bfactor=B_iso_or_equiv;
            
            residue.atoms.push_back(atom);

            asym_prev=asym_id;
            seq_id_prev=seq_id;
            ins_code_prev=pdbx_PDB_ins_code;
        }

        /* clean up */
        clear_line_vec(line_vec);
        lines[l].clear();
    }
    lines.clear();

    deepClean(atom);
    if (residue.atoms.size())
    {
        chain.residues.push_back(residue);
        deepClean(residue);
    }
    if (chain.residues.size())
    {
        pep.chains.push_back(chain);
        deepClean(chain);
    }
                
    if (pdbid.size()==0)
    {
        cerr<<"ERROR: no PDB ID in "<<infile<<endl;
        return -1;
    }

    /* consolidate the list of chains */
    vector<string> chainID_vec;
    for (c=0;c<pep.chains.size();c++)
    {
        asym_id=pep.chains[c].asym_id;
        if (find(chainID_vec.begin(), chainID_vec.end(), asym_id) != 
            chainID_vec.end()) continue;
        chainID_vec.push_back(asym_id);
    }
    /* TODO: map modified residues to standard residue */
    /* TODO: check molecule type */

    /* write macromolecules */
    ofstream fout;
    string filename;
    int serial;
    string resSeq;
    stringstream fout_buf;
    size_t c0;
    for (c0=0;c0<chainID_vec.size();c0++)
    {
        asym_id=chainID_vec[c0];
        serial=0;
        for (c=0;c<pep.chains.size();c++)
        {
            if (asym_id!=pep.chains[c].asym_id) continue;
        /*
COLUMNS        DATA  TYPE    FIELD        DEFINITION
-------------------------------------------------------------------------------------
 1 -  6        Record name   "ATOM  "
 7 - 11        Integer       serial       Atom  serial number.
13 - 16        Atom          name         Atom name.
17             Character     altLoc       Alternate location indicator.
18 - 20        Residue name  resName      Residue name.
22             Character     chainID      Chain identifier.
23 - 26        Integer       resSeq       Residue sequence number.
27             AChar         iCode        Code for insertion of residues.
31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
55 - 60        Real(6.2)     occupancy    Occupancy.
61 - 66        Real(6.2)     tempFactor   Temperature  factor.
77 - 78        LString(2)    element      Element symbol, right-justified.
79 - 80        LString(2)    charge       Charge  on the atom.
         */
            for (r=0;r<pep.chains[c].residues.size();r++)
            {
                if (pep.chains[c].residues[r].het>=3) continue;
                buf<<setw(4)<<pep.chains[c].residues[r].resi
                    <<pep.chains[c].residues[r].icode;
                resSeq=buf.str().substr(0,5);
                buf.str(string());
                for (a=0;a<pep.chains[c].residues[r].atoms.size();a++)
                {
                    serial++;
                    fout_buf<<"ATOM  "<<right<<setw(5)<<serial%100000<<' '
                        <<pep.chains[c].residues[r].atoms[a].name<<' '
                        <<pep.chains[c].residues[r].resn.substr(0,3)
                        <<' '<<pep.chains[c].asym_id[0]<<resSeq<<"   "
                        <<writeDouble(pep.chains[c].residues[r].atoms[a].xyz[0],8,3)
                        <<writeDouble(pep.chains[c].residues[r].atoms[a].xyz[1],8,3)
                        <<writeDouble(pep.chains[c].residues[r].atoms[a].xyz[2],8,3)
                        <<"  1.00"
                        <<writeDouble(pep.chains[c].residues[r].atoms[a].bfactor,6,2)
                        <<"          "
                        <<setw(2)<<pep.chains[c].residues[r].atoms[a].element<<'\n';
                }
            }
        }
        if (serial)
        {
            fout_buf<<"TER"<<endl;
            filename=pdbid+asym_id+".pdb";
            cout<<filename<<endl;
            fout.open(filename.c_str());
            fout<<fout_buf.str()<<flush;
            fout.close();
        }
        fout_buf.str(string());
    }


    /* clean up */
    filename.clear();
    _citation_title.clear();
    _citation_pdbx_database_id_PubMed.clear();
    group_PDB.clear();
    type_symbol.clear();
    atom_id.clear();
    comp_id.clear();
    asym_id.clear();
    pdbx_PDB_ins_code.clear();
    alt_id.clear();
    pdbx_PDB_model_num.clear();
    asym_prev.clear();
    ins_code_prev.clear();
    vector<pair<string,string> >().swap(dbref_vec);
    pdbx_db_accession.clear();
    deepClean(pep);
    line.clear();
    vector<string>().swap(lines);
    vector<string>().swap(chainID_vec);
    return 0;
}

int main(int argc,char **argv)
{
    string infile ="";
    string pdbid  ="";

    for (int a=1;a<argc;a++)
    {
        if (infile.size()==0)
            infile=argv[a];
        else if (pdbid.size()==0)
            pdbid=argv[a];
        else
        {
            cerr<<"ERROR: unknown option "<<argv[a]<<endl;
            return 1;
        }
    }

    if (infile.size()==0)
    {
        cerr<<docstring;
        return 1;
    }

    cif2pdb(infile,pdbid);

    /* clean up */
    string ().swap(infile);
    string ().swap(pdbid);
    return 0;
}

/* main END */
