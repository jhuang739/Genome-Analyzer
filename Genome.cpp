#include "provided.h"
#include <string>
#include <vector>
#include <iostream>
#include <istream>
#include <fstream>
using namespace std;

class GenomeImpl
{
public:
    GenomeImpl(const string& nm, const string& sequence);
    static bool load(istream& genomeSource, vector<Genome>& genomes);
    int length() const;
    string name() const;
    bool extract(int position, int length, string& fragment) const;
private:
    string m_name;
    string genome;
};

GenomeImpl::GenomeImpl(const string& nm, const string& sequence)
{
    m_name = nm;
    genome = sequence;
}

bool GenomeImpl::load(istream& genomeSource, vector<Genome>& genomes)
{
    if (!genomeSource)              // Did opening the file fail?
    {
        cerr << "Error: Cannot open data.txt!" << endl;
        return false;
    }
    
    genomes.clear();
    
    string s;                           // getline
    bool nameLine = false;              // is this a genome name line
    string name = "";                   // genome name
    string genome = "";                 // genome sequence

    while (getline(genomeSource, s))
    {
        for(int i = 0; i < s.size(); i++)   // for each char in string
        {
            char c = s[i];
            
            if(i == 0 && c == '>')          // check for name line indicator
            {
                nameLine = true;
                
                if(!name.empty() && !genome.empty())        // if prior genome data ended
                {
                    Genome current = Genome(name, genome);  // create new genome
                    genomes.push_back(current);
                    name = "";                              // reset genome name
                    genome = "";                            // reset genome sequence
                }
                else if(!name.empty())                      // consequetive name line check
                    return false;
                continue;
            }
            
            if(i == 0 && c != '>')                          // if first char isn't >, not name line
            {
                nameLine = false;
            }
            
            if(c == '\n')
                continue;
            
            if(nameLine)
            {
                name += c;
            }
            else
            {
                if(name.empty())        // cannot fill genome with no genome name
                    return false;
                if(c != 'A' && c != 'C' && c != 'T' && c != 'G' && c != 'N')
                    return false;       // check for valid genome character
                
                genome += c;
            }
        }
    }
    
    Genome current = Genome(name, genome);      // add last genome's information
    genomes.push_back(current);
    return true;
}

int GenomeImpl::length() const
{
    return (int)genome.size();
}

string GenomeImpl::name() const
{
    return m_name;
}

bool GenomeImpl::extract(int position, int length, string& fragment) const
{
    string local = "";              // temp container for desired extracted sequence
    for(int i = position; i < position + length; i++)
    {
        if(i >= genome.length())    // check for attempt to extract past end of genome
            return false;
        local += genome[i];
    }
    
    fragment = local;
    return true;
}

//******************** Genome functions ************************************

// These functions simply delegate to GenomeImpl's functions.
// You probably don't want to change any of this code.

Genome::Genome(const string& nm, const string& sequence)
{
    m_impl = new GenomeImpl(nm, sequence);
}

Genome::~Genome()
{
    delete m_impl;
}

Genome::Genome(const Genome& other)
{
    m_impl = new GenomeImpl(*other.m_impl);
}

Genome& Genome::operator=(const Genome& rhs)
{
    GenomeImpl* newImpl = new GenomeImpl(*rhs.m_impl);
    delete m_impl;
    m_impl = newImpl;
    return *this;
}

bool Genome::load(istream& genomeSource, vector<Genome>& genomes)
{
    return GenomeImpl::load(genomeSource, genomes);
}

int Genome::length() const
{
    return m_impl->length();
}

string Genome::name() const
{
    return m_impl->name();
}

bool Genome::extract(int position, int length, string& fragment) const
{
    return m_impl->extract(position, length, fragment);
}
