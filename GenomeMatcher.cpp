#include "provided.h"
#include "Trie.h"
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <fstream>
using namespace std;

class GenomeMatcherImpl
{
public:
    GenomeMatcherImpl(int minSearchLength);
    void addGenome(const Genome& genome);
    int minimumSearchLength() const;
    bool findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const;
    bool findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const;
private:
    int minSearchLen;
    Trie<Genome> genomes;           // store genomes
    Trie<DNAMatch> trie;            // store fragments of minSearchLen mapping to DNAMatch
    
    static bool compare(const GenomeMatch& a, const GenomeMatch& b);
};

GenomeMatcherImpl::GenomeMatcherImpl(int minSearchLength)
: trie(), genomes(), minSearchLen(minSearchLength)
{}

void GenomeMatcherImpl::addGenome(const Genome& genome)
{
    genomes.insert(genome.name(), genome);      // store genome sequence in genomes trie
    
    for(int i = 0; i <= genome.length()-minSearchLen; i++)
    {
        string fragment;
        genome.extract(i, minSearchLen, fragment);  // find fragment of DNA in genome
        DNAMatch value;
        value.genomeName = genome.name();
        value.length = genome.length();
        value.position = i;
        trie.insert(fragment, value);               // insert into search trie
    }
}

int GenomeMatcherImpl::minimumSearchLength() const
{
    return minSearchLen;
}

bool GenomeMatcherImpl::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
    matches.clear();
    
    if(fragment.length() < minimumLength)
        return false;
    if(minimumLength < minimumSearchLength())
        return false;
    
    // find hits across genomes, prefix/snip of minSearchLen of fragment
    vector<DNAMatch> trieSearch = trie.find(fragment.substr(0, minimumSearchLength()), exactMatchOnly);
    unordered_map<string, DNAMatch> recorder;       // record info of genome matching fragments
    
    for(int i = 0; i < trieSearch.size(); i++)      // go through each distinct hit across genomes
    {
        Genome currentGenome = genomes.find(trieSearch[i].genomeName, true)[0]; // find exact genome
        string comp = "";
        int pos = trieSearch[i].position;
        int len = 0;
        currentGenome.extract(pos, (int)fragment.length(), comp);       // extract complete fragment from hit
        bool snip = false;
        bool push = false;
        
        for(int j = 0; j < fragment.length(); j++)  // go through each char in fragment
        {
            if(fragment[j] == comp[j])     // triesearch[i] genome at pos+j
            {
                len++;
                if(j == minimumLength-1)    // if fragments match, push to result
                    push = true;
                continue;
            }
            else if(!snip && !exactMatchOnly)
            {
                len++;
                if(j == minimumLength-1)    // if snip and fragment match, push to result
                    push = true;
                snip = true;
            }
            else
                break;
        }
        
        trieSearch[i].length = len;
        
        if(recorder.find(trieSearch[i].genomeName) != recorder.end())
        {
            // only record fragment of greatest length, earliest position tiebreaker
            if(trieSearch[i].length > recorder[trieSearch[i].genomeName].length || ((trieSearch[i].length == recorder[trieSearch[i].genomeName].length) && trieSearch[i].position < recorder[trieSearch[i].genomeName].position))
                recorder.erase(trieSearch[i].genomeName);
            else
                continue;
        }

        if(push)        // record match
        {
            pair<string, DNAMatch> inserted(trieSearch[i].genomeName, trieSearch[i]);
            recorder.insert(inserted);
        }
    }
    
    // add to matches result information
    for(unordered_map<string, DNAMatch>::iterator it = recorder.begin(); it != recorder.end(); it++)
        matches.push_back(it->second);
    
    if(matches.empty())     // no matches found
        return false;
    return true;
}

bool GenomeMatcherImpl::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
    results.clear();
    
    if(fragmentMatchLength < minimumSearchLength())
        return false;
    
    int numQueries = query.length()/fragmentMatchLength;
    unordered_map<string, int> recorder;        // records number of matches for each genome
    vector<DNAMatch> matches;                   // store sequence matches
    
    for(int i = 0; i < numQueries; i++)
    {
        matches.clear();
        int index = i * fragmentMatchLength;    // fragments start index increment by fragmentMatchLength
        string sequence = "";
        query.extract(index, fragmentMatchLength, sequence);    // extract sequence from query
        findGenomesWithThisDNA(sequence, fragmentMatchLength, exactMatchOnly, matches); // find matches
        for(int j = 0; j < matches.size(); j++)
        {
            if(recorder.find(matches[j].genomeName) == recorder.end())      // if not initialized in map
            {
                pair<string, int> inserted(matches[j].genomeName, 0);
                recorder.insert(inserted);
            }
            recorder[matches[j].genomeName]++;      // increment number of matching sequences
        }
    }
    
    // for each genome in GenomeMatcherImpl, see if percentage match meets threshold
    for(unordered_map<string, int>::iterator it = recorder.begin(); it != recorder.end(); it++)
    {
        double percentage = (double)(it->second)/numQueries;
        percentage *= 100;
        if(percentage >= matchPercentThreshold)
        {
            GenomeMatch toBeInserted;               // if meets threshold, push to result vector
            toBeInserted.genomeName = it->first;
            toBeInserted.percentMatch = percentage;
            results.push_back(toBeInserted);
        }
    }
    
    // sort results using <algorithm> sort
    sort(results.begin(), results.end(), &GenomeMatcherImpl::compare);
    
    if(results.empty())     // no related genomes
        return false;
    return true;
}

// compare function for sort
bool GenomeMatcherImpl::compare(const GenomeMatch& a, const GenomeMatch& b)
{
    if(a.percentMatch == b.percentMatch)            // tiebreaker by name
    {
        if(a.genomeName > b.genomeName)
            return false;
        else
            return true;
    }
    return a.percentMatch > b.percentMatch;         // descending order
}

//******************** GenomeMatcher functions ********************************

// These functions simply delegate to GenomeMatcherImpl's functions.
// You probably don't want to change any of this code.

GenomeMatcher::GenomeMatcher(int minSearchLength)
{
    m_impl = new GenomeMatcherImpl(minSearchLength);
}

GenomeMatcher::~GenomeMatcher()
{
    delete m_impl;
}

void GenomeMatcher::addGenome(const Genome& genome)
{
    m_impl->addGenome(genome);
}

int GenomeMatcher::minimumSearchLength() const
{
    return m_impl->minimumSearchLength();
}

bool GenomeMatcher::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
    return m_impl->findGenomesWithThisDNA(fragment, minimumLength, exactMatchOnly, matches);
}

bool GenomeMatcher::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
    return m_impl->findRelatedGenomes(query, fragmentMatchLength, exactMatchOnly, matchPercentThreshold, results);
}
