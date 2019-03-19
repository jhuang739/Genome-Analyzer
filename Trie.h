#ifndef TRIE_INCLUDED
#define TRIE_INCLUDED

#include <string>
#include <vector>

template<typename ValueType>
class Trie
{
public:
    Trie();
    ~Trie();
    void reset();
    void insert(const std::string& key, const ValueType& value);
    std::vector<ValueType> find(const std::string& key, bool exactMatchOnly) const;
    
    // C++11 syntax for preventing copying and assignment
    Trie(const Trie&) = delete;
    Trie& operator=(const Trie&) = delete;
private:
    struct Node                 // node struct
    {
        char label;
        std::vector<ValueType> values;
        std::vector<Node*> children;
    };

    Node* root;
    
    // helper to delete nodes
    void freeNodes(Node* current);
    
    // helper to find nodes
    void findNodes(const std::string& key, bool exactMatchOnly, Node* current, std::vector<ValueType>& result) const;
};

// Implementations

template<typename ValueType>    // constructor
Trie<ValueType>::Trie()
{
    root = new Node();
}

template<typename ValueType>    // destructor
Trie<ValueType>::~Trie()
{
    freeNodes(root);
}

template<typename ValueType>
void Trie<ValueType>::freeNodes(Node* current)    // helper for destruction
{
    if(current->children.size() == 0)   // finished
        return;
    
    for(int i = 0; i < current->children.size(); i++)
    {
        if(current->children[i] != nullptr)
            freeNodes(current->children[i]);    // postorder traversal
        delete current->children[i];
    }
}

template<typename ValueType>
void Trie<ValueType>::reset()       // reset implementation
{
    freeNodes(root);
    root = new Node();
}

template<typename ValueType>
void Trie<ValueType>::insert(const std::string& key, const ValueType& value)        // insert implementation
{
    Node* curr = root;
    bool found = true;
    
    for(int i = 0; i < key.length(); i++)           // iterate through each char in key
    {
        found = false;
        for(int j = 0; j < curr->children.size(); j++)           // iterate through children
        {
            if(curr->children[j]->label == key[i])
            {
                curr = curr->children[j];
                found = true;
                break;
            }
        }
        
        if(!found)          // no existing node in trie with char label of key[i]
        {
            Node* newNode = new Node();             // insert node into children
            newNode->label = key[i];
            
            curr->children.push_back(newNode);
            curr = curr->children[curr->children.size()-1];
        }
    }
    
    curr->values.push_back(value);          // add value to node
}

template<typename ValueType>
std::vector<ValueType> Trie<ValueType>::find(const std::string& key, bool exactMatchOnly) const
{
    std::vector<ValueType> result;
    findNodes(key, exactMatchOnly, root, result);       // call recursive helper function
    return result;
}

template<typename ValueType>
void Trie<ValueType>::findNodes(const std::string& key, bool exactMatchOnly, Node* current, std::vector<ValueType>& result) const              // recursive helper function for find
{
    // base case pseudocode:
    // if key length is zero, push values of node into result
    
    if(key.length() == 0)
    {
        for(int i = 0; i < current->values.size(); i++)
            result.push_back(current->values[i]);
        return;
    }
    
    // simplifying step pseudocode:
    // if first character of key matches a label in children, traverse into child
    // if exactMatchOnly is false and not at root node, traverse into child and set exactMatchOnly to true
    
    for(int i = 0; i < current->children.size(); i++)
    {
        if(key[0] == current->children[i]->label)
            findNodes(key.substr(1, key.length()-1), exactMatchOnly, current->children[i], result);
        else if(!exactMatchOnly && current != root)
            findNodes(key.substr(1, key.length()-1), true, current->children[i], result);
    }
}

#endif // TRIE_INCLUDED


