// Copyright (c) 2005, Google Inc.
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following disclaimer
// in the documentation and/or other materials provided with the
// distribution.
//     * Neither the name of Google Inc. nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// ---
// Author: Craig Silverstein

#include <stdio.h>
#include <time.h>              // for silly random-number-seed generator
#include <math.h>              // for sqrt()
#include <map>
#include <set>
#include <iterator>            // for insert_iterator
#include <iostream>
#include <iomanip>             // for setprecision()
#include <string>


#include <mm/cache_map.hpp>
#include <mm/cache_set.hpp>
#include <mm/hash_fun.hpp>

using mm::cache_map;
using mm::cache_set;
using mm::hash;

using std::map;
using std::pair;
using std::string;
using std::insert_iterator;
using std::allocator;
using std::equal_to;
using std::ostream;

static const int DEFAULT_SIZE = 100000;

#define CHECK(cond)  do {                                  \
  if (!(cond)) {                                           \
    std::cout << std::endl;                                \
    std::cout << ">>>> at line " << __LINE__ << std::endl; \
    std::cout << ">>>> Test failed: " #cond << std::endl;  \
    std::exit( 1 );                                        \
  }                                                        \
} while (0)

const char *words[] = { "Baffin",        // in 'words'
                        "Boffin",        // not in
                        "baffin",        // not in
                        "genial",        // last word in
                        // "Aarhus",        // first word alphabetically
                        "Zurich",        // last word alphabetically
};

const char *nwords[] = { "Boffin",
                         "baffin"
                       };

// Apparently identity is not stl-standard, so we define our own
template<class Value>
struct Identity {
    Value& operator()(Value& v) const { return v; }
    const Value& operator()(const Value& v) const { return v; }
};

// Let us log the pairs that make up a hash_map
template<class P1, class P2>
ostream& operator<<(ostream& s, const pair<P1, P2>& p) {
    s << "pair(" << p.first << ", " << p.second << ")";
    return s;
}

// Comparison function for 'char *'
struct strcmp_fnc
{
    bool operator()(const char* s1, const char* s2) const {
        return ((s1 == 0 && s2 == 0) ||
                (s1 && s2 && *s1 == *s2 && strcmp(s1, s2) == 0));
    }
};


template <class K, class V, class H, class C, class DF>
static void set_empty_key( cache_map<K,V,H,C,DF> *ht, K key )
{
    ht->set_empty_key( key );
}

template <class V, class H, class C, class DF>
static void set_empty_key( cache_set<V,H,C,DF> *ht, V val )
{
    ht->set_empty_key( val );
}

template <class K, class V, class H, class C, class DF>
static void insert( cache_map<K,V,H,C,DF> *ht, K val )
{
    ht->insert( pair<K,V>( val, V() ) );
}

template <class V, class H, class C, class DF>
static void insert( cache_set<V,H,C,DF> *ht, V val )
{
    ht->insert( val );
}

template <class HT, class Iterator>
static void insert( HT *ht, Iterator begin, Iterator end )
{
    ht->insert(begin, end);
}

// For hashtable's and hash_set's, the iterator insert works fine (and
// is used). But for the hash_map's, the iterator insert expects the
// iterators to point to pair's. So by looping over and calling insert
// on each element individually, the code below automatically expands
// into inserting a pair.
template <class K, class V, class H, class C, class DF, class Iterator>
static void insert(cache_map<K,V,H,C,DF> *ht, Iterator begin, Iterator end)
{
    while (begin != end) {
        insert(ht, *begin);
        ++begin;
    }
}

template <class V, class H, class C, class DF, class Iterator>
static void insert(cache_set<V,H,C,DF> *ht, Iterator begin, Iterator end)
{
    while ( begin != end ) {
        insert( ht, *begin );
        ++begin;
    }
}

// A version of insert that uses the insert_iterator.  But insert_iterator
// isn't defined for the low level hashtable classes, so we just punt to insert.

template <class K, class V, class H, class C, class DF>
static void iterator_insert( cache_map<K,V,H,C,DF>* , K key,
                             insert_iterator<cache_map<K,V,H,C> >* ii)
{
    *(*ii)++ = pair<K,V>( key, V() );
}

template <class V, class H, class C, class DF>
static void iterator_insert( cache_set<V,H,C,DF>* , V val,
                             insert_iterator<cache_set<V,H,C> >* ii)
{
    *(*ii)++ = val;
}

static char* read_line(FILE* fp, char* line, int linesize)
{
    if ( fgets(line, linesize, fp) == NULL )
        return NULL;
    // normalize windows files :-(
    const int linelen = strlen(line);
    if ( linelen >= 2 && line[linelen-2] == '\r' && line[linelen-1] == '\n' ) {
        line[linelen-2] = '\n';
        line[linelen-1] = '\0';
    }
    return line;
}

// Performs tests where the hashtable's value type is assumed to be int.
template <class htint>
void test_int()
{
    htint y(1000);
    htint z(64);
    set_empty_key(&y, 0xefefef);
    set_empty_key(&z, 0xefefef);
    
    CHECK(y.empty());
    insert(&y, 1);
    CHECK(!y.empty());
    insert(&y, 11);
    insert(&y, 111);
    insert(&y, 1111);
    insert(&y, 11111);
    insert(&y, 111111);
    insert(&y, 1111111);     // 1M, more or less
    insert(&y, 11111111);
    insert(&y, 111111111);
    insert(&y, 1111111111);  // 1B, more or less

    for ( int i = 0; i < 32; ++i )
        insert(&z, i);
    // test the second half of the insert with an insert_iterator
    insert_iterator<htint> insert_iter(z, z.begin());
    for ( int i = 32; i < 64; ++i )
        iterator_insert(&z, i, &insert_iter);
    
    for ( typename htint::const_iterator it = y.begin(); it != y.end(); ++it )
        std::cout << "y: " << *it << "\n";
    z.insert(y.begin(), y.end());
    swap( y, z );
  
    for ( typename htint::iterator it = y.begin(); it != y.end(); ++it )
        std::cout << "y+z: " << *it << "\n";
    std::cout << "z has " << z.bucket_count() << " buckets\n";
    std::cout << "y has " << y.bucket_count() << " buckets\n";
    std::cout << "z size: " << z.size() << "\n";

    for (int i = 0; i < 64; ++i) {
        // This test fails because 'y' has 64 buckets and inserting
        // more integers will overwrite existing ones..
        ///// CHECK(y.find(i) != y.end());
    }

    CHECK(z.size() == 10);
    z.erase(11111);
    CHECK(z.size() == 9);
    insert(&z, 11111);                  // should retake deleted value
    CHECK(z.size() == 10);
    // Do the delete/insert again.  Last time we probably resized; this time no
    z.erase(11111);
    insert(&z, 11111);                  // should retake deleted value
    CHECK(z.size() == 10);
    
    z.erase(-11111);                    // shouldn't do anything
    CHECK(z.size() == 10);
    z.erase(1);
    CHECK(z.size() == 9);
    typename htint::iterator itdel = z.find(1111);
    
    z.erase(itdel);
    
    CHECK(z.size() == 8);
    
    itdel = z.find(2222);               // should be end()
    z.erase(itdel);                     // shouldn't do anything
    
    CHECK(z.size() == 8);
    for ( typename htint::const_iterator it = z.begin(); it != z.end(); ++it )
        std::cout << "y: " << *it << "\n";
    for ( typename htint::const_iterator it = z.begin(); it != z.end(); ++it )
        std::cout << "y: " << *it << "\n";
    std::cout << "That's " << z.size() << " elements\n";
    std::cout << "Distance " << mm::distance( z.begin(), z.end() ) << " elements\n";
    z.erase( z.begin(), z.end() );
    std::cout << "That's " << z.size() << " elements\n";
    CHECK(z.empty());

    y.clear();
    CHECK(y.empty());
    std::cout << "y has " << y.bucket_count() << " buckets\n";
}

// Performs tests where the hashtable's value type is assumed to be char*.
// The read_write parameters specifies whether the read/write tests
// should be performed. Note that densehashtable::write_metdata is not
// implemented, so we only do the read/write tests for the
// sparsehashtable varieties.
template <class ht>
void test_charptr()
{
    ht w( DEFAULT_SIZE );
    set_empty_key(&w, (char*) NULL);
    insert(&w, const_cast<char **>(nwords),
           const_cast<char **>(nwords) + sizeof(nwords) / sizeof(*nwords));
    std::cout << "w has " << w.size() << " items\n";
    CHECK(w.size() == 2);
    CHECK(w == w);
    
    ht x( DEFAULT_SIZE );
    set_empty_key(&x, (char*) NULL);

    map<string, int> counts;
    // Hash the dictionary
    {
        const char* file = "words";
        FILE *fp = fopen( file, "rb" );
    
        assert( fp != NULL );
        
        char line[1024];
        while ( read_line( fp, line, sizeof(line) ) )
        {
            string s( line );
            s = s.substr( 0, s.length() -1  );
            
            insert( &x, strdup( s.c_str() ) );
            
            counts[ s ] = 0;
        }
        std::cout << "Read " << x.size() << " words from " << file << "\n";
        
        fclose(fp);
        
        for ( char **word = const_cast<char **>(words);
              word < ( const_cast<char **>(words)
                       + sizeof(words) / sizeof(*words) );
              ++word )
        {
            if ( x.find( *word ) == x.end() )
            {
                CHECK( w.find(*word) != w.end() );
            } else
            {
                CHECK( w.find(*word) == w.end() );
            }
        }
    }
    
    CHECK( counts.size() == ( x.size() + x.num_collisions() ) );
}

template <class ht>
void test_string()
{
    ht w( DEFAULT_SIZE );
    set_empty_key( &w, string("-*- empty key -*-") );
    const int N = sizeof(nwords) / sizeof(*nwords);
    string* nwords1 = new string[ N ];

    std::cout << "N: " << N << std::endl;
    
    for (int i = 0; i < N; ++i)
//        insert( &w, std::pair<string,int>( * nwords[i], i ) );
        nwords1[ i ] = nwords[ i ];
    
    insert( &w, nwords1, nwords1 + N );
    
    delete[] nwords1;
    
    std::cout << "w has " << w.size() << " items: " << std::endl;
    for ( typename ht::iterator it = w.begin(); it != w.end(); ++it )
    {
        std::cout << *it << std::endl;
    }
    
    
    CHECK( w.size() == 2 );
    CHECK( w == w );
    
    ht x( DEFAULT_SIZE );
    set_empty_key(&x, string("-*- empty key -*-"));

    std::set<string> counts;
    // Hash the dictionary
    {
        const char* file = "words";
        FILE *fp = fopen(file, "r");
        if ( fp == NULL )
        {
            std::cout << "Can't open " << file << ", skipping dictionary hash..." << std::endl;
        }
        else
        {
            char line[1024];
            while ( fgets(line, sizeof(line), fp) ) {
                string s( line );
                s = s.substr( 0, s.length() -1  );
                
                insert( &x, s );
                counts.insert( s );
            }
            fclose( fp );
            
            std::cout << "Read " << x.size() << " words from " << file << std::endl;
            std::cout << "Collisions: " << x.num_collisions() << std::endl;
            
            for ( char **word = const_cast<char **>(words);
                  word < ( const_cast<char **>(words)
                           + sizeof(words) / sizeof(*words) );
                  ++word )
            {
                std::cout << "checking: '" << *word << "':.. ";
                
                if ( x.find( *word ) == x.end() )
                {
                    std::cout << " NOT found in x " << std::endl;
                    
                    CHECK( w.find( *word ) != w.end() );
                }
                else
                {
                    std::cout << " found in x : "
                              << *( x.find( *word ) )
                              << std::endl;
                    if ( w.find( *word ) != w.end() )
                    {
                        std::cout << "w.size(): " << w.size() << std::endl;
                        
                        std::cout << "Error: found "
                                  << *( w.find( *word ) )
                                  << " in w" << std::endl;
                    }
                    
                    CHECK( w.find( *word ) == w.end() );
                }
            }
        }
    }

    std::cout << "counts.size(): " << counts.size() << std::endl;
    std::cout << "x.size(): " << x.size() << std::endl;
    std::cout << "x.max_size(): " << x.max_size() << std::endl;
    CHECK( counts.size() == ( x.size() + x.num_collisions() ) );  
}

// The read_write parameters specifies whether the read/write tests
// should be performed. Note that densehashtable::write_metdata is not
// implemented, so we only do the read/write tests for the
// sparsehashtable varieties.
template<class ht, class htstr, class htint>
void test()
{
    test_int<htint>();
    test_string<htstr>();
    test_charptr<ht>();
}

void print_bin( size_t n )
{
    const int hash_size = 29;
    
    for ( int i = hash_size-1; i >= 0; i--  )
        std::cout << ( ( n >> i ) & 0x0001ul );
}

void test_hash( const string& s1, const string& s2 )
{
    mm::hash<string> hasher;
    size_t h1 = hasher( s1 );
    size_t h2 = hasher( s2 );
    
    std::cout << "hash( " << s1 << " ): 0x" << std::hex << h1 << std::endl;
    std::cout << "hash( " << s2 << " ): 0x" << std::hex << h2 << std::endl;
    std::cout << "Diff: "; print_bin( h1 ^ h2 ); std::cout << std::endl;

    std::cout << std::endl;
}

int main( int argc, char **argv )
{
   /*
    std::cout << "\n\nTEST WITH CACHE_MAP\n\n";
    test < cache_map<char *, int, mm::hash<char *>, strcmp_fnc>,
           cache_map< string, int >,
           cache_map<int, int>
         >();
*/
    std::cout << "\n\nTEST WITH CACHE_SET\n\n";
    test < cache_set<char *, mm::hash<char *>, strcmp_fnc>,
           cache_set<string>,
           cache_set<int>
         >();

    std::cout << "\nAll tests pass.\n";

    std::cout << std::endl;
    
    return 0;
}
