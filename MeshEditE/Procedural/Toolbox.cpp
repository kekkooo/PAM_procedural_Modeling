//
//  Toolbox.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 14/05/15.
//  Copyright (c) 2015 J. Andreas Bærentzen. All rights reserved.
//

#include "Toolbox.h"

#include <string>
#include <fstream>
#include <streambuf>
#include <iostream>


#include <GEL/HMesh/obj_load.h>

#include "rapidjson/document.h"


using namespace std;
using namespace HMesh;


namespace Procedural {
    
    Toolbox::Toolbox(){
        unsigned seed   = chrono::system_clock::now().time_since_epoch().count();
        randomizer.seed( seed );
        rand_max =  static_cast<float>( randomizer.max( ));
    }
    
    Toolbox& Toolbox::getToolboxInstance(){
        static Toolbox instance;
        return instance;
    }
    
    /// Loads a toolbox from a JSON configuration file
    void Toolbox::fromJson( std::string path ){
        
        std::cout << path;
        std::ifstream t( path );
        assert( t.good() );
        std::string json( (std::istreambuf_iterator<char>(t)),
                          std::istreambuf_iterator<char>());
        
        rapidjson::Document d;
        d.Parse( json.c_str() );
        rapidjson::Value &tb = d["toolbox"];
        assert( tb.IsArray() );
        
        for( rapidjson::SizeType i = 0; i < tb.Size(); ++i ){
            assert( tb[i].HasMember( "filename" ));
            assert( tb[i].HasMember( "config" ));
            assert( tb[i].HasMember( "type" ));
            assert( tb[i].HasMember( "probability" ));
            assert( tb[i].HasMember( "no_pieces" ));
            assert( tb[i].HasMember( "no_glueings" ));
            
            string mFilename    = tb[i]["filename"].GetString();
            string mConfig      = tb[i]["config"].GetString();
            Moduletype mType    = tb[i]["type"].GetInt();
            double mProbability = tb[i]["probability"].GetDouble();
            int    mNoPieces    = tb[i]["no_pieces"].GetInt();
            int    mNoGlueings  = tb[i]["no_glueings"].GetInt();
            
            // get the name of the module
            char delimiter = '/';
            int index = mFilename.length() - 1;
            while( mFilename[index] != delimiter ){ --index; assert( index >= 0 );  }
            string mName = mFilename.substr( index, mFilename.length() + 1 - index -3  );

            ModuleInfo* mInfo = new ModuleInfo;
            mInfo->m                    = new Module( mFilename, mType);
            mInfo->m->no_of_glueings    = mNoGlueings;
            mInfo->no_glueings          = mNoGlueings;
            mInfo->no_pieces            = mNoPieces;
            mInfo->probability          = mProbability;
            mInfo->name                 = mName;
            this->modules.push_back( mInfo );
            total_pieces += mNoPieces;
        }
    }
    
    void Toolbox::clear(){
        total_pieces = 0;
        modules.clear();
    }
    
    void Toolbox::addModule( std::string path, size_t no_pieces ){
        
        Manifold* manifold = new Manifold();
        obj_load( path, *manifold );
        // here I should load configuration
        // no_glueings
        // type
        size_t no_glueings = 1;
        Module* module = new Module( *manifold, 0, 1 );
        ModuleInfo mi;
        mi.m = module;
        mi.no_pieces = no_pieces;
        mi.no_glueings = no_glueings;

        total_pieces += no_pieces;
    }
    
    void Toolbox::updateProbabilities(){
        float tot_pieces_f = static_cast<float>(total_pieces);
        
        for( size_t i = 0; i < modules.size(); ++i){
            float no_pieces_f = static_cast<float>( modules[i]->no_pieces );
            modules[i]->probability = no_pieces_f / tot_pieces_f;
        }
    }
    
    Module& Toolbox::getNext(){
        
        assert( this->hasNext() );
        
        size_t index;
        bool done = false;
        
        // how do I choose which will be the next module?
#warning a lot of things to do here!

        while( !done ){
            // choose the index
            do{
                index = randomizer( ) % modules.size();
            }while( modules[index]->no_pieces <= 0 );
            
            size_t attempt = randomizer() % total_pieces;
            done = ( attempt <= modules[index]->no_pieces );
        }
        
        modules[index]->no_pieces   -= 1;
        total_pieces                -= 1;
        last_used_module            = index;
        used_module                 = true;
        
        cout << " picking module :  " << modules[index]->name << endl;
        
        return *(modules[index]->m);
    }
    
    bool Toolbox::hasNext() const {
#warning a lot of things to do here!
        return ( total_pieces > 0 );
    }
    
    void Toolbox::undoLast(){
        assert( used_module );
        assert(modules.size() > last_used_module );
        
        modules[last_used_module]->no_pieces += 1;
        total_pieces                         += 1;
        used_module                          = false;
    }
    
    void Toolbox::print() const{
        for ( const auto& m : modules ) {
            cout << " ##### MODULE : " << m->name <<  "######" << endl
                 << m->no_glueings << " glueings " << endl
                 << m->no_pieces   << " available" << endl << endl;
        }
    }
    
    
    
}