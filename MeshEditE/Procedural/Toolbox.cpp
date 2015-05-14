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


using namespace HMesh;


namespace Procedural {
    
    Toolbox::Toolbox(){}
    
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
        rapidjson::Value &toolbox = d["toolbox"];
        
        
        
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
    }
    
    void Toolbox::updateProbabilities(){
        float tot_pieces_f = static_cast<float>(total_pieces);
        
        for( size_t i = 0; i < modules.size(); ++i){
            float no_pieces_f = static_cast<float>( modules[i].no_pieces );
            modules[i].probability = no_pieces_f / tot_pieces_f;
        }
    }
    
    const Module& Toolbox::getNext(){
        
        // a lot of things to do here!
        int selected = 0;
        
        modules[0].no_pieces -= 1;
        total_pieces -= 1;
        return *(modules[0].m);
    }
    
    
    
}