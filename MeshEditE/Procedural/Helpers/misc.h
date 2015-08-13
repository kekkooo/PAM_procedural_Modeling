//
//  misc.h
//  MeshEditE
//
//  Created by Francesco Usai on 13/08/15.
//  Copyright (c) 2015 J. Andreas Bærentzen. All rights reserved.
//

#ifndef MeshEditE_misc_h
#define MeshEditE_misc_h

#include <string>

#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <sys/stat.h>


namespace Procedural{
    namespace Helpers{
        namespace Misc{
            
inline std::string get_filename_stem( std::string full_path ){
//    std::string full_path( "/Users/francescousai/Documents/Dottorato/Conferenze/CGI_PG2015/Shapes/Modules/1plus_90_sim2.obj" );

    std::string::reverse_iterator revit = full_path.rbegin();
    size_t point_pos = full_path.length() - 1;
    size_t slash_pos = point_pos;
    while (( revit != full_path.rend()) && (*revit != '/' ) ) {
        
        std::cout << *revit;
        if( *revit == '.' ){ point_pos = slash_pos; }
        
        --slash_pos;
        ++revit;
    }
    assert( slash_pos < point_pos );
    std::string file_stem = full_path.substr( slash_pos + 1, point_pos - slash_pos - 1 );
    
//    std::cout << "last slash at " << slash_pos << " ### point at " << point_pos << std::endl
//    << "filename stem is : " << file_stem << std::endl;

    return file_stem;
}
            
inline void new_folder( std::string base_path, std::string folder_name ){
    std::string full_path = base_path + "/" + folder_name;
    mkdir( full_path.c_str(), 0775 );
}
            
}}}

#endif
