# this works for gcc
# compatibility for other compilers should go here
function( target_link_libraries_whole_archive target library )
    target_link_libraries( ${target} "-Wl,--whole-archive" ${library} "-Wl,--no-whole-archive" )
endfunction( target_link_libraries_whole_archive )