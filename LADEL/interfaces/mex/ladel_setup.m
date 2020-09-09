function ladel_setup()
% Script to compile the mex-interface of LADEL

warning('off','MATLAB:mex:GccVersion_link')
% store the path from which this function is called
current_path = pwd;

% Change to directory containing this function
this_path = fileparts(mfilename('fullpath'));
cd(this_path);

details = 0 ;	    % 1 if details of each command are to be printed

try
    % ispc does not appear in MATLAB 5.3
    pc = ispc ;
    mac = ismac ;
catch                                                                       
    % if ispc fails, assume we are on a Windows PC if it's not unix
    pc = ~isunix ;
    mac = 0 ;
end

flags = '-silent -g' ;
is64 = ~isempty (strfind (computer, '64')) ;
if (is64)
    % 64-bit MATLAB
    flags = [flags ' -largeArrayDims'];
end

include = '-I. -I./include -I../../include -I../../amd/Include'; 

ladel_src = { ...
    'ladel_col_counts.c', ...
    'ladel_copy.c', ...
    'ladel_debug_print.c', ...
    'ladel_etree.c', ...
    'ladel_global.c', ...
    'ladel.c', ...
    'ladel_ldl_numeric.c', ...
    'ladel_ldl_symbolic.c', ...
    'ladel_matvec.c', ...
    'ladel_matmat.c', ...
    'ladel_submatrix.c', ...
    'ladel_add.c', ...
    'ladel_pattern.c', ...
    'ladel_permutation.c', ...
    'ladel_postorder.c', ...
    'ladel_rank1_mod.c', ...
    'ladel_row_mod.c', ...
    'ladel_scale.c', ...
    'ladel_transpose.c', ...
    'ladel_upper_diag.c', ...
    };

ladel_src_path = '../../src';
for i = 1:length(ladel_src)
    ladel_src{i} = [ladel_src_path '/' ladel_src{i}] ;
end

ladel_mex_util_src = ['ladel_mex_util.c']; 

amd_src = { ...
    'amd_1.c', ...
    'amd_2.c', ...
    'amd_aat.c', ...
    'amd_control.c', ...
    'amd_defaults.c', ...
    'amd_dump.c', ...
    'amd_global.c', ...
    'amd_info.c', ...
    'amd_order.c', ...
    'amd_postorder.c', ...
    'amd_post_tree.c', ...
    'amd_preprocess.c', ...
    'amd_valid.c', ...
    'SuiteSparse_config.c', ...
    };

amd_src_path = '../../amd/Source';
for i = 1:length(amd_src)
    amd_src{i} = [amd_src_path '/' amd_src{i}];
end

if (pc)
    obj_extension = '.obj' ;
else
    obj_extension = '.o' ;
end
obj = '' ;

cflags = ' CFLAGS="\$CFLAGS -std=c99 -fPIC -DMATLAB -DDAMD -DSIMPLE_COL_COUNTS -O3 -DDLONG"';
flags = [flags cflags];
source = [amd_src ladel_src ladel_mex_util_src];
kk = 0;
for f = source
    ff = f {1} ;
    slash = strfind (ff, '/') ;
    if (isempty (slash))
        slash = 1 ;
    else
        slash = slash (end) + 1 ;
    end
    o = ff (slash:end) ;
    
    % fprintf ('%s\n', o) ;
    o = [o(1:end-2) obj_extension] ;
    obj = [obj  ' ' o] ;					            %#ok
    s = sprintf ('mex %s %s -c %s', flags, include, ff) ;
    s = [s ' ' o] ;
    kk = do_cmd (s, kk, details) ;
end

s = sprintf ('mex %s %s %s.c', flags, include, 'ladel_mex') ;
s = [s obj] ;	
kk = do_cmd (s, kk, details) ;

 % clean up
s = ['delete ' obj] ;
do_cmd (s, kk, details) ;
fprintf ('\nLADEL successfully compiled\n') ;

end


function kk = do_cmd (s, kk, details)
 %DO_CMD: evaluate a command, and either print it or print a "."
if (details)
    fprintf ('%s\n', s) ;
else
    if (mod (kk, 60) == 0)
	fprintf ('\n') ;
    end
    kk = kk + 1 ;
    fprintf ('.') ;
end
eval (s) ;
end