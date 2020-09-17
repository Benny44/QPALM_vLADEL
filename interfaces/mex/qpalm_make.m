function qpalm_make (solver)
%QPALM_MAKE compiles the QPALM mexFunction
%
% Example:
%   qpalm_make
%   qpalm_make('ladel')

%
% QPALM relies either on LADEL or on CHOLMOD, AMD, COLAMD and optionally CCOLAMD, CAMD, and METIS.
warning('off','MATLAB:mex:GccVersion_link')

close all

% store the path from which this function is called
current_path = pwd;

% Change to directory containing this function
this_path = fileparts(mfilename('fullpath'));
cd(this_path);

!git submodule update --init
!git submodule update --recursive --remote

%% For the matlab version of QPALM we need the following mex functions.
% PWAlinesearch_setup;

% cholmod_path = fullfile(this_path, '../../suitesparse/CHOLMOD/MATLAB');
% addpath(cholmod_path);
% cholmod_path = fileparts(which('cholmod_make.m'));
% cd(cholmod_path);
% cholmod_make;

% cd(this_path);

%% Determine the linear systems solver and specify the include paths
fprintf ('Compiling QPALM mex function\n')

details = 0 ;	    % 1 if details of each command are to be printed

v = version ;
try
    % ispc does not appear in MATLAB 5.3
    pc = ispc ;
    mac = ismac ;
catch                                                                       %#ok
    % if ispc fails, assume we are on a Windows PC if it's not unix
    pc = ~isunix ;
    mac = 0 ;
end

flags = '' ;
is64 = ~isempty (strfind (computer, '64')) ;
if (is64)
    % 64-bit MATLAB
    flags = '-largeArrayDims' ;
end

% MATLAB 8.3.0 now has a -silent option to keep 'mex' from burbling too much
if (~verLessThan ('matlab', '8.3.0'))
    flags = ['-silent ' flags] ;
end

if nargin < 1
    solver = 'ladel';


elseif strcmp(solver, 'ladel')
    fprintf('Using LADEL as the linear systems solver.\n');
elseif strcmp(solver, 'cholmod');
    error('CHOLMOD is not available in this version of QPALM\n');
else
    error('Unrecognised linear systems solver. Use ''ladel'' \n');
end

if strcmp(solver, 'ladel')
    include = '-I. -I../../include -I../../LADEL/include -I../../LADEL/amd/Include -I../../LADEL/interfaces/mex/include';
    flags = [flags ' -DUSE_LADEL'];
end

 %---------------------------------------------------------------------------
 % BLAS option
 %---------------------------------------------------------------------------

 % This is exceedingly ugly.  The MATLAB mex command needs to be told where to
 % find the LAPACK and BLAS libraries, which is a real portability nightmare.

if (pc)
    if (verLessThan ('matlab', '7.5'))
        lapack = 'libmwlapack.lib' ;
    else
        lapack = 'libmwlapack.lib libmwblas.lib' ;
    end
else
    if (verLessThan ('matlab', '7.5'))
        lapack = '-lmwlapack' ;
    else
        lapack = '-lmwlapack -lmwblas' ;
    end
end

if (is64 && ~verLessThan ('matlab', '7.8'))
    % versions 7.8 and later on 64-bit platforms use a 64-bit BLAS
%     fprintf ('with 64-bit BLAS\n') ;
    flags = [flags ' -DBLAS64'] ;
end

if (~(pc || mac))
    % for POSIX timing routine
    lapack = [lapack ' -lrt'] ;
end

%% List all the source files
qpalm_src = { ...
    'qpalm',...
    'solver_interface',...
    'lin_alg',...
    'linesearch',...
    'newton',...
    'scaling',...
    'termination',...
    'util',...
    'validate',...
    'nonconvex',...
    'iteration',...
    };
    
qpalm_src_path = '../../src';
for i = 1:length (qpalm_src)
    qpalm_src {i} = [qpalm_src_path '/' qpalm_src{i}] ;
end

if strcmp(solver, 'ladel')
    
    ladel_src = { ...
        'ladel_col_counts', ...
        'ladel_copy', ...
        'ladel_debug_print', ...
        'ladel_etree', ...
        'ladel_global', ...
        'ladel', ...
        'ladel_ldl_numeric', ...
        'ladel_ldl_symbolic', ...
        'ladel_matvec', ...
        'ladel_matmat', ...
        'ladel_submatrix', ...
        'ladel_add', ...
        'ladel_pattern', ...
        'ladel_permutation', ...
        'ladel_postorder', ...
        'ladel_rank1_mod', ...
        'ladel_row_mod', ...
        'ladel_scale', ...
        'ladel_transpose', ...
        'ladel_upper_diag', ...
        };

    ladel_src_path = '../../LADEL/src';
    for i = 1:length(ladel_src)
        ladel_src{i} = [ladel_src_path '/' ladel_src{i}] ;
    end

    ladel_mex_util_src = ['../../LADEL/interfaces/mex/ladel_mex_util']; 

    amd_src = { ...
        'amd_1', ...
        'amd_2', ...
        'amd_aat', ...
        'amd_control', ...
        'amd_defaults', ...
        'amd_dump', ...
        'amd_global', ...
        'amd_info', ...
        'amd_order', ...
        'amd_postorder', ...
        'amd_post_tree', ...
        'amd_preprocess', ...
        'amd_valid', ...
        'SuiteSparse_config', ...
        };

    amd_src_path = '../../LADEL/amd/Source';
    for i = 1:length(amd_src)
        amd_src{i} = [amd_src_path '/' amd_src{i}];
    end
    
    source = [amd_src ladel_src ladel_mex_util_src qpalm_src] ;
    
end

%% Compile everything
if (pc)
    obj_extension = '.obj' ;
else
    obj_extension = '.o' ;
end
obj = '' ;

kk = 0 ;
cflags = 'CFLAGS="\$CFLAGS -std=c99 -fPIC -DMATLAB -DDAMD -O3 -DPROFILING -DPRINTING"';
flags = [cflags ' ' flags];

for f = source
    ff = f {1} ;
    if (strcmp(solver, 'cholmod') && isequal (ff, [metis_path '/GKlib/util']))
        % special case, since a file with the same name also exists in libmetis
        copyfile ([ff '.c'], 'GKlib_util.c', 'f') ;
        ff = 'GKlib_util' ;
        o = 'GKlib_util' ;
    elseif (strcmp(solver, 'cholmod') && isequal (ff, [metis_path '/GKlib/graph']))
        % special case, since a file with the same name also exist in libmetis
        copyfile ([ff '.c'], 'GKlib_graph.c', 'f') ;
        ff = 'GKlib_graph' ;
        o = 'GKlib_graph' ;
    elseif (isequal (ff, [qpalm_src_path '/util']))
        % special case, since a file with the same name also exists in
        % qpalm
        copyfile ([ff '.c'], 'qpalm_util.c', 'f') ;
        ff = 'qpalm_util' ;
        o = 'qpalm_util' ;
    else
        slash = strfind (ff, '/') ;
        if (isempty (slash))
            slash = 1 ;
        else
            slash = slash (end) + 1 ;
        end
        o = ff (slash:end) ;
    end
    % fprintf ('%s\n', o) ;
    o = [o obj_extension] ;
    obj = [obj  ' ' o] ;					            
    s = sprintf ('mex %s -DDLONG -O %s -c %s.c', flags, include, ff) ;
    s = [s obj ' ' lapack] ;
    kk = do_cmd (s, kk, details) ;
end

s = sprintf ('mex %s -DDLONG -O %s %s.c', flags, include, 'qpalm_mex') ;
s = [s obj ' ' lapack] ;	
kk = do_cmd (s, kk, details) ;

 % clean up
s = ['delete ' obj] ;
do_cmd (s, kk, details) ;
fprintf ('\nQPALM successfully compiled\n') ;

% remove the renamed METIS files, if they exist
if (exist ('GKlib_util.c', 'file'))
    delete ('GKlib_util.c') ;
end
if (exist ('GKlib_graph.c', 'file'))
    delete ('GKlib_graph.c') ;
end
if (exist ('qpalm_util.c', 'file'))
    delete ('qpalm_util.c') ;
end

%change to qpalm main path
% cd(current_path);
qpalm_path = fullfile(this_path, '..');
cd(qpalm_path);

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
