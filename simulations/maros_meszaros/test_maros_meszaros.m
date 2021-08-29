myFolder = 'maros/meszaros/maros_meszaros_data_old'; %test set to run
myFile = 'output/maros_meszaros'; %file to save results in
myFile = 'output/MM_1e6'; %file to save results in

run_from_scratch = false;

if run_from_scratch
    options.qpalm_matlab = false;
    options.qpalm_c = true;
    options.osqp = true;
    options.qpoases = false;
    options.gurobi = true;

    options.SCALING_ITER = 10;
    options.VERBOSE = false;
    options.TIME_LIMIT = 3600; 
    options.MAXITER = 1e8;

    options.EPS_ABS = 1e-3;

    [results, options] = run_test_set(myFolder, options);

    save(myFile, 'results', 'options');

else
    load(myFile);
end

postprocess_maros_meszaros(results, options);


