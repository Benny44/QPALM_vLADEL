close all;
clear;

myFolder = 'cutest_set/nonconvexQPs'; %test set to run
myFile = 'output/cutest'; %file to save results in

run_from_scratch = false;

if run_from_scratch
    options.qpalm_c = false;
    options.ipopt = false;

    options.SCALING_ITER=10;
    options.EPS_ABS=1e-6;
    options.VERBOSE = true;
    options.TIME_LIMIT = 3600;
    options.MAXITER = 1e9;
    options.NONCONVEX = true;

    [results, options] = run_test_set(myFolder, options);

    save(myFile, 'results', 'options');

else
    load(myFile);
end

postprocess_cutest(results, options);

