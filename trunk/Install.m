function Install()
% On the Matlab prompt, change into the D-MiRaGe directory and run Install.
  dir = pwd;
  addpath(genpath(dir));
  disp(sprintf('\n%s\n','MatrixClust directory added to Matlab search path.'));
  disp(['Click <a href="matlab: doc MatrixClust;">doc MatrixClust</a> for ' ...
        'more information about this toolbox.']);
  disp(sprintf('\n%s','To permanently add MatrixClust to the Matlab search path,'));
  disp('create a "startup.m" file in your home directory containing');
  disp(['the line "addpath(genpath('  dir  '))".']);
  disp(sprintf('\n%s\n',['Click <a href="matlab: doc startup;">doc startup</a> for more ' ...
        'information about startup options.'])); 

