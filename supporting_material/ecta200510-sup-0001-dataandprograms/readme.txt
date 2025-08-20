All numerical results in the paper were computed with fortran program "main.f90" in this folder, which in turn calls on subroutines contained in modules in the other files of this folder. 

Compilation uses the Intel® Fortran Compiler 19.1.3.311, along with the IMSL Numerical Libraries Version 7.01 with options

/nologo /O2 /QxHost /Qip /fpp /I"C:\Program Files (x86)\VNI\imsl\fnl701\Intel64\include\dll" /I"C:\Program Files (x86)\VNI\imsl\fnl701\Intel64\include\static" /assume:nosource_include /assume:nocc_omp /Qopenmp /standard-semantics /real-size:64 /names:uppercase /module:"x64\Release\\" /object:"x64\Release\\" /Fd"x64\Release\vc150.pdb" /check:none /libs:static /threads /Qmkl:sequential /c

The US_Light data are in the folder "data" and are read in the module "locmod.f90". The path to these files needs to be adjusted as necessary in that module.
