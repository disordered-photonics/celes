In order to execute the CELES code, it is necessary to have the following installed:



1. The CUDA Toolkit. 

You can check the verison that you need by running the command "gpuDevice" in MATLAB and look for ToolkitVersion in the output.



2. A C++ compiler which is suppoerted by MATLAB in combination with the given CUDA Toolkit version.

Linux users:
The built-in gcc compiler works fine.

Windows users:
One combination that works is MATLAB 2016b, CUDA Toolkit 7.5 and MS Visual Studio 2013 (not 2015).
With MATLAB 2017a: CUDA Toolkit 8 and MS Visual Studio 2013 .


-> MS Visual Studio 2013 can be downloaded from:
https://www.microsoft.com/en-us/download/details.aspx?id=48138