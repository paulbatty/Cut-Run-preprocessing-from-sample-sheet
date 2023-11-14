Creating environments to perform analysis is a good idea as in this way by specifying particular versions of a package then they should be 
reproducible by other people who install the same set of packages as you.

The environment necessary to carry out analyses with this script can be installed locally using conda with the 'base_env.yml' file.

Install Conda. Then generate the enviroment with the following command (replace 'environment' with the name of your enviroment':

conda env create -f environment.yml

To activate your enviroment run the following command (again replacing 'my_env' with your environment name).

conda activate my_env

If you want to use Visual Studio Code to carry out your analyses you can select your enviroment by activating the Python interpreter (Ctrl+Shift+P)
and then choosing your environment.


