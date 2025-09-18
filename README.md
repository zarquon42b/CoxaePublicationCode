## This README contains instructions to install R-packages to run the publication code
## It runs on Ubuntu >=20.04 (which is also possible to install on Windows using [WSL](https://zarquon42b.github.io/2020/02/22/WSL-RvtkStatismoUpdate/))

Install required R-packages and underlying libraries. We assume that a working R-version is already installed on the machine.
Open a terminal and execute the following code:
	
	sudo apt-add-repository ppa:zarquon42/statismo-develop
	sudo apt update
	sudo apt install statismo-dev


Once all packages are installed, go to R and install the required packages (from CRAN & GitHub)
Start R and rund the following code:
	
```
install.packages("devtools")
install.packages("Morpho")
devtools::install_github("zarquon42b/mesheR")
devtools::install_github("zarquon42b/RvtkStatismo",ref="develop")

```
	
	
