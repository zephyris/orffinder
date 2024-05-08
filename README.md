# OrfFinder

A simple Python and tBLASTn-based open reading frame finding tool.

## Installation

OrfFinder requires `python`, and is easiest to install using `pip` and `git`. It also requires `ncbi-blast+` to be installed and available to the system command line.

On Linux, we recommend installing using your normal package manager, for example:
``` shell
sudo apt update
sudo apt install ncbi-blast+ python3 python3-pip git
```

On Windows, download and install `ncbi-blast+` using the installer, available at [https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/].

You can install `python`, `pip` and `git` using `winget` in either Command Prompt or PowerShell:
``` shell
winget install Python.Python.3.0; winget install Git.Git
python3 -m pip install -update pip
```

Install `orffinder` using `pip` and `git`:
``` shell
pip install git+https://github.com/zephyris/orffinder
```

To reinstall and upgrade use `pip` and `git`:
``` shell
pip install --upgrade --force-reinstall git+https://github.com/zephyris/orffinder
```


To uninstall use `pip`
``` shell
pip uninstall orffinder
```

## Usage
