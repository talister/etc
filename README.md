# etc
This is a generalized Exposure Time Calculator for imaging and spectroscopic instruments for the optical-NIR regime (thermal emission is not included). It also contains a growing collection of "optical components" data in `etc/data/comp` for:
* atmospheric transmission
* mirror reflectivities (protected and enhanced Aluminium, gold and silver coatings)
* CCD, CMOS, NIR detector QE curves (23 so far)
* various optical glasses (BK7, CaF2, FK5, fused silica etc)
* sky background/atmospheric emission data
* filter transmissions for LCO and selected ones from ESO, ING, Gemini, SOAR

## Usage

## Development Guide - Getting started

Before installing any dependencies or writing code, it's a great idea to create a
virtual environment. We primarily use regular `virtualenv`s to manage virtual
environments, although `conda` *may* work. Using virtualenv, you can run the
following to create and activate a new environment:

```
>> python -m venv /path/to/new/virtual/environment
>> source /path/to/new/virtual/environment/bin/activate
```

(On Windows, invoke the venv command as follows:

```
>> C:\>Python310\python -m venv C:\path\to\myenv
>> C:\> C:\path\to\myenv\Scripts\activate.bat (cmd.exe)
or
>> PS C:\> C:\path\to\myenv\Scripts\Activate.ps1 (PowerShell)
```

If you have conda installed locally, you can run the following to
create and activate a new environment.

```
>> conda create env -n <env_name> python=3.10
>> conda activate <env_name>
```

Once you have created a new environment, you can install this project for local
development using the following commands:

```
>> pip install -e .
>> pip install -e .'[dev]'
>> pre-commit install
```
