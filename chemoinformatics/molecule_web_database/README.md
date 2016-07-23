# Simple web and database of molecules

This is one of my credit's projects.
It's just a simple web connected to database, where you can list, filter, search, order/purchase and so on molecules. I am heavily using
the `jQuery` and sometimes the `AJAX` calls.

# Requirements

I am using the `conda` package manager and it's great: http://conda.pydata.org/docs/

I strongly recommend you to create new environment
using `conda` and install these SW's (later I will add here the `conda`'s installed packages list output):

`Python 3.5+`

[Django](https://www.djangoproject.com/) 1.9.2+

`RDKit` stuff:

[RDKit](http://www.rdkit.org/) - `conda` makes the installation easy, because normally it's a nightmare
(compilation, libraries, system variables etc.). Install from http://anaconda.org/rdkit/rdkit

`PostgreSQL` with `RDKit cartridge` - install from http://anaconda.org/rdkit/rdkit-postgresql
