# Simple web and database of molecules

This is one of my credit's projects.
It's just a simple web connected to database, where you can list, filter, search, order/purchase and so on molecules. I am heavily using
the `jQuery` and sometimes the `AJAX` calls.

# Requirements

I am using the `conda` package manager and it's great: http://conda.pydata.org/docs/

I strongly recommend you to create new environment
using `conda` and [requirements.txt](https://github.com/gorgitko/bioinformatics-chemoinformatics/blob/master/chemoinformatics/molecule_web_database/requirements.txt) or install these SW's manually with `conda install`:

`Python 3.5+`

[Django](https://www.djangoproject.com/) 1.9.2+

`RDKit` stuff:

[RDKit](http://www.rdkit.org/) - `conda` makes the installation easy -- normally it's a nightmare
(compilation, libraries, system variables etc.). Install from http://anaconda.org/rdkit/rdkit

`PostgreSQL` with [`RDKit cartridge`](http://www.rdkit.org/docs/Cartridge.html) - install from http://anaconda.org/rdkit/rdkit-postgresql. PostgreSQL will be installed in `/path/to/your/conda/env/bin`. Than somewhere create working directory for PostgreSQL and start it with `/path/to/your/conda/env/bin/postgres -D /path/to/your/postgresql/working/directory &`. Create new user `moldb` and set its privileges, set password etc. Google will help you :) After that you can use script [`reset_db.sh`](https://github.com/gorgitko/bioinformatics-chemoinformatics/blob/master/chemoinformatics/molecule_web_database/reset_db.sh) for clean DB init or [`init_db.sh`](https://github.com/gorgitko/bioinformatics-chemoinformatics/blob/master/chemoinformatics/molecule_web_database/init_db.sh) to reset DB and insert some initial data with the script [`insert_data.py`](https://github.com/gorgitko/bioinformatics-chemoinformatics/blob/master/chemoinformatics/molecule_web_database/insert_data.py)
