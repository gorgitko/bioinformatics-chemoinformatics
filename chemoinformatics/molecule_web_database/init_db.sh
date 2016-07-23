#!/bin/bash

dropdb moldb;

psql -d template1 -c "CREATE DATABASE "moldb"
  WITH OWNER "moldb"
  ENCODING 'UTF8'
  LC_COLLATE = 'en_US.UTF-8'
  LC_CTYPE = 'en_US.UTF-8';"
psql -d moldb -c "CREATE EXTENSION rdkit;"

rm -r ./moldb/migrations/__pycache__
rm ./moldb/migrations/*.py
python manage.py makemigrations moldb
python manage.py migrate

python insert_data.py