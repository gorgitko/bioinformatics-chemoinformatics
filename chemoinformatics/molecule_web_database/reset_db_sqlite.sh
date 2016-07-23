#!/bin/bash

rm moldb.sqlite3
rm -r moldb/migrations/*

python manage.py makemigrations
python manage.py makemigrations moldb
python manage.py migrate