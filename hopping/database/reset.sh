#!/bin/bash
set -e  # halts execution on error

database=thin_films.db
init_sql=init.sql

if [ -f "$database" ]; then
    rm -f $database
    echo "Deleted database $database! >:)"
    sqlite3 $database < $init_sql
    echo "Succesfully created new database $database from $init_sql! :D"
else
    sqlite3 $database < init.sql
    echo "Database $database didn't exist! Don't worry, I created one for you from $init_sql :)"
fi
