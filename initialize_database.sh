#!/bin/bash

source env.sh

python -c '
from immuno.ui import db
db.create_all()
'
