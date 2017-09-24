#/bin/bash

find ./ -type f -exec sed -i -e 's/wellington_20/wellington_16/g' {} \;
find ./ -type f -exec sed -i -e 's/hint_20/hint_16/g' {} \;
