#!/bin/sh

#  clean_hist.sh
#  
#
#  Created by Ingo on 28/09/2017.
#

git checkout --orphan new_branch
git add -A
git commit -am "clean commit"
git branch -D master
git branch -m master
git push -f origin master
