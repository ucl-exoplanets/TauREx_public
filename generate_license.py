#! /usr/bin/python

################################################
#script generate_license
#
# script generating license file
#
# Input: user input
#
#
# Output: - license.dat
#
#
# Modification History:
#   - v1.0 : first definition, Ingo Waldmann, Oct 2014
#
################################################

import sys
sys.path.append('./classes')
import license
from license import *


print 'TauREx License file generator'
print''


FULL = False; USER_bool = False; HOST_bool = False; DATE_bool = False;

FULL = raw_input('Full access (y/n)? [n]: ')
if FULL == '' or FULL == 'n':
    FULL = False
else:
    FULL = True
    
if FULL:
    pass
else:
    USER_input = raw_input('Require Username (y/n)? [y]: ')
    if USER_input == '' or USER_input == 'y':
        USER_bool = True
        USER = raw_input('Enter Username: ')
        
    HOST_input = raw_input('Require Machine Hostname (y/n)? [y]: ')
    if HOST_input == '' or HOST_input == 'y':
        HOST_bool = True
        HOST = raw_input('Enter Hostname: ')  
        
    DATE_input = raw_input('Require Date limit (y/n)? [y]: ')
    if DATE_input == '' or DATE_input == 'y':
        DATE_bool = True
        DATE = raw_input('Enter Date [YEAR-MONTH-DAY]: ')
    
licenseob = license_manager()    

if FULL:
    licenseob.generate_license_file(FULLACCESS=True)
elif (USER_bool is True) and (DATE_bool is False) and (HOST_bool is False):
    licenseob.generate_license_file(FULLACCESS=False,USER=USER)
    
elif (USER_bool is False) and (DATE_bool is True) and (HOST_bool is False):
    licenseob.generate_license_file(FULLACCESS=False,DATE=DATE)
    
elif (USER_bool is False) and (DATE_bool is False) and (HOST_bool is True):
    licenseob.generate_license_file(FULLACCESS=False,HOSTNAME=HOST)
    
elif (USER_bool is True) and (DATE_bool is True) and (HOST_bool is False):
    licenseob.generate_license_file(FULLACCESS=False,USER=USER,DATE=DATE)   

elif (USER_bool is False) and (DATE_bool is True) and (HOST_bool is True):
    licenseob.generate_license_file(FULLACCESS=False,DATE=DATE,HOSTNAME=HOST)

elif (USER_bool is True) and (DATE_bool is False) and (HOST_bool is True):
    licenseob.generate_license_file(FULLACCESS=False,USER=USER,HOSTNAME=HOST)   
    
elif (USER_bool is True) and (DATE_bool is True) and (HOST_bool is True):
    licenseob.generate_license_file(FULLACCESS=False,USER=USER,DATE=DATE,HOSTNAME=HOST)

print 'done'
    
        
  