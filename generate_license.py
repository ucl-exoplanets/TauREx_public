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


from classes.license import *


print 'TauREx License file generator'
print''

FULL = raw_input('Full access (y/n)? [n]: ')
if FULL == '':
    FULL = False
    
if FULL:
    pass
else:
    USER_bool = raw_input('Require Username (y/n)? [y]: ')
    if USER_bool == '' or USER_bool == 'y':
        USER_bool = True
        USER = raw_input('Enter Username: ')
        
    HOST_bool = raw_input('Require Machine Hostname (y/n)? [y]: ')
    if HOST_bool == '' or HOST_bool == 'y':
        HOST_bool = True
        HOST = raw_input('Enter Hostname: ')  
        
    DATE_bool = raw_input('Require Date limit (y/n)? [y]: ')
    if DATE_bool == '' or DATE_bool == 'y':
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
    
        
  