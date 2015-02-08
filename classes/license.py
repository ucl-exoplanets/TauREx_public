#! /usr/bin/python

################################################
#class license
#
# License manager for taurex code. It can generate new license keys and 
# validate existing. Two forms will be implemented. Personalised time dependent 
# and time independent licenses. The keys are encrypted using base64 and stored 
# in an ascii file which the license manager needs to read in.  
#
# Input: -license.dat
#
#
# Output: -  none
#
#
# Modification History:
#   - v1.0 : first definition, Ingo Waldmann, Oct 2014
#
################################################


import numpy as np
import datetime,getpass,json,base64, socket, random,logging




class license_manager(object):
    def __init__(self):
        
        logging.info('Loading License manager')
        
        self.local_user   = getpass.getuser()
        self.local_host   = socket.getfqdn()
        self.current_time = datetime.datetime.now()
        self.seed         = 1567 #encryption seed
        
        
#         self.generate_license_file(FULLACCESS=False,USER='ingowaldmann',DATE='2014-10-10')
#         self.generate_license_file(FULLACCESS=True)


        
    def run(self):
        #run license check
        self.read_license_file('./')
        self.verify()
    
    def print_error(self):
        #print error message
        logging.error('Invalid or expired license file.')
        logging.info('To obtain a new license please contact Ingo Waldmann (ingo@star.ucl.ac.uk)')

        
        
    
    def verify(self):
        #checking if license is valid or not
                
        valid    = False #final tag
        validmin = False #minimal usage
        validuser= False
        validtime= False
        validhost= False
        
        if self.license['FULLACCESS']:
            valid = True; validmin = True; validuser = True; validtime = True; validhost = True;
        else:
            validmin = self.license['MINIMALACCESS']
            
            #checking for date requirement 
            if self.license['REQUIRE DATE']:
                license_date = datetime.datetime.strptime(self.license['DATE'], "%Y-%m-%d")
                if license_date > self.current_time:
                    validtime = True 
                    diff = (license_date - self.current_time).days
                    if int(diff) < 10:
                        logging.warn('Warning: License expires in '+str(diff)+' day(s)')
                else:
                    validtime = False        
            else:
                validtime = True
                
            #checking for user requirement
            if self.license['REQUIRE USER']:
                if self.license['USER'] == self.local_user:
                    validuser = True
                else:
                    validuser = False
            else:
                validuser = True
            
            #checking for hostname requirement
            if self.license['REQUIRE HOSTNAME']:
                if self.license['HOSTNAME'] == self.local_host:
                    validhost = True
                else:
                    validhost = False
            else:
                validhost = True


        #verify flags
        if validuser and validtime and validhost:
            valid = True 
        
        #end program if valid = False
        if valid:
            logging.info('Valid TauREx license found')
        else:
            self.print_error()
            exit()
        
    
    def read_license_file(self,PATH):
        #reading in license file
        
        #read file
        try:
            with open(PATH+'license.dat','r') as infile:
                content = infile.readlines()
        except IOError:
            logging.error('No license.dat file found in base directory.')
            exit()
            
        #decoding
        licensekey = content[2][:-1]
        try:
            decoded = base64.b64decode(licensekey)
        except TypeError:
            self.print_error()
            exit()
        
        #derandomising        
        random.seed(self.seed)
        randomidx = random.sample(xrange(len(decoded)),len(decoded))
        derandomidx = np.argsort(randomidx)
        de_rand = [decoded[i] for i in derandomidx]
        de_randstr = ''.join(de_rand)
           
        #converting to dictionary
        self.license = eval(de_randstr,{'false': False, 'true': True, 'null': None})

        
        
    def generate_license_file(self,FULLACCESS= False, MINIMALACCESS=False, DATE= None, USER = None, HOSTNAME= None):
        #generating license file 
        #
        #FULLACCESS:     no restriction to date or user or usage
        #MINIMALACCESS:  only use forward models but no retrieval 
        #TIMEVAR:        license times out at given DATE
        #DATE:           date on which licence expires. Format: YEAR:MONTH:DAY eg. 2014:10:08
        #USER:           if set only specified system user (getpass.getuser()) can use license
        #############################
        
        #if FULLACCESS = False, either User or Date requirement needed
        if FULLACCESS is False:
            if (USER is False) and (DATE is False): 
                print 'Error: Either User or Date requirement needed if FULLACCESS = False' 
                exit()
                        
        #generating dictionary
        taurex_license = {}
        taurex_license['TAUREX_LICENSE_VERSION'] = '1.0'
        taurex_license['TIMESTAMP'] = str(datetime.datetime.now())
        if FULLACCESS:
            taurex_license['FULLACCESS'] = True
        else:
            taurex_license['FULLACCESS'] = False
            taurex_license['MINIMALACCESS'] = MINIMALACCESS
            if USER is None:
                taurex_license['REQUIRE USER'] = False
            else:
                taurex_license['REQUIRE USER'] = True
                taurex_license['USER'] = USER
            if DATE is None:
                taurex_license['REQUIRE DATE'] = False
            else:
                taurex_license['REQUIRE DATE'] = True
                taurex_license['DATE'] = DATE
            if HOSTNAME is None:
                taurex_license['REQUIRE HOSTNAME'] = False
            else:
                taurex_license['REQUIRE HOSTNAME'] = True
                taurex_license['HOSTNAME'] = HOSTNAME
                
        #parsing dictionary into string 
        dicdump = json.dumps(taurex_license)
        random.seed(self.seed)
        randomidx = random.sample(xrange(len(dicdump)),len(dicdump))
    
        dicrand = [dicdump[i] for i in randomidx]
        dicrandstr = ''.join(dicrand)
           
        encoded = base64.b64encode(dicrandstr)

        #setting up file
        header = 'Begin of TauREx License\n####################################\n'
        footer = '####################################\nEnd of License'
        
        with open('license.dat','w') as outfile:
            outfile.write(header)
            outfile.write(encoded+'\n')
            outfile.write(footer)
            
            
       
        
            
            
    

