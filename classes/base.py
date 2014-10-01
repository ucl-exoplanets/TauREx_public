################################################
#class base
#
#Class defining common class methods from which other classes
#are derived. 
#
#
# Modification History:
#   - v1.0 : first definition, Ingo Waldmann, Apr 2013 
#
################################################


class base(object):

    #basic class methods and overloading. These are common functions that will be 
    #inherited by other classes 
    
    
    def list(self,name=None):
        if name is None:
            return dir(self)[2:-1]
        else:
            lst = dir(self)
            return filter(lambda k: name in k, lst)
        
    def __getattribute__(self,name):
        return object.__getattribute__(self, name)
    
    def __getitem__(self,name):
        return self.__dict__[name] 

    def reset(self,params):
    #allows to reset the original instance to reflect changes in the data instance
    #this avoids an initialisation of a separate instance.
        self.__init__(params)