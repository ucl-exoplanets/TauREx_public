################################################
#class base
#
#Class defining common class methods from which other classes
#are derived. 
#
#
# Modification History:
#   - v1.0 : first definition, Ingo Waldmann, Sept 2014 
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

    # def __getattribute__(self,name):
    #     return object.__getattribute__(self, name)
    #
    # def __getitem__(self,name):
    #     return self.__dict__[name]

    def reset(self,params, **kwargs):
    #allows to reset the original instance to reflect changes in the data instance
    #this avoids an initialisation of a separate instance.
        self.__init__(params,**kwargs)
        
    
    def set_model(self, INPUT = None):
        #loads emission/transmission model pointer into fitting class
        if INPUT == None: 
            self.model = None
            self.__MODEL_ID__ = None
        else:
            if INPUT.__ID__ == 'transmission':
                model = INPUT.cpath_integral
            elif INPUT.__ID__ == 'emission':
                model = INPUT.path_integral
                
            self.model    = model
            self.__MODEL_ID__ = INPUT.__ID__
            