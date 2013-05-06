"""
 param_dict.py
 class to handle parameter dictionaries
 
 author: Nick Hand 
 contributions: Toby Marriage, Sudeep Das, Jake VanderPlas
 contact: nhand@berkeley.edu
 creation date: 05/01/2013
"""

import os
import string
import re

def ask_for(key):
    """
    @brief ask for the input key/value
    """
    s = raw_input("param_dict: enter value for '%s': " %key)
    try:
        val = eval(s)
    except NameError:
        # allow people to enter unquoted strings
        val = s
    return val

class param_dict( object ):
    """
    @brief param_dict is a borg class.  Any instance points to the
    common dictionary __shared_state.  Any modification to a single
    param_container instance modifies the state of all 
    param_container instances.
    """
    
    # create a shared parameter dictionary
    __shared_state = {}
    
    def __init__(self,ask=False):
        self.__dict__ = self.__shared_state
        self.ask = ask

    def __getitem__(self,item):
        if item not in self.keys():
            if self.ask:
                print "param_dict: parameter '%s' not found" %item
                val = ask_for(key)
                print "param_dict: setting '%s' = %s" % (item,repr(val))
                setattr(self,item,val)
            else:
                return None
        return getattr(self,item)

    def __setitem__(self,item,val):
        return setattr(self,item,val)

    def update_param_dict(self,D):
        self.__dict__.update(D)

    def clear_param_dict(self):
        keys = self.__dict__.keys()
        x = self.ask
        for k in keys: del self.__dict__[k]
        self.ask = x

    def get_keys(self):
        return self.__dict__.keys()
    keys = get_keys
    
    def dump(self):
        """
        @brief print all parameters and types
        """
        print "parameters"

        lines = [ (p,str(v),type(v).__name__) \
                      for p,v in self.__dict__.iteritems() ]

        if len(lines)==0:
            print ' [no parameters]'

        else:
            L0 = max([len(l[0]) for l in lines])
            L1 = max([len(l[1]) for l in lines])

            for p,v,t in lines:
                print ' %s %s %s' % (p.ljust(L0),v.ljust(L1),t)
                
    def load(self, filename, clear_current=False):
        """
        @brief Fill the package variable params with the parameters specified
        in the filename.

        if clear_current is True, then first empty current parameter settings

        filename has one parameter per line, space delimited:
        name val [type] #comment

        if type is not specified, string type is assumed

        For example:
          #start params.dict
          dir = "/local/tmp/"              
          pi  = 3.14
          num_entries = 3                       
          file1 = "$(dir)/myfile1.txt"          # dollar sign indicates variable
                                                # substitution
          file2 =  dir + "myfile2.txt"
          #end params.dat
        """
        if clear_current:
            D = {}
        else:
            D = self.__dict__
        linecount = 0
        for line in open(filename):
            linecount += 1
            line = ' '.join(line.split('#')[:1]).split('\\')
            line = ' '.join(line)
            exec(line)
            line = line.split('=')
            line = [x.strip() for x in line]

            if len(line)==0 or line[0]=='': 
                continue
            if len(line)==1:
                raise ValueError, \
                    "Must specify value for parameter %s on line %i" \
                    % (line[0], linecount)
            elif len(line) != 2:
                raise ValueError, "cannot understand line %i of %s" \
                                    % (linecount,filename)     
            if (line[0][0]>='0' and line[0][0]<='9'):
                raise ValueError, "invalid variable name %s" % line[0]

            # check for variables in the value
            if '$' in line[1]:
                line[1] = replace_vars(line[1],D)
             
            print line[1]  
            D[line[0]] = eval(line[1])

        if clear_current:
            self.clear_param_dict()
            self.update_param_dict(D)

    def dump_to_file(self,filename,mode = 'w'):
        f = open(filename, mode)
        keys = self.keys()
        keys.sort()
        for key in keys:
            f.write("%s = %s\n" % (key,repr(self[key])))
        f.close()
        
    def unify(self, D):
        self.clear_param_dict()
        self.update_param_dict(D.__dict__)
    update = update_param_dict    
    
        
def replace_vars(s,D):
    """
    @brief given a string s and a dictionary of variables D, replace all variable
    names with the value from the dict D

    variable names in s are denoted by '$(' at the beginning and ')' at
    the end, or '$' at the beginning and a non-variable character at the
    end.  Variables must be valid python variable names, that is they
    consist of only alphanumeric characters (A-Z,a-z,0-9) and underscores,
    and cannot start with a number.

    example:
    >> D = {'my_var1' : 'abc',
            'my_var2' : '123' }
    >> s = "I know my $(my_var1)s and $my_var2's"
    >> print replace_vars(s,D)

    I know my abcs and 123's
    """
    s_in = str(s) 
    s_out = ''

    while True:
        i = s_in.find('$')
        if i==-1:
            s_out += s_in
            break

        s_out += s_in[:i]
        s_in = s_in[i+1:]

        if len(s_in)==0:
            raise ValueError, "trailing $"

        elif s_in[0] == '(':
            i = s_in.find(')')
            if i==-1:
                raise ValueError, "unmatched '('"
            var = s_in[1:i]

            s_in = s_in[i+1:]
            try:
                s_out += str(D[var])
            except:
                s_out += os.environ[var]

        else:
            var = ''
            i = 0
            while True:
                if i>=len(s_in):
                    break
                s = s_in[i]
                if (s >= 'a' and s <= 'z') \
                        or (s >= 'A' and s <= 'Z') \
                        or (s >= '0' and s <= '9') \
                        or s=='_':
                    var += s
                    i += 1
                else:
                    break
            s_in = s_in[i:]
            s_out += str(D[var])
    return s_out
