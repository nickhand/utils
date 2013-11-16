import re

class bibitem(dict):
    """
    A class designed to emulate a bibitem in a bibtex file
    """
    #---------------------------------------------------------------------------
    def __init__(self, bib_type, bib_id, bib_info):
        """
        Initialize the bibitem as a dictionary according to the input dictionary
        
        Parameters
        ----------
        bib_type : str
            the type of entry for the bibitem, i.e., article, etc
        bib_id : str
            the identification key used for this bibitem
        bib_info : dict
            the dictionary of key and value pairs that holds the bibitem info
        """ 
        # compile some regexes
        self.white = re.compile(r"[\n|\s]+")
        self.nl = re.compile(r"[\n]")
        self.token_re = re.compile(r"([^\s\"#%'(){}@,=]+|\n|@|\"|{|}|=|,)")

        self.id = bib_id
        self.type = bib_type
        
        # initialize hidden variables
        self._author_count = 0
        self._authors = ()
        self._parsed_names = ()
        
        # set the dictionary values
        for k, v in bib_info.items():
            dict.__setitem__(self, k, v)
         
        # find the author key   
        self.author_key = None
        if 'Author' in self.keys():
            self.author_key = 'Author'
        elif 'author' in self.keys():
            self.author_key = 'author'
            
        # parse the author list
        if self.author_key is not None: self._parse_authors()
    #end __init__
    
    #---------------------------------------------------------------------------
    def add(self, key, val):
        """
        Add a field to the bibitem
        """
        self[key] = "{%s}" %str(val)
    #end add
    
    #---------------------------------------------------------------------------
    def delete(self, key, val):
        """
        Delete a field from the bibitem
        """
        if key in self.keys():
            del self[key]
    #end delete
    
    #---------------------------------------------------------------------------
    def keys(self):
        """
        Return the bibitem field items
        """
        return sorted([k for k in self], key=str.lower)
    #end keys
    fields = keys
    #---------------------------------------------------------------------------
    def __setitem__(self, key, val):
        """
        Bultin python setitem function. Re-parses the authors if we are updating
        the author key
        """
        dict.__setitem__(self, key, val)
        if key == self.author_key:
            self._parse_authors()
    #end __setitem__
        
    #---------------------------------------------------------------------------
    def _get_author_count(self):
        return len(self._authors)
 
    author_count = property(_get_author_count)
    
    #---------------------------------------------------------------------------
    def _get_authors(self):
        return self._authors   
    
    def _set_authors(self, val):
        if self.author_key is None: raise KeyError("no author field for this bibitem")
        if not isinstance(val, (tuple, list)):
            val = (val, )
        elif isinstance(val, list):
            val = tuple(val)
        self[self.author_key] = val
        
    authors = property(_get_authors, _set_authors)
    
    #---------------------------------------------------------------------------
    def _get_parsed_names(self):
        return self._parsed_names
        
    parsed_names = property(_get_parsed_names)
    
    #---------------------------------------------------------------------------
    def _authors_to_bibtex(self):
        authors = []
        for a in self._authors:
            if 'given' in a.keys():
                authors.append("%s, %s" %(a['family'], a['given']))
            else:
                authors.append("%s" %a['family'])
        ret = ' and '.join(authors)
        return "{%s}" %ret
    #---------------------------------------------------------------------------
    def _parse_authors(self):
        """
        Internal function to parse the authors bibtex input
        """
        if self.author_key is None: raise KeyError("no author field for this bibitem")
        self._authors = self[self.author_key]
        self._parsed_names = ()
        
        for author in self._authors:
            family = author['family']
            last_name = family[family.find("{")+1:family.rfind("}")]
            if 'given' in author.keys():
                first_name = author['given']
                self._parsed_names += (first_name + " " + last_name, )
            else:
                self._parsed_names += (last_name, )
    #end _parse_authors
    
    #---------------------------------------------------------------------------
    def __str__(self):
        """
        The string representation of the bibitem
        """     
        copy = self.copy()
        if self.author_key is not None:
            copy[self.author_key] = self._authors_to_bibtex()
            
        val = "@%s{%s, \n" %(self.type, self.id)
        val += ",\n".join("%s = %s" %(k, copy[k]) for k in sorted(self.keys(), key=str.lower))
        val += "}\n"
        return val
    #end __str__
    
    #---------------------------------------------------------------------------
    def __repr__(self):
        """
        The builtin representation function
        """
        return "bibitem of type %s with key %s" %(self.type, self.id)
    #end __repr__
    
    #---------------------------------------------------------------------------
    def truncate_authors(self, N=1, ending='others'):
        """
        Truncate the author list to N authors, with the option of appending
        an ending to the author list
        """
        if self.author_count > 0:
            if self.author_count > N: 
                self.authors = self.authors[:N]
            
            if ending is not None:
                self[self.author_key] += ({'family' :ending}, )
    #end truncate_authors
    
    #---------------------------------------------------------------------------
    def info(self):
        """
        Print out the info for this bibitem
        """
        title = self.get('Title', self.get('title', None))
        print "key: %s" %self.id
        if title is not None: print "title: %s" %title
        print "authors: %s" %str(self.parsed_names)
    #end info
#endclass bibitem
#-------------------------------------------------------------------------------
class bibfile(object):
    """
    A class designed to emualte a bibtex file
    """
    def _tokenize(self, ignore_white=False) :
        """
        Return a token iterator
        """        
        for item in self.token_re.finditer(self.data):
            i = item.group(0)
            if ignore_white:
                if self.white.match(i) :
                    continue
            if self.nl.match(i) :
                self.line += 1
            yield i            
    #end _tokenize
    
    #---------------------------------------------------------------------------
    def __init__(self, filename):
        """
        Initialize the bibfile, given the bibtex file name
        
        Parameters
        ----------
        filename : str
            the name of the bibtex file to initialize
        """
        
        self.data = open(filename, 'r').readlines()    
        self.items = {}        
        self.line = 0

    #end __init__
    
    #---------------------------------------------------------------------------
    def __iter__(self):
        return self.items.__iter__()
        
    def iteritems(self):
        return self.items.iteritems()
        
    def items(self):
        return self.items.items()
        
    def __getitem__(self, key):
        return self.items.__getitem__(key)
        
    def __setitem__(self, key, val):
        return self.items.__setitem__(key, val)
    #---------------------------------------------------------------------------
    def keys(self):
        """
        Return the alphabetically sorted keys of the entries in the file
        """
        return sorted(self.items.keys(), key=str.lower)
    #end keys
    
    #---------------------------------------------------------------------------
    @property
    def size(self):
        """
        Define a size attribute that returns the number of bibitems in the file
        """
        return len(self.items)
    #end size
    
    #---------------------------------------------------------------------------
    def __repr__(self):
        return "bibfile of size %d" %self.size
    #end __repr__

    #---------------------------------------------------------------------------
    def __len__(self):
        """
        The builtin length function returns the size of the bib file
        """
        return self.size
    #end __len__
    
    #---------------------------------------------------------------------------
    def __str__(self):
        """
        Return the joined strings of each bibitem in the bib file
        """
        return "\n".join(str(self.items[k]) for k in self.keys())
    #end __str__
    
    #---------------------------------------------------------------------------
    def write(self, filename):
        """
        Write out the bib file to the specified filename
        """
        out = open(filename, 'w')
        out.write("%s" %self)
        out.close()
    #end write
    
    #---------------------------------------------------------------------------
    def load(self):
        """
        Load the data from the bib file
        """
        in_field       = False  # keep track of whether we found a entry field
        unclosed_brackets = 0      # keep track of number of unclosed brackets     
        
        # loop over each line in the file
        for line in self.data:
            self.line += 1
            if line.isspace(): continue  # skip white space
                
            # now in an entry
            if "@" in line and not unclosed_brackets:
                
                # find the entry type and id
                entry_type = re.findall("@(\w+){", line)[0].strip()
                entry_id = re.findall("{(.*?),", line)[0].strip()
                
                in_entry = True
                info = {}
                
                # count the unclosed brackets
                unclosed_brackets += self._check_brackets(line)
                
            elif "@" in line and unclosed_brackets:
                raise ValueError("unclosed brackets on line %d" %(self.line-1))
            elif unclosed_brackets:
                
                # allow for one word and then = sign, and remove empty strings and space
                if unclosed_brackets == 1:
                    line = line.split('=', 1) #filter(None, re.split(r'(\s*\b[a-zA-Z0-9\-]+\b)\s*=', line))
                else:
                    line = [line]
                line = [x.strip() for x in line] 
                
                # skip empty lines
                if len(line) == 0 or line[0] == '': 
                    continue
                    
                # update the bracket count
                unclosed_brackets += self._check_brackets(line)    
                
                # if line just contains spaces and onethen skip to next line
                if len(line) == 1 and not unclosed_brackets:
                    self.items[entry_id] = bibitem(entry_type, entry_id, info)
                    in_field = False
                    continue
                
                # continue if its just a single comma
                if len(line) == 1 and re.match("\s*[,]+\s*", line[0]):
                    continue
                
                # if line has two components, we just read a new field value
                if len(line) == 2:
                    
                    current_line = line[1]                    
                    field = line[0]
                    in_field = True
                # only one component, then we are reading a value separated on multiple lines
                else:
                    if not in_field:
                        raise ValueError("error finding the field value of line %d" %self.line)
                    
                    # join the current line with spaces
                    current_line = " ".join([current_line, line[0]])
                
                # skip ahead if we are in the middle of reading a multiline field value
                if unclosed_brackets > 1:
                    continue
                
                # check for comma at end of line
                if current_line.rstrip()[-1] != "," and unclosed_brackets:
                    value = current_line[:]
                else:
                    value = current_line[:-1]
                    
                # parse author list
                if field.lower() == 'author':
                    value = self._parse_authors(value)
                
                # store the field and value in the info dict
                info[field] = value
                in_field = False
                
                # make the bibitem
                if not unclosed_brackets:
                    self.items[entry_id] = bibitem(entry_type, entry_id, info)
            else:
                raise ValueError("parsing error on line %d" %self.line) 
            
        if unclosed_brackets != 0:
            raise ValueError("unclosed brackets found when parsing file; "
                            "last field read was the '%s' field" %field)
    #end load
    
    #---------------------------------------------------------------------------
    def _check_brackets(self, line):
        
        count = 0
        if not isinstance(line, list):
            line = [line]
        for x in line:
            for char in x:
                if char == "}":
                    count -= 1
                elif char == "{":
                    count += 1
        return count
    #---------------------------------------------------------------------------
    def remove(self, key):
        """
        Remove the entry from the bib file that corresponds to the input key
        """
        if self.contains(key):
            del self.items[key]
    #end remove
    
    #---------------------------------------------------------------------------
    def contains(self, key):
        """
        Return whether the bibfile contains the input key
        """
        return key in self.keys()
    #end contains
    
    #---------------------------------------------------------------------------
    def find(self, key):
        """
        Return the bibitem corresponding to key
        """
        if self.contains(key):
            return self.items[key]
        else:
            return None
    #end find
    
    #---------------------------------------------------------------------------
    def _parse_authors(self, authors):
        authors = authors[authors.find("{")+1 : authors.rfind("}")]
        res = []        
        authors = authors.split(' and ')
        for author in authors:
            _author = author.split(',')
            family = _author[0].strip().rstrip()
            rec = {'family':family}
            try :
                given = _author[1].strip().rstrip()
                rec['given'] = given
            except IndexError:
                pass
            res.append(rec)
        return res
    #end _parse_authors
#endclass bibfile
