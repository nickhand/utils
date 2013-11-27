"""
Generate a LaTeX authors list (to stdout) given a first author
(or multiple lead authors).  Also the associated affiliation 
addresses are generated.

Usage: authors.py leadAuthorLastName [otherLeadLastNames...]
"""
import string, sys


class Author(object):
    """
    Contains a single author.  Knows how to number his/her affiliation.
    Can render his/her name in LaTeX with appropriate affil number.

    WARNING: Could (should?) be made capable of handling a list of
    affiliations for one individual.
    """

    def __init__(self, name, addresses):
        self.name = name

        # Remove non-printing characters from the last name, for sorting
        raw_last_name = name.split()[-1]
        identity = string.maketrans("","") # an identity translation table
        BAD_CHAR = '{}"\''
        self.last_name = raw_last_name.translate( identity, BAD_CHAR )

        self.addresses = addresses

    def __cmp__(self, other):
        """Sorting Author objects should compare by last name first.
        To do: figure out how to get D\"{u}nner to sort as Dunner.
        """
        c1 = cmp(self.last_name, other.last_name)
        if c1 != 0: return c1
        return cmp(self.name, other.name)

    def __affil_mark(self, affiliations):
        ERROR_AFFIL = "unknown"

        marks=[]
        for address in self.addresses:
            address = address.strip()
            try:
                marks.append( "%d"%affiliations.affil_num(address) )
            except KeyError:
                sys.stderr.write( "Unknown affiliation: %s\n"%address)
                return ERROR_AFFIL
        # marks.sort()
        return ",".join(marks)

    def latex(self, affiliations, initials_only = False, last=False, add_email=None):
        """
        affiliations is an Affiliation object.
        Set last = True for the last author in a list (supresses the comma)
        """
        names = self.name.split()
        if names[-1]=="Dunner": names[-1]=r'D\"{u}nner'
        if names[-1]=="Dunner-Planella": names[-1]=r'D\"{u}nner-Planella'
        ret = ""

        # Replace whole names with initials
        dash = False
        if initials_only:
           for i,n in enumerate(names[:-1]):
                if n[0] == "-":
                    names[i] = n[0:2] + "."
                    dash = True
                else:
                    names[i] = n[0] + "."
        for i in range(0, len(names)):
            if i > 0 and names[i][0] != "-":
                ret = ret + "~" + names[i]
            else:
                ret = ret + names[i]
        
        if add_email is None:
            ret = r"\author{%s}" %ret + "\n"
        else:
            ret = r"\author{%s} \email{%s}" %(ret, add_email) + "\n"
        
        
        for ad in self.addresses:
            affil = affiliations.affil_text(ad.strip())
            ret += r"\affiliation{%s}" %affil + "\n"
        return ret


class AuthorList(list):
    """
    Holds a list of authors (by subclassing list).
    Knows how to render the whole author list in LaTeX.
    """

    def __init__(self, *args, **kwargs):
        self.initials_only = False
        list.__init__(self, *args, **kwargs)
        

    def latex(self, affiliations, addEtAl = False, first_author=None, first_email=None):
        lines=[]
        lines.append(75*r"%")
        lines.append(
            "%% WARNING: This LaTeX block was generated automatically by %s"%
            __file__)
        lines.append("% Do not change by hand: your changes will be lost.\n")

        for a in self[:-1]:
            add_email = None
            if a.last_name == first_author and first_email is not None:
                add_email = first_email
            lines.append( a.latex( affiliations, 
                                   initials_only=self.initials_only, add_email=add_email ) )
        
        if addEtAl:
            lines.append( self[-1].latex( affiliations, last=False,
                                      initials_only=self.initials_only) )
            lines.append(" et al.") 
        else:
            lines.append( self[-1].latex( affiliations, last=True,
                                      initials_only=self.initials_only) )
        #lines.append("}")

        #lines.append( affiliations.latex() )
        lines.append("\n% End auto-generated block")
        lines.append(75*r"%")
        return "\n".join(lines)

    def __add__(self, other):
        """Without this method, al1+al2 returns a plain list."""
        out = AuthorList()        
        out.extend(self)
        out.extend(other)
        return out


class Affiliations(dict):
    """
    Contains a group of affiliations, indexed by code name.
    Subclasses dictionary but adds a "numbering" feature.
    The first affiliation used gets a persistent number starting at 1.
    By "used" we mean that its number was checked with affil_num().
    
    Note that numbering is reset explicitly with reset_numbering()
    or implicitly each time a new affiliation is added.
    """

    def affil_num(self, key):
        """Numbers the affiliations in the order they were looked up.
        use Affiliations.reset_numbering() to start over"""

        if key not in self.keys():
            raise KeyError("key %s not known"%key)
        try:
            return self.ordering[key]
        except KeyError:
            self.anum += 1
            self.ordering[key] = self.anum
            return self.anum
    
    def affil_text(self, key):
        """Returns the affilition text"""

        if key not in self.keys():
            raise KeyError("key %s not known"%key)
        try:
            return self[key]
        except KeyError:
            pass

    def reset_numbering(self):
        self.anum = 0
        self.ordering={}

    def __setitem__(self, key, value):
        dict.__setitem__(self, key, value)
        self.reset_numbering()

    def latex(self):
        lines=[]

        ordered_affil={}
        for k in self.keys():
            ordered_affil[self.affil_num(k)] = self[k]

        for codenum,val in ordered_affil.iteritems():
            lines.append(r"\altaffiltext{%d}{%s}" %
                         (codenum, val) )
        return "\n".join(lines)
            



def parse_addresses(string):
    inst = Affiliations()
    lines = string.split("\n\n")
    for l in lines:
        if ":" not in l: continue
        key, address = l.split(":", 1)
        inst[key.strip()] = address.strip()
    return inst


def parse_names(string):
    auth = AuthorList()
    lines = string.split("\n")
    for l in lines:
        if ":" not in l: continue
        name, addressset = l.split(":", 1)
        addresses = addressset.split(",")
        auth.append( Author(name.strip(), addresses) )

    auth.sort()
    return auth
