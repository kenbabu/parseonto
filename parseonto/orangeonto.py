from itertools import chain
from collections import defaultdict
from StringIO import StringIO
import itertools
import re

CONVERSION_DICTIONARIS=[

]


BUILTIN_OBO_OBJECTS = [
"""[Typedef]
id: is_a
name: is_a
range: OBO:TERM_OR_TYPE
domain: OBO:TERM_OR_TYPE
definition: The basic subclassing relationship [OBO:defs]"""
,
"""[Typedef]
id: disjoint_from
name: disjoint_from
range: OBO:TERM
domain: OBO:TERM
definition: Indicates that two classes are disjoint [OBO:defs]"""
,
"""[Typedef]
id: instance_of
name: instance_of
range: OBO:TERM
domain: OBO:INSTANCE
definition: Indicates the type of an instance [OBO:defs]"""
,
"""[Typedef]
id: inverse_of
name: inverse_of
range: OBO:TYPE
domain: OBO:TYPE
definition: Indicates that one relationship type is the inverse of another [OBO:defs]"""
,
"""[Typedef]
id: union_of
name: union_of
range: OBO:TERM
domain: OBO:TERM
definition: Indicates that a term is the union of several others [OBO:defs]"""
,
"""[Typedef]
id: intersection_of
name: intersection_of
range: OBO:TERM
domain: OBO:TERM
definition: Indicates that a term is the intersection of several others [OBO:defs]"""
]
    
def _split_and_strip(string, sep):
    head, tail = string.split(sep, 1)
    return head.rstrip(" "), tail.lstrip(" ")


class OBOObject(object):
    """ Represents a generic OBO object (e.g. Term, Typedef, Instance, ...)
    Example::
        >>> term = OBOObject(stanza_type="Term", id="FOO:001", name="bar")
    """
    def __init__(self, stanza_type="Term", **kwargs):
        """ Init from keyword arguments.
        Example::
            >>> term = OBOObject(stanza_type="Term", id="FOO:001", name="bar", def_="Example definition { modifier=frob } ! Comment")
            >>> term = OBOObject(stanza_type="Term", id="FOO:001", name="bar", def_=("Example definition", [("modifier", "frob")], "Comment"))
            >>> term = OBOObject(stanza_type="Term", id="FOO:001", name="bar", def_=("Example definition", [("modifier", "frob")])) # without the comment
            >>> term = OBOObject(stanza_type="Term", id="FOO:001", name="bar", def_=("Example definition",)) # without the modifiers and comment
        """
        self.stanza_type = stanza_type
        
        self.modifiers = []
        self.comments = []
        self.tag_values = []
        self.values = {}
        
        sorted_tags = sorted(kwargs.iteritems(), key=lambda key_val: chr(1) if key_val[0] == "id" else key_val[0])
        for tag, value in sorted_tags:
            if isinstance(value, basestring):
                tag, value, modifiers, comment = self.parse_tag_value(self.name_demangle(tag), value)
            elif isinstance(value, tuple):
                tag, value, modifiers, comment = ((self.name_demangle(tag),) + value + (None, None))[:4]
            self.add_tag(tag, value, modifiers, comment)
        
        self.related = set()
#        self.related_to = set()
            
    @property
    def is_annonymous(self):
        value = self.get_value("is_annonymous")
        return bool(value)
    
    def name_mangle(self, tag):
        """ Mangle tag name if it conflicts with python keyword
        Example::
            >>> term.name_mangle("def"), term.name_mangle("class")
            ('def_', 'class_')
        """
        if tag in ["def", "class", "in", "not"]:
            return tag + "_"
        else:
            return tag
        
    def name_demangle(self, tag):
        """ Reverse of name_mangle
        """
        if tag in ["def_", "class_", "in_", "not_"]:
            return tag[:-1]
        else:
            return tag
        
    def add_tag(self, tag, value, modifiers=None, comment=None):
        """ Add `tag`, `value` pair to the object with optional modifiers and
        comment.
        Example::
            >>> term = OBOObject("Term")
            >>> term.add_tag("id", "FOO:002", comment="This is an id")
            >>> print term
            [Term]
            id: FOO:002 ! This is an id
             
        """
        tag = intern(tag) # a small speed and memory benefit
        self.tag_values.append((tag, value))
        self.modifiers.append(modifiers)
        self.comments.append(comment)
        self.values.setdefault(tag, []).append(value)
        
        #  TODO: fix multiple tags grouping
        if hasattr(self, tag):
            if isinstance(getattr(self, tag), list):
                getattr(self, tag).append(value)
            else:
                setattr(self, tag, [getattr(self, tag)] + [value])
        else:
            setattr(self, self.name_mangle(tag), value)
            
    def update(self, other):
        """ Update the term with tag value pairs from `other` 
        (a OBOObject instance). The tag value pairs are appended
        to the end except for the `id` tag.
        """ 
        for (tag, value), modifiers, comment in zip(other.tag_values, other.modifiers, other.comments):
            if tag != "id":
                self.add_tag(tag, value, modifiers, comment)
        
    def get_value(self, tag, group=True):
        if group:
            pairs = [pair for pair in self.tag_values if pair[0] == tag]
            return pairs
        else:
            tag = self.name_mangle(tag)
            if tag in self.__dict__:
                return self.__dict__[tag]
            else:
                raise ValueError("No value for tag: %s" % tag)
        
    def tag_count(self):
        """ Retrun the number of tags in this object
        """
        return len(self.tag_values)
    
    def tags(self):
        """ Retrun an iterator over the (tag, value) pairs.
        """
        for i in range(self.tag_count()):
            yield self.tag_values[i] + (self.modifiers[i], self.comments[i])
        
    def format_single_tag(self, index):
        """Return a formated string representing index-th tag pair value
        Example::
            >>> term = OBOObject("Term", id="FOO:001", name="bar", def_="Example definition {modifier=frob} ! Comment")
            >>> term.format_single_tag(0)
            'id: FOO:001'
            >>> term.format_single_tag(1)
            'def: Example definition { modifier=frob } ! Comment'
        """
        tag, value = self.tag_values[index]
        modifiers = self.modifiers[index]
        comment = self.comments[index]
        res = ["%s: %s" % (tag, value)]
        if modifiers:
            res.append("{ %s }" % modifiers)
        if comment:
            res.append("! " + comment)
        return " ".join(res)
    
    def format_stanza(self):
        """ Return a string stanza representation of this object 
        """
        stanza = ["[%s]" % self.stanza_type]
        for i in range(self.tag_count()):
            stanza.append(self.format_single_tag(i))
        return "\n".join(stanza)
            
    @classmethod     
    def parse_stanza(cls, stanza):
        r""" Parse and return an OBOObject instance from a single stanza.
        Example::
            >>> term = OBOObject.parse_stanza("[Term]\nid: FOO:001\nname:bar")
            >>> print term.id, term.name
            FOO:001 bar
            
        """
        lines = stanza.splitlines()
        stanza_type = lines[0].strip("[]")
        tag_values = []
        for line in lines[1:]:
            if ":" in line:
                tag_values.append(cls.parse_tag_value(line))
        
        obo = OBOObject(stanza_type)
        for i, (tag, value, modifiers, comment) in enumerate(tag_values):
#            print tag, value, modifiers, comment
            obo.add_tag(tag, value, modifiers, comment)
        return obo
    
        
    @classmethod
    def parse_tag_value_1(cls, tag_value_pair, *args):
        """ Parse and return a four-tuple containing a tag, value, a list of modifier pairs, comment.
        If no modifiers or comments are present the corresponding entries will be None.
        
        Example::
            >>> OBOObject.parse_tag_value("foo: bar {modifier=frob} ! Comment")
            ('foo', 'bar', 'modifier=frob', 'Comment')
            >>> OBOObject.parse_tag_value("foo: bar")
            ('foo', 'bar', None, None)
            >>> #  Can also pass tag, value pair already split   
            >>> OBOObject.parse_tag_value("foo", "bar {modifier=frob} ! Comment")
            ('foo', 'bar', 'modifier=frob', 'Comment')
        """
        if args and ":" not in tag_value_pair:
            tag, rest = tag_value_pair, args[0]
        else:
            tag, rest = _split_and_strip(tag_value_pair, ":")
        value, modifiers, comment = None, None, None
        
        if "{" in rest:
            value, rest = _split_and_strip(rest, "{",)
            modifiers, rest = _split_and_strip(rest, "}")
        if "!" in rest:
            if value is None:
                value, comment = _split_and_strip(rest, "!")
            else:
                _, comment = _split_and_strip(rest, "!")
        if value is None:
            value = rest
            
        if modifiers is not None:
            modifiers = modifiers #TODO: split modifiers in a list
            
        return tag, value, modifiers, comment
    
    _RE_TAG_VALUE = re.compile(r"^(?P<tag>.+?[^\\])\s*:\s*(?P<value>.+?)\s*(?P<modifiers>[^\\]{.+?[^\\]})?\s*(?P<comment>[^\\]!.*)?$")
    _RE_VALUE = re.compile(r"^\s*(?P<value>.+?)\s*(?P<modifiers>[^\\]{.+?[^\\]})?\s*(?P<comment>[^\\]!.*)?$")
    
    @classmethod
    def parse_tag_value(cls, tag_value_pair, arg=None):
        """ Parse and return a four-tuple containing a tag, value, a list of modifier pairs, comment.
        If no modifiers or comments are present the corresponding entries will be None.
        
        Example::
            >>> OBOObject.parse_tag_value("foo: bar {modifier=frob} ! Comment")
            ('foo', 'bar', 'modifier=frob', 'Comment')
            >>> OBOObject.parse_tag_value("foo: bar")
            ('foo', 'bar', None, None)
            >>> #  Can also pass tag, value pair already split   
            >>> OBOObject.parse_tag_value("foo", "bar {modifier=frob} ! Comment")
            ('foo', 'bar', 'modifier=frob', 'Comment')
            
        .. warning: This function assumes comment an modifiers are prefixed
            with a whitespace i.e. 'tag: bla! comment' will be parsed incorrectly!
        """
        if arg is not None: # tag_value_pair is actually a tag only
            tag = tag_value_pair
            value, modifiers, comment =  cls._RE_VALUE.findall(arg)[0]
        else:
            tag, value, modifiers, comment = cls._RE_TAG_VALUE.findall(tag_value_pair)[0]
        none_if_empyt = lambda val: None if not val.strip() else val.strip()
        modifiers = modifiers.strip(" {}")
        comment = comment.lstrip(" !")
        return (none_if_empyt(tag), none_if_empyt(value),
                none_if_empyt(modifiers), none_if_empyt(comment))
         
    def related_objects(self):
        """ Return a list of tuple pairs where the first element is relationship (typedef id)
        is and the second object id whom the relationship applies to.
        """
        result = [(type_id, id) for type_id in ["is_a"] for id in self.values.get(type_id, [])] ##TODO add other defined Typedef ids
        result = result + [tuple(r.split(None, 1)) for r in self.values.get("relationship", [])]
        return result

    def __repr__(self):
        """ Return a string representation of the object in OBO format
        """
        return self.format_stanza()

    def __iter__(self):
        """ Iterates over sub terms
        """
        for type_id, id in self.related_objects():
            yield (type_id, id)
        
        
class Term(OBOObject):
    def __init__(self, *args, **kwargs):
        OBOObject.__init__(self, "Term", *args, **kwargs)

class Typedef(OBOObject):
    def __init__(self, *args, **kwargs):
        OBOObject.__init__(self, "Typedef", *args, **kwargs)

class Instance(OBOObject):
    def __init__(self, *args, **kwargs):
        OBOObject.__init__(self, "Instance", *args, **kwargs)


import re

class OBOParser(object):
    r""" A simple parser for .obo files (inspired by xml.dom.pulldom)
    
    Example::
        >>> from StringIO import StringIO
        >>> file = StringIO("header_tag: header_value\n[Term]\nid: FOO { modifier=bar } ! comment\n\n") 
        >>> parser = OBOParser(file)
        >>> for event, value in parser:
        ...     print event, value
        ...     
        HEADER_TAG ['header_tag', 'header_value']
        START_STANZA Term
        TAG_VALUE ('id', 'FOO', 'modifier=bar', 'comment')
        CLOSE_STANZA None
                
    """
    def __init__(self, file):
        self.file = file
        
    def parse(self, progress_callback=None):
        """ Parse the file and yield parse events.
        
        .. TODO: List events and values
        """
        data = self.file.read()
        header = data[: data.index("\n[")]
        body = data[data.index("\n[") + 1:]
        for line in header.splitlines():
            if line.strip():
                yield "HEADER_TAG", line.split(": ", 1)
                
        current = None
        #  For speed make these functions local
        startswith = str.startswith
        endswith = str.endswith
        parse_tag_value = OBOObject.parse_tag_value
        
        for line in body.splitlines():
#            line = line.strip()
            if startswith(line, "[") and endswith(line, "]"):
                yield "START_STANZA", line.strip("[]")
                current = line
            elif startswith(line, "!"):
                yield "COMMENT", line[1:]
            elif line:
                yield "TAG_VALUE", parse_tag_value(line)
            else: #  empty line is the end of a term
                yield "CLOSE_STANZA", None
                current = None
        if current is not None:
            yield "CLOSE_STANZA", None
    
    def __iter__(self):
        """ Iterate over parse events (same as parse())
        """
        return self.parse() 
        
        
class OBOOntology(object):
    """ The main ontology object.
    """
    
    BUILTINS = BUILTIN_OBO_OBJECTS
    def __init__(self, file=None):
        """ Init an ontology instance from a file like object (.obo format)
        
        """
        self.objects = []
        self.header_tags = []
        self.id2term = {}
        self.alt2id = {}
        self._resolved_imports = []
        self._invalid_cache_flag = False
        self._related_to = {}
        
        # First load the built in OBO objects
        builtins = StringIO("\n" + "\n\n".join(self.BUILTINS) + "\n")
        self.load(builtins)
        if file:
            self.load(file)
#     Original function - to revert remove '_OLD'   
    def add_object_OLD(self, object):
        """ Add OBOObject instance to this ontology.
        """
        if object.id in self.id2term:
            raise ValueError("OBOObject with id: %s already in the ontology" % object.id)
        self.objects.append(object)
        self.id2term[object.id] = object
        self._invalid_cache_flag = True
        
    def add_object(self, object):
#     	print "Adding object"
    	self.objects.append(object)
        self.id2term[object.id] = object
        self._invalid_cache_flag = True
    
        
    def add_header_tag(self, tag, value):
        """ Add header tag, value pair to this ontology
        """
        self.header_tags.append((tag, value))
    
    def load(self, file, progress_callback=None):
        """ Load terms from a file.
        """
        if isinstance(file, basestring):
            file = open(file, "rb")
        parser = OBOParser(file)
        current = None
        for event, value in parser.parse(progress_callback=progress_callback):
            if event == "TAG_VALUE":
                current.add_tag(*value)
            elif event == "START_STANZA":
                current = OBOObject(value)
            elif event == "CLOSE_STANZA":
                self.add_object(current)
                current = None
            elif event == "HEADER_TAG":
                self.add_header_tag(*value)
            elif event != "COMMENT":
                raise Exception("Parse Error! Unknown parse event {0}".format(event))
            
        imports = [value for tag, value, in self.header_tags if tag == "import"]
        
        while imports:
            url = imports.pop(0)
            if uri not in self._resolved_imports:
                imported = self.parse_file(open(url, "rb"))
                ontology.update(imported)
                self._resolved_imports.append(uri)
                
    def dump(self, file):
        """ Dump the contents of the ontology to a .obo `file`.
        """
        if isinstance(file, basestring):
            file = open(file, "wb")
            
        for key, value in self.header_tags:
            file.write(key + ": " + value + "\n")
            
        # Skip the builtins
        for object in self.objects[len(self.BUILTINS):]:
            file.write("\n")
            file.write(object.format_stanza())
            file.write("\n")
    
    def update(self, other):
        """ Update this ontology with the terms from another. 
        """
        for term in other:
            if term.id in self:
                if not term.is_annonymous:
                    self.term(term.id).update(term)
                else: #  Do nothing
                    pass 
            else:
                self.add_object(term)
        self._invalid_cache_flag = True
        
    def _cache_validate(self, force=False):
        """ Update the relations cache if `self._invalid_cache` flag is set. 
        """
        if self._invalid_cache_flag or force:
            self._cache_relations()
            
    def _cache_relations(self):
        """ Collect all relations from parent to a child and store it in
        `self._related_to` member.
        
        """
        related_to = defaultdict(list)
        for obj in self.objects:
            for rel_type, id in self.related_terms(obj):
                term = self.term(id)
                related_to[term].append((rel_type, obj))
                
        self._related_to = related_to
        self._invalid_cache_flag = False
        
    def term(self, id):
        """ Return the OBOObject associated with this id.
        """
        if isinstance(id, basestring):
            if id in self.id2term:
                return self.id2term[id]
            elif id in self.alt2id:
                return self.id2term[self.alt2id[id]]
            else:
                raise ValueError("Unknown term id: %r" % id)
                pass
        elif isinstance(id, OBOObject):
            return id
        
    def terms(self):
        """ Return all `Term` instances in the ontology.
        """
        return [obj for obj in self.objects if obj.stanza_type == "Term"]
    
    def term_by_name(self, name):
        """ Return the term with name ``name``. 
        """
        terms = [t for t in self.terms() if t.name == name]
        if len(terms) != 1:
            raise ValueError("Unknown term name: %r" % name)
        return terms[0]
    
    def typedefs(self):
        """ Return all `Typedef` instances in the ontology.
        """
        return [obj for obj in self.objects if obj.stanza_type == "Typedef"]
    
    def instances(self):
        """ Return all `Instance` instances in the ontology.
        """
        return [obj for obj in self.objects if obj.stanza_type == "Instance"]
        
    def related_terms(self, term):
        """ Return a list of (rel_type, term_id) tuples where rel_type is
        relationship type (e.g. 'is_a', 'has_part', ...) and term_id is the
        id of the term in the relationship.
        
        """
        term = self.term(term) if not isinstance(term, OBOObject) else term
        related = [(tag, value) for tag in ["is_a"] for value in term.values.get(tag, [])] #TODO: add other typedef ids
        relationships = term.values.get("relationship", [])
        for rel in relationships:
            related.append(tuple(rel.split(None, 1)))
        return related
        
    def edge_types(self):
        """ Return a list of all edge types in the ontology
        """
        return [obj.id for obj in self.objects if obj.stanza_type == "Typedef"]
    
    def parent_edges(self, term):
        """ Return a list of (rel_type, parent_term) tuples 
        """
        term = self.term(term)
        parents = []
        for rel_type, parent in self.related_terms(term):
            parents.append((rel_type, self.term(parent)))
        return parents
        
    def child_edges(self, term):
        """ Return a list of (rel_type, source_term) tuples
        """
        self._cache_validate()
        term = self.term(term)
        return self._related_to.get(term, [])
        
        
    def super_terms(self, term):
        """ Return a set of all super terms of `term` up to the most general one.
        """
        terms = self.parent_terms(term)
        visited = set()
        queue = set(terms)
        while queue:
            term = queue.pop()
            visited.add(term)
            queue.update(self.parent_terms(term) - visited)
        return visited
    
    def sub_terms(self, term):
        """ Return a set of all sub terms for `term`.
        """
        terms = self.child_terms(term)
        visited = set()
        queue = set(terms)
        while queue:
            term = queue.pop()
            visited.add(term)
            queue.update(self.child_terms(term) - visited)
        return visited
    
    def child_terms(self, term):
        """ Return a set of all child terms for this `term`.
        """
        self._cache_validate()
        term = self.term(term)
        children = []
        for rel_type, term in self.child_edges(term):
            children.append(term)
        return set(children)
        
    def parent_terms(self, term):
        """ Return a set of all parent terms for this `term`
        """
        term = self.term(term)
        parents = []
        for rel_type, id in self.parent_edges(term): #term.related_objects():
            parents.append(self.term(id))
        return set(parents)
    # Ken
    # Find the root of a given ontology
    def root(self):
        try:
            return [t for t in self.terms() if self.parent_terms(t) == set([]) and self.child_terms(t)!=set([])][0]
        except:
            raise("Can't find the root of the ontology. Please check the intergrity of the file!")

    #  Check if  a term is the root term
    def is_root(self, term):
        return self.root() == self.term(term)

    def relations(self):
        """ Return a list of all relations in the ontology.
        """
        relations = []
        for obj in self.objects:
            for type_id, id in  obj.related:
                target_term = self.term(id)
            relations.append((obj, type_id, target_term))
        return relations
    
    def __len__(self):
        return len(self.objects)
    
    def __iter__(self):
        return iter(self.objects)
    
    def __contains__(self, obj):
        if isinstance(obj, basestring):
            return obj in self.id2term
        else:
            return obj in self.objects
    
    def __getitem__(self, key):
        return self.id2term[key]
    
    def has_key(self, key):
        return self.id2term.has_key(key)

    # Ken
    def isRoot(self, term):
        return  self.parent_terms(term)== set([])

    # Calculate the topological characteristics of a term

    # Kenneth added function
    def topology(self, term):
        if self.is_root(term):
            return  1.0
        num_children_parents= [self.child_terms(i).__len__() for i in list(self.climb_bf(term))]
        return  1.0/reduce(lambda x,y: x * y, num_children_parents)

    def mca(self, term, other):
        anc_term1= [i.id for i in self.super_terms(term)]
        anc_term2= [i.id for i in self.super_terms(other)]

        set(anc_term1).intersection(set(anc_term2))

        # print sorted([i.id for i in anc_term1])
        # print sorted([i.id for i in anc_term2])
        return (self.topology(term), self.topology(other))

        # return zip([i.id for i in anc_term1], [i.id for i in anc_term2])

    def getGlobalPairwiseSimilarity(self, term):
        # for t in self.terms():
        #     print self.mca(t.id, term)
        return [self.mca(i.id, term) for i in self.terms() ]


        # return len(self.terms())

    def climb_bf(self, term):
        queue = list(self.parent_terms(term))
        while queue:
            term = queue.pop(0)
            queue.extend(self.parent_terms(term))
            yield  term
    def traverse_bf(self, term):
        """ BF traverse of the ontology down from term.
        """
        queue = list(self.child_terms(term))
        while queue:
            term = queue.pop(0)
            queue.extend(self.child_terms(term))
            yield term
                    
    def traverse_df(self, term, depth=1e30):
        """ DF traverse of the ontology down from term.
        """
        if depth >= 1:
            for child in self.child_terms(term):
                yield child
                for t in self.traverse_df(child, depth-1):
                    yield t
        
    
    def to_network(self, terms=None):
        """ Return an Orange.network.Network instance constructed from
        this ontology.
        
        """
        edge_types = self.edge_types()
        terms = self.terms()
        from Orange.orng import orngNetwork
        import orange
        
        network = orngNetwork.Network(len(terms), True, len(edge_types))
        network.objects = dict([(term.id, i) for i, term in enumerate(terms)])
        
        edges = defaultdict(set)
        for term in self.terms():
            related = self.related_terms(term)
            for relType, relTerm in related:
                edges[(term.id, relTerm)].add(relType)
                
        edgeitems = edges.items()
        for (src, dst), eTypes in edgeitems:
            network[src, dst] = [1 if e in eTypes else 0 for e in edge_types]
            
        domain = orange.Domain([orange.StringVariable("id"),
                                orange.StringVariable("name"),
                                orange.StringVariable("def"),
                                ], False)
        
        items = orange.ExampleTable(domain)
        for term in terms:
            ex = orange.Example(domain, [term.id, term.name, term.values.get("def", [""])[0]])
            items.append(ex)
        
        relationships = set([", ".join(sorted(eTypes)) for (_, _), eTypes in edgeitems])
        domain = orange.Domain([orange.FloatVariable("u"),
                                orange.FloatVariable("v"),
                                orange.EnumVariable("relationship", values=list(edge_types))
                                ], False)
        
        id2index = dict([(term.id, i + 1) for i, term in enumerate(terms)])
        links = orange.ExampleTable(domain)
        for (src, dst), eTypes in edgeitems:
            ex = orange.Example(domain, [id2index[src], id2index[dst], eTypes.pop()])
            links.append(ex)
            
        network.items = items
        network.links = links
        network.optimization = None
        return network
    
    def to_networkx(self, terms=None):
        """ Return a NetworkX graph of this ontology
        """
        import networkx
        graph = networkx.Graph()
        
        edge_types = self.edge_types()
        
        edge_colors = {"is_a": "red"}
        
        if terms is None:
            terms = self.terms()
        else:
            terms = [self.term(term) for term in terms]
            super_terms = [self.super_terms(term) for term in terms]
            terms = reduce(set.union, super_terms, set(terms))
            
        for term in terms:
            graph.add_node(term.id, name=term.name)
            
        for term in terms:
            for rel_type, rel_term in self.related_terms(term):
                rel_term = self.term(rel_term)
                if rel_term in terms:
                    graph.add_edge(term.id, rel_term.id, label=rel_type, color=edge_colors.get(rel_type, "blue"))
                    
        return graph
    
    def to_graphviz(self, terms=None):
        """ Return an pygraphviz.AGraph representation of the ontology in.
        If `terms` is not `None` it must be a list of terms in the ontology.
        The graph will in this case contain only the super graph of those
        terms.  
        
        """
        import pygraphviz as pgv
        graph = pgv.AGraph(directed=True, name="ontology")
        
        edge_types = self.edge_types()
        
        edge_colors = {"is_a": "red"}
        
        if terms is None:
            terms = self.terms()
        else:
            terms = [self.term(term) for term in terms]
            super_terms = [self.super_terms(term) for term in terms]
            terms = reduce(set.union, super_terms, set(terms))
            
        for term in terms:
            graph.add_node(term.id, name=term.name)
            
        for term in terms:
            for rel_type, rel_term in self.related_terms(term):
                rel_term = self.term(rel_term)
                if rel_term in terms:
                    graph.add_edge(term.id, rel_term.id, label=rel_type, color=edge_colors.get(rel_type, "blue"))
                    
        return graph
    
    
def load(file):
    """ Load an ontology from a .obo file
    """
    return OBOOntology(file)
    
    
def foundry_ontologies():
    """ List ontologies available from the OBOFoundry website
    (`http://www.obofoundry.org/`_) 
    Example::
        >>> foundry_ontologies()
        [('Biological process', 'http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/genomic-proteomic/gene_ontology_edit.obo'), ...
    
    """
    import urllib2, re
    stream = urllib2.urlopen("http://www.obofoundry.org/")
    text = stream.read()
    pattern = r'<td class=".+?">\s*<a href=".+?">(.+?)</a>\s*</td>\s*<td class=".+?">.*?</td>\s*<td class=".+?">.*?</td>\s*?<td class=".+?">\s*<a href="(.+?obo)">.+?</a>'
    return re.findall(pattern, text)
    
    
if __name__ == "__main__":
    import doctest
    stanza = '''[Term]
id: FOO:001
name: bar
'''
    from StringIO import StringIO
    seinfeld = StringIO("""
[Typedef]
id: parent

[Typedef]
id: child
inverse_of: parent ! not actually used yet

[Term]
id: 001
name: George

[Term]
id: 002
name: Estelle
relationship: parent 001 ! George

[Term]
id: 003
name: Frank
relationship: parent 001 ! George

""") # TODO: fill the ontology with all characters
    term = OBOObject.parse_stanza(stanza)
    
    seinfeld = OBOOntology(seinfeld)
    print seinfeld.child_edges("001")
    
    doctest.testmod(extraglobs={"stanza": stanza, "term": term}, optionflags=doctest.ELLIPSIS)


