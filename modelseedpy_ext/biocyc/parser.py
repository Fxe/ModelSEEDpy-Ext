import lxml.etree as et
from enum import Enum


class XmlParseMethod(Enum):
    UNIQUE = 'singleton'
    UNIQUE_E = 'singleton element only'
    LIST = 'list append'


def _parse_taxonomic_range(elem, parser):
    for action, elem in parser:
        tag = elem.tag
        # print('_parse_taxonomic_range', action, tag)
        if action == 'start' and tag == 'Pathway':
            # p, b = _capture_single(elem, parser, 'Pathway')
            # print(p, b)
            pass
        elif action == 'end' and tag == 'taxonomic-range':
            return


def _parse_reaction_ordering(elem, parser):
    for action, elem in parser:
        tag = elem.tag
        # print('_parse_reaction_ordering', action, tag)
        if action == 'start' and tag == 'Pathway':
            # p, b = _capture_single(elem, parser, 'Pathway')
            # print(p, b)
            pass
        elif action == 'end' and tag == 'reaction-ordering':
            return


def _parse_reaction_layout(elem, parser):
    for action, elem in parser:
        tag = elem.tag
        # print('_parse_reaction_ordering', action, tag)
        if action == 'start' and tag == 'Pathway':
            # p, b = _capture_single(elem, parser, 'Pathway')
            # print(p, b)
            pass
        elif action == 'end' and tag == 'reaction-layout':
            return


def _parse_reaction_list(elem, parser):
    res = []
    for action, elem in parser:
        tag = elem.tag
        # print('_parse_reaction_list', action, tag)
        if action == 'start' and tag == 'Reaction':
            p, b = _capture_single(elem, parser, 'Reaction')
            res.append((p, b))
        elif action == 'end' and tag == 'reaction-list':
            return res


def _parse_reaction_species(elem, parser):
    for action, elem in parser:
        tag = elem.tag
        # print('_parse_reaction_species', action, tag)
        if action == 'start' and tag == 'Pathway':
            # p, b = _capture_single(elem, parser, 'Pathway')
            # print(p, b)
            pass
        elif action == 'end' and tag == 'species':
            return


def _parse_parent(elem, parser):
    p = None
    b = None
    for action, elem in parser:
        tag = elem.tag
        print('_parse_parent', action, tag)
        if action == 'start' and tag == 'Pathway':
            p, b = _capture_single(elem, parser, 'Pathway')
        elif action == 'end' and tag == 'parent':
            return p, b


class BiocycPathwayParser:

    def __init__(self):
        self.tag_capture = {
            ('start', 'parent'): _parse_parent,
            ('start', 'taxonomic-range'): _parse_taxonomic_range,
            ('start', 'reaction-ordering'): _parse_reaction_ordering,
            ('start', 'species'): _parse_reaction_species,
            ('start', 'reaction-layout'): _parse_reaction_layout,
            ('start', 'reaction-list'): _parse_reaction_list,
        }

    def _parse(self, elem, parser):
        res = {
            'parent': []
        }
        for action, elem in parser:
            tag = elem.tag
            t = (action, tag)
            # print('_parse', action, tag)
            if t in self.tag_capture:
                if tag == 'parent':
                    res['parent'].append(self.tag_capture[t](elem, parser))
                else:
                    res[t[1]] = self.tag_capture[t](elem, parser)
        return res

    def parse(self, fh):
        parser = et.iterparse(fh, events=("end", "start"))
        for action, elem in parser:
            tag = elem.tag
            if action == "start" and tag == 'Pathway':
                return self._parse(elem, parser)


def _capture_single(elem, parser, end_tag):
    b = None
    p = {}
    for action, elem in parser:
        tag = elem.tag
        p = dict(elem.attrib)
        b = elem.text
        if action == 'end' and tag == end_tag:
            return p, b
        else:
            raise Exception('invalid single')


class BiocycParserSingle:

    def __init__(self, single_tag, end_tag):
        self.single_tag = single_tag
        self.end_tag = end_tag

    def parse(self, elem, parser):
        p = None
        b = None
        for action, elem in parser:
            tag = elem.tag
            # print('BiocycParserSingle', self.single_tag, action, tag)
            if action == 'start' and tag == self.single_tag:
                p, b = _capture_single(elem, parser, self.single_tag)
            elif action == 'end' and tag == self.end_tag:
                return p, b
            else:
                raise Exception(f'bad traversal expected: [{self.single_tag}], found [{tag}]')


class BiocycParseBase:

    def __init__(self, end_tag):
        self.end_tag = end_tag
        self.tag_capture = {}

    def _parse(self, elem, parser):
        body = dict(elem.attrib)
        res = {}
        for action, elem in parser:
            t = (action, elem.tag)
            # print('BiocycParseEnzymaticReaction', t)
            if t in self.tag_capture:
                p, method = self.tag_capture[t]
                d = p.parse(elem, parser)
                attr = None
                if attr is None:
                    attr = t[1]
                if method == XmlParseMethod.UNIQUE:
                    if attr in res:
                        raise Exception(f'non unique capture {attr}')
                    else:
                        res[attr] = d
                elif method == XmlParseMethod.UNIQUE_E:
                    if attr in res:
                        raise Exception(f'non unique capture {attr}')
                    else:
                        res[attr] = (dict(elem.attrib), elem.text)
                elif method == XmlParseMethod.LIST:
                    if attr not in res:
                        res[attr] = []
                    res[attr].append(d)
                else:
                    raise Exception(f'method not implemented {method}')
            elif action == 'end' and elem.tag == self.end_tag:
                return body, res
            else:
                raise Exception(f'bad traversal expected: {self.tag_capture.keys()}, found {t}')


class BiocycParseEnzymaticReaction(BiocycParseBase):

    def __init__(self):
        super().__init__('Enzymatic-Reaction')
        self.tag_capture = {
            ('start', 'enzyme'): (BiocycParserSingle('Protein', 'enzyme'), XmlParseMethod.UNIQUE),
            ('start', 'reaction'): (BiocycParserSingle('Reaction', 'reaction'), XmlParseMethod.UNIQUE),
            ('start', 'common-name'): (BiocycParserSingle('common-name', 'common-name'), XmlParseMethod.UNIQUE_E),
        }


class BiocycParseEnzymaticReactions(BiocycParseBase):

    def __init__(self):
        super().__init__('enzymatic-reaction')
        self.tag_capture = {
            ('start', 'Enzymatic-Reaction'): (BiocycParseEnzymaticReaction(), XmlParseMethod.LIST),
            ('end', 'Enzymatic-Reaction'): None,
        }


class BiocycParserDblink:

    def __init__(self):
        self.end_tag = 'dblink'
        self.capture = {'dblink-db', 'dblink-oid', 'dblink-relationship', 'dblink-URL'}

    def parse(self, elem, parser):
        dblink = {}
        for action, elem in parser:
            tag = elem.tag
            # print('BiocycParserDblink', action, tag, tag in self.capture)
            if action == 'start' and tag in self.capture:
                dblink[tag] = elem.text
            elif action == 'end' and tag in self.capture:
                pass
            elif action == 'end' and tag == self.end_tag:
                return dblink
            else:
                raise Exception(f'bad traversal expected: [{self.capture}], found [{tag}] {action}')


class BiocycParserPass:

    def __init__(self, end_tag):
        self.end_tag = end_tag

    def parse(self, elem, parser):
        for action, elem in parser:
            tag = elem.tag
            if action == 'end' and tag == self.end_tag:
                return None


class BiocycReactionParser:

    def __init__(self):
        self.tag_capture = {
            ('start', 'parent'): BiocycParserSingle('Reaction', 'parent'),
            ('start', 'left'): BiocycParserSingle('Compound', 'left'),
            ('start', 'right'): BiocycParserSingle('Compound', 'right'),
            ('start', 'dblink'): BiocycParserDblink(),
            ('start', 'enzymatic-reaction'): BiocycParseEnzymaticReactions(),
            # ('start', 'reaction-ordering'): _parse_reaction_ordering,
            # ('start', 'species'): _parse_reaction_species,
            # ('start', 'reaction-layout'): _parse_reaction_layout,
            # ('start', 'reaction-list'): _parse_reaction_list,
        }

    def _parse(self, elem, parser):
        res = {
            'parent': [],
            'left': [],
            'right': [],
            'dblink': [],
        }
        for action, elem in parser:
            tag = elem.tag
            t = (action, tag)
            # print('_parse', action, tag)
            if t in self.tag_capture:
                if tag == 'parent':
                    res['parent'].append(self.tag_capture[t].parse(elem, parser))
                elif tag == 'left':
                    res['left'].append(self.tag_capture[t].parse(elem, parser))
                elif tag == 'right':
                    res['right'].append(self.tag_capture[t].parse(elem, parser))
                elif tag == 'dblink':
                    res['dblink'].append(self.tag_capture[t].parse(elem, parser))
                else:
                    res[t[1]] = self.tag_capture[t].parse(elem, parser)
        return res

    def parse(self, fh):
        parser = et.iterparse(fh, events=("end", "start"))
        for action, elem in parser:
            tag = elem.tag
            if action == "start" and tag == 'Reaction':
                return self._parse(elem, parser)


class BiocycProteinParser:

    def __init__(self):
        self.end_tag = 'Protein'
        self.tag_capture = {
            ('start', 'gene'): (BiocycParserSingle('Gene', 'gene'), XmlParseMethod.UNIQUE),
            ('start', 'parent'): (BiocycParserSingle('Protein', 'parent'), XmlParseMethod.UNIQUE),
            ('start', 'synonym'): (BiocycParserSingle('synonym', 'synonym'), XmlParseMethod.UNIQUE),
            ('start', 'common-name'): (BiocycParserSingle('common-name', 'common-name'), XmlParseMethod.UNIQUE),
            ('start', 'molecular-weight-seq'): (
                BiocycParserSingle('molecular-weight-seq', 'molecular-weight-seq'), XmlParseMethod.UNIQUE),
            ('start', 'species'): (BiocycParserPass('species'), XmlParseMethod.UNIQUE),
            ('start', 'catalyzes'): (BiocycParserPass('catalyzes'), XmlParseMethod.UNIQUE),
            ('start', 'credits'): (BiocycParserPass('credits'), XmlParseMethod.UNIQUE),
            ('start', 'comment'): (BiocycParserPass('comment'), XmlParseMethod.UNIQUE),
            ('start', 'dblink'): (BiocycParserDblink(), XmlParseMethod.LIST),
        }

    def _parse(self, elem, parser):
        res = {}
        for action, elem in parser:
            t = (action, elem.tag)
            print('BiocycProteinParser', t)
            if t in self.tag_capture:
                p, method = self.tag_capture[t]
                d = p.parse(elem, parser)
                attr = None
                if attr == None:
                    attr = t[1]
                if method == XmlParseMethod.UNIQUE:
                    if attr in res:
                        raise Exception(f'non unique capture {attr}')
                    else:
                        res[attr] = d
                elif method == XmlParseMethod.UNIQUE_E:
                    if attr in res:
                        raise Exception(f'non unique capture {attr}')
                    else:
                        res[attr] = (dict(elem.attrib), elem.text)
                elif method == XmlParseMethod.LIST:
                    if attr not in res:
                        res[attr] = []
                    res[attr].append(d)
                else:
                    raise Exception(f'method not implemented {method}')
            elif action == 'end' and elem.tag == self.end_tag:
                return res
            else:
                raise Exception(f'bad traversal expected: {self.tag_capture.keys()}, found {t}')

    def parse(self, fh):
        parser = et.iterparse(fh, events=("end", "start"))
        for action, elem in parser:
            tag = elem.tag
            if action == "start" and tag == 'Protein':
                return self._parse(elem, parser)


class BiocycParseTag:

    def __init__(self, tag):
        self.tag = tag

    def parse(self, elem, parser):
        attr = dict(elem.attrib)
        text = elem.text
        for action, elem in parser:
            if action == 'end' and elem.tag == self.tag:
                return attr, text
            else:
                raise Exception(f'Element [{self.tag}] contains sub element [{elem.tag}]')


class BiocycXmlParser(BiocycParseBase):

    def __init__(self, tag):
        super().__init__(tag)
        self.tag = tag

    def parse(self, h):
        parser = et.iterparse(h, events=("end", "start"))
        for action, elem in parser:
            tag = elem.tag
            if action == "start" and elem.tag == self.tag:
                return super()._parse(elem, parser)


class BiocycCompoundParser(BiocycXmlParser):

    def __init__(self):
        super().__init__('Compound')
        self.tag_capture = {
            ('start', 'parent'): (BiocycParserSingle('Compound', 'parent'), XmlParseMethod.LIST),
            ('start', 'cml'): (BiocycParserPass('cml'), XmlParseMethod.UNIQUE),
            ('start', 'appears-in-right-side-of'): (
                BiocycParserPass('appears-in-right-side-of'), XmlParseMethod.UNIQUE),
            ('start', 'inchi'): (BiocycParseTag('inchi'), XmlParseMethod.UNIQUE),
            ('start', 'synonym'): (BiocycParseTag('synonym'), XmlParseMethod.LIST),
            ('start', 'molecular-weight'): (BiocycParseTag('molecular-weight'), XmlParseMethod.UNIQUE),
            ('start', 'monoisotopic-mw'): (BiocycParseTag('monoisotopic-mw'), XmlParseMethod.UNIQUE),
            ('start', 'dblink'): (BiocycParserDblink(), XmlParseMethod.LIST),
            ('start', 'inchi-key'): (BiocycParserPass('inchi-key'), XmlParseMethod.UNIQUE),
            ('start', 'common-name'): (BiocycParseTag('common-name'), XmlParseMethod.UNIQUE),
            ('start', 'gibbs-0'): (BiocycParserPass('gibbs-0'), XmlParseMethod.UNIQUE),
            ('start', 'credits'): (BiocycParserPass('credits'), XmlParseMethod.UNIQUE),
            ('start', 'comment'): (BiocycParserPass('comment'), XmlParseMethod.UNIQUE),
            ('start', 'appears-in-left-side-of'): (BiocycParserPass('appears-in-left-side-of'), XmlParseMethod.UNIQUE),
        }
