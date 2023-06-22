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


class BiocycParseReactionEnzymaticReaction(BiocycParseBase):

    def __init__(self):
        super().__init__('Enzymatic-Reaction')
        self.tag_capture = {
            ('start', 'enzyme'): (BiocycParseSubElements('enzyme', ['Protein'], True), XmlParseMethod.UNIQUE),
            ('start', 'reaction'): (BiocycParseSubElements('reaction', ['Reaction'], True), XmlParseMethod.UNIQUE),
            ('start', 'common-name'): (BiocycParseTag('common-name'), XmlParseMethod.UNIQUE),
            ('start', 'synonym'): (BiocycParseTag('synonym'), XmlParseMethod.UNIQUE),
        }

    def parse(self, elem, parser):
        return self._parse(elem, parser)


class BiocycParseReactionEnzymaticReactions(BiocycParseBase):

    def __init__(self):
        super().__init__('enzymatic-reaction')
        self.tag_capture = {
            ('start', 'Enzymatic-Reaction'): (BiocycParseReactionEnzymaticReaction(), XmlParseMethod.LIST),
        }

    def parse(self, elem, parser):
        return self._parse(elem, parser)


class BiocycParseSubElements:

    def __init__(self, end_tag, tags, single: bool):
        self.tags = tags
        self.end_tag = end_tag
        self.single = single

    def parse(self, elem, parser):
        attr = dict(elem.attrib)
        text = elem.text
        sub_elements = []
        for action, elem in parser:
            if action == 'end' and elem.tag == self.end_tag:
                return attr, text, sub_elements
            elif action == 'start' and elem.tag in self.tags:
                attr, text = BiocycParseTag(elem.tag).parse(elem, parser)
                if self.single and len(sub_elements) > 0:
                    raise Exception('multiple sub elements')
                else:
                    sub_elements.append((attr, text))
            elif action == 'end' and elem.tag in self.tags:
                pass
            else:
                raise Exception(f'Bad element [{elem.tag}] expected either {self.tags}')


class BiocycParseReactionLeftRight:

    def parse(self, elem, parser):
        species = None
        compartment = None
        coefficient = None
        if elem.tag != 'left' and elem.tag != 'right':
            raise Exception(f'Invalid element [{elem.tag}] expected either [left, right]')

        for action, elem in parser:
            if action == 'end' and elem.tag == 'left' or elem.tag == 'right':
                return species, compartment, coefficient
            elif action == 'start' and elem.tag in ['Compound', 'Protein']:
                if species is None:
                    attr, text = BiocycParseTag(elem.tag).parse(elem, parser)
                    species = attr
                else:
                    raise Exception(f'Found multiple species')
            elif action == 'start' and elem.tag == 'compartment':
                if compartment is None:
                    attr, text, sub_elements = BiocycParseSubElements('compartment', ['cco'], True).parse(elem, parser)
                    compartment = sub_elements[0]
                else:
                    raise Exception(f'Found multiple compartment')
            elif action == 'start' and elem.tag == 'coefficient':
                if coefficient is None:
                    attr, text = BiocycParseTag(elem.tag).parse(elem, parser)
                    coefficient = text
                else:
                    raise Exception(f'Found multiple coefficient')
            elif action == 'end':
                pass
            else:
                raise Exception(f'Bad element [{elem.tag}]')


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


class BiocycReactionParser(BiocycXmlParser):

    def __init__(self):
        super().__init__('Reaction')
        self.tag_capture = {
            ('start', 'parent'): (BiocycParserSingle('Reaction', 'parent'), XmlParseMethod.LIST),
            ('start', 'left'): (BiocycParseReactionLeftRight(), XmlParseMethod.LIST),
            ('start', 'right'): (BiocycParseReactionLeftRight(), XmlParseMethod.LIST),
            ('start', 'reaction-ordering'): (BiocycParserPass('reaction-ordering'), XmlParseMethod.UNIQUE),
            ('start', 'enzymes-not-used'): (BiocycParserPass('enzymes-not-used'), XmlParseMethod.UNIQUE),
            ('start', 'enzymatic-reaction'): (BiocycParseReactionEnzymaticReactions(), XmlParseMethod.UNIQUE),
            ('start', 'dblink'): (BiocycParserDblink(), XmlParseMethod.LIST),
            ('start', 'ec-number'): (BiocycParserPass('ec-number'), XmlParseMethod.UNIQUE),
            ('start', 'reaction-list'): (BiocycParserPass('reaction-list'), XmlParseMethod.UNIQUE),
            ('start', 'reaction-direction'): (BiocycParserPass('reaction-direction'), XmlParseMethod.UNIQUE),
            ('start', 'orphan'): (BiocycParserPass('orphan'), XmlParseMethod.UNIQUE),
            ('start', 'physiologically-relevant'): (
                BiocycParserPass('physiologically-relevant'), XmlParseMethod.UNIQUE),
            ('start', 'credits'): (BiocycParserPass('credits'), XmlParseMethod.UNIQUE),
            ('start', 'in-pathway'): (BiocycParserPass('in-pathway'), XmlParseMethod.UNIQUE),
            ('start', 'gibbs-0'): (BiocycParserPass('gibbs-0'), XmlParseMethod.UNIQUE),
            ('start', 'comment'): (BiocycParserPass('comment'), XmlParseMethod.UNIQUE),
            ('start', 'citation'): (BiocycParserPass('citation'), XmlParseMethod.LIST),
        }
