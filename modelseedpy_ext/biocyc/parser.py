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
                        raise Exception(f'{self.end_tag} non unique capture {attr}')
                    else:
                        res[attr] = d
                elif method == XmlParseMethod.UNIQUE_E:
                    if attr in res:
                        raise Exception(f'{self.end_tag} non unique capture {attr}')
                    else:
                        res[attr] = (dict(elem.attrib), elem.text)
                elif method == XmlParseMethod.LIST:
                    if attr not in res:
                        res[attr] = []
                    res[attr].append(d)
                else:
                    raise Exception(f'{self.end_tag} method not implemented {method}')
            elif action == 'end' and elem.tag == self.end_tag:
                return body, res
            else:
                raise Exception(f'{self.end_tag} bad traversal expected: {self.tag_capture.keys()}, found {t}')


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
        count = 0
        for action, elem in parser:
            tag = elem.tag
            if action == 'end' and tag == self.end_tag:
                if count == 0:
                    return None
                else:
                    count -= 1
            elif action == 'start' and tag == self.end_tag:
                count += 1


class BiocycParseReactionEnzymaticReaction(BiocycParseBase):

    def __init__(self):
        super().__init__('Enzymatic-Reaction')
        self.tag_capture = {
            ('start', 'enzyme'): (BiocycParseSubElements('enzyme', ['Protein'], True), XmlParseMethod.UNIQUE),
            ('start', 'reaction'): (BiocycParseSubElements('reaction', ['Reaction'], True), XmlParseMethod.UNIQUE),
            ('start', 'common-name'): (BiocycParseTag('common-name'), XmlParseMethod.UNIQUE),
            ('start', 'synonym'): (BiocycParseTag('synonym'), XmlParseMethod.LIST),
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
                    raise Exception(f'{self.end_tag} multiple sub elements')
                else:
                    sub_elements.append((attr, text))
            elif action == 'end' and elem.tag in self.tags:
                pass
            else:
                raise Exception(f'Bad element [{elem.tag}] expected either {self.tags}')


class BiocycParseReactionLeftRight:

    def __init__(self):
        self.capture = {'coefficient', 'name-slot', 'n-name', 'n-1-name'}

    def parse(self, elem, parser):
        species = None
        compartment = None
        capture_data = {}
        if elem.tag != 'left' and elem.tag != 'right':
            raise Exception(f'Invalid element [{elem.tag}] expected either [left, right]')

        for action, elem in parser:
            if action == 'end' and elem.tag == 'left' or elem.tag == 'right':
                return {
                    'species': species,
                    'compartment': compartment,
                    'attributes': capture_data
                }
            elif action == 'start' and elem.tag in ['Compound', 'Protein', 'RNA']:
                if species is None:
                    attr, text = BiocycParseTag(elem.tag).parse(elem, parser)
                    species = attr
                else:
                    raise Exception(f'Found multiple species')
            elif action == 'start' and elem.tag in self.capture:
                if elem.tag not in capture_data:
                    attr, text = BiocycParseTag(elem.tag).parse(elem, parser)
                    capture_data[elem.tag] = text
                else:
                    raise Exception(f'Found multiple {elem.tag}')
            elif action == 'start' and elem.tag == 'compartment':
                if compartment is None:
                    attr, text, sub_elements = BiocycParseSubElements('compartment', ['cco'], True).parse(elem, parser)
                    compartment = sub_elements[0]
                else:
                    raise Exception(f'Found multiple compartment')
            elif action == 'end':
                pass
            else:
                raise Exception(f'Bad element [{elem.tag}]')


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
            if action == "start" and elem.tag == self.tag:
                return super()._parse(elem, parser)


class BiocycCompoundParser(BiocycXmlParser):

    def __init__(self):
        super().__init__('Compound')
        self.tag_capture = {
            ('start', 'subclass'): (BiocycParseSubElements('subclass',
                                                           ['Compound', 'Protein'], True), XmlParseMethod.LIST),
            ('start', 'instance'): (BiocycParseSubElements('instance',
                                                           ['Compound', 'Protein'], True), XmlParseMethod.LIST),
            ('start', 'parent'): (BiocycParserSingle('Compound', 'parent'), XmlParseMethod.LIST),
            ('start', 'cml'): (BiocycParserPass('cml'), XmlParseMethod.UNIQUE),  # want
            ('start', 'appears-in-left-side-of'): (
                BiocycParserPass('appears-in-left-side-of'), XmlParseMethod.UNIQUE),
            ('start', 'appears-in-right-side-of'): (
                BiocycParserPass('appears-in-right-side-of'), XmlParseMethod.UNIQUE),
            ('start', 'inchi'): (BiocycParseTag('inchi'), XmlParseMethod.UNIQUE),
            ('start', 'synonym'): (BiocycParseTag('synonym'), XmlParseMethod.LIST),
            ('start', 'molecular-weight'): (BiocycParseTag('molecular-weight'), XmlParseMethod.UNIQUE),
            ('start', 'monoisotopic-mw'): (BiocycParseTag('monoisotopic-mw'), XmlParseMethod.UNIQUE),
            ('start', 'dblink'): (BiocycParserDblink(), XmlParseMethod.LIST),
            ('start', 'inchi-key'): (BiocycParseTag('inchi-key'), XmlParseMethod.UNIQUE),

            ('start', 'abbrev-name'): (BiocycParseTag('abbrev-name'), XmlParseMethod.UNIQUE),
            ('start', 'common-name'): (BiocycParseTag('common-name'), XmlParseMethod.UNIQUE),
            ('start', 'n-1-name'): (BiocycParseTag('n-1-name'), XmlParseMethod.UNIQUE),
            ('start', 'n-name'): (BiocycParseTag('n-name'), XmlParseMethod.UNIQUE),
            ('start', 'n-plus-1-name'): (BiocycParseTag('n-plus-1-name'), XmlParseMethod.UNIQUE),
            ('start', 'systematic-name'): (BiocycParseTag('systematic-name'), XmlParseMethod.UNIQUE),
            ('start', 'pka1'): (BiocycParseTag('pka1'), XmlParseMethod.UNIQUE),
            ('start', 'pka2'): (BiocycParseTag('pka2'), XmlParseMethod.UNIQUE),
            ('start', 'pka3'): (BiocycParseTag('pka3'), XmlParseMethod.UNIQUE),

            ('start', 'gibbs-0'): (BiocycParserPass('gibbs-0'), XmlParseMethod.UNIQUE),
            ('start', 'regulates'): (BiocycParserPass('regulates'), XmlParseMethod.UNIQUE),  # want
            ('start', 'species'): (BiocycParserPass('species'), XmlParseMethod.LIST),
            ('start', 'cofactors-of'): (BiocycParserPass('cofactors-of'), XmlParseMethod.UNIQUE),
            ('start', 'credits'): (BiocycParserPass('credits'), XmlParseMethod.UNIQUE),
            ('start', 'comment'): (BiocycParserPass('comment'), XmlParseMethod.LIST),
            ('start', 'citation'): (BiocycParserPass('citation'), XmlParseMethod.LIST),
            ('start', 'component'): (BiocycParserPass('component'), XmlParseMethod.LIST),
            ('start', 'component-of'): (BiocycParseSubElements('component-of',
                                                               ['Compound', 'Protein'], False), XmlParseMethod.UNIQUE),

        }


class BiocycReactionParser(BiocycXmlParser):

    def __init__(self):
        super().__init__('Reaction')
        self.tag_capture = {
            ('start', 'subclass'): (BiocycParserSingle('Reaction', 'subclass'), XmlParseMethod.LIST),
            ('start', 'parent'): (BiocycParserSingle('Reaction', 'parent'), XmlParseMethod.LIST),
            ('start', 'instance'): (BiocycParserSingle('Reaction', 'instance'), XmlParseMethod.LIST),

            ('start', 'left'): (BiocycParseReactionLeftRight(), XmlParseMethod.LIST),
            ('start', 'right'): (BiocycParseReactionLeftRight(), XmlParseMethod.LIST),
            ('start', 'reaction-ordering'): (BiocycParserPass('reaction-ordering'), XmlParseMethod.UNIQUE),
            ('start', 'enzymes-not-used'): (BiocycParserPass('enzymes-not-used'), XmlParseMethod.UNIQUE),
            ('start', 'enzymatic-reaction'): (BiocycParseReactionEnzymaticReactions(), XmlParseMethod.UNIQUE),
            ('start', 'dblink'): (BiocycParserDblink(), XmlParseMethod.LIST),

            ('start', 'common-name'): (BiocycParseTag('common-name'), XmlParseMethod.UNIQUE),
            ('start', 'systematic-name'): (BiocycParseTag('systematic-name'), XmlParseMethod.UNIQUE),
            ('start', 'synonym'): (BiocycParseTag('synonym'), XmlParseMethod.LIST),
            ('start', 'spontaneous'): (BiocycParseTag('spontaneous'), XmlParseMethod.UNIQUE),

            ('start', 'std-reduction-potential'): (BiocycParseTag('std-reduction-potential'), XmlParseMethod.UNIQUE),
            ('start', 'regulated-by'): (BiocycParserPass('regulated-by'), XmlParseMethod.UNIQUE),
            ('start', 'signal'): (BiocycParserPass('signal'), XmlParseMethod.LIST),

            ('start', 'ec-number'): (BiocycParserPass('ec-number'), XmlParseMethod.LIST),
            ('start', 'reaction-list'): (BiocycParserPass('reaction-list'), XmlParseMethod.UNIQUE),
            ('start', 'reaction-direction'): (BiocycParserPass('reaction-direction'), XmlParseMethod.UNIQUE),
            ('start', 'orphan'): (BiocycParserPass('orphan'), XmlParseMethod.UNIQUE),
            ('start', 'physiologically-relevant'): (
                BiocycParserPass('physiologically-relevant'), XmlParseMethod.UNIQUE),
            ('start', 'credits'): (BiocycParserPass('credits'), XmlParseMethod.UNIQUE),
            ('start', 'in-pathway'): (BiocycParserPass('in-pathway'), XmlParseMethod.UNIQUE),
            ('start', 'gibbs-0'): (BiocycParserPass('gibbs-0'), XmlParseMethod.UNIQUE),
            ('start', 'comment'): (BiocycParseTag('comment'), XmlParseMethod.UNIQUE),
            ('start', 'citation'): (BiocycParserPass('citation'), XmlParseMethod.LIST),
        }


class BiocycPathwayParser(BiocycXmlParser):

    def __init__(self):
        super().__init__('Pathway')
        self.tag_capture = {
            ('start', 'parent'): (BiocycParseSubElements('parent',
                                                         ['Pathway'], True), XmlParseMethod.LIST),
            ('start', 'instance'): (BiocycParseSubElements('instance',
                                                           ['Pathway'], True), XmlParseMethod.LIST),
            ('start', 'subclass'): (BiocycParseSubElements('subclass',
                                                           ['Pathway'], True), XmlParseMethod.LIST),
            ('start', 'super-pathways'): (BiocycParseSubElements('super-pathways',
                                                                 ['Pathway'], False), XmlParseMethod.UNIQUE),
            ('start', 'in-pathway'): (BiocycParseSubElements('in-pathway',
                                                             ['Pathway'], False), XmlParseMethod.UNIQUE),
            ('start', 'sub-pathways'): (BiocycParseSubElements('sub-pathways',
                                                               ['Pathway'], False), XmlParseMethod.UNIQUE),
            ('start', 'common-name'): (BiocycParseTag('common-name'), XmlParseMethod.UNIQUE),
            ('start', 'abbrev-name'): (BiocycParseTag('abbrev-name'), XmlParseMethod.UNIQUE),
            ('start', 'comment'): (BiocycParseTag('comment'), XmlParseMethod.LIST),
            ('start', 'synonym'): (BiocycParseTag('synonym'), XmlParseMethod.LIST),

            ('start', 'taxonomic-range'): (BiocycParserPass('taxonomic-range'), XmlParseMethod.UNIQUE),
            ('start', 'pathway-link'): (BiocycParserPass('pathway-link'), XmlParseMethod.LIST),  # want
            ('start', 'credits'): (BiocycParserPass('credits'), XmlParseMethod.UNIQUE),
            ('start', 'citation'): (BiocycParserPass('citation'), XmlParseMethod.LIST),
            ('start', 'evidence'): (BiocycParserPass('evidence'), XmlParseMethod.LIST),
            ('start', 'reaction-ordering'): (BiocycParserPass('reaction-ordering'), XmlParseMethod.LIST),
            ('start', 'species'): (BiocycParserPass('species'), XmlParseMethod.LIST),
            ('start', 'enzymes-not-used'): (BiocycParserPass('enzymes-not-used'), XmlParseMethod.UNIQUE),
            ('start', 'reaction-layout'): (BiocycParserPass('reaction-layout'), XmlParseMethod.LIST),
            ('start', 'dblink'): (BiocycParserDblink(), XmlParseMethod.LIST),
            ('start', 'reaction-list'): (
                BiocycParseSubElements('reaction-list', ['Reaction', 'Pathway'], False), XmlParseMethod.UNIQUE),
            ('start', 'hypothetical-reactions'): (
                BiocycParseSubElements('hypothetical-reactions', ['Reaction'], False), XmlParseMethod.UNIQUE),

        }


class BiocycProteinParser(BiocycXmlParser):

    def __init__(self):
        super().__init__('Protein')
        self.tag_capture = {
            ('start', 'cml'): (BiocycParserPass('cml'), XmlParseMethod.UNIQUE),  # want
            ('start', 'regulates'): (BiocycParserPass('regulates'), XmlParseMethod.UNIQUE),  # want

            ('start', 'parent'): (BiocycParseSubElements('parent', ['Protein', 'Compound'], True), XmlParseMethod.LIST),
            ('start', 'subclass'): (BiocycParserSingle('Protein', 'subclass'), XmlParseMethod.LIST),
            ('start', 'instance'): (BiocycParserSingle('Protein', 'instance'), XmlParseMethod.LIST),
            ('start', 'modified-form'): (BiocycParseSubElements('modified-form', ['Protein'], False), XmlParseMethod.UNIQUE),
            ('start', 'unmodified-form'): (
                BiocycParseSubElements('unmodified-form', ['Protein'], True), XmlParseMethod.UNIQUE),
            ('start', 'location'): (BiocycParserSingle('cco', 'location'), XmlParseMethod.LIST),
            ('start', 'has-go-term'): (BiocycParserPass('has-go-term'), XmlParseMethod.LIST),  # want
            ('start', 'cofactors-of'): (BiocycParserPass('cofactors-of'), XmlParseMethod.UNIQUE),  # want
            ('start', 'has-feature'): (BiocycParserPass('has-feature'), XmlParseMethod.LIST),
            ('start', 'gene'): (BiocycParserSingle('Gene', 'gene'), XmlParseMethod.UNIQUE),
            ('start', 'species'): (BiocycParserPass('species'), XmlParseMethod.UNIQUE),

            ('start', 'coding-segment'): (BiocycParserPass('coding-segment'), XmlParseMethod.LIST),

            ('start', 'component-of'): (BiocycParserSingle('Protein', 'component-of'), XmlParseMethod.UNIQUE),
            ('start', 'component'): (BiocycParserPass('component'), XmlParseMethod.LIST),  # want Protein/coefficient

            ('start', 'common-name'): (BiocycParseTag('common-name'), XmlParseMethod.UNIQUE),
            ('start', 'abbrev-name'): (BiocycParseTag('abbrev-name'), XmlParseMethod.UNIQUE),

            ('start', 'n-plus-1-name'): (BiocycParseTag('n-plus-1-name'), XmlParseMethod.UNIQUE),
            ('start', 'n-name'): (BiocycParseTag('n-name'), XmlParseMethod.UNIQUE),
            ('start', 'n-1-name'): (BiocycParseTag('n-1-name'), XmlParseMethod.UNIQUE),
            ('start', 'inchi-key'): (BiocycParseTag('inchi-key'), XmlParseMethod.UNIQUE),
            ('start', 'inchi'): (BiocycParseTag('inchi'), XmlParseMethod.UNIQUE),
            ('start', 'monoisotopic-mw'): (BiocycParseTag('monoisotopic-mw'), XmlParseMethod.UNIQUE),

            ('start', 'synonym'): (BiocycParseTag('synonym'), XmlParseMethod.LIST),
            ('start', 'molecular-weight'): (BiocycParseTag('molecular-weight'), XmlParseMethod.UNIQUE),
            ('start', 'molecular-weight-seq'): (BiocycParserPass('molecular-weight-seq'), XmlParseMethod.UNIQUE),
            ('start', 'accession-1'): (BiocycParseTag('accession-1'), XmlParseMethod.UNIQUE),

            ('start', 'catalyzes'): (BiocycParserPass('catalyzes'), XmlParseMethod.UNIQUE),
            ('start', 'symmetry'): (BiocycParserPass('symmetry'), XmlParseMethod.LIST),
            ('start', 'dna-footprint-size'): (BiocycParserPass('dna-footprint-size'), XmlParseMethod.UNIQUE),
            ('start', 'consensus-sequence'): (BiocycParserPass('consensus-sequence'), XmlParseMethod.LIST),
            ('start', 'gene'): (BiocycParserPass('gene'), XmlParseMethod.LIST),
            ('start', 'molecular-weight-exp'): (BiocycParserPass('molecular-weight-exp'), XmlParseMethod.LIST),
            ('start', 'credits'): (BiocycParserPass('credits'), XmlParseMethod.UNIQUE),
            ('start', 'evidence'): (BiocycParserPass('evidence'), XmlParseMethod.LIST),
            ('start', 'comment'): (BiocycParseTag('comment'), XmlParseMethod.LIST),
            ('start', 'citation'): (BiocycParserPass('citation'), XmlParseMethod.LIST),

            ('start', 'appears-in-right-side-of'): (BiocycParserPass('appears-in-right-side-of'), XmlParseMethod.UNIQUE),
            ('start', 'appears-in-left-side-of'): (BiocycParserPass('appears-in-left-side-of'), XmlParseMethod.UNIQUE),

            ('start', 'expression-mechanism'): (BiocycParserPass('expression-mechanism'), XmlParseMethod.UNIQUE),
            ('start', 'isozyme-sequence-similarity'): (BiocycParserPass('isozyme-sequence-similarity'),
                                                       XmlParseMethod.LIST),
            ('start', 'pi'): (BiocycParserPass('pi'), XmlParseMethod.LIST),

            ('start', 'dblink'): (BiocycParserDblink(), XmlParseMethod.LIST),
        }

