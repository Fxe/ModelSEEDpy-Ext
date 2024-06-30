import requests
import xml.etree.ElementTree as ET


def parse_link_set_db(tree):
    res = set()
    for el in tree:
        if el.tag == "Link":
            for o in el:
                if o.tag == "Id":
                    res.add(int(o.text))
    return res


def parse_link_set(tree):
    link_set = {}
    for el in tree:
        if el.tag == "LinkSetDb":
            link_set["link_set_db"] = parse_link_set_db(el)
    return link_set


def get_link_set(res):
    link_set = None
    tree = ET.fromstring(res.content)
    for el in tree:
        if el.tag == "LinkSet":
            link_set = parse_link_set(el)
    return link_set


def parse_document_summary_set(tree):
    documents = []
    for el in tree:
        if el.tag == "DocumentSummary":
            documents.append(parse_document_summary(el))
        else:
            print("DocumentSummary?", el, el.tag, el.text, el.attrib)
    return documents


def parse_document_summary(tree):
    parse_function = {
        "RsUid": int,
        "GbUid": int,
        "AssemblyAccession": str,
        "LastMajorReleaseAccession": str,
        "LatestAccession": str,
        "ChainId": int,
        "AssemblyName": str,
        "UCSCName": str,
        "EnsemblName": str,
        "Taxid": int,
        "Organism": str,
        "SpeciesTaxid": int,
        "SpeciesName": str,
        "AssemblyType": str,
        "AssemblyStatus": str,
        "AssemblyStatusSort": str,
        "WGS": str,
        "Coverage": float,
        "ReleaseLevel": str,
        "ReleaseType": str,
        "SubmitterOrganization": str,
        "RefSeq_category": str,
        "FromType": str,
        "ContigN50": str,
        "ScaffoldN50": str,
        "FtpPath_GenBank": str,
        "FtpPath_RefSeq": str,
        "FtpPath_Assembly_rpt": str,
        "FtpPath_Stats_rpt": str,
        "FtpPath_Regions_rpt": str,
        "BioSampleAccn": str,
        "BioSampleId": str,
        "AsmReleaseDate_GenBank": str,
        "AsmReleaseDate_RefSeq": str,
        "SeqReleaseDate": str,
        "AsmUpdateDate": str,
        "SubmissionDate": str,
        "LastUpdateDate": str,
        "Primary": str,
    }
    doc = {}
    for el in tree:
        if el.tag in parse_function:
            try:
                doc[el.tag] = parse_function[el.tag](el.text)
            except Exception as ex:
                print("error!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!", el.tag, el.text, ex)
        elif el.tag == "PropertyList":
            property_list = []
            for e in el:
                property_list.append(e.text)
            doc["PropertyList"] = property_list
        elif el.tag == "Synonym":
            synonym = {}
            for e in el:
                synonym[e.tag] = e.text
            doc["Synonym"] = synonym
        elif el.tag == "Meta":
            pass
        else:
            # print(el, el.tag, el.text, el.attrib, el)
            pass
    return doc


def parse_doc_sum(t):
    doc_sum = {}
    for el in t:
        if el.tag == "Id":
            doc_sum["id"] = el.text
        elif el.tag == "Item":
            if "Name" in el.attrib and "Type" in el.attrib:
                if el.attrib["Type"] == "String":
                    doc_sum[el.attrib["Name"]] = str(el.text)
                elif el.attrib["Type"] == "Integer":
                    doc_sum[el.attrib["Name"]] = int(el.text)
                else:
                    print("attrb?", el, el.attrib, el.text)
            else:
                print("element?", el, el.attrib, el.text)
        else:
            print("tag?", el, el.attrib, el.text)
    return doc_sum


class NcbiEutils:

    def __init__(self, api_key=None):
        self.base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        self.api_key = api_key
        pass

    @staticmethod
    def build_url_params(p):
        return "&".join(map(lambda x: str(x[0]) + "=" + str(x[1]), p.items()))

    def esearch(self, db, term, retmode="xml", retmax=10, retstart=0):
        params = {
            "retmode": retmode,
            "db": db,
            "term": term,
            "retmax": retmax,
            "retstart": retstart,
        }
        if self.api_key:
            params["api_key"] = self.api_key

        res = requests.get(
            self.base_url + "/esearch.fcgi?" + self.build_url_params(params)
        )

        tree = ET.fromstring(res.content)

        id_list = set()
        ret_max = None
        ret_start = None
        count = None
        query_translation = None
        for o in tree:
            if o.tag == "IdList":
                id_list = {i.text for i in o}
            elif o.tag == "Count":
                count = int(o.text)
            elif o.tag == "RetStart":
                ret_start = int(o.text)
            elif o.tag == "RetMax":
                ret_max = int(o.text)
            elif o.tag == "QueryTranslation":
                query_translation = o.text
            else:
                print(o.tag, o.text)

        return {
            "ids": id_list,
            "count": count,
            "ret_max": ret_max,
            "ret_start": ret_start,
            "query_translation": query_translation,
        }

    def esummary(self, db, id, retmode="xml"):
        params = {"retmode": retmode, "db": db, "id": id}
        if self.api_key:
            params["api_key"] = self.api_key

        res = requests.get(
            self.base_url + "/esummary.fcgi?" + self.build_url_params(params)
        )
        tree = ET.fromstring(res.content)

        docs = []
        for el in tree:
            if el.tag == "DocSum":
                docs.append(parse_doc_sum(el))
            elif el.tag == "DocumentSummarySet":
                docs.append(parse_document_summary_set(el))
            else:
                print(el.tag)
        return docs

    def elink(self, dbfrom, db, ids):
        params = {"dbfrom": dbfrom, "db": db, "id": ids}
        if self.api_key:
            params["api_key"] = self.api_key
        res = requests.get(
            self.base_url + "/elink.fcgi?" + self.build_url_params(params)
        )
        return res
