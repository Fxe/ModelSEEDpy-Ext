from cobrakbase.kbaseapi import KBaseAPI


class KBaseRE(KBaseAPI):

    def __init__(self, re, token=None, dev=False, config=None):
        super().__init__(token, dev, config)

    def get_from_ws(self, id_or_ref, workspace=None):
        super().get_from_ws(id_or_ref, workspace)
