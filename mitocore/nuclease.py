import json
from abc import ABCMeta, abstractmethod
from typing import Dict, Sequence, NoReturn, Any

from mitocore.molecular import GeneticSequence, Locus, Strand
from mitocore.search import Site, SearchMatchHandler


class CasSite(Site, json.JSONEncoder):
    aspects: Dict[str, any]

    def __init__(self, of: 'Cas', on: GeneticSequence, at: Locus):
        super().__init__(on, at)
        self.of = of
        self.aspects = {
            'spacer': self.spacer
        }

    @property
    def spacer(self) -> GeneticSequence:
        return self.on[self.at + self.of.spacer]

    def repr(self) -> str:
        return '<{cas} site {pam} {spacer}>'.format_map({
            'cas': self.of.name,
            'pam': self.on[self.at].forward,
            'spacer': self.spacer.forward
        })


class CasSiteEncoder(json.JSONEncoder):

    def default(self, obj: CasSite) -> object:
        dct = {
            '__cas_site__': True,
            'of': obj.of.name,
            'at': obj.at.position,
            'strand': str(obj.at.strand),
        }
        return dct


class CasModel(metaclass=ABCMeta):
    ...


class CasLocalModel(CasModel, metaclass=ABCMeta):
    """ Model, describing aspects of single CasSite, independent of another sites
    """

    @abstractmethod
    def evaluate(self, site: CasSite) -> NoReturn:
        ...


class CasGlobalModel(CasModel, metaclass=ABCMeta):
    """ Model, describing sequence-wise aspects of single CasSite"""

    @abstractmethod
    def evaluate(self, site: CasSite, counterparts: Sequence[CasSite]) -> NoReturn:
        ...


class Cas(SearchMatchHandler):
    name: str
    pam: GeneticSequence
    spacer: Locus
    cuts: Sequence[Locus]
    models: Sequence[CasLocalModel]

    def __init__(self,
                 name: str,
                 pam: GeneticSequence,
                 spacer: Locus,
                 cuts: Sequence[Locus] = (),
                 models: Sequence[CasLocalModel] = ()):
        self.name = name
        self.pam = pam
        self.spacer = spacer
        self.cuts = cuts
        self.models = models

    def __call__(self, sequence: GeneticSequence, locus: Locus) -> CasSite:
        site = CasSite(self, sequence, locus)
        for m in self.models:
            m.evaluate(site)
        return site

    @classmethod
    def json_hook(cls, dct: Dict) -> Any:
        if '__cas__' in dct.keys():
            return Cas(
                name=dct['name'],
                pam=GeneticSequence(dct['pam']),
                spacer=dct['spacer'],
            )
        elif '__locus__' in dct.keys():
            return Locus(
                position=dct['position'],
                length=dct['length'],
                strand=Strand(dct['strand'])
            )
        else:
            return dct

    @classmethod
    def from_json(cls, json_: str) -> Sequence['Cas']:
        return json.loads(json_, object_hook=Cas.json_hook)
