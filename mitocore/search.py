from abc import ABCMeta, abstractmethod
from typing import List, Tuple, Dict, Union, Optional, NamedTuple, Any

from mitocore.molecular import Strand, Locus, GeneticSequence


class Site:
    on: GeneticSequence
    at: Locus

    def __init__(self, on: GeneticSequence, at: Locus):
        self.on = on
        self.at = at

    def __len__(self):
        return len(self.at)

    @property
    def sequence(self) -> GeneticSequence:
        if self.at.length is not None:
            return self.on[self.at]
        else:
            return None

    def __repr__(self) -> str:
        return '<Site {}>'.format(self.sequence.forward)


class SearchMatchHandler(metaclass=ABCMeta):

    @abstractmethod
    def __call__(self, sequence: GeneticSequence, locus: Locus) -> Any:
        ...


class SearchTree(NamedTuple):
    results: List[Tuple[Strand, Union[SearchMatchHandler, None]]]
    branches: Dict[str, 'SearchTree']


class SiteSearch:
    _search_tree: SearchTree

    def __init__(self):
        self._search_tree = SearchTree([], {})

    def add_query(self, query: GeneticSequence, handler: Optional[SearchMatchHandler] = None):
        for instance in query.instances:
            for strand in Strand:
                base = self._search_tree
                for n in instance.forward if strand == Strand.forward else instance.reverse:
                    if n not in base.branches.keys():
                        base.branches[n] = SearchTree([], {})
                    base = base.branches[n]
                base.results.append((strand, handler))

    def search_in(self, sequence: GeneticSequence) -> List[Any]:
        results = []
        for i in range(len(sequence)):
            base = self._search_tree
            depth = 0
            while sequence.forward[(i + depth) % len(sequence)] in base.branches.keys():
                base = base.branches[sequence.forward[(i + depth) % len(sequence)]]
                depth += 1
                for strand, handler in base.results:
                    locus = Locus(i if strand == Strand.forward else (i + depth), depth, strand)
                    results.append(Site(sequence, locus) if handler is None else handler(sequence, locus))
        return results
