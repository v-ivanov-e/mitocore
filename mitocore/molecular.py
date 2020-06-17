from abc import ABCMeta, abstractmethod
from enum import Enum
from itertools import product
from typing import NamedTuple, Sized, Dict, Sequence, Optional, Tuple, Union


class Strand(Enum):
    forward = '+'
    reverse = '-'


class Locus(NamedTuple):
    position: int
    length: int
    strand: Strand

    def __add__(self, other: 'Locus') -> 'Locus':
        strand = self.strand if other.strand == Strand.forward else \
            (Strand.forward if self.strand == Strand.reverse else Strand.reverse)
        position = self.position + other.position
        length = other.length
        return Locus(position, length, strand)

    @classmethod
    def json_hook(cls, dct: Dict) -> Union[Dict, 'Locus']:
        if '__locus__' in dct.keys():
            return Locus(
                position=dct['position'],
                length=dct['length'],
                strand=Strand(dct['strand'])
            )
        else:
            return dct


class AbstractGeneticSequence(Sized, metaclass=ABCMeta):

    @abstractmethod
    def __len__(self) -> int:
        ...

    @abstractmethod
    def __getitem__(self, item: Locus) -> 'AbstractGeneticSequence':
        ...

    _complements_table: Dict[str, str] = {
        'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
        'S': 'W', 'W': 'S',
        'Y': 'R', 'R': 'Y',
        'K': 'M', 'M': 'K',
        'V': 'B', 'B': 'V',
        'D': 'H', 'H': 'D',
        'N': 'N'
    }

    _inclusions_table: Dict[str, Sequence[str]] = {
        'S': ['C', 'G', 'S'], 'W': ['A', 'T', 'W'],
        'Y': ['C', 'T', 'Y'], 'R': ['A', 'G', 'R'],
        'K': ['G', 'T', 'K'], 'M': ['A', 'C', 'K'],
        'B': ['C', 'G', 'T', 'S', 'Y', 'K', 'B'],
        'D': ['A', 'G', 'T', 'R', 'D'],
        'H': ['A', 'C', 'T', 'M', 'H'],
        'V': ['A', 'C', 'G', 'S', 'R', 'M', 'V'],
        'N': ['A', 'C', 'G', 'T', 'S', 'W', 'Y', 'R', 'K', 'M', 'N'],
    }

    @classmethod
    def generate_complement(cls, sequence: str, _skip: str = '-') -> str:
        table = cls._complements_table
        return ''.join(table[n] if n in table.keys() else _skip for n in sequence)

    @classmethod
    def generate_instances(cls, sequence: str) -> Sequence[str]:
        table = cls._inclusions_table
        return [''.join(instance) for instance in product(*[table[n] if n in table.keys() else tuple(n)
                                                            for n in sequence])]


class GeneticSequence(AbstractGeneticSequence):
    forward: str
    reverse: str

    def __init__(self, strand: str, reverse_complement_strand: Optional[str] = None):
        self.forward = strand
        self.reverse = reverse_complement_strand

    @property
    def reverse(self):
        if self._reverse is None:
            self._reverse = self.generate_complement(self.forward[::-1])
        return self._reverse

    @reverse.setter
    def reverse(self, value: str):
        self._reverse = value

    @property
    def instances(self) -> Sequence['GeneticSequence']:
        return [GeneticSequence(i) for i in self.generate_instances(self.forward)]

    def __len__(self) -> int:
        return len(self.forward)

    def __neg__(self):
        return GeneticSequence(self.reverse, self.forward)

    def __add__(self, other: 'GeneticSequence') -> 'GeneticSequence':
        return GeneticSequence(self.forward + other.forward)

    def __getitem__(self, item: Locus) -> 'GeneticSequence':
        if item.strand == Strand.forward:
            start = item.position % len(self)
            stop = start + item.length
            return GeneticSequence(
                self.forward[start:min(stop, len(self))] +
                self.forward[0:max(0, stop - len(self))]
            )
        else:
            start = (len(self) - item.position) % len(self)
            stop = start + item.length
            return GeneticSequence(
                self.reverse[start:min(stop, len(self))] +
                self.reverse[0:max(0, stop - len(self))]
            )

    @classmethod
    def from_fasta(cls,
                   fasta: str,
                   _sequence_comment: Optional[str] = '>',
                   _fasta_comment: Optional[str] = '!') -> Sequence[Tuple[str, 'GeneticSequence']]:
        com, seq = None, None
        result = []
        for line in fasta.split('\n'):
            if line.startswith(_sequence_comment):
                if com is not None:
                    result.append((com, GeneticSequence(seq)))
                com, seq = line[len(_sequence_comment):], ''
            elif line.startswith(_fasta_comment):
                continue
            else:
                seq += line
        if com is not None:
            result.append((com, GeneticSequence(seq)))
        return result
