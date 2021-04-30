from collections import Iterable, OrderedDict
from typing import Iterator, List, Optional, Dict

try:
    from ete3 import PhyloTree
except:
    pass


GAP_CHARS = "X-"


class Record(Iterable):
    """
    Very lightweight wrapper around fasta sequence files.
    """

    def __init__(self, title: str, sequence: str, seq_len: int = 60):
        self.title = title
        self.sequence = sequence
        self.seq_len = seq_len

    def __getitem__(self, item):
        return Record(self.title, self.sequence[item], self.seq_len)

    def as_string(self, pad_to_len: int = -1):
        header = ">" + self.title if self.title else None
        seqs = []
        counter = 0
        tracker = ""

        sequence = self.sequence

        if pad_to_len > len(sequence):
            sequence = sequence + ('-' * (pad_to_len - len(sequence)))

        for f in sequence:
            if counter == self.seq_len:
                counter = 0
                seqs.append(tracker)
                tracker = ""
            counter += 1
            tracker += f

        if len(tracker) > 0:
            seqs.append(tracker)

        if header:
            return header + "\n" + "\n".join(seqs)
        else:
            return "\n".join(seqs)

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return self.as_string()

    def __add__(self, other):
        return Record(self.title, self.sequence + other.sequence, self.seq_len)

    def __len__(self):
        return len(self.sequence)

    def save_to(self, filename, write_mode='w'):
        with open(filename, write_mode) as f:
            f.write(str(self))

    def __iter__(self) -> Iterator[chr]:
        for f in self.sequence:
            yield f

    def strip_gaps(self, include_unknowns: bool = False):
        return Record(self.title, self.sequence.replace(GAP_CHARS if include_unknowns else '-', ''), self.seq_len)

    def gap_percentage(self) -> float:
        total_gaps = 0
        for char in self.sequence:
            if char in GAP_CHARS:
                total_gaps += 1

        return total_gaps / len(self.sequence)

    def pairwise_comparison(self, other: 'Record', leniency=0) -> float:
        max_match = len(self.sequence)
        max_all = max(max_match, len(other.sequence))
        min_all = min(max_match, len(other.sequence))
        len_diff = abs(len(self.sequence) - len(other.sequence))
        first_pass_matches = 0
        for first, second in zip(self.sequence, other.sequence):
            if first == second:
                first_pass_matches += 1

        if first_pass_matches == max_match and len_diff == 0:
            return 1.0

        curr_misses = (max_match - first_pass_matches) + len_diff
        curr_score = 1 - (curr_misses / max_all)

        window_tracker = 1
        matches = 0.0
        while window_tracker < leniency+1:
            i = 0
            for first, second in zip(self.sequence, other.sequence):
                if first == second:
                    matches += 1.0
                else:
                    curr_window = 1
                    while curr_window <= window_tracker:
                        lookahead = i + curr_window
                        lookbehind = i - curr_window
                        if lookahead < min_all:
                            if first == other.sequence[lookahead]:
                                matches += 1.0
                                break
                        if lookbehind >= 0:
                            if first == other.sequence[lookbehind]:
                                matches += 1.0
                                break

                        curr_window += 1

                i += 1

            round_misses = (float(max_match) - matches) + float(len_diff)
            round_score = 1 - (round_misses / float(max_all))
            if round_score > curr_score:
                curr_score = round_score

            window_tracker += 1

        return curr_score

    def __copy__(self):
        return Record(self.title, self.sequence, self.seq_len)

    def copy(self) -> 'Record':
        return self.__copy__()


def phylogenetic_sort(records: List[Record], tree: 'PhyloTree', important_taxon: str = None) -> List[Record]:
    to_return = order_records(records, list(tree.iter_leaf_names()))

    if important_taxon:
        to_return = [[r for r in records if r.title.strip() == important_taxon.strip()][0]] + to_return

    return to_return


def order_records(records: List[Record], ordering: List[str]) -> List[Record]:
    title2index = {t.strip(): i for i, t in enumerate(ordering)}

    return list(sorted(records, key=lambda r: (title2index.get(r.title.strip(), len(records)+1), r.title)))


def parse_record(fasta_str: str) -> Record:
    lines = fasta_str.splitlines()
    header = None
    was_header_checked = False
    seq = ""
    seq_len = None

    for line in lines:
        if not was_header_checked:
            header = line[1:].lstrip() if line.startswith('>') else None
            if header:
                was_header_checked = True
                continue
        seq += line

        if seq_len is None:
            seq_len = len(line)

    return Record(header, seq, seq_len)


def parse_records(fastas: str) -> List[Record]:
    return list(generate_records(fastas))  # For backwards compat


def generate_records(fastas: str):  # Designed to handle large files
    fastas = fastas.strip()

    if len(fastas) > 0 and fastas[0] != '>':
        yield parse_record(fastas)

    else:
        curr_seq = None
        for line in fastas.splitlines():
            if len(line.strip()) == 0:
                continue

            if line.startswith(">") and curr_seq is not None:
                yield parse_record("\n".join(curr_seq))
                curr_seq = None

            if curr_seq is None:
                curr_seq = []

            curr_seq.append(line)

        if curr_seq is not None:
            yield parse_record("\n".join(curr_seq))


def _read_records_generator(filename: str, _generate_records: bool = False, generate_chunk_size: int = 100):
    counter = 0
    lines = ""
    with open(filename, 'r') as f:
        for line in f:
            if '>' in line:
                if counter > generate_chunk_size:
                    yield from generate_records(lines)
                    lines = ""
                    counter = 0
                counter += 1
            lines += line
    if len(lines) > 0:
        yield from generate_records(lines)


def read_records(filename: str, _generate_records: bool = False, generate_chunk_size: int = 100) -> List[Record]:
    if not _generate_records:
        with open(filename, 'r') as f:
            x = parse_records(f.read())
            return x
    else:
        return _read_records_generator(filename, _generate_records, generate_chunk_size)


def write_records(filename: str, records: List[Record], pad: bool = False):
    if isinstance(records, Record):
        records = [records]

    max_len = -1 if not pad else max([len(r.sequence) for r in records])

    with open(filename, 'w') as f:
        f.writelines([record.as_string(max_len) + "\n" for record in records])


def records2dict(records: List[Record]) -> Dict[str, List[Record]]:
    record_dict = dict()
    for r in records:
        if r not in record_dict:
            record_dict[r.title] = [r]
        else:
            record_dict[r.title].append(r)
    return record_dict


def remove_paralog_groups(filename: str, outfile: str):
    recs = records2dict(read_records(filename))

    pruned = [rs[0] for rs in recs.values() if len(rs) == 1]

    write_records(outfile, pruned)


def pairwise_comparison(fasta1: str, fasta2: str, leniency: int = 0) -> float:
    record1 = parse_record(fasta1)
    record2 = parse_record(fasta2)
    return record1.pairwise_comparison(record2, leniency)


def prune_nonsimilar_paralogs(fasta: str) -> Optional[str]:
    records = parse_records(fasta)
    grouped = {}
    for record in records:
        if record.title not in grouped:
            grouped[record.title] = []
        grouped[record.title].append(record)

    unambiguous_records = []
    ambiguous_records = []
    for group in grouped.values():
        if len(group) == 1:
            unambiguous_records.append(group[0])
        elif len(group) > 1:
            ambiguous_records.append(group)

    if len(unambiguous_records) == 0:
        return None

    if len(ambiguous_records) == 0:
        return "\n".join([str(x) for x in unambiguous_records])

    disambiguated_records = []
    for records in ambiguous_records:
        max_similarity = -1
        best_match = None

        for record in records:
            similarities = []
            for unambiguous_record in unambiguous_records:
                similarities.append(record.pairwise_comparison(unambiguous_record, leniency=1))
            avg_similarity = sum(similarities) / len(similarities)

            if not best_match or avg_similarity > max_similarity:
                max_similarity = avg_similarity
                best_match = record

        disambiguated_records.append(best_match)

    return "\n".join([str(x) for x in unambiguous_records + disambiguated_records])


def gap_percentage(fasta: str) -> float:
    return parse_record(fasta).gap_percentage()
