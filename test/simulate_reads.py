#!/usr/bin/env python3
"""Simulate the hybseq test reads from the committed test references.

test/data/ holds real data: Drosophila melanogaster chromosome 4 (ref_a,
GCF_000001215.4, 1,348,131 bp), Drosophila simulans w501 chromosome 4 (ref_b,
Prin_Dsim_3.1, 1,146,867 bp) and the FlyBase r6.54 chromosome-4 annotation
(ref_a.gtf). The reads are simulated from those references instead of being
committed, so they cannot drift away from them again: placeholder reads that no
longer match the references are the failure this script removes.

    python3 test/simulate_reads.py                 # (re)generate, 5000 reads/sample
    python3 test/simulate_reads.py --reads 1000    # quicker, lower coverage
    python3 test/simulate_reads.py --check         # verify the existing fastq

Only the Python standard library is used (no pysam, no containers, no network)
and the seed is fixed, so the same references and the same --reads always give
the same reads.
"""

import argparse
import bisect
import gzip
import math
import os
import random
import re
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data")
FASTQ_DIR = os.path.join(DATA_DIR, "fastq")

DEFAULT_SEED = 20261015
DEFAULT_READS = 5000          # read pairs per paired-end sample, and reads per
                              # single-end sample (see --reads)
DEFAULT_READ_LENGTH = 150
DEFAULT_INSERT_MEAN = 350
DEFAULT_INSERT_SD = 50
DEFAULT_RNA_INSERT_MEAN = 300
DEFAULT_RNA_INSERT_SD = 50
DEFAULT_ERROR_RATE = 0.005
DEFAULT_CHECK_READS = 50
DEFAULT_PAIR_CHECKS = 10
# A mixed sample must show every source in the checked reads; a real 50:50 mix
# clears this by a wide margin, while reads from a single source score 0.
MIN_SOURCE_SHARE = 0.25

# Exact-match window used by the --check mapper: long enough to be unique in a
# 1.3 Mbp chromosome, short enough that a read still contains a clean one at the
# default error rate.
WINDOW = 25

# How many alternative placements a repeat can offer before --check gives up
# trying to fit a mate pair into one of them.
MAX_ANCHORS = 8

# Sample plan, kept in step with test/samplesheet.csv. "sources" lists
# (kind, reference, weight) with kind "genome" (uniform over the contigs) or
# "transcript" (weighted by spliced transcript length, so RNA reads can span
# exon-exon junctions). Only ref_a has an annotation, so both RNA samples are
# ref_a-derived and differ in layout alone.
SAMPLES = [
    {
        "sample": "dna_a",
        "assay": "dnaseq",
        "layout": "pe",
        "files": ["dna_a_R1.fq.gz", "dna_a_R2.fq.gz"],
        "sources": [("genome", "ref_a", 1.0)],
    },
    {
        "sample": "dna_b",
        "assay": "dnaseq",
        "layout": "pe",
        "files": ["dna_b_R1.fq.gz", "dna_b_R2.fq.gz"],
        "sources": [("genome", "ref_b", 1.0)],
    },
    {
        "sample": "dna_hybrid",
        "assay": "dnaseq",
        "layout": "pe",
        "files": ["dna_hybrid_R1.fq.gz", "dna_hybrid_R2.fq.gz"],
        "sources": [("genome", "ref_a", 0.5), ("genome", "ref_b", 0.5)],
    },
    {
        "sample": "dna_c",
        "assay": "dnaseq",
        "layout": "se",
        "files": ["dna_c.fq.gz"],
        "sources": [("genome", "ref_a", 0.5), ("genome", "ref_b", 0.5)],
    },
    {
        "sample": "rna_a",
        "assay": "rnaseq",
        "layout": "se",
        "files": ["rna_a.fastq.gz"],
        "sources": [("transcript", "ref_a", 1.0)],
    },
    {
        "sample": "rna_b",
        "assay": "rnaseq",
        "layout": "pe",
        "files": ["rna_b_R1.fq.gz", "rna_b_R2.fq.gz"],
        "sources": [("transcript", "ref_a", 1.0)],
    },
]

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_TRANSITIONS = {"A": "G", "G": "A", "C": "T", "T": "C"}
_TRANSVERSIONS = {"A": ("C", "T"), "G": ("C", "T"), "C": ("A", "G"), "T": ("A", "G")}
_SEPARATOR = "N" * 10


def die(message):
    sys.stderr.write("simulate_reads.py: error: %s\n" % message)
    sys.exit(2)


def revcomp(seq):
    return seq.translate(_COMPLEMENT)[::-1]


# ---------------------------------------------------------------- references


def load_fasta(path):
    """Return [(contig name, upper-case sequence)]."""
    contigs = []
    name = None
    chunks = []
    try:
        handle = open(path)
    except OSError as exc:
        die("cannot read %s (%s)" % (path, exc))
    with handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    contigs.append((name, "".join(chunks).upper()))
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        contigs.append((name, "".join(chunks).upper()))
    if not contigs:
        die("no sequences in %s" % path)
    return contigs


def n_runs(seq):
    """[(start, end)] of the N runs in `seq`, 0-based, end exclusive."""
    runs = []
    start = seq.find("N")
    while start != -1:
        end = start
        while end < len(seq) and seq[end] == "N":
            end += 1
        runs.append((start, end))
        start = seq.find("N", end)
    return runs


def spans_n(runs, start, end):
    return any(run_start < end and run_end > start for run_start, run_end in runs)


def load_genomes(ref_a, ref_b):
    """{"ref_a": {...}, "ref_b": {...}} with contigs and N runs per contig."""
    genomes = {}
    for name, path in (("ref_a", ref_a), ("ref_b", ref_b)):
        contigs = load_fasta(path)
        cumulative = []
        total = 0
        for _, seq in contigs:
            total += len(seq)
            cumulative.append(total)
        genomes[name] = {
            "contigs": contigs,
            "cumulative": cumulative,
            "total": total,
            "nruns": {contig: n_runs(seq) for contig, seq in contigs},
        }
    return genomes


def load_transcriptome(path, contigs, minimum_length):
    """Spliced transcript sequences for the annotated reference.

    For a '-' transcript the exons are concatenated in genomic order and the
    result is reverse-complemented, so every sequence is in transcript
    orientation and a fragment can span exon-exon junctions. Transcripts shorter
    than `minimum_length` cannot host a fragment and are left out of the
    simulation (they stay in `flat`, which --check searches).
    """
    sequences = dict(contigs)
    exons = defaultdict(list)
    strand = {}
    contig_of = {}
    try:
        handle = open(path)
    except OSError as exc:
        die("cannot read %s (%s)" % (path, exc))
    with handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            match = re.search(r'transcript_id "([^"]+)"', fields[8])
            if not match:
                continue
            if fields[0] not in sequences:
                die("GTF contig %s is not in the reference FASTA" % fields[0])
            transcript = match.group(1)
            exons[transcript].append((int(fields[3]), int(fields[4])))
            strand[transcript] = fields[6]
            contig_of[transcript] = fields[0]
    transcripts = []
    for transcript in sorted(exons):
        seq = "".join(sequences[contig_of[transcript]][start - 1:end]
                      for start, end in sorted(exons[transcript]))
        if strand[transcript] == "-":
            seq = revcomp(seq)
        transcripts.append((transcript, seq))
    if not transcripts:
        die("no exon features with a transcript_id in %s" % path)
    eligible = [(name, seq) for name, seq in transcripts
                if len(seq) >= minimum_length]
    if not eligible:
        die("no transcript in %s is at least %d bp long"
            % (path, minimum_length))
    cumulative = []
    total = 0
    for _, seq in eligible:
        total += len(seq)
        cumulative.append(total)
    return {
        "eligible": eligible,
        "cumulative": cumulative,
        "total": total,
        "num_transcripts": len(transcripts),
        # One searchable copy with N separators: no match can run across a
        # transcript boundary, and a flat position maps back through `spans`.
        "flat": _SEPARATOR.join(seq for _, seq in transcripts),
        "spans": transcript_spans(transcripts),
    }


def transcript_spans(transcripts):
    spans = []
    cursor = 0
    for name, seq in transcripts:
        spans.append((cursor, cursor + len(seq), name))
        cursor += len(seq) + len(_SEPARATOR)
    return spans


# --------------------------------------------------------------- simulation


def poisson(rng, lam):
    """Knuth's sampler: the error count per read is Binomial(L, rate), which is
    Poisson(L * rate) at the rates used here."""
    limit = math.exp(-lam)
    count = 0
    product = 1.0
    while True:
        count += 1
        product *= rng.random()
        if product <= limit:
            return count - 1


def substitute(rng, base):
    if rng.random() < 0.7:            # transition bias, as in real data
        return _TRANSITIONS[base]
    return rng.choice(_TRANSVERSIONS[base])


def mutate(rng, seq, rate):
    """Apply substitution errors; return (sequence, quality) as bytes.

    Correct bases are Q40 ('I'), error bases Q20 ('5'), Phred+33.
    """
    if rate <= 0:
        return seq.encode(), b"I" * len(seq)
    count = poisson(rng, rate * len(seq))
    if count <= 0:
        return seq.encode(), b"I" * len(seq)
    bases = list(seq)
    quality = bytearray(b"I" * len(seq))
    for position in rng.sample(range(len(seq)), min(count, len(seq))):
        bases[position] = substitute(rng, bases[position])
        quality[position] = ord("5")
    return "".join(bases).encode(), bytes(quality)


def draw_insert(rng, mean, sd, minimum):
    """Fragment length from a normal distribution truncated at 3 SD.

    Libraries are size-selected, so the simulator does not emit the rare
    multi-sigma outliers a plain normal would produce (and which --check would
    rightly call implausible).
    """
    low = max(minimum, int(round(mean - 3 * sd)))
    high = max(low, int(round(mean + 3 * sd)))
    for _ in range(100):
        length = int(rng.gauss(mean, sd))
        if low <= length <= high:
            return length
    return max(low, min(int(mean), high))


def weighted_pick(rng, cumulative, total):
    """Index drawn proportionally to the lengths behind `cumulative`."""
    target = rng.random() * total
    return min(bisect.bisect_right(cumulative, target), len(cumulative) - 1)


def pick_source(rng, sources):
    target = rng.random()
    accumulated = 0.0
    for kind, reference, weight in sources:
        accumulated += weight
        if target < accumulated:
            return kind, reference
    return sources[-1][0], sources[-1][1]


def simulate_dna_pair(rng, genome, args, paired):
    """(r1, r2) from a random place in `genome`; r2 is None for single-end."""
    while True:
        index = weighted_pick(rng, genome["cumulative"], genome["total"])
        contig, seq = genome["contigs"][index]
        insert = draw_insert(rng, args.insert_mean, args.insert_sd,
                             args.read_length + 1)
        if insert > len(seq):
            continue
        start = rng.randrange(0, len(seq) - insert + 1)
        if spans_n(genome["nruns"][contig], start, start + insert):
            continue
        fragment = seq[start:start + insert]
        if rng.random() < 0.5:
            fragment = revcomp(fragment)
        r1 = fragment[:args.read_length]
        r2 = revcomp(fragment[-args.read_length:]) if paired else None
        return r1, r2


def simulate_rna_pair(rng, transcriptome, args, paired):
    """(r1, r2) from a spliced transcript, in transcript orientation."""
    while True:
        index = weighted_pick(rng, transcriptome["cumulative"],
                              transcriptome["total"])
        seq = transcriptome["eligible"][index][1]
        insert = draw_insert(rng, args.rna_insert_mean, args.rna_insert_sd,
                             args.read_length + 1)
        if insert > len(seq):
            continue
        start = rng.randrange(0, len(seq) - insert + 1)
        fragment = seq[start:start + insert]
        r1 = fragment[:args.read_length]
        r2 = revcomp(fragment[-args.read_length:]) if paired else None
        return r1, r2


def simulate_pair(rng, sources, genomes, transcriptome, args, paired):
    kind, reference = pick_source(rng, sources)
    if kind == "genome":
        return simulate_dna_pair(rng, genomes[reference], args, paired)
    return simulate_rna_pair(rng, transcriptome, args, paired)


# ------------------------------------------------------------------ writing


def write_record(handle, name, seq, quality):
    handle.write(b"@" + name.encode() + b"\n" + seq + b"\n+\n" + quality + b"\n")


def generate(args):
    genomes = load_genomes(args.ref_a, args.ref_b)
    transcriptome = load_transcriptome(args.gtf, genomes["ref_a"]["contigs"],
                                       args.read_length + 1)
    os.makedirs(args.outdir, exist_ok=True)
    print("generating %d read pairs (%d reads) per paired-end sample and "
          "%d reads per single-end sample into %s"
          % (args.reads, 2 * args.reads, args.reads, args.outdir))
    written = 0
    for entry in SAMPLES:
        paired = entry["layout"] == "pe"
        rng = random.Random("%s:%s" % (args.seed, entry["sample"]))
        paths = [os.path.join(args.outdir, name) for name in entry["files"]]
        handles = [gzip.GzipFile(filename=path, mode="wb", compresslevel=6,
                                 mtime=0) for path in paths]
        try:
            for index in range(args.reads):
                r1, r2 = simulate_pair(rng, entry["sources"], genomes,
                                       transcriptome, args, paired)
                name = "%s_%d" % (entry["sample"], index + 1)
                seq, quality = mutate(rng, r1, args.error_rate)
                write_record(handles[0], name + "/1" if paired else name,
                             seq, quality)
                written += 1
                if paired:
                    seq, quality = mutate(rng, r2, args.error_rate)
                    write_record(handles[1], name + "/2", seq, quality)
                    written += 1
        finally:
            for handle in handles:
                handle.close()
        print("  %-44s %s" % (", ".join(entry["files"]),
                              "paired-end" if paired else "single-end"))
    print("wrote %d reads" % written)


# -------------------------------------------------------------------- check


def offsets(length):
    return range(0, max(length - WINDOW, 0) + 1, WINDOW)


def mismatches(query, seq, start):
    return sum(1 for a, b in zip(query, seq[start:start + len(query)]) if a != b)


def find_approx(query, seq, accept):
    """Approximate match of `query` in `seq` (forward strand only).

    Candidate positions come from exact WINDOW-mers of the query, so reads with
    a few simulated errors are found without building a k-mer index (which would
    cost hundreds of MB on a large reference). Returns (mismatches, start) for
    the first hit of at most `accept` mismatches, otherwise the best hit found,
    or None.
    """
    best = None
    length = len(query)
    for offset in offsets(length):
        window = query[offset:offset + WINDOW]
        position = seq.find(window)
        while position != -1:
            start = position - offset
            if start >= 0 and start + length <= len(seq):
                distance = mismatches(query, seq, start)
                if distance <= accept:
                    return (distance, start)
                if best is None or distance < best[0]:
                    best = (distance, start)
            position = seq.find(window, position + 1)
    return best


def footprint(read, seq, accept):
    """Best match of `read` in `seq` on either strand.

    Returns (mismatches, start, strand), where `start` is the 0-based position of
    the read's footprint on the forward strand, or None.
    """
    best = None
    for strand, query in (("+", read), ("-", revcomp(read))):
        hit = find_approx(query, seq, accept)
        if hit is None:
            continue
        if hit[0] <= accept:
            return (hit[0], hit[1], strand)
        if best is None or hit[0] < best[0]:
            best = (hit[0], hit[1], strand)
    return best


def source_footprint(read, kind, reference, genomes, transcriptome, accept):
    """(mismatches, sequence name, start, strand) of the best hit in one source."""
    if kind == "transcript":
        hit = footprint(read, transcriptome["flat"], accept)
        return None if hit is None else (hit[0], "transcriptome", hit[1], hit[2])
    best = None
    for contig, seq in genomes[reference]["contigs"]:
        hit = footprint(read, seq, accept)
        if hit is None:
            continue
        if hit[0] <= accept:
            return (hit[0], contig, hit[1], hit[2])
        if best is None or hit[0] < best[0]:
            best = (hit[0], contig, hit[1], hit[2])
    return best


def sequence_of(kind, reference, name, genomes, transcriptome):
    if kind == "transcript":
        return transcriptome["flat"]
    for contig, seq in genomes[reference]["contigs"]:
        if contig == name:
            return seq
    return None


def best_source(read, entry, strict, genomes, transcriptome):
    """First source matching within `strict`, else the closest one."""
    best = None
    for kind, reference, _ in entry["sources"]:
        hit = source_footprint(read, kind, reference, genomes, transcriptome,
                               strict)
        if hit is None:
            continue
        if best is None or hit[0] < best[1][0]:
            best = ((kind, reference), hit)
        if hit[0] <= strict:
            break
    return best


def within_insert_window(anchor, query, seq, mean, sd, read_length, loose,
                         downstream):
    """Does the mate's forward-strand `query` sit where the anchor says it must?

    The anchor fixes the fragment's place and strand, so its mate can only
    appear inside the window spanned by the plausible insert lengths. Searching
    that window alone is cheap and immune to repeats elsewhere in the assembly.
    `downstream` is True when the mate is one insert length ahead of the anchor
    (the anchor is R1), False when it is behind it (the anchor is revcomp(R2)).
    """
    low = max(read_length + 1, mean - 3 * sd)
    high = mean + 3 * sd
    sign = 1 if (anchor[3] == "+") == downstream else -1
    offset_a = sign * (low - read_length)
    offset_b = sign * (high - read_length)
    window_start = max(0, anchor[2] + min(offset_a, offset_b))
    window_end = min(len(seq), anchor[2] + max(offset_a, offset_b) + read_length)
    if window_end - window_start < read_length:
        return False
    hit = find_approx(query, seq[window_start:window_end], loose)
    return hit is not None and hit[0] <= loose


def find_all(query, seq, accept, limit):
    """Up to `limit` placements of `query` in `seq` with at most `accept` errors."""
    hits = []
    seen = set()
    length = len(query)
    for offset in offsets(length):
        window = query[offset:offset + WINDOW]
        position = seq.find(window)
        while position != -1:
            start = position - offset
            if start >= 0 and start + length <= len(seq) and start not in seen:
                seen.add(start)
                distance = mismatches(query, seq, start)
                if distance <= accept:
                    hits.append((distance, start))
                    if len(hits) >= limit:
                        return hits
            position = seq.find(window, position + 1)
    return hits


def source_candidates(read, reference, genomes, accept, limit):
    """Best-first placements of `read` in one genome, repeated elements included."""
    hits = []
    for contig, seq in genomes[reference]["contigs"]:
        for strand, query in (("+", read), ("-", revcomp(read))):
            for distance, start in find_all(query, seq, accept, limit):
                hits.append((distance, contig, start, strand))
    hits.sort()
    return hits[:limit]


def mate_in_window(anchor, forward_query, reverse_query, downstream, reference,
                   mean, sd, args, genomes, transcriptome, loose):
    """Is the mate of `anchor` where the fragment geometry puts it?"""
    query = forward_query if anchor[3] == "+" else reverse_query
    seq = sequence_of("genome", reference, anchor[1], genomes, transcriptome)
    return within_insert_window(anchor, query, seq, mean, sd, args.read_length,
                                loose, downstream)


def genome_pair_consistent(r1, r2, reference, mean, sd, args, genomes,
                           transcriptome, loose):
    """Anchor each mate in turn and check the other against the expected locus.

    The unique placement is tried first; only if that fails are all placements
    collected, because a read inside a repeated element sits in several places
    and the first one found need not be the fragment's own.
    """
    folded = revcomp(r2)
    for anchor_read, forward_query, reverse_query, downstream in (
            (r1, folded, r2, True), (folded, r1, revcomp(r1), False)):
        geometry = (forward_query, reverse_query, downstream, reference, mean, sd,
                    args, genomes, transcriptome, loose)
        best = source_footprint(anchor_read, "genome", reference, genomes,
                                transcriptome, loose)
        if best is not None and best[0] <= loose and \
                mate_in_window(best, *geometry):
            return True
        for candidate in source_candidates(anchor_read, reference, genomes, loose,
                                           MAX_ANCHORS):
            if candidate != best and mate_in_window(candidate, *geometry):
                return True
    return False


def transcript_pair_consistent(r1, r2, transcriptome, mean, sd, args, loose):
    """Both mates must sit in one transcript, a plausible distance apart.

    Searching transcript by transcript: a fragment comes from a single
    transcript, and isoforms that share exons make a flat anchor ambiguous.
    """
    folded = revcomp(r2)
    expected = mean - args.read_length
    for _, seq in transcriptome["eligible"]:
        first = footprint(r1, seq, loose)
        if first is None or first[0] > loose:
            continue
        second = footprint(folded, seq, loose)
        if second is None or second[0] > loose:
            continue
        if first[2] != second[2]:
            continue
        if abs(abs(first[1] - second[1]) - expected) > 3 * sd + 10:
            continue
        return True
    return False


def pair_consistent(r1, r2, entry, mean, sd, args, genomes, transcriptome, loose):
    """Is this pair consistent with any of the sample's sources?"""
    for kind, reference, _ in entry["sources"]:
        if kind == "transcript":
            if transcript_pair_consistent(r1, r2, transcriptome, mean, sd, args,
                                          loose):
                return True
        elif genome_pair_consistent(r1, r2, reference, mean, sd, args, genomes,
                                    transcriptome, loose):
            return True
    return False


def iter_fastq(path):
    with gzip.open(path, "rt") as handle:
        while True:
            name = handle.readline()
            if not name:
                return
            seq = handle.readline().rstrip("\n")
            handle.readline()
            quality = handle.readline().rstrip("\n")
            yield name[1:].rstrip("\n"), seq, quality


def check_sample(entry, args, genomes, transcriptome):
    """Verify one sample's fastq; return (ok, summary)."""
    paths = [os.path.join(args.outdir, name) for name in entry["files"]]
    for path in paths:
        if not os.path.exists(path):
            die("missing %s; generate the test reads first" % path)
    records = [list(iter_fastq(path)) for path in paths]
    counts = [len(block) for block in records]
    references = [reference for _, reference, _ in entry["sources"]]

    if len(set(counts)) != 1:
        return False, "record counts differ across mates: %s" % counts
    if counts[0] != args.reads:
        return False, ("holds %d reads, expected %d (regenerate with --reads %d)"
                       % (counts[0], args.reads, args.reads))
    for index in range(counts[0]):
        first = records[0][index]
        if len(first[1]) != args.read_length or len(first[2]) != args.read_length:
            return False, "read %d has length %d, expected %d" % (
                index + 1, len(first[1]), args.read_length)
        if entry["layout"] != "pe":
            continue
        second = records[1][index]
        if first[0].rsplit("/", 1)[0] != second[0].rsplit("/", 1)[0] or \
                not first[0].endswith("/1") or not second[0].endswith("/2"):
            return False, "mates do not pair at record %d: %s / %s" % (
                index + 1, first[0], second[0])

    # Reads carry Binomial(L, rate) errors, so accept a 5-sigma bound; the strict
    # bound only decides which source is credited when a sample mixes two.
    errors = args.error_rate * args.read_length
    strict = max(2, int(math.ceil(errors)) + 2)
    loose = max(strict, int(math.ceil(errors + 5 * math.sqrt(errors))))
    step = max(1, counts[0] // max(args.check_reads, 1))
    tally = defaultdict(int)
    worst = 0
    checked = 0
    for index in range(0, counts[0], step):
        winner = best_source(records[0][index][1], entry, strict, genomes,
                             transcriptome)
        if winner is None or winner[1][0] > loose:
            distance = "none" if winner is None else winner[1][0]
            return False, ("read %d matches none of %s (best mismatches: %s)"
                           % (index + 1, "/".join(references), distance))
        tally[winner[0][1]] += 1
        worst = max(worst, winner[1][0])
        checked += 1

    if checked == 0:
        return False, "no reads in %s" % ", ".join(entry["files"])
    if len(entry["sources"]) > 1:
        for reference in references:
            share = tally[reference] / float(checked)
            if share < MIN_SOURCE_SHARE:
                return False, ("%s contributes only %.0f%% of the checked reads"
                               % (reference, share * 100))

    pairs_checked = 0
    pairs_bad = []
    if entry["layout"] == "pe":
        mean = args.rna_insert_mean if entry["assay"] == "rnaseq" \
            else args.insert_mean
        sd = args.rna_insert_sd if entry["assay"] == "rnaseq" else args.insert_sd
        for index in range(0, counts[0], max(1, counts[0] // args.pair_checks)):
            if pair_consistent(records[0][index][1], records[1][index][1], entry,
                               mean, sd, args, genomes, transcriptome, loose):
                pairs_checked += 1
            else:
                pairs_bad.append(index + 1)
        if pairs_bad:
            return False, ("mate geometry is wrong for pair(s) %s (the mates "
                           "should sit one insert length apart, on one strand)"
                           % ", ".join(str(number) for number in pairs_bad[:5]))

    summary = "%s: %d %s; %d checked (%s); max mismatches %d" % (
        ", ".join(entry["files"]), counts[0],
        "pairs" if entry["layout"] == "pe" else "reads", checked,
        ", ".join("%s=%d" % (reference, tally[reference])
                  for reference in references), worst)
    if entry["layout"] == "pe":
        summary += "; %d pairs geometry-checked" % pairs_checked
    return True, summary


def check(args):
    genomes = load_genomes(args.ref_a, args.ref_b)
    transcriptome = load_transcriptome(args.gtf, genomes["ref_a"]["contigs"],
                                       args.read_length + 1)
    failed = False
    for entry in SAMPLES:
        ok, summary = check_sample(entry, args, genomes, transcriptome)
        print("[%s] %s" % ("ok" if ok else "FAIL", summary))
        failed = failed or not ok
    if failed:
        die("the test reads do not match the references")
    print("all test reads match the references")


# --------------------------------------------------------------------- main


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Simulate (or verify) the hybseq test reads from "
                    "test/data/ref_a.fa, ref_b.fa and ref_a.gtf.")
    parser.add_argument("--reads", type=int, default=DEFAULT_READS,
                        help="read pairs per paired-end sample, and reads per "
                             "single-end sample (default: %(default)s)")
    parser.add_argument("--read-length", type=int, default=DEFAULT_READ_LENGTH,
                        help="read length in bp (default: %(default)s)")
    parser.add_argument("--insert-mean", type=int, default=DEFAULT_INSERT_MEAN,
                        help="mean DNA fragment length (default: %(default)s)")
    parser.add_argument("--insert-sd", type=int, default=DEFAULT_INSERT_SD,
                        help="DNA fragment length SD (default: %(default)s)")
    parser.add_argument("--rna-insert-mean", type=int,
                        default=DEFAULT_RNA_INSERT_MEAN,
                        help="mean RNA fragment length (default: %(default)s)")
    parser.add_argument("--rna-insert-sd", type=int,
                        default=DEFAULT_RNA_INSERT_SD,
                        help="RNA fragment length SD (default: %(default)s)")
    parser.add_argument("--error-rate", type=float, default=DEFAULT_ERROR_RATE,
                        help="per-base substitution rate (default: %(default)s)")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED,
                        help="master seed (default: %(default)s)")
    parser.add_argument("--ref-a", default=os.path.join(DATA_DIR, "ref_a.fa"),
                        help="reference FASTA (default: %(default)s)")
    parser.add_argument("--ref-b", default=os.path.join(DATA_DIR, "ref_b.fa"),
                        help="second reference FASTA (default: %(default)s)")
    parser.add_argument("--gtf", default=os.path.join(DATA_DIR, "ref_a.gtf"),
                        help="annotation for the RNA samples (default: "
                             "%(default)s)")
    parser.add_argument("--outdir", default=FASTQ_DIR,
                        help="where the fastq live (default: %(default)s)")
    parser.add_argument("--check", action="store_true",
                        help="verify the fastq in --outdir instead of writing "
                             "them; uses the same --reads/--read-length")
    parser.add_argument("--check-reads", type=int, default=DEFAULT_CHECK_READS,
                        help="reads verified per file with --check "
                             "(default: %(default)s)")
    parser.add_argument("--pair-checks", type=int, default=DEFAULT_PAIR_CHECKS,
                        help="mate pairs geometry-checked per paired-end file "
                             "with --check (default: %(default)s)")
    args = parser.parse_args(argv)
    if args.reads < 1:
        die("--reads must be at least 1")
    if args.read_length <= WINDOW:
        die("--read-length must be greater than %d" % WINDOW)
    if args.error_rate < 0:
        die("--error-rate must not be negative")
    return args


def main(argv):
    args = parse_args(argv)
    if args.check:
        check(args)
    else:
        generate(args)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
