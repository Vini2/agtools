#!/usr/bin/env python3

import os
import sys
import warnings
from collections.abc import Mapping

from rsbio_seq import SeqReader, SeqReaderIndexed, build_index

from agtools.core._io import is_bgzf, is_gzip, open_text
from agtools.log_config import logger

# Extensions rsbio-seq recognises when inferring the file format.
_FASTA_EXTENSIONS = (".fasta", ".fa", ".fna")


def _looks_like_fasta(file_path: str) -> bool:
    """Check whether rsbio-seq will recognise this path as FASTA."""
    name = str(file_path).lower()
    if name.endswith(".gz"):
        name = name[:-3]
    return name.endswith(_FASTA_EXTENSIONS)


class FaiIndex(Mapping):
    """
    A read-only, memory-light mapping of sequence name to file offset.

    This behaves like the plain ``dict`` that :class:`FastaParser` used to
    build, but it does not keep one Python string and one dict slot per
    record. On an assembly with a few million contigs that difference is
    hundreds of megabytes.

    Instead:

    * iteration, lengths and record counts stream straight off the ``.fai``
      file and retain nothing;
    * the name-to-offset table is materialised only if something actually
      subscripts the mapping or tests membership -- and then only once.

    Parameters
    ----------
    fai_path : str
        Path to the ``.fai`` index file.
    key_for : callable
        Maps a FASTA record ID to the name ``agtools`` looks records up by.
    identity_keys : bool
        True when ``key_for`` is the identity, in which case a lookup key is
        already a FASTA record ID and no name table is ever needed.
    """

    __slots__ = ("_fai_path", "_key_for", "_identity", "_offsets", "_count", "_names")

    def __init__(self, fai_path, key_for, identity_keys):
        self._fai_path = fai_path
        self._key_for = key_for
        self._identity = identity_keys
        self._offsets = None
        self._names = None
        self._count = None

    def _records(self):
        """Stream ``(key, offset)`` pairs off the .fai without retaining them."""
        with open(self._fai_path) as handle:
            for line in handle:
                if not line.strip():
                    continue
                fields = line.split("\t")
                yield self._key_for(fields[0]), int(fields[2])

    def _offset_table(self) -> dict:
        """Materialise the full offset table, once, on first subscript."""
        if self._offsets is None:
            self._offsets = {key: offset for key, offset in self._records()}
        return self._offsets

    def record_id(self, key: str) -> str:
        """
        Resolve a lookup key back to the FASTA record ID it came from.

        For most assemblers the two are the same and this is free. Only
        assemblers whose keys are derived from the record ID (myloasm) pay
        for a name table here.
        """
        if self._identity:
            return key
        return self._key_table().get(key, key)

    def _key_table(self) -> dict:
        """Name table for assemblers whose keys differ from record IDs."""
        if self._names is None:
            names = {}
            with open(self._fai_path) as handle:
                for line in handle:
                    if not line.strip():
                        continue
                    record_id = line.split("\t", 1)[0]
                    names[self._key_for(record_id)] = record_id
            self._names = names
        return self._names

    def __getitem__(self, key):
        return self._offset_table()[key]

    def __contains__(self, key):
        if self._identity:
            return key in self._offset_table()
        return key in self._key_table()

    def __iter__(self):
        if self._identity:
            seen = None
        else:
            seen = set()
        for key, _ in self._records():
            if seen is not None:
                if key in seen:
                    continue
                seen.add(key)
            yield key

    def __len__(self):
        if self._count is None:
            self._count = sum(1 for _ in self._records())
        return self._count

    def iter_lengths(self):
        """
        Stream ``(key, sequence length)`` pairs straight from the .fai.

        Lengths are recorded in the index, so this answers "how long is
        every contig" without reading a single base of sequence.
        """
        with open(self._fai_path) as handle:
            for line in handle:
                if not line.strip():
                    continue
                fields = line.split("\t")
                yield self._key_for(fields[0]), int(fields[1])


class FastaParser:
    """
    A lightweight FASTA parser with on-demand sequence retrieval.

    Sequences are fetched lazily, so the whole FASTA file is never held in
    memory. Two retrieval strategies are used, picked automatically:

    * **Indexed (fast path).** For plain-text and BGZF-compressed FASTA
      files, lookups are served by `rsbio-seq
      <https://github.com/anuradhawick/rsbio-seq>`_ through a samtools-style
      ``.fai`` index (plus a ``.gzi`` index for BGZF). The index is built
      once on first use and reused on later runs, so repeated analyses over
      the same assembly skip the indexing pass entirely. Record names,
      offsets and lengths stay on the Rust side rather than being copied
      into Python dictionaries.
    * **Scanning (fallback).** For ordinary gzip files -- which cannot be
      randomly accessed -- and for FASTA files with unrecognised extensions,
      a byte-offset table is built by scanning the headers, and sequences
      are read by seeking. Recompressing such files with ``bgzip`` instead
      of ``gzip`` moves them onto the fast path.

    Attributes
    ----------
    file_path : str
        Path to the FASTA file (plain, gzip or BGZF).
    assembler : str
        Assembler used to get the GFA file.
    mapping : dict
        Name mapping of contigs in graph and FASTA file (MEGAHIT).
    index : Mapping
        Mapping of sequence name to byte offset of the record in the file.
        On the indexed path this is a lazy :class:`FaiIndex` view rather
        than a materialised dict.
    gzipped : bool
        True if the file is gzip-compressed (ordinary gzip or BGZF).
    bgzf : bool
        True if the file is BGZF-compressed and therefore randomly accessible.
    indexed : bool
        True if lookups are served from an ``.fai`` index via rsbio-seq.

    Methods
    -------
    get_sequence(seq_name)
        Retrieve a DNA sequence by sequence name.
    get_index(seq_name)
        Retrieve the file offset of the DNA sequence by sequence name.
    iter_sequences()
        Stream every record in file order.
    iter_lengths()
        Stream every record's length without reading sequence data.

    Examples
    --------
    >>> from agtools.core.fasta_parser import FastaParser
    >>> parser = FastaParser("contigs.fasta")
    >>> parser.indexed
    True
    """

    def __init__(self, file_path, assembler="general", mapping=None):
        self.file_path = str(file_path)
        self.assembler = assembler
        self.mapping = mapping  # for MEGAHIT
        self.index = {}
        self.gzipped = is_gzip(self.file_path)
        self.bgzf = self.gzipped and is_bgzf(self.file_path)
        self.indexed = False

        self._reader = None

        if self.gzipped and not self.bgzf:
            # Ordinary gzip is a single compressed stream: there is no way to
            # jump to a record without decompressing everything before it, so
            # no index can help and lookups stay linear in file size.
            logger.warning(
                f"{self.file_path} is gzip-compressed, which cannot be indexed "
                "for random access. Sequence lookups will be slow. Recompress "
                "with bgzip (BGZF) for indexed access: "
                f"gunzip -c {self.file_path} | bgzip > {self.file_path}.tmp && "
                f"mv {self.file_path}.tmp {self.file_path}"
            )
        elif _looks_like_fasta(self.file_path):
            self.indexed = self._build_rsbio_index()

        if not self.indexed:
            self._build_scan_index()

    @property
    def reader(self):
        """
        The indexed rsbio-seq reader, opened on first use.

        Opening it loads the index on the Rust side, which for a few million
        contigs costs a few hundred megabytes. Most of what ``agtools`` does
        with a contigs file -- building a graph, counting contigs, computing
        lengths, N50 or membership -- never needs a single base of sequence,
        so the reader is not opened until a sequence is actually requested.
        """
        if self._reader is None:
            fai_path, gzi_path = self._index_paths()
            self._reader = (
                SeqReaderIndexed(self.file_path, fai_path, gzi_path)
                if self.bgzf
                else SeqReaderIndexed(self.file_path, fai_path)
            )
        return self._reader

    @property
    def _keys_are_record_ids(self) -> bool:
        """True when a lookup key is just the FASTA record ID."""
        return self.assembler != "myloasm"

    def _key_for(self, record_id: str) -> str:
        """Derive the lookup key ``agtools`` uses for a FASTA record ID."""
        if self.assembler == "myloasm":
            return record_id.split("_")[0]
        return record_id.split()[0] if record_id else record_id

    def lookup_key(self, seq_name: str) -> str:
        """
        Translate a graph-side name into the name used in the FASTA file.

        Only MEGAHIT needs a translation; for every other assembler the two
        names are the same.

        Parameters
        ----------
        seq_name : str
            Name as it appears in the assembly graph.

        Returns
        -------
        str
            Name as it appears in the FASTA file.
        """
        return self.mapping[seq_name] if self.assembler == "megahit" else seq_name

    def _index_paths(self) -> tuple[str, str]:
        """Return the ``.fai`` and ``.gzi`` companion index paths."""
        return f"{self.file_path}.fai", f"{self.file_path}.gzi"

    def _index_is_current(self) -> bool:
        """Check whether an existing index is present and not stale."""
        fai_path, gzi_path = self._index_paths()
        if not os.path.exists(fai_path):
            return False
        if self.bgzf and not os.path.exists(gzi_path):
            return False

        source_mtime = os.path.getmtime(self.file_path)
        if os.path.getmtime(fai_path) < source_mtime:
            return False
        if self.bgzf and os.path.getmtime(gzi_path) < source_mtime:
            return False

        return True

    def _build_rsbio_index(self) -> bool:
        """
        Set up indexed access through rsbio-seq.

        Builds the ``.fai`` (and ``.gzi`` for BGZF) index if it is missing or
        older than the FASTA file, then opens an indexed reader over it. No
        per-record Python objects are created here -- the index stays on the
        Rust side and :class:`FaiIndex` reads from it on demand.

        Returns
        -------
        bool
            True if indexed access is available, False if the caller should
            fall back to scanning (e.g. duplicate record IDs, or a directory
            that cannot be written to).
        """
        fai_path, gzi_path = self._index_paths()

        try:
            if not self._index_is_current():
                build_index(self.file_path)
            if not os.path.exists(fai_path):
                raise OSError(f"{fai_path} was not created")
            if self.bgzf and not os.path.exists(gzi_path):
                raise OSError(f"{gzi_path} was not created")

            self.index = FaiIndex(fai_path, self._key_for, self._keys_are_record_ids)
        except (OSError, ValueError, RuntimeError, IndexError) as error:
            warnings.warn(
                f"Could not use an index for {self.file_path} ({error}). "
                "Falling back to scanning the file.",
                RuntimeWarning,
            )
            self._reader = None
            self.index = {}
            return False

        return True

    def _build_scan_index(self):
        """
        Build a byte-offset table by scanning the FASTA headers.

        Used for ordinary gzip files, which cannot be indexed for random
        access, and for files whose extension rsbio-seq does not recognise.
        Names are interned so they are shared with the graph-side names
        rather than duplicated.
        """
        self.index = {}
        with open_text(self.file_path, "rt") as handle:
            pos = handle.tell()
            line = handle.readline()
            while line:
                if line.startswith(">"):
                    self.index[sys.intern(self._key_for(line[1:].strip()))] = pos
                pos = handle.tell()
                line = handle.readline()

    def get_sequence(self, seq_name: str) -> str:
        """
        Retrieve a DNA sequence by sequence name.

        Parameters
        ----------
        seq_name : str
            The sequence name to fetch (matching the FASTA header without '>').

        Returns
        -------
        str
            The DNA sequence corresponding to the given sequence name, or an
            empty string if the sequence is not present.

        Raises
        ------
        RuntimeWarning
            If the sequence is not found in the contigs FASTA file

        Examples
        --------
        >>> seq = parser.get_sequence("contig_1")
        >>> len(seq)
        1500
        >>> seq[:10]
        'TGGCTCTTCA'
        """
        seq_name = self.lookup_key(seq_name)

        if self.indexed:
            # Membership is answered by the Rust index, so the hot path never
            # builds a Python-side name table.
            record_id = self.index.record_id(seq_name)
            reader = self.reader
            if record_id not in reader:
                self._warn_missing(seq_name)
                return ""
            return reader[record_id].seq

        if seq_name not in self.index:
            self._warn_missing(seq_name)
            return ""

        seq_lines = []

        with open_text(self.file_path, "rt") as handle:
            handle.seek(self.index[seq_name])

            handle.readline()  # skip header line
            for line in handle:
                if line.startswith(">"):
                    break
                seq_lines.append(line.strip())

        return "".join(seq_lines)

    @staticmethod
    def _warn_missing(seq_name: str):
        warnings.warn(
            f"The sequence {seq_name} is not found in the contigs FASTA file",
            RuntimeWarning,
        )

    def get_index(self, seq_name: str) -> int:
        """
        Retrieve the file offset of the DNA sequence by sequence name.

        Parameters
        ----------
        seq_name : str
            The sequence name to fetch (matching the FASTA header without '>').

        Returns
        -------
        int
            Byte offset of the record in the FASTA file. For indexed files
            this is the offset of the first sequence base, as recorded in the
            ``.fai`` index; for scanned files it is the offset of the header
            line.

        Raises
        ------
        KeyError
            If the sequence is not found in the index

        Examples
        --------
        >>> parser.get_index("contig_1")
        8487228
        """
        try:
            return self.index[seq_name]
        except KeyError:
            raise KeyError(f"{seq_name} not found in the index") from None

    def iter_sequences(self):
        """
        Iterate over every record in the FASTA file in file order.

        This is a single sequential pass handled by rsbio-seq, which is
        faster than looking each record up by name when the whole file is
        needed, and keeps memory flat regardless of assembly size.

        Yields
        ------
        tuple[str, str]
            Pairs of ``(sequence name, sequence)``. The name follows the same
            convention as :attr:`index`.

        Examples
        --------
        >>> total = sum(len(seq) for _, seq in parser.iter_sequences())
        """
        for record in SeqReader(self.file_path):
            yield self._key_for(record.id), record.seq

    def iter_lengths(self):
        """
        Iterate over every record's length, in file order.

        On the indexed path lengths come straight from the ``.fai``, so no
        sequence data is read at all -- counting the bases in a 1 GB assembly
        becomes a scan of a few megabytes of index. On the scanning path this
        falls back to reading sequences.

        Yields
        ------
        tuple[str, int]
            Pairs of ``(sequence name, sequence length)``.

        Examples
        --------
        >>> total_bases = sum(length for _, length in parser.iter_lengths())
        """
        if self.indexed:
            yield from self.index.iter_lengths()
        else:
            for name, sequence in self.iter_sequences():
                yield name, len(sequence)

    def has_cheap_lengths(self) -> bool:
        """
        Whether :meth:`iter_lengths` avoids reading sequence data.

        True when an ``.fai`` index is in use and lookup keys are the record
        IDs, so lengths can be streamed off the index without ambiguity.
        """
        return self.indexed and self._keys_are_record_ids
