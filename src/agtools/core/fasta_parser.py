import gzip
import warnings

from Bio.Seq import Seq


class FastaParser:
    """
    A minimal, lightweight FASTA parser with on-demand sequence retrieval.

    This parser builds an index mapping of sequence names to byte
    offsets in the file, allowing sequences to be fetched lazily
    without loading the entire FASTA file into memory.
    Works with both plain-text FASTA and gzip-compressed FASTA (.gz).

    Attributes
    ----------
    file_path : str
        Path to the FASTA file (plain or gzipped).
    assembler : str
        Assembler used to get the GFA file.
    mapping : dict
        Name mapping of contigs in graph and FASTA file (MEGAHIT).
    index : dict
        Mapping of sequence name to file offset for the header line.
    gzipped : bool
        True if the file is gzip-compressed.

    Methods
    -------
    get_sequence(seq_name)
        Retrieve a DNA sequence by sequence name.
    get_index(seq_name)
        Retrieve the file pointer of the DNA sequence by sequence name.

    Examples
    --------
    >>> from agtools.core.fasta_parser import FastaParser
    >>> parser = FastaParser("contigs.fasta")
    """

    def __init__(self, file_path, assembler="general", mapping=None):
        self.file_path = file_path
        self.assembler = assembler
        self.mapping = mapping  # for MEGAHIT
        self.index = {}
        self.gzipped = str(file_path).endswith(".gz")
        self._build_index()

    def _open(self, mode="rt"):
        """
        Open the FASTA file in text mode, supporting gzip if needed.

        Parameters
        ----------
        mode : str, optional
            File mode, default is 'rt' (read text).

        Returns
        -------
        file object
            An open file handle.
        """
        if self.gzipped:
            return gzip.open(self.file_path, mode)
        return open(self.file_path, mode)

    def _build_index(self):
        """
        Build an index mapping sequence names to byte offsets.

        For each header line starting with '>', store the current file position.
        This allows seeking to the start of a sequence later.
        """
        with self._open("rt") as f:
            pos = f.tell() if not self.gzipped else f.fileobj.tell()
            line = f.readline()
            while line:
                if line.startswith(">"):
                    if self.assembler == "myloasm":
                        seq_name = line[1:].strip().split("_")[0]
                    else:
                        seq_name = line[1:].strip().split()[0]
                    self.index[seq_name] = pos
                pos = f.tell() if not self.gzipped else f.fileobj.tell()
                line = f.readline()

    def get_sequence(self, seq_name: str) -> Seq:
        """
        Retrieve a DNA sequence by sequence name.

        Parameters
        ----------
        seq_name : str
            The sequence name to fetch (matching the FASTA header without '>').

        Returns
        -------
        Bio.Seq.Seq
            The DNA sequence corresponding to the given sequence name.

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
        Seq('TGGCTCTTCA')
        """
        seq_name = self.mapping[seq_name] if self.assembler == "megahit" else seq_name
        if seq_name not in self.index:
            warnings.warn(
                f"The sequence {seq_name} is not found in the contigs FASTA file",
                RuntimeWarning,
            )
            return ""

        seq_lines = []

        with self._open("rt") as f:
            if not self.gzipped:
                f.seek(self.index[seq_name])
            else:
                # For gzip, use fileobj.seek
                f.fileobj.seek(self.index[seq_name])

            f.readline()  # skip header line
            for line in f:
                if line.startswith(">"):
                    break
                seq_lines.append(line.strip())

        return Seq("".join(seq_lines))

    def get_index(self, seq_name: str) -> int:
        """
        Retrieve the file pointer of the DNA sequence by sequence name.

        Parameters
        ----------
        seq_name : str
            The sequence name to fetch (matching the FASTA header without '>').

        Returns
        -------
        int
            The DNA sequence as a string, or None if the name is not found.

        Raises
        ------
        KeyError
            If the sequence is not found in the index

        Examples
        --------
        >>> parser.get_index("contig_1")
        8487228
        """
        if seq_name in self.index:
            return self.index[seq_name]
        else:
            raise KeyError(f"{seq_name} not found in the index")
